use anyhow::{anyhow, bail, Result};
use bio::bio_types::strand;
use clap::Parser;
use env_logger::Env;
use log::{debug, info, log_enabled, warn, Level};
use memopair::utils::{modtype, motif, motif::MotifLike, strand::Strand};
use std::{
    collections::btree_map::Keys, f32::NAN, fs::File, path::Path, process::Output, time::Instant,
};

mod cli;
mod data;
mod fasta_reader;
mod pileup;
mod sequence;

fn main() {
    let args = cli::Cli::parse();
    // Set up logging level
    match args.verbosity {
        cli::LogLevel::silent => {
            env_logger::Builder::from_env(Env::default().default_filter_or("off")).init();
        }
        cli::LogLevel::normal => {
            env_logger::Builder::from_env(Env::default().default_filter_or("info")).init();
        }
        cli::LogLevel::verbose => {
            env_logger::Builder::from_env(Env::default().default_filter_or("debug")).init();
        }
    }

    // Create output directory
    info!("Running motif methylation state");
    let out_path = Path::new(&args.out);
    match out_path.exists() {
        true => {
            panic!("Output directory already exists");
        }
        false => match std::fs::create_dir(out_path) {
            Ok(_) => info!("Created output directory"),
            Err(e) => panic!("Could not create output directory: {}", e),
        },
    }

    // Run the main function
    match memopair(&args) {
        Ok(_) => info!("Finished running motif methylation state"),
        Err(e) => panic!("Error running motif methylation state: {}", e),
    }
}

fn memopair(args: &cli::Cli) -> Result<(), anyhow::Error> {
    let global_timer = Instant::now();
    let motifs = match &args.motifs {
        Some(motifs) => parse_motif_pair_strings(motifs.clone())?,
        None => bail!("No motifs provided"),
    };
    let reference_file = Path::new(&args.reference);
    let reference = fasta_reader::read_fasta_file(reference_file)
        .map_err(|e| anyhow::anyhow!("Error reading reference file: {}", e))?;
    info!("Loaded {} reference records", reference.len());

    let pileup_file = File::open(&args.pileup)
        .map_err(|e| anyhow::anyhow!("Could not open pileup file: {} ({})", args.pileup, e))?;
    let mut pileup_reader = pileup::PileupChunkReader::new(pileup_file, args.min_cov);

    info!("Processing pileup file: {}", args.pileup);
    loop {
        info!("Processing a batch");
        let timer = Instant::now();
        let chunks = pileup_reader.load_n_chunks(1);
        match chunks {
            Some(chunks) => {
                info!("Loaded batch {:?}", timer.elapsed());
                let mut builder = data::GenomeWorkSpaceBuilder::new();
                for chunk in chunks {
                    let contig_id = &chunk.reference;
                    info!("Processing contig: {}", contig_id);
                    debug!("Adding contig to workspace");
                    builder.add_contig(reference.get(contig_id).unwrap().clone());
                    debug!("Adding records to contig");
                    builder.push_records(chunk);
                }
                let genome_work_space = builder.build();

                for (refenrece_id, contig) in genome_work_space.contigs.into_iter() {
                    motif_methylation_pattern(&contig, &motifs, &args.out)?;
                }
            }
            None => {
                info!("Contig did not contain any records");
            }
        }

        info!("Finished batch in {:?}", timer.elapsed());

        if pileup_reader.eof_reached {
            break;
        }
    }
    info!("Finished processing in {:?}", global_timer.elapsed());
    Ok(())
}

fn motif_methylation_pattern(
    contig: &sequence::Contig,
    motifs: &Vec<motif::MotifPair>,
    out: &str,
) -> Result<(), anyhow::Error> {
    let out_path = format!("{}/{}.tsv", out, contig.reference);
    let mut record_writer = MotifPairRecordWriter::new(&out_path)?;
    record_writer.write_header()?;

    for motif in motifs {
        debug!("Processing motif pair: {:?}", motif);
        debug!("Processing forward strand");
        let mod_position_shift =
            motif.reverse.reverse_complement().unwrap().position - motif.forward.position;

        // Process forward strand
        let fwd_indices = match contig.find_motif_indeces(&motif.forward) {
            Some(i) => i,
            None => continue,
        };
        for index in fwd_indices {
            let record_1 =
                match contig
                    .records
                    .get(&(index, Strand::Positive, motif.forward.mod_type))
                {
                    Some(r) => r,
                    None => continue,
                };
            let record_2 = match contig.records.get(&(
                index + mod_position_shift as usize,
                Strand::Negative,
                motif.reverse.mod_type,
            )) {
                Some(r) => r,
                None => continue,
            };
            record_writer.write_record(motif, record_1, record_2)?;
        }
        // If motif pair is palindromic, the reverse is captured in the reverse complement of the forward motif
        if motif.is_palindromic {
            debug!("Skipping reverse strand, as motifs are palindromic");
            continue;
        }

        // Process reverse strand
        debug!("Processing reverse strand");
        let rev_indices = match contig.find_complement_motif_indeces(&motif.forward) {
            Some(i) => i,
            None => continue,
        };
        for index in rev_indices {
            let key1 = (index, Strand::Negative, motif.forward.mod_type);
            let record_1 = match contig.records.get(&key1) {
                Some(r) => r,
                None => continue,
            };
            let key2 = (
                index - mod_position_shift as usize,
                Strand::Positive,
                motif.reverse.mod_type,
            );
            let record_2 = match contig.records.get(&key2) {
                Some(r) => r,
                None => continue,
            };
            record_writer.write_record(motif, record_1, record_2)?;
        }
    }
    record_writer.flush()?;
    Ok(())
}

#[derive(Debug)]
struct MotifPairRecordWriter {
    csv_writer: csv::Writer<File>,
}

impl MotifPairRecordWriter {
    pub fn new(out_path: &str) -> Result<Self, anyhow::Error> {
        let csv_writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_path(out_path)?;
        Ok(Self { csv_writer })
    }

    pub fn write_header(&mut self) -> Result<(), anyhow::Error> {
        self.csv_writer.write_record(&[
            "contig_id",
            "motif_start_position",
            "strand",
            "motif_sequence",
            "motif_mod_position",
            "mod_type_1",
            "position_1",
            "n_mod_1",
            "n_nomod_1",
            "n_diff_1",
            "motif_mod_position_2",
            "mod_type_2",
            "position_2",
            "n_mod_2",
            "n_nomod_2",
            "n_diff_2",
            "methylation_difference",
            "odds_ratio",
            "classification",
        ])?;
        Ok(())
    }

    pub fn write_record(
        &mut self,
        motif_pair: &motif::MotifPair,
        record_1: &pileup::PileupRecord,
        record_2: &pileup::PileupRecord,
    ) -> Result<(), anyhow::Error> {
        let start_position = record_1.position as u32 - motif_pair.forward.position as u32;
        let n_nomod_1 = record_1.n_valid_cov - record_1.n_mod;
        let mean_mod_1 = record_1.n_mod as f64 / record_1.n_valid_cov as f64;
        let n_nomod_2 = record_2.n_valid_cov - record_2.n_mod;
        let mean_mod_2 = record_2.n_mod as f64 / record_2.n_valid_cov as f64;
        let motif_2_mod_pos = motif_pair.reverse.reverse_complement().unwrap().position;

        let methylation_diff = mean_mod_1 - mean_mod_2;
        let abs_methylation_diff = methylation_diff.abs();
        let odds_1 = mean_mod_1 / (1.0 - mean_mod_1);
        let odds_2 = mean_mod_2 / (1.0 - mean_mod_2);
        let odds_ratio: f64;
        if odds_2 == 0.0 || odds_1 == 0.0 {
            odds_ratio = f64::NAN;
        } else {
            odds_ratio = odds_1 / odds_2;
        }
        let classification = match abs_methylation_diff {
            x if x > 0.5 => "differential",
            x if x > 0.1 => "moderately differential",
            _ => "non-differential",
        };
        self.csv_writer.write_record(&[
            record_1.reference.clone(),
            start_position.to_string(),
            record_1.strand.to_string(),
            motif_pair.forward.sequence_string(),
            motif_pair.forward.position.to_string(),
            motif_pair.forward.mod_type.to_string().to_string(),
            record_1.position.to_string(),
            record_1.n_mod.to_string(),
            n_nomod_1.to_string(),
            record_1.n_diff.to_string(),
            motif_2_mod_pos.to_string(),
            motif_pair.reverse.mod_type.to_string().to_string(),
            record_2.position.to_string(),
            record_2.n_mod.to_string(),
            n_nomod_2.to_string(),
            record_2.n_diff.to_string(),
            abs_methylation_diff.to_string(),
            odds_ratio.to_string(),
            classification.to_string(),
        ])?;
        Ok(())
    }
    pub fn flush(&mut self) -> Result<(), anyhow::Error> {
        self.csv_writer.flush()?;
        Ok(())
    }
}

fn parse_motif_pair_string(motif_pair_string: String) -> Result<motif::MotifPair, anyhow::Error> {
    let parts: Vec<&str> = motif_pair_string.split('_').collect();
    if parts.len() != 5 {
        bail!("Invalid motif pair string: {}", motif_pair_string);
    }
    let sequence_1 = parts[0];
    let mod_type_1 = parts[1];
    let position_1 = parts[2].parse::<u8>()?;
    let motif_1 = motif::Motif::new(sequence_1, mod_type_1, position_1)?;
    let sequence_2 = motif_1.reverse_complement_sequence();
    let mod_type_2 = parts[3];
    let position_2 = parts[4].parse::<u8>()?;
    let position_2 = sequence_2.len() as u8 - position_2 - 1;
    let motif_2 = motif::Motif::new(sequence_2.as_str(), mod_type_2, position_2)?;
    let pair = motif::MotifPair::new(motif_1, motif_2)?;
    Ok(pair)
}

fn parse_motif_pair_strings(
    motif_pair_strings: Vec<String>,
) -> Result<Vec<motif::MotifPair>, anyhow::Error> {
    motif_pair_strings
        .into_iter()
        .map(|s| parse_motif_pair_string(s))
        .collect()
}
