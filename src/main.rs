use anyhow::Ok;
use data::GenomeWorkspace;
use motif_methylation_state::utils::strand::Strand;
use motif_methylation_state::utils::{iupac, modtype, motif};
use anyhow::{Result, bail};
use env_logger::Env;
use log::{debug, info, log_enabled, warn, Level};
use std::fs::File;
use std::path::Path;
use clap::Parser;
use anyhow::anyhow;
use polars::prelude::*;
use polars::df;
use std::time::Instant;
mod pileup;
mod cli;
mod fasta_reader;
mod sequence;
mod data;
fn main() {
    env_logger::Builder::from_env(
        Env::default().default_filter_or("debug")
    ).init();
    info!("Running motif methylation state");
    
    std::fs::create_dir("motif_complement_methylation")
    .map_err(|e| anyhow!("Error creating output directory: {}", e));
    let args = cli::Cli::parse();
    motif_methylation_state(&args).unwrap();
}


fn motif_methylation_state(
    args: &cli::Cli,
) -> Result<(), anyhow::Error> {
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
                    builder.add_contig(reference.get(contig_id).unwrap());
                    builder.push_records(chunk);
                }
                let genome_work_space = builder.build();

                for (refenrece_id, contig) in genome_work_space.contigs.into_iter() {
                    motif_methylation_pattern(&contig, &motifs)?;
                    
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

struct MotifMethylationState {
    contig: sequence::Contig,
    start_position: u32,
    motif: Vec<motif::Motif>,
    pos_1: u32,
    mod_type_1: modtype::ModType,
    n_mod_1: u32,
    n_nomod_1: u32,
    n_diff_1: u32,
    pos_2: u32,
    mod_type_2: modtype::ModType,
    motif_2: Vec<motif::Motif>,
    n_mod_2: u32,
    n_nomod_2: u32,
    n_diff_2: u32,
}

fn motif_methylation_pattern(
    contig: &sequence::Contig,
    motifs: &Vec<motif::MotifPair>,
) -> Result<(), anyhow::Error> {
    let out_path = format!("motif_complement_methylation/{}.tsv", contig.reference);
    let mut csv_writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(out_path)?;
    csv_writer.write_record(&[
        "contig",
        "start_position",
        "motif",
        "mod_position",
        "mod_type",
        "position",
        "n_mod",
        "n_nomod",
        "n_diff",
        "mod_position_2",
        "mod_type_2",
        "position_2",
        "n_mod_2",
        "n_nomod_2",
        "n_diff_2",
    ])?;

    for motif in motifs {
        debug!("Processing motif pair: {:?}", motif);        // motif 1
        let fwd_indices = contig.find_motif_indeces(motif.forward.clone());
        debug!("Found {} forward indices for motif 1", fwd_indices.len());
        let rev_indices = contig.find_complement_motif_indeces(motif.forward.clone());
        debug!("Found {} reverse indices for motif 1", rev_indices.len());
        let modtype = motif.forward.mod_type;
        let fwd_records = contig.get_records(fwd_indices, Strand::Positive,modtype);
        let rev_records = contig.get_records(rev_indices, Strand::Negative,modtype);
        let records = {
            match (&fwd_records, &rev_records) {
                (Some(fwd), Some(rev)) => Some(fwd.iter().chain(rev.iter())),
                (Some(fwd), None) => Some(fwd.iter().chain([].iter())),
                (None, Some(rev)) => Some(rev.iter().chain([].iter())),
                (None, None) => None,
            }
        };
        if records.is_none() {
            info!("No records found for motif 1 of pair: {:?}", motif);
            continue;
        }

        // motif 2
        let fwd_indices_2 = contig.find_motif_indeces(motif.reverse.clone());
        let rev_indices_2 = contig.find_complement_motif_indeces(motif.reverse.clone());
        let modtype_2 = motif.reverse.mod_type;
        let fwd_records_2 = contig.get_records(fwd_indices_2, Strand::Positive,modtype_2);
        let rev_records_2 = contig.get_records(rev_indices_2, Strand::Negative,modtype_2);
        let records_2 = {
            match (&fwd_records_2, &rev_records_2) {
                (Some(fwd), Some(rev)) => Some(rev.iter().chain(fwd.iter())),
                (Some(fwd), None) => Some(fwd.iter().chain([].iter())),
                (None, Some(rev)) => Some(rev.iter().chain([].iter())),
                (None, None) => None,
            }
        };
        if records_2.is_none() {
            info!("No records found for motif 2 of pair: {:?}", motif);
            continue;
        }

        // combine records
        let zip_iter = records.unwrap().zip(records_2.unwrap());
        for (record, record_2) in zip_iter {
            let start_position = record.position as u32 - motif.forward.position as u32;
            let n_nomod = record.n_valid_cov - record.n_mod;
            let n_nomod_2 = record_2.n_valid_cov - record_2.n_mod;

            csv_writer.write_record(
                &[
                    record.reference.clone(),
                    start_position.to_string(),
                    motif.forward.as_string(),
                    motif.forward.position.to_string(),
                    motif.forward.mod_type.to_string().to_string(),
                    record.position.to_string(),
                    record.n_mod.to_string(),
                    n_nomod.to_string(),
                    record.n_diff.to_string(),
                    motif.reverse.position.to_string(),
                    motif.reverse.mod_type.to_string().to_string(),
                    record_2.position.to_string(),
                    record_2.n_mod.to_string(),
                    n_nomod_2.to_string(),
                    record_2.n_diff.to_string(),
                ],
            )?;

        }
    }
    csv_writer.flush()?;
    Ok(())
}

fn parse_motif_string(motif_string: String) -> Result<motif::Motif, anyhow::Error> {
    let parts: Vec<&str> = motif_string.split('_').collect();
    if parts.len() != 3 {
        bail!("Invalid motif string: {}", motif_string);
    }
    let sequence = parts[0];
    let mod_type = parts[1];
    let position = parts[2].parse::<u8>()?;
    motif::Motif::new(sequence, mod_type, position)
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

fn parse_motif_strings(motif_strings: Vec<String>) -> Result<Vec<motif::Motif>, anyhow::Error> {
    motif_strings.into_iter().map(|s| parse_motif_string(s)).collect()
}
fn parse_motif_pair_strings(motif_pair_strings: Vec<String>) -> Result<Vec<motif::MotifPair>, anyhow::Error> {
    motif_pair_strings.into_iter().map(|s| parse_motif_pair_string(s)).collect()
}









