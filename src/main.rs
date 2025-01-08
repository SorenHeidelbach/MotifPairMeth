use anyhow::Ok;
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
    let args = cli::Cli::parse();
    motif_methylation_state(&args).unwrap();
}


fn motif_methylation_state(
    args: &cli::Cli,
) -> Result<(), anyhow::Error> {
    let global_timer = Instant::now();
    let motifs = match &args.motifs {
        Some(motifs) => parse_motif_strings(motifs.clone())?,
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
                builder.add_contigs(&reference);
                for chunk in chunks {
                    let contig_id = &chunk.reference;
                    info!("Processing contig: {}", contig_id);
                    if !builder.contigs.contains_key(contig_id) {
                        warn!("Contig not found in reference: {}", contig_id);
                        continue;
                    }
                    builder.push_records(chunk);
                }
                let genome_work_space = builder.build();
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
    motifs: &Vec<motif::Motif>,
) -> Result<DataFrame, anyhow::Error> {
    let mut data = Vec::new();

    for motif in motifs {
        let fwd_indices = contig.find_motif_indeces(motif.clone());
        let rev_indices = contig.find_complement_motif_indeces(motif.clone());
        let modtype = motif.mod_type;
        let fwd_records = contig.get_records(fwd_indices, Strand::Positive,modtype);
        let rev_records = contig.get_records(rev_indices, Strand::Negative,modtype);


    }
    let df = DataFrame::new(data).map_err(|e| anyhow!("Error creating DataFrame: {}", e))?;
    Ok(df)
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

fn parse_motif_strings(motif_strings: Vec<String>) -> Result<Vec<motif::Motif>, anyhow::Error> {
    motif_strings.into_iter().map(|s| parse_motif_string(s)).collect()
}










