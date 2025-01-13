
// src/cli.rs
use clap::{Parser, ValueEnum};
/// A CLI tool that processes a file with optional numeric parameters.
#[derive(Parser, Debug)]
#[command(name = "my_cli", version, about = "An example CLI")]
pub struct Cli {
    #[arg(
        value_name = "REFERENCE", 
        help = "File path to the fasta file with references"
    )]
    pub reference: String,

    #[arg(
        value_name = "PILEUP",
        help = "File path to the pileup file with methylation data"
    )]
    pub pileup: String,

    #[arg(
        value_name = "MOTIFS",
        help = "Comeplement motif pairs in the format: 'MOTIF_TYPE1_POS1_TYPE2_POS2', e.g. 'ACGT_a_0_m_3' or 'CCWGG_4mC_0_5mC_3'"
    )]
    pub motifs: Option<Vec<String>>,

    #[arg(
        long,
        short,
        default_value = "motif_methylation_state",
        value_name = "OUT",
        help = "Output file path"
    )]
    pub out: String,

    #[arg(
        long, 
        default_value = "5", 
        help = "Minimum coverage required to consider a position"
    )]
    pub min_cov: u32,
    
    #[arg(
        long, 
        short,
        default_value = "5", 
        help = "Number of threads to use"
    )]
    pub threads: u32,

    #[arg(
        long, 
        default_value = "100", 
        help = "Number of contigs to load and process at once"
    )]
    pub batch_size: u32,

    #[arg(
        value_enum,
        long,
        default_value = "normal",
        value_name = "VERBOSITY",
        help = "Verbosity level"
    )]
    pub verbosity: LogLevel,
}

#[derive(ValueEnum, Clone, Debug)]
pub enum LogLevel {
    verbose,
    normal,
    silent
}