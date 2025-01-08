// src/cli.rs
use clap::Parser;

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
        help = "Motifs to analyze in the format: 'ACGT_a_0'"

    )]
    pub motifs: Option<Vec<String>>,

    #[arg(
        long, 
        default_value = "5", 
        help = "Minimum coverage required to consider a position"
    )]
    pub min_cov: u32,
    
    #[arg(
        long, 
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
        long, 
        global=true, 
        help="Verbose mode"
    )]
    pub verbose: bool,
}
