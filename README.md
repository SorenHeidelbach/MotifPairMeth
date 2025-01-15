**Me**thylation **mo**tif **pair**

A small utility CLI for getting methylation of motif pairs. 
```
# Usage
Usage: motif_methylation_state [OPTIONS] <REFERENCE> <PILEUP> [MOTIFS]...

Arguments:
  <REFERENCE>  File path to the fasta file with references
  <PILEUP>     File path to the pileup file with methylation data
  [MOTIFS]...  Comeplement motif pairs in the format: 'MOTIF_TYPE1_POS1_TYPE2_POS2', e.g. 'ACGT_a_0_m_3' or 'CCWGG_4mC_0_5mC_3'

Options:
  -o, --out <OUT>                Output file path [default: motif_methylation_state]
      --min-cov <MIN_COV>        Minimum coverage required to consider a position [default: 5]
      --verbosity <VERBOSITY>    Verbosity level [default: normal] [possible values: verbose, normal, silent]
  -h, --help                     Print help
  -V, --version                  Print version
```
