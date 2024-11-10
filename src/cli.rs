use std::path::PathBuf;

use clap::Parser;

#[derive(Parser)]
pub struct Args {
    /// Input PAF between input_ref_fa and input_query_fa. Requires cg tag.
    #[arg(short = 'i', long)]
    pub paf: PathBuf,

    /// Input reference fasta file.
    #[arg(short = 'r', long)]
    pub ref_fa: PathBuf,

    /// Input query fasta file.
    #[arg(short = 'q', long)]
    pub query_fa: PathBuf,

    /// Reference misassembly bed.
    #[arg(long)]
    pub ref_misasm_bed: PathBuf,

    /// Regions to filter alignments in reference
    #[arg(long)]
    pub roi_bed: Option<PathBuf>,

    /// Output fasta file.
    #[arg(short, long)]
    pub output_fa: Option<PathBuf>,

    /// Output BED file with misassemblies.
    #[arg(short = 'b', long)]
    pub output_bed: Option<PathBuf>,
}
