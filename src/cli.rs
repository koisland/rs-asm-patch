use std::path::PathBuf;

use clap::Parser;

#[derive(Debug, Parser)]
pub struct Args {
    /// Input PAF between input_ref_fa and input_query_fa. Requires cg tag.
    #[arg(short = 'i', long, required = true)]
    pub paf: PathBuf,

    /// Input reference fasta file.
    #[arg(short = 'r', long, required = true)]
    pub ref_fa: PathBuf,

    /// Input query fasta file.
    #[arg(short = 'q', long, required = true)]
    pub query_fa: PathBuf,

    /// Reference misassembly bed.
    #[arg(long, required = true)]
    pub ref_misasm_bed: PathBuf,

    /// Query misassembly bed.
    #[arg(long, required = true)]
    pub qry_misasm_bed: PathBuf,

    /// Merge misassemblies by provided bp.
    #[arg(long)]
    pub bp_merge_misasm: Option<u32>,

    /// Extend patched region on both ends by provided bp.
    #[arg(long)]
    pub bp_extend_patch: Option<u32>,

    /// Regions to filter alignments in reference.
    #[arg(long)]
    pub ref_roi_bed: Option<PathBuf>,

    /// Output fasta file.
    #[arg(short, long)]
    pub output_fa: Option<PathBuf>,

    /// Output BED file with misassemblies.
    #[arg(short = 'b', long)]
    pub output_bed: Option<PathBuf>,

    /// Log level.
    #[arg(long, default_value = "Info")]
    pub log_level: String,
}
