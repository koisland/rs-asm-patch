use std::collections::HashMap;

use clap::Parser;
use coitrees::{COITree, Interval};
use log::LevelFilter;
use simple_logger::SimpleLogger;

type RegionIntervals = HashMap<String, Vec<Interval<Option<String>>>>;
type RegionIntervalTrees = HashMap<String, COITree<Option<String>, usize>>;

mod cli;
mod concensus;
mod io;

use cli::Args;

fn main() -> eyre::Result<()> {
    SimpleLogger::new().with_level(LevelFilter::Debug).init()?;
    let args = Args::parse();

    let paf_records = io::read_paf(args.paf)?;
    let roi_records = io::read_bed(args.roi_bed, |_| None)?;
    let misasm_records = io::read_bed(Some(args.ref_misasm_bed), |rec| Some(rec.to_owned()))?;

    let new_ctgs = concensus::get_concensus(&paf_records, roi_records, misasm_records);

    // Write fasta.

    Ok(())
}
