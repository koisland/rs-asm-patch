use std::{collections::HashMap, fmt::Display, str::FromStr};

use clap::Parser;
use coitrees::{COITree, Interval, IntervalTree};
use itertools::Itertools;
use log::LevelFilter;
use simple_logger::SimpleLogger;

type RegionIntervals = HashMap<String, Vec<Interval<Option<String>>>>;

struct RegionIntervalTrees(HashMap<String, COITree<Option<String>, usize>>);

impl Display for RegionIntervalTrees {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for (region, itvs) in self.0.iter() {
            writeln!(
                f,
                "{region}: [{}]",
                itvs.iter()
                    .map(|itv| format!(
                        "({}, {}, {:?})",
                        itv.first,
                        itv.last,
                        itv.metadata.as_ref()
                    ))
                    .join(",")
            )?;
        }
        Ok(())
    }
}

mod cli;
mod concensus;
mod io;

use cli::Args;

fn main() -> eyre::Result<()> {
    let args = Args::parse();
    let log_level = LevelFilter::from_str(&args.log_level)?;
    SimpleLogger::new().with_level(log_level).init()?;

    let paf_records = io::read_paf(args.paf)?;
    let ref_roi_records = io::read_bed(args.ref_roi_bed, |_| None)?;
    let ref_misasm_records = io::read_bed(Some(args.ref_misasm_bed), |rec| Some(rec.to_owned()))?;
    let qry_misasm_records = io::read_bed(Some(args.qry_misasm_bed), |rec| Some(rec.to_owned()))?;

    let new_ctgs = concensus::get_concensus(
        &paf_records,
        ref_roi_records,
        ref_misasm_records,
        qry_misasm_records,
    )?;

    // Write fasta.
    //  // Update boundary coordinates of first and last contigs.
    //  if let Some(ctg) = formatted_rows.get_mut(0) {
    //     ctg.start = 0;
    // }
    // let last_idx = formatted_rows.len() - 1;
    // if let Some(ctg) = formatted_rows.get_mut(last_idx) {
        
    // }
    // println!("{formatted_rows:?}")
    Ok(())
}
