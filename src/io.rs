use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

use coitrees::{COITree, Interval, IntervalTree};
use itertools::Itertools;
use paf::{PafRecord, Reader};

use super::{RegionIntervalTrees, RegionIntervals};

pub fn read_bed(
    bed: Option<impl AsRef<Path>>,
    metadata_fn: impl Fn(&str) -> Option<String>,
) -> eyre::Result<RegionIntervalTrees> {
    let mut intervals: RegionIntervals = HashMap::new();
    let mut trees: RegionIntervalTrees = RegionIntervalTrees(HashMap::new());

    let Some(bed) = bed else {
        return Ok(trees);
    };
    let bed_fh = File::open(bed)?;
    let bed_reader = BufReader::new(bed_fh);

    for line in bed_reader.lines() {
        let line = line?;

        let Some((name, start, stop, other_cols)) = line.splitn(4, '\t').collect_tuple() else {
            log::error!("Invalid line: {line}");
            continue;
        };
        let (first, last) = (start.parse::<i32>()?, stop.parse::<i32>()?);

        intervals
            .entry(name.to_owned())
            .and_modify(|intervals| {
                intervals.push(Interval {
                    first,
                    last,
                    metadata: metadata_fn(other_cols),
                })
            })
            .or_default();
    }
    for (roi, intervals) in intervals.into_iter() {
        trees.0.entry(roi).or_insert(COITree::new(&intervals));
    }
    Ok(trees)
}

pub fn read_paf(paf: impl AsRef<Path>) -> eyre::Result<Vec<PafRecord>> {
    Ok(Reader::from_path(&paf)?
        .records()
        .flatten()
        .filter(|r| r.mapping_quality() == 60 && r.strand() == '+')
        .collect_vec())
}
