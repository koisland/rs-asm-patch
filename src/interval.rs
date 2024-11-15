use std::{collections::HashMap, fmt::Display};

use coitrees::{COITree, Interval, IntervalTree};
use itertools::Itertools;

pub type RegionIntervals = HashMap<String, Vec<Interval<Option<String>>>>;

pub struct RegionIntervalTrees(pub HashMap<String, COITree<Option<String>, usize>>);

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

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ContigType {
    Target,
    Query,
}

#[derive(Debug, Clone)]
pub struct Contig {
    pub name: String,
    pub category: ContigType,
    pub start: u32,
    pub stop: u32,
    pub full_len: u32,
}

pub fn get_overlapping_intervals(
    start: i32,
    stop: i32,
    misasm_itree: &RegionIntervalTrees,
    name: &str,
) -> Option<Vec<(i32, i32, Option<String>)>> {
    let misasms = misasm_itree.0.get(name)?;
    let mut overlapping_itvs = vec![];
    misasms.query(start, stop, |n| {
        overlapping_itvs.push((n.first, n.last, n.metadata.clone()));
    });
    Some(overlapping_itvs)
}
