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
    Spacer,
}

#[derive(Debug)]
pub struct Contig {
    pub name: String,
    pub category: ContigType,
    pub start: u32,
    pub stop: u32,
}

impl From<ContigType> for Contig {
    fn from(category: ContigType) -> Self {
        Contig {
            name: String::new(),
            category,
            start: 0,
            stop: 0,
        }
    }
}

pub fn get_largest_overlapping_interval(
    start: i32,
    stop: i32,
    misasm_itree: &RegionIntervalTrees,
    name: &str,
) -> Option<coitrees::IntervalNode<Option<String>, usize>> {
    let misasms = misasm_itree.0.get(name)?;
    let mut overlapping_itvs = vec![];
    misasms.query(start, stop, |n| {
        overlapping_itvs.push(n.clone());
    });
    overlapping_itvs
        .into_iter()
        .max_by(|a, b| a.len().cmp(&b.len()))
}
