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
    pub category: Option<ContigType>,
    pub start: u32,
    pub stop: u32,
    pub full_len: u32,
}

impl From<Option<ContigType>> for Contig {
    fn from(category: Option<ContigType>) -> Self {
        Contig {
            name: String::new(),
            category,
            start: 0,
            stop: 0,
            full_len: 0,
        }
    }
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

pub fn get_largest_overlapping_interval(
    start: i32,
    stop: i32,
    misasm_itree: &RegionIntervalTrees,
    name: &str,
) -> Option<(i32, i32, Option<String>)> {
    get_overlapping_intervals(start, stop, misasm_itree, name)?
        .into_iter()
        .max_by(|a, b| {
            let len_a = a.1 - a.0;
            let len_b = b.1 - b.0;
            len_a.cmp(&len_b)
        })
}
