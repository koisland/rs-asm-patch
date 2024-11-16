use coitrees::{COITree, Interval, IntervalTree};
use std::collections::HashMap;

pub type RegionIntervals = HashMap<String, Vec<Interval<Option<String>>>>;
pub type RegionIntervalTrees = HashMap<String, COITree<Option<String>, usize>>;

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
    itree: &COITree<Option<String>, usize>,
) -> Vec<(i32, i32, Option<String>)> {
    let mut overlapping_itvs = vec![];
    itree.query(start, stop, |n| {
        overlapping_itvs.push((n.first, n.last, n.metadata.to_owned()));
    });
    overlapping_itvs
}

pub fn in_roi(
    start: i32,
    stop: i32,
    roi_intervals: Option<&COITree<Option<String>, usize>>,
) -> bool {
    if let Some(itvs) = roi_intervals.as_ref() {
        itvs.query_count(start, stop) > 0
    } else {
        true
    }
}
