use coitrees::{COITree, Interval, IntervalTree};
use impg::paf::Strand;
use itertools::Itertools;
use std::collections::{HashMap, VecDeque};

pub type RegionIntervals<T> = HashMap<String, Vec<Interval<T>>>;
pub type RegionIntervalTrees<T> = HashMap<String, COITree<T, usize>>;
pub type ContigInfo = (ContigType, String, Strand);

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ContigType {
    Target,
    Query,
}

#[allow(unused)]
/// Get overlapping intervals in a [`COITree`]. Simply wraps [`COITree::query`].
///
/// # Arguments
/// * `start`: Start position.
/// * `stop`: Stop position.
/// * `itree`: Interval tree to check for overlaps.
///
/// # Returns
/// * [`Vec`] of intervals as a tuple.
///     * [`COITree`]'s `IntervalNode` doesn't impl [`Debug`].
pub fn get_overlapping_intervals<T>(
    start: i32,
    stop: i32,
    itree: &COITree<T, usize>,
) -> Vec<(i32, i32, T)>
where
    T: Clone,
{
    let mut overlapping_itvs = vec![];
    itree.query(start, stop, |n| {
        overlapping_itvs.push((n.first, n.last, n.metadata.clone()));
    });
    overlapping_itvs
}

/// Merge intervals in a [`COITree`]. Includes book-ended intervals.
///
/// # Arguments
/// * `intervals`: Intervals to merge. Elements are cloned.
/// * `data_reducer`: Function to reduce metadata.
/// * `data_finalizer`: Function to apply some final operation on intervals.
///
/// # Returns
/// * Merged overlapping intervals.
pub fn merge_overlapping_intervals<I, T>(
    intervals: I,
    data_reducer: impl Fn(&T, &T) -> T,
    data_finalizer: Option<impl Fn(Interval<T>) -> Interval<T>>,
) -> Vec<Interval<T>>
where
    I: Iterator<Item = Interval<T>>,
    T: Clone,
{
    // let mut merged: Vec<Interval<T>> = Vec::with_capacity(intervals.len());
    let mut merged: Vec<Interval<T>> = Vec::new();
    let mut intervals: VecDeque<Interval<T>> = intervals
        .into_iter()
        .sorted_by(|a, b| a.first.cmp(&b.first))
        .collect();
    while !intervals.is_empty() {
        let Some(itv_1) = intervals.pop_front() else {
            unreachable!()
        };
        let Some(itv_2) = intervals.pop_front() else {
            merged.push(itv_1);
            break;
        };
        // (if) First case:
        // 1-2
        //     3-4
        // (else) Second case:
        // 1-2
        //   2-3
        // (else) Third case:
        // 1-2
        // 1-2
        if itv_1.last < itv_2.first {
            merged.push(itv_1);
            intervals.push_front(itv_2);
        } else {
            let new_data = data_reducer(&itv_1.metadata, &itv_2.metadata);
            let merged_interval = Interval::new(itv_1.first, itv_2.last, new_data);
            intervals.push_front(merged_interval);
        }
    }
    // Apply finalizer function
    if let Some(data_finalizer) = data_finalizer {
        merged.into_iter().map(data_finalizer).collect_vec()
    } else {
        merged
    }
}

pub fn in_roi<T: Clone>(start: i32, stop: i32, roi_intervals: Option<&COITree<T, usize>>) -> bool {
    if let Some(itvs) = roi_intervals.as_ref() {
        itvs.query_count(start, stop) > 0
    } else {
        true
    }
}

#[cfg(test)]
mod tests {
    use std::fmt::Debug;

    use coitrees::Interval;

    use super::merge_overlapping_intervals;

    fn reduce_to_a<'a>(a: &usize, _b: &usize) -> usize {
        *a
    }

    fn noop(a: Interval<usize>) -> Interval<usize> {
        a
    }

    fn assert_itvs_equal<T: Clone + PartialEq + Debug>(
        itvs_1: &[Interval<T>],
        itvs_2: &[Interval<T>],
    ) {
        itertools::assert_equal(
            itvs_1
                .iter()
                .map(|itv| (itv.first, itv.last, itv.metadata.clone())),
            itvs_2
                .iter()
                .map(|itv| (itv.first, itv.last, itv.metadata.clone())),
        );
    }

    #[test]
    fn test_no_merge_intervals() {
        let itvs = vec![
            Interval::new(1, 2, 1),
            Interval::new(3, 5, 2),
            Interval::new(6, 9, 3),
        ];
        let merged_itvs =
            merge_overlapping_intervals(itvs.clone().into_iter(), reduce_to_a, Some(noop));
        assert_itvs_equal(&itvs, &merged_itvs);
    }

    #[test]
    fn test_merge_intervals_single() {
        let itvs = vec![
            Interval::new(1, 3, 1),
            Interval::new(3, 5, 2),
            Interval::new(6, 9, 3),
        ];
        let merged_itvs = merge_overlapping_intervals(itvs.into_iter(), reduce_to_a, Some(noop));
        let exp_itvs = vec![Interval::new(1, 5, 1), Interval::new(6, 9, 3)];

        assert_itvs_equal(&exp_itvs, &merged_itvs);
    }

    #[test]
    fn test_merge_intervals_multiple() {
        let itvs = vec![
            Interval::new(1, 3, 1),
            Interval::new(6, 9, 3),
            Interval::new(3, 6, 2),
        ];
        let merged_itvs = merge_overlapping_intervals(itvs.into_iter(), reduce_to_a, Some(noop));
        let exp_itvs = vec![Interval::new(1, 9, 1)];
        assert_itvs_equal(&exp_itvs, &merged_itvs);
    }
}
