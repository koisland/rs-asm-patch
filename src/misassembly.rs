use std::{fmt::Debug, path::Path};

use coitrees::{Interval, IntervalTree};

use crate::{
    interval::{merge_overlapping_intervals, RegionIntervalTrees},
    io,
};

/// Read misassemblies and optionally merge intervals by some bp distance.
pub fn read_misassemblies(
    bedfile: impl AsRef<Path> + Debug,
    bp_merge_misasm: Option<u32>,
) -> eyre::Result<RegionIntervalTrees> {
    // Get bps to merge misasm.
    // Modify the misassemblies and then convert back.
    let misasm_interval_fn = |start, stop, other_cols: &str| {
        let (start, stop) = bp_merge_misasm
            .map(|bp| (start - bp as i32, stop + bp as i32))
            .unwrap_or_else(|| (start, stop));

        Interval::new(start, stop, Some(other_cols.to_owned()))
    };
    // Read misassemblies.
    // Always required.
    let ref_misasm_records = io::read_bed(Some(&bedfile), misasm_interval_fn)?.unwrap();

    // Merge overlaps.
    // Then reverse bp modification above.
    if let Some(bp_merge_misasm) = bp_merge_misasm {
        log::info!("Merging overlapping intervals in {bedfile:?} by {bp_merge_misasm} bp.");
        Ok(ref_misasm_records
            .into_iter()
            .map(|(name, itree)| {
                let corrected_intervals = merge_overlapping_intervals(
                    itree.iter(),
                    |a, _b| a.clone(),
                    Some(|itv: Interval<Option<String>>| {
                        Interval::new(
                            itv.first + bp_merge_misasm as i32,
                            itv.last - bp_merge_misasm as i32,
                            itv.metadata,
                        )
                    }),
                );
                (name, corrected_intervals)
            })
            .collect())
    } else {
        Ok(ref_misasm_records)
    }
}

#[cfg(test)]
mod tests {
    use coitrees::{Interval, IntervalTree};
    use itertools::Itertools;

    use crate::{interval::RegionIntervalTrees, io};

    use super::read_misassemblies;

    fn assert_region_itree_equal(
        rgn_itree_1: RegionIntervalTrees,
        rgn_itree_2: RegionIntervalTrees,
    ) {
        for ((rgn_1, rgn_1_itree), (rgn_2, rgn_2_itree)) in rgn_itree_1
            .into_iter()
            .sorted_by(|a, b| a.0.cmp(&b.0))
            .zip(rgn_itree_2.into_iter().sorted_by(|a, b| a.0.cmp(&b.0)))
        {
            assert_eq!(rgn_1, rgn_2);
            assert_eq!(
                rgn_1_itree
                    .iter()
                    .map(|itv| (itv.first, itv.last, itv.metadata))
                    .collect_vec(),
                rgn_2_itree
                    .iter()
                    .map(|itv| (itv.first, itv.last, itv.metadata))
                    .collect_vec(),
            )
        }
    }

    #[test]
    fn test_read_misassemblies() {
        let wd = std::env::current_dir().unwrap();
        let infile = wd.join("data/test/merge_overlaps/input.bed");
        let outfile = wd.join("data/test/merge_overlaps/output.bed");
        // Merge intervals with 100,000 bp distance.
        let input_misassemblies = read_misassemblies(infile, Some(100_000)).unwrap();
        let expected_misassemblies = io::read_bed(Some(outfile), |st, end, other| {
            Interval::new(st, end, Some(other.to_owned()))
        })
        .unwrap()
        .unwrap();

        assert_region_itree_equal(input_misassemblies, expected_misassemblies)
    }
}
