use std::{collections::HashMap, fmt::Debug};

use coitrees::{GenericInterval, Interval};
use eyre::bail;
use impg::{
    impg::{AdjustedInterval, CigarOp},
    paf::{PafRecord, Strand},
};
use itertools::Itertools;

use crate::interval::{
    in_roi, merge_overlapping_intervals, ContigInfo, ContigType, RegionIntervalTrees,
    RegionIntervals,
};

fn merge_collapse_rows_by_rle_id<T: Clone + std::cmp::PartialEq>(
    itvs: Vec<(usize, Interval<T>)>,
) -> Vec<Interval<T>> {
    let mut final_itvs = vec![];

    // Merge overlapping intervals based on rle id.
    for (_, itv_grp) in &itvs
        .into_iter()
        .sorted_by(|a, b| a.0.cmp(&b.0))
        .chunk_by(|a| a.1.metadata.clone())
    {
        final_itvs.extend(merge_overlapping_intervals(
            itv_grp.into_iter().map(|itv| itv.1),
            |a, _b| a.clone(),
            Some(|a| a),
        ));
    }
    final_itvs
}

// Interval subset sum problem. Try to maximize the intervals and get as close as possible to target_len.
// Allow overshoot from target len.
// Based on https://stackoverflow.com/a/57526728
fn dp_itv_subset(itvs: &[AdjustedInterval], target_len: usize) -> Option<Vec<AdjustedInterval>> {
    let mut vals: Vec<Option<Vec<AdjustedInterval>>> =
        vec![None::<Vec<AdjustedInterval>>; target_len + 1];
    vals[0] = Some(vec![]);

    let mut best_val_idx = None;
    let mut best_val = None;

    // Iterate thru all intervals.
    for (qitv, cgs, ritv) in itvs {
        let itv_len = if qitv.first() < qitv.last() {
            qitv.len() as usize
        } else {
            (qitv.first() - qitv.last()) as usize
        };
        // Iterate in reverse thru all possible target lengths.
        for i in (0..target_len.saturating_sub(1)).rev() {
            let new_val_idx = i + itv_len;
            if vals[i].is_some() {
                // Update entry.
                if new_val_idx <= target_len && vals[new_val_idx].is_none() {
                    let mut new_vec = vec![];
                    new_vec.extend(vals[i].clone().unwrap());
                    new_vec.push((*qitv, cgs.to_vec(), *ritv));
                    vals[new_val_idx] = Some(new_vec)
                }
                // First, check if overshoot value.
                // And if no best value and new_val_idx close, update best value.
                if new_val_idx > target_len
                    && (best_val_idx.is_none() || Some(new_val_idx) < best_val_idx)
                {
                    best_val_idx = Some(new_val_idx);
                    let mut new_vec = vec![];
                    new_vec.extend(vals[i].clone().unwrap());
                    new_vec.push((*qitv, cgs.to_vec(), *ritv));
                    best_val = Some(new_vec)
                }
            }
        }
    }

    vals[target_len].take().or(best_val)
}

pub fn blast_identity(cigar: &[CigarOp]) -> eyre::Result<f32> {
    let (mut matches, mut mismatches, mut indels) = (0, 0, 0);
    for cg in cigar.iter() {
        match cg.op() {
            '=' | 'M' => matches += cg.len(),
            'X' => mismatches += cg.len(),
            'I' | 'D' => indels += cg.len(),
            _ => bail!("Invalid cigar op. {cg:?}"),
        }
    }
    let aln_len = matches + mismatches + indels;
    let mismatches_gaps = mismatches + indels;
    Ok((aln_len - mismatches_gaps) as f32 / aln_len as f32)
}

pub fn get_concensus<T: Clone + Debug>(
    paf_rows: &[PafRecord],
    paf_file: &str,
    ref_roi_itree: Option<RegionIntervalTrees<T>>,
    ref_misasm_itree: RegionIntervalTrees<T>,
    _qry_misasm_itree: RegionIntervalTrees<T>,
    bp_extend_patch: Option<u32>,
) -> eyre::Result<RegionIntervals<ContigInfo>> {
    let mut final_rows = HashMap::new();
    let bp_extend_patch = bp_extend_patch.unwrap_or(0) as i32;
    if bp_extend_patch != 0 {
        log::info!("Extending patched regions by {bp_extend_patch} bp.")
    }

    let impg = impg::impg::Impg::from_paf_records(paf_rows, paf_file)
        .map_err(|err| eyre::Error::msg(format!("{err:?}")))?;

    for (rname, ref_rgn_misassemblies) in ref_misasm_itree.iter() {
        let grp_roi_intervals: Option<&coitrees::BasicCOITree<T, usize>> =
            if let Some(itvs) = ref_roi_itree.as_ref() {
                itvs.get(rname)
            } else {
                None
            };
        let mut new_ctgs: Vec<Interval<ContigInfo>> = vec![];
        let Some((rid, rlen)) = impg
            .seq_index
            .get_id(rname)
            .and_then(|rid| impg.seq_index.get_len_from_id(rid).map(|rlen| (rid, rlen)))
        else {
            continue;
        };
        let mut rstart = 0;

        for misasm_itv in ref_rgn_misassemblies
            .into_iter()
            .sorted_by(|a, b| a.first.cmp(&b.first))
        {
            let (mstart, mstop, mtype) = (misasm_itv.first, misasm_itv.last, misasm_itv.metadata);
            let mlen = TryInto::<usize>::try_into(mstop - mstart)?;
            let is_roi = in_roi(mstart, mstop, grp_roi_intervals);
            if !is_roi {
                continue;
            }

            // Add more to misassembled region to extend patched region.
            let (mstart, mstop) = (
                mstart.saturating_sub(bp_extend_patch).clamp(0, i32::MAX),
                mstop + bp_extend_patch,
            );
            let correct_coords = (rstart, mstart.saturating_sub(1).clamp(0, i32::MAX));
            rstart = mstop;
            // Ignore if null interval.
            if correct_coords.0 == correct_coords.1 {
                continue;
            }
            new_ctgs.push(Interval::new(
                correct_coords.0,
                correct_coords.1,
                (ContigType::Target, rname.to_owned(), Strand::Forward),
            ));

            log::debug!("{mtype:?} at {rname}:{mstart}-{mstop} of {mlen} bp.");
            // Liftover coordinates.
            let liftover = impg.query(rid, mstart, mstop);

            // TODO: Remove and minmax.
            // Groupby rid.
            for (_, ritv_grps) in &liftover
                .into_iter()
                // Omit input interval.
                .filter(|(itv, _, _)| (itv.first, itv.last, itv.metadata) != (mstart, mstop, rid))
                .sorted_by(|(itv1, _, _), (itv2, _, _)| {
                    impg.seq_index
                        .get_name(itv1.metadata)
                        .cmp(&impg.seq_index.get_name(itv2.metadata))
                })
                .chunk_by(|(itv, _, _)| itv.metadata)
            {
                let itvs = ritv_grps.into_iter().collect_vec();
                let Some(optimal_qitvs) = dp_itv_subset(&itvs, mlen) else {
                    log::debug!("\tUnable to get a subset of {:?} whose length is greater than or equal to {mlen:?}.", itvs.iter().map(|a| a.0).collect_vec());
                    continue;
                };

                for (qitv, cg, _ritv) in optimal_qitvs
                    .into_iter()
                    .sorted_by(|a, b| a.2.first.cmp(&b.2.first))
                {
                    let Some(qname) = impg.seq_index.get_name(qitv.metadata) else {
                        bail!("Invalid query index. ({})", qitv.metadata)
                    };
                    let itv_diff = qitv.last - qitv.first;
                    let strand = if itv_diff > 0 {
                        Strand::Forward
                    } else {
                        Strand::Reverse
                    };
                    if itv_diff == 0 {
                        continue;
                    }
                    // Correct coordinates.
                    let (qstart, qstop) = if strand == Strand::Reverse {
                        (qitv.last, qitv.first)
                    } else {
                        (qitv.first, qitv.last)
                    };
                    let identity = blast_identity(&cg)?;
                    log::debug!("\tLifted to {qname}:{qstart}-{qstop} ({strand:?}) with identity {identity}.");

                    // TODO: Use ritv to sort after optimal found
                    new_ctgs.push(Interval::new(
                        qstart,
                        qstop,
                        (ContigType::Query, qname.to_owned(), strand),
                    ));
                }
            }
        }

        // Add remainder of reference contig.
        new_ctgs.push(Interval::new(
            rstart,
            rlen.try_into()?,
            (ContigType::Target, rname.to_owned(), Strand::Forward),
        ));

        let mut rle_id = 0;
        let mut collapsed_rows: Vec<(usize, Interval<ContigInfo>)> = vec![];
        // Group by ctg_name and ctg_type and find min start and max stop coordinate.
        for (_, ctg_grp_rows) in &new_ctgs.iter().chunk_by(|ctg| &ctg.metadata) {
            // Store id to sort later.
            rle_id += 1;

            let rows: Vec<&Interval<ContigInfo>> = ctg_grp_rows.collect_vec();
            let (mut min_start, mut max_stop) = (i32::MAX, 0);
            for ctg in rows.iter() {
                if ctg.first < min_start {
                    min_start = ctg.first
                }
                if ctg.last > max_stop {
                    max_stop = ctg.last
                }
            }
            if min_start > max_stop {
                continue;
            }

            // Create collapsed row.
            let new_ctg = Interval::new(min_start, max_stop, rows[0].metadata.clone());
            collapsed_rows.push((rle_id, new_ctg));
        }
        let collapsed_rows = merge_collapse_rows_by_rle_id(collapsed_rows);
        final_rows.entry(rname.to_owned()).or_insert(collapsed_rows);
    }

    Ok(final_rows)
}

#[cfg(test)]
mod tests {
    use coitrees::{GenericInterval, Interval};
    use impg::impg::AdjustedInterval;
    use itertools::Itertools;

    use super::dp_itv_subset;

    fn test_intervals_eq(itvs_1: &[AdjustedInterval], itvs_2: &[AdjustedInterval]) {
        itertools::assert_equal(
            itvs_1
                .into_iter()
                .sorted_by(|a, b| a.0.first.cmp(&b.0.first))
                .map(|a| (a.0.first, a.0.last, a.0.metadata)),
            itvs_2
                .into_iter()
                .sorted_by(|a, b| a.0.first.cmp(&b.0.first))
                .map(|a| (a.0.first, a.0.last, a.0.metadata)),
        );
    }

    #[test]
    fn test_dp_itv_subset() {
        let target_itv = Interval::new(16995741, 19569078, 1);
        let target_len = target_itv.len().try_into().unwrap();
        let null_itv = Interval::new(1, 1, 1);

        let itvs = [
            (Interval::new(0, 396846, 1), vec![], null_itv),
            (Interval::new(58, 17030, 2), vec![], null_itv),
            (Interval::new(88882357, 91081164, 3), vec![], null_itv),
        ];

        let res = dp_itv_subset(&itvs, target_len).unwrap();
        test_intervals_eq(
            &res,
            &vec![
                (
                    Interval {
                        first: 0,
                        last: 396846,
                        metadata: 1,
                    },
                    vec![],
                    null_itv,
                ),
                (
                    Interval {
                        first: 88882357,
                        last: 91081164,
                        metadata: 3,
                    },
                    vec![],
                    null_itv,
                ),
            ],
        );
    }

    #[test]
    fn test_dp_itv_subset_overshoot() {
        let target_itv = Interval::new(0, 80, 1);
        let target_len = target_itv.len().try_into().unwrap();
        let null_itv = Interval::new(1, 1, 1);
        let itvs = [
            (Interval::new(0, 50, 1), vec![], null_itv),
            (Interval::new(0, 40, 2), vec![], null_itv),
        ];
        // target itv of 80, itvs of 90 allowed.
        let res = dp_itv_subset(&itvs, target_len).unwrap();
        test_intervals_eq(&res, &itvs);
    }

    #[test]
    fn test_dp_itv_subset_none() {
        let target_itv = Interval::new(73876415, 74812999, 1);
        let target_len = target_itv.len().try_into().unwrap();
        let null_itv = Interval::new(1, 1, 1);

        let itvs = [
            (Interval::new(0, 62934, 1), vec![], null_itv),
            (Interval::new(0, 132609, 2), vec![], null_itv),
            (Interval::new(71603268, 71691116, 3), vec![], null_itv),
            (Interval::new(116337563, 116828997, 4), vec![], null_itv),
        ];
        let res = dp_itv_subset(&itvs, target_len);
        assert!(res.is_none())
    }
}
