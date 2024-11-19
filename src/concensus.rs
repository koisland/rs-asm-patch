use std::{collections::HashMap, fmt::Debug};

use coitrees::{COITree, GenericInterval, Interval, IntervalTree};
use eyre::bail;
use impg::paf::{PafRecord, Strand};
use itertools::Itertools;

use crate::interval::{
    in_roi, merge_overlapping_intervals, ContigType, RegionIntervalTrees, RegionIntervals,
};

fn merge_collapse_rows_by_rle_id(
    itvs: Vec<(usize, Interval<(ContigType, String)>)>,
) -> Vec<Interval<(ContigType, String)>> {
    let mut ref_itvs = vec![];
    let mut qry_itvs = vec![];
    for (i, itv) in itvs {
        let (ctype, _name) = itv.metadata();
        let itv = Interval::new(itv.first, itv.last, (i, itv.metadata.clone()));
        match ctype {
            ContigType::Target => {
                ref_itvs.push(itv);
            }
            ContigType::Query => {
                qry_itvs.push(itv);
            }
        }
    }
    let ref_itree: COITree<(usize, (ContigType, String)), usize> = COITree::new(&ref_itvs);
    let qry_itree: COITree<(usize, (ContigType, String)), usize> = COITree::new(&qry_itvs);
    // Merge overlapping intervals based on rle id.
    merge_overlapping_intervals(ref_itree.iter(), |a, _b| a.clone(), Some(|a| a))
        .into_iter()
        .chain(merge_overlapping_intervals(
            qry_itree.iter(),
            |a, _b| a.clone(),
            Some(|a| a),
        ))
        .sorted_by(|a, b| a.metadata.0.cmp(&b.metadata.0))
        .map(|a| Interval::new(a.first, a.last, a.metadata.1))
        .collect()
}

pub fn get_concensus<T: Clone + Debug>(
    paf_rows: &[PafRecord],
    paf_file: &str,
    ref_roi_itree: Option<RegionIntervalTrees<T>>,
    ref_misasm_itree: RegionIntervalTrees<T>,
    _qry_misasm_itree: RegionIntervalTrees<T>,
    bp_extend_patch: Option<u32>,
) -> eyre::Result<RegionIntervals<(ContigType, String)>> {
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
        let mut new_ctgs: Vec<Interval<(ContigType, String)>> = vec![];
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
                (ContigType::Target, rname.to_owned()),
            ));

            log::debug!("{mtype:?} at {rname}:{mstart}-{mstop}.");
            // Liftover coordinates.
            let liftover = impg.query(rid, mstart, mstop);

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
                for (qitv, _, _) in ritv_grps {
                    let (Some(qname), Some(qlen)) = (
                        impg.seq_index.get_name(qitv.metadata),
                        impg.seq_index.get_len_from_id(qitv.metadata),
                    ) else {
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
                    log::debug!(
                        "\tLifted to {qname}:{}-{} ({strand:?}).",
                        qitv.first,
                        qitv.last,
                    );
                    // Correct coordinates.
                    let (qstart, qstop) = if strand == Strand::Reverse {
                        let qlen = TryInto::<i32>::try_into(qlen)?;
                        let (adj_qstart, adj_qstop) = (qlen - qitv.first, qlen - qitv.last);
                        log::debug!("\tAdjusted coordinates to {qname}:{adj_qstart}-{adj_qstop}");
                        (adj_qstart, adj_qstop)
                    } else {
                        (qitv.first, qitv.last)
                    };
                    new_ctgs.push(Interval::new(
                        qstart,
                        qstop,
                        (ContigType::Query, qname.to_owned()),
                    ));
                }
            }
        }

        // Add remainder of reference contig.
        new_ctgs.push(Interval::new(
            rstart,
            rlen.try_into()?,
            (ContigType::Target, rname.to_owned()),
        ));

        let mut rle_id = 0;
        let mut collapsed_rows: Vec<(usize, Interval<(ContigType, String)>)> = vec![];
        // Group by ctg_name and ctg_type and find min start and max stop coordinate.
        for (_, ctg_grp_rows) in &new_ctgs.iter().chunk_by(|ctg| &ctg.metadata) {
            // Store id to sort later.
            rle_id += 1;

            let rows: Vec<&Interval<(ContigType, String)>> = ctg_grp_rows.collect_vec();
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

        final_rows
            .entry(rname.to_owned())
            .or_insert_with(|| merge_collapse_rows_by_rle_id(collapsed_rows));
    }

    Ok(final_rows)
}
