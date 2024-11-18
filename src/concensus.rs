use std::{collections::HashMap, fmt::Debug};

use coitrees::{COITree, GenericInterval, Interval, IntervalTree};
use impg::paf::PafRecord;
use itertools::Itertools;

use crate::interval::{
    get_overlapping_intervals, in_roi, merge_overlapping_intervals, ContigType,
    RegionIntervalTrees, RegionIntervals,
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

    for (tname, pafs) in &paf_rows
        .iter()
        .sorted_by(|a, b| a.target_name.cmp(&b.target_name))
        .chunk_by(|c| c.target_name.as_str())
    {
        let grp_roi_intervals: Option<&coitrees::BasicCOITree<T, usize>> =
            if let Some(itvs) = ref_roi_itree.as_ref() {
                itvs.get(tname)
            } else {
                None
            };
        let mut new_ctgs: Vec<Interval<(ContigType, String)>> = vec![];

        let mut paf_recs = pafs
            .filter(|rec| {
                if let Some(itvs) = ref_roi_itree.as_ref() {
                    itvs.contains_key(&rec.target_name)
                } else {
                    true
                }
            })
            .sorted_by(|a, b| a.target_start.cmp(&b.target_start))
            .peekable();

        let Some(target_id) = impg.seq_index.get_id(tname) else {
            continue;
        };
        let mut tstart = 0;
        while let Some(paf_rec) = paf_recs.next() {
            let Some(query_id) = impg.seq_index.get_id(&paf_rec.query_name) else {
                continue;
            };

            // Get misassemblies
            let Some(ref_misassemblies) = ref_misasm_itree.get(&paf_rec.target_name) else {
                continue;
            };
            let ref_rgn_misassemblies = get_overlapping_intervals(
                paf_rec.target_start as i32,
                paf_rec.target_end as i32,
                ref_misassemblies,
            );
            let qname = paf_rec.query_name.to_owned();

            for (mstart, mstop, mtype) in ref_rgn_misassemblies
                .into_iter()
                .sorted_by(|a, b| a.0.cmp(&b.0))
            {
                let is_roi = in_roi(mstart, mstop, grp_roi_intervals);
                if !is_roi {
                    continue;
                }

                // Add more to misassembled region to extend patched region.
                let (mstart, mstop) = (
                    mstart.saturating_sub(bp_extend_patch).clamp(0, i32::MAX),
                    mstop + bp_extend_patch,
                );
                let correct_coords = (tstart, mstart.saturating_sub(1).clamp(0, i32::MAX));
                tstart = mstop;
                // Ignore if null interval.
                if correct_coords.0 == correct_coords.1 {
                    continue;
                }
                new_ctgs.push(Interval::new(
                    correct_coords.0,
                    correct_coords.1,
                    (ContigType::Target, tname.to_owned()),
                ));
                log::debug!("{mtype:?} at {tname}:{mstart}-{mstop}.");
                let liftover = impg.query(target_id, mstart, mstop);

                for (interval, _, _) in liftover
                    .iter()
                    .filter(|(itv, _, _)| itv.metadata == query_id)
                {
                    log::debug!(
                        "\tLifted to {}:{}-{}.",
                        impg.seq_index.get_name(interval.metadata).unwrap(),
                        interval.first,
                        interval.last
                    )
                }

                // Collapse intervals first and last position.
                let adjusted_interval = liftover
                    .iter()
                    .filter(|(itv, _, _)| itv.metadata == query_id)
                    .fold(Interval::new(i32::MAX, 0, query_id), |a, b| {
                        Interval::new(a.first.min(b.0.first), a.last.max(b.0.last), query_id)
                    });

                if adjusted_interval.last == 0 {
                    log::debug!(
                        "\tUnable to liftover {tname}:{mstart}-{mstop} ({mtype:?}) to {}.",
                        &qname
                    );
                    continue;
                };
                log::debug!(
                    "\tReplacing sequence with {}:{}-{}.",
                    &qname,
                    adjusted_interval.first,
                    adjusted_interval.last
                );

                new_ctgs.push(Interval::new(
                    adjusted_interval.first,
                    adjusted_interval.last,
                    (ContigType::Query, qname.clone()),
                ));
            }

            // Check gaps in alignment.
            let Some(_next_paf_rec) = paf_recs
                .peek()
                .filter(|rec| rec.query_name == paf_rec.query_name)
            else {
                // Add remainder of reference contig.
                new_ctgs.push(Interval::new(
                    tstart,
                    paf_rec.target_length.try_into()?,
                    (ContigType::Target, paf_rec.target_name.to_owned()),
                ));
                break;
            };

            //     // TODO: contig joining.

            //     // Case 1: Same qry ctg.
            //     // qry:     --------
            //     // alns:    ---  ---
            //     // ref:     ---??---

            //     // Case 2: Different qry ctgs
            //     // qry:     ---
            //     // qry:          ---
            //     // alns:    ---  ---
            //     // ref:     ---??---
            //     let gap_start = paf_rec.target_end;
            //     let gap_stop = next_paf_rec.target_start;

            //     let overlap_cnt = grp_roi_intervals
            //         .map(|itvs| itvs.query_count(gap_start as i32, gap_stop as i32))
            //         .unwrap_or(0);

            //     if overlap_cnt == 0 {
            //         continue;
            //     }
            //     let ref_aln_gap_misassemblies =
            //         get_overlapping_intervals(gap_start as i32, gap_stop as i32, ref_misassemblies);
            //     let qry_aln_gap_miassemblies = qry_misasm_itree
            //         .get(&paf_rec.query_name)
            //         .map_or_else(Vec::new, |itvs| {
            //             get_overlapping_intervals(gap_start as i32, gap_stop as i32, itvs)
            //         });

            //     let ref_aln_gap_misassembled = !ref_aln_gap_misassemblies.is_empty();
            //     let qry_aln_gap_misassembled = !qry_aln_gap_miassemblies.is_empty();
            //     if ref_aln_gap_misassembled || qry_aln_gap_misassembled {
            //         log::debug!(
            //             "Misassemblies between {tname}:{gap_start}-{gap_stop} with no alignment with {}",
            //             &qname
            //         );
            //         log::debug!("\tTarget: {:?}", ref_aln_gap_misassemblies,);
            //         log::debug!("\tQuery: {:?}", qry_aln_gap_miassemblies,);
            //     }

            //     match (ref_aln_gap_misassembled, qry_aln_gap_misassembled) {
            //         // Fix misassembly in target with query.
            //         (true, false) => {
            //             let qry_start = paf_rec.query_end;
            //             let qry_stop = next_paf_rec.query_start;

            //             log::debug!(
            //                 "Replacing {tname}:{gap_start}-{gap_stop} with {}:{qry_start}-{qry_stop}.",
            //                 &qname
            //             );
            //             new_ctgs.push(Interval::new(
            //                 qry_start.try_into()?,
            //                 qry_stop.try_into()?,
            //                 (ContigType::Query, qname.clone()),
            //             ));
            //         }
            //         // Misassembly in target and query or no information.
            //         // Cannot repair. Keep original sequence.
            //         (_, _) => {
            //             new_ctgs.push(Interval::new(
            //                 gap_start.try_into()?,
            //                 gap_stop.try_into()?,
            //                 (ContigType::Target, tname.to_owned()),
            //             ));
            //         }
            //     }
        }

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
            .entry(tname.to_owned())
            .or_insert_with(|| merge_collapse_rows_by_rle_id(collapsed_rows));
    }

    Ok(final_rows)
}
