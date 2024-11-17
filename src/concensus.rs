use std::{collections::HashMap, fmt::Debug};

use coitrees::{Interval, IntervalTree};
use itertools::Itertools;
use paf::PafRecord;

use crate::{
    interval::{
        get_overlapping_intervals, in_roi, ContigType, RegionIntervalTrees, RegionIntervals,
    },
    liftover::Liftover,
};

pub fn get_concensus<T: Clone + Debug>(
    paf_rows: &[PafRecord],
    ref_roi_itree: Option<RegionIntervalTrees<T>>,
    ref_misasm_itree: RegionIntervalTrees<T>,
    qry_misasm_itree: RegionIntervalTrees<T>,
    bp_extend_patch: Option<u32>,
) -> eyre::Result<RegionIntervals<(ContigType, String)>> {
    let mut final_rows = HashMap::new();
    let bp_extend_patch = bp_extend_patch.unwrap_or(0) as i32;
    if bp_extend_patch != 0 {
        log::info!("Extending patched regions by {bp_extend_patch} bp.")
    }

    for (tname, pafs) in &paf_rows
        .iter()
        .sorted_by(|a, b| a.target_name().cmp(b.target_name()))
        .chunk_by(|c| c.target_name())
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
                    itvs.contains_key(rec.target_name())
                } else {
                    true
                }
            })
            .sorted_by(|a, b| a.target_start().cmp(&b.target_start()))
            .peekable();

        let mut tstart = 0;
        while let Some(paf_rec) = paf_recs.next() {
            let liftover_itree = Liftover::new(paf_rec, grp_roi_intervals)?;
            // Get misassemblies
            let Some(ref_misassemblies) = ref_misasm_itree.get(paf_rec.target_name()) else {
                continue;
            };
            let ref_rgn_misassemblies = get_overlapping_intervals(
                paf_rec.target_start() as i32,
                paf_rec.target_end() as i32,
                ref_misassemblies,
            );
            let qname = paf_rec.query_name().to_owned();

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
                let Some((qry_start, qry_stop, qry_ident)) = liftover_itree.query(mstart, mstop)
                else {
                    log::debug!(
                        "Unable to liftover {tname}:{mstart}-{mstop} ({mtype:?}) to {}.",
                        &qname
                    );
                    continue;
                };
                log::debug!("Replacing {tname}:{mstart}-{mstop} ({mtype:?}) with {}:{qry_start}-{qry_stop} with identity {qry_ident}.", &qname);

                // TODO: use identity?
                new_ctgs.push(Interval::new(
                    qry_start,
                    qry_stop,
                    (ContigType::Query, qname.clone()),
                ));
            }

            // Check gaps in alignment.
            let Some(next_paf_rec) = paf_recs
                .peek()
                .filter(|rec| rec.query_name() == paf_rec.query_name())
            else {
                // Add remainder of reference contig.
                new_ctgs.push(Interval::new(
                    tstart,
                    paf_rec.target_len().try_into()?,
                    (ContigType::Target, paf_rec.target_name().to_owned()),
                ));
                break;
            };

            // Case 1: Same qry ctg.
            // qry:     --------
            // alns:    ---  ---
            // ref:     ---??---

            // Case 2: Different qry ctgs
            // qry:     ---
            // qry:          ---
            // alns:    ---  ---
            // ref:     ---??---
            let gap_start = paf_rec.target_end();
            let gap_stop = next_paf_rec.target_start();

            let overlap_cnt = grp_roi_intervals
                .map(|itvs| itvs.query_count(gap_start as i32, gap_stop as i32))
                .unwrap_or(0);

            if overlap_cnt == 0 {
                continue;
            }
            let ref_aln_gap_misassemblies =
                get_overlapping_intervals(gap_start as i32, gap_stop as i32, ref_misassemblies);
            let qry_aln_gap_miassemblies = qry_misasm_itree
                .get(paf_rec.query_name())
                .map_or_else(Vec::new, |itvs| {
                    get_overlapping_intervals(gap_start as i32, gap_stop as i32, itvs)
                });

            let ref_aln_gap_misassembled = !ref_aln_gap_misassemblies.is_empty();
            let qry_aln_gap_misassembled = !qry_aln_gap_miassemblies.is_empty();
            if ref_aln_gap_misassembled || qry_aln_gap_misassembled {
                log::debug!(
                    "Misassemblies between {tname}:{gap_start}-{gap_stop} with no alignment with {}",
                    &qname
                );
                log::debug!("\tTarget: {:?}", ref_aln_gap_misassemblies,);
                log::debug!("\tQuery: {:?}", qry_aln_gap_miassemblies,);
            }

            match (ref_aln_gap_misassembled, qry_aln_gap_misassembled) {
                // Fix misassembly in target with query.
                (true, false) => {
                    let qry_start = paf_rec.query_end();
                    let qry_stop = next_paf_rec.query_start();

                    log::debug!(
                        "Replacing {tname}:{gap_start}-{gap_stop} with {}:{qry_start}-{qry_stop}.",
                        &qname
                    );
                    new_ctgs.push(Interval::new(
                        qry_start.try_into()?,
                        qry_stop.try_into()?,
                        (ContigType::Query, qname.clone()),
                    ));
                }
                // Misassembly in target and query or no information.
                // Cannot repair. Keep original sequence.
                (_, _) => {
                    new_ctgs.push(Interval::new(
                        gap_start.try_into()?,
                        gap_stop.try_into()?,
                        (ContigType::Target, tname.to_owned()),
                    ));
                }
            }
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

        final_rows.entry(tname.to_owned()).or_insert(
            collapsed_rows
                .into_iter()
                .sorted_by(|a, b| a.0.cmp(&b.0))
                .map(|c| c.1)
                .collect_vec(),
        );
    }

    Ok(final_rows)
}
