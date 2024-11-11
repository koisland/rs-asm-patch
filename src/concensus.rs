use std::{collections::HashMap, str::FromStr};

use coitrees::{GenericInterval, IntervalTree};
use itertools::Itertools;
use paf::PafRecord;

use crate::{
    cigar::{get_cg_ops, CigarOp},
    interval::{get_largest_overlapping_interval, Contig, ContigType, RegionIntervalTrees},
    misassembly::{get_misassembly_from_itv, Misassembly},
};

pub fn get_concensus(
    paf_rows: &[PafRecord],
    ref_roi_itree: RegionIntervalTrees,
    ref_misasm_itree: RegionIntervalTrees,
    qry_misasm_itree: RegionIntervalTrees,
) -> eyre::Result<HashMap<String, Vec<Contig>>> {
    let mut final_rows = HashMap::new();
    for (tname, pafs) in &paf_rows
        .iter()
        .sorted_by(|a, b| a.target_name().cmp(b.target_name()))
        .chunk_by(|c| c.target_name())
    {
        let grp_roi_intervals = ref_roi_itree.0.get(tname);
        let mut new_ctgs: Vec<Contig> = vec![];
        let mut paf_recs = pafs
            .sorted_by(|a, b| a.target_start().cmp(&b.target_start()))
            .peekable();
        while let Some(paf_rec) = paf_recs.next() {
            let Some(cg) = paf_rec.cg() else {
                continue;
            };
            let cg_ops = get_cg_ops(cg)?;
            let mut target_bp_accounted = 0;
            let mut query_bp_accounted = 0;
            let mut _query_bp_deleted = 0;
            let mut _query_bp_inserted = 0;
            for (bp, op) in cg_ops.into_iter() {
                let target_start = paf_rec.target_start() + target_bp_accounted;
                let target_stop = paf_rec.target_start() + target_bp_accounted + bp;
                let query_start = paf_rec.query_start() + query_bp_accounted;
                let query_stop = paf_rec.query_start() + query_bp_accounted + bp;

                let overlap_cnt = grp_roi_intervals
                    .map(|itvs| itvs.query_count(target_start as i32, target_stop as i32))
                    .unwrap_or(0);

                // Skip alignments not within region of interest
                // Only take intervals overlapping ROIs.
                if grp_roi_intervals.is_some() && overlap_cnt == 0 {
                    target_bp_accounted += bp;
                    query_bp_accounted += bp;
                    continue;
                }

                match op {
                    CigarOp::Match => {
                        new_ctgs.push(Contig {
                            name: paf_rec.target_name().to_string(),
                            category: ContigType::Target,
                            start: target_start,
                            stop: target_stop,
                        });

                        target_bp_accounted += bp;
                        query_bp_accounted += bp;
                    }
                    // Deletion in query with respect to reference
                    CigarOp::Deletion => {
                        let largest_target_misassembly = get_largest_overlapping_interval(
                            target_start as i32,
                            target_stop as i32,
                            &ref_misasm_itree,
                            paf_rec.target_name(),
                        );
                        if largest_target_misassembly.is_some() {
                            log::debug!(
                                "Deletion in query with respect to reference ({}, {}):",
                                target_start,
                                target_stop
                            );
                            log::debug!(
                                "\tTarget name: {}, {:?}",
                                paf_rec.target_name(),
                                largest_target_misassembly.as_ref().map(|m| (
                                    m.first,
                                    m.last,
                                    m.metadata()
                                )),
                            );
                        }
                        // Added sequence in target due to misassemblies associated with drop in read coverage.
                        // Remove sequence.
                        if let Some(Ok(Misassembly::MISJOIN | Misassembly::GAP)) =
                            largest_target_misassembly
                                .and_then(|itv| itv.metadata.map(|s| Misassembly::from_str(&s)))
                        {
                            new_ctgs.push(Contig::from(ContigType::Spacer));
                        }
                        _query_bp_deleted += bp;
                        target_bp_accounted += bp;
                    }
                    CigarOp::Insertion => {
                        let largest_target_misassembly = get_largest_overlapping_interval(
                            target_start as i32,
                            target_stop as i32,
                            &ref_misasm_itree,
                            paf_rec.target_name(),
                        );
                        if largest_target_misassembly.is_some() {
                            log::debug!(
                                "Insertion in query with respect to reference ({}, {}):",
                                target_start,
                                target_stop
                            );
                            log::debug!(
                                "\tTarget name: {}, {:?}",
                                paf_rec.target_name(),
                                largest_target_misassembly.as_ref().map(|m| (
                                    m.first,
                                    m.last,
                                    m.metadata()
                                )),
                            );
                        }
                        // Deleted sequence in target as insertion in other assembly
                        // Add sequence from query.
                        if let Some(Ok(
                            Misassembly::MISJOIN | Misassembly::GAP | Misassembly::ERROR,
                        )) = largest_target_misassembly
                            .and_then(|itv| itv.metadata.map(|s| Misassembly::from_str(&s)))
                        {
                            new_ctgs.push(Contig {
                                name: paf_rec.query_name().to_owned(),
                                category: ContigType::Query,
                                start: query_start,
                                stop: query_stop,
                            });
                        }
                        _query_bp_inserted += bp;
                        query_bp_accounted += bp;
                    }
                    _ => {
                        target_bp_accounted += bp;
                        query_bp_accounted += bp;
                    }
                }
            }

            // Check gaps in alignment.
            let Some(next_paf_rec) = paf_recs.peek() else {
                continue;
            };
            let gap_start = paf_rec.target_end();
            let gap_stop = next_paf_rec.target_start();

            let ref_gap_misassembly = get_largest_overlapping_interval(
                gap_start as i32,
                gap_stop as i32,
                &ref_misasm_itree,
                paf_rec.target_name(),
            );
            let qry_gap_misassembly = get_largest_overlapping_interval(
                gap_start as i32,
                gap_stop as i32,
                &qry_misasm_itree,
                paf_rec.query_name(),
            );
            let ref_gap_mtype = get_misassembly_from_itv(ref_gap_misassembly);
            let qry_gap_mtype = get_misassembly_from_itv(qry_gap_misassembly);
            if ref_gap_mtype.is_some() || qry_gap_mtype.is_some() {
                log::debug!(
                    "Misassemblies between alignments at ({}, {}):",
                    gap_start,
                    gap_stop
                );
                log::debug!(
                    "\tTarget name: {}, {:?}",
                    paf_rec.target_name(),
                    ref_gap_mtype,
                );
                log::debug!(
                    "\tQuery name: {}, {:?}",
                    paf_rec.query_name(),
                    qry_gap_mtype,
                );
            }

            match (ref_gap_mtype, qry_gap_mtype) {
                // Fix misassembly in target with query.
                (Some(_), None) => {
                    new_ctgs.push(Contig {
                        name: paf_rec.query_name().to_owned(),
                        category: ContigType::Query,
                        start: paf_rec.query_end(),
                        stop: next_paf_rec.query_start(),
                    });
                }
                // Misassembly in target and query or no information.
                // Cannot repair. Keep original sequence.
                (Some(_), Some(_)) | (None, None) | (None, Some(_)) => {
                    new_ctgs.push(Contig {
                        name: paf_rec.target_name().to_owned(),
                        category: ContigType::Target,
                        start: gap_start,
                        stop: gap_stop,
                    });
                }
            }
        }

        let mut rle_id = 0;
        let mut collapsed_rows: Vec<(usize, Contig)> = vec![];
        // Group by ctg_name and ctg_type and find min start and max stop coordinate.
        for ((_ctg_name, ctg_categ), ctg_grp_rows) in
            &new_ctgs.iter().chunk_by(|ctg| (&ctg.name, &ctg.category))
        {
            // Store id to sort later.
            rle_id += 1;

            if *ctg_categ == ContigType::Spacer {
                continue;
            }

            let rows: Vec<&Contig> = ctg_grp_rows.collect_vec();
            let (mut min_start, mut max_stop) = (u32::MAX, 0);
            for ctg in rows.iter() {
                if ctg.start < min_start {
                    min_start = ctg.start
                }
                if ctg.stop > max_stop {
                    max_stop = ctg.stop
                }
            }
            // Create collapsed row.
            let new_ctg = Contig {
                name: rows[0].name.clone(),
                category: rows[0].category.clone(),
                start: min_start,
                stop: max_stop,
            };
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
