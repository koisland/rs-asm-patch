use std::{collections::HashMap, str::FromStr};

use coitrees::{GenericInterval, IntervalNode, IntervalTree};
use eyre::bail;
use itertools::Itertools;
use paf::PafRecord;

use crate::RegionIntervalTrees;

#[derive(Debug, Clone, PartialEq, Eq)]
enum ContigType {
    Target,
    Query,
    Spacer,
}

#[derive(Debug)]
pub struct Contig {
    name: String,
    category: ContigType,
    start: u32,
    stop: u32,
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

#[allow(clippy::upper_case_acronyms)]
#[derive(Debug, PartialEq, Eq)]
pub enum Misassembly {
    MISJOIN,
    GAP,
    COLLAPSE,
    ERROR,
    HET,
}

impl FromStr for Misassembly {
    type Err = eyre::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(match s {
            "MISJOIN" => Misassembly::MISJOIN,
            "GAP" => Misassembly::GAP,
            "COLLAPSE" | "COLLAPSE_VAR" => Misassembly::COLLAPSE,
            "ERROR" => Misassembly::ERROR,
            "HET" => Misassembly::HET,
            _ => bail!("Invalid misassembly type. ({s})"),
        })
    }
}

#[derive(Debug)]
enum CigarOp {
    Match,
    Mismatch,
    Insertion,
    Deletion,
}

impl TryFrom<char> for CigarOp {
    type Error = eyre::Error;

    fn try_from(value: char) -> Result<Self, Self::Error> {
        match value {
            '=' => Ok(CigarOp::Match),
            'X' => Ok(CigarOp::Mismatch),
            'I' => Ok(CigarOp::Insertion),
            'D' => Ok(CigarOp::Deletion),
            _ => bail!("Invalid cigar operation ({value})."),
        }
    }
}

fn get_cg_ops(cg: &str) -> eyre::Result<Vec<(u32, CigarOp)>> {
    let mut all_cg_ops = vec![];
    let mut bp_elem = None;
    let mut op_elem = None;
    for (is_digit, mut cg_ops) in &cg.chars().chunk_by(|c| c.is_ascii_digit()) {
        if is_digit {
            bp_elem = Some(cg_ops.join("").parse::<u32>()?);
        } else if let Some(cg_op) = cg_ops.next() {
            op_elem = Some(CigarOp::try_from(cg_op)?)
        } else {
            unreachable!()
        }
        if bp_elem.is_some() && op_elem.is_some() {
            all_cg_ops.push((bp_elem.unwrap(), op_elem.unwrap()));
            bp_elem = None;
            op_elem = None;
        }
    }
    Ok(all_cg_ops)
}

fn get_largest_overlapping_interval(
    start: i32,
    stop: i32,
    misasm_itree: &RegionIntervalTrees,
    name: &str,
) -> Option<coitrees::IntervalNode<Option<String>, usize>> {
    let Some(misasms) = misasm_itree.0.get(name) else {
        return None;
    };
    let mut overlapping_itvs = vec![];
    misasms.query(start, stop, |n| {
        overlapping_itvs.push(n.clone());
    });
    overlapping_itvs
        .into_iter()
        .max_by(|a, b| a.len().cmp(&b.len()))
}

fn get_misassembly_from_itv(
    itv: Option<IntervalNode<Option<String>, usize>>,
) -> Option<Misassembly> {
    itv.as_ref()
        .and_then(|m| {
            m.metadata()
                .as_ref()
                .and_then(|m| Misassembly::from_str(m).ok())
        })
        // Filter HETs that aren't serious misassemblies.
        .filter(|m| *m != Misassembly::HET)
}

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
                let query_stop = paf_rec.query_end() + query_bp_accounted + bp;

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
                        // Added sequence in target due to misassemblies associated with drop in read coverage.
                        // Remove sequence.
                        if let Some(Ok(misassembly)) = largest_target_misassembly
                            .and_then(|itv| itv.metadata.map(|s| Misassembly::from_str(&s)))
                        {
                            if let Misassembly::MISJOIN | Misassembly::GAP = misassembly {
                                new_ctgs.push(Contig::from(ContigType::Spacer));
                            }
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
                        // Deleted sequence in target as insertion in other assembly
                        // Add sequence from query.
                        if let Some(Ok(misassembly)) = largest_target_misassembly
                            .and_then(|itv| itv.metadata.map(|s| Misassembly::from_str(&s)))
                        {
                            if let Misassembly::MISJOIN | Misassembly::GAP | Misassembly::ERROR =
                                misassembly
                            {
                                new_ctgs.push(Contig {
                                    name: paf_rec.query_name().to_owned(),
                                    category: ContigType::Query,
                                    start: query_start,
                                    stop: query_stop,
                                });
                            }
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
