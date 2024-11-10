use std::str::FromStr;

use coitrees::IntervalTree;
use eyre::bail;
use itertools::Itertools;
use paf::PafRecord;

use crate::RegionIntervalTrees;

#[derive(Debug)]
enum ContigType {
    Target,
    Query,
    Spacer,
}

#[derive(Debug)]
struct Contig {
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

#[derive(Debug)]
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
        } else if let Some(cg_op) = cg_ops.nth(0) {
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
    aln: &PafRecord,
) -> Option<coitrees::IntervalNode<Option<String>, usize>> {
    let Some(misasms) = misasm_itree.get(aln.target_name()) else {
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

pub fn get_concensus(
    paf_rows: &[PafRecord],
    roi_itree: RegionIntervalTrees,
    misasm_itree: RegionIntervalTrees,
) -> eyre::Result<()> {
    for (tname, pafs) in &paf_rows
        .iter()
        .sorted_by(|a, b| a.target_name().cmp(b.target_name()))
        .chunk_by(|c| c.target_name())
    {
        let grp_roi_intervals = roi_itree.get(tname);
        let mut new_ctgs: Vec<Contig> = vec![];
        let mut paf_recs = pafs.peekable();
        while let Some(paf_rec) = paf_recs.next() {
            let Some(cg) = paf_rec.cg() else {
                continue;
            };
            let cg_ops = get_cg_ops(&cg)?;
            let mut target_bp_accounted = 0;
            let mut query_bp_accounted = 0;
            let mut query_bp_deleted = 0;
            let mut query_bp_inserted = 0;
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
                if overlap_cnt == 0 {
                    target_bp_accounted += bp;
                    query_bp_accounted += bp;
                    continue;
                }

                match op {
                    CigarOp::Match => {
                        new_ctgs.push(
                            Contig {
                                name: paf_rec.target_name().to_string(),
                                category: ContigType::Target,
                                start: target_start,
                                stop: target_stop,
                            }
                        );

                        target_bp_accounted += bp;
                        query_bp_accounted += bp;
                    }
                    // Deletion in query with respect to reference
                    CigarOp::Deletion => {
                        let largest_target_misassembly = get_largest_overlapping_interval(
                            target_start as i32,
                            target_stop as i32,
                            &misasm_itree,
                            paf_rec,
                        );
                        // Added sequence in target due to misassemblies associated with drop in read coverage.
                        // Remove sequence.
                        if let Some(Ok(misassembly)) = largest_target_misassembly
                            .and_then(|itv| itv.metadata.map(|s| Misassembly::from_str(&s)))
                        {
                            if let Misassembly::MISJOIN | Misassembly::GAP = misassembly {
                                new_ctgs.push(
                                    Contig::from(ContigType::Spacer)
                                );
                            }
                        }
                        query_bp_deleted += bp;
                        target_bp_accounted += bp;
                    }
                    CigarOp::Insertion => {
                        let largest_target_misassembly = get_largest_overlapping_interval(
                            target_start as i32,
                            target_stop as i32,
                            &misasm_itree,
                            paf_rec,
                        );
                        // Deleted sequence in target as insertion in other assembly
                        // Add sequence from query.
                        if let Some(Ok(misassembly)) = largest_target_misassembly
                            .and_then(|itv| itv.metadata.map(|s| Misassembly::from_str(&s)))
                        {
                            if let Misassembly::MISJOIN | Misassembly::GAP | Misassembly::ERROR =
                                misassembly
                            {
                                new_ctgs.push(
                                    Contig {
                                        name: paf_rec.query_name().to_owned(),
                                        category: ContigType::Query,
                                        start: query_start,
                                        stop: query_stop,
                                    }
                                );
                            }
                        }
                        query_bp_inserted += bp;
                        query_bp_accounted += bp;
                    }
                    _ => {
                        target_bp_accounted += bp;
                        query_bp_accounted += bp;
                    },
                }
            }

            // Check gaps in alignment.
            let Some(next_paf_rec) = paf_recs.peek() else {
                continue;
            };
            let aln_gap_start = paf_rec.target_end();
            let aln_gap_stop = next_paf_rec.target_start();

        }
    }

    Ok(())
}
