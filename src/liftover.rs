use coitrees::{COITree, Interval, IntervalTree};
use eyre::bail;
use paf::PafRecord;

use crate::{
    cigar::{get_cg_ops, CigarOp},
    interval::in_roi,
};

pub struct Liftover {
    intervals: COITree<(i32, i32, CigarOp), usize>,
}

impl Liftover {
    pub fn new<T: Clone>(
        paf_rec: &PafRecord,
        filter_intervals: Option<&COITree<T, usize>>,
    ) -> eyre::Result<Self> {
        let Some(cg) = paf_rec.cg() else {
            bail!("No cigar in {paf_rec:?}.")
        };
        let cg_ops = get_cg_ops(cg)?;
        let mut target_bp_accounted = 0;
        let mut query_bp_accounted = 0;
        let mut _query_bp_deleted = 0;
        let mut _query_bp_inserted = 0;

        let mut liftover_intervals = vec![];

        for (bp, op) in cg_ops.into_iter() {
            let target_start = paf_rec.target_start() + target_bp_accounted;
            let target_stop = paf_rec.target_start() + target_bp_accounted + bp;
            let query_start = paf_rec.query_start() + query_bp_accounted;
            let query_stop = paf_rec.query_start() + query_bp_accounted + bp;

            // Skip alignments not within region of interest
            // Only take intervals overlapping ROIs.
            let is_roi = in_roi(target_start as i32, target_stop as i32, filter_intervals);
            if !is_roi {
                target_bp_accounted += bp;
                query_bp_accounted += bp;
                continue;
            }

            match op {
                CigarOp::Match | CigarOp::Mismatch => {
                    target_bp_accounted += bp;
                    query_bp_accounted += bp;
                }
                // Deletion in query with respect to reference
                CigarOp::Deletion => {
                    _query_bp_deleted += bp;
                    target_bp_accounted += bp;
                }
                CigarOp::Insertion => {
                    _query_bp_inserted += bp;
                    query_bp_accounted += bp;
                }
            }

            let interval = Interval::new(
                target_start as i32,
                target_stop as i32,
                (query_start as i32, query_stop as i32, op),
            );
            liftover_intervals.push(interval);
        }

        Ok(Liftover {
            intervals: COITree::new(&liftover_intervals),
        })
    }

    pub fn query(&self, start: i32, stop: i32) -> Option<(i32, i32, f32)> {
        let mut qry_coords = (i32::MAX, 0);
        let mut matches = 0;
        let mut mismatches = 0;
        let mut indels = 0;
        self.intervals.query(start, stop, |itv| {
            let start_diff = itv.first - start;
            let end_diff = itv.last - stop;
            let (qry_start, qry_stop, cg_op) = &itv.metadata;
            let bp = qry_stop - qry_start;
            let qry_start = *qry_start - start_diff;
            let qry_stop = *qry_stop - end_diff;
            match cg_op {
                CigarOp::Match => matches += bp,
                CigarOp::Mismatch => mismatches += bp,
                CigarOp::Insertion | CigarOp::Deletion => indels += bp,
            }
            if qry_start.is_negative() || qry_stop.is_negative() {
                return;
            }
            qry_coords.0 = std::cmp::min(qry_start, qry_coords.0);
            qry_coords.1 = std::cmp::max(qry_stop, qry_coords.1);
        });
        if qry_coords == (i32::MAX, 0) {
            return None;
        }
        // blast ident
        // https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity
        let aln_len = matches + mismatches + indels;
        let mismatches_gaps = mismatches + indels;
        let identity = (aln_len - mismatches_gaps) as f32 / aln_len as f32;
        Some((qry_coords.0, qry_coords.1, identity))
    }
}
