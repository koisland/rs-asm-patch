use std::{
    fs::File,
    io::{stdout, BufWriter, Write},
    path::PathBuf,
    str::FromStr,
};

use clap::Parser;
use coitrees::IntervalTree;
use log::LevelFilter;
use simple_logger::SimpleLogger;

mod cigar;
mod cli;
mod concensus;
mod interval;
mod io;
mod liftover;
use cli::Args;

#[allow(unused)]
fn debug_args() -> eyre::Result<Args> {
    Ok(Args {
        paf: PathBuf::from_str("data/mPanPan1_trim.paf")?,
        ref_fa: PathBuf::from_str("data/mPanPan1_merged_dedup_asm.fa.gz")?,
        query_fa: PathBuf::from_str("data/mPanPan1_merged_dedup_asm_query.fa.gz")?,
        ref_misasm_bed: PathBuf::from_str("data/mPanPan1_cen_misassemblies.bed")?,
        qry_misasm_bed: PathBuf::from_str("data/mPanPan1_cen_misassemblies_query.bed")?,
        ref_roi_bed: Some(PathBuf::from_str("data/mPanPan1.matpat.ALR.500kbp.bed")?),
        output_fa: Some(PathBuf::from_str("/dev/null")?),
        output_bed: Some(PathBuf::from_str("test.bed")?),
        log_level: String::from("Debug"),
    })
}

fn main() -> eyre::Result<()> {
    let args = Args::parse();
    // let args = debug_args()?;
    let log_level = LevelFilter::from_str(&args.log_level)?;
    SimpleLogger::new().with_level(log_level).init()?;

    log::info!("Using params: {args:#?}");

    // TODO: Add option to replace with minimap2-rs
    // TODO: Remove overlapping alignments.
    let paf_records = io::read_paf(args.paf)?;

    // Read ref roi bed.
    let ref_roi_records = io::read_bed(args.ref_roi_bed, |_| None)?;
    if let Some(ref_itvs) = ref_roi_records.as_ref() {
        log::info!(
            "Patching {} provided reference regions.",
            ref_itvs.values().map(|itvs| itvs.len()).sum::<usize>()
        );
    } else {
        log::info!("Patching all regions since no reference region bed provided.");
    }

    // Read misassemblies.
    // Always required.
    let ref_misasm_records =
        io::read_bed(Some(args.ref_misasm_bed), |rec| Some(rec.to_owned()))?.unwrap();
    let qry_misasm_records =
        io::read_bed(Some(args.qry_misasm_bed), |rec| Some(rec.to_owned()))?.unwrap();

    let mut new_ctgs = concensus::get_concensus(
        &paf_records,
        ref_roi_records,
        ref_misasm_records,
        qry_misasm_records,
    )?;
    log::info!(
        "Generated {:?} consensus sequences.",
        new_ctgs
            .iter()
            .flat_map(|(_, new_ctgs)| (new_ctgs.len() > 1).then_some(1))
            .sum::<usize>()
    );

    let mut ref_fh = io::FastaReaderHandle::new(args.ref_fa)?;
    let mut qry_fh = io::FastaReaderHandle::new(&args.query_fa)?;

    let output_fa: Box<dyn Write> =
        if let Some(outfile) = args.output_fa.filter(|fpath| *fpath != PathBuf::from("-")) {
            Box::new(BufWriter::new(File::create(outfile)?))
        } else {
            Box::new(BufWriter::new(stdout().lock()))
        };

    let output_bed = if let Some(output_bed) = args.output_bed {
        Some(BufWriter::new(File::create(output_bed)?))
    } else {
        None
    };

    io::update_contig_boundaries(&mut new_ctgs, &ref_fh, &qry_fh)?;
    log::info!("Writing concensus fasta...");
    io::write_consensus_fa(new_ctgs, &mut ref_fh, &mut qry_fh, output_fa, output_bed)?;

    Ok(())
}
