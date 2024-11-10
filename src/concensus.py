import argparse
import sys
from itertools import islice
from collections import defaultdict, deque
from typing import NamedTuple, Self, TextIO
from enum import StrEnum, auto

import cigar
import polars as pl
import pysam
from intervaltree import Interval, IntervalTree

DEF_PAF_COLS = [
    ("qname", pl.String),
    ("qlen", pl.Int64),
    ("qst", pl.Int64),
    ("qend", pl.Int64),
    ("ort", pl.String),
    ("tname", pl.String),
    ("tlen", pl.Int64),
    ("tst", pl.Int64),
    ("tend", pl.Int64),
    ("matches", pl.Int64),
    ("aln_blk_len", pl.Int64),
    ("q", pl.Int64),
]


class ContigType(StrEnum):
    Target = auto()
    Query = auto()
    Spacer = auto()


class Contig(NamedTuple):
    name: str
    category: ContigType
    start: int = 0
    end: int = 0

    @classmethod
    def spacer(cls) -> Self:
        return Contig(name="", category=ContigType.Spacer)


def read_roi_bed(roi_bed: str | TextIO | None) -> defaultdict[str, IntervalTree]:
    roi_intervals = defaultdict(IntervalTree)
    if roi_bed:
        df = pl.read_csv(
            roi_bed,
            has_header=False,
            separator="\t",
            new_columns=["name", "start", "stop"],
        )
        for row in df.iter_rows():
            roi_intervals[row[0]].add(Interval(row[1], row[2]))

    return roi_intervals

def read_paf(paf: str | TextIO) -> pl.DataFrame:
    paf_rows = []
    desired_tags = {"cg"}
    for line in paf:
        line: bytes
        paf_rec = {}
        for i, elem in enumerate(line.decode().strip().split(), 0):
            try:
                val, _ = DEF_PAF_COLS[i]
                paf_rec[val] = elem
            except IndexError:
                tag, _, tag_val = elem.split(":")
                if tag not in desired_tags:
                    continue
                paf_rec[tag] = tag_val

        paf_rows.append(paf_rec)

    return (
        pl.DataFrame(
            paf_rows, schema={"cg": pl.String} | dict(DEF_PAF_COLS)
        )
        .filter((pl.col("q") == 60) & (pl.col("ort") == "+"))
        .sort(by=["tst"])
    )

def create_intervaltrees(
    file: str | None,
) -> tuple[defaultdict[str, IntervalTree], defaultdict[str, IntervalTree]]:
    contig_intervals = defaultdict(IntervalTree)
    misasm_intervals = defaultdict(IntervalTree)
    if not file:
        return contig_intervals, misasm_intervals

    df_query_misasm = pl.read_csv(
        file,
        has_header=False,
        separator="\t",
        new_columns=["name", "start", "stop", "misassembly"],
    )
    for row in df_query_misasm.iter_rows():
        name: str = row[0]
        misassembly_type = row[3]
        contig_intervals[name].add(Interval(row[1], row[2], misassembly_type))
        misasm_intervals[misassembly_type].add(Interval(row[1], row[2], name))

    return contig_intervals, misasm_intervals

# https://docs.python.org/3/library/itertools.html
def sliding_window(iterable, n):
    iterator = iter(iterable)
    window = deque(islice(iterator, n - 1), maxlen=n)
    if len(window) < n:
        for _ in range(n - len(window)):
            window.append(None)
        yield tuple(window)
        return
    for x in iterator:
        window.append(x)
        yield tuple(window)
    

def main():
    ap = argparse.ArgumentParser(
        description="Generate a concensus sequence from an alignment between two contigs."
    )
    ap.add_argument(
        "-i",
        "--input_paf",
        type=argparse.FileType("rb"),
        required=True,
        help="Alignment PAF.",
    )
    ap.add_argument(
        "-r",
        "--input_ref_fa",
        type=str,
        required=True,
        help="Reference fasta sequence.",
    )
    ap.add_argument(
        "-q", "--input_query_fa", type=str, required=True, help="Query fasta sequence."
    )
    ap.add_argument(
        "-o",
        "--output_fa",
        default=sys.stdout,
        type=argparse.FileType("wt"),
        help="Output fasta. Can be merged with the --merge flag.",
    )
    ap.add_argument(
        "--input_ref_misasm_bed",
        type=argparse.FileType("rb"),
        default=None,
        help="Reference misassembly bed.",
    )
    ap.add_argument(
        "--input_query_misasm_bed",
        type=argparse.FileType("rb"),
        default=None,
        help="Query misassembly bed.",
    )
    ap.add_argument(
        "-b",
        "--output_bed",
        default=None,
        type=argparse.FileType("wt"),
        help="Output bedfile.",
    )
    ap.add_argument(
        "--roi_bed",
        type=argparse.FileType("rb"),
        default=None,
        help="Regions to filter alignments in reference.",
    )
    ap.add_argument("--merge", action="store_true", help="Merge output_fa sequences.")
    ap.add_argument("--verbose", action="store_true", help="Verbose logging output.")

    args = ap.parse_args()

    roi_intervals = read_roi_bed(args.roi_bed)
    df_paf = read_paf(args.input_paf)

    target_ctg_misassemblies, target_misassemblies = create_intervaltrees(
        args.input_ref_misasm_bed
    )
    query_ctg_misassemblies, query_misassemblies = create_intervaltrees(
        args.input_query_misasm_bed
    )
    # TODO: chr10_mat and h1tg000019l
    for grp, df_grp in df_paf.group_by(["tname"], maintain_order=True):
        grp = grp[0]
        alns = df_grp.rows(named=True)
        grp_roi_intervals = roi_intervals[grp]
        new_ctgs: list[Contig] = []
        for aln, next_aln in sliding_window(alns, 2):
            cg = cigar.Cigar(aln["cg"])
            target_bp_accounted = 0
            query_bp_accounted = 0
            query_bp_deleted = 0
            query_bp_inserted = 0

            for bp, op in cg.items():
                target_interval = Interval(
                    aln["tst"] + target_bp_accounted, aln["tst"] + target_bp_accounted + bp
                )
                query_interval = Interval(
                    aln["qst"] + query_bp_accounted, aln["qst"] + query_bp_accounted + bp
                )

                # Skip alignments not within region of interest
                # Only take intervals overlapping ROIs.
                if not grp_roi_intervals.overlaps(target_interval):
                    # print(grp, grp_roi_intervals, target_interval)
                    target_bp_accounted += bp
                    query_bp_accounted += bp
                    continue

                new_ctg = None
                if op == "=":
                    new_ctg = Contig(
                        aln["tname"],
                        ContigType.Target,
                        aln["tst"] + target_bp_accounted,
                        aln["tst"] + target_bp_accounted + bp,
                    )
                    target_bp_accounted += bp
                    query_bp_accounted += bp
                # Deletion in query with respect to reference.
                elif op == "D":
                    overlapping_target_misasms = target_ctg_misassemblies[
                        aln["tname"]
                    ].overlap(target_interval)
                    overlapping_query_misasms = query_ctg_misassemblies[
                        aln["qname"]
                    ].overlap(query_interval)

                    if args.verbose:
                        print(
                            f"# Deletion in query ({query_interval}) relative to target ({target_interval}):",
                            file=sys.stderr,
                        )
                        print(
                            "Target name:",
                            aln["tname"],
                            overlapping_target_misasms,
                            file=sys.stderr,
                        )
                        print(
                            "Query name:",
                            aln["qname"],
                            overlapping_query_misasms,
                            file=sys.stderr,
                        )

                    if overlapping_target_misasms and overlapping_query_misasms:
                        largest_target_misassembly: Interval = max(
                            overlapping_target_misasms, key=lambda x: x.begin - x.end
                        )
                        largest_query_misassembly: Interval = max(
                            overlapping_query_misasms, key=lambda x: x.begin - x.end
                        )
                        match (
                            largest_target_misassembly.data,
                            largest_query_misassembly.data,
                        ):
                            # IF was just gap, could conclude that target had sequence added.
                            # BUT because MISJOIN, deletion in query explained.
                            # TWO OPTIONS: Add spacer deleting non-mapping seq or leave alone.
                            case ("GAP", "MISJOIN"):
                                new_ctg = Contig.spacer()
                            case ("COLLAPSE", "MISJOIN"):
                                pass
                            case _:
                                pass
                    elif overlapping_target_misasms:
                        largest_misassembly: Interval = max(
                            overlapping_target_misasms, key=lambda x: x.length()
                        )
                        match largest_misassembly.data:
                            # If misjoin or gap in query, sequence in target is added. Remove it from target.
                            case "MISJOIN" | "GAP":
                                new_ctg = Contig.spacer()
                            case "COLLAPSE":
                                pass
                            case "COLLAPSE_VAR":
                                pass
                            case _:
                                pass
                        pass
                    elif overlapping_query_misasms:
                        largest_misassembly: Interval = max(
                            overlapping_query_misasms, key=lambda x: x.length()
                        )

                        pass

                    query_bp_deleted += bp
                    target_bp_accounted += bp
                # Insertion in query with respect to reference.
                elif op == "I":
                    overlapping_target_misasms = target_ctg_misassemblies[
                        aln["tname"]
                    ].overlap(target_interval)
                    overlapping_query_misasms = query_ctg_misassemblies[
                        aln["qname"]
                    ].overlap(query_interval)

                    if args.verbose:
                        print(
                            f"# Insertion in query ({query_interval}) relative to target ({target_interval}):",
                            file=sys.stderr,
                        )
                        print(
                            "Target name:",
                            aln["tname"],
                            overlapping_target_misasms,
                            file=sys.stderr,
                        )
                        print(
                            "Query name:",
                            aln["qname"],
                            overlapping_query_misasms,
                            file=sys.stderr,
                        )
                    elif overlapping_target_misasms:
                        largest_misassembly: Interval = max(
                            overlapping_target_misasms, key=lambda x: x.length()
                        )
                        match largest_misassembly.data:
                            # If misjoin or gap, sequence in target is missing. Add it from query.
                            case "MISJOIN" | "GAP" | "ERROR":
                                new_ctg = Contig(
                                    aln["qname"],
                                    ContigType.Query,
                                    aln["qst"] + query_bp_accounted,
                                    aln["qst"] + query_bp_accounted + bp,
                                )
                            case "COLLAPSE":
                                pass
                            case "COLLAPSE_VAR":
                                pass
                            case _:
                                pass

                    elif overlapping_query_misasms:
                        largest_misassembly = max(
                            overlapping_query_misasms, key=lambda x: x.begin - x.end
                        )
                        match largest_misassembly.data:
                            case "MISJOIN" | "GAP":
                                pass
                            case "COLLAPSE":
                                pass
                            case "COLLAPSE_VAR":
                                pass
                            case _:
                                pass

                    query_bp_inserted += bp
                    query_bp_accounted += bp
                else:
                    target_bp_accounted += bp
                    query_bp_accounted += bp

                # If new ctg
                if new_ctg:
                    new_ctgs.append(new_ctg)

            if not next_aln:
                continue

            gap_interval = Interval(aln["tend"], next_aln["tst"])
            # TODO: Does this need to be specific for either.
            overlapping_target_misasms: set[Interval] = target_ctg_misassemblies[
                aln["tname"]
            ].overlap(gap_interval)
            overlapping_query_misasms: set[Interval] = query_ctg_misassemblies[
                aln["qname"]
            ].overlap(gap_interval)
            if args.verbose:
                print("# Misassembly between alignments:", gap_interval, file=sys.stderr)
                print(
                    "Target name:",
                    aln["tname"],
                    overlapping_target_misasms,
                    file=sys.stderr,
                )
                print(
                    "Query name:",
                    aln["qname"],
                    overlapping_query_misasms,
                    file=sys.stderr,
                )
            try:
                misasm_in_query = next(
                    ms for ms in overlapping_query_misasms if ms.data
                )
            except StopIteration:
                misasm_in_query = None

            try:
                misasm_in_target = next(
                    ms for ms in overlapping_target_misasms
                )
            except StopIteration:
                misasm_in_target = None

            original_ctg = Contig(
                aln["tname"],
                ContigType.Target,
                gap_interval.begin,
                gap_interval.end,
            )

            # misasm in target and query. Cannot repair. Keep original sequence.
            if misasm_in_query and misasm_in_target:
                replacement_interval = original_ctg
            # Fix misasm in target with query
            elif not misasm_in_query and misasm_in_target:
                replacement_start = aln["qend"]
                replacement_end = next_aln["qst"]

                replacement_interval = Contig(
                    aln["qname"],
                    ContigType.Query,
                    # Adjust query coordinates by what was inserted/deleted.
                    replacement_start,
                    replacement_end,
                )
            else:
                # Unknown case.
                replacement_interval = original_ctg

            new_ctgs.append(replacement_interval)

        df_new_ctgs = pl.DataFrame(
            [tuple(ctg) for ctg in new_ctgs],
            orient="row",
            schema={
                "ctg_name": pl.String,
                "ctg_type": pl.String,
                "st": pl.Int64,
                "end": pl.Int64,
            },
        )

        # with pl.Config(fmt_str_lengths=1000, tbl_width_chars=1000):
        #     breakpoint()
        
        df_new_ctgs = (
            df_new_ctgs
            # Group and aggregate to min max coords
            .with_columns(row_grp_id=(pl.col("ctg_name") + pl.col("ctg_type")).rle_id())
            .group_by(["row_grp_id"])
            .agg(
                pl.col("ctg_name").first(),
                pl.col("ctg_type").first(),
                pl.col("st").min(),
                pl.col("end").max(),
            )
            .sort(by="row_grp_id")
            # Remove group agg spacers.
            .filter(pl.col("ctg_type") != "spacer")
            .with_row_index()
            # Then add boundary coordinates.
            .with_columns(
                st=pl.when(pl.col("index") == pl.col("index").first())
                .then(pl.lit(0))
                .otherwise(pl.col("st")),
                end=pl.when(pl.col("index") == pl.col("index").last())
                .then(None)
                .otherwise(pl.col("end")),
            )
            .drop("index", "row_grp_id")
        )

        # with pl.Config(fmt_str_lengths=1000, tbl_width_chars=1000):
        #     breakpoint()

        if df_new_ctgs.is_empty():
            return

        main_ctg = df_new_ctgs["ctg_name"][0]
        if args.merge:
            args.output_fa.write(f">{main_ctg}\n")

        with (
            pysam.FastaFile(args.input_ref_fa) as ref_fa,
            pysam.FastaFile(args.input_query_fa) as query_fa,
        ):
            for it in df_new_ctgs.iter_rows(named=True):
                st, end, ctg_name, ctg_type = it["st"], it["end"], it["ctg_name"], it["ctg_type"]

                if ctg_type == "query":
                    fa_file = query_fa
                elif ctg_type == "target":
                    fa_file = ref_fa
                else:
                    fa_file = ref_fa

                # TODO: Some cases where contig cut short.
                end = (
                    end
                    if end
                    else fa_file.get_reference_length(ctg_name)
                )
                if st > end or st == end:
                    continue

                try:
                    seq = fa_file.fetch(ctg_name, start=st, end=end)
                except ValueError:
                    print(f"Failed to extract: {ctg_name}:{st}-{end}", file=sys.stderr)
                    continue
                header = f'>{ctg_name}:{st}-{end}\n'
                if not args.merge:
                    args.output_fa.write(header)
                    args.output_fa.write(seq)
                    args.output_fa.write("\n")
                else:
                    args.output_fa.write(seq)

                if args.output_bed:
                    args.output_bed.write(
                        f"{ctg_name}\t{st}\t{end}\t{ctg_type}\t{main_ctg}\n"
                    )

            if args.merge:
                args.output_fa.write("\n")


if __name__ == "__main__":
    raise SystemExit(main())
