import pprint
from typing import Iterable
import polars as pl
from intervaltree import Interval
from itertools import combinations
tests = [
    (
        Interval(8365339, 24983596, "chr17_pat_hsa18"),
        [
            (0, 409620, "h1tg000053l"),
            (73334462, 77196229, "h2tg000001l"),
            (5229808, 9992334, "h2tg000040l"),
            (4155690, 5229705, "h2tg000040l"),
            (0, 4155563, "h2tg000040l"),
            (8135, 1306939, "h2tg000041l"),
            (0, 489040, "h2tg000042l"),
        ]
    ),
    (
        Interval(73876415, 74812999, "chr5_pat_hsa6"),
        [
            (71670703, 71691116, "h1tg000017l"),
            (71653942, 71674651, "h1tg000017l"),
            (71620033, 71641130, "h1tg000017l"),
            (71603268, 71624365, "h1tg000017l"),
            (71634779, 71655500, "h1tg000017l"),
            (116426891, 116828997, "h1tg000022l"),
            (116337563, 116428864, "h1tg000022l"),
            (0, 62934, "h1tg000091l"),
            (0, 132609, "h2tg000058l"),
        ]
    ),
    (
        Interval(16995741, 19569078, "chr15_mat_hsa14"),
        [
            (89314215, 91081164, "h1tg000013l"),
            (89261880, 89312209, "h1tg000013l"),
            (89245435, 89259589, "h1tg000013l"),
            (88882357, 89243620, "h1tg000013l"),
            (58, 527, "h1tg000092l"),
            (808, 17030, "h1tg000092l"),
            (361208, 396846, "h2tg000049l"),
            (0, 354655, "h2tg000049l"),
        ]
    ),
    (
        Interval(57475940, 57519664, "chr10_pat_hsa12"),
        [(81525222, 81571660, "h2tg000006l")]
    )
]


# https://stackoverflow.com/a/57526728
def dp_itv_subset(itvs: Iterable[Interval], target_len: int):
    # Generate dp matrix with lists.
    vals = [None] * (target_len + 1)
    vals[0] = ()

    # Keep track of best.
    best_v = None
    best_k = 0

    # Iterate thru all intervals.
    for itv in itvs:
        itv_len = itv.length()
        # Iterate in reverse thru all possible target lengths.
        for i in range(target_len - 1, -1, -1):
            new_val = i + itv_len
            if vals[i] is not None:
                # print(best_k, best_v)
                # print(i, new_val)
                # print(vals[i])
                if new_val <= target_len and vals[new_val] is None:
                    vals[new_val] = (*vals[i], itv_len)
                # First, check if overshoot value.
                # And if no best value and new_val closer.
                if new_val > target_len and (best_v is None or new_val < best_v):
                    best_v = new_val
                    best_k = (*vals[i], itv_len)

    if vals[target_len] is not None:
        return vals[target_len]
    else:
        return best_k

for target_itv, all_itvs in tests:
    df = pl.DataFrame(all_itvs, orient="row", schema=["st", "end", "name"])
    itvs = [
        Interval(itv["st"], itv["end"], itv["name"])
        for itv in df.group_by(["name"]).agg(pl.col("st").min(), pl.col("end").max()).iter_rows(named=True)
    ]
    res = dp_itv_subset(
        itvs, target_itv.length(),
    )
    print(target_itv, target_itv.length())
    pprint.pprint(itvs)
    print(res)
