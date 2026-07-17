#!/usr/bin/env python3
"""merge band-aid WC-FST-new columns into the diploid aggregated tables.

for each diploid 2-pop arm a separate WC-only arm (wcfst_<arm>) re-simulated the
same seeds and computed only pixy Weir-Cockerham FST. this joins fst_wc_new,
fst_wc_new_multi, fst_wc_old from data/aggregated/wcfst_<arm>.tsv into
data/aggregated/dxy_<arm>.tsv on (replicate, chromosome, window_pos_1), leaving
other columns untouched. fst_wc_old is overwritten with the band-aid value (old
pixy rerun at --fst_maf_filter 0). idempotent.

usage:  python3 analysis/merge_wcfst.py                  # all 5 diploid arms
        python3 analysis/merge_wcfst.py dxy_dip_2pop ...  # specific targets
"""
from __future__ import annotations

import csv
import os
import sys

AGG = os.path.join(os.path.dirname(__file__), "..", "data", "aggregated")
KEY = ("replicate", "chromosome", "window_pos_1")
# fst_wc_old re-merged too: band-aid reruns old pixy at --fst_maf_filter 0,
# overwriting the original maf-0.05-filtered (upward-biased) fst_wc_old
WC_COLS = ("fst_wc_new", "fst_wc_new_multi", "fst_wc_old")

ARMS = [
    "dxy_dip_2pop",
    "dxy_dip_2pop_miss10",
    "dxy_dip_2pop_miss25",
    "dxy_dip_2pop_miss50",
    "dxy_dip_2pop_miss75",
]


def read_tsv(path: str) -> tuple[list[str], list[dict[str, str]]]:
    with open(path, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        rows = list(r)
        return list(r.fieldnames or []), rows


def merge_arm(target_arm: str) -> None:
    dxy_path = os.path.join(AGG, f"{target_arm}.tsv")
    wc_arm = target_arm.replace("dxy_", "wcfst_", 1)
    wc_path = os.path.join(AGG, f"{wc_arm}.tsv")
    if not os.path.isfile(dxy_path):
        raise SystemExit(f"ERROR: missing target {dxy_path}")
    if not os.path.isfile(wc_path):
        raise SystemExit(f"ERROR: missing WC source {wc_path}")

    dxy_cols, dxy_rows = read_tsv(dxy_path)
    _, wc_rows = read_tsv(wc_path)

    wc_lookup = {tuple(row[k] for k in KEY): row for row in wc_rows}

    n_matched = 0
    n_missing = 0
    for row in dxy_rows:
        k = tuple(row[k_] for k_ in KEY)
        src = wc_lookup.get(k)
        if src is None:
            n_missing += 1
            for c in WC_COLS:
                row[c] = "NA"
        else:
            n_matched += 1
            for c in WC_COLS:
                row[c] = src.get(c, "NA")

    # existing columns, WC cols appended once
    out_cols = list(dxy_cols)
    for c in WC_COLS:
        if c not in out_cols:
            out_cols.append(c)

    tmp = dxy_path + ".tmp"
    with open(tmp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=out_cols, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(dxy_rows)
    os.replace(tmp, dxy_path)

    print(f"[merge] {target_arm}: {n_matched} windows matched, "
          f"{n_missing} unmatched (NA); cols -> {out_cols[-2:]}")


def main(argv: list[str]) -> None:
    targets = argv[1:] or ARMS
    for arm in targets:
        merge_arm(arm)


if __name__ == "__main__":
    main(sys.argv)
