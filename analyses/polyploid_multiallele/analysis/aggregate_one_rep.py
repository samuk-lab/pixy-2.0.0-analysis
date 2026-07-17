#!/usr/bin/env python3
"""aggregate one replicate's pixy / vcftools output into a wide per-window tsv row.

the runner (02_run_arm.sh) knows which files exist for this replicate and passes
them via repeated --column-spec flags. joins on (chromosome, window_pos_1) and
emits one row per window, one column per spec, in declared order.

usage:
    aggregate_one_rep.py --replicate REP_ID
        --column-spec LABEL:FILE:KIND:VALUE_COL_CANDIDATES [--column-spec ...]
        [--header-only]

column-spec:
  LABEL   output column name (pi_new, fst_wc_new_multi, vcftools_pi, ...)
  FILE    path to the pixy/vcftools file for this replicate
  KIND    pixy | vcftools_pi | vcftools_fst
  VALUE_COL_CANDIDATES  comma-separated column names, first present wins

keys: pixy = (chromosome, window_pos_1), end window_pos_2;
      vcftools = (CHROM, BIN_START), end BIN_END.
value cols: vcftools_pi -> PI, vcftools_fst -> WEIGHTED_FST, pixy -> avg_<stat>
(avg_pi, avg_dxy, avg_wc_fst, avg_hudson_fst, avg_wattersons_theta, avg_tajima_d).

missing file -> NA for every key other files contribute (warns to stderr).
present file missing its value column -> exit 2 (schema mismatch, fail loud).

output header: replicate  chromosome  window_pos_1  window_pos_2  <label1> ...
one row per (chromosome, window_pos_1) any declared file contains.
--header-only emits just the header line.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from typing import Dict, Iterable, List, Optional, Tuple

KEY = Tuple[str, str]  # (chromosome, window_pos_1)


def _select_value_col(fields: Iterable[str], candidates: List[str], path: str) -> str:
    fields_set = set(fields)
    for c in candidates:
        if c in fields_set:
            return c
    sys.stderr.write(
        f"[aggregate_one_rep.py] ERROR: {path} has none of value cols "
        f"{candidates}; header was {sorted(fields_set)}\n"
    )
    sys.exit(2)


def _read_pixy(path: str, value_candidates: List[str]) -> Tuple[Dict[KEY, str], Dict[KEY, str]]:
    """returns (values, win2)."""
    if not os.path.exists(path):
        sys.stderr.write(f"[aggregate_one_rep.py] WARNING: missing pixy file {path}\n")
        return {}, {}
    values: Dict[KEY, str] = {}
    win2: Dict[KEY, str] = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        fields = reader.fieldnames or []
        for required in ("chromosome", "window_pos_1"):
            if required not in fields:
                sys.stderr.write(
                    f"[aggregate_one_rep.py] ERROR: {path} missing required column "
                    f"{required!r}; header was {fields}\n"
                )
                sys.exit(2)
        value_col = _select_value_col(fields, value_candidates, path)
        for row in reader:
            key = (row["chromosome"], row["window_pos_1"])
            values[key] = row[value_col]
            win2.setdefault(key, row.get("window_pos_2", "NA"))
    return values, win2


def _read_vcftools(path: str, value_candidates: List[str]) -> Tuple[Dict[KEY, str], Dict[KEY, str]]:
    if not os.path.exists(path):
        sys.stderr.write(f"[aggregate_one_rep.py] WARNING: missing vcftools file {path}\n")
        return {}, {}
    values: Dict[KEY, str] = {}
    win2: Dict[KEY, str] = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        fields = reader.fieldnames or []
        for required in ("CHROM", "BIN_START"):
            if required not in fields:
                sys.stderr.write(
                    f"[aggregate_one_rep.py] ERROR: {path} missing required column "
                    f"{required!r}; header was {fields}\n"
                )
                sys.exit(2)
        value_col = _select_value_col(fields, value_candidates, path)
        for row in reader:
            key = (row["CHROM"], row["BIN_START"])
            values[key] = row[value_col]
            win2.setdefault(key, row.get("BIN_END", "NA"))
    return values, win2


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--replicate", default="")
    p.add_argument(
        "--column-spec",
        action="append",
        default=[],
        help=(
            "Repeatable. Format: LABEL:FILE:KIND:VALUE_COL_CANDIDATES. "
            "KIND is one of pixy, vcftools_pi, vcftools_fst. "
            "VALUE_COL_CANDIDATES is a comma-separated list."
        ),
    )
    p.add_argument("--header-only", action="store_true")
    args = p.parse_args()

    # parse specs into (label, path, kind, candidates)
    specs: List[Tuple[str, str, str, List[str]]] = []
    for raw in args.column_spec:
        # split on first + last colons so a windows C:/... path in FILE survives
        first_colon = raw.find(":")
        last_colon = raw.rfind(":")
        second_to_last = raw.rfind(":", 0, last_colon)
        if first_colon == -1 or second_to_last == -1 or first_colon == second_to_last:
            sys.stderr.write(f"[aggregate_one_rep.py] ERROR: bad --column-spec {raw!r}\n")
            sys.exit(2)
        label = raw[:first_colon]
        path = raw[first_colon + 1 : second_to_last]
        kind = raw[second_to_last + 1 : last_colon]
        cands_str = raw[last_colon + 1 :]
        if kind not in ("pixy", "vcftools_pi", "vcftools_fst"):
            sys.stderr.write(
                f"[aggregate_one_rep.py] ERROR: bad kind {kind!r} in spec {raw!r}\n"
            )
            sys.exit(2)
        cands = [c.strip() for c in cands_str.split(",") if c.strip()]
        if not cands:
            sys.stderr.write(
                f"[aggregate_one_rep.py] ERROR: no value-col candidates in spec {raw!r}\n"
            )
            sys.exit(2)
        specs.append((label, path, kind, cands))

    if args.header_only:
        cols = ["replicate", "chromosome", "window_pos_1", "window_pos_2"] + [s[0] for s in specs]
        sys.stdout.write("\t".join(cols) + "\n")
        return

    if not specs:
        sys.stderr.write("[aggregate_one_rep.py] ERROR: no --column-spec given\n")
        sys.exit(2)

    # read every declared file
    column_values: List[Dict[KEY, str]] = []
    win2_master: Dict[KEY, str] = {}
    for label, path, kind, cands in specs:
        if kind == "pixy":
            vals, win2 = _read_pixy(path, cands)
        elif kind == "vcftools_pi":
            vals, win2 = _read_vcftools(path, cands)
        elif kind == "vcftools_fst":
            vals, win2 = _read_vcftools(path, cands)
        else:  # unreachable
            vals, win2 = {}, {}
        column_values.append(vals)
        for k, v in win2.items():
            win2_master.setdefault(k, v)

    all_keys: set = set()
    for d in column_values:
        all_keys.update(d.keys())

    if not all_keys:
        sys.stderr.write(
            f"[aggregate_one_rep.py] WARNING: replicate {args.replicate}: "
            f"no usable windows in any input file; emitting nothing\n"
        )
        return

    # warn on columns missing windows the union has
    for (label, path, _kind, _cands), d in zip(specs, column_values):
        missing = [k for k in all_keys if k not in d]
        if missing:
            sys.stderr.write(
                f"[aggregate_one_rep.py] WARNING rep {args.replicate}: column "
                f"{label!r} (from {path}) missing {len(missing)} window(s); "
                f"writing NA for those\n"
            )

    def sort_key(k: KEY):
        try:
            return (k[0], int(k[1]))
        except ValueError:
            return (k[0], k[1])

    for key in sorted(all_keys, key=sort_key):
        chrom, win1 = key
        win2 = win2_master.get(key, "NA")
        row = [args.replicate, chrom, win1, win2] + [d.get(key, "NA") for d in column_values]
        sys.stdout.write("\t".join(row) + "\n")


if __name__ == "__main__":
    main()
