#!/usr/bin/env python3
"""
reproduce figure 2 (tajima's d / watterson's theta vs ploidy) under the jc69 sim
model, pixy's multiallelic estimator vs the mutation-count (eta) estimator.

per ploidy cell: simulate N_OBS neutral 10kb windows at figure-2 params
(Ne=1720600, mu=5.49e-9, r=mu, sample_size=10), record per window
  - tajimaD_multi, thetaw_multi  : pixy.calc (the 'new_multi' columns)
  - tajimaD_eta,   thetaw_eta     : eta_tajima
multi values are the validation anchor: per-ploidy means should reproduce the
stored thetaw_tajimaD_summary.tsv 'new_multi' rows.

writes one row per window for aggregate_eta.py (per-ploidy means, 2.5/97.5% quantiles).
"""

from __future__ import annotations

import argparse
import glob
import os
import shutil
import subprocess
import sys
import tempfile
import warnings

import allel
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from eta_tajima import eta_stats  # noqa: E402
from pixy.calc import calc_tajima_d, calc_watterson_theta  # noqa: E402

VCFSIM_PY = os.environ.get("VCFSIM_PY", "/home/ksamuk/miniconda3/envs/vcfsim/bin/python")
VCFSIM_SRC = os.environ.get("VCFSIM_SRC", "/mnt/f/Dropbox/02_Projects/vcfsim")
TMPROOT = os.environ.get("ROUTEB_TMP", "/tmp")


def simulate(outdir, prefix, seed, reps, length, ploidy, ne, mu, recomb, samp):
    cmd = [
        VCFSIM_PY, "-m", "vcfsim",
        "--seed", str(seed), "--replicates", str(reps),
        "--sequence_length", str(length), "--ploidy", str(ploidy),
        "--Ne", str(ne), "--mu", f"{mu:.10g}",
        "--recombination_rate", f"{recomb:.10g}", "--sample_size", str(samp),
        "--percent_missing_sites", "0", "--percent_missing_genotypes", "0",
        "--chromosome", "chr1", "--output_file", os.path.join(outdir, prefix),
        "--population_mode", "1",
    ]
    cwd = VCFSIM_SRC if (VCFSIM_SRC and os.path.isdir(VCFSIM_SRC)) else None
    subprocess.run(cmd, cwd=cwd, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def read_gt(path, ploidy):
    # allel.read_vcf defaults to diploid and silently truncates polyploid genotypes
    # to the first 2 alleles; numbers= forces true ploidy
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cs = allel.read_vcf(path, fields=["calldata/GT"],
                            numbers={"calldata/GT": ploidy})
    if cs is None or "calldata/GT" not in cs:
        return None
    return allel.GenotypeArray(cs["calldata/GT"])


def _num(x):
    return float(x) if isinstance(x, (int, float)) and np.isfinite(x) else float("nan")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ploidy", type=int, required=True)
    ap.add_argument("--n_obs", type=int, default=1500)
    ap.add_argument("--length", type=int, default=10000)
    ap.add_argument("--ne", type=int, default=1720600)
    ap.add_argument("--mu", type=float, default=5.49e-9)
    ap.add_argument("--recomb", type=float, default=5.49e-9)  # r = mu
    ap.add_argument("--samp", type=int, default=10)
    ap.add_argument("--seedbase", type=int, default=1)
    ap.add_argument("--raw", required=True)
    a = ap.parse_args()

    work = tempfile.mkdtemp(prefix="eta_", dir=TMPROOT)
    rows = []
    try:
        obsdir = os.path.join(work, "obs")
        os.makedirs(obsdir)
        simulate(obsdir, "obs", a.seedbase, a.n_obs, a.length, a.ploidy,
                 a.ne, a.mu, a.recomb, a.samp)
        for of in sorted(glob.glob(os.path.join(obsdir, "*.vcf"))):
            gt = read_gt(of, a.ploidy)
            if gt is None:
                continue
            td = calc_tajima_d(gt)          # pixy multiallelic
            tw = calc_watterson_theta(gt)
            es = eta_stats(gt)              # mutation-count
            if es is None:
                continue
            rows.append((
                a.ploidy,
                _num(td.tajima_d), _num(tw.avg_theta),
                _num(es["tajimaD_eta"]), _num(es["thetaw_eta"]),
            ))
    finally:
        shutil.rmtree(work, ignore_errors=True)

    with open(a.raw, "a") as fh:
        for ploidy, dm, twm, de, twe in rows:
            fh.write(f"{ploidy}\t{dm:.6g}\t{twm:.6g}\t{de:.6g}\t{twe:.6g}\n")

    d_multi = np.array([r[1] for r in rows], float)
    d_eta = np.array([r[3] for r in rows], float)
    d_multi = d_multi[np.isfinite(d_multi)]
    d_eta = d_eta[np.isfinite(d_eta)]
    print(f"ploidy={a.ploidy} n={d_multi.size}  "
          f"D_multi mean={d_multi.mean():+.4f}  D_eta mean={d_eta.mean():+.4f}")


if __name__ == "__main__":
    main()
