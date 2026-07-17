#!/usr/bin/env python3
"""
figure-2 comparison table from the per-window eta dumps.

reads results/eta_raw/cell_*.raw.tsv (ploidy  D_multi  thetaw_multi  D_eta  thetaw_eta),
prints per-ploidy mean and 2.5/97.5% quantiles of tajima's d and watterson's theta
for the multiallelic vs mutation-count (eta) estimator. stored figure-2 'new_multi' /
'new' means (0% missingness) and theory targets shown for validation.
"""

from __future__ import annotations

import glob
import sys
from collections import defaultdict

import numpy as np

# stored figure-2 references at 0% missingness (analysis/thetaw_tajimaD_summary.tsv)
REF = {
    2: dict(D_multi=0.0720, D_bi=-0.0332, tw_multi=0.035383, tw_theory=0.034906),
    4: dict(D_multi=0.1344, D_bi=+0.0100, tw_multi=0.034717, tw_theory=0.034528),
    6: dict(D_multi=0.1602, D_bi=+0.0243, tw_multi=0.034402, tw_theory=0.034303),
    8: dict(D_multi=0.1748, D_bi=+0.0309, tw_multi=0.034203, tw_theory=0.034142),
}


def q(a, p):
    return float(np.quantile(a, p)) if a.size else float("nan")


def main() -> None:
    pat = sys.argv[1] if len(sys.argv) > 1 else "results/eta_raw/cell_*.raw.tsv"
    data = defaultdict(lambda: defaultdict(list))
    for f in glob.glob(pat):
        with open(f) as fh:
            for ln in fh:
                p = ln.rstrip("\n").split("\t")
                if len(p) < 5:
                    continue
                ploidy = int(p[0])
                for key, val in zip(("Dm", "twm", "De", "twe"), p[1:5]):
                    try:
                        v = float(val)
                    except ValueError:
                        continue
                    if np.isfinite(v):
                        data[ploidy][key].append(v)

    print("=== Tajima's D ===")
    print(f"{'ploidy':>6} {'n':>5} | {'D_multi(re-sim)':>22} {'[ref]':>8} | "
          f"{'D_eta(corrected)':>22} | {'D_biallelic[ref]':>16}")
    for ploidy in sorted(data):
        dm = np.array(data[ploidy]["Dm"], float)
        de = np.array(data[ploidy]["De"], float)
        r = REF.get(ploidy, {})
        print(f"{ploidy:>6} {dm.size:>5} | "
              f"{dm.mean():+.4f} [{q(dm,0.025):+.3f},{q(dm,0.975):+.3f}] "
              f"{r.get('D_multi', float('nan')):>+8.4f} | "
              f"{de.mean():+.4f} [{q(de,0.025):+.3f},{q(de,0.975):+.3f}] | "
              f"{r.get('D_bi', float('nan')):>+16.4f}")

    print("\n=== Watterson's theta  (theory target in brackets) ===")
    print(f"{'ploidy':>6} {'n':>5} | {'thetaw_multi(re-sim)':>22} {'[ref]':>9} | "
          f"{'thetaw_eta(corrected)':>22} {'[theory]':>9}")
    for ploidy in sorted(data):
        twm = np.array(data[ploidy]["twm"], float)
        twe = np.array(data[ploidy]["twe"], float)
        r = REF.get(ploidy, {})
        print(f"{ploidy:>6} {twm.size:>5} | "
              f"{twm.mean():.5f} [{q(twm,0.025):.4f},{q(twm,0.975):.4f}] "
              f"{r.get('tw_multi', float('nan')):>9.5f} | "
              f"{twe.mean():.5f} [{q(twe,0.025):.4f},{q(twe,0.975):.4f}] "
              f"{r.get('tw_theory', float('nan')):>9.5f}")


if __name__ == "__main__":
    main()
