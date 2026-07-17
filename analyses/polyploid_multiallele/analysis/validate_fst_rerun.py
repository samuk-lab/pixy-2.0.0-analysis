#!/usr/bin/env python3
"""validate the Figure-4 FST rerun at n = 2000 (2026-07-15).

produces the numbers in finite_sites_fst_derivation.md section 5.4 and
HANDOFF_fst_multiallelic_2026-07-14.md. run from this directory:
    python validate_fst_rerun.py

reads ../data/aggregated/dxy_*_2pop_theta*.tsv (rerun, pixy 02d6d91) and
../data/aggregated/_pre_fstfix_2026-07-14/ (pre-fix backup).

three checks:
  1. new_multi_fstbi reproduces the pre-fix new_multi column exactly -> the old
     mislabelled column was the "read multiallelic in, filter FST to biallelic"
     estimand, not garbage.
  2. pooled FST vs each estimator's own expectation, bootstrap SE over reps. the
     ratio test from the derivation's section 5.4 is NOT reproduced: it assumes
     bi == E_inf exactly (false), so fails at n = 2000 for unrelated reasons.
  3. E[FST]_finite = 1 - E[pi_w]/E[dxy] against the data.

pooling is ratio-of-sums (Bhatia 2013): FST = sum(num)/sum(den); averaging
per-window or per-rep ratios is Jensen-biased downward.
"""

import csv
import math
import os
import random
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
NEW_DIR = os.path.join(HERE, "..", "data", "aggregated")
OLD_DIR = os.path.join(NEW_DIR, "_pre_fstfix_2026-07-14")

# model constants -- must match config/sim_params.tsv for the dxy_*_theta arms
NE = 1720600.0
SPLIT_TIME = 1e6
ALPHA = 2.0
TAU = SPLIT_TIME / (2 * NE)

MU = {"005": 7.267e-10, "010": 1.4534e-9, "025": 3.6335e-9,
      "050": 7.267e-9, "100": 1.4534e-8}

PLOIDIES = ["dip", "tet", "hex", "oct"]
THETAS = ["005", "010", "025", "050", "100"]

NBOOT = 400
SEED = 12345


def theta_of(tag):
    return 4 * NE * MU[tag]


def fst_infinite():
    """Slatkin (1991): ratio of mean coalescence times, theta-free."""
    return (TAU + (ALPHA - 1) * (1 - math.exp(-TAU))) / (TAU + ALPHA)


def L_dxy(theta):
    return math.exp(-(4.0 / 3.0) * theta * TAU) / (1 + ALPHA * (4.0 / 3.0) * theta)


def L_pi_w(theta):
    k = 1 + (4.0 / 3.0) * theta
    return ((1 - math.exp(-TAU * k)) / k
            + math.exp(-TAU * k) / (1 + ALPHA * (4.0 / 3.0) * theta))


def fst_finite(theta):
    """1 - E[pi_w]/E[dxy] under Jukes-Cantor; see derivation section 3."""
    return 1 - (1 - L_pi_w(theta)) / (1 - L_dxy(theta))


def e_pi_w(theta):
    return 0.75 * (1 - L_pi_w(theta))


def e_dxy(theta):
    return 0.75 * (1 - L_dxy(theta))


def fnum(s):
    if s in ("NA", "", "nan", None):
        return None
    try:
        return float(s)
    except ValueError:
        return None


def load(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def arm_path(base, p, t):
    return os.path.join(base, "dxy_%s_2pop_theta%s.tsv" % (p, t))


def arms_present(base):
    out = []
    for p in PLOIDIES:
        for t in THETAS:
            path = arm_path(base, p, t)
            if os.path.exists(path):
                out.append((p, t, path))
    return out


def pool(rows, suffix):
    """ratio-of-sums FST over all rows for one variant suffix."""
    n = d = 0.0
    for r in rows:
        a = fnum(r.get("fst_hudson_num_" + suffix))
        b = fnum(r.get("fst_hudson_den_" + suffix))
        if a is None or b is None:
            continue
        n += a
        d += b
    return (n / d) if d else float("nan")


def mean_col(rows, col):
    vals = [fnum(r.get(col)) for r in rows]
    vals = [v for v in vals if v is not None]
    return (sum(vals) / len(vals)) if vals else float("nan")


def boot_pooled(rows, suffixes, nboot=NBOOT, seed=SEED):
    """bootstrap over replicates -- the independent simulation unit.

    windows within a replicate share a genealogy, so the replicate (not the
    window) is what gets resampled.
    """
    byrep = defaultdict(list)
    for r in rows:
        byrep[r["replicate"]].append(r)
    reps = sorted(byrep)
    sums = {}
    for s in suffixes:
        sums[s] = {}
        for rep in reps:
            n = d = 0.0
            for r in byrep[rep]:
                a = fnum(r.get("fst_hudson_num_" + s))
                b = fnum(r.get("fst_hudson_den_" + s))
                if a is None or b is None:
                    continue
                n += a
                d += b
            sums[s][rep] = (n, d)
    rng = random.Random(seed)
    draws = defaultdict(list)
    nrep = len(reps)
    for _ in range(nboot):
        pick = [reps[rng.randrange(nrep)] for _ in range(nrep)]
        for s in suffixes:
            n = sum(sums[s][p][0] for p in pick)
            d = sum(sums[s][p][1] for p in pick)
            draws[s].append(n / d if d else float("nan"))
    out = {}
    for k, v in draws.items():
        m = sum(v) / len(v)
        var = sum((x - m) ** 2 for x in v) / (len(v) - 1)
        out[k] = (m, math.sqrt(var))
    return out


def check1():
    print("=" * 78)
    print("CHECK 1: new_multi_fstbi (rerun) == new_multi (pre-fix backup), per window")
    print("=" * 78)
    print("%-8s %-6s %8s %12s %14s %14s" %
          ("ploidy", "theta", "n_cmp", "n_mismatch", "max_abs_diff", "max_rel_diff"))
    allok = True
    for p, t, newpath in arms_present(NEW_DIR):
        oldpath = arm_path(OLD_DIR, p, t)
        if not os.path.exists(oldpath):
            continue
        oldmap = {(r["replicate"], r["chromosome"], r["window_pos_1"]): r
                  for r in load(oldpath)}
        ncmp = nmis = 0
        maxabs = maxrel = 0.0
        for r in load(newpath):
            o = oldmap.get((r["replicate"], r["chromosome"], r["window_pos_1"]))
            if o is None:
                continue
            a = fnum(r.get("fst_hudson_new_multi_fstbi"))
            b = fnum(o.get("fst_hudson_new_multi"))
            if a is None or b is None:
                continue
            ncmp += 1
            diff = abs(a - b)
            rel = diff / abs(b) if b else diff
            maxabs = max(maxabs, diff)
            maxrel = max(maxrel, rel)
            if rel > 1e-9 and diff > 1e-12:
                nmis += 1
        if nmis:
            allok = False
        print("%-8s %-6s %8d %12d %14.3e %14.3e" % (p, t, ncmp, nmis, maxabs, maxrel))
    print("\nVERDICT check 1:",
          "PASS - exact reproduction (differences are float noise)" if allok
          else "FAIL - real mismatches found")
    return allok


def check2():
    inf_fst = fst_infinite()
    print()
    print("=" * 78)
    print("CHECK 2: pooled FST at n=2000 vs each estimator's own expectation")
    print("=" * 78)
    print("E[FST]_infinite = %.9f  (tau=%.6f, alpha=%g)\n" % (inf_fst, TAU, ALPHA))
    hdr = ("%-6s %-6s %19s %19s" %
           ("ploid", "theta", "bi - E_inf", "multi - E_fin"))
    print(hdr)
    print("-" * len(hdr))
    for p, t, path in arms_present(NEW_DIR):
        rows = load(path)
        th = theta_of(t)
        b = boot_pooled(rows, ["new", "new_multi"])
        bi_m, bi_se = b["new"]
        mu_m, mu_se = b["new_multi"]
        print("%-6s %-6s %9.6f+-%-8.6f %9.6f+-%-8.6f" %
              (p, t, bi_m - inf_fst, bi_se, mu_m - fst_finite(th), mu_se))
    print("\nbi - E_inf   : biallelic minus infinite-sites expectation.")
    print("               Expect a small negative sag at high theta (residual")
    print("               homoplasy), growing with ploidy.")
    print("multi - E_fin: multiallelic minus finite-sites expectation. This is the")
    print("               derivation's out-of-sample test -- expect ~0 at every cell.")


def check3():
    print()
    print("=" * 78)
    print("CHECK 3: E[FST]_finite = 1 - E[pi_w]/E[dxy] against data")
    print("=" * 78)
    print("NOTE: the pi_* columns hold ONE population's pi, not the mean of pop1/pop2 --")
    print("      aggregate_one_rep.py keys on (chromosome, window_pos_1) and the last row")
    print("      wins, while pixy emits one pi row per population. The two populations are")
    print("      exchangeable so this stays unbiased for E[pi_w], but it uses half the data.\n")
    hdr = ("%-6s %-6s %10s %10s %10s %10s | %10s %10s %10s" %
           ("ploid", "theta", "pi_multi", "E[pi_w]", "dxy_multi", "E[dxy]",
            "1-pi/dxy", "E_fin", "diff"))
    print(hdr)
    print("-" * len(hdr))
    for p, t, path in arms_present(NEW_DIR):
        rows = load(path)
        th = theta_of(t)
        pim = mean_col(rows, "pi_new_multi")
        dxym = mean_col(rows, "dxy_new_multi")
        implied = 1 - pim / dxym
        efin = fst_finite(th)
        print("%-6s %-6s %10.6f %10.6f %10.6f %10.6f | %10.6f %10.6f %10.6f" %
              (p, t, pim, e_pi_w(th), dxym, e_dxy(th), implied, efin, implied - efin))


if __name__ == "__main__":
    have = arms_present(NEW_DIR)
    if not have:
        raise SystemExit("No arms found under %s -- sync the rerun first." % NEW_DIR)
    print("arms found: %d (%s)\n" %
          (len(have), ", ".join(sorted(set(p for p, _, _ in have)))))
    check1()
    check2()
    check3()
