#!/usr/bin/env python
"""is the eta-substituted tajima's d variance term calibrated at true polyploidy?

the mutation-count (eta = sum(k-1)) estimator has two externally validated arms
and one unvalidated:

  * theta_W  -> matches Tajima 1996 closed-form E(s*)/a1 to ~1e-5 across ploidy
                (analysis/thetaw_tajimaD_summary.tsv)
  * pi       -> matches tskit's unbiased pairwise pi to ~1e-16
  * variance -> calc_tajima_d_stdev computes e1*eta + e2*eta*(eta-1), but the
                Tajima 1989 constants are derived for S segregating sites under
                infinite sites. substituting eta treats two mutations sharing a
                site as two independent units. never externally checked.

the D floor in neutral_tajimaD_reference.tsv can't check it: that floor comes from
eta_stats, the same variance machinery under test, so "pipeline D matches the floor"
is an implementation check not a calibration one.

the one prior direct test (results/run.log, sd(D) at r=0) is invalid: it came from
analyze.py, which called allel.read_vcf without numbers=, silently truncating every
polyploid genotype to its first 2 alleles. every "polyploid" cell there was diploid --
visible as a total absence of the ploidy trend the correct re-sim shows (D_old flat
~-0.01 across 2n-8n vs +0.069 -> +0.180 in eta_validate.py). see handoff.md.

design: drive msprime directly to a GenotypeArray. no VCF is written, so the truncating
reader can't apply and the ~0.01 vcfsim-VCF D residual (handoff.md) is out of the loop.

two recombination regimes. r=0 (one non-recombining tree per locus) is the only regime
where Tajima's analytic variance is derived and sd(D)=1 is the right target -- the clean
test of the substitution. r=mu matches the figure 2/3 arms and measures how far
recombination alone moves sd(D); it can't discriminate between estimators, since the
analytic variance is wrong there for any tajima's d.

the pre-existing figures (sd ~0.85 at r=0, ~0.19 at r=mu, from results/run.log) came from
the ploidy-truncated run above and at a different theta (0.025 vs the figure 2/3 regime's
4*NE*mu = 0.0378). this script supersedes them: sd ~0.88 at r=0, ~0.16 at r=mu, at true
polyploidy and the figures' own parameters.

two mutation models per ploidy:

  * infinite  -> discrete_genome=False. no homoplasy, every k=2, so eta == S and D_new is
                 arithmetically identical to D_old. positive control: if sd(D) is not ~1
                 here the harness (or Tajima's variance at this theta) is the problem, not
                 the eta substitution.
  * jc69      -> discrete_genome=True + JC69. recurrent mutation -> multiallelic sites.
                 the actual test: how far does sd(D_new) drift from the infinite-sites
                 control as multiallelic fraction rises with n?

reading the result: compare sd(D_new) against sd(D_old) in the jc69 arm and against the
infinite-sites control. the control absorbs whatever finite-theta departure from sd=1
Tajima's variance already has, isolating the eta effect.

env: pixy-py311-test (pixy + scikit-allel + msprime), PYTHONPATH -> pixy checkout.
override reps with N_REPS.
"""

from __future__ import annotations

import os
import sys
from collections import Counter

import allel
import msprime
import numpy as np
from allel import GenotypeArray

from pixy.calc import _harmonic_sum, calc_tajima_d_stdev

# matched to config/sim_params.tsv tajimaD arms (the Figure 2/3 regime)
NE = 1_720_600
MU = 5.49e-9
WINDOW_LEN = 10_000
N_IND = 10
PLOIDIES = [2, 4, 6, 8]
N_REPS = int(os.environ.get("N_REPS", "2000"))
SEED = 20260715

# recombination regimes. r=0 is where Tajima's variance is derived and sd(D)=1 is the right
# target. r=mu matches the figure 2/3 arms, only measures how far recombination moves sd(D)
RECOMB = {"r0": 0.0, "rmu": MU}


def both_stats(gt: GenotypeArray):
    """tajima's d by site counting (master) and mutation counting (branch).

    both share pixy's _harmonic_sum / calc_tajima_d_stdev and the same pi, so the
    only difference is what goes into the per-n counter. None if the window has no
    usable variation.
    """
    ac = gt.count_alleles()
    if int(np.count_nonzero(ac.sum(axis=1))) == 0:
        return None

    raw_pi = float(np.sum(allel.mean_pairwise_difference(ac=ac, fill=0)))
    variant = ac[ac[:, 1:].sum(axis=1) != 0]
    if variant.shape[0] == 0:
        return None

    n_per = variant.sum(axis=1)
    k_per = np.count_nonzero(variant, axis=1)

    # pixy master: one tally per variant row, keyed by n. counts monomorphic-for-alt
    # rows (k==1) as segregating -- what the old figure numbers did
    site_counts: Counter = Counter(int(x) for x in n_per.tolist())
    # branch: sum of (k-1) over rows with k>=2
    eta_counts: Counter = Counter()
    for n_i, k_i in zip(n_per.tolist(), k_per.tolist(), strict=True):
        if k_i >= 2:
            eta_counts[int(n_i)] += int(k_i) - 1

    # denominator inputs, shared by both arms, so the only difference between D_old and D_new
    # is the counting currency (sites vs mutations). denominator is pixy's post-#160 form: one
    # evaluation at the pooled mean n (see calc_tajima_d_stdev). with no missing data every
    # site has the same n, so this equals a per-class evaluation; only matters under ragged
    # missingness, which this script doesn't simulate
    per_site_n = ac.sum(axis=1)
    total_allele_count = int(per_site_n[per_site_n > 0].sum())
    num_sites = int(np.count_nonzero(per_site_n))

    def theta_and_d(counter):
        th = 0.0
        for n_i, c in counter.items():
            a1 = float(_harmonic_sum(int(n_i)))
            if a1 > 0:
                th += c / a1
        sd = calc_tajima_d_stdev(
            total_allele_count=total_allele_count,
            num_sites=num_sites,
            num_mutations=int(sum(counter.values())),
        )
        d = (raw_pi - th) / sd if sd > 0 else float("nan")
        return th, d

    th_old, d_old = theta_and_d(site_counts)
    th_new, d_new = theta_and_d(eta_counts)

    n_var = int(variant.shape[0])
    return {
        "D_old": d_old,
        "D_new": d_new,
        "thetaW_old": th_old,
        "thetaW_new": th_new,
        "n_var_sites": float(n_var),
        # fraction of variant sites with >2 alleles
        "multi_frac": float(np.count_nonzero(k_per > 2) / n_var),
        "max_k": float(k_per.max()),
    }


def run_cell(ploidy, model, n_reps, recomb_label="r0"):
    """one (ploidy x mutation model x recombination) cell, returns per-rep rows."""
    rng = np.random.default_rng(
        SEED + 1000 * ploidy + (0 if model == "infinite" else 7)
        + (0 if recomb_label == "r0" else 13)
    )
    n_chrom = N_IND * ploidy

    # ploidy=1 with population_size=2*NE reproduces the vcfsim coalescent timescale
    # (theta = 4*NE*mu, constant across organism ploidy; only n_chrom varies). see
    # calibrate_neutral_tajimaD.py. recombination_rate=0 -> one tree per locus
    anc = msprime.sim_ancestry(
        samples=n_chrom,
        ploidy=1,
        population_size=2 * NE,
        sequence_length=WINDOW_LEN,
        recombination_rate=RECOMB[recomb_label],
        num_replicates=n_reps,
        random_seed=int(rng.integers(1, 2**31)),
    )

    rows = []
    for ts in anc:
        if model == "infinite":
            mts = msprime.sim_mutations(
                ts, rate=MU, discrete_genome=False,
                random_seed=int(rng.integers(1, 2**31)),
            )
        else:
            mts = msprime.sim_mutations(
                ts, rate=MU, model=msprime.JC69(), discrete_genome=True,
                random_seed=int(rng.integers(1, 2**31)),
            )
        if mts.num_sites == 0:
            continue
        mat = mts.genotype_matrix().astype(np.int8)
        gt = GenotypeArray(mat[:, :, None])
        s = both_stats(gt)
        if s is None or not np.isfinite(s["D_old"]) or not np.isfinite(s["D_new"]):
            continue
        rows.append(s)
    return rows


def summarize(ploidy, model, rows, recomb_label="r0"):
    d_old = np.array([r["D_old"] for r in rows])
    d_new = np.array([r["D_new"] for r in rows])
    mf = np.array([r["multi_frac"] for r in rows])
    mk = np.array([r["max_k"] for r in rows])
    n = len(rows)
    # se(sd) ~ sd/sqrt(2(n-1)) for an approximately normal sample
    se_sd_old = d_old.std(ddof=1) / np.sqrt(2 * (n - 1))
    se_sd_new = d_new.std(ddof=1) / np.sqrt(2 * (n - 1))
    print(
        f"ploidy={ploidy} n_chrom={N_IND*ploidy} model={model} recomb={recomb_label} reps={n}\n"
        f"  multiallelic frac of variant sites: {mf.mean():.4f}  (mean max k = {mk.mean():.2f})\n"
        f"  D_old : mean={d_old.mean():+.4f}  sd={d_old.std(ddof=1):.4f} +/- {se_sd_old:.4f}\n"
        f"  D_new : mean={d_new.mean():+.4f}  sd={d_new.std(ddof=1):.4f} +/- {se_sd_new:.4f}",
        flush=True,
    )
    return {
        "ploidy": ploidy, "n_chrom": N_IND * ploidy, "model": model,
        "recomb": recomb_label, "reps": n,
        "multi_frac": mf.mean(), "mean_max_k": mk.mean(),
        "mean_D_old": d_old.mean(), "sd_D_old": d_old.std(ddof=1), "se_sd_old": se_sd_old,
        "mean_D_new": d_new.mean(), "sd_D_new": d_new.std(ddof=1), "se_sd_new": se_sd_new,
    }


def main():
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "results", "variance_sweep.tsv")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    summaries = []
    for recomb_label in ["r0", "rmu"]:
        for model in ["infinite", "jc69"]:
            for ploidy in PLOIDIES:
                rows = run_cell(ploidy, model, N_REPS, recomb_label)
                if not rows:
                    print(f"ploidy={ploidy} model={model} recomb={recomb_label}: no usable reps",
                          file=sys.stderr)
                    continue
                summaries.append(summarize(ploidy, model, rows, recomb_label))
                print(flush=True)

    cols = list(summaries[0].keys())
    with open(out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for s in summaries:
            fh.write("\t".join(str(s[c]) for c in cols) + "\n")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
