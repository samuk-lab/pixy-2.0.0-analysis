#!/usr/bin/env python
"""finite-sites neutral expectation of the multiallelic-corrected Tajima's D,
per (ploidy x missingness) cell of the Figure-3 sweep.

the reference line for D is the value the corrected estimator attains on
finite-sites (homoplasy-generating) neutral data, not the infinite-sites ideal
of 0. the mutation-count eta = sum(k-1) estimator removes the site-vs-mutation
counting bias, but under JC69 homoplasy eta undercounts recurrent mutations,
leaving a small negative floor (~ -0.01 to -0.04, deepening with missingness as
the effective sample size drops). see ../../tajima_multiallelic_variance/handoff.md.

per cell, replays vcfsim's post-processing on neutral msprime windows:
  * percent_missing_genotypes -> per site, k_ind = round(pct/100 * n_ind) whole
    individuals (all `ploidy` gene copies) set missing. count fixed per site, so
    the called size drops to a constant (n_ind - k_ind) * ploidy gene copies.
    moves E[D] (vcfsim SimulatorClass: hap_idx = s*ploidy .. (s+1)*ploidy = -1).
  * percent_missing_sites -> whole sites dropped ~Bernoulli(pct); fewer sites per
    window shrinks |D| toward 0.

windows are finite-sites JC69 msprime (discrete genome, recurrent mutation ->
multiallelic sites), matched to the vcfsim window (10 kb, r = mu, Ne, mu, 10
individuals x ploidy gene copies), scored with eta_tajima.eta_stats. per-cell
mean of D_eta is the `theoretical` reference consumed by 07_thetaw_tajimaD.R.

genealogies are drawn ONCE per ploidy and shared across the 5 missingness cells,
so the missingness shape is clean but the absolute level is one MC draw. default
N_WINDOWS = 80000 (MCSE ~5e-4), at which the floor converges onto the pipeline's
own neutral mean.

modes:
  * python calibrate_neutral_tajimaD.py           -> all 4 ploidies serial,
    writes merged neutral_tajimaD_reference.tsv
  * python calibrate_neutral_tajimaD.py <ploidy>  -> one ploidy, writes partial
    neutral_tajimaD_reference.p<ploidy>.tsv (parallel; see run_calibrate_*.sh)
  seeded per-ploidy (SEED + 1000*ploidy), so serial and parallel agree.

env: pixy + scikit-allel + msprime. override the window count with N_WINDOWS.
"""
from __future__ import annotations

import os
import sys

import allel
import msprime
import numpy as np

sys.path.insert(0, "../../tajima_multiallelic_variance")
from eta_tajima import eta_stats  # noqa: E402  (project's corrected eta-D)

# matched to config/sim_params.tsv tajimaD arms
NE = 1_720_600
MU = 5.49e-9
WINDOW_LEN = 10_000          # 100 kb sequence / 10 windows per replicate
N_IND = 10                   # diploid-equivalent samples per population (sample_size)
PLOIDIES = [2, 4, 6, 8]      # -> n_chrom = 20, 40, 60, 80 sampled gene copies
N_WINDOWS = int(os.environ.get("N_WINDOWS", "80000"))
SEED = 20260628

# Figure-3 missingness cells: (miss_pct, percent_missing_sites,
# percent_missing_genotypes) as split in config/sim_params.tsv
CELLS = [
    (0.00,  0,  0),
    (0.10,  5,  5),
    (0.25, 12, 13),
    (0.50, 25, 25),
    (0.75, 37, 38),
]

HEADER = "ploidy\tmiss_pct\tn_chrom\tn_chrom_eff\tmean_D\tse_D\tn_windows\n"


def apply_vcfsim_missing(mat, ploidy, n_ind, pct_site, pct_geno, rng):
    """replay vcfsim's missing-data post-processing on a window's allele matrix.

    mat is (n_sites, n_chrom) int allele indices, columns grouped by individual
    (gene copies s*ploidy .. (s+1)*ploidy). returns a copy with dropped sites
    removed and masked genotypes set to -1 (allel's missing sentinel).
    """
    # missing sites: drop variant sites ~Bernoulli(pct_site)
    if pct_site > 0 and mat.shape[0] > 0:
        keep = rng.random(mat.shape[0]) >= (pct_site / 100.0)
        mat = mat[keep]
    if mat.shape[0] == 0:
        return mat

    # missing genotypes: fixed k_ind whole individuals per site
    k_ind = round((pct_geno / 100.0) * n_ind)  # round() = vcfsim's rounding
    if k_ind > 0:
        mat = mat.copy()
        n_sites = mat.shape[0]
        rand = rng.random((n_sites, n_ind))
        kth = min(k_ind, n_ind - 1)
        miss_ind = np.argpartition(rand, kth, axis=1)[:, :k_ind]  # (n_sites, k_ind)
        for j in range(ploidy):
            cols = miss_ind * ploidy + j
            np.put_along_axis(mat, cols, -1, axis=1)
    return mat


def run_ploidy(ploidy, n_windows=N_WINDOWS):
    """finite-sites D floor for one ploidy across all missingness cells.

    seeded per ploidy (SEED + 1000*ploidy) so serial and parallel agree.
    genealogies are simulated ONCE and every missingness mask applied to the
    same draws (clean missingness shape).
    """
    rng = np.random.default_rng(SEED + 1000 * ploidy)
    n_chrom = N_IND * ploidy

    # ploidy=1 keeps exactly n_chrom haploid gene-copy columns (per-individual
    # missingness masking stays valid). use population_size=2*NE, NOT NE, to match
    # the vcfsim timescale: msprime's top-level ploidy scales the coalescent (rate
    # 1/(ploidy*N)); vcfsim leaves it at 2 (it sets ploidy only in the SampleSet,
    # which controls sample count not timescale), so its timescale is 2*NE ->
    # theta=4*NE*mu for every organism ploidy. ploidy=1,population_size=2*NE
    # reproduces that; timescale is constant across ploidy, only n_chrom varies.
    anc = msprime.sim_ancestry(
        samples=n_chrom,
        ploidy=1,
        population_size=2 * NE,
        sequence_length=WINDOW_LEN,
        recombination_rate=MU,
        num_replicates=n_windows,
        random_seed=int(rng.integers(1, 2**31)),
    )
    mats = []
    for ts in anc:
        # finite-sites JC69: discrete genome + recurrent mutation -> multiallelic sites
        mts = msprime.sim_mutations(
            ts, rate=MU, model=msprime.JC69(), discrete_genome=True,
            random_seed=int(rng.integers(1, 2**31)),
        )
        if mts.num_sites == 0:
            continue
        mats.append(mts.genotype_matrix().astype(np.int8))  # (n_sites, n_chrom)

    rows = []
    for miss_pct, pct_site, pct_geno in CELLS:
        mask_rng = np.random.default_rng(SEED + ploidy)  # same masks reproducibly
        ds = []
        for mat in mats:
            m = apply_vcfsim_missing(mat, ploidy, N_IND, pct_site, pct_geno, mask_rng)
            if m.shape[0] == 0:
                continue
            gt = allel.GenotypeArray(m[:, :, None])  # ploidy-1 GT, -1 = missing
            s = eta_stats(gt)
            if s is not None and np.isfinite(s["tajimaD_eta"]):
                ds.append(s["tajimaD_eta"])
        d = np.asarray(ds)
        k_ind = round((pct_geno / 100.0) * N_IND)
        n_chrom_eff = (N_IND - k_ind) * ploidy
        se = d.std(ddof=1) / np.sqrt(len(d)) if len(d) > 1 else float("nan")
        rows.append((ploidy, miss_pct, n_chrom, n_chrom_eff, d.mean(), se, len(d)))
        print(f"ploidy {ploidy} miss {miss_pct:.0%} (n_eff={n_chrom_eff}): "
              f"mean D_eta = {d.mean():.4f} +/- {se:.4f}  (n_windows={len(d)})",
              flush=True)
    return rows


def write_rows(path, rows):
    with open(path, "w") as fh:
        fh.write(HEADER)
        for ploidy, miss_pct, n_chrom, n_eff, m, se, nw in rows:
            fh.write(f"{ploidy}\t{miss_pct}\t{n_chrom}\t{n_eff}\t{m:.6f}\t{se:.6f}\t{nw}\n")


def main():
    if len(sys.argv) > 1:
        # single-ploidy worker -> partial file (for parallel runs)
        ploidy = int(sys.argv[1])
        rows = run_ploidy(ploidy)
        out = f"neutral_tajimaD_reference.p{ploidy}.tsv"
        write_rows(out, rows)
        print(f"wrote {out}", flush=True)
    else:
        # all ploidies (serial) -> merged reference
        rows = []
        for ploidy in PLOIDIES:
            rows.extend(run_ploidy(ploidy))
        write_rows("neutral_tajimaD_reference.tsv", rows)
        print("wrote neutral_tajimaD_reference.tsv", flush=True)


if __name__ == "__main__":
    main()
