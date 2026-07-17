# tajima_multiallelic_variance → supplementary note, Figure 3 reference line

Theory and validation for the mutation-count Watterson's θ / Tajima's D estimator that
pixy 2.x uses in multiallelic-aware mode, and the source of the finite-sites Tajima's D
expectation drawn as the reference line in Figure 3.

## Overview

In multiallelic-aware mode π counts allelic pairwise differences — a site with three or
four alleles contributes its full heterozygosity — while classic Watterson's θ<sub>W</sub>
counts segregating sites. Under finite-sites JC69 mutation the two diverge, and Tajima's
`D = (π − θ_W) / sd` is biased positive, with the bias growing as more lineages are
sampled. This is the positive D seen across ploidies in Figure 3. The fix generalises
θ<sub>W</sub> to count mutations, `M = Σ (k − 1)` (Tajima 1996's parsimony count `s*`);
every biallelic site has `k − 1 = 1`, so the result matches stock pixy on biallelic data.

`variance_theory.md` is the full write-up; `finite_sites_expectations_writeup.md` covers
the finite-sites expectations.

## Scripts

| Script | Purpose |
|--------|---------|
| `eta_tajima.py` | Mutation-count θ<sub>W</sub> + D (`eta_stats`). Imported by `../polyploid_multiallele/analysis/calibrate_neutral_tajimaD.py` to build the Figure 3 D reference line |
| `variance_validate.py` | Variance-term validation (r=0 and r=μ, with an infinite-sites control) → `results/variance_sweep.tsv` |
| `eta_validate.py` | Re-simulation driver for the ploidy sweep (ploidy-aware VCF reader) |
| `aggregate_eta.py` | Per-ploidy means and quantiles from the re-simulation output |
| `tajima1996_expectations.py` | Analytic JC69 E(s), E(s\*), E(π), checked against Tajima 1996 Table 1 |
| `cluster/eta_array.sbatch` | 4-ploidy SLURM array for the re-simulation |

## Run the variance validation

```sh
PYTHONPATH=/path/to/pixy N_REPS=2000 python variance_validate.py
```

The infinite-sites arm is a positive control: η ≡ S there, so `D_old == D_new` must
hold exactly. Compare `sd(D_new)` in the JC69 arm against that control (≈0.88 at r=0),
not against 1.0 — this is an estimator correction, not a calibrated test statistic.

> `allel.read_vcf` defaults to diploid and silently truncates polyploid genotypes;
> always pass `numbers={"calldata/GT": ploidy}`. An earlier round of this analysis was
> retracted over that bug, and its Route A/B and self-normalisation scripts have been
> dropped from this copy.

## Dependencies

Python with `scikit-allel`, `msprime`, `numpy`, and `pandas`, plus an importable copy
of the pixy branch under test (`PYTHONPATH`). The per-replicate raw cells under `data/`
(~4 GB) are not tracked; regenerate them with `cluster/eta_array.sbatch`.
