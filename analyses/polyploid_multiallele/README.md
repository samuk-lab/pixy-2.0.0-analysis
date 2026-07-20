# polyploid_multiallele → Figures 3, 4, S2

Simulation sweeps testing pixy 2.x estimators across ploidy (2n, 4n, 6n, 8n),
missingness (0–75%), and θ, using vcfsim, pixy, old pixy, and vcftools. The summary
tables land in `analysis/` and feed the figure scripts in `../../figures/`.

## Overview

Each "arm" is a ploidy × condition cell defined in `config/sim_params.tsv`. An arm
simulates many replicate VCFs, runs the estimators on each replicate, and aggregates
them into a per-arm table; the summary scripts then reduce those tables to what the
figures read. Figure 3 shows estimator bias against finite-sites expectations across
missingness; Figure 4 compares the biallelic and multiallelic-aware estimators across
θ; Figure S2 is a diploid cross-version check of old pixy against 2.0.0.

| Summary table | Figure |
|---|---|
| `missingness_summary.tsv`, `thetaw_tajimaD_summary.tsv`, `per_rep_*.tsv`, `neutral_tajimaD_reference.tsv` | Figure 3 |
| `multiallelic_frac_summary.tsv`, `dxy_fst_theta_summary.tsv` | Figure 4 |
| `per_seed_diploid.tsv` | Figure S2 |

## Scripts

Numbered in execution order; `.sbatch`/`.sh` run under SLURM, the R and Python
scripts run natively.

| Script | Purpose |
|--------|---------|
| `00_create_conda_envs.sh` | Build envs in `envs/` (`pixy`, `old_pixy`, `vcfsim`, `vcftools`) |
| `01_run_all.sh` | Submit `02_run_arm.sh` as an array, one job per arm in `config/sim_params.tsv` |
| `02_run_arm.sh` | For one arm: simulate VCFs, run pixy (+ old pixy, + vcftools), join each replicate with `analysis/aggregate_one_rep.py` |
| `03_concat_and_analyze.sh` | Concatenate replicates into per-arm TSVs under `data/aggregated/` |
| `04_run_R_analyses.sh` | Run the four summary scripts below that build the figure tables |
| `analysis/03_missingness_sweep.R` | → `missingness_summary.tsv` (Figures 3, S2) |
| `analysis/07_thetaw_tajimaD.R` | → `thetaw_tajimaD_summary.tsv` (Figure 3) |
| `analysis/08_multiallelic_frac.R` | → `multiallelic_frac_summary.tsv` (Figure 4) |
| `analysis/11_dxy_fst_theta_sweep.R` | → `dxy_fst_theta_summary.tsv` (Figure 4) |
| `analysis/per_seed_diploid.R` | → `per_seed_diploid.tsv` (Figure S2); WC-F<sub>ST</sub> for new pixy backfilled by `analysis/merge_wcfst.py` |
| `analysis/calibrate_neutral_tajimaD.py` | Finite-sites Tajima's D reference line for Figure 3 (uses `../tajima_multiallelic_variance/`) |
| `analysis/04_ploidy_grid.R` | → `ploidy_grid_summary.tsv`, the provenance of a hardcoded expected-d<sub>xy</sub> constant in `Inline_statistical_tests.R` |
| `update_pixy_env.sbatch` | Repoint the pixy env to the `biallelic-fst` build for the F<sub>ST</sub> rerun; checked by `analysis/validate_fst_rerun.py` |

## Notes

Old pixy runs with `--fst_maf_filter 0` (see `02_run_arm.sh`) so its Weir–Cockerham
F<sub>ST</sub> is comparable to new pixy; the 0.95.01 default of 0.05 biases it upward
otherwise. The FST rerun produced the genuine multiallelic Hudson-F<sub>ST</sub>
columns in Figure 4; the finite-sites E[F<sub>ST</sub>] derivation is in
`analysis/finite_sites_fst_derivation.md`.

## Dependencies

Conda envs in `envs/`: `pixy`, `old_pixy`, `vcfsim`, `vcftools`. The aggregated
per-arm TSVs (~380 MB) are not tracked; regenerate them with steps `01`–`03`.
