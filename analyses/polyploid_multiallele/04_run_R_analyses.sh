#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G
#SBATCH --job-name=pixy_R_analyses
#SBATCH -p batch
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x.%j.out

# regenerate the manuscript summary tables. runs the four figure-feeding summary
# scripts against data/aggregated/, writing the *_summary.tsv files that
# ../../../figures/ consume, plus sanity-check pdfs in analysis/figs/.
#
#   sbatch 04_run_R_analyses.sh        # under slurm
#   ./04_run_R_analyses.sh             # locally
#
# Rscript from PATH by default; override with RSCRIPT=/path/to/Rscript
#
# qc / diagnostic scripts (01_baseline_check, 02_diagnose_residual, 04_ploidy_grid,
# 09_vcftools_compare) live in analysis/diagnostics/ and are not run here.
# 02_diagnose_residual.R needs a raw per-rep pixy output the runner deletes by
# default (comment out the final rm -rf "$rep_work" in 02_run_arm.sh, run one
# replicate, then invoke it with the saved rep_<seed>_pi.txt).

set -euo pipefail

cd "$(dirname "$(readlink -f "$0")")/analysis"

RSCRIPT="${RSCRIPT:-Rscript}"

ts() { date '+%F %T'; }

echo "[$(ts)] R analyses: starting in $(pwd)"
echo "[$(ts)] Rscript: $($RSCRIPT --version 2>&1 | head -1)"

for script in 03_missingness_sweep.R \
              07_thetaw_tajimaD.R \
              08_multiallelic_frac.R \
              11_dxy_fst_theta_sweep.R; do
    echo
    echo "[$(ts)] === Running $script ==="
    "$RSCRIPT" "$script"
    echo "[$(ts)] === Done: $script ==="
done

echo
echo "[$(ts)] All R analyses complete. Figures:"
ls figs/*.pdf
