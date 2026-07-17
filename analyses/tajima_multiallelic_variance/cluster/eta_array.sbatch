#!/bin/bash
#SBATCH --job-name=etafig2
#SBATCH --array=1-4
#SBATCH --partition=epyc
#SBATCH --time=03:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/eta_%A_%a.out
#SBATCH --error=logs/eta_%A_%a.err
#
# re-simulate figure-2 arms (one ploidy per task), pixy multiallelic tajima's d /
# watterson theta vs the eta estimator
#   sbatch --array=1-4 cluster/eta_array.sbatch
#   override: sbatch --export=ALL,N_OBS=1500 ...
set -eo pipefail

PROJECT_DIR="${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/.." && pwd)}}"
cd "$PROJECT_DIR"

ENVS=/bigdata/samuklab/ksamuk/.conda/envs
export VCFSIM_PY="$ENVS/vcfsim/bin/python"
export VCFSIM_SRC=""
export ROUTEB_TMP="${SLURM_TMPDIR:-/tmp}"
PIXY_PY="$ENVS/pixy/bin/python"

case "$SLURM_ARRAY_TASK_ID" in
  1) P=2 ;;
  2) P=4 ;;
  3) P=6 ;;
  4) P=8 ;;
  *) echo "bad task id $SLURM_ARRAY_TASK_ID"; exit 1 ;;
esac

mkdir -p results/eta_raw logs
RAW="results/eta_raw/cell_${P}.raw.tsv"
: > "$RAW"

echo "[task $SLURM_ARRAY_TASK_ID] ploidy=$P n_obs=${N_OBS:-1500}"
"$PIXY_PY" eta_validate.py --ploidy "$P" --n_obs "${N_OBS:-1500}" \
    --seedbase "$((P * 1000003))" --raw "$RAW"
echo "[task $SLURM_ARRAY_TASK_ID] done -> $RAW"
