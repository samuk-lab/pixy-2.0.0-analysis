#!/usr/bin/env bash
# local (wsl) driver for 02_run_arm.sh — loops the array tasks in parallel, no slurm
#
# usage: analysis/run_local_arm.sh ARM [N_ARRAY_TASKS] [PARALLEL] [N_CORES]
#   ARM            arm name in config/sim_params.tsv
#   N_ARRAY_TASKS  tasks to split n_replicates into (default 100)
#   PARALLEL       tasks to run at once (default 15)
#   N_CORES        --n_cores passed to pixy per task (default 1)
#
# resolves the project dir relative to this script
set -euo pipefail

ARM="${1:?usage: run_local_arm.sh ARM [N_ARRAY_TASKS] [PARALLEL] [N_CORES]}"
NT="${2:-100}"
PAR="${3:-15}"
NC="${4:-1}"

PROJ="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJ"
export PATH="$HOME/miniconda3/bin:$PATH"

echo "[$(date '+%F %T')] local run: arm=$ARM tasks=$NT parallel=$PAR n_cores=$NC"
rm -rf "data/_work/${ARM}" 2>/dev/null || true
rm -f data/"${ARM}".part_* 2>/dev/null || true
mkdir -p logs

seq 1 "$NT" | xargs -P "$PAR" -I{} bash -c '
  ARM_NAME="'"$ARM"'" N_ARRAY_TASKS="'"$NT"'" SLURM_ARRAY_TASK_ID="$1" N_CORES="'"$NC"'" \
    bash 02_run_arm.sh > "logs/'"$ARM"'.task_$1.log" 2>&1
' _ {}

n_parts=$(ls data/"${ARM}".part_*.tsv 2>/dev/null | wc -l)
echo "[$(date '+%F %T')] done $ARM: ${n_parts} partials"
