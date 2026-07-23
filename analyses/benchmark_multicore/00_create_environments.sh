#!/usr/bin/env bash
# build the three conda envs via sbatch (run before 05_submit_all_waves.sh)
# --chain also submits the pipeline launcher with afterok deps

set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
LOG_DIR="${PROJECT_DIR}/logs"
CHAIN=false
[[ "${1:-}" == "--chain" ]] && CHAIN=true

mkdir -p "${LOG_DIR}"

# find sbatch, fall back to hpcc path
SBATCH="$(command -v sbatch 2>/dev/null || \
          echo /opt/linux/rocky/8.x/x86_64/pkgs/slurm/24.11.1/bin/sbatch)"
[[ -x "${SBATCH}" ]] || { echo "sbatch not found at ${SBATCH}" >&2; exit 1; }

submit_env() {
  local name="$1" yml="${PROJECT_DIR}/envs/$2"
  [[ -f "${yml}" ]] || { echo "Missing ${yml}" >&2; exit 1; }
  "${SBATCH}" \
    --job-name="env_${name}" \
    --partition=short \
    --ntasks=1 --cpus-per-task=1 \
    --mem=16G --time=01:00:00 \
    --output="${LOG_DIR}/env_${name}_%j.out" \
    --error="${LOG_DIR}/env_${name}_%j.err" \
    --wrap="source ~/.bashrc
            # clear the cached conda-forge index first: a stale index predating a
            # freshly published pixy release once made the solve fail *after* the old
            # env was already removed, destroying it (see RERUN_HANDOFF run log).
            mamba clean --index-cache -y 2>/dev/null || true
            mamba env remove -n '${name}' -y 2>/dev/null || true
            mamba env create -f '${yml}'" \
  | awk '{print $NF}'
}

echo "=== Creating conda environments ==="
vcfsim_jid=$(submit_env vcfsim vcfsim.yml);       echo "  vcfsim    → job ${vcfsim_jid}"
pixy_old_jid=$(submit_env old_pixy pixy_old.yml); echo "  old_pixy  → job ${pixy_old_jid}"
pixy_new_jid=$(submit_env new_pixy pixy_new.yml); echo "  new_pixy  → job ${pixy_new_jid}"

all_jids="${vcfsim_jid}:${pixy_old_jid}:${pixy_new_jid}"
echo ""
echo "All env jobs submitted. Monitor with:"
echo "  squeue -u \$USER -j ${vcfsim_jid},${pixy_old_jid},${pixy_new_jid}"

if [[ "${CHAIN}" == "true" ]]; then
  echo ""
  echo "=== Chaining pipeline launcher (afterok:${all_jids}) ==="
  launch_jid=$("${SBATCH}" \
    --job-name=launch_benchmark_v2 \
    --partition=short \
    --ntasks=1 --cpus-per-task=1 \
    --mem=1G --time=00:15:00 \
    --dependency="afterok:${all_jids}" \
    --output="${LOG_DIR}/launch_benchmark_%j.out" \
    --error="${LOG_DIR}/launch_benchmark_%j.err" \
    --wrap="source ~/.bashrc
            export PROJECT_DIR='${PROJECT_DIR}'
            cd '${PROJECT_DIR}'
            bash 05_submit_all_waves.sh" \
    | awk '{print $NF}')
  echo "  Pipeline launcher → job ${launch_jid}"
  echo ""
  echo "Pipeline will start automatically once all envs are ready."
else
  echo ""
  echo "When all three env jobs complete successfully, run:"
  echo "  bash ${PROJECT_DIR}/05_submit_all_waves.sh"
  echo ""
  echo "Or rerun with --chain to submit the pipeline launcher automatically:"
  echo "  bash ${PROJECT_DIR}/00_create_environments.sh --chain"
fi
