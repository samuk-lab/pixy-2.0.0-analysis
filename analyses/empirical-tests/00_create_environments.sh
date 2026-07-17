#!/usr/bin/env bash
# submit sbatch jobs to build the three conda envs for empirical-tests
# --chain also submits 01_pick_regions.sh afterok on the fetch env
set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
LOG_DIR="${PROJECT_DIR}/logs"
mkdir -p "${LOG_DIR}"

CHAIN=false
[[ "${1:-}" == "--chain" ]] && CHAIN=true

# find sbatch, fall back to hpcc path
SBATCH="$(command -v sbatch 2>/dev/null || \
          echo /opt/linux/rocky/8.x/x86_64/pkgs/slurm/24.11.1/bin/sbatch)"
[[ -x "${SBATCH}" ]] || { echo "sbatch not found at ${SBATCH}" >&2; exit 1; }

submit_env() {
  # args: env name, yml filename under envs/
  local name="$1" yml="${PROJECT_DIR}/envs/$2"
  [[ -f "${yml}" ]] || { echo "Missing ${yml}" >&2; exit 1; }
  "${SBATCH}" \
    --job-name="env_${name}" \
    --partition=intel \
    --ntasks=1 --cpus-per-task=2 \
    --mem=16G --time=01:30:00 \
    --output="${LOG_DIR}/env_${name}_%j.out" \
    --error="${LOG_DIR}/env_${name}_%j.err" \
    --wrap="set -eo pipefail
            set +u; source ~/.bashrc; set -u
            mamba env remove -n '${name}' -y 2>/dev/null || true
            mamba env create -f '${yml}' -n '${name}'
            # pixy env: bump dask for numpy 2.x (scikit-allel pins old dask)
            if [ '${name}' = 'empirical_tests_pixy' ]; then
                mamba run -n '${name}' pip install 'dask[array]>=2024.1.0' --upgrade --quiet
            fi" \
  | awk '{print $NF}'
}

echo "=== Creating conda environments ==="
fetch_jid=$(submit_env empirical_tests_fetch fetch.yml)
echo "  empirical_tests_fetch  -> job ${fetch_jid}"
gatk_jid=$(submit_env empirical_tests_gatk  gatk.yml)
echo "  empirical_tests_gatk   -> job ${gatk_jid}"
pixy_jid=$(submit_env empirical_tests_pixy  pixy.yml)
echo "  empirical_tests_pixy   -> job ${pixy_jid}"

all_jids="${fetch_jid}:${gatk_jid}:${pixy_jid}"
echo ""
echo "Monitor: squeue -u \$USER -j ${fetch_jid},${gatk_jid},${pixy_jid}"

if [[ "${CHAIN}" == "true" ]]; then
  echo ""
  echo "=== Chaining 01_pick_regions.sh (afterok:${fetch_jid}) ==="
  pick_jid=$("${SBATCH}" \
    --job-name=pick_regions \
    --partition=short \
    --ntasks=1 --cpus-per-task=4 \
    --mem=8G --time=01:00:00 \
    --dependency="afterok:${fetch_jid}" \
    --output="${LOG_DIR}/pick_regions_%j.out" \
    --error="${LOG_DIR}/pick_regions_%j.err" \
    --wrap="set +u; source ~/.bashrc; set -u
            cd '${PROJECT_DIR}'
            bash 01_pick_regions.sh" \
    | awk '{print $NF}')
  echo "  pick_regions -> job ${pick_jid}"
else
  echo ""
  echo "When fetch env is ready, run:"
  echo "  bash ${PROJECT_DIR}/01_pick_regions.sh"
fi
