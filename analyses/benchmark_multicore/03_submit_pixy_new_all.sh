#!/usr/bin/env bash
# submit pixy 2.0.0 array jobs across all (stat x cores) combos
# called from 05_submit_all_waves.sh or standalone
# progress -> stderr; one job id per line -> stdout for dep chaining
# vars: DEP, STATS, CORES, ARRAY_SPEC

set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}}"
LOG_DIR="${PROJECT_DIR}/logs"

STATS="${STATS:-pi dxy fst}"
CORES="${CORES:-1 2 4 8 16}"
ARRAY_SPEC="${ARRAY_SPEC:-1-100%25}"
DEP="${DEP:-}"

mkdir -p "${LOG_DIR}"

script="${PROJECT_DIR}/03_pixy_new_array.sbatch"
[[ -f "${script}" ]] || { echo "Missing ${script}" >&2; exit 1; }

for ncores in ${CORES}; do
  for stat in ${STATS}; do
    run_tag="cores_${ncores}"
    dep_args=()
    [[ -n "${DEP}" ]] && dep_args=(--dependency="${DEP}")

    # mem ceilings, not measured; re-tune from seff. peak rss in all_cells_long.tsv
    case "${ncores}" in
      1|2)  mem_mb=8192  ;;   #  8G
      4|8)  mem_mb=12288 ;;   # 12G
      *)    mem_mb=16384 ;;   # 16G (16 cores)
    esac

    jid=$(sbatch \
      --array="${ARRAY_SPEC}" \
      --job-name="pixy_new_${stat}_${run_tag}" \
      --cpus-per-task="${ncores}" \
      --mem="${mem_mb}M" \
      "${dep_args[@]+"${dep_args[@]}"}" \
      --export=ALL,PROJECT_DIR="${PROJECT_DIR}",STAT="${stat}",N_CORES="${ncores}",RUN_TAG="${run_tag}" \
      "${script}" | awk '{print $NF}')

    echo "  [submit] pixy_new ${stat} cores=${ncores} → job ${jid}" >&2
    echo "${jid}"   # stdout: one job id per line
  done
done
