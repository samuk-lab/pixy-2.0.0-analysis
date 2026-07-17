#!/usr/bin/env bash
# submit pixy 2.0.0 fst-only array jobs across all core counts
# wraps 03_submit_pixy_new_all.sh with STATS=fst
# run after 03a patches new_pixy. 04_aggregate de-dupes by (version,stat,cores,seed)
# keeping last row, so reruns overwrite per seed
# vars: DEP, CORES, ARRAY_SPEC

set -euo pipefail

# slurm not on PATH in non-interactive ssh
export PATH="/opt/linux/rocky/8.x/x86_64/pkgs/slurm/24.11.1/bin:${PATH}"

PROJECT_DIR="${PROJECT_DIR:-$(cd "$(dirname "$0")" && pwd)}"

CORES="${CORES:-1 2 4 8 16}"
ARRAY_SPEC="${ARRAY_SPEC:-1-100%25}"
DEP="${DEP:-}"

RESULTS_DIR="${PROJECT_DIR}/data/results"

# clear new-pixy fst tsvs for a clean run; leave old-pixy fst untouched
echo "=== Clearing existing new-pixy FST result TSVs ===" >&2
for ncores in ${CORES}; do
  f="${RESULTS_DIR}/pixy_new_fst_10Mb_cores_${ncores}.tsv"
  if [[ -f "${f}" ]]; then
    rm -f "${f}"
    echo "  removed $(basename "${f}")" >&2
  fi
done

echo "=== Submitting pixy 2.0.0 FST benchmark (GitHub head) ===" >&2
echo "    Cores : ${CORES}" >&2
echo "    Array : ${ARRAY_SPEC}" >&2
[[ -n "${DEP}" ]] && echo "    Dep   : ${DEP}" >&2

STATS="fst" \
CORES="${CORES}" \
ARRAY_SPEC="${ARRAY_SPEC}" \
DEP="${DEP}" \
PROJECT_DIR="${PROJECT_DIR}" \
bash "${PROJECT_DIR}/03_submit_pixy_new_all.sh"
