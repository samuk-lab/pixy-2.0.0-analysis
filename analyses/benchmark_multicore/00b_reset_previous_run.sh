#!/usr/bin/env bash
# wipe every product of a previous benchmark run so the next one starts clean.
# the per-cell summaries are appended to under flock with a header-if-absent
# guard, so rerunning on top of them silently doubles every cell -- this must be
# run before resubmitting, not after.
#
# simulated VCFs are KEPT (they are expensive and vcfsim is deterministic per
# seed), except the two affected by the seed-collision bug fixed in
# 01_simulate_vcfs.sbatch. 01 skips seeds whose .vcf.gz + .tbi already exist, so
# deleting them here is what makes them regenerate.
#
# dry run by default; pass --yes to actually delete.

set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
DATA_DIR="${PROJECT_DIR}/data"
VCF_DIR="${DATA_DIR}/vcfs/dm_10Mb"

# seeds whose VCFs were byte-identical (11 was built from 1's simulation)
COLLIDED_SEEDS=(1 11)

apply=0
[[ "${1:-}" == "--yes" ]] && apply=1

rm_path() {
  local p="$1"
  [[ -e "${p}" ]] || return 0
  if (( apply )); then
    rm -rf "${p}"
    echo "  removed  ${p#"${PROJECT_DIR}/"}"
  else
    echo "  would rm ${p#"${PROJECT_DIR}/"}"
  fi
}

echo "== results and per-cell summaries"
rm_path "${DATA_DIR}/results"

echo "== pixy working output"
for d in pixy_old_out pixy_new_out pixy_mem_out pixy_mem_old_out zarr zarr_mem; do
  rm_path "${DATA_DIR}/${d}"
done

echo "== logs"
rm_path "${PROJECT_DIR}/logs"

echo "== simulated VCFs (kept, except the collided seeds)"
if [[ -d "${VCF_DIR}" ]]; then
  for seed in "${COLLIDED_SEEDS[@]}"; do
    rm_path "${VCF_DIR}/dm_sim_vcf_seed_${seed}.vcf.gz"
    rm_path "${VCF_DIR}/dm_sim_vcf_seed_${seed}.vcf.gz.tbi"
  done
  # unrenamed intermediates from tasks that died mid-simulation, and the
  # per-task work tree; both are regenerable and neither should be reused
  rm_path "${VCF_DIR}/_work"
  while IFS= read -r stray; do
    rm_path "${stray}"
  done < <(find "${VCF_DIR}" -maxdepth 1 -name '*.vcf' -print)

  kept=$(find "${VCF_DIR}" -maxdepth 1 -name '*.vcf.gz' | wc -l)
  echo "  ${kept} VCFs currently present (100 expected after 01 backfills)"
else
  echo "  ${VCF_DIR#"${PROJECT_DIR}/"} does not exist; 01 will simulate all seeds"
fi

if (( apply )); then
  echo
  echo "Reset complete. Next: sbatch 01_simulate_vcfs.sbatch (backfills the deleted seeds)."
else
  echo
  echo "Dry run — nothing deleted. Re-run with --yes to apply."
fi
