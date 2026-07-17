#!/usr/bin/env bash
#SBATCH --job-name=run_pixy
#SBATCH --partition=intel
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --array=1-4
#SBATCH --output=logs/run_pixy_%A_%a.out
#SBATCH --error=logs/run_pixy_%A_%a.err
# array job: 2 arms x 2 variants = 4 tasks (biallelic / +multiallelic)
# each task runs pixy pi+dxy then fst; outputs under data/<arm>/pixy_out/<variant>/
# window 10 kb from arms.tsv; fst wc for diploid, hudson for ploidy > 2
# (pixy 2.0.0 refuses wc on non-diploid)
set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}}"

# shellcheck disable=SC1091
source "${PROJECT_DIR}/lib/common.sh"
activate_env empirical_tests_pixy
require_tool pixy
require_tool bcftools

TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"
case "${TASK_ID}" in
    1) ARM_ID=anopheles; VARIANT=biallelic ;;
    2) ARM_ID=anopheles; VARIANT=multi     ;;
    3) ARM_ID=arenosa;   VARIANT=biallelic ;;
    4) ARM_ID=arenosa;   VARIANT=multi     ;;
    *) die "Unknown SLURM_ARRAY_TASK_ID=${TASK_ID}; expected 1..4" ;;
esac

WINDOW_SIZE=$(arm_config_lookup "${ARM_ID}" window_size)
PLOIDY=$(arm_config_lookup "${ARM_ID}" ploidy)
if [[ "${PLOIDY}" -eq 2 ]]; then
    FST_TYPE=wc
else
    FST_TYPE=hudson
fi
DATASET_DIR="${PROJECT_DIR}/data/${ARM_ID}"
VCF="${DATASET_DIR}/vcfs/${ARM_ID}.all_sites.vcf.gz"
POPS_SRC="${PROJECT_DIR}/config/pops_${ARM_ID}.tsv"
PIXY_POPS="${DATASET_DIR}/pops_${ARM_ID}.pixy.tsv"
OUT_DIR="${DATASET_DIR}/pixy_out/${VARIANT}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

[[ -s "${VCF}"      ]] || die "Missing all-sites VCF: ${VCF} - run 0{3,5}_call_*.sbatch first"
[[ -s "${VCF}.tbi"  ]] || die "Missing tabix index: ${VCF}.tbi"
[[ -s "${POPS_SRC}" ]] || die "Missing pops file: ${POPS_SRC}"

# emit pixy 2-column pops file (sample_id<TAB>population)
awk 'BEGIN{FS=OFS="\t"} /^#/||/^sample_id/{next} {print $1, $2}' "${POPS_SRC}" > "${PIXY_POPS}"
N_SAMPLES=$(wc -l < "${PIXY_POPS}")
info "[${ARM_ID}/${VARIANT}] pops: ${PIXY_POPS} (${N_SAMPLES} samples, window=${WINDOW_SIZE}, ploidy=${PLOIDY}, fst_type=${FST_TYPE})"

mkdir -p "${OUT_DIR}"

EXTRA=()
if [[ "${VARIANT}" == "multi" ]]; then
    EXTRA+=(--include_multiallelic_snps)
fi

PREFIX="${ARM_ID}_${VARIANT}"

##########
# pi + dxy
##########
info "[${ARM_ID}/${VARIANT}] pixy pi dxy"
pixy --stats pi dxy \
    --vcf          "${VCF}" \
    --populations  "${PIXY_POPS}" \
    --window_size  "${WINDOW_SIZE}" \
    --n_cores      "${THREADS}" \
    --output_folder "${OUT_DIR}" \
    --output_prefix "${PREFIX}_pidxy" \
    "${EXTRA[@]}"

##########
# fst (wc diploid, hudson ploidy > 2)
##########
info "[${ARM_ID}/${VARIANT}] pixy fst (${FST_TYPE})"
pixy --stats fst \
    --vcf          "${VCF}" \
    --populations  "${PIXY_POPS}" \
    --window_size  "${WINDOW_SIZE}" \
    --n_cores      "${THREADS}" \
    --fst_type     "${FST_TYPE}" \
    --output_folder "${OUT_DIR}" \
    --output_prefix "${PREFIX}_fst" \
    "${EXTRA[@]}"

info "[${ARM_ID}/${VARIANT}] done. Outputs in ${OUT_DIR}/"
ls -la "${OUT_DIR}/"
