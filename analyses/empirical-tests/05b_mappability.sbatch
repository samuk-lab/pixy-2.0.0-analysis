#!/usr/bin/env bash
#SBATCH --job-name=mappability
#SBATCH --partition=intel
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=1-2
#SBATCH --output=logs/mappability_%A_%a.out
#SBATCH --error=logs/mappability_%A_%a.err
# genmap uniquely-mappable mask per arm (K=100, E=1, keep mappability==1)
# array: 1=anopheles 2=arenosa; writes config/regions/<arm>.mappable.bed
# consumed by 05c_filter.sbatch as a site-level mask (variant + invariant alike)
set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}}"

# shellcheck disable=SC1091
source "${PROJECT_DIR}/lib/common.sh"
activate_env empirical_tests_mappability
require_tool genmap
require_tool bedtools
require_tool awk

# mappability params (short-read population genomics defaults)
KMER="${KMER:-100}"
ERRORS="${ERRORS:-1}"
MAP_MIN="${MAP_MIN:-1.0}"

TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"
case "${TASK_ID}" in
    1) ARM_ID=anopheles ;;
    2) ARM_ID=arenosa   ;;
    *) die "Unknown SLURM_ARRAY_TASK_ID=${TASK_ID}; expected 1..2" ;;
esac

REF_NAME=$(arm_config_lookup "${ARM_ID}" reference_name)
CHROM=$(arm_config_lookup    "${ARM_ID}" chrom)
DATASET_DIR="${PROJECT_DIR}/data/${ARM_ID}"
REF_FASTA="${DATASET_DIR}/reference/${REF_NAME}.${CHROM}.fa"
OUT_DIR="${PROJECT_DIR}/config/regions"
WORK="${DATASET_DIR}/mappability"
IDX_DIR="${WORK}/genmap_index"
MAP_PREFIX="${WORK}/${ARM_ID}.K${KMER}_E${ERRORS}"
MASK_BED="${OUT_DIR}/${ARM_ID}.mappable.bed"

[[ -s "${REF_FASTA}"     ]] || die "Missing reference: ${REF_FASTA} - run 01_pick_regions.sh"
[[ -s "${REF_FASTA}.fai" ]] || die "Missing reference index: ${REF_FASTA}.fai"
mkdir -p "${OUT_DIR}" "${WORK}"

##########
# genmap index (dir must not pre-exist)
##########
if [[ -d "${IDX_DIR}" ]]; then
    info "[${ARM_ID}] genmap index exists, reusing: ${IDX_DIR}"
else
    info "[${ARM_ID}] genmap index -> ${IDX_DIR}"
    genmap index -F "${REF_FASTA}" -I "${IDX_DIR}"
fi

##########
# genmap map -> per-base mappability bedgraph
##########
info "[${ARM_ID}] genmap map -K ${KMER} -E ${ERRORS} (bedgraph)"
genmap map -K "${KMER}" -E "${ERRORS}" \
    -I "${IDX_DIR}" \
    -O "${MAP_PREFIX}" \
    -T "${SLURM_CPUS_PER_TASK:-8}" \
    -bg

##########
# threshold mappability >= MAP_MIN, merge to mask bed
##########
info "[${ARM_ID}] threshold mappability >= ${MAP_MIN} -> ${MASK_BED}"
awk -v m="${MAP_MIN}" 'BEGIN{FS=OFS="\t"} $4 >= m {print $1, $2, $3}' "${MAP_PREFIX}.bedgraph" \
    | bedtools merge -i - > "${MASK_BED}.tmp"
mv "${MASK_BED}.tmp" "${MASK_BED}"

n_bp=$(awk 'BEGIN{FS="\t"; s=0} {s += $3 - $2} END{print s}' "${MASK_BED}")
info "[${ARM_ID}] mappable bases (K=${KMER},E=${ERRORS},>=${MAP_MIN}): ${n_bp} on ${CHROM}"
info "[${ARM_ID}] mappability mask done: ${MASK_BED}"
