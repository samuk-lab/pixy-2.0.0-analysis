#!/usr/bin/env bash
#SBATCH --job-name=call_arenosa
#SBATCH --partition=intel
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/call_arenosa_%j.out
#SBATCH --error=logs/call_arenosa_%j.err
# gatk joint-call, 16 arenosa on LR999451.1, autotetraploid, raw all-sites vcf
# mirrors 03 but ploidy=4; 05c_filter finalizes the filtered all-sites vcf
set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}}"
ARM_ID=arenosa
PLOIDY=4
DATASET_DIR="${PROJECT_DIR}/data/${ARM_ID}"
TARGETS_BED="${PROJECT_DIR}/config/regions/${ARM_ID}.targets.bed"
BAMS_LIST="${DATASET_DIR}/bams.list"

# shellcheck disable=SC1091
source "${PROJECT_DIR}/lib/common.sh"
activate_env empirical_tests_gatk
require_tool gatk
require_tool bcftools
require_tool samtools
require_tool parallel

REF_NAME=$(arm_config_lookup "${ARM_ID}" reference_name)
CHROM=$(arm_config_lookup    "${ARM_ID}" chrom)
REF_FASTA="${DATASET_DIR}/reference/${REF_NAME}.${CHROM}.fa"
DICT="${REF_FASTA%.fa}.dict"

[[ -s "${BAMS_LIST}"   ]] || die "Missing bams.list: ${BAMS_LIST} — run 04_fetch_arenosa.sbatch first"
[[ -s "${TARGETS_BED}" ]] || die "Missing targets BED: ${TARGETS_BED} — run 01_pick_regions.sh"
[[ -s "${REF_FASTA}"   ]] || die "Missing reference: ${REF_FASTA} — run 01_pick_regions.sh"
[[ -s "${REF_FASTA}.fai" && -s "${DICT}" ]] || die "Reference index/dict missing for ${REF_FASTA}"

missing=0
while IFS= read -r bam_rel; do
    bam_abs="${PROJECT_DIR}/${bam_rel}"
    if [[ ! -s "${bam_abs}" || ! -s "${bam_abs}.bai" ]]; then
        printf 'missing BAM or .bai: %s\n' "${bam_abs}" >&2
        missing=$((missing + 1))
    fi
done < "${BAMS_LIST}"
[[ "${missing}" -eq 0 ]] || die "${missing} BAM(s) missing — rerun 04_fetch_arenosa.sbatch"

THREADS="${SLURM_CPUS_PER_TASK:-16}"
N_SAMPLES=$(wc -l < "${BAMS_LIST}")

ABS_BAMS_LIST="${DATASET_DIR}/bams.abspath.list"
awk -v root="${PROJECT_DIR}/" '{print root $0}' "${BAMS_LIST}" > "${ABS_BAMS_LIST}"

VCF_DIR="${DATASET_DIR}/vcfs"
GVCF_DIR="${VCF_DIR}/gvcf"
TMP_DIR="${VCF_DIR}/tmp"
mkdir -p "${VCF_DIR}" "${GVCF_DIR}" "${TMP_DIR}"

RAW_VCF="${VCF_DIR}/${ARM_ID}.gatk.p${PLOIDY}.raw.vcf.gz"
COMBINED_GVCF="${VCF_DIR}/${ARM_ID}.combined.g.vcf.gz"

JAVA_OPTS_PER_SAMPLE="-Xmx4g -XX:ParallelGCThreads=1"

call_one_sample() {
    local bam_abs="$1"
    local sample
    sample=$(samtools view -H "${bam_abs}" | awk '/^@RG/{for(i=1;i<=NF;i++)if($i~/^SM:/){sub("SM:","",$i);print $i;exit}}')
    [[ -n "${sample}" ]] || { echo "[hc] no SM tag in ${bam_abs}" >&2; return 1; }
    local gvcf="${GVCF_DIR}/${sample}.g.vcf.gz"
    if [[ -s "${gvcf}" && -s "${gvcf}.tbi" ]]; then
        echo "[hc] ${sample}: GVCF exists, skipping" >&2
        return 0
    fi
    rm -f "${gvcf}"{,.tbi}
    gatk --java-options "${JAVA_OPTS_PER_SAMPLE}" HaplotypeCaller \
        -R "${REF_FASTA}" \
        -I "${bam_abs}" \
        -L "${TARGETS_BED}" \
        -O "${gvcf}" \
        --sample-ploidy "${PLOIDY}" \
        -ERC GVCF \
        --tmp-dir "${TMP_DIR}" \
        >&2
}
export -f call_one_sample
export GVCF_DIR TMP_DIR REF_FASTA TARGETS_BED PLOIDY JAVA_OPTS_PER_SAMPLE

info "[${ARM_ID}] HaplotypeCaller x ${N_SAMPLES} samples (${THREADS} workers, ploidy=${PLOIDY})"
parallel --jobs "${THREADS}" --halt now,fail=1 call_one_sample :::: "${ABS_BAMS_LIST}"

n_gvcfs=$(find "${GVCF_DIR}" -maxdepth 1 -name '*.g.vcf.gz' | wc -l)
[[ "${n_gvcfs}" -eq "${N_SAMPLES}" ]] || die "[${ARM_ID}] expected ${N_SAMPLES} GVCFs, found ${n_gvcfs}"
info "[${ARM_ID}] HaplotypeCaller done: ${n_gvcfs} GVCFs"

GVCF_ARGS=()
for g in "${GVCF_DIR}"/*.g.vcf.gz; do
    GVCF_ARGS+=(-V "${g}")
done

info "[${ARM_ID}] CombineGVCFs -> ${COMBINED_GVCF}"
rm -f "${COMBINED_GVCF}"{,.tbi}
gatk --java-options "-Xmx32g -XX:ParallelGCThreads=2" CombineGVCFs \
    -R "${REF_FASTA}" \
    "${GVCF_ARGS[@]}" \
    -O "${COMBINED_GVCF}" \
    --tmp-dir "${TMP_DIR}"

info "[${ARM_ID}] GenotypeGVCFs --sample-ploidy ${PLOIDY} -all-sites -> ${RAW_VCF}"
rm -f "${RAW_VCF}"{,.tbi}
gatk --java-options "-Xmx32g -XX:ParallelGCThreads=2" GenotypeGVCFs \
    -R "${REF_FASTA}" \
    -V "${COMBINED_GVCF}" \
    -L "${TARGETS_BED}" \
    -O "${RAW_VCF}" \
    --sample-ploidy "${PLOIDY}" \
    -all-sites \
    --tmp-dir "${TMP_DIR}"
[[ -s "${RAW_VCF}.tbi" ]] || bcftools index -t "${RAW_VCF}"

# drop combined gvcf; keep the raw vcf (05c_filter finalizes the all-sites vcf).
# NB 05c defers --trim-alt-alleles until after FORMAT/PL is stripped, so the
# ploidy-4 PL-length crash in bcftools 1.23.1 is still avoided
rm -f "${COMBINED_GVCF}"{,.tbi}

info "[${ARM_ID}] raw all-sites VCF: ${RAW_VCF}"
info "[${ARM_ID}] call done (run 05b_mappability + 05c_filter next)"
