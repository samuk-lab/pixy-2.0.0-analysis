#!/usr/bin/env bash
#SBATCH --job-name=fetch_arenosa
#SBATCH --partition=intel
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --array=1-16%8
#SBATCH --output=logs/fetch_arenosa_%A_%a.out
#SBATCH --error=logs/fetch_arenosa_%A_%a.err
# arenosa arm: ena fastq -> fastp -> bwa-mem2 -> sorted dedup bam, sliced to LR999451.1
# 16 samples (8 SPI + 8 TRE) from PRJNA484107 (Monnahan et al. 2019)
# bwa-mem2 indexes the full ~150 Mb ref; array throttled to 8 (%8), ~60 GB/task
set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}}"
ARM_ID=arenosa
DATASET_DIR="${PROJECT_DIR}/data/${ARM_ID}"
TARGETS_BED="${PROJECT_DIR}/config/regions/${ARM_ID}.targets.bed"
POPS_FILE="${PROJECT_DIR}/config/pops_${ARM_ID}.tsv"

# shellcheck disable=SC1091
source "${PROJECT_DIR}/lib/common.sh"
activate_env empirical_tests_fetch
require_tool samtools
require_tool bgzip
require_tool curl
require_tool bwa-mem2
require_tool fastp

[[ -s "${TARGETS_BED}" ]] || die "Missing targets BED: ${TARGETS_BED}. Run 01_pick_regions.sh first."
[[ -s "${POPS_FILE}"   ]] || die "Missing pops file: ${POPS_FILE}"

REF_NAME=$(arm_config_lookup "${ARM_ID}" reference_name)
# bwa-mem2 indexes the full fasta from 01_pick_regions.sh
REF_FULL="${DATASET_DIR}/sources/${REF_NAME}.full.fa"
[[ -s "${REF_FULL}.fai" ]] || die "Missing full reference: ${REF_FULL}.fai (run 01_pick_regions.sh)"

mkdir -p "${DATASET_DIR}/bams" "${DATASET_DIR}/fastq" "${DATASET_DIR}/tmp" "${PROJECT_DIR}/logs"

TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"

write_manifests() {
    {
        printf "sample_id\tpopulation\talignment\n"
        awk 'BEGIN{FS=OFS="\t"} /^#/||/^sample_id/{next} {printf "%s\t%s\tdata/'"${ARM_ID}"'/bams/%s.bam\n",$1,$2,$1}' "${POPS_FILE}"
    } > "${DATASET_DIR}/samples.tsv"
    awk 'BEGIN{FS="\t"} /^#/||/^sample_id/{next} {printf "data/'"${ARM_ID}"'/bams/%s.bam\n",$1}' "${POPS_FILE}" \
        > "${DATASET_DIR}/bams.list"
}

build_bwa_index() {
    if [[ -s "${REF_FULL}.bwt.2bit.64" ]]; then
        info "[${ARM_ID}] BWA-MEM2 index already present"
        return
    fi
    info "[${ARM_ID}] building BWA-MEM2 index on full reference (~3 min for ~150 Mb)"
    bwa-mem2 index "${REF_FULL}"
}

if [[ "${TASK_ID}" -eq 1 ]]; then
    write_manifests
    build_bwa_index
else
    for _ in $(seq 1 360); do
        if [[ -s "${DATASET_DIR}/samples.tsv" && -s "${REF_FULL}.bwt.2bit.64" ]]; then
            break
        fi
        sleep 10
    done
fi

read -r SAMPLE POP RUN_ACC FASTQ_FTP < <(
    awk -v t="${TASK_ID}" '
        BEGIN{FS="\t"}
        /^#/ || /^sample_id/ { next }
        { n++; if (n == t) { print $1, $2, $3, $4; exit } }
    ' "${POPS_FILE}"
)
[[ -n "${SAMPLE:-}" ]] || die "Task ${TASK_ID}: no sample on that row of ${POPS_FILE}"
info "Task ${TASK_ID}: ${SAMPLE} (${POP}) run=${RUN_ACC}"

OUT_BAM="${DATASET_DIR}/bams/${SAMPLE}.bam"
if [[ -s "${OUT_BAM}" && -s "${OUT_BAM}.bai" ]]; then
    info "[${SAMPLE}] BAM already present, skipping"
    exit 0
fi

# resolve paired fastq urls; PRJNA484107 gives 3 per run, pick the _1/_2 pair
IFS=';' read -ra FQ_URLS <<< "${FASTQ_FTP}"
if [[ "${#FQ_URLS[@]}" -lt 2 ]]; then
    info "[${SAMPLE}] resolving FASTQ paths via ENA filereport API"
    REPORT=$(curl -fsSL "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${RUN_ACC}&result=read_run&fields=fastq_ftp&format=tsv")
    FASTQ_FTP=$(printf '%s\n' "${REPORT}" | awk 'NR==2{print $1}')
    IFS=';' read -ra FQ_URLS <<< "${FASTQ_FTP}"
fi
[[ "${#FQ_URLS[@]}" -ge 2 ]] || die "[${SAMPLE}] could not resolve paired FASTQ URLs"

FQ1_URL=""
FQ2_URL=""
for url in "${FQ_URLS[@]}"; do
    case "${url}" in
        *_1.fastq.gz|*_1.fq.gz) FQ1_URL="${url}" ;;
        *_2.fastq.gz|*_2.fq.gz) FQ2_URL="${url}" ;;
    esac
done
if [[ -z "${FQ1_URL}" || -z "${FQ2_URL}" ]]; then
    if [[ "${#FQ_URLS[@]}" -eq 2 ]]; then
        FQ1_URL="${FQ_URLS[0]}"; FQ2_URL="${FQ_URLS[1]}"
    else
        die "[${SAMPLE}] could not identify paired _1/_2 FASTQ URLs in: ${FASTQ_FTP}"
    fi
fi
[[ "${FQ1_URL}" == http* ]] || FQ1_URL="https://${FQ1_URL}"
[[ "${FQ2_URL}" == http* ]] || FQ2_URL="https://${FQ2_URL}"
info "[${SAMPLE}] paired URLs: R1=${FQ1_URL##*/} R2=${FQ2_URL##*/}"

WORK="${DATASET_DIR}/tmp/${SAMPLE}"
mkdir -p "${WORK}"
FQ1="${WORK}/${SAMPLE}.R1.fq.gz"
FQ2="${WORK}/${SAMPLE}.R2.fq.gz"

robust_download() {
    local url="$1" dest="$2" max="${3:-5}"
    local n=0
    while (( n < max )); do
        (( n++ )) || true
        info "[${SAMPLE}] download attempt ${n}/${max}: ${url##*/}"
        curl -fSL --retry 3 --retry-delay 10 --retry-connrefused \
            -C - -o "${dest}.tmp" "${url}" || {
            info "[${SAMPLE}] curl failed on attempt ${n}; will retry"
            rm -f "${dest}.tmp"
            sleep $(( n * 30 ))
            continue
        }
        if gzip -t "${dest}.tmp" 2>/dev/null; then
            mv -f "${dest}.tmp" "${dest}"
            return 0
        fi
        info "[${SAMPLE}] gzip integrity check failed on attempt ${n} — deleting and retrying"
        rm -f "${dest}.tmp"
        sleep $(( n * 30 ))
    done
    die "[${SAMPLE}] download failed after ${max} attempts: ${url}"
}

info "[${SAMPLE}] downloading FASTQs"
{ [[ -s "${FQ1}" ]] && gzip -t "${FQ1}" 2>/dev/null; } || { rm -f "${FQ1}"; robust_download "${FQ1_URL}" "${FQ1}"; }
{ [[ -s "${FQ2}" ]] && gzip -t "${FQ2}" 2>/dev/null; } || { rm -f "${FQ2}"; robust_download "${FQ2_URL}" "${FQ2}"; }

info "[${SAMPLE}] fastp adapter/quality trim"
TRIM1="${WORK}/${SAMPLE}.R1.trim.fq.gz"
TRIM2="${WORK}/${SAMPLE}.R2.trim.fq.gz"
fastp \
    --in1 "${FQ1}" --in2 "${FQ2}" \
    --out1 "${TRIM1}" --out2 "${TRIM2}" \
    --json "${WORK}/fastp.json" --html "${WORK}/fastp.html" \
    --thread 4 --detect_adapter_for_pe

info "[${SAMPLE}] BWA-MEM2 align + fixmate + sort"
RAW_BAM="${WORK}/${SAMPLE}.sorted.bam"
bwa-mem2 mem -t 8 \
    -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SAMPLE}\tPL:ILLUMINA" \
    "${REF_FULL}" "${TRIM1}" "${TRIM2}" \
    | samtools fixmate -m -O bam - - \
    | samtools sort -@ 4 -m 1G -O bam -o "${RAW_BAM}" -

info "[${SAMPLE}] markdup"
MARKED_BAM="${WORK}/${SAMPLE}.markdup.bam"
samtools markdup -@ 4 --mode s "${RAW_BAM}" "${MARKED_BAM}"
samtools index -@ 4 "${MARKED_BAM}"

info "[${SAMPLE}] slicing to target regions"
REGION_ARGS=()
while IFS=$'\t' read -r chrom start end _; do
    [[ -z "${chrom}" || "${chrom}" == \#* ]] && continue
    REGION_ARGS+=("${chrom}:$((start + 1))-${end}")
done < "${TARGETS_BED}"

TMP_BAM="${OUT_BAM}.tmp.$$"
samtools view -@ 4 -b -o "${TMP_BAM}" "${MARKED_BAM}" "${REGION_ARGS[@]}"
mv "${TMP_BAM}" "${OUT_BAM}"
samtools index -@ 4 "${OUT_BAM}"

rm -f "${FQ1}" "${FQ2}" "${TRIM1}" "${TRIM2}" "${RAW_BAM}" "${MARKED_BAM}" "${MARKED_BAM}.bai"

info "[${SAMPLE}] done"
