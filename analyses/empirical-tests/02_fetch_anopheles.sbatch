#!/usr/bin/env bash
#SBATCH --job-name=fetch_anopheles
#SBATCH --partition=intel
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --array=1-16%8
#SBATCH --output=logs/fetch_anopheles_%A_%a.out
#SBATCH --error=logs/fetch_anopheles_%A_%a.err
# stream whole-chrom bam slices (chr3R) from ena Ag1000G phase 3
# for the 16 samples in config/pops_anopheles.tsv
set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}}"
ARM_ID=anopheles
DATASET_DIR="${PROJECT_DIR}/data/${ARM_ID}"
TARGETS_BED="${PROJECT_DIR}/config/regions/${ARM_ID}.targets.bed"
POPS_FILE="${PROJECT_DIR}/config/pops_${ARM_ID}.tsv"

# shellcheck disable=SC1091
source "${PROJECT_DIR}/lib/common.sh"
activate_env empirical_tests_fetch
require_tool samtools
require_tool curl

[[ -s "${TARGETS_BED}" ]] || die "Missing targets BED: ${TARGETS_BED}. Run 01_pick_regions.sh first."
[[ -s "${POPS_FILE}"   ]] || die "Missing pops file: ${POPS_FILE}"

# full-genome fasta decodes the streamed bams
REF_NAME=$(arm_config_lookup "${ARM_ID}" reference_name)
REF_FULL="${DATASET_DIR}/sources/${REF_NAME}.full.fa"
[[ -s "${REF_FULL}.fai" ]] || die "Missing full reference: ${REF_FULL}.fai (run 01_pick_regions.sh)"

mkdir -p "${DATASET_DIR}/bams" "${PROJECT_DIR}/logs"

TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"

# task 1 writes samples.tsv + bams.list; other tasks wait for it
write_manifests() {
    {
        printf "sample_id\tpopulation\talignment\n"
        awk 'BEGIN{FS=OFS="\t"} /^#/||/^sample_id/{next} {printf "%s\t%s\tdata/'"${ARM_ID}"'/bams/%s.bam\n",$1,$2,$1}' "${POPS_FILE}"
    } > "${DATASET_DIR}/samples.tsv"
    awk 'BEGIN{FS="\t"} /^#/||/^sample_id/{next} {printf "data/'"${ARM_ID}"'/bams/%s.bam\n",$1}' "${POPS_FILE}" \
        > "${DATASET_DIR}/bams.list"
}

if [[ "${TASK_ID}" -eq 1 ]]; then
    write_manifests
else
    for _ in $(seq 1 60); do
        [[ -s "${DATASET_DIR}/samples.tsv" ]] && break
        sleep 5
    done
fi

read -r SAMPLE POP SOURCE < <(
    awk -v t="${TASK_ID}" '
        /^#/ || /^sample_id/ { next }
        { n++; if (n == t) { print $1, $2, $3; exit } }
    ' "${POPS_FILE}"
)
[[ -n "${SAMPLE:-}" ]] || die "Task ${TASK_ID}: no sample on that row of ${POPS_FILE}"
info "Task ${TASK_ID}: ${SAMPLE} (${POP}) from ${SOURCE}"

OUT_BAM="${DATASET_DIR}/bams/${SAMPLE}.bam"
if [[ -s "${OUT_BAM}" && -s "${OUT_BAM}.bai" ]]; then
    info "[${SAMPLE}] BAM already present, skipping"
    exit 0
fi

# resolve ENA:<accession> to a bam url via the ena filereport api.
# prefer <sample_id>.bam (original Ag3.0 aln) over re-analyses like <sample_id>x.bam
resolve_ena_bam_url() {
    local accession="$1"
    local sample_id="$2"
    local report
    report=$(curl -fSL "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${accession}&result=analysis&fields=submitted_ftp,generated_ftp&format=tsv") \
        || die "ENA filereport API failed for ${accession}"
    python3 - "${report}" "${sample_id}" <<'PY'
import sys
lines = sys.argv[1].strip().splitlines()
sid = sys.argv[2]
if len(lines) < 2:
    sys.exit("ENA filereport returned no rows")
header = lines[0].split("\t")
bam_urls, bai_urls = [], []
for line in lines[1:]:
    row = dict(zip(header, line.split("\t")))
    for col in ("submitted_ftp", "generated_ftp"):
        for url in row.get(col, "").split(";"):
            if url.endswith(".bam"):
                bam_urls.append(url)
            elif url.endswith(".bam.bai"):
                bai_urls.append(url)
if not bam_urls:
    sys.exit(f"no .bam URL in ENA filereport for {sid}")
preferred = [u for u in bam_urls if u.endswith("/" + sid + ".bam")]
bam = preferred[0] if preferred else bam_urls[0]
bai = next((u for u in bai_urls if u == bam + ".bai"), "")
if not bai:
    sys.exit(f"no matching .bam.bai for {bam}")
if not bam.startswith("http"): bam = "https://" + bam
if not bai.startswith("http"): bai = "https://" + bai
print(bam)
print(bai)
PY
}

BAI_URL=""
if [[ "${SOURCE}" == ENA:* ]]; then
    accession="${SOURCE#ENA:}"
    mapfile -t _RESOLVED < <(resolve_ena_bam_url "${accession}" "${SAMPLE}")
    BAM_URL="${_RESOLVED[0]:-}"
    BAI_URL="${_RESOLVED[1]:-}"
    [[ -n "${BAM_URL}" && -n "${BAI_URL}" ]] || die "[${SAMPLE}] failed to resolve BAM/BAI URLs from ${accession}"
else
    BAM_URL="${SOURCE}"
    BAI_URL="${SOURCE}.bai"
fi
info "[${SAMPLE}] BAM URL: ${BAM_URL}"
info "[${SAMPLE}] BAI URL: ${BAI_URL}"

# build region args from targets bed (single row: whole chrom)
REGION_ARGS=()
while IFS=$'\t' read -r chrom start end _; do
    [[ -z "${chrom}" || "${chrom}" == \#* ]] && continue
    REGION_ARGS+=("${chrom}:$((start + 1))-${end}")
done < "${TARGETS_BED}"

# pre-fetch the .bai locally; htslib remote-index auto-discovery is flaky
LOCAL_BAI="${OUT_BAM}.remote.bai"
info "[${SAMPLE}] downloading remote BAM index to ${LOCAL_BAI}"
curl -fSL --retry 5 --retry-delay 10 --retry-connrefused \
    -o "${LOCAL_BAI}.tmp" "${BAI_URL}"
mv "${LOCAL_BAI}.tmp" "${LOCAL_BAI}"

info "[${SAMPLE}] streaming BAM slice (regions: ${REGION_ARGS[*]})"
TMP_BAM="${OUT_BAM}.tmp.$$"
samtools view -@ 4 -b \
    -T "${REF_FULL}" \
    -X \
    -o "${TMP_BAM}" \
    "${BAM_URL}" "${LOCAL_BAI}" \
    "${REGION_ARGS[@]}"
rm -f "${LOCAL_BAI}"

# force SM tag to sample_id so joint-call lines up
samtools view -H "${TMP_BAM}" | sed -E "s/(SM:)[^\t]+/\1${SAMPLE}/" \
    | samtools reheader - "${TMP_BAM}" > "${OUT_BAM}.tmp2"
mv "${OUT_BAM}.tmp2" "${OUT_BAM}"
rm -f "${TMP_BAM}"
samtools index -@ 4 "${OUT_BAM}"

info "[${SAMPLE}] done"
