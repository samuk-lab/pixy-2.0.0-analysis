#!/usr/bin/env bash
# per arm: write whole-chrom targets bed + extract single-chrom reference
# run from empirical_tests_fetch; .dict step shells out to gatk env
set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
ARMS_TSV="${PROJECT_DIR}/config/arms.tsv"
OUT_DIR="${PROJECT_DIR}/config/regions"
FORCE=0
[[ "${1:-}" == "--force" ]] && FORCE=1

[[ -s "${ARMS_TSV}" ]] || { echo "Missing ${ARMS_TSV}" >&2; exit 1; }
mkdir -p "${OUT_DIR}"

# shellcheck disable=SC1091
source "${PROJECT_DIR}/lib/common.sh"
activate_env empirical_tests_fetch
require_tool curl
require_tool samtools
require_tool bgzip

# vectorbase AgamP4 headers are AgamP4_2R; ena Ag1000G bams use bare 2R.
# strip the AgamP4_ prefix so bam and targets bed agree. arenosa
# (GCA_905216605.1) already uses bare accessions (LR999451.1), no rewrite
strip_chrom_prefix_anopheles() {
    sed -E 's/^>AgamP4_([0-9]+[LR]|X)\b/>\1/'
}

while IFS=$'\t' read -r arm_id species ploidy refname refurl chrom chrom_size \
                        pop_a pop_b pops_file n_per_pop window_size rng_seed; do
    [[ "${arm_id}" == "arm_id" || "${arm_id}" == \#* ]] && continue

    info "[${arm_id}] target chromosome: ${chrom} (${chrom_size} bp)"

    ##########
    # targets bed
    ##########
    bed="${OUT_DIR}/${arm_id}.targets.bed"
    if [[ -s "${bed}" && "${FORCE}" -eq 0 ]]; then
        info "[${arm_id}] targets BED exists, skipping: ${bed}"
    else
        printf '%s\t0\t%s\n' "${chrom}" "${chrom_size}" > "${bed}"
        info "[${arm_id}] wrote ${bed}"
    fi

    ##########
    # single-chrom reference fasta
    ##########
    ref_dir="${PROJECT_DIR}/data/${arm_id}/reference"
    src_dir="${PROJECT_DIR}/data/${arm_id}/sources"
    mkdir -p "${ref_dir}" "${src_dir}"

    # downloaded full ref; ncbi is .fna.gz, vectorbase is uncompressed .fasta
    full_basename=$(basename "${refurl}")
    full_cache="${src_dir}/${full_basename}"
    if [[ ! -s "${full_cache}" ]]; then
        info "[${arm_id}] downloading reference: ${refurl}"
        curl -fSL --retry 5 --retry-delay 10 -o "${full_cache}.tmp" "${refurl}"
        mv "${full_cache}.tmp" "${full_cache}"
    else
        info "[${arm_id}] reference cached: ${full_cache}"
    fi

    # decompress full fasta to a bare-chrom working copy for faidx extract
    full_fa="${src_dir}/${refname}.full.fa"
    if [[ ! -s "${full_fa}" || ! -s "${full_fa}.fai" ]]; then
        info "[${arm_id}] preparing full reference (${refname}) for chromosome extract"
        case "${full_basename}" in
            *.gz) gunzip_cmd="gunzip -c" ;;
            *)    gunzip_cmd="cat" ;;
        esac
        case "${arm_id}" in
            anopheles)
                ${gunzip_cmd} "${full_cache}" | strip_chrom_prefix_anopheles > "${full_fa}.tmp"
                ;;
            *)
                ${gunzip_cmd} "${full_cache}" > "${full_fa}.tmp"
                ;;
        esac
        mv "${full_fa}.tmp" "${full_fa}"
        samtools faidx "${full_fa}"
    fi

    # target chrom must be in the fai
    if ! awk -v c="${chrom}" 'BEGIN{FS="\t"} $1 == c {found=1} END{exit found?0:1}' \
              "${full_fa}.fai"; then
        echo "ERROR: chromosome '${chrom}' not in ${full_fa}.fai. Available chroms:" >&2
        cut -f1 "${full_fa}.fai" | head -20 >&2
        exit 1
    fi

    ##########
    # extract single chromosome
    ##########
    subset_bgz="${ref_dir}/${refname}.${chrom}.fa.gz"
    subset_fa="${ref_dir}/${refname}.${chrom}.fa"
    if [[ -s "${subset_bgz}" && -s "${subset_bgz}.fai" && -s "${subset_bgz}.gzi" \
          && -s "${subset_fa}" && -s "${subset_fa}.fai" \
          && "${FORCE}" -eq 0 ]]; then
        info "[${arm_id}] single-chrom reference already built: ${subset_bgz}"
    else
        info "[${arm_id}] extracting ${chrom} into ${subset_bgz}"
        samtools faidx "${full_fa}" "${chrom}" | bgzip -c > "${subset_bgz}.tmp"
        mv "${subset_bgz}.tmp" "${subset_bgz}"
        samtools faidx "${subset_bgz}"
        samtools faidx "${full_fa}" "${chrom}" > "${subset_fa}.tmp"
        mv "${subset_fa}.tmp" "${subset_fa}"
        samtools faidx "${subset_fa}"
    fi

    ##########
    # sequence dictionary (.dict) for gatk
    ##########
    dict="${subset_fa%.fa}.dict"
    if [[ -s "${dict}" && "${FORCE}" -eq 0 ]]; then
        info "[${arm_id}] sequence dictionary present: ${dict}"
    else
        info "[${arm_id}] creating sequence dictionary: ${dict}"
        if command -v gatk >/dev/null 2>&1; then
            gatk CreateSequenceDictionary -R "${subset_fa}" -O "${dict}"
        else
            # shell out to gatk env if gatk not on path
            mamba run -n empirical_tests_gatk gatk CreateSequenceDictionary \
                -R "${subset_fa}" -O "${dict}"
        fi
    fi

    info "[${arm_id}] done: targets=${bed} ref=${subset_fa}"
done < "${ARMS_TSV}"

info "[01_pick_regions] all arms done"
ls -1 "${OUT_DIR}"
