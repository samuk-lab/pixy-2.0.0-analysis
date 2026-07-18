#!/usr/bin/env bash
#SBATCH --job-name=count_multi
#SBATCH --partition=short
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=logs/count_multi_%j.out
#SBATCH --error=logs/count_multi_%j.err
# count multiallelic SNP sites (>= 3 alleles) per 10 kb window, directly from the
# canonical filtered all_sites VCF (05c output). independent of pixy's pi/fst
# accounting; used for the empirical uplift ~ multiallelic-density Spearman rho
# reported in the Results (see figures/Inline_statistical_tests.R, section 7-8).
set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}}"
# shellcheck disable=SC1091
source "${PROJECT_DIR}/lib/common.sh"
activate_env empirical_tests_gatk
require_tool bcftools

OUT="${PROJECT_DIR}/data/aggregated"
mkdir -p "${OUT}"

for arm in anopheles arenosa; do
    vcf="${PROJECT_DIR}/data/${arm}/vcfs/${arm}.all_sites.vcf.gz"
    [[ -s "${vcf}" ]] || die "Missing filtered VCF: ${vcf} - run 05c_filter.sbatch"
    # -m3: >= 3 alleles (multiallelic); -v snps: SNPs only. bin to pixy 10 kb window start.
    bcftools view -H -m3 -v snps "${vcf}" \
        | gawk 'BEGIN{OFS="\t"} {w=int(($2-1)/10000)*10000+1; c[w]++} END{for(x in c) print x, c[x]}' \
        | sort -k1,1n > "${OUT}/${arm}.multi_per_window.tsv"
    info "[${arm}] $(wc -l < "${OUT}/${arm}.multi_per_window.tsv") windows with multiallelic SNPs"
done
