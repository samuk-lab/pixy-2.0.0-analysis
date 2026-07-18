#!/usr/bin/env bash
#SBATCH --job-name=filter_vcf
#SBATCH --partition=intel
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --array=1-2
#SBATCH --output=logs/filter_vcf_%A_%a.out
#SBATCH --error=logs/filter_vcf_%A_%a.err
# build a filtered all-sites vcf from the raw gatk vcf, pixy-safe:
#   1. genotype DP below MIN_DP -> missing (denominator stays correct)
#   2. drop indels + SNPs within SNP_GAP bp of an indel (near-indel misalignment
#      is a top source of spurious multiallelic SNPs); keep SNP + invariant only
#   3. gatk hard filters + ExcessHet on variant sites; failing SITES dropped (invariants kept)
#   4. mappability mask from 05b (variant + invariant alike)
#   5. site depth band: drop sites with INFO/DP outside [qlo,qhi] quantiles
#   6. max-missing: each pop must have >= MIN_PRESENT called genotypes
#   7. strip FORMAT annots + trim unused alt alleles (as 03/05)
# steps 2,4-6 are symmetric over variant + invariant sites, so pixy denominators
# stay honest. array: 1=anopheles(2n) 2=arenosa(4n). reuses *.gatk.p*.raw.vcf.gz
# (INFO DP/QD/FS/MQ/SOR/ExcessHet + FORMAT DP), NOT a re-call. this is the SOLE producer
# of the canonical <arm>.all_sites.vcf.gz that 06_run_pixy consumes (03/05 emit raw only)
set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}}"

# shellcheck disable=SC1091
source "${PROJECT_DIR}/lib/common.sh"
activate_env empirical_tests_gatk
require_tool gatk
require_tool bcftools
require_tool gawk

# thresholds (gatk best-practice snp hard filters; DP genotype floor)
MIN_DP="${MIN_DP:-5}"
QD_MIN="${QD_MIN:-2.0}"
FS_MAX="${FS_MAX:-60.0}"
MQ_MIN="${MQ_MIN:-40.0}"
MQRS_MIN="${MQRS_MIN:--12.5}"
RPRS_MIN="${RPRS_MIN:--8.0}"
SOR_MAX="${SOR_MAX:-3.0}"
EXCESSHET_MAX="${EXCESSHET_MAX:-54.69}"
# indel handling: drop indels and SNPs within SNP_GAP bp of an indel
SNP_GAP="${SNP_GAP:-5}"
# site depth band quantiles (on INFO/DP, total site depth) and max-missing
DP_QLO="${DP_QLO:-0.05}"
DP_QHI="${DP_QHI:-0.95}"
MIN_PRESENT="${MIN_PRESENT:-6}"

TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"
case "${TASK_ID}" in
    1) ARM_ID=anopheles; PLOIDY=2 ;;
    2) ARM_ID=arenosa;   PLOIDY=4 ;;
    *) die "Unknown SLURM_ARRAY_TASK_ID=${TASK_ID}; expected 1..2" ;;
esac

REF_NAME=$(arm_config_lookup "${ARM_ID}" reference_name)
CHROM=$(arm_config_lookup    "${ARM_ID}" chrom)
N_PER_POP=$(arm_config_lookup "${ARM_ID}" n_per_pop)
DATASET_DIR="${PROJECT_DIR}/data/${ARM_ID}"
REF_FASTA="${DATASET_DIR}/reference/${REF_NAME}.${CHROM}.fa"
VCF_DIR="${DATASET_DIR}/vcfs"
TMP_DIR="${VCF_DIR}/tmp"
RAW_VCF="${VCF_DIR}/${ARM_ID}.gatk.p${PLOIDY}.raw.vcf.gz"
MASK_BED="${PROJECT_DIR}/config/regions/${ARM_ID}.mappable.bed"
POPS_SRC="${PROJECT_DIR}/config/pops_${ARM_ID}.tsv"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

GT_VCF="${VCF_DIR}/${ARM_ID}.gtfilt.vcf.gz"
SNP_VCF="${VCF_DIR}/${ARM_ID}.snp.vcf.gz"
FLAG_VCF="${VCF_DIR}/${ARM_ID}.hardflag.vcf.gz"
SITEPASS_VCF="${VCF_DIR}/${ARM_ID}.sitepass.vcf.gz"
MASK_VCF="${VCF_DIR}/${ARM_ID}.masked.vcf.gz"
DPBAND_VCF="${VCF_DIR}/${ARM_ID}.dpband.vcf.gz"
MISS_VCF="${VCF_DIR}/${ARM_ID}.maxmiss.vcf.gz"
NOPL_VCF="${VCF_DIR}/${ARM_ID}.filt.no_pl.vcf.gz"
OUT_VCF="${VCF_DIR}/${ARM_ID}.all_sites.vcf.gz"
GROUPS_TSV="${VCF_DIR}/${ARM_ID}.groups.tsv"

[[ -s "${RAW_VCF}"   ]] || die "Missing raw VCF: ${RAW_VCF} - rerun 0{3,5}_call (raw kept?)"
[[ -s "${REF_FASTA}" ]] || die "Missing reference: ${REF_FASTA}"
[[ -s "${MASK_BED}"  ]] || die "Missing mappability mask: ${MASK_BED} - run 05b_mappability.sbatch"
[[ -s "${POPS_SRC}"  ]] || die "Missing pops file: ${POPS_SRC}"
[[ -s "${RAW_VCF}.tbi" ]] || bcftools index -t "${RAW_VCF}"
mkdir -p "${TMP_DIR}"

# max-missing threshold: F_MISSING per pop must be <= (n - MIN_PRESENT)/n
MAXMISS=$(gawk -v n="${N_PER_POP}" -v k="${MIN_PRESENT}" 'BEGIN{printf "%.6f", (n-k)/n}')
# 2-col sample<TAB>pop for +fill-tags; distinct pop names for the filter expr
awk 'BEGIN{FS=OFS="\t"} /^#/||/^sample_id/{next} {print $1, $2}' "${POPS_SRC}" > "${GROUPS_TSV}"
mapfile -t POPS < <(cut -f2 "${GROUPS_TSV}" | sort -u)
[[ "${#POPS[@]}" -ge 1 ]] || die "[${ARM_ID}] no populations parsed from ${POPS_SRC}"

##########
# 1. genotype depth floor -> missing (DP absent = left as-is)
##########
info "[${ARM_ID}] genotypes with DP<${MIN_DP} -> missing"
rm -f "${GT_VCF}"{,.tbi}
bcftools filter -S . -e "FMT/DP<${MIN_DP}" "${RAW_VCF}" -Oz -o "${GT_VCF}"
bcftools index -t "${GT_VCF}"

##########
# 2. drop indels + near-indel SNPs; keep SNP + invariant only
# --SnpGap needs indels present, so this runs before VariantFiltration; TYPE
# membership drops any record carrying an indel/mnp/other allele (incl. mixed)
##########
info "[${ARM_ID}] flag SNPs within ${SNP_GAP} bp of an indel, drop indels + flagged SNPs"
rm -f "${SNP_VCF}"{,.tbi}
bcftools filter --SnpGap "${SNP_GAP}" "${GT_VCF}" -Oz \
  | bcftools view --threads "${THREADS}" \
      -e 'TYPE="indel" || TYPE="mnp" || TYPE="other" || TYPE="bnd" || FILTER="SnpGap"' \
      -Oz -o "${SNP_VCF}"
bcftools index -t "${SNP_VCF}"

##########
# 3. gatk hard filters + ExcessHet (variant sites); invariants lack annots -> stay PASS
##########
info "[${ARM_ID}] VariantFiltration (QD/FS/MQ/MQRankSum/ReadPosRankSum/SOR/ExcessHet)"
rm -f "${FLAG_VCF}"{,.tbi}
gatk --java-options "-Xmx24g" VariantFiltration \
    -R "${REF_FASTA}" \
    -V "${SNP_VCF}" \
    --filter-name QD        --filter-expression "QD < ${QD_MIN}" \
    --filter-name FS        --filter-expression "FS > ${FS_MAX}" \
    --filter-name MQ        --filter-expression "MQ < ${MQ_MIN}" \
    --filter-name MQRS      --filter-expression "MQRankSum < ${MQRS_MIN}" \
    --filter-name RPRS      --filter-expression "ReadPosRankSum < ${RPRS_MIN}" \
    --filter-name SOR       --filter-expression "SOR > ${SOR_MAX}" \
    --filter-name ExcessHet --filter-expression "ExcessHet > ${EXCESSHET_MAX}" \
    -O "${FLAG_VCF}" \
    --tmp-dir "${TMP_DIR}"

# keep PASS + unfiltered/invariant (-f PASS,. never drops invariant sites)
info "[${ARM_ID}] drop hard-fail sites (keep PASS + invariant)"
rm -f "${SITEPASS_VCF}"{,.tbi}
bcftools view -f 'PASS,.' "${FLAG_VCF}" -Oz -o "${SITEPASS_VCF}"
bcftools index -t "${SITEPASS_VCF}"

##########
# 4. mappability mask (symmetric over variant + invariant)
##########
info "[${ARM_ID}] apply mappability mask: ${MASK_BED}"
rm -f "${MASK_VCF}"{,.tbi}
bcftools view --threads "${THREADS}" -T "${MASK_BED}" "${SITEPASS_VCF}" -Oz -o "${MASK_VCF}"
bcftools index -t "${MASK_VCF}"

##########
# 5. site depth band: INFO/DP quantiles (one-pass integer histogram)
##########
info "[${ARM_ID}] computing INFO/DP quantiles (qlo=${DP_QLO}, qhi=${DP_QHI})"
read -r DP_LO DP_HI < <(
    bcftools query -f '%INFO/DP\n' "${MASK_VCF}" \
      | gawk -v qlo="${DP_QLO}" -v qhi="${DP_QHI}" '
          $1!="." && $1!="" { c[$1]++; n++ }
          END{
              if(n==0){ print "NA NA"; exit }
              m=asorti(c, keys, "@ind_num_asc")
              cum=0; lo=""; hi=""
              for(i=1;i<=m;i++){
                  cum+=c[keys[i]]
                  if(lo=="" && cum>=qlo*n) lo=keys[i]
                  if(hi=="" && cum>=qhi*n) hi=keys[i]
              }
              printf "%s %s\n", lo, hi
          }'
)
[[ "${DP_LO}" != "NA" && -n "${DP_LO}" && -n "${DP_HI}" ]] || die "[${ARM_ID}] failed to compute DP quantiles"
info "[${ARM_ID}] depth band: keep ${DP_LO} <= INFO/DP <= ${DP_HI}"
rm -f "${DPBAND_VCF}"{,.tbi}
bcftools view -e "INFO/DP<${DP_LO} || INFO/DP>${DP_HI}" "${MASK_VCF}" -Oz -o "${DPBAND_VCF}"
bcftools index -t "${DPBAND_VCF}"

##########
# 6. max-missing: each pop >= MIN_PRESENT called genotypes (per-pop F_MISSING)
##########
info "[${ARM_ID}] max-missing: each pop >= ${MIN_PRESENT}/${N_PER_POP} present (F_MISSING <= ${MAXMISS})"
TAGGED_VCF="${VCF_DIR}/${ARM_ID}.ftags.vcf.gz"
rm -f "${TAGGED_VCF}"{,.tbi}
bcftools +fill-tags "${DPBAND_VCF}" -Oz -o "${TAGGED_VCF}" -- -S "${GROUPS_TSV}" -t 'F_MISSING'
bcftools index -t "${TAGGED_VCF}"
# exclude a site if any population exceeds the missing threshold
MISS_EXPR=""
for p in "${POPS[@]}"; do
    [[ -n "${MISS_EXPR}" ]] && MISS_EXPR+=" || "
    MISS_EXPR+="F_MISSING_${p}>${MAXMISS}"
done
info "[${ARM_ID}] max-missing exclude expr: ${MISS_EXPR}"
rm -f "${MISS_VCF}"{,.tbi}
bcftools view -e "${MISS_EXPR}" "${TAGGED_VCF}" -Oz -o "${MISS_VCF}"
bcftools index -t "${MISS_VCF}"

##########
# 7. strip FORMAT annots + trim unused alt alleles (mirrors 03/05)
##########
info "[${ARM_ID}] strip FORMAT annots + trim-alt-alleles"
rm -f "${NOPL_VCF}"{,.tbi}
bcftools annotate \
    -x 'FORMAT/PL,FORMAT/AD,FORMAT/GQ,FORMAT/PGT,FORMAT/PID,FORMAT/PS,FORMAT/MIN_DP,FORMAT/RGQ,FORMAT/SB' \
    "${MISS_VCF}" -Oz -o "${NOPL_VCF}"
bcftools index -t "${NOPL_VCF}"

rm -f "${OUT_VCF}"{,.tbi}
bcftools view --threads "${THREADS}" --trim-alt-alleles "${NOPL_VCF}" -Oz -o "${OUT_VCF}"
bcftools index -t "${OUT_VCF}"

# drop intermediates, keep the filtered all-sites vcf
rm -f "${GT_VCF}"{,.tbi} "${SNP_VCF}"{,.tbi} "${FLAG_VCF}"{,.tbi} "${SITEPASS_VCF}"{,.tbi} \
      "${MASK_VCF}"{,.tbi} "${DPBAND_VCF}"{,.tbi} "${TAGGED_VCF}"{,.tbi} \
      "${MISS_VCF}"{,.tbi} "${NOPL_VCF}"{,.tbi} "${GROUPS_TSV}"

n_raw=$(bcftools index -n "${RAW_VCF}")
n_out=$(bcftools index -n "${OUT_VCF}")
info "[${ARM_ID}] records raw=${n_raw} filtered=${n_out} (DP band ${DP_LO}-${DP_HI}, maxmiss ${MAXMISS})"
info "[${ARM_ID}] canonical filtered all-sites VCF: ${OUT_VCF}"
