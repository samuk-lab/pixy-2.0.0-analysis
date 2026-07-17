#!/usr/bin/env bash
# run pixy 2.0.0 on each 10 Mb vcf for one stat at one core count
# STAT, N_CORES, RUN_TAG passed via --export. one seed per array task
# window 25 kb = chunk 25 kb: 10 Mb / 25 kb = 400 chunks -> 25 chunks/core at
# 16 cores, near-linear scaling. prior benchmark used 100 kb (6 chunks/core), capping speedup
# env: new_pixy (envs/pixy_new.yml)
#SBATCH --job-name=pixy_new
#SBATCH --output=logs/pixy_new_%x_%A_%a.out
#SBATCH --error=logs/pixy_new_%x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --mem=8G
# overridden per-arm by 03_submit_pixy_new_all.sh; re-tune from seff
#SBATCH --time=0-00:30:00
#SBATCH --array=1-100%25

set -euo pipefail
IFS=$'\n\t'

##########
# parameters
##########
PROJECT_DIR="${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}}"
VCF_DIR="${PROJECT_DIR}/data/vcfs/dm_10Mb"
LOG_DIR="${PROJECT_DIR}/logs"
WINDOW=25000
CHUNK=25000

STAT="${STAT:-pi}"
N_CORES="${N_CORES:-1}"
RUN_TAG="${RUN_TAG:-cores_${N_CORES}}"

case "${STAT}" in
  pi)       POP_MAP="${PROJECT_DIR}/config/pop_map_1pop.tsv" ;;
  dxy|fst)  POP_MAP="${PROJECT_DIR}/config/pop_map_2pop.tsv" ;;
  *) echo "Unknown STAT '${STAT}'; must be pi, dxy, or fst" >&2; exit 1 ;;
esac

case "${STAT}" in
  pi)  STAT_COL="Pi";  VAL_COL="avg_pi"     ;;
  dxy) STAT_COL="Dxy"; VAL_COL="avg_dxy"    ;;
  fst) STAT_COL="Fst"; VAL_COL="avg_wc_fst" ;;
esac

WORK_DIR="${PROJECT_DIR}/data/pixy_new_out/10Mb_${STAT}/${RUN_TAG}"
SUMMARY="${PROJECT_DIR}/data/results/pixy_new_${STAT}_10Mb_${RUN_TAG}.tsv"
LOCKFILE="${SUMMARY}.lock"

##########
# env
##########
if command -v conda >/dev/null 2>&1; then
  eval "$(conda shell.bash hook)"
  conda activate new_pixy
fi
command -v pixy          >/dev/null 2>&1 || { echo "pixy not on PATH" >&2; exit 1; }
command -v /usr/bin/time >/dev/null 2>&1 || { echo "/usr/bin/time not on PATH" >&2; exit 1; }
command -v flock         >/dev/null 2>&1 || { echo "flock not on PATH" >&2; exit 1; }

mkdir -p "${WORK_DIR}" "${LOG_DIR}" "$(dirname "${SUMMARY}")"
flock -x "${LOCKFILE}" bash -c "
  [[ -f '${SUMMARY}' ]] || printf 'Seed\tCores\t${STAT_COL}\tElapsed_us\tElapsed_s\tUser_s\tSys_s\tRSS_kb\tStatus\n' > '${SUMMARY}'
"

##########
# per-task
##########
seed="${SLURM_ARRAY_TASK_ID}"
vcf="${VCF_DIR}/dm_sim_vcf_seed_${seed}.vcf.gz"

if [[ ! -f "${vcf}" ]]; then
  flock -x "${LOCKFILE}" \
    printf "%s\t%s\tNA\t0\t0\t0\t0\t0\tMISSING_VCF\n" "${seed}" "${N_CORES}" >> "${SUMMARY}"
  exit 0
fi

prefix="pixy_new_${STAT}_${seed}"
out_file="${WORK_DIR}/${prefix}_${STAT}.txt"
timing="${WORK_DIR}/timing_${seed}.txt"
run_log="${WORK_DIR}/pixy_${seed}.log"

status="OK"; val="NA"; rc=0
set +e
/usr/bin/time -f "%e %U %S %M" -o "${timing}" \
  pixy --stats "${STAT}" \
       --vcf "${vcf}" \
       --output_folder "${WORK_DIR}" \
       --output_prefix "${prefix}" \
       --populations "${POP_MAP}" \
       --window_size "${WINDOW}" \
       --chunk_size "${CHUNK}" \
       --n_cores "${N_CORES}" \
  > "${run_log}" 2>&1
rc=$?
set -e

elapsed_s=0; user_s=0; sys_s=0; elapsed_us=0
if [[ -s "${timing}" ]]; then
  IFS=' ' read -r elapsed_s user_s sys_s rss_kb < "${timing}" || true
  elapsed_us=$(awk -v s="${elapsed_s}" 'BEGIN{printf "%.0f", s*1000000}')
fi

if [[ -s "${out_file}" ]]; then
  val=$(awk -F'\t' -v col="${VAL_COL}" '
    NR==1 { for(i=1;i<=NF;i++) if($i==col||(col=="avg_wc_fst"&&($i=="wc_fst"||$i=="fst"))) c=i }
    NR>1 && c { s+=$c; n++ }
    END { if(n>0) printf "%.15f", s/n; else print "NA" }' "${out_file}")
fi

if   [[ ${rc} -ne 0 ]];      then status="PIXY_FAIL_${rc}"
elif [[ ! -s "${out_file}" ]]; then status="NO_OUTPUT"
elif [[ -z "${val}" || "${val}" == "NA" ]]; then status="NO_VALUE"
fi

flock -x "${LOCKFILE}" \
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${seed}" "${N_CORES}" "${val:-NA}" \
    "${elapsed_us}" "${elapsed_s}" "${user_s}" "${sys_s}" "${rss_kb:-0}" "${status}" >> "${SUMMARY}"

echo "Seed ${seed} | stat=${STAT} cores=${N_CORES} | val=${val} | status=${status} | elapsed=${elapsed_s}s"

rm -f "${out_file}" "${timing}" "${run_log}"
