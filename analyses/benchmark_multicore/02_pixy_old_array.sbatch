#!/usr/bin/env bash
# run pixy 0.95.01 (zarr, single-core) on each 10 Mb vcf for one stat
# STAT set via --export: pi, dxy, or fst. one seed per array task. window 25 kb
# zarr build is timed on purpose (part of the 0.95.01 workflow)
# if --bypass_filtration is absent in 0.95.01, drop it from the pixy call below
# env: old_pixy (envs/pixy_old.yml)
#SBATCH --job-name=pixy_old
#SBATCH --output=logs/pixy_old_%x_%A_%a.out
#SBATCH --error=logs/pixy_old_%x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
# re-tune from seff after first run
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

STAT="${STAT:-pi}"

case "${STAT}" in
  pi)       POP_MAP="${PROJECT_DIR}/config/pop_map_1pop.tsv" ;;
  dxy|fst)  POP_MAP="${PROJECT_DIR}/config/pop_map_2pop.tsv" ;;
  *) echo "Unknown STAT '${STAT}'; must be pi, dxy, or fst" >&2; exit 1 ;;
esac

case "${STAT}" in
  pi)  STAT_COL="Pi";  OUT_FILE="pixy_output_pi.txt";  VAL_COL="avg_pi"     ;;
  dxy) STAT_COL="Dxy"; OUT_FILE="pixy_output_dxy.txt"; VAL_COL="avg_dxy"    ;;
  fst) STAT_COL="Fst"; OUT_FILE="pixy_output_fst.txt"; VAL_COL="avg_wc_fst" ;;
esac

WORK_BASE="${PROJECT_DIR}/data/pixy_old_out/10Mb_${STAT}"
ZARR_BASE="${PROJECT_DIR}/data/zarr/10Mb_${STAT}"
SUMMARY="${PROJECT_DIR}/data/results/pixy_old_${STAT}_10Mb.tsv"
LOCKFILE="${SUMMARY}.lock"

##########
# env
##########
if command -v conda >/dev/null 2>&1; then
  eval "$(conda shell.bash hook)"
  conda activate old_pixy
fi
command -v pixy          >/dev/null 2>&1 || { echo "pixy not on PATH" >&2; exit 1; }
command -v /usr/bin/time >/dev/null 2>&1 || { echo "/usr/bin/time not on PATH" >&2; exit 1; }
command -v flock         >/dev/null 2>&1 || { echo "flock not on PATH" >&2; exit 1; }

mkdir -p "${WORK_BASE}" "${ZARR_BASE}" "${LOG_DIR}" "$(dirname "${SUMMARY}")"
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
    printf "%s\t0.95.01\tNA\t0\t0\t0\t0\t0\tMISSING_VCF\n" "${seed}" >> "${SUMMARY}"
  exit 0
fi

RUN_DIR="${WORK_BASE}/seed_${seed}"; mkdir -p "${RUN_DIR}"
ZARR_SEED="${ZARR_BASE}/seed_${seed}"; mkdir -p "${ZARR_SEED}"
out_dst="${RUN_DIR}/${OUT_FILE}"
timing="${RUN_DIR}/timing.txt"
run_log="${RUN_DIR}/pixy_run.log"

# 0.95.01 writes to cwd; run from per-seed dir
pushd "${RUN_DIR}" >/dev/null
ln -sf "${vcf}" ./input.vcf.gz
[[ -f "${vcf}.tbi" ]] && ln -sf "${vcf}.tbi" ./input.vcf.gz.tbi

status="OK"; val="NA"; rc=0
set +e
/usr/bin/time -f "%e %U %S %M" -o "${timing}" \
  pixy --stats "${STAT}" \
       --vcf "./input.vcf.gz" \
       --populations "${POP_MAP}" \
       --window_size "${WINDOW}" \
       --zarr_path "${ZARR_SEED}" \
       --bypass_filtration yes \
  > "${run_log}" 2>&1
rc=$?
set -e

popd >/dev/null

elapsed_s=0; user_s=0; sys_s=0; elapsed_us=0
if [[ -s "${timing}" ]]; then
  IFS=' ' read -r elapsed_s user_s sys_s rss_kb < "${timing}" || true
  elapsed_us=$(awk -v s="${elapsed_s}" 'BEGIN{printf "%.0f", s*1000000}')
fi

if [[ -s "${out_dst}" ]]; then
  val=$(awk -F'\t' -v col="${VAL_COL}" '
    NR==1 { for(i=1;i<=NF;i++) if($i==col||(col=="avg_wc_fst"&&($i=="wc_fst"||$i=="fst"))) c=i }
    NR>1 && c { s+=$c; n++ }
    END { if(n>0) printf "%.15f", s/n; else print "NA" }' "${out_dst}")
fi

if   [[ ${rc} -ne 0 ]];       then status="PIXY_FAIL_${rc}"
elif [[ ! -s "${out_dst}" ]];  then status="NO_OUTPUT"
elif [[ -z "${val}" || "${val}" == "NA" ]]; then status="NO_VALUE"
fi

# cores col = "0.95.01" to distinguish from pixy_new 1-core
flock -x "${LOCKFILE}" \
  printf "%s\t0.95.01\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${seed}" "${val:-NA}" \
    "${elapsed_us}" "${elapsed_s}" "${user_s}" "${sys_s}" "${rss_kb:-0}" "${status}" >> "${SUMMARY}"

echo "Seed ${seed} | stat=${STAT} | val=${val} | status=${status} | elapsed=${elapsed_s}s"

# clean up per-seed dir and zarr
rm -rf "${RUN_DIR}" "${ZARR_SEED}"
