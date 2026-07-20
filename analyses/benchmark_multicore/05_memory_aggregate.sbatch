#!/usr/bin/env bash
# aggregate (process-tree) peak memory for pixy 2.2.3 at one core count
# /usr/bin/time %M in 02/03 reports the largest single process, not the sum over
# concurrent workers. this samples the whole process tree so parent + workers are
# counted simultaneously. reports both quantities so they can be compared directly.
# STAT, N_CORES, RUN_TAG passed via --export. one seed per array task
# env: new_pixy (envs/pixy_new.yml)
#SBATCH --job-name=pixy_mem
#SBATCH --output=logs/pixy_mem_%x_%A_%a.out
#SBATCH --error=logs/pixy_mem_%x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --time=0-00:30:00
#SBATCH --array=1-10

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
SAMPLE_INTERVAL=0.05

STAT="${STAT:-pi}"
N_CORES="${N_CORES:-1}"
RUN_TAG="${RUN_TAG:-cores_${N_CORES}}"

case "${STAT}" in
  pi)       POP_MAP="${PROJECT_DIR}/config/pop_map_1pop.tsv" ;;
  dxy|fst)  POP_MAP="${PROJECT_DIR}/config/pop_map_2pop.tsv" ;;
  *) echo "Unknown STAT '${STAT}'; must be pi, dxy, or fst" >&2; exit 1 ;;
esac

WORK_DIR="${PROJECT_DIR}/data/pixy_mem_out/10Mb_${STAT}/${RUN_TAG}"
SUMMARY="${PROJECT_DIR}/data/results/pixy_mem_${STAT}_10Mb_${RUN_TAG}.tsv"
LOCKFILE="${SUMMARY}.lock"

##########
# env
##########
if command -v conda >/dev/null 2>&1; then
  eval "$(conda shell.bash hook)"
  conda activate new_pixy
fi
command -v pixy   >/dev/null 2>&1 || { echo "pixy not on PATH" >&2; exit 1; }
command -v setsid >/dev/null 2>&1 || { echo "setsid not on PATH" >&2; exit 1; }
command -v flock  >/dev/null 2>&1 || { echo "flock not on PATH" >&2; exit 1; }

mkdir -p "${WORK_DIR}" "${LOG_DIR}" "$(dirname "${SUMMARY}")"
flock -x "${LOCKFILE}" bash -c "
  [[ -f '${SUMMARY}' ]] || printf 'Seed\tCores\tPixy_version\tTree_peak_kb\tProc_peak_kb\tN_proc_at_peak\tElapsed_s\tStatus\n' > '${SUMMARY}'
"

##########
# per-task
##########
seed="${SLURM_ARRAY_TASK_ID}"
vcf="${VCF_DIR}/dm_sim_vcf_seed_${seed}.vcf.gz"
pixy_version="$(pixy --version 2>&1 | tr -d '\n')"

if [[ ! -f "${vcf}" ]]; then
  flock -x "${LOCKFILE}" \
    printf "%s\t%s\t%s\tNA\tNA\tNA\t0\tMISSING_VCF\n" "${seed}" "${N_CORES}" "${pixy_version}" >> "${SUMMARY}"
  exit 0
fi

prefix="pixy_mem_${STAT}_${seed}"
out_file="${WORK_DIR}/${prefix}_${STAT}.txt"
run_log="${WORK_DIR}/pixy_${seed}.log"

# own process group so the sampler can find every worker by pgid
start_s=$(date +%s.%N)
setsid pixy --stats "${STAT}" \
            --vcf "${vcf}" \
            --output_folder "${WORK_DIR}" \
            --output_prefix "${prefix}" \
            --populations "${POP_MAP}" \
            --window_size "${WINDOW}" \
            --chunk_size "${CHUNK}" \
            --n_cores "${N_CORES}" \
            > "${run_log}" 2>&1 &
pixy_pid=$!
pgid="${pixy_pid}"

tree_peak=0; proc_peak=0; n_at_peak=0
while kill -0 "${pixy_pid}" 2>/dev/null; do
  read -r sum maxrss n < <(
    ps -eo pgid=,rss= 2>/dev/null |
    awk -v g="${pgid}" '$1==g { s+=$2; if($2>m) m=$2; n++ } END { print s+0, m+0, n+0 }'
  )
  if (( sum > tree_peak )); then tree_peak=${sum}; n_at_peak=${n}; fi
  if (( maxrss > proc_peak )); then proc_peak=${maxrss}; fi
  sleep "${SAMPLE_INTERVAL}"
done

set +e
wait "${pixy_pid}"
rc=$?
set -e
elapsed_s=$(awk -v a="${start_s}" -v b="$(date +%s.%N)" 'BEGIN{printf "%.2f", b-a}')

status="OK"
if   [[ ${rc} -ne 0 ]];        then status="PIXY_FAIL_${rc}"
elif [[ ! -s "${out_file}" ]]; then status="NO_OUTPUT"
elif (( tree_peak == 0 ));     then status="NO_SAMPLES"
fi

flock -x "${LOCKFILE}" \
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${seed}" "${N_CORES}" "${pixy_version}" \
    "${tree_peak}" "${proc_peak}" "${n_at_peak}" "${elapsed_s}" "${status}" >> "${SUMMARY}"

echo "Seed ${seed} | stat=${STAT} cores=${N_CORES} | tree_peak=${tree_peak}kB proc_peak=${proc_peak}kB (${n_at_peak} procs) | status=${status}"

rm -f "${out_file}" "${run_log}"
