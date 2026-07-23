#!/usr/bin/env bash
# aggregate (process-tree) peak memory for the legacy pixy 0.95.01 baseline.
#
# the companion 05_memory_aggregate.sbatch does this for pixy 2.2.3. both are
# needed: comparing a process-tree peak for the new version against the
# `/usr/bin/time %M` figure recorded for the old one would compare a sum over
# concurrent workers to a single-process maximum, which is the very conflation
# this benchmark exists to remove. 0.95.01 is single-core, so its tree peak and
# process peak should agree closely -- but that has to be measured, not assumed,
# and the zarr build is part of the timed workflow.
#
# STAT passed via --export. one seed per array task.
# env: old_pixy (envs/pixy_old.yml -- pixy 0.95.01, self-reports "0.95.0")
#SBATCH --job-name=pixy_mem_old
#SBATCH --output=logs/pixy_mem_old_%x_%A_%a.out
#SBATCH --error=logs/pixy_mem_old_%x_%A_%a.err
#SBATCH --ntasks=1 --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=intel
#SBATCH --time=0-01:45:00
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
SAMPLE_INTERVAL=0.05

STAT="${STAT:-pi}"

case "${STAT}" in
  pi)       POP_MAP="${PROJECT_DIR}/config/pop_map_1pop.tsv" ;;
  dxy|fst)  POP_MAP="${PROJECT_DIR}/config/pop_map_2pop.tsv" ;;
  *) echo "Unknown STAT '${STAT}'; must be pi, dxy, or fst" >&2; exit 1 ;;
esac

case "${STAT}" in
  pi)  OUT_FILE="pixy_output_pi.txt"  ;;
  dxy) OUT_FILE="pixy_output_dxy.txt" ;;
  fst) OUT_FILE="pixy_output_fst.txt" ;;
esac

WORK_BASE="${PROJECT_DIR}/data/pixy_mem_old_out/10Mb_${STAT}"
ZARR_BASE="${PROJECT_DIR}/data/zarr_mem/10Mb_${STAT}"
SUMMARY="${PROJECT_DIR}/data/results/pixy_mem_old_${STAT}_10Mb.tsv"
LOCKFILE="${SUMMARY}.lock"

##########
# env
##########
if command -v conda >/dev/null 2>&1; then
  eval "$(conda shell.bash hook)"
  conda activate old_pixy
fi
command -v pixy    >/dev/null 2>&1 || { echo "pixy not on PATH" >&2; exit 1; }
command -v python3 >/dev/null 2>&1 || { echo "python3 not on PATH" >&2; exit 1; }
command -v flock   >/dev/null 2>&1 || { echo "flock not on PATH" >&2; exit 1; }
[[ -f "${PROJECT_DIR}/sample_tree_rss.py" ]] || { echo "missing sample_tree_rss.py" >&2; exit 1; }

mkdir -p "${WORK_BASE}" "${ZARR_BASE}" "${LOG_DIR}" "$(dirname "${SUMMARY}")"
flock -x "${LOCKFILE}" bash -c "
  [[ -f '${SUMMARY}' ]] || printf 'Seed\tCores\tPixy_version\tTree_peak_kb\tProc_peak_kb\tN_proc_at_peak\tN_samples\tElapsed_s\tStatus\n' > '${SUMMARY}'
"

##########
# per-task
##########
seed="${SLURM_ARRAY_TASK_ID}"
vcf="${VCF_DIR}/dm_sim_vcf_seed_${seed}.vcf.gz"
# 0.95.01 self-reports "0.95.0"; record the pinned package version instead
pixy_version="0.95.01"

if [[ ! -f "${vcf}" ]]; then
  flock -x "${LOCKFILE}" \
    printf "%s\t1\t%s\tNA\tNA\tNA\tNA\t0\tMISSING_VCF\n" "${seed}" "${pixy_version}" >> "${SUMMARY}"
  exit 0
fi

RUN_DIR="${WORK_BASE}/seed_${seed}"; mkdir -p "${RUN_DIR}"
ZARR_SEED="${ZARR_BASE}/seed_${seed}"; mkdir -p "${ZARR_SEED}"
sample_out="${RUN_DIR}/tree_rss.txt"
run_log="${RUN_DIR}/pixy_run.log"

# 0.95.01 writes to cwd; run from per-seed dir
pushd "${RUN_DIR}" >/dev/null
ln -sf "${vcf}" ./input.vcf.gz
[[ -f "${vcf}.tbi" ]] && ln -sf "${vcf}.tbi" ./input.vcf.gz.tbi

start_s=$(date +%s.%N)
pixy --stats "${STAT}" \
     --vcf "./input.vcf.gz" \
     --populations "${POP_MAP}" \
     --window_size "${WINDOW}" \
     --zarr_path "${ZARR_SEED}" \
     --bypass_filtration yes \
     > "${run_log}" 2>&1 &
pixy_pid=$!

python3 "${PROJECT_DIR}/sample_tree_rss.py" \
        "${pixy_pid}" "${SAMPLE_INTERVAL}" "${sample_out}" &
sampler_pid=$!

set +e
wait "${pixy_pid}"
rc=$?
set -e
wait "${sampler_pid}" 2>/dev/null || true
popd >/dev/null

elapsed_s=$(awk -v a="${start_s}" -v b="$(date +%s.%N)" 'BEGIN{printf "%.2f", b-a}')

tree_peak=0; proc_peak=0; n_at_peak=0; n_samples=0
if [[ -s "${sample_out}" ]]; then
  IFS=$'\t' read -r tree_peak proc_peak n_at_peak n_samples < "${sample_out}"
fi

status="OK"
if   [[ ${rc} -ne 0 ]];                     then status="PIXY_FAIL_${rc}"
elif [[ ! -s "${RUN_DIR}/${OUT_FILE}" ]];   then status="NO_OUTPUT"
elif (( n_samples == 0 ));                  then status="NO_SAMPLES"
fi

flock -x "${LOCKFILE}" \
  printf "%s\t1\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${seed}" "${pixy_version}" \
    "${tree_peak}" "${proc_peak}" "${n_at_peak}" "${n_samples}" "${elapsed_s}" "${status}" >> "${SUMMARY}"

echo "Seed ${seed} | stat=${STAT} | 0.95.01 | tree_peak=${tree_peak}kB proc_peak=${proc_peak}kB (${n_at_peak} procs, ${n_samples} samples) | status=${status}"

# zarr stores are large; keep the summary, drop the intermediates
rm -rf "${ZARR_SEED}" "${RUN_DIR}"
