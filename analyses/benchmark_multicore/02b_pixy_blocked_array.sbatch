#!/usr/bin/env bash
# blocked timing benchmark -- supersedes 02_pixy_old_array.sbatch and
# 03_pixy_new_array.sbatch for Figure 2a. those two remain valid but collect
# each (version x cores) cell in a separate array at a separate time, so the
# panel-a ratio median(0.95.01) / median(2.2.3 @ n) divides two numbers measured
# under different cluster load. that is unfixable by scheduling: /bigdata is
# shared and other users cannot be kept off it.
#
# instead, one array task = one (stat, seed) BLOCK, running all six measurements
# back to back on one node:
#
#     0.95.01 (1 core) -> 2.2.3 @ 1 -> @2 -> @4 -> @8 -> @16
#
# so every ratio is formed from runs taken minutes apart under the same ambient
# load, and ambient load largely cancels. two further protections:
#
#   1. the vcf is copied to node-local scratch and every run reads it from
#      there, with all pixy scratch (zarr, outputs) local too. only the summary
#      append touches /bigdata. the shared filesystem is out of the timed path.
#   2. run order is rotated by seed so no run type sits permanently first, where
#      it would pay the cold-cache cost for the rest of the block.
#
# output is written to the SAME per-cell summary TSVs as 02/03, so
# 04_aggregate_summaries.sh and Figure_2_Performance.R need no changes. the
# pairing is recoverable by joining cells on Seed.
#
# one array of 300 tasks = 3 stats x 100 seeds, interleaved so the stats stay
# mixed in time. resubmit failures with an explicit --array list.
# env: old_pixy (envs/pixy_old.yml), new_pixy (envs/pixy_new.yml)
#SBATCH --job-name=pixy_block
#SBATCH --output=logs/pixy_block_%A_%a.out
#SBATCH --error=logs/pixy_block_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --partition=intel
# worst observed block is ~16 min (old 280s + new 236/248/86/43/39s)
#SBATCH --time=0-01:30:00
# 24 x 16 cpus = 384 = the samuklab group cap
#SBATCH --array=1-300%24

set -euo pipefail
IFS=$'\n\t'

##########
# parameters
##########
PROJECT_DIR="${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}}"
VCF_DIR="${PROJECT_DIR}/data/vcfs/dm_10Mb"
LOG_DIR="${PROJECT_DIR}/logs"
RESULTS_DIR="${PROJECT_DIR}/data/results"
# node-local disk, keeping /bigdata out of the timed path. Use /tmp: it is
# node-local, has a sticky bit (safe in a shared dir), and is writable by the job
# on every compute node. NOT /scratch/$USER -- that per-user dir exists root-owned
# on compute nodes so a subdir create is denied there, and /scratch top-level has
# no sticky bit. The node probe showed /tmp at 91 GB free / 206 MB/s, ample for a
# 42 MB VCF copy plus zarr, which each task cleans up on exit.
LOCAL_BASE="${LOCAL_BASE:-/tmp}"
WINDOW=25000
CHUNK=25000

CORES=(1 2 4 8 16)
STATS=(pi dxy fst)

# task id -> (stat, seed), interleaved: 1=pi/1, 2=dxy/1, 3=fst/1, 4=pi/2, ...
task="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID not set (submit via sbatch with --array)}"
STAT="${STATS[$(( (task - 1) % 3 ))]}"
seed=$(( (task - 1) / 3 + 1 ))

case "${STAT}" in
  pi)       POP_MAP="${PROJECT_DIR}/config/pop_map_1pop.tsv" ;;
  dxy|fst)  POP_MAP="${PROJECT_DIR}/config/pop_map_2pop.tsv" ;;
esac

case "${STAT}" in
  pi)  STAT_COL="Pi";  OLD_OUT="pixy_output_pi.txt";  VAL_COL="avg_pi"     ;;
  dxy) STAT_COL="Dxy"; OLD_OUT="pixy_output_dxy.txt"; VAL_COL="avg_dxy"    ;;
  fst) STAT_COL="Fst"; OLD_OUT="pixy_output_fst.txt"; VAL_COL="avg_wc_fst" ;;
esac

OLD_SUMMARY="${RESULTS_DIR}/pixy_old_${STAT}_10Mb.tsv"

new_summary() { echo "${RESULTS_DIR}/pixy_new_${STAT}_10Mb_cores_${1}.tsv"; }

##########
# env
##########
# ~/.bashrc references MODULESHOME with no default and aborts under set -u
activate_env() {
    set +u
    eval "$(conda shell.bash hook)"
    conda activate "$1"
    set -u
    command -v pixy >/dev/null 2>&1 || { echo "pixy not on PATH in $1" >&2; exit 1; }
}

command -v conda         >/dev/null 2>&1 || { echo "conda not on PATH" >&2; exit 1; }
command -v /usr/bin/time >/dev/null 2>&1 || { echo "/usr/bin/time not on PATH" >&2; exit 1; }
command -v flock         >/dev/null 2>&1 || { echo "flock not on PATH" >&2; exit 1; }

mkdir -p "${LOG_DIR}" "${RESULTS_DIR}"

# headers, once per summary file
init_header() {
  local f="$1"
  flock -x "${f}.lock" bash -c "
    [[ -f '${f}' ]] || printf 'Seed\tCores\t${STAT_COL}\tElapsed_us\tElapsed_s\tUser_s\tSys_s\tRSS_kb\tStatus\n' > '${f}'
  "
}
init_header "${OLD_SUMMARY}"
for n in "${CORES[@]}"; do init_header "$(new_summary "${n}")"; done

##########
# block scratch on node-local disk
##########
vcf="${VCF_DIR}/dm_sim_vcf_seed_${seed}.vcf.gz"

if [[ ! -f "${vcf}" ]]; then
  flock -x "${OLD_SUMMARY}.lock" \
    printf "%s\t0.95.01\tNA\t0\t0\t0\t0\t0\tMISSING_VCF\n" "${seed}" >> "${OLD_SUMMARY}"
  for n in "${CORES[@]}"; do
    s="$(new_summary "${n}")"
    flock -x "${s}.lock" \
      printf "%s\t%s\tNA\t0\t0\t0\t0\t0\tMISSING_VCF\n" "${seed}" "${n}" >> "${s}"
  done
  exit 0
fi

mkdir -p "${LOCAL_BASE}" || { echo "cannot create ${LOCAL_BASE} on ${HOSTNAME}" >&2; exit 1; }
BLOCK_DIR="${LOCAL_BASE}/pixy_block_${SLURM_JOB_ID:-$$}_${STAT}_${seed}"
# nothing purges /scratch (it holds files here from weeks ago), so clean up on the
# way out for signals too, not just a normal exit -- a walltime kill would
# otherwise strand the vcf copy and zarr store on the node
trap 'rm -rf "${BLOCK_DIR}"' EXIT INT TERM
mkdir -p "${BLOCK_DIR}"

cp "${vcf}" "${BLOCK_DIR}/input.vcf.gz"
[[ -f "${vcf}.tbi" ]] && cp "${vcf}.tbi" "${BLOCK_DIR}/input.vcf.gz.tbi"
LOCAL_VCF="${BLOCK_DIR}/input.vcf.gz"
# equalise page-cache state across the six runs
cat "${LOCAL_VCF}" > /dev/null

##########
# one measurement
##########
extract_val() {
  awk -F'\t' -v col="${VAL_COL}" '
    NR==1 { for(i=1;i<=NF;i++) if($i==col||(col=="avg_wc_fst"&&($i=="wc_fst"||$i=="fst"))) c=i }
    NR>1 && c { s+=$c; n++ }
    END { if(n>0) printf "%.15f", s/n; else print "NA" }' "$1"
}

# $1=timing file -> sets elapsed_s user_s sys_s rss_kb elapsed_us
read_timing() {
  elapsed_s=0; user_s=0; sys_s=0; rss_kb=0; elapsed_us=0
  if [[ -s "$1" ]]; then
    IFS=' ' read -r elapsed_s user_s sys_s rss_kb < "$1" || true
    elapsed_us=$(awk -v s="${elapsed_s}" 'BEGIN{printf "%.0f", s*1000000}')
  fi
}

run_old() {
  local run_dir="${BLOCK_DIR}/old" zarr_dir="${BLOCK_DIR}/zarr_old"
  mkdir -p "${run_dir}" "${zarr_dir}"

  activate_env old_pixy

  # 0.95.01 writes to cwd
  pushd "${run_dir}" >/dev/null
  ln -sf "${LOCAL_VCF}" ./input.vcf.gz
  [[ -f "${LOCAL_VCF}.tbi" ]] && ln -sf "${LOCAL_VCF}.tbi" ./input.vcf.gz.tbi

  local rc=0
  set +e
  /usr/bin/time -f "%e %U %S %M" -o "${run_dir}/timing.txt" \
    pixy --stats "${STAT}" \
         --vcf "./input.vcf.gz" \
         --populations "${POP_MAP}" \
         --window_size "${WINDOW}" \
         --zarr_path "${zarr_dir}" \
         --bypass_filtration yes \
    > "${run_dir}/pixy_run.log" 2>&1
  rc=$?
  set -e
  popd >/dev/null

  read_timing "${run_dir}/timing.txt"

  local val="NA" status="OK"
  [[ -s "${run_dir}/${OLD_OUT}" ]] && val=$(extract_val "${run_dir}/${OLD_OUT}")
  if   [[ ${rc} -ne 0 ]];                     then status="PIXY_FAIL_${rc}"
  elif [[ ! -s "${run_dir}/${OLD_OUT}" ]];    then status="NO_OUTPUT"
  elif [[ -z "${val}" || "${val}" == "NA" ]]; then status="NO_VALUE"
  fi

  flock -x "${OLD_SUMMARY}.lock" \
    printf "%s\t0.95.01\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "${seed}" "${val:-NA}" "${elapsed_us}" "${elapsed_s}" \
      "${user_s}" "${sys_s}" "${rss_kb:-0}" "${status}" >> "${OLD_SUMMARY}"

  echo "  [block ${STAT}/${seed}] 0.95.01        val=${val} elapsed=${elapsed_s}s status=${status}"
  rm -rf "${run_dir}" "${zarr_dir}"
}

run_new() {
  local n="$1"
  local run_dir="${BLOCK_DIR}/new_${n}" prefix="pixy_new_${STAT}_${seed}"
  mkdir -p "${run_dir}"

  activate_env new_pixy

  local rc=0
  set +e
  /usr/bin/time -f "%e %U %S %M" -o "${run_dir}/timing.txt" \
    pixy --stats "${STAT}" \
         --vcf "${LOCAL_VCF}" \
         --output_folder "${run_dir}" \
         --output_prefix "${prefix}" \
         --populations "${POP_MAP}" \
         --window_size "${WINDOW}" \
         --chunk_size "${CHUNK}" \
         --n_cores "${n}" \
    > "${run_dir}/pixy_run.log" 2>&1
  rc=$?
  set -e

  read_timing "${run_dir}/timing.txt"

  local out_file="${run_dir}/${prefix}_${STAT}.txt"
  local val="NA" status="OK"
  [[ -s "${out_file}" ]] && val=$(extract_val "${out_file}")
  if   [[ ${rc} -ne 0 ]];                     then status="PIXY_FAIL_${rc}"
  elif [[ ! -s "${out_file}" ]];              then status="NO_OUTPUT"
  elif [[ -z "${val}" || "${val}" == "NA" ]]; then status="NO_VALUE"
  fi

  local s; s="$(new_summary "${n}")"
  flock -x "${s}.lock" \
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "${seed}" "${n}" "${val:-NA}" "${elapsed_us}" "${elapsed_s}" \
      "${user_s}" "${sys_s}" "${rss_kb:-0}" "${status}" >> "${s}"

  echo "  [block ${STAT}/${seed}] 2.2.3 @${n}       val=${val} elapsed=${elapsed_s}s status=${status}"
  rm -rf "${run_dir}"
}

##########
# run the block, order rotated by seed
##########
ORDER=(old 1 2 4 8 16)
shift_by=$(( (seed - 1) % ${#ORDER[@]} ))
ROTATED=()
for ((i=0; i<${#ORDER[@]}; i++)); do
  ROTATED+=( "${ORDER[$(( (i + shift_by) % ${#ORDER[@]} ))]}" )
done

echo "[$(date '+%F %T')] block stat=${STAT} seed=${seed} on ${HOSTNAME}"
# IFS is \n\t here, so join the order explicitly for the log line
echo "  order: $(IFS=' '; echo "${ROTATED[*]}")"

for r in "${ROTATED[@]}"; do
  if [[ "${r}" == "old" ]]; then run_old; else run_new "${r}"; fi
done

echo "[$(date '+%F %T')] block stat=${STAT} seed=${seed} complete"
