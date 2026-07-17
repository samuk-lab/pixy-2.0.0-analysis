#!/bin/bash
# submit the polyploid / multiallelic expansion to slurm
# submits one 02_run_arm.sh array per arm in config/sim_params.tsv
# (filtered by --include / --exclude), then a single
# 03_concat_and_analyze.sh job with afterok on all of them.
#
# usage:
#
#   ./01_run_all.sh                              # all non-legacy arms
#   ./01_run_all.sh --include miss               # regex match on arm name
#   ./01_run_all.sh --include 'hex|oct|unequal|mixed|theta'
#   ./01_run_all.sh --exclude theta              # everything except that regex
#   ./01_run_all.sh --include-legacy             # also rerun the legacy 10k-rep arms
#   ./01_run_all.sh --dry-run                    # print sbatch lines, do not submit
#   ./01_run_all.sh --clean                      # wipe old partials for selected arms
#   ./01_run_all.sh -h | --help                  # show this help
#
# each arm is submitted as
#
#   sbatch --job-name="run_${ARM}" \
#          --array=1-${N_ARRAY_TASKS}%${CONCURRENT} \
#          --export=ALL,ARM_NAME=${ARM},N_ARRAY_TASKS=${N_ARRAY_TASKS} \
#          02_run_arm.sh
#
# N_ARRAY_TASKS defaults to 100, CONCURRENT to 50; override from the env:
#
#   N_ARRAY_TASKS=200 CONCURRENT=100 ./01_run_all.sh --include miss
#
#
#
# PARTITION / TIMELIMIT override the sbatch defaults for the arm arrays and the
# concat job. DEPENDENCY gates the arm arrays only (concat keeps its afterok):
#
#   PARTITION=short TIMELIMIT=01:45:00 DEPENDENCY=afterok:12345 \
#       ./01_run_all.sh --include theta
#
#
#

set -euo pipefail

usage() {
    sed -n '3,30p' "$0" | sed 's/^# \{0,1\}//'
}

INCLUDE=""
EXCLUDE=""
INCLUDE_LEGACY=0
DRY_RUN=0
CLEAN=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --include)        INCLUDE="$2"; shift 2 ;;
        --exclude)        EXCLUDE="$2"; shift 2 ;;
        --include-legacy) INCLUDE_LEGACY=1; shift ;;
        --dry-run)        DRY_RUN=1; shift ;;
        --clean)          CLEAN=1; shift ;;
        -h|--help)        usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

cd "$(dirname "$(readlink -f "$0")")"

PARAMS_TSV="${PARAMS_TSV:-config/sim_params.tsv}"
N_ARRAY_TASKS="${N_ARRAY_TASKS:-100}"
CONCURRENT="${CONCURRENT:-50}"

# optional partition / walltime override (e.g. PARTITION=short TIMELIMIT=01:45:00).
# applied to both the arm arrays and the concat job; empty = script defaults.
PARTITION="${PARTITION:-}"
TIMELIMIT="${TIMELIMIT:-}"
SBATCH_EXTRA=()
[[ -n "$PARTITION" ]] && SBATCH_EXTRA+=(--partition "$PARTITION")
[[ -n "$TIMELIMIT" ]] && SBATCH_EXTRA+=(--time "$TIMELIMIT")

# optional gate for the arm arrays (e.g. DEPENDENCY=afterok:12345 to hold them
# behind an env-setup job). sbatch does NOT read a SBATCH_DEPENDENCY env var, so
# pass it as an explicit flag.
# kept out of SBATCH_EXTRA: that array also goes to the concat job, and a second
# --dependency would override concat's afterok on the arms, letting it aggregate
# before they finish.
DEPENDENCY="${DEPENDENCY:-}"
ARM_SBATCH_EXTRA=()
[[ -n "$DEPENDENCY" ]] && ARM_SBATCH_EXTRA+=(--dependency "$DEPENDENCY")

if (( ! DRY_RUN )) && ! command -v sbatch >/dev/null 2>&1; then
    echo "[run_all] ERROR: sbatch not on PATH. Use --dry-run to inspect what would be submitted." >&2
    exit 1
fi

mkdir -p logs data

# discover arms via awk: one arm name per line, filtered by include/exclude
# regex and the legacy flag
ARMS=$(awk -F'\t' -v inc_re="$INCLUDE" -v exc_re="$EXCLUDE" -v keep_legacy="$INCLUDE_LEGACY" '
    NR==1 {
        for (i=1; i<=NF; i++) col[$i]=i
        next
    }
    {
        arm  = $(col["arm"])
        notes = $(col["notes"])
        if (keep_legacy == 0 && notes ~ /do not rerun/) next
        if (inc_re != "" && arm !~ inc_re) next
        if (exc_re != "" && arm ~  exc_re) next
        print arm
    }
' "$PARAMS_TSV")

if [[ -z "$ARMS" ]]; then
    echo "[run_all] No arms matched the filter. Inspect $PARAMS_TSV or your --include / --exclude." >&2
    exit 1
fi

N_ARMS=$(echo "$ARMS" | wc -l | tr -d ' ')
echo "[run_all] Selected ${N_ARMS} arm(s):"
echo "$ARMS" | sed 's/^/    /'
echo "  N_ARRAY_TASKS=${N_ARRAY_TASKS}  CONCURRENT=${CONCURRENT}"
echo

if (( CLEAN )); then
    for ARM in $ARMS; do
        rm -f  "data/${ARM}.part_"*.tsv "data/${ARM}.part_0001.tsv.header"
        rm -rf "data/_work/${ARM}"
        rm -f  "data/aggregated/${ARM}.tsv"
    done
    echo "[run_all] --clean: wiped partials for selected arms"
fi

ARM_JOB_IDS=()
for ARM in $ARMS; do
    if (( DRY_RUN )); then
        echo "  [DRY] sbatch --job-name=run_${ARM} --array=1-${N_ARRAY_TASKS}%${CONCURRENT} --export=ALL,ARM_NAME=${ARM},N_ARRAY_TASKS=${N_ARRAY_TASKS} ${SBATCH_EXTRA[*]-} ${ARM_SBATCH_EXTRA[*]-} 02_run_arm.sh"
        continue
    fi
    JOB_ID=$(sbatch --parsable \
        --job-name="run_${ARM}" \
        --array=1-"${N_ARRAY_TASKS}"%"${CONCURRENT}" \
        --export=ALL,ARM_NAME="${ARM}",N_ARRAY_TASKS="${N_ARRAY_TASKS}" \
        --chdir "$(pwd)" \
        ${SBATCH_EXTRA[@]+"${SBATCH_EXTRA[@]}"} \
        ${ARM_SBATCH_EXTRA[@]+"${ARM_SBATCH_EXTRA[@]}"} \
        02_run_arm.sh)
    ARM_JOB_IDS+=("$JOB_ID")
    printf "  %-40s array %s\n" "$ARM" "$JOB_ID"
done

if (( DRY_RUN )); then
    echo "[run_all] dry run: 03_concat_and_analyze.sh would follow with afterok on the above"
    exit 0
fi

# build dependency list for concat
DEP=$(IFS=':'; echo "${ARM_JOB_IDS[*]}")

# pass the same ARMS list to concat so it doesn't rediscover
CONCAT_JOB=$(sbatch --parsable \
    --dependency=afterok:"$DEP" \
    --export=ALL,ARMS="$(echo "$ARMS" | tr '\n' ' ')" \
    --chdir "$(pwd)" \
    ${SBATCH_EXTRA[@]+"${SBATCH_EXTRA[@]}"} \
    03_concat_and_analyze.sh)
echo
echo "[run_all] concat + aggregate  ${CONCAT_JOB}  (afterok on ${#ARM_JOB_IDS[@]} arm arrays)"
echo
echo "[run_all] Track progress with:"
echo "    squeue -u \$USER"
echo "    sacct -j $(IFS=,; echo "${ARM_JOB_IDS[*]}"),${CONCAT_JOB}"
echo
echo "[run_all] After ${CONCAT_JOB} completes, aggregated TSVs live in data/aggregated/."
