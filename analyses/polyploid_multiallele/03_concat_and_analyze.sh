#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name="pixy_concat"
#SBATCH -p batch
#SBATCH --output=logs/%x.%j.out

# final concat step. discovers arms from config/sim_params.tsv (or env var ARMS,
# space-separated, set by 01_run_all.sh). per arm:
#   1. read the per-task header task 1 emitted to data/${ARM}.part_0001.tsv.header
#   2. concat every data/${ARM}.part_NNNN.tsv into data/aggregated/${ARM}.tsv,
#      header prepended
#   3. report row counts
#
#   sbatch 03_concat_and_analyze.sh                   # all arms in sim_params.tsv
#   ARMS="pi_tet_1pop_miss25" sbatch 03_concat_and_analyze.sh
#
# or invoked from 01_run_all.sh, which sets ARMS to the wave it submitted and
# adds --dependency=afterok on the array jobs.

set -euo pipefail

PARAMS_TSV="${PARAMS_TSV:-config/sim_params.tsv}"

# if ARMS unset, discover from sim_params.tsv, skipping legacy rows
if [[ -z "${ARMS:-}" ]]; then
    ARMS=$(awk -F'\t' '
        NR==1 {
            for (i=1; i<=NF; i++) col[$i]=i
            next
        }
        {
            notes = $(col["notes"])
            if (notes ~ /do not rerun/) next
            print $(col["arm"])
        }
    ' "$PARAMS_TSV")
fi

mkdir -p data/aggregated

echo "[$(date '+%F %T')] concat: building per-arm TSVs for: $ARMS"

for ARM in $ARMS; do
    RESULTS="data/aggregated/${ARM}.tsv"
    HEADER_FILE="data/${ARM}.part_0001.tsv.header"

    shopt -s nullglob
    PARTS=( data/${ARM}.part_*.tsv )
    shopt -u nullglob

    if [[ ${#PARTS[@]} -eq 0 ]]; then
        echo "[concat] WARNING: no partials for ${ARM} (expected data/${ARM}.part_*.tsv); skipping" >&2
        continue
    fi

    if [[ ! -f "$HEADER_FILE" ]]; then
        echo "[concat] WARNING: no header for ${ARM} (expected $HEADER_FILE); skipping" >&2
        continue
    fi

    cat "$HEADER_FILE" > "$RESULTS"
    # sort partials so concat order matches array task order
    printf '%s\n' "${PARTS[@]}" | sort | xargs cat >> "$RESULTS"

    n_rows=$(( $(wc -l < "$RESULTS") - 1 ))
    echo "[concat] ${ARM}: ${#PARTS[@]} partials -> ${RESULTS} (${n_rows} rows)"
done

echo
echo "[$(date '+%F %T')] concat complete. Aggregated TSVs in data/aggregated/"
