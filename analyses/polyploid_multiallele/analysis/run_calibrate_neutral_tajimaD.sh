#!/usr/bin/env bash
# regenerate the neutral tajima D reference: run the 4 ploidies of
# calibrate_neutral_tajimaD.py in parallel, merge partials into
# neutral_tajimaD_reference.tsv. serial mode (no arg) gives an identical file.
#
# env:
#   PYTHON      python with msprime + scikit-allel + pixy importable
#               (default python3; set PYTHONPATH to the pixy checkout if needed)
#   N_WINDOWS   windows per ploidy (default 80000; ~5e-4 MCSE on the floor)
#
# run from analyses/polyploid_multiallele/analysis/ ; then re-run
# 07_thetaw_tajimaD.R and figures/Inline_statistical_tests.R
set -euo pipefail

PYTHON="${PYTHON:-python3}"
export N_WINDOWS="${N_WINDOWS:-80000}"

echo "Regenerating D reference: N_WINDOWS=${N_WINDOWS}, python=${PYTHON}"
pids=()
for p in 2 4 6 8; do
    "$PYTHON" calibrate_neutral_tajimaD.py "$p" > "calib_p${p}.log" 2>&1 &
    pids+=("$!")
    echo "  launched ploidy $p (pid $!)"
done

fail=0
for pid in "${pids[@]}"; do
    wait "$pid" || fail=1
done
if [[ "$fail" -ne 0 ]]; then
    echo "ERROR: a ploidy worker failed; see calib_p*.log" >&2
    exit 1
fi

# merge partials (header once, then data rows in ploidy order)
head -1 neutral_tajimaD_reference.p2.tsv > neutral_tajimaD_reference.tsv
for p in 2 4 6 8; do
    tail -n +2 "neutral_tajimaD_reference.p${p}.tsv" >> neutral_tajimaD_reference.tsv
done
rm -f neutral_tajimaD_reference.p*.tsv
echo "wrote neutral_tajimaD_reference.tsv"
cat neutral_tajimaD_reference.tsv
