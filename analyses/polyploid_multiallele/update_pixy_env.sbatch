#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=1G
#SBATCH -p short
#SBATCH --time=1:00:00
#SBATCH --job-name=update_pixy_env
#SBATCH --output=logs/%x.%j.out

# repoint the pixy conda env at the commit pinned in envs/pixy.yml (02d6d91:
# multiallelic Hudson FST + --fst_biallelic, unreleased as of 2026-07-14), then
# smoke-test one replicate before any production array starts.
#
# reinstalls only the pixy package, leaving htslib/samtools/python untouched.
# --force-reinstall needed because pip caches by URL and would treat the installed
# 2.1.2 as satisfying the requirement.
#
# gate the rerun via dependency so a bad env can't produce 40,000 garbage reps:
#   JOB=$(sbatch --parsable update_pixy_env.sbatch)
#   SBATCH_DEPENDENCY=afterok:$JOB PARTITION=short TIMELIMIT=01:45:00 \
#       ./01_run_all.sh --include 'dxy_.*_theta'

set -euo pipefail

# do NOT cd "$(dirname "$(readlink -f "$0")")" here — under sbatch $0 is Slurm's
# spool copy, not the repo file, so it resolves to an unwritable spool dir. sbatch
# already starts in the submission directory.
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
mkdir -p logs

PIN=02d6d91
SMOKE_ARM=dxy_dip_2pop_theta100   # highest theta => most multiallelic sites

set +u
eval "$(conda shell.bash hook)"
conda activate pixy
set -u

echo "[$(date '+%F %T')] before:"
pip show pixy 2>/dev/null | grep -E '^Version' || echo "  (pixy not found)"

pip install --force-reinstall --no-deps "git+https://github.com/ksamuk/pixy@${PIN}"

echo "[$(date '+%F %T')] after:"
pip show pixy 2>/dev/null | grep -E '^Version'

# 1. new cli surface exists
echo "=== --fst_biallelic present? ==="
if pixy --help 2>&1 | grep -q -- --fst_biallelic; then
    echo "OK: --fst_biallelic found"
else
    echo "FAIL: --fst_biallelic missing — env did not update" >&2
    exit 1
fi

# 2. multiallelic code path exists
echo "=== _hudson_fst present in pixy.calc? ==="
python - <<'PY'
import inspect

from pixy import calc

src = inspect.getsource(calc)
assert "_hudson_fst" in src, "FAIL: _hudson_fst not in pixy.calc"
print("OK: _hudson_fst present in pixy.calc")
PY

# 3. one real replicate end-to-end
# N_ARRAY_TASKS = n_replicates => exactly one rep in task 1
echo "=== smoke: one replicate of ${SMOKE_ARM} ==="
N_REPS=$(awk -F'\t' -v arm="$SMOKE_ARM" '
    NR==1 { for (i=1; i<=NF; i++) if ($i == "n_replicates") c=i; next }
    $1 == arm { print $c; exit }
' config/sim_params.tsv)

ARM_NAME="$SMOKE_ARM" N_ARRAY_TASKS="$N_REPS" SLURM_ARRAY_TASK_ID=1 \
    bash 02_run_arm.sh

SMOKE_TSV="data/${SMOKE_ARM}.part_0001.tsv"
SMOKE_HDR="${SMOKE_TSV}.header"
echo "--- header ---"
cat "$SMOKE_HDR"
echo "--- first row ---"
head -1 "$SMOKE_TSV"
echo "--- row count (expect 10 windows) ---"
wc -l < "$SMOKE_TSV"

# 4. the three fst columns are present and behave
# the pre-fix aggregated tsv is passed for comparison, not as a pass/fail gate:
# its env was patched in place over the project, so a mismatch is ambiguous (fix
# vs accumulated drift) and is only reported.
BACKUP="data/aggregated/_pre_fstfix_2026-07-14/${SMOKE_ARM}.tsv"
python - "$SMOKE_HDR" "$SMOKE_TSV" "$BACKUP" <<'PY'
import sys

hdr = open(sys.argv[1]).read().strip().split("\t")
rows = [line.rstrip("\n").split("\t") for line in open(sys.argv[2]) if line.strip()]

need = [
    "dxy_new", "dxy_new_multi",
    "pi_new", "pi_new_multi",
    "fst_hudson_new", "fst_hudson_num_new", "fst_hudson_den_new",
    "fst_hudson_new_multi", "fst_hudson_num_new_multi", "fst_hudson_den_new_multi",
    "fst_hudson_new_multi_fstbi", "fst_hudson_num_new_multi_fstbi",
    "fst_hudson_den_new_multi_fstbi",
]
missing = [c for c in need if c not in hdr]
if missing:
    sys.exit(f"FAIL: missing columns: {missing}\ngot: {hdr}")
print(f"OK: all {len(need)} expected columns present")


def pooled(prefix):
    ni, di = hdr.index(f"fst_hudson_num_{prefix}"), hdr.index(f"fst_hudson_den_{prefix}")
    num = sum(float(r[ni]) for r in rows if r[ni] not in ("NA", ""))
    den = sum(float(r[di]) for r in rows if r[di] not in ("NA", ""))
    return num / den if den else float("nan")


f_new = pooled("new")
f_multi = pooled("new_multi")
f_fstbi = pooled("new_multi_fstbi")
print(f"pooled FST  new={f_new:.6f}  new_multi={f_multi:.6f}  new_multi_fstbi={f_fstbi:.6f}")

# new_multi must now genuinely differ from biallelic (the whole point);
# new_multi_fstbi should sit close to new (differs only by sites that pass the
# read filter but not is_biallelic())
if abs(f_multi - f_new) < 1e-9:
    sys.exit("FAIL: new_multi FST identical to new — multiallelic path not active")
print(f"OK: new_multi differs from new by {100*(f_multi-f_new)/f_new:+.2f}%")
print(f"NOTE: new_multi_fstbi vs new: {100*(f_fstbi-f_new)/f_new:+.3f}% "
      f"(expected small; this is the old mislabelled new_multi estimand)")

# informational: same seed, so compare against the pre-fix run.
# prediction: the new new_multi_fstbi column reproduces the old new_multi column
# (both = read multiallelic in, then filter fst to biallelic), while the new
# new_multi should NOT (it is now genuinely multiallelic).
try:
    bak = open(sys.argv[3])
except OSError as e:
    print(f"NOTE: no pre-fix backup to compare ({e}); skipping")
    sys.exit(0)

bhdr = bak.readline().rstrip("\n").split("\t")
brows = {}
for line in bak:
    f = line.rstrip("\n").split("\t")
    brows[(f[bhdr.index("replicate")], f[bhdr.index("window_pos_1")])] = f

rep = rows[0][hdr.index("replicate")]
print(f"--- vs pre-fix backup, replicate {rep} ---")


def bpooled(prefix):
    ni, di = bhdr.index(f"fst_hudson_num_{prefix}"), bhdr.index(f"fst_hudson_den_{prefix}")
    num = den = 0.0
    n = 0
    for r in rows:
        key = (r[hdr.index("replicate")], r[hdr.index("window_pos_1")])
        b = brows.get(key)
        if b is None or b[ni] in ("NA", ""):
            continue
        num += float(b[ni])
        den += float(b[di])
        n += 1
    return (num / den if den else float("nan")), n


ob_new, n_matched = bpooled("new")
ob_multi, _ = bpooled("new_multi")
if n_matched == 0:
    print("NOTE: replicate not found in backup; cannot compare")
    sys.exit(0)

print(f"matched {n_matched} windows")
print(f"  old new       = {ob_new:.9f}   new new       = {f_new:.9f}  "
      f"delta={100*(f_new-ob_new)/ob_new:+.4f}%")
print(f"  old new_multi = {ob_multi:.9f}   new fstbi     = {f_fstbi:.9f}  "
      f"delta={100*(f_fstbi-ob_multi)/ob_multi:+.4f}%   <- expect ~0 if diagnosis right")
print(f"  old new_multi = {ob_multi:.9f}   new new_multi = {f_multi:.9f}  "
      f"delta={100*(f_multi-ob_multi)/ob_multi:+.4f}%   <- expect NON-zero (the fix)")
PY

# the production array's task 1 rewrites these, but leave no half-arm behind
rm -f "$SMOKE_TSV" "$SMOKE_HDR"
rm -rf "data/_work/${SMOKE_ARM}"

echo "[$(date '+%F %T')] update_pixy_env + smoke complete"
