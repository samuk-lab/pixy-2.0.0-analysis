#!/bin/bash -l
# one-time conda env setup for the pipeline
# create the four conda envs from envs/*.yml; run on a login/interactive node
#
#   ./00_create_conda_envs.sh             # create envs that don't exist
#   ./00_create_conda_envs.sh --update    # refresh existing envs from yml
#   ./00_create_conda_envs.sh --force     # remove + recreate every env
#   ./00_create_conda_envs.sh -h | --help # show this help
#
# envs:
#   vcfsim    - simulation scripts; vcfsim from bioconda
#   pixy      - current pixy from conda-forge, htslib/samtools from bioconda
#   old_pixy  - pixy 0.95.01 from conda-forge, htslib/samtools from bioconda
#   vcftools  - vcftools comparator (pi / Weir-Cockerham FST) + bcftools for
#               header/concat in the mixed-ploidy arm
#
# then smoke-tests each env and applies the upstream-source patches below
#
#

set -euo pipefail

ENVS=(vcfsim pixy old_pixy vcftools)
ENV_DIR=envs

usage() {
    sed -n '3,20p' "$0" | sed 's/^# \{0,1\}//'
}

MODE=create
for arg in "$@"; do
    case "$arg" in
        --force)  MODE=force ;;
        --update) MODE=update ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown argument: $arg" >&2; usage >&2; exit 2 ;;
    esac
done

# run from the script's dir so envs/*.yml resolves
cd "$(dirname "$(readlink -f "$0")")"

if ! command -v conda >/dev/null 2>&1; then
    echo "[00_create_envs] ERROR: conda not on PATH." >&2
    echo "  Load conda first, e.g. 'module load anaconda' or source your" >&2
    echo "  miniconda/conda init script, then re-run." >&2
    exit 1
fi

# enable conda activate / env list / run inside scripts
eval "$(conda shell.bash hook)"

# clear the cached package index up front: a stale index predating a freshly
# published pin can make a solve fail *after* the old env was removed (--force),
# which once destroyed several envs during the 2.2.3 rerun. a clean index also
# makes a from-scratch build reproduce the pinned versions in envs/*.yml.
conda clean --index-cache -y >/dev/null 2>&1 || true

env_exists() {
    conda env list | awk 'NF && $1 !~ /^#/ {print $1}' | grep -qxF "$1"
}

for env in "${ENVS[@]}"; do
    yml="${ENV_DIR}/${env}.yml"
    if [[ ! -f "$yml" ]]; then
        echo "[00_create_envs] ERROR: $yml not found (expected in $ENV_DIR/)" >&2
        exit 1
    fi

    if env_exists "$env"; then
        case "$MODE" in
            create)
                echo "[00_create_envs] $env: already exists; skipping"
                echo "                 (use --update to refresh from $yml, --force to recreate)"
                continue
                ;;
            force)
                echo "[00_create_envs] $env: removing existing env (--force)"
                conda env remove -n "$env" -y
                ;;
            update)
                echo "[00_create_envs] $env: updating from $yml"
                conda env update -n "$env" -f "$yml" --prune
                continue
                ;;
        esac
    fi

    echo "[00_create_envs] $env: creating from $yml"
    conda env create -f "$yml"
done

echo
echo "[00_create_envs] Sanity checks:"

# vcfsim: check imports, show cli if present
conda run -n vcfsim python - <<'PY'
import msprime
import vcfsim
print("  vcfsim    OK  (imports: vcfsim, msprime)")
PY
if conda run -n vcfsim vcfsim --help >/dev/null 2>&1; then
    echo "             vcfsim CLI callable"
else
    echo "             note: vcfsim CLI not detected, but Python import succeeded"
fi

# vcftools env: expose both vcftools and bcftools
if conda run -n vcftools vcftools --version >/dev/null 2>&1; then
    vcftools_version=$(conda run -n vcftools vcftools --version 2>&1 | head -1)
    bcftools_version=$(conda run -n vcftools bcftools --version 2>&1 | head -1)
    printf "  %-10s OK  (%s; %s)\n" "vcftools" "$vcftools_version" "$bcftools_version"
else
    echo "  vcftools  FAILED -- vcftools not callable in env" >&2
    exit 1
fi

# pixy envs: pixy and samtools callable from the env
for env in pixy old_pixy; do
    pixy_version=$(conda run -n "$env" pixy --version 2>&1 | head -1 || true)
    samtools_version=$(conda run -n "$env" samtools --version 2>&1 | head -1 || true)

    if [[ -z "$pixy_version" ]]; then
        echo "  $env FAILED -- pixy not callable in env" >&2
        exit 1
    fi
    if [[ -z "$samtools_version" ]]; then
        echo "  $env FAILED -- samtools not callable in env" >&2
        exit 1
    fi

    printf "  %-10s OK  (%s; %s)\n" "$env" "$pixy_version" "$samtools_version"
done

##########
# patches applied to installed env source files
##########
# one upstream bug worked around, a literal find/replace in one env file.
# apply_patch is idempotent: skips if already patched, warns if neither old nor
# new pattern matches (upstream changed the line).
#
#   1. vcfsim SimulatorClass.row_changes() slices a pandas Series with .values,
#      read-only in newer numpy/pandas. fix: add .copy().
#
# two old_pixy patches were dropped when the legacy pin moved 0.93.1 -> 0.95.01:
# gzipped-vcf chromosome detection (0.95.01 selects gunzip -c via cat_prog) and
# the --interval_start default (0.95.01 already defaults to min(pos_array)).
# Both were fixed upstream between the two releases.

apply_patch() {
    # apply_patch <env> <path-glob-under-env-prefix> <label> <old> <new>
    OLD="$4" NEW="$5" PATH_GLOB="$2" LABEL="$3" ENV_NAME="$1" \
        python3 - <<'PATCH_PY'
import glob, os, subprocess, sys
env = os.environ["ENV_NAME"]
prefix = subprocess.check_output(
    ["conda", "run", "-n", env, "python", "-c", "import sys; print(sys.prefix)"],
    text=True,
).strip()
label = os.environ["LABEL"]
paths = glob.glob(os.path.join(prefix, os.environ["PATH_GLOB"]))
if not paths:
    print(f"  {label}: WARNING - file not found under {prefix}", file=sys.stderr)
    sys.exit(0)
path = paths[0]
old, new = os.environ["OLD"], os.environ["NEW"]
with open(path) as f:
    txt = f.read()
if new in txt:
    print(f"  {label}: already applied")
elif old not in txt:
    print(f"  {label}: WARNING - expected pattern not found in {path}", file=sys.stderr)
else:
    n = txt.count(old)
    with open(path, "w") as f:
        f.write(txt.replace(old, new))
    print(f"  {label}: applied ({n} replacement(s) in {path})")
PATCH_PY
}

echo "[00_create_envs] Applying upstream-source patches..."

apply_patch vcfsim "lib/python*/site-packages/vcfsim/SimulatorClass.py" \
    "vcfsim read-only .values slice" \
    "altlist = row[col_start:col_end+1].values" \
    "altlist = row[col_start:col_end+1].values.copy()"

echo
echo "[00_create_envs] All set. Next step:"
echo "    ./01_run_all.sh"
