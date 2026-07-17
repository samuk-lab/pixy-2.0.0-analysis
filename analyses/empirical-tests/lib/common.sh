#!/usr/bin/env bash
# shared helpers sourced by the per-arm sbatch scripts

info() {
    printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"
}

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

require_tool() {
    command -v "$1" >/dev/null 2>&1 || die "Required tool not on PATH: $1"
}

# activate a named conda env, mamba or conda
activate_env() {
    local env_name="$1"
    [[ -n "${env_name}" ]] || die "activate_env: env name is required"
    # hpcc ~/.bashrc breaks under set -euo pipefail (unbound MODULESHOME,
    # nonzero module load); relax around the source
    if [[ -f ~/.bashrc ]]; then
        set +eu
        # shellcheck disable=SC1091
        source ~/.bashrc
        set -eu
    fi
    set +u
    local _hook=""
    if command -v mamba >/dev/null 2>&1; then
        _hook="$(mamba shell hook --shell bash 2>/dev/null)" || true
    fi
    if [[ -n "${_hook}" ]]; then
        eval "${_hook}"
        mamba activate "${env_name}"
    elif command -v conda >/dev/null 2>&1; then
        eval "$(conda shell.bash hook)"
        conda activate "${env_name}"
    else
        set -u
        die "Neither mamba nor conda found on PATH"
    fi
    set -u
}

# look up a field from config/arms.tsv for an arm_id
arm_config_lookup() {
    local arm_id="$1"
    local field="$2"
    local arms_tsv="${PROJECT_DIR:?PROJECT_DIR not set}/config/arms.tsv"
    [[ -s "${arms_tsv}" ]] || die "arms.tsv not found at ${arms_tsv}"
    python3 - "${arms_tsv}" "${arm_id}" "${field}" <<'PY'
import csv, sys
arms_tsv, arm_id, field = sys.argv[1:4]
with open(arms_tsv) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        if row["arm_id"] == arm_id:
            if field not in row:
                sys.exit(f"arms.tsv has no field '{field}'")
            print(row[field])
            sys.exit(0)
sys.exit(f"arm_id '{arm_id}' not found in {arms_tsv}")
PY
}
