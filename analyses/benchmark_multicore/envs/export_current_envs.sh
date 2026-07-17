#!/usr/bin/env bash
# conda-export each named env to <env>.lock.yml next to this script

set -euo pipefail

ENV_DIR="$(cd "$(dirname "$0")" && pwd)"

for env in vcfsim pixy1 pixy2; do
  if conda env list | awk '{print $1}' | grep -qx "${env}"; then
    out="${ENV_DIR}/${env}.lock.yml"
    echo "Exporting ${env} -> ${out}"
    conda env export -n "${env}" --no-builds > "${out}"
  else
    echo "WARN: conda env '${env}' not found; skipping." >&2
  fi
done
