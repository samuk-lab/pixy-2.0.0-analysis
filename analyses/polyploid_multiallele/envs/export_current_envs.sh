#!/bin/bash
# snapshot pinned package versions from each conda env into *.lock.yml
# run on the machine where the pipeline actually executes

set -euo pipefail

cd "$(dirname "$0")"

for env in pixy old_pixy vcfsim; do
    if conda env list | awk '{print $1}' | grep -qx "$env"; then
        conda env export -n "$env" > "${env}.lock.yml"
        echo "Wrote ${env}.lock.yml"
    else
        echo "Skipping $env — env not found"
    fi
done
