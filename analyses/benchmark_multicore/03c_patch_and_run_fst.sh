#!/usr/bin/env bash
# one-shot: patch new_pixy with github-head pixy, then submit the fst multicore benchmark
# array jobs carry afterok on the patch job
# aggregate after with 04_aggregate_summaries.sh

set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-$(cd "$(dirname "$0")" && pwd)}"

patch_jid=$(PROJECT_DIR="${PROJECT_DIR}" bash "${PROJECT_DIR}/03a_patch_pixy_env.sh")

DEP="afterok:${patch_jid}" \
PROJECT_DIR="${PROJECT_DIR}" \
bash "${PROJECT_DIR}/03b_submit_pixy_new_fst_only.sh"
