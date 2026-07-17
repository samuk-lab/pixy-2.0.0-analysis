#!/usr/bin/env bash
# sbatch job to swap conda pixy in new_pixy for github head (ksamuk/pixy default branch)
# run before 03b_submit_pixy_new_fst_only.sh to benchmark unreleased code
# pip --no-deps --force-reinstall: only pixy is replaced, conda deps untouched
# prints patch job id to stdout for dep chaining

set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
LOG_DIR="${PROJECT_DIR}/logs"
GITHUB_URL="https://github.com/ksamuk/pixy.git"
ENV_NAME="new_pixy"

mkdir -p "${LOG_DIR}"

SBATCH="$(command -v sbatch 2>/dev/null || \
          echo /opt/linux/rocky/8.x/x86_64/pkgs/slurm/24.11.1/bin/sbatch)"
[[ -x "${SBATCH}" ]] || { echo "sbatch not found" >&2; exit 1; }

patch_jid=$("${SBATCH}" \
  --job-name="patch_pixy_git" \
  --partition=short \
  --ntasks=1 --cpus-per-task=1 \
  --mem=8G --time=00:20:00 \
  --output="${LOG_DIR}/patch_pixy_%j.out" \
  --error="${LOG_DIR}/patch_pixy_%j.err" \
  --wrap="source ~/.bashrc
          echo 'Before patch:'; conda run -n '${ENV_NAME}' pixy --version
          conda run -n '${ENV_NAME}' pip install --no-deps --force-reinstall 'git+${GITHUB_URL}'
          echo 'After patch:';  conda run -n '${ENV_NAME}' pixy --version" \
  | awk '{print $NF}')

echo "  [patch] ${ENV_NAME} ← GitHub head → job ${patch_jid}" >&2
echo "${patch_jid}"
