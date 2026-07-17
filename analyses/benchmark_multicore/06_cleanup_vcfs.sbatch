#!/usr/bin/env bash
# delete each seed's vcf (.vcf.gz + .tbi) once all pixy jobs are done
# submitted by 05_submit_all_waves.sh with afterok on every pixy array job,
# so a vcf is removed only after all arms have read it
# safe to run manually for an already-processed seed
#SBATCH --job-name=cleanup_vcfs
#SBATCH --output=logs/cleanup_%A_%a.out
#SBATCH --error=logs/cleanup_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=0-00:05:00
#SBATCH --array=1-100%25

set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}}"
VCF_DIR="${PROJECT_DIR}/data/vcfs/dm_10Mb"

seed="${SLURM_ARRAY_TASK_ID}"
vcf="${VCF_DIR}/dm_sim_vcf_seed_${seed}.vcf.gz"

if [[ -f "${vcf}" ]]; then
  rm -f "${vcf}" "${vcf}.tbi"
  echo "Deleted VCF for seed ${seed}"
else
  echo "VCF for seed ${seed} already gone; nothing to do"
fi
