#!/bin/bash -l
# simulate diploid 10 Mb vcfs with vcfsim, one per array task
# params match prior benchmark (d. melanogaster Ne/mu); output bgzipped + tabix'd
# env: vcfsim (envs/vcfsim.yml)
#SBATCH --job-name=vcfsim_10Mb
#SBATCH --output=logs/vcfsim_10Mb_%A_%a.out
#SBATCH --error=logs/vcfsim_10Mb_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
# re-tune from seff after first run
#SBATCH --partition=intel
#SBATCH --time=0-08:00:00
#SBATCH --array=1-100%25

set -euo pipefail
umask 002

##########
# parameters
##########
PROJECT_DIR="${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}}"
OUT_DIR="${OUT_DIR:-${PROJECT_DIR}/data/vcfs/dm_10Mb}"
LOG_DIR="${LOG_DIR:-${PROJECT_DIR}/logs}"

CHROM=1
REPS=1
SEQ_LEN=10000000        # 10 Mb
PLOIDY=2
NE=1720600              # drosophila melanogaster
MU=0.00000000549
MISSING_SITES=0
MISSING_GENOS=0
SAMPLE_SIZE=20

##########
# env
##########
if command -v conda >/dev/null 2>&1; then
  eval "$(conda shell.bash hook)"
  conda activate vcfsim
fi
module load htslib 2>/dev/null || true
command -v vcfsim >/dev/null 2>&1 || { echo "vcfsim not on PATH" >&2; exit 1; }
command -v bgzip  >/dev/null 2>&1 || { echo "bgzip not on PATH"  >&2; exit 1; }
command -v tabix  >/dev/null 2>&1 || { echo "tabix not on PATH"  >&2; exit 1; }

mkdir -p "${OUT_DIR}" "${LOG_DIR}"
cd "${OUT_DIR}"

seed="${SLURM_ARRAY_TASK_ID}"
prefix="dm_sim_vcf_seed_${seed}"

if [[ -f "${prefix}.vcf.gz" && -f "${prefix}.vcf.gz.tbi" ]]; then
  echo "Seed ${seed} already exists; skipping."
  exit 0
fi

start_time=$(date +%s)
echo "[$(date)] Starting seed=${seed} on ${HOSTNAME}"

# each task simulates in its own directory so it cannot see another task's
# intermediates. previously every task ran in OUT_DIR and normalised with
# `ls "${prefix}"*.vcf | head -n1`; seed 11's glob matched seed 1's unrenamed
# intermediate and seed 11 was bgzipped from seed 1's simulation.
work_dir="${OUT_DIR}/_work/seed_${seed}"
rm -rf "${work_dir}"
mkdir -p "${work_dir}"
cd "${work_dir}"

if ! vcfsim \
     --chromosome                 "${CHROM}" \
     --replicates                 "${REPS}" \
     --seed                       "${seed}" \
     --sequence_length            "${SEQ_LEN}" \
     --ploidy                     "${PLOIDY}" \
     --Ne                         "${NE}" \
     --mu                         "${MU}" \
     --percent_missing_sites      "${MISSING_SITES}" \
     --percent_missing_genotypes  "${MISSING_GENOS}" \
     --output_file                "${prefix}" \
     --sample_size                "${SAMPLE_SIZE}"
then
  echo "vcfsim failed for seed ${seed}" >&2
  exit 1
fi

# vcfsim appends the seed to --output_file (vcfsim/__main__.py:107), so the name
# is deterministic; match it exactly rather than globbing
generated="${prefix}${seed}.vcf"
if [[ ! -f "${generated}" ]]; then
  echo "Seed ${seed}: expected ${generated}, found: $(ls -1 ./*.vcf 2>/dev/null | tr '\n' ' ')" >&2
  exit 2
fi
mv "${generated}" "${prefix}.vcf"

bgzip -f "${prefix}.vcf"
tabix -f -p vcf "${prefix}.vcf.gz"
mv -f "${prefix}.vcf.gz" "${prefix}.vcf.gz.tbi" "${OUT_DIR}/"
cd "${OUT_DIR}"
rm -rf "${work_dir}"

elapsed=$(( $(date +%s) - start_time ))
printf "Seed %4d done in %02d:%02d:%02d\n" \
  "${seed}" "$((elapsed/3600))" "$(((elapsed%3600)/60))" "$((elapsed%60))"
