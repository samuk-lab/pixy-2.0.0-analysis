#!/usr/bin/env bash
#SBATCH --job-name=compare
#SBATCH --partition=short
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=logs/compare_%j.out
#SBATCH --error=logs/compare_%j.err
# per arm: join biallelic + multi pixy outputs (pi, dxy, fst) on
# (chromosome, window_pos_1) into data/aggregated/<arm>.joined.tsv
# fst col is avg_wc_fst (wc) or avg_hudson_fst (hudson); handles both
set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}}"
ARMS_TSV="${PROJECT_DIR}/config/arms.tsv"
OUT_DIR="${PROJECT_DIR}/data/aggregated"
mkdir -p "${OUT_DIR}"

# shellcheck disable=SC1091
source "${PROJECT_DIR}/lib/common.sh"
activate_env empirical_tests_pixy   # python + pandas here
require_tool python3

while IFS=$'\t' read -r arm_id rest; do
    [[ "${arm_id}" == "arm_id" || "${arm_id}" == \#* ]] && continue

    info "[${arm_id}] joining biallelic + multi pixy outputs"
    DATASET_DIR="${PROJECT_DIR}/data/${arm_id}"
    BI_DIR="${DATASET_DIR}/pixy_out/biallelic"
    MU_DIR="${DATASET_DIR}/pixy_out/multi"
    OUT_TSV="${OUT_DIR}/${arm_id}.joined.tsv"

    for f in \
        "${BI_DIR}/${arm_id}_biallelic_pidxy_pi.txt"  \
        "${BI_DIR}/${arm_id}_biallelic_pidxy_dxy.txt" \
        "${BI_DIR}/${arm_id}_biallelic_fst_fst.txt"   \
        "${MU_DIR}/${arm_id}_multi_pidxy_pi.txt"      \
        "${MU_DIR}/${arm_id}_multi_pidxy_dxy.txt"     \
        "${MU_DIR}/${arm_id}_multi_fst_fst.txt"       \
    ; do
        [[ -s "${f}" ]] || die "[${arm_id}] missing pixy output: ${f}"
    done

    python3 - "${arm_id}" "${BI_DIR}" "${MU_DIR}" "${OUT_TSV}" <<'PY'
import sys
import pandas as pd

arm_id, bi_dir, mu_dir, out_tsv = sys.argv[1:5]

def load(prefix_dir, variant, stat):
    # pixy cols: pop1[, pop2], chromosome, window_pos_1, window_pos_2, then stat cols
    path = f"{prefix_dir}/{arm_id}_{variant}_{'pidxy' if stat in ('pi','dxy') else 'fst'}_{stat}.txt"
    df = pd.read_csv(path, sep="\t")
    return df

def flatten_col(c, tag):
    # flatten the 2-level MultiIndex from pivot_table: keep index cols verbatim,
    # rename (value, pop) pairs to <value>_<pop>_<tag>
    if isinstance(c, str):
        return c
    if c[1] == "":
        return c[0]
    return f"{c[0]}_{c[1]}_{tag}"

def join_pi(bi, mu):
    # pi has one row per (pop, window); pivot wide on pop, then merge bi+mu
    def wide(df, tag):
        w = df.pivot_table(index=["chromosome","window_pos_1","window_pos_2"],
                           columns="pop",
                           values=["avg_pi","no_sites"],
                           aggfunc="first").reset_index()
        w.columns = [flatten_col(c, tag) for c in w.columns]
        return w
    return wide(bi, "biallelic"), wide(mu, "multi")

# pi
pi_bi = load(bi_dir, "biallelic", "pi")
pi_mu = load(mu_dir, "multi",     "pi")
pi_bi_w, pi_mu_w = join_pi(pi_bi, pi_mu)
pi = pd.merge(pi_bi_w, pi_mu_w, on=["chromosome","window_pos_1","window_pos_2"], how="outer")

# dxy — already one row per window with pop1/pop2 cols
dxy_bi = load(bi_dir, "biallelic", "dxy").rename(columns={"avg_dxy":"dxy_biallelic","no_sites":"n_sites_dxy_biallelic"})
dxy_mu = load(mu_dir, "multi",     "dxy").rename(columns={"avg_dxy":"dxy_multi",    "no_sites":"n_sites_dxy_multi"})
dxy = pd.merge(
    dxy_bi[["pop1","pop2","chromosome","window_pos_1","window_pos_2","dxy_biallelic","n_sites_dxy_biallelic"]],
    dxy_mu[["pop1","pop2","chromosome","window_pos_1","window_pos_2","dxy_multi","n_sites_dxy_multi"]],
    on=["pop1","pop2","chromosome","window_pos_1","window_pos_2"], how="outer",
)

# fst — col is avg_wc_fst or avg_hudson_fst depending on --fst_type
def load_fst(prefix_dir, variant, tag):
    df = load(prefix_dir, variant, "fst")
    fst_col = next((c for c in ("avg_wc_fst","avg_hudson_fst") if c in df.columns), None)
    if fst_col is None:
        raise SystemExit(f"no FST estimator column found in {variant} fst output; got {list(df.columns)}")
    return df.rename(columns={fst_col: f"fst_{tag}", "no_snps": f"n_snps_fst_{tag}"})

fst_bi = load_fst(bi_dir, "biallelic", "biallelic")
fst_mu = load_fst(mu_dir, "multi",     "multi")
fst_keep_bi = [c for c in ("pop1","pop2","chromosome","window_pos_1","window_pos_2","fst_biallelic","n_snps_fst_biallelic") if c in fst_bi.columns]
fst_keep_mu = [c for c in ("pop1","pop2","chromosome","window_pos_1","window_pos_2","fst_multi","n_snps_fst_multi") if c in fst_mu.columns]
fst = pd.merge(
    fst_bi[fst_keep_bi],
    fst_mu[fst_keep_mu],
    on=["pop1","pop2","chromosome","window_pos_1","window_pos_2"], how="outer",
)

# merge all on (chromosome, window_pos_1, window_pos_2)
out = pi.merge(dxy, on=["chromosome","window_pos_1","window_pos_2"], how="outer") \
        .merge(fst, on=["pop1","pop2","chromosome","window_pos_1","window_pos_2"], how="outer")

# reorder cols
front = ["chromosome","window_pos_1","window_pos_2","pop1","pop2"]
out = out[[c for c in front if c in out.columns] +
          [c for c in out.columns if c not in front]]
out.sort_values(["chromosome","window_pos_1"], inplace=True)
out.to_csv(out_tsv, sep="\t", index=False, na_rep="NA")
print(f"wrote {out_tsv}: {len(out)} rows, {len(out.columns)} cols", file=sys.stderr)
PY

    info "[${arm_id}] -> ${OUT_TSV}"
    head -3 "${OUT_TSV}" || true
done < "${ARMS_TSV}"

info "[07_compare] done. Aggregated TSVs:"
ls -la "${OUT_DIR}/"
