#!/usr/bin/env bash
# concat every per-cell summary tsv in data/results/ into one long df
# adds pixy_version, statistic, n_cores; de-dupes by (version,stat,cores,seed)
# keeping last row (drops retry rows)
# writes data/results/all_cells_long.tsv

set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
RESULTS_DIR="${PROJECT_DIR}/data/results"
OUT="${RESULTS_DIR}/all_cells_long.tsv"

[[ -d "${RESULTS_DIR}" ]] || { echo "No ${RESULTS_DIR}" >&2; exit 1; }

shopt -s nullglob
files=( "${RESULTS_DIR}"/pixy_*.tsv )
mapfile -t files < <(printf "%s\n" "${files[@]}" | grep -v '/all_cells_long.tsv$')

if [[ ${#files[@]} -eq 0 ]]; then
  echo "No summary TSVs found in ${RESULTS_DIR}" >&2
  exit 1
fi

printf "pixy_version\tstatistic\tn_cores\tseed\tvalue\telapsed_us\telapsed_s\tuser_s\tsys_s\trss_kb\tstatus\n" > "${OUT}"

emit_rows() {
  # $1=file $2=pixy_version $3=statistic $4=default_cores
  local f="$1" ver="$2" stat="$3" dc="$4"
  awk -v ver="${ver}" -v stat="${stat}" -v dc="${dc}" '
    BEGIN { FS=OFS="\t" }
    NR==1 {
      for (i=1;i<=NF;i++) {
        if      ($i=="Seed")                    iSeed=i
        else if ($i=="Cores" || $i=="N_Cores")  iCores=i
        else if ($i=="Pi"||$i=="Dxy"||$i=="Fst") iVal=i
        else if ($i=="Elapsed_us")              iEus=i
        else if ($i=="Elapsed_s")               iEs=i
        else if ($i=="User_s")                  iUs=i
        else if ($i=="Sys_s")                   iSs=i
        else if ($i=="RSS_kb")                  iRss=i
        else if ($i=="Status")                  iStat=i
      }
      next
    }
    {
      seed  = (iSeed  ? $iSeed  : "NA")
      cores = (iCores ? $iCores : dc)
      val   = (iVal   ? $iVal   : "NA")
      print ver, stat, cores, seed, val,
            (iEus?$iEus:""), (iEs?$iEs:""), (iUs?$iUs:""), (iSs?$iSs:""),
            (iRss?$iRss:""), (iStat?$iStat:"")
    }
  ' "${f}"
}

for f in "${files[@]}"; do
  base=$(basename "${f}" .tsv)
  case "${base}" in
    # 0.95.01 single core (cores col already "0.95.01" in tsv)
    pixy_old_pi_10Mb)  emit_rows "${f}" "0.95.01" pi  "0.95.01" ;;
    pixy_old_dxy_10Mb) emit_rows "${f}" "0.95.01" dxy "0.95.01" ;;
    pixy_old_fst_10Mb) emit_rows "${f}" "0.95.01" fst "0.95.01" ;;
    # 2.2.3 any core count
    pixy_new_pi_10Mb_cores_*)  c="${base##*_}"; emit_rows "${f}" "2.2.3" pi  "${c}" ;;
    pixy_new_dxy_10Mb_cores_*) c="${base##*_}"; emit_rows "${f}" "2.2.3" dxy "${c}" ;;
    pixy_new_fst_10Mb_cores_*) c="${base##*_}"; emit_rows "${f}" "2.2.3" fst "${c}" ;;
    all_cells_long) ;;
    *) echo "WARN: unrecognized file ${base}" >&2 ;;
  esac
done >> "${OUT}"

# dedup: keep last row per (version, stat, cores, seed)
{
  head -1 "${OUT}"
  tail -n +2 "${OUT}" | awk -F'\t' '
    { line[NR] = $0; key[NR] = $1 FS $2 FS $3 FS $4; last[key[NR]] = NR }
    END { for (i = 1; i <= NR; i++) if (last[key[i]] == i) print line[i] }
  '
} > "${OUT}.dedup"
mv -f "${OUT}.dedup" "${OUT}"

n=$(wc -l < "${OUT}")
echo "Wrote ${OUT} ($((n-1)) data rows from ${#files[@]} cell files)"
