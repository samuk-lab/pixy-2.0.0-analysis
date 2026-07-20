#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=1G
#SBATCH -p batch
#SBATCH --output=logs/%x.%A_%a.out

# unified per-vcf, per-replicate pipeline. reads one arm row from
# config/sim_params.tsv, then for each replicate in this task's seed chunk:
#   1. simulate a vcf (single-chrom or per-chromosome mixed-ploidy)
#   2. bgzip + tabix
#   3. run pixy 2.0 (new + new_multi), all declared stats in one call
#   4. run pixy 0.95.01 if comparators has old_pixy (pi/dxy only; old pixy has
#      no fst-alone, watterson_theta, or tajima_d)
#   5. run vcftools if comparators has vcftools (windowed pi for 1-pop,
#      windowed Weir FST for 2-pop)
#   6. aggregate across variants into one row per window, appended to
#      data/${ARM_NAME}.part_NNNN.tsv
#   7. delete the vcf, index, per-rep work dir, and zarr scratch
#
# no vcfs survive past their statistics (cluster disk); re-runs reproduce the
# vcf from (arm, seed).
#
# submitted by 01_run_all.sh (defaults N_ARRAY_TASKS=100, CONCURRENT=50);
# manual submission can override.

set -euo pipefail

# required env vars
ARM_NAME="${ARM_NAME:?ARM_NAME must be set via sbatch --export=ARM_NAME=...}"
N_ARRAY_TASKS="${N_ARRAY_TASKS:-100}"
PARAMS_TSV="${PARAMS_TSV:-config/sim_params.tsv}"
WINDOW_SIZE="${WINDOW_SIZE:-10000}"
AGGREGATOR="${AGGREGATOR:-analysis/aggregate_one_rep.py}"

TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID not set (submit via sbatch with --array)}"
TASK_ID_PADDED=$(printf '%04d' "$TASK_ID")
N_CORES="${N_CORES:-${SLURM_CPUS_PER_TASK:-8}}"

echo "[$(date '+%F %T')] 02_run_arm.sh: arm=${ARM_NAME} task=${TASK_ID}/${N_ARRAY_TASKS} n_cores=${N_CORES}"

##########
# parse arm row
##########

declare -A ARM
while IFS=$'\t' read -r key value; do
    ARM[$key]="$value"
done < <(awk -F'\t' -v arm="$ARM_NAME" '
    NR==1 { for (i=1; i<=NF; i++) header[i]=$i; next }
    $1 == arm {
        for (i=1; i<=NF; i++) printf "%s\t%s\n", header[i], $i
        found=1; exit
    }
    END { if (!found) { print "ERROR: arm \""arm"\" not found in '"$PARAMS_TSV"'" > "/dev/stderr"; exit 2 } }
' "$PARAMS_TSV")

# bind locals from the arm row
STATS="${ARM[stats]}"
N_POPULATIONS="${ARM[n_populations]}"
SAMPLE_SIZE="${ARM[sample_size]}"
SAMPLE_SPLIT="${ARM[sample_split]}"
PLOIDY_CHR1="${ARM[ploidy_chr1]}"
PLOIDY_CHR2="${ARM[ploidy_chr2]}"
SEQ_LENGTH="${ARM[sequence_length]}"
NE="${ARM[Ne]}"
MU="${ARM[mu]}"
PCT_MISS_SITES="${ARM[percent_missing_sites]}"
PCT_MISS_GENO="${ARM[percent_missing_genotypes]}"
N_REPS="${ARM[n_replicates]}"
SEED_START_BASE="${ARM[seed_start]}"
SPLIT_TIME="${ARM[split_time]}"
FST_TYPE="${ARM[fst_type]}"
COMPARATORS="${ARM[comparators]}"
# recombination_rate added 2026-06-23; defaults to 0 for older config files
# that lack it (single-tree-per-rep). set to mu (r/mu = 1) for all non-legacy
# missingness-sweep, baseline, and theta-sweep arms; legacy 10k-rep "do not
# rerun" rows and pi-only theta-sweep arms stay at 0.
RECOMBINATION_RATE="${ARM[recombination_rate]:-0}"

##########
# compute this task's seed range
##########

if (( N_REPS % N_ARRAY_TASKS != 0 )); then
    echo "ERROR: n_replicates ($N_REPS) must be a multiple of N_ARRAY_TASKS ($N_ARRAY_TASKS)" >&2
    exit 1
fi
REPS_PER_TASK=$(( N_REPS / N_ARRAY_TASKS ))
TASK_SEED_START=$(( SEED_START_BASE + (TASK_ID - 1) * REPS_PER_TASK ))
TASK_SEED_END=$(( SEED_START_BASE + TASK_ID * REPS_PER_TASK - 1 ))
echo "[$(date '+%F %T')] task ${TASK_ID}: seeds ${TASK_SEED_START}..${TASK_SEED_END} (${REPS_PER_TASK} reps)"

##########
# per-task work area + partial output
##########

TASK_WORK="data/_work/${ARM_NAME}/task_${TASK_ID_PADDED}"
TASK_RESULTS="data/${ARM_NAME}.part_${TASK_ID_PADDED}.tsv"
mkdir -p "$TASK_WORK" "$(dirname "$TASK_RESULTS")" logs
: > "$TASK_RESULTS"   # truncate any previous attempt

##########
# build pops file (once per task)
##########

POPS_FILE="${TASK_WORK}/pops.tsv"
: > "$POPS_FILE"
if (( N_POPULATIONS == 1 )); then
    for ((i=1; i<=SAMPLE_SIZE; i++)); do
        printf "tsk_%d\tpop1\n" "$i" >> "$POPS_FILE"
    done
else
    if [[ "$SAMPLE_SPLIT" == "NA" || -z "$SAMPLE_SPLIT" ]]; then
        # default even split
        n_pop1=$(( SAMPLE_SIZE / 2 ))
        n_pop2=$(( SAMPLE_SIZE - n_pop1 ))
    else
        IFS=',' read -r n_pop1 n_pop2 <<< "$SAMPLE_SPLIT"
    fi
    if (( n_pop1 + n_pop2 != SAMPLE_SIZE )); then
        echo "ERROR: sample_split ($SAMPLE_SPLIT) does not sum to sample_size ($SAMPLE_SIZE)" >&2
        exit 1
    fi
    idx=1
    for ((j=0; j<n_pop1; j++)); do
        printf "tsk_%d\tpop1\n" "$idx" >> "$POPS_FILE"; idx=$((idx+1))
    done
    for ((j=0; j<n_pop2; j++)); do
        printf "tsk_%d\tpop2\n" "$idx" >> "$POPS_FILE"; idx=$((idx+1))
    done
fi

##########
# decide which variants run
##########

VARIANTS=(new new_multi)
# old pixy does pi / dxy only; stat-filtering happens downstream when we pick
# which old-pixy outputs to invoke. here we just decide whether to invoke it.
RUN_OLD_PIXY=0
RUN_VCFTOOLS=0
case ",${COMPARATORS}," in *,old_pixy,*) RUN_OLD_PIXY=1 ;; esac
case ",${COMPARATORS}," in *,vcftools,*) RUN_VCFTOOLS=1 ;; esac
echo "[$(date '+%F %T')] task ${TASK_ID}: stats=${STATS} variants=new,new_multi old=${RUN_OLD_PIXY} vcftools=${RUN_VCFTOOLS}"

##########
# translate stat names to pixy cli
##########
# config/sim_params.tsv uses short names; pixy uses long ones.
#   thetaw   -> watterson_theta
#   tajimaD  -> tajima_d
PIXY_STATS=$(echo "$STATS" | sed 's/thetaw/watterson_theta/g; s/tajimaD/tajima_d/g')
# stats list with fst stripped, for the no-fst-type sub-call below
PIXY_STATS_NO_FST=$(echo "$PIXY_STATS" | sed 's/,fst//g; s/fst,//g; s/^fst$//')
PIXY_STATS_HAS_FST=0
case ",${STATS}," in *,fst,*) PIXY_STATS_HAS_FST=1 ;; esac

# fst flavours: NA or comma-list of {wc, hudson}
FST_FLAVOURS=()
if (( PIXY_STATS_HAS_FST )); then
    if [[ "$FST_TYPE" == "NA" || -z "$FST_TYPE" ]]; then
        FST_FLAVOURS=(wc)
    else
        IFS=',' read -ra FST_FLAVOURS <<< "$FST_TYPE"
    fi
fi

# new_multi_fstbi: a third, fst-only variant that reads multiallelic sites in but
# excludes them from fst (--include_multiallelic_snps --fst_biallelic). this is the
# estimand the pre-2026-07 new_multi fst column held, when the fst path re-filtered
# to biallelic regardless of the flag. its pi/dxy would duplicate new_multi, so it
# emits none. only meaningful for hudson: WC is biallelic-only anyway (flag is a no-op).
RUN_FST_BIALLELIC=0
if (( PIXY_STATS_HAS_FST )); then
    for flav in "${FST_FLAVOURS[@]}"; do
        [[ "$flav" == "hudson" ]] && RUN_FST_BIALLELIC=1
    done
fi

##########
# header emission
##########
# per-arm header is the same across tasks; only task 1 emits it (to
# ${TASK_RESULTS}.header), which 03_concat_and_analyze.sh prepends to the partials.

# build the column-spec list once, reused per replicate.
# format: LABEL:FILE:KIND:VALUE_COL_CANDIDATES; FILE substituted per replicate.
build_column_specs() {
    local seed="$1"
    local rep_work="$2"
    local specs=()

    # pixy new / new_multi non-fst stats
    add_pixy_nofst() {
        local variant="$1"   # new or new_multi
        local subdir="$2"    # full directory
        local _stats
        IFS=',' read -ra _stats <<< "$PIXY_STATS_NO_FST"
        for s in "${_stats[@]}"; do
            [[ -z "$s" ]] && continue
            local label=""
            local cands=""
            case "$s" in
                pi)
                    label="pi_${variant}"
                    cands="avg_pi"
                    ;;
                dxy)
                    label="dxy_${variant}"
                    cands="avg_dxy"
                    ;;
                watterson_theta)
                    label="thetaw_${variant}"
                    cands="avg_wattersons_theta,avg_watterson_theta,avg_thetaw"
                    ;;
                tajima_d)
                    label="tajimaD_${variant}"
                    cands="tajima_d,avg_tajima_d,avg_tajimaD,avg_tajimas_d"
                    ;;
                *)
                    continue
                    ;;
            esac
            specs+=("${label}:${subdir}/rep${seed}_${s}.txt:pixy:${cands}")
        done
    }

    # pixy fst.
    # one fst flavour -> folded into the main pixy call (all 5 stats), output in
    # the main pixy_<variant> dir. two flavours -> one extra call per flavour in
    # its own ${base}_fst_<flav> subdir (two-call path).
    # for hudson, also capture per-window num/den so the analysis side can pool
    # ΣN/ΣD across windows (ratio-of-sums) not Jensen-biased mean-of-window-fsts.
    add_pixy_fst() {
        local variant="$1"
        local subdir_base="$2"   # ${rep_work}/pixy_${variant}
        if (( ${#FST_FLAVOURS[@]} == 1 )); then
            local flav="${FST_FLAVOURS[0]}"
            local label="fst_${flav}_${variant}"
            local cand="avg_${flav}_fst"
            local file="${subdir_base}/rep${seed}_fst.txt"
            specs+=("${label}:${file}:pixy:${cand}")
            if [[ "$flav" == "hudson" ]]; then
                specs+=("fst_hudson_num_${variant}:${file}:pixy:hudson_fst_num")
                specs+=("fst_hudson_den_${variant}:${file}:pixy:hudson_fst_den")
            fi
        else
            for flav in "${FST_FLAVOURS[@]}"; do
                local label="fst_${flav}_${variant}"
                local cand="avg_${flav}_fst"
                local file="${subdir_base}_fst_${flav}/rep${seed}_fst.txt"
                specs+=("${label}:${file}:pixy:${cand}")
                if [[ "$flav" == "hudson" ]]; then
                    specs+=("fst_hudson_num_${variant}:${file}:pixy:hudson_fst_num")
                    specs+=("fst_hudson_den_${variant}:${file}:pixy:hudson_fst_den")
                fi
            done
        fi
    }

    # the fst-only new_multi_fstbi variant. always one pixy call (hudson), output
    # directly in its own dir; the single/two-flavour subdir split does not apply.
    add_pixy_fst_biallelic() {
        local file="${rep_work}/pixy_new_multi_fstbi/rep${seed}_fst.txt"
        specs+=("fst_hudson_new_multi_fstbi:${file}:pixy:avg_hudson_fst")
        specs+=("fst_hudson_num_new_multi_fstbi:${file}:pixy:hudson_fst_num")
        specs+=("fst_hudson_den_new_multi_fstbi:${file}:pixy:hudson_fst_den")
    }

    add_pixy_nofst "new" "${rep_work}/pixy_new"
    add_pixy_nofst "new_multi" "${rep_work}/pixy_new_multi"
    (( PIXY_STATS_HAS_FST )) && add_pixy_fst "new" "${rep_work}/pixy_new"
    (( PIXY_STATS_HAS_FST )) && add_pixy_fst "new_multi" "${rep_work}/pixy_new_multi"
    (( RUN_FST_BIALLELIC )) && add_pixy_fst_biallelic

    if (( RUN_OLD_PIXY )); then
        # old pixy 0.95 emits pi / dxy / wc-fst (col avg_wc_fst)
        case "$STATS" in
            *pi*)
                specs+=("pi_old:${rep_work}/pixy_old/rep${seed}_pi.txt:pixy:avg_pi")
                ;;
        esac
        case "$STATS" in
            *dxy*)
                specs+=("dxy_old:${rep_work}/pixy_old/rep${seed}_dxy.txt:pixy:avg_dxy")
                ;;
        esac
        case "$STATS" in
            *fst*)
                specs+=("fst_wc_old:${rep_work}/pixy_old/rep${seed}_fst.txt:pixy:avg_wc_fst")
                ;;
        esac
    fi

    if (( RUN_VCFTOOLS )); then
        case "$STATS" in
            *pi*)
                specs+=("pi_vcftools:${rep_work}/vcftools/rep${seed}.windowed.pi:vcftools_pi:PI")
                ;;
        esac
        case "$STATS" in
            *fst*)
                specs+=("fst_wc_vcftools:${rep_work}/vcftools/rep${seed}.windowed.weir.fst:vcftools_fst:WEIGHTED_FST")
                ;;
        esac
    fi

    printf '%s\n' "${specs[@]}"
}

##########
# activate pixy env (once)
##########
# set +u around activation: conda's compiler activate.d scripts (binutils/gcc/gxx)
# read toolchain vars without defaults, tripping set -u. harmless on the cluster
# (no compiler envs); needed for local runs where pixy pulls them in.
set +u
eval "$(conda shell.bash hook)"
conda activate pixy
set -u

# emit header from task 1, using an arbitrary seed (header doesn't depend on seed)
if (( TASK_ID == 1 )); then
    HEADER_SPECS=$(build_column_specs "1" "${TASK_WORK}/rep_1")
    HEADER_ARGS=()
    while IFS= read -r spec; do
        [[ -z "$spec" ]] && continue
        HEADER_ARGS+=(--column-spec "$spec")
    done <<< "$HEADER_SPECS"
    python3 "$AGGREGATOR" --header-only "${HEADER_ARGS[@]}" > "${TASK_RESULTS}.header"
fi

##########
# per-replicate loop
##########

run_one_replicate() {
    local seed=$1
    local rep_work="${TASK_WORK}/rep_${seed}"
    mkdir -p "$rep_work"

    # simulate
    local vcf_path="${rep_work}/myvcf${seed}.vcf"
    if [[ "$PLOIDY_CHR2" == "NA" ]]; then
        # single-chromosome
        local sim_args=(
            --chromosome 1
            --replicates 1
            --seed "$seed"
            --sequence_length "$SEQ_LENGTH"
            --ploidy "$PLOIDY_CHR1"
            --Ne "$NE"
            --mu "$MU"
            --percent_missing_sites "$PCT_MISS_SITES"
            --percent_missing_genotypes "$PCT_MISS_GENO"
            --output_file "${rep_work}/myvcf"
            --sample_size "$SAMPLE_SIZE"
            --recombination_rate "$RECOMBINATION_RATE"
        )
        if (( N_POPULATIONS == 2 )); then
            sim_args+=(--population_mode 2 --div_time "$SPLIT_TIME")
        fi
        if ! conda run -n vcfsim vcfsim "${sim_args[@]}" >/dev/null 2>&1; then
            echo "[task ${TASK_ID} rep ${seed}] WARNING: vcfsim failed; skipping replicate" >&2
            rm -rf "$rep_work"; return 0
        fi
    else
        # mixed-ploidy across chromosomes 1 and 2
        local chrom_vcfs=()
        for chrom in 1 2; do
            local pvar="PLOIDY_CHR${chrom}"
            local ploidy=${!pvar}
            local prefix="${rep_work}/myvcf_chr${chrom}"
            if ! conda run -n vcfsim vcfsim \
                    --chromosome "$chrom" \
                    --replicates 1 \
                    --seed "$seed" \
                    --sequence_length "$SEQ_LENGTH" \
                    --ploidy "$ploidy" \
                    --Ne "$NE" \
                    --mu "$MU" \
                    --percent_missing_sites "$PCT_MISS_SITES" \
                    --percent_missing_genotypes "$PCT_MISS_GENO" \
                    --output_file "$prefix" \
                    --sample_size "$SAMPLE_SIZE" \
                    --recombination_rate "$RECOMBINATION_RATE" \
                    >/dev/null 2>&1; then
                echo "[task ${TASK_ID} rep ${seed}] WARNING: vcfsim chrom ${chrom} failed; skipping" >&2
                rm -rf "$rep_work"; return 0
            fi
            # vcfsim writes ${prefix}${seed}.vcf
            chrom_vcfs+=("${prefix}${seed}.vcf")
        done
        # bgzip + tabix each, then concat
        for v in "${chrom_vcfs[@]}"; do
            conda run -n vcftools bgzip -f "$v"
            conda run -n vcftools tabix -f "${v}.gz"
        done
        if ! conda run -n vcftools bcftools concat \
                -O v -o "$vcf_path" \
                "${chrom_vcfs[0]}.gz" "${chrom_vcfs[1]}.gz" \
                2>/dev/null; then
            echo "[task ${TASK_ID} rep ${seed}] WARNING: bcftools concat failed; skipping" >&2
            rm -rf "$rep_work"; return 0
        fi
        rm -f "${chrom_vcfs[0]}.gz" "${chrom_vcfs[0]}.gz.tbi"
        rm -f "${chrom_vcfs[1]}.gz" "${chrom_vcfs[1]}.gz.tbi"
    fi

    # compress + index
    if [[ -f "$vcf_path" ]]; then
        bgzip -f "$vcf_path"
    fi
    local vcf_gz="${vcf_path}.gz"
    if [[ ! -f "$vcf_gz" ]]; then
        echo "[task ${TASK_ID} rep ${seed}] WARNING: VCF not produced after sim; skipping" >&2
        rm -rf "$rep_work"; return 0
    fi
    tabix -f "$vcf_gz"

    # pixy 2.0: new (biallelic) and new_multi (--include_multiallelic_snps).
    # one pixy call covers all 5 stats (pi, dxy, fst, watterson_theta, tajima_d)
    # when the arm requests at most one fst flavour. pixy accepts one --fst_type
    # per call; if both wc and hudson are requested (no current arm does, but the
    # code stays generic) split: one call for non-fst stats + one per fst flavour.
    run_pixy_new() {
        local variant="$1"   # "new" or "new_multi"
        local out_base="${rep_work}/pixy_${variant}"
        local extra=()
        [[ "$variant" == "new_multi" ]] && extra=(--include_multiallelic_snps)

        if (( PIXY_STATS_HAS_FST )) && (( ${#FST_FLAVOURS[@]} == 1 )); then
            mkdir -p "$out_base"
            local flav="${FST_FLAVOURS[0]}"
            local stats_list="$PIXY_STATS_NO_FST"
            stats_list="${stats_list:+${stats_list},}fst"
            pixy --stats $(echo "$stats_list" | tr ',' ' ') \
                 --vcf "$vcf_gz" --populations "$POPS_FILE" \
                 --window_size "$WINDOW_SIZE" --n_cores "$N_CORES" \
                 --output_folder "$out_base" \
                 --output_prefix "rep${seed}" \
                 --fst_type "$flav" --fst_components \
                 "${extra[@]}" \
                 >/dev/null 2>&1 \
              || echo "[task ${TASK_ID} rep ${seed}] WARNING: pixy ${variant} (all stats) failed" >&2
            return
        fi

        if [[ -n "$PIXY_STATS_NO_FST" ]]; then
            mkdir -p "$out_base"
            pixy --stats $(echo "$PIXY_STATS_NO_FST" | tr ',' ' ') \
                 --vcf "$vcf_gz" --populations "$POPS_FILE" \
                 --window_size "$WINDOW_SIZE" --n_cores "$N_CORES" \
                 --output_folder "$out_base" \
                 --output_prefix "rep${seed}" \
                 "${extra[@]}" \
                 >/dev/null 2>&1 \
              || echo "[task ${TASK_ID} rep ${seed}] WARNING: pixy ${variant} (no-fst) failed" >&2
        fi

        if (( PIXY_STATS_HAS_FST )); then
            for flav in "${FST_FLAVOURS[@]}"; do
                local out_fst="${out_base}_fst_${flav}"
                mkdir -p "$out_fst"
                pixy --stats fst \
                     --vcf "$vcf_gz" --populations "$POPS_FILE" \
                     --window_size "$WINDOW_SIZE" --n_cores "$N_CORES" \
                     --output_folder "$out_fst" \
                     --output_prefix "rep${seed}" \
                     --fst_type "$flav" \
                     --fst_components \
                     "${extra[@]}" \
                     >/dev/null 2>&1 \
                  || echo "[task ${TASK_ID} rep ${seed}] WARNING: pixy ${variant} fst_${flav} failed" >&2
            done
        fi
    }

    # hudson fst over multiallelic-read-in but biallelic-filtered sites. always a
    # single fst-only call, so unlike run_pixy_new it needs no flavour-count split.
    run_pixy_fst_biallelic() {
        local out_base="${rep_work}/pixy_new_multi_fstbi"
        mkdir -p "$out_base"
        pixy --stats fst \
             --vcf "$vcf_gz" --populations "$POPS_FILE" \
             --window_size "$WINDOW_SIZE" --n_cores "$N_CORES" \
             --output_folder "$out_base" \
             --output_prefix "rep${seed}" \
             --fst_type hudson --fst_components \
             --include_multiallelic_snps --fst_biallelic \
             >/dev/null 2>&1 \
          || echo "[task ${TASK_ID} rep ${seed}] WARNING: pixy new_multi_fstbi failed" >&2
    }

    run_pixy_new new
    run_pixy_new new_multi
    (( RUN_FST_BIALLELIC )) && run_pixy_fst_biallelic

    # pixy 0.95.01 (old) — only for diploid arms with pi or dxy.
    # old pixy defaults --fst_maf_filter to 0.05, dropping MAF<=5% snps before
    # WC-FST and biasing it upward (~+0.02). force 0 so old WC-FST uses the same
    # site set as new pixy / vcftools (FigS2 agreement).
    if (( RUN_OLD_PIXY )); then
        mkdir -p "${rep_work}/pixy_old"
        local zarr_dir="${rep_work}/_zarr"
        local old_stats=""
        case "$STATS" in *pi*)  old_stats="${old_stats}${old_stats:+,}pi"  ;; esac
        case "$STATS" in *dxy*) old_stats="${old_stats}${old_stats:+,}dxy" ;; esac
        case "$STATS" in *fst*) old_stats="${old_stats}${old_stats:+,}fst" ;; esac
        if [[ -n "$old_stats" ]]; then
            for s in $(echo "$old_stats" | tr ',' ' '); do
                rm -rf "$zarr_dir"
                conda run -n old_pixy pixy --stats "$s" \
                    --vcf "$vcf_gz" \
                    --zarr_path "$zarr_dir" \
                    --populations "$POPS_FILE" \
                    --window_size "$WINDOW_SIZE" \
                    --chromosomes 1 \
                    --variant_filter_expression 'GT>=0' \
                    --invariant_filter_expression 'GT>=0' \
                    --bypass_filtration yes \
                    --fst_maf_filter 0 \
                    --outfile_prefix "${rep_work}/pixy_old/rep${seed}" \
                    >/dev/null 2>&1 \
                  || echo "[task ${TASK_ID} rep ${seed}] WARNING: old pixy ${s} failed" >&2
            done
            rm -rf "$zarr_dir"
        fi
    fi

    # vcftools
    if (( RUN_VCFTOOLS )); then
        mkdir -p "${rep_work}/vcftools"
        local vt_prefix="${rep_work}/vcftools/rep${seed}"
        case "$STATS" in
            *pi*)
                conda run -n vcftools vcftools \
                    --gzvcf "$vcf_gz" --window-pi "$WINDOW_SIZE" \
                    --out "$vt_prefix" \
                    >/dev/null 2>&1 \
                  || echo "[task ${TASK_ID} rep ${seed}] WARNING: vcftools --window-pi failed" >&2
                ;;
        esac
        case "$STATS" in
            *fst*)
                # vcftools needs separate sample-list files per pop
                local pop1_list="${rep_work}/vcftools/pop1.txt"
                local pop2_list="${rep_work}/vcftools/pop2.txt"
                awk -F'\t' '$2=="pop1" {print $1}' "$POPS_FILE" > "$pop1_list"
                awk -F'\t' '$2=="pop2" {print $1}' "$POPS_FILE" > "$pop2_list"
                conda run -n vcftools vcftools \
                    --gzvcf "$vcf_gz" \
                    --weir-fst-pop "$pop1_list" \
                    --weir-fst-pop "$pop2_list" \
                    --fst-window-size "$WINDOW_SIZE" \
                    --out "$vt_prefix" \
                    >/dev/null 2>&1 \
                  || echo "[task ${TASK_ID} rep ${seed}] WARNING: vcftools --weir-fst-pop failed" >&2
                ;;
        esac
    fi

    # aggregate
    local spec_args=()
    while IFS= read -r spec; do
        [[ -z "$spec" ]] && continue
        spec_args+=(--column-spec "$spec")
    done < <(build_column_specs "$seed" "$rep_work")
    python3 "$AGGREGATOR" --replicate "$seed" "${spec_args[@]}" \
        >> "$TASK_RESULTS" \
      || echo "[task ${TASK_ID} rep ${seed}] WARNING: aggregator failed" >&2

    # delete
    rm -rf "$rep_work"
}

REP_COUNT=0
for ((seed=TASK_SEED_START; seed<=TASK_SEED_END; seed++)); do
    run_one_replicate "$seed"
    REP_COUNT=$((REP_COUNT + 1))
done

echo "[$(date '+%F %T')] task ${TASK_ID} complete: ${REP_COUNT} replicates -> ${TASK_RESULTS}"
