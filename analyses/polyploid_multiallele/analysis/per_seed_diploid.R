#!/usr/bin/env Rscript
# per-vcf (per-seed) values for Figure S2: collapse each replicate's windows to
# one value per (stat, estimator, miss_pct, replicate) via mean-of-windows -- the
# same per-rep pooling 03_missingness_sweep.R uses for pi/dxy/WC-FST. diploid only.
#   pi     from pi_dip_1pop[_miss*]  : new, old, vcftools
#   dxy    from dxy_dip_2pop[_miss*] : new, old
#   fst_wc from dxy_dip_2pop[_miss*] : new (band-aid), old, vcftools
# -> analysis/per_seed_diploid.tsv (stat, estimator, miss_pct, replicate, value)
# feeds the per-seed scatter matrix in Figure_S2_Diploid_OldNew.R

suppressPackageStartupMessages({
    library(tidyverse)
})

AGG <- "../data/aggregated"

miss_of <- function(arm) {
    case_when(
        str_detect(arm, "miss10") ~ 0.10,
        str_detect(arm, "miss25") ~ 0.25,
        str_detect(arm, "miss50") ~ 0.50,
        str_detect(arm, "miss75") ~ 0.75,
        TRUE                      ~ 0.00
    )
}

# per-rep mean-of-windows for a set of estimator columns in one arm's table
per_seed_arm <- function(arm, stat, cols, prefix) {
    f <- file.path(AGG, paste0(arm, ".tsv"))
    if (!file.exists(f)) { warning("missing: ", f); return(tibble()) }
    df <- read_tsv(f, show_col_types = FALSE)
    have <- intersect(cols, names(df))
    if (length(have) == 0) { warning("no ", stat, " cols in ", arm); return(tibble()) }
    df |>
        select(replicate, all_of(have)) |>
        pivot_longer(all_of(have), names_to = "estimator", values_to = "value") |>
        mutate(estimator = str_remove(estimator, prefix)) |>
        group_by(estimator, replicate) |>
        summarise(value = mean(value, na.rm = TRUE), .groups = "drop") |>
        mutate(stat = stat, miss_pct = miss_of(arm))
}

miss_suffix <- c("", "_miss10", "_miss25", "_miss50", "_miss75")
pi_arms  <- paste0("pi_dip_1pop",  miss_suffix)
dxy_arms <- paste0("dxy_dip_2pop", miss_suffix)

rows <- bind_rows(
    map(pi_arms,  ~per_seed_arm(.x, "pi",
                                c("pi_new", "pi_old", "pi_vcftools"), "^pi_")),
    map(dxy_arms, ~per_seed_arm(.x, "dxy",
                                c("dxy_new", "dxy_old"), "^dxy_")),
    map(dxy_arms, ~per_seed_arm(.x, "fst_wc",
                                c("fst_wc_new", "fst_wc_old", "fst_wc_vcftools"),
                                "^fst_wc_"))
)

if (nrow(rows) == 0) stop("No per-seed rows produced — check aggregated tables.")

rows <- rows |>
    select(stat, estimator, miss_pct, replicate, value) |>
    arrange(stat, estimator, miss_pct, replicate)

write_tsv(rows, "per_seed_diploid.tsv")

cat("Wrote analysis/per_seed_diploid.tsv\n")
rows |>
    count(stat, estimator, miss_pct) |>
    print(n = Inf)
