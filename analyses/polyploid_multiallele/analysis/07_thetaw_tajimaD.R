#!/usr/bin/env Rscript
# polyploid extension of the missingness-aware Watterson's theta + Tajima's D
# (Bailey, Stevison & Samuk 2025 validated these on diploids). across every 1-pop
# arm with thetaw + tajimaD in stats (missingness sweep + ploidy grid), plot:
#   - mean(thetaw_new_multi) vs theoretical 4*Ne*mu
#   - mean(tajimaD_new_multi) vs the neutral reference
# over ploidy x missingness.
#
# inputs: ../data/aggregated/{pi_dip_1pop_miss*, pi_tet_1pop, ...}.tsv with cols
#         thetaw_new, thetaw_new_multi, tajimaD_new, tajimaD_new_multi.
# -> figs/07_thetaw_tajimaD.pdf

suppressPackageStartupMessages({
    library(tidyverse)
})

fig_dir <- "figs"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

sim_params <- read_tsv("../config/sim_params.tsv", show_col_types = FALSE)

# Finite-sites JC69 expectation for Watterson's theta on n sampled gene copies.
#   E[theta_W] = (1 - P(monomorphic)) / a_n
#   P(monomorphic) = 4 * (theta/3)^{(n)} / (4*theta/3)^{(n)}
# where x^{(n)} is the rising factorial gamma(x+n)/gamma(x). Uses lgamma for
# numerical stability across the full ploidy x sample-size grid.
e_thetaW_jc69 <- function(theta, n) {
    if (n < 2) return(NA_real_)
    a_n <- sum(1 / seq_len(n - 1))
    log_p_mono <- log(4) +
        (lgamma(theta / 3 + n)     - lgamma(theta / 3)) -
        (lgamma(4 * theta / 3 + n) - lgamma(4 * theta / 3))
    (1 - exp(log_p_mono)) / a_n
}

# Finite-sites JC69 expectation for the MUTATION-count (eta = sum k-1) estimator,
# E(s*)/a_n, from Tajima (1996, Genetics 143:1457) Eq 15: E(s*) = 3 - 4*p_ijk,
# where p_ijk = P(a particular nucleotide is absent) = (theta)^{(n)} / (4*theta/3)^{(n)}.
# This is the analytic target for pixy's CORRECTED multiallelic theta_W (branch
# multiallelic-mutation-count-theta-d); the site-count E(s)/a_n above is the target
# for the old site-count theta_W and the biallelic estimator.
e_thetaStar_jc69 <- function(theta, n) {
    if (n < 2) return(NA_real_)
    a_n <- sum(1 / seq_len(n - 1))
    log_p_ijk <- (lgamma(4 * theta / 3) + lgamma(theta + n)) -
                 (lgamma(4 * theta / 3 + n) + lgamma(theta))
    (3 - 4 * exp(log_p_ijk)) / a_n
}

# Calibrated finite-sites neutral E[D] per (ploidy x missingness) cell (finite-
# sites JC69 coalescent matched to the vcfsim window, with vcfsim's missing-data
# masking replayed on each window; see calibrate_neutral_tajimaD.py).
# Tajima's D is an n-normalised per-window contrast, so E[D] != 0 under finite
# sites AND moves with missingness: genotype-missing lowers the effective per-site
# sample size, and site-missing thins each window. The reference therefore varies
# by missingness, not just ploidy -- joined on (ploidy, miss_pct) below.
neutral_ref <- read_tsv("neutral_tajimaD_reference.tsv", show_col_types = FALSE) |>
    select(ploidy, miss_pct, theoretical_tajimaD = mean_D, n_chrom_eff)

# 1-pop arms with thetaw + tajimaD in their stats
tw_arms <- sim_params |>
    filter(n_populations == 1, str_detect(stats, "thetaw"), str_detect(stats, "tajimaD")) |>
    mutate(
        ploidy = ploidy_chr1,
        miss_pct = case_when(
            str_detect(arm, "miss10") ~ 0.10,
            str_detect(arm, "miss25") ~ 0.25,
            str_detect(arm, "miss50") ~ 0.50,
            str_detect(arm, "miss75") ~ 0.75,
            TRUE                      ~ 0.00
        ),
        theta_nom = 4 * Ne * mu,
        # n = number of sampled gene copies = diploid-equivalent samples x ploidy
        n_chrom = sample_size * ploidy
    ) |>
    left_join(neutral_ref, by = c("ploidy", "miss_pct")) |>
    mutate(
        # vcfsim genotype-missingness drops whole individuals, lowering the
        # effective per-site sample size to a constant n_chrom_eff (computed by the
        # calibration). Watterson's theta = eta / a_1(n), so its finite-sites JC69
        # expectation must be evaluated at n_eff, NOT the full n_chrom -- that is
        # what makes the theta_W reference slope with missingness (matching the
        # data) instead of sitting flat. (site-missingness only thins the per-window
        # site count, leaving the per-site expectation unchanged.) Falls back to the
        # full n_chrom for any cell not present in the calibration reference.
        n_eff = dplyr::coalesce(n_chrom_eff, n_chrom),
        theoretical_thetaw    = mapply(e_thetaW_jc69, theta_nom, n_eff),
        theoretical_thetaStar = mapply(e_thetaStar_jc69, theta_nom, n_eff)
    )

load_arm <- function(arm_row) {
    f <- file.path("../data/aggregated", paste0(arm_row$arm, ".tsv"))
    if (!file.exists(f)) {
        warning("missing: ", f)
        return(tibble())
    }
    df <- read_tsv(f, show_col_types = FALSE)
    cols <- intersect(names(df), c(
        "thetaw_new", "thetaw_new_multi",
        "tajimaD_new", "tajimaD_new_multi"
    ))
    if (length(cols) == 0) return(tibble())
    df |>
        select(replicate, all_of(cols)) |>
        pivot_longer(all_of(cols), names_to = "col", values_to = "value") |>
        mutate(
            stat      = if_else(str_starts(col, "thetaw"), "thetaw", "tajimaD"),
            estimator = str_remove(col, "^(thetaw|tajimaD)_"),
            arm       = arm_row$arm,
            ploidy    = arm_row$ploidy,
            miss_pct  = arm_row$miss_pct,
            # theta_W theoretical is estimator-specific: the corrected multiallelic
            # estimator (new_multi) counts mutations -> E(s*)/a_n (Tajima 1996 Eq 15);
            # the biallelic estimator (new) counts sites -> E(s)/a_n (Eq 6). Tajima's
            # D reference is the calibrated neutral mean (~0 here). This column is
            # consumed by figures/Figure_3_Polyploid.R.
            theoretical = case_when(
                !str_starts(col, "thetaw") ~ arm_row$theoretical_tajimaD,
                str_detect(col, "_multi")  ~ arm_row$theoretical_thetaStar,
                TRUE                       ~ arm_row$theoretical_thetaw
            )
        )
}

all <- bind_rows(lapply(seq_len(nrow(tw_arms)), function(i) load_arm(tw_arms[i, ])))

if (nrow(all) == 0) {
    cat("No Phase F-relevant aggregated TSVs found yet.\n")
    quit(status = 0)
}

# Per-replicate values (mean-of-windows per rep) for the Figure 3 violins.
per_rep <- all |>
    group_by(stat, ploidy, miss_pct, estimator, replicate) |>
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
write_tsv(per_rep, "per_rep_thetaw_tajimaD.tsv")

summary_tbl <- all |>
    group_by(stat, ploidy, miss_pct, estimator) |>
    summarise(
        mean = mean(value, na.rm = TRUE),
        q025 = quantile(value, 0.025, na.rm = TRUE),
        q975 = quantile(value, 0.975, na.rm = TRUE),
        sd_value = sd(value, na.rm = TRUE),
        n_reps   = sum(!is.na(value)),
        theoretical = first(theoretical),
        .groups = "drop"
    ) |>
    mutate(se_mean = sd_value / sqrt(n_reps),
           ci_lo   = mean - 1.96 * se_mean,
           ci_hi   = mean + 1.96 * se_mean)

write_tsv(summary_tbl, "thetaw_tajimaD_summary.tsv")

p <- summary_tbl |>
    ggplot(aes(miss_pct, mean, color = factor(ploidy), shape = estimator)) +
    geom_line(aes(group = interaction(ploidy, estimator))) +
    geom_pointrange(aes(ymin = q025, ymax = q975),
                    position = position_dodge(width = 0.02), size = 0.3) +
    geom_hline(aes(yintercept = theoretical),
               linetype = "dashed", linewidth = 0.4) +
    facet_wrap(~ stat, scales = "free_y", labeller = label_both) +
    scale_x_continuous(labels = scales::percent) +
    scale_color_viridis_d() +
    labs(x = "missingness fraction",
         y = "estimate (mean +/- 95%)",
         title = "Watterson's theta and Tajima's D on polyploid VCFs",
         subtitle = paste(
             "Dashed = finite-sites JC69 expectation (Tajima 1996): E(s*)/a_n for",
             "multiallelic, E(s)/a_n for biallelic; D ~ 0 under the null"),
         color = "ploidy", shape = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position = "top")

ggsave(file.path(fig_dir, "07_thetaw_tajimaD.pdf"), p, width = 9, height = 6)

cat("Wrote analysis/thetaw_tajimaD_summary.tsv and figs/07_thetaw_tajimaD.pdf\n")
