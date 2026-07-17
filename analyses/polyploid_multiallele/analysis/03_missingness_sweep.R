#!/usr/bin/env Rscript
# missingness x ploidy sweep -> missingness_summary.tsv + figs/03_missingness_sweep.pdf
# one row per (stat, ploidy, miss_pct, estimator) for stat in {pi, dxy, fst_hudson, fst_wc}
# fst_wc (new/old/vcftools) feeds the diploid old-vs-new Figure S2
#
# fst_hudson pools --fst_components: per-rep value = Σ_w num_w / Σ_w den_w, then a
# grand pool Σ_r num_total / Σ_r den_total for the central estimate. see
# 11_dxy_fst_theta_sweep.R for the derivation.
# hudson theoretical: structured coalescent, alpha = Ne_anc/Ne = 2 (vcfsim ancestor)
#   tau = T / (2 Ne);  E[FST] = (tau + (alpha-1)(1 - exp(-tau))) / (tau + alpha)

suppressPackageStartupMessages({
    library(tidyverse)
})

fig_dir <- "figs"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

sim_params <- read_tsv("../config/sim_params.tsv", show_col_types = FALSE)

# missingness arms + their miss=0 references (2n baseline, 6n/8n, 2n/4n anchors)
phase_b_arms <- sim_params |>
    filter(str_detect(arm, "_miss") | arm %in% c(
        "pi_dip_1pop", "pi_tet_1pop", "pi_hex_1pop", "pi_oct_1pop",
        "dxy_dip_2pop", "dxy_tet_2pop", "dxy_hex_2pop", "dxy_oct_2pop",
        "pi_dip_1pop_tw", "pi_tet_1pop_tw"
    )) |>
    mutate(
        miss_pct = case_when(
            str_detect(arm, "miss10") ~ 0.10,
            str_detect(arm, "miss25") ~ 0.25,
            str_detect(arm, "miss50") ~ 0.50,
            str_detect(arm, "miss75") ~ 0.75,
            TRUE                      ~ 0.00
        ),
        theta_nom = 4 * Ne * mu,
        tau       = if_else(!is.na(split_time), split_time / (2 * Ne), NA_real_),
        ploidy    = ploidy_chr1
    )

# theoretical expectations per (stat, demography)
theoretical_for <- function(stat, theta_nom, mu, split_time, tau, alpha = 2) {
    switch(stat,
        pi  = theta_nom / (1 + (4/3) * theta_nom),
        dxy = (3/4) * (1 - exp(-8 * mu * split_time / 3) /
                              (1 + (8/3) * theta_nom)),
        fst_hudson = (tau + (alpha - 1) * (1 - exp(-tau))) / (tau + alpha),
        fst_wc     = (tau + (alpha - 1) * (1 - exp(-tau))) / (tau + alpha),
        NA_real_
    )
}

# per-arm loader; returns rows for each (stat, estimator) present in the arm
load_arm <- function(arm_row) {
    f <- file.path("../data/aggregated", paste0(arm_row$arm, ".tsv"))
    if (!file.exists(f)) {
        warning("missing: ", f)
        return(tibble())
    }
    df <- read_tsv(f, show_col_types = FALSE)
    stats_in_arm <- str_split(arm_row$stats, ",")[[1]]
    out <- list()

    # pi: window-mean across reps (mean-of-window-pis unbiased)
    if ("pi" %in% stats_in_arm) {
        pi_cols <- intersect(names(df),
                             c("pi_old", "pi_new", "pi_new_multi", "pi_vcftools"))
        if (length(pi_cols) > 0) {
            out$pi <- df |>
                select(replicate, all_of(pi_cols)) |>
                pivot_longer(all_of(pi_cols),
                             names_to = "estimator", values_to = "value") |>
                mutate(estimator = str_remove(estimator, "^pi_"),
                       stat = "pi",
                       theoretical = theoretical_for("pi",
                                                     arm_row$theta_nom,
                                                     arm_row$mu,
                                                     arm_row$split_time,
                                                     arm_row$tau))
        }
    }

    # dxy: window-mean across reps
    if ("dxy" %in% stats_in_arm) {
        dxy_cols <- intersect(names(df),
                              c("dxy_old", "dxy_new", "dxy_new_multi"))
        if (length(dxy_cols) > 0) {
            out$dxy <- df |>
                select(replicate, all_of(dxy_cols)) |>
                pivot_longer(all_of(dxy_cols),
                             names_to = "estimator", values_to = "value") |>
                mutate(estimator = str_remove(estimator, "^dxy_"),
                       stat = "dxy",
                       theoretical = theoretical_for("dxy",
                                                     arm_row$theta_nom,
                                                     arm_row$mu,
                                                     arm_row$split_time,
                                                     arm_row$tau))
        }
    }

    # fst_hudson: components-pool per rep, grand-pool across reps
    if ("fst" %in% stats_in_arm) {
        pool_fst <- function(estimator_suffix) {
            num_col <- paste0("fst_hudson_num_", estimator_suffix)
            den_col <- paste0("fst_hudson_den_", estimator_suffix)
            if (!(num_col %in% names(df) && den_col %in% names(df))) return(NULL)
            df |>
                select(replicate,
                       num = all_of(num_col),
                       den = all_of(den_col)) |>
                filter(!is.na(num), !is.na(den)) |>
                group_by(replicate) |>
                summarise(num_total = sum(num),
                          den_total = sum(den),
                          value     = num_total / den_total,
                          .groups   = "drop") |>
                mutate(estimator = estimator_suffix)
        }
        fst_pooled <- bind_rows(pool_fst("new"), pool_fst("new_multi"))
        if (nrow(fst_pooled) > 0) {
            out$fst_hudson <- fst_pooled |>
                mutate(stat = "fst_hudson",
                       theoretical = theoretical_for("fst_hudson",
                                                     arm_row$theta_nom,
                                                     arm_row$mu,
                                                     arm_row$split_time,
                                                     arm_row$tau))
        }
    }

    # fst_wc: per-window Weir-Cockerham FST (new/new_multi/old/vcftools).
    # old pixy and vcftools emit only per-window WC FST (no num/den), so pool
    # mean-of-windows for ALL wc estimators for a fair cross-method comparison
    # (Figure S2) -- differs from the hudson ratio-of-sums.
    if ("fst" %in% stats_in_arm) {
        wc_cols <- intersect(names(df),
                             c("fst_wc_new", "fst_wc_new_multi",
                               "fst_wc_old", "fst_wc_vcftools"))
        if (length(wc_cols) > 0) {
            out$fst_wc <- df |>
                select(replicate, all_of(wc_cols)) |>
                pivot_longer(all_of(wc_cols),
                             names_to = "estimator", values_to = "value") |>
                mutate(estimator = str_remove(estimator, "^fst_wc_"),
                       stat = "fst_wc",
                       theoretical = theoretical_for("fst_wc",
                                                     arm_row$theta_nom,
                                                     arm_row$mu,
                                                     arm_row$split_time,
                                                     arm_row$tau))
        }
    }

    bind_rows(out) |>
        mutate(arm = arm_row$arm,
               ploidy = arm_row$ploidy,
               miss_pct = arm_row$miss_pct)
}

all <- bind_rows(lapply(seq_len(nrow(phase_b_arms)),
                        function(i) load_arm(phase_b_arms[i, ])))

if (nrow(all) == 0) {
    cat("No Phase B aggregated TSVs found yet. Run ./01_run_all.sh --include miss first.\n")
    quit(status = 0)
}

# per-rep values for the Figure 3 violins: mean-of-windows per rep for
# pi/dxy/fst_wc; fst_hudson is already one pooled value per rep (per-rep mean is
# a no-op). same pooling the summary quantiles use.
per_rep <- all |>
    group_by(stat, ploidy, miss_pct, estimator, replicate) |>
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
write_tsv(per_rep, "per_rep_missingness.tsv")

# per-(stat, ploidy, miss_pct, estimator) summary.
# pi/dxy: mean-of-window-values across reps; q025/q975 from those values too.
# fst_hudson: grand-pool ratio-of-sums for mean; q025/q975 from per-rep pooled values.
mean_pi_dxy <- all |>
    filter(stat %in% c("pi", "dxy", "fst_wc")) |>
    group_by(stat, ploidy, miss_pct, estimator) |>
    summarise(mean = mean(value, na.rm = TRUE),
              n_windows = n(),
              .groups = "drop")

mean_fst <- all |>
    filter(stat == "fst_hudson") |>
    group_by(stat, ploidy, miss_pct, estimator) |>
    summarise(mean = sum(num_total, na.rm = TRUE) /
                     sum(den_total, na.rm = TRUE),
              n_windows = n(),
              .groups = "drop")

summary_tbl <- all |>
    group_by(stat, ploidy, miss_pct, estimator) |>
    summarise(q025        = quantile(value, 0.025, na.rm = TRUE),
              q975        = quantile(value, 0.975, na.rm = TRUE),
              sd_value    = sd(value, na.rm = TRUE),
              n_reps      = sum(!is.na(value)),
              theoretical = first(theoretical),
              .groups     = "drop") |>
    left_join(bind_rows(mean_pi_dxy, mean_fst),
              by = c("stat", "ploidy", "miss_pct", "estimator")) |>
    mutate(se_mean = sd_value / sqrt(n_reps),
           ci_lo   = mean - 1.96 * se_mean,
           ci_hi   = mean + 1.96 * se_mean) |>
    relocate(mean, n_windows, .before = q025)

write_tsv(summary_tbl, "missingness_summary.tsv")

p <- summary_tbl |>
    ggplot(aes(miss_pct, mean, color = estimator)) +
    geom_line() +
    geom_pointrange(aes(ymin = q025, ymax = q975),
                    position = position_dodge(width = 0.02),
                    size = 0.3) +
    geom_hline(aes(yintercept = theoretical),
               linetype = "dashed", linewidth = 0.4) +
    facet_grid(stat ~ ploidy, scales = "free_y",
               labeller = labeller(ploidy = function(x) paste0("ploidy ", x))) +
    scale_x_continuous(labels = scales::percent) +
    scale_color_brewer(palette = "Dark2") +
    labs(x = "missingness fraction",
         y = "estimate (mean +/- 95% interval)",
         title = "Effect of missing data on pi, dxy, and Hudson FST by ploidy",
         subtitle = "Dashed = theoretical expectation",
         color = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position = "top")

ggsave(file.path(fig_dir, "03_missingness_sweep.pdf"), p,
       width = 9, height = 7)

cat("Wrote analysis/missingness_summary.tsv and figs/03_missingness_sweep.pdf\n")
