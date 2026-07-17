#!/usr/bin/env Rscript
# dxy + Hudson FST multiallelic sweep across theta. mirrors 08_multiallelic_frac.R
# for the 2-pop arms dxy_{dip,tet,hex,oct}_2pop_theta{005,010,025,050,100}.
# -> dxy_fst_theta_summary.tsv, one row per (arm, ploidy, theta_nominal, stat,
#    estimator) for stat in {dxy, fst_hudson}.
#
# theoretical expectations:
#
# dxy (full JC69 with ancestral polymorphism; vcfsim doubles the ancestral pop to
# 2*Ne, hence the (8/3)*theta_nom denominator coefficient not (4/3)*theta_nom):
#     E[dxy] = (3/4) * (1 - exp(-8*mu*T/3) / (1 + (8/3) * theta_nom))
#
# Hudson FST (Bhatia 2013 ratio-of-averages, what pixy computes) has TWO estimands:
# the biallelic and multiallelic estimators target different ones. which applies
# depends on the estimator (see the per-estimator join in load_arm); attaching the
# wrong one is how the original mislabelling bug arose.
#
# 1. infinite-sites (Slatkin 1991) -> target for the BIALLELIC estimator. ratio of
#    coalescent times, ancestral pop Ne_anc = alpha * Ne:
#     tau   = T / (2 * Ne)
#     T_W   = 2 * Ne * (1 + (alpha - 1) * exp(-tau))
#     T_B   = 2 * Ne * (tau + alpha)
#     E[FST] = 1 - T_W / T_B = (tau + (alpha-1)(1 - exp(-tau))) / (tau + alpha)
#    vcfsim doubles the ancestral pop, so alpha = 2. the naive 1 - exp(-T/(2*Ne))
#    is alpha = 0 and over-predicts FST by ~6% absolute.
#    theta-free, and the biallelic estimator tracks it within 0.2% -- because the
#    biallelic filter conditions away visible homoplasy, not because mutation is
#    absent. the filter is imperfect: back/parallel mutation keeps some sites
#    biallelic, so the estimator sags -0.2 to -0.4% at theta = 0.1 (significant at
#    n = 2000, growing with ploidy).
#
# 2. finite-sites (Jukes & Cantor 1969) -> target for the MULTIALLELIC estimator,
#    which measures observed dissimilarity and saturates. both terms (3/4)(1 - L(beta)),
#    L = Laplace transform of the pairwise coalescence time:
#     E[FST]_finite = 1 - E[pi_w] / E[dxy]
#    only mildly theta-dependent (-1.6% at theta = 0.1); the saturation largely
#    cancels between numerator and denominator. derivation + 4 checks:
#    analysis/finite_sites_fst_derivation.md. validated at n = 2000 within one
#    bootstrap SE at every theta x ploidy.
#
# theoretical_dxy needs no split -- already finite-sites, and dxy has one estimand.
# also writes figs/11_dxy_fst_theta.pdf as a sanity-check plot.

suppressPackageStartupMessages({
    library(tidyverse)
})

fig_dir <- "figs"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

sim_params <- read_tsv("../config/sim_params.tsv", show_col_types = FALSE)

# the theta sweep ran at 2n/4n/6n/8n, but the multiallelic analysis (Figure 4)
# reports 2n and 8n only (decision 2026-07-15). the 4n/6n arms are simulated and
# aggregated, just not carried here. Figure 3's ploidy x missingness grid is a
# different arm set and still uses all four ploidies.
FIG4_PLOIDIES <- c("diploid", "octoploid")

theta_arms <- sim_params |>
    filter(str_detect(arm, "^dxy_.*_theta")) |>
    mutate(
        ploidy = case_when(
            str_detect(arm, "_dip_") ~ "diploid",
            str_detect(arm, "_tet_") ~ "tetraploid",
            str_detect(arm, "_hex_") ~ "hexaploid",
            str_detect(arm, "_oct_") ~ "octoploid"
        ),
        theta_nominal = case_when(
            str_detect(arm, "theta005") ~ 0.005,
            str_detect(arm, "theta010") ~ 0.010,
            str_detect(arm, "theta025") ~ 0.025,
            str_detect(arm, "theta050") ~ 0.050,
            str_detect(arm, "theta100") ~ 0.100
        ),
        theta_nom        = 4 * Ne * mu,
        theoretical_dxy  = (3/4) * (1 - exp(-8 * mu * split_time / 3) /
                                        (1 + (8/3) * theta_nom)),
        fst_tau          = split_time / (2 * Ne),
        fst_alpha        = 2,
        # infinite-sites (Slatkin 1991) -> BIALLELIC estimator
        theoretical_fst_infinite = (fst_tau + (fst_alpha - 1) *
                                      (1 - exp(-fst_tau))) /
                                   (fst_tau + fst_alpha),
        # finite-sites (JC69) -> MULTIALLELIC estimator
        # as theta -> 0 this converges to theoretical_fst_infinite
        fs_k             = 1 + (4/3) * theta_nom,
        fs_L_dxy         = exp(-(4/3) * theta_nom * fst_tau) /
                           (1 + fst_alpha * (4/3) * theta_nom),
        fs_L_pi_w        = (1 - exp(-fst_tau * fs_k)) / fs_k +
                           exp(-fst_tau * fs_k) /
                           (1 + fst_alpha * (4/3) * theta_nom),
        theoretical_fst_finite   = 1 - (1 - fs_L_pi_w) / (1 - fs_L_dxy)
    ) |>
    filter(ploidy %in% FIG4_PLOIDIES)

load_arm <- function(arm_row) {
    f <- file.path("../data/aggregated", paste0(arm_row$arm, ".tsv"))
    if (!file.exists(f)) {
        warning("missing: ", f)
        return(tibble())
    }
    df <- read_tsv(f, show_col_types = FALSE)

    dxy_cols <- intersect(names(df),
                          c("dxy_new", "dxy_new_multi"))
    fst_cols <- intersect(names(df),
                          c("fst_hudson_new", "fst_hudson_new_multi"))

    out <- list()

    # dxy: per-window values, mean()'d later
    if (length(dxy_cols) > 0) {
        out$dxy <- df |>
            select(replicate, all_of(dxy_cols)) |>
            pivot_longer(all_of(dxy_cols),
                         names_to = "estimator", values_to = "value") |>
            mutate(
                estimator     = str_remove(estimator, "^dxy_"),
                stat          = "dxy",
                arm           = arm_row$arm,
                ploidy        = arm_row$ploidy,
                theta_nominal = arm_row$theta_nominal,
                theoretical   = arm_row$theoretical_dxy
            )
    }

    # fst: pool ΣN/ΣD within each replicate, not mean-of-windows.
    # pixy's per-window avg_hudson_fst is itself ΣN_sites/ΣD_sites within the
    # window (unbiased); averaging per-window FSTs across windows is Jensen-biased
    # downward. --fst_components gives hudson_fst_num/den, so per-rep pooled FST is
    # the exact ratio-of-sums pooled_r = Σ_w num_w / Σ_w den_w. carrying
    # num_total/den_total per rep lets the summary further pool across reps for a
    # Jensen-unbiased grand mean: Σ_r Σ_w num / Σ_r Σ_w den.
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
            summarise(
                num_total = sum(num),
                den_total = sum(den),
                value     = num_total / den_total,
                .groups   = "drop"
            ) |>
            mutate(estimator = estimator_suffix)
    }

    fst_pooled <- bind_rows(
        pool_fst("new"),
        pool_fst("new_multi"),
        pool_fst("new_multi_fstbi")
    )

    # expectation is PER ESTIMATOR, not per arm -- the two target different
    # estimands (see header). new_multi_fstbi reads multiallelic sites in but
    # filters FST back to biallelic, so it targets the coalescent estimand like
    # new does; reproduces the pre-2026-07 new_multi column exactly (n = 2000).
    if (nrow(fst_pooled) > 0) {
        out$fst <- fst_pooled |>
            mutate(
                stat          = "fst_hudson",
                arm           = arm_row$arm,
                ploidy        = arm_row$ploidy,
                theta_nominal = arm_row$theta_nominal,
                theoretical   = if_else(
                    estimator == "new_multi",
                    arm_row$theoretical_fst_finite,
                    arm_row$theoretical_fst_infinite
                )
            )
    }

    bind_rows(out)
}

all <- bind_rows(lapply(seq_len(nrow(theta_arms)),
                        function(i) load_arm(theta_arms[i, ])))

if (nrow(all) == 0) {
    cat("No Phase G ext aggregated TSVs found yet. Run the dxy_*_theta arms first.\n")
    quit(status = 0)
}

# dxy summary: mean of per-window values across all reps.
# fst summary: per-rep value is the pooled ratio-of-sums within a rep; the central
# estimate further pools across reps via Σ num_total / Σ den_total (Jensen-unbiased
# grand mean). per-rep value still feeds q025/q975 for the rep-to-rep CI.
mean_dxy <- all |>
    filter(stat == "dxy") |>
    group_by(arm, ploidy, theta_nominal, stat, estimator) |>
    summarise(mean = mean(value, na.rm = TRUE), .groups = "drop")

mean_fst <- all |>
    filter(stat == "fst_hudson") |>
    group_by(arm, ploidy, theta_nominal, stat, estimator) |>
    summarise(mean = sum(num_total, na.rm = TRUE) /
                     sum(den_total, na.rm = TRUE),
              .groups = "drop")

summary_tbl <- all |>
    group_by(arm, ploidy, theta_nominal, stat, estimator) |>
    summarise(
        q025        = quantile(value, 0.025, na.rm = TRUE),
        q975        = quantile(value, 0.975, na.rm = TRUE),
        theoretical = first(theoretical),
        .groups     = "drop"
    ) |>
    left_join(bind_rows(mean_dxy, mean_fst),
              by = c("arm", "ploidy", "theta_nominal", "stat", "estimator")) |>
    relocate(mean, .before = q025) |>
    mutate(
        bias     = mean - theoretical,
        bias_pct = 100 * bias / theoretical
    )

write_tsv(summary_tbl, "dxy_fst_theta_summary.tsv")

p <- summary_tbl |>
    ggplot(aes(theta_nominal, bias_pct,
               color = estimator, linetype = ploidy, shape = ploidy)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
    geom_line() +
    geom_pointrange(aes(ymin = 100 * (q025 - theoretical) / theoretical,
                        ymax = 100 * (q975 - theoretical) / theoretical),
                    position = position_dodge(width = 0.005), size = 0.3) +
    facet_wrap(~ stat, scales = "free_y") +
    scale_x_log10() +
    scale_color_brewer(palette = "Dark2") +
    labs(x = expression(theta == 4 * N[e] * mu ~ "(log scale)"),
         y = "bias relative to theory (%)",
         title = "Bias of dxy and Hudson FST vs theta",
         subtitle = "biallelic-only `new` vs multiallelic-aware `new_multi`",
         color = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position = "top")

ggsave(file.path(fig_dir, "11_dxy_fst_theta.pdf"), p, width = 9, height = 5)

cat("Wrote analysis/dxy_fst_theta_summary.tsv and figs/11_dxy_fst_theta.pdf\n")
