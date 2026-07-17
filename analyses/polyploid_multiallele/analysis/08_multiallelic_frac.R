#!/usr/bin/env Rscript
# multiallelic-fraction sweep across theta: bias of biallelic-only pi vs
# multiallelic-aware pi as theta = 4*Ne*mu grows (more multiallelic sites).
# relates to Sopniewski et al. 2024 (excluding multiallelic sites biases het).
#
# inputs: ../data/aggregated/pi_dip_1pop_theta{005,010,025,050,100}.tsv
# -> figs/08_multiallelic_frac.pdf

suppressPackageStartupMessages({
    library(tidyverse)
})

fig_dir <- "figs"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

sim_params <- read_tsv("../config/sim_params.tsv", show_col_types = FALSE)

# the theta sweep ran at 2n/4n/6n/8n, but the multiallelic analysis (Figure 4)
# reports 2n and 8n only (decision 2026-07-15). the 4n/6n arms are simulated and
# aggregated, just not carried here. Figure 3's ploidy x missingness grid is a
# different arm set (pi_*_miss*, dxy_*_miss*) and still uses all four ploidies.
FIG4_PLOIDIES <- c("diploid", "octoploid")

# match the pi 1-pop theta arms only. "_theta" alone also matches dxy_*_2pop_theta*,
# which since 2026-07-14 carry pi in --stats -- but that pi is one population's
# diversity within the 2-pop split (E[pi_w] = 0.138 at theta = 0.1), not the
# panmictic pi this sweep is about (E[pi] = 0.088). the loose filter would double
# every cell and put two points per theta in Figure 4's pi row. keep the anchor.
theta_arms <- sim_params |>
    filter(str_detect(arm, "^pi_.*_1pop_theta")) |>
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
        theta_nom      = 4 * Ne * mu,
        theoretical_pi = theta_nom / (1 + (4/3) * theta_nom)
    ) |>
    filter(ploidy %in% FIG4_PLOIDIES)

load_arm <- function(arm_row) {
    f <- file.path("../data/aggregated", paste0(arm_row$arm, ".tsv"))
    if (!file.exists(f)) {
        warning("missing: ", f)
        return(tibble())
    }
    df <- read_tsv(f, show_col_types = FALSE)
    cols <- intersect(names(df), c("pi_new", "pi_new_multi", "pi_old", "pi_vcftools"))
    if (length(cols) == 0) return(tibble())
    df |>
        select(replicate, all_of(cols)) |>
        pivot_longer(all_of(cols), names_to = "estimator", values_to = "value") |>
        mutate(
            estimator     = str_remove(estimator, "^pi_"),
            arm           = arm_row$arm,
            ploidy        = arm_row$ploidy,
            theta_nominal = arm_row$theta_nominal,
            theoretical   = arm_row$theoretical_pi
        )
}

all <- bind_rows(lapply(seq_len(nrow(theta_arms)), function(i) load_arm(theta_arms[i, ])))

if (nrow(all) == 0) {
    cat("No Phase G aggregated TSVs found yet.\n")
    quit(status = 0)
}

summary_tbl <- all |>
    group_by(arm, ploidy, theta_nominal, estimator) |>
    summarise(
        mean   = mean(value, na.rm = TRUE),
        q025   = quantile(value, 0.025, na.rm = TRUE),
        q975   = quantile(value, 0.975, na.rm = TRUE),
        theoretical = first(theoretical),
        .groups = "drop"
    ) |>
    mutate(
        bias        = mean - theoretical,
        bias_pct    = 100 * bias / theoretical
    )

# one row per (ploidy, theta, estimator). a second arm matching the filter would
# add a duplicate row (two points per theta in Figure 4's pi row). fail loud.
stopifnot(
    nrow(summary_tbl) ==
        nrow(distinct(summary_tbl, ploidy, theta_nominal, estimator))
)

write_tsv(summary_tbl, "multiallelic_frac_summary.tsv")

p <- summary_tbl |>
    ggplot(aes(theta_nominal, bias_pct, color = estimator, linetype = ploidy, shape = ploidy)) +
    geom_line() +
    geom_pointrange(aes(ymin = 100 * (q025 - theoretical) / theoretical,
                        ymax = 100 * (q975 - theoretical) / theoretical),
                    position = position_dodge(width = 0.005), size = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
    scale_x_log10() +
    scale_color_brewer(palette = "Dark2") +
    labs(x = expression(theta == 4 * N[e] * mu ~ "(log scale)"),
         y = "bias relative to theory (%)",
         title = "Bias of pi vs theta (multiallelic site fraction proxy)",
         subtitle = paste("biallelic-only `new` should drift below zero as theta increases;",
                          "`new_multi` should track zero"),
         color = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position = "top")

ggsave(file.path(fig_dir, "08_multiallelic_frac.pdf"), p, width = 8, height = 5)

cat("Wrote analysis/multiallelic_frac_summary.tsv and figs/08_multiallelic_frac.pdf\n")
