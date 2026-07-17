#!/usr/bin/env Rscript
# pi and dxy across the ploidy grid {2,4,6,8}, 1-pop and 2-pop, zero missingness.
# reads the diploid/tetraploid arms + pi_hex_1pop, pi_oct_1pop, dxy_hex_2pop,
# dxy_oct_2pop.
# -> figs/04_ploidy_grid.pdf (mean +/- 95% interval per estimator, theoretical dashed)

suppressPackageStartupMessages({
    library(tidyverse)
})

fig_dir <- "figs"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

sim_params <- read_tsv("../config/sim_params.tsv", show_col_types = FALSE)

ploidy_arms <- sim_params |>
    filter(arm %in% c(
        "pi_dip_1pop", "pi_tet_1pop", "pi_hex_1pop", "pi_oct_1pop",
        "dxy_dip_2pop", "dxy_tet_2pop", "dxy_hex_2pop", "dxy_oct_2pop"
    )) |>
    mutate(
        stat   = sapply(strsplit(stats, ","), `[`, 1),
        ploidy = ploidy_chr1,
        theta_nom   = 4 * Ne * mu,
        theoretical = if_else(
            stat == "pi",
            theta_nom / (1 + (4/3) * theta_nom),
            (3/4) * (1 - exp(-8 * mu * split_time / 3) / (1 + (8/3) * theta_nom))
        )
    )

load_arm <- function(arm_row) {
    f <- file.path("../data/aggregated", paste0(arm_row$arm, ".tsv"))
    if (!file.exists(f)) {
        warning("missing: ", f)
        return(tibble())
    }
    df <- read_tsv(f, show_col_types = FALSE)
    # old + new schema both work for this stat (pi or dxy)
    value_cols <- intersect(names(df), paste0(arm_row$stat, c("_old", "_new", "_new_multi")))
    if (length(value_cols) == 0) return(tibble())
    df |>
        select(replicate, all_of(value_cols)) |>
        pivot_longer(all_of(value_cols), names_to = "estimator", values_to = "value") |>
        mutate(
            estimator   = str_remove(estimator, paste0("^", arm_row$stat, "_")),
            arm         = arm_row$arm,
            stat        = arm_row$stat,
            ploidy      = arm_row$ploidy,
            theoretical = arm_row$theoretical
        )
}

all <- bind_rows(lapply(seq_len(nrow(ploidy_arms)), function(i) load_arm(ploidy_arms[i, ])))

if (nrow(all) == 0) {
    cat("No ploidy-arm aggregated TSVs found yet.\n")
    quit(status = 0)
}

summary_tbl <- all |>
    group_by(stat, ploidy, estimator) |>
    summarise(
        n_windows = n(),
        mean      = mean(value, na.rm = TRUE),
        q025      = quantile(value, 0.025, na.rm = TRUE),
        q975      = quantile(value, 0.975, na.rm = TRUE),
        theoretical = first(theoretical),
        .groups   = "drop"
    )

write_tsv(summary_tbl, "ploidy_grid_summary.tsv")

p <- summary_tbl |>
    ggplot(aes(factor(ploidy), mean, color = estimator)) +
    geom_pointrange(aes(ymin = q025, ymax = q975),
                    position = position_dodge(width = 0.5),
                    size = 0.4) +
    geom_hline(aes(yintercept = theoretical),
               linetype = "dashed", linewidth = 0.4) +
    facet_wrap(~ stat, scales = "free_y") +
    scale_color_brewer(palette = "Dark2") +
    labs(x = "ploidy", y = "estimate (mean +/- 95%)",
         title = "pi and dxy across the ploidy grid (no missing data)",
         subtitle = "Dashed = theoretical expectation",
         color = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position = "top")

ggsave(file.path(fig_dir, "04_ploidy_grid.pdf"), p,
       width = 8, height = 5)

cat("Wrote analysis/ploidy_grid_summary.tsv and figs/04_ploidy_grid.pdf\n")
