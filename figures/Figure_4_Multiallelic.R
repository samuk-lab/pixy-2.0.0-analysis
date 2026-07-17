# figure 4 — biallelic vs multiallelic-aware pi/dxy/fst across a theta sweep
# 2n and 8n only; renders pi-only if dxy_fst_theta_summary.tsv is missing

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(patchwork)
})

script_dir <- tryCatch(
  dirname(rstudioapi::getActiveDocumentContext()$path),
  error = function(e) getwd()
)
here <- function(...) file.path(script_dir, ...)

source(here("figure_theme.R"))

POLY_DIR <- here("..", "analyses", "polyploid_multiallele", "analysis")

fig_dir  <- here("figs")
fig_path <- file.path(fig_dir, "Figure4_multiallelic.pdf")

KEEP_ESTIMATORS <- c("new", "new_multi")
KEEP_PLOIDY     <- c("diploid", "octoploid")

##########
# load
##########
pi_tbl <- read_tsv(file.path(POLY_DIR, "multiallelic_frac_summary.tsv"),
                   show_col_types = FALSE) |>
  mutate(stat = "pi") |>
  select(arm, ploidy, theta_nominal, stat, estimator,
         mean, q025, q975, theoretical, bias_pct)

dxy_fst_path <- file.path(POLY_DIR, "dxy_fst_theta_summary.tsv")
if (file.exists(dxy_fst_path)) {
  dxy_fst_tbl <- read_tsv(dxy_fst_path, show_col_types = FALSE) |>
    select(arm, ploidy, theta_nominal, stat, estimator,
           mean, q025, q975, theoretical, bias_pct)
} else {
  message("dxy_fst_theta_summary.tsv not found; rendering pi only. ",
          "Run analysis/11_dxy_fst_theta_sweep.R after the dxy/fst arms ",
          "complete on the cluster.")
  dxy_fst_tbl <- tibble()
}

dat <- bind_rows(pi_tbl, dxy_fst_tbl) |>
  filter(estimator %in% KEEP_ESTIMATORS,
         ploidy    %in% KEEP_PLOIDY) |>
  mutate(
    ploidy    = factor(ploidy,
                       levels = c("diploid", "octoploid"),
                       labels = c("2n", "8n")),
    estimator = factor(estimator, levels = KEEP_ESTIMATORS),
    stat      = factor(stat, levels = c("pi", "dxy", "fst_hudson"))
  )

stat_labels <- c(
  pi         = "pi",
  dxy        = "d[xy]",
  fst_hudson = "F[ST]"
)

##########
# theoretical reference per (stat, ploidy, theta)
##########
# one reference line per panel = finite-sites expectation in all three rows.
# pi/dxy: E[pi] under JC finite sites (0.088 at theta = 0.1, not 0.100).
# fst: E[FST]_finite = 1 - E[pi_w]/E[dxy], verified against these data 2026-07-15.
# biallelic fst sits ~1.6% high here because it recovers the coalescent fst
# (flat 0.237, Slatkin 1991) — a different estimand, not a broken estimator.
# no single line flatters both; figure uses the same estimand as the pi/dxy rows.
# analysis/finite_sites_fst_derivation.md carries both; summary TSV keeps both
# expectations for the inline stats
th_fst <- dat |>
  filter(stat == "fst_hudson", estimator == "new_multi") |>
  distinct(stat, ploidy, theta_nominal, theoretical)

th <- dat |>
  filter(stat != "fst_hudson") |>
  distinct(stat, ploidy, theta_nominal, theoretical) |>
  bind_rows(th_fst) |>
  filter(!is.na(theoretical)) |>
  arrange(stat, ploidy, theta_nominal)

# one expectation per (stat, ploidy, theta); a second row makes geom_line zigzag
stopifnot(nrow(th) == nrow(distinct(th, stat, ploidy, theta_nominal)))

# short estimator labels for the colour legend
estimator_labels_short <- c(new = "biallelic", new_multi = "multiallelic")

##########
# panel a: mean estimate vs theta
##########
p_a <- ggplot(dat,
              aes(theta_nominal, mean,
                  colour = estimator, group = estimator)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2.0) +
  geom_line(data = th,
            aes(x = theta_nominal, y = theoretical, linetype = "theoretical"),
            colour = "grey20", linewidth = 0.4, inherit.aes = FALSE) +
  # fix the fst row to 0.223 .. 0.247: the ~0.004 divergence is ~16% of panel
  # height and legible, while 1 bootstrap SE (~5e-4) reads as noise. fixed not
  # free so both ploidy columns share a scale
  geom_blank(data = tidyr::expand_grid(
               stat          = factor("fst_hudson", levels = levels(dat$stat)),
               ploidy        = factor(levels(dat$ploidy),
                                      levels = levels(dat$ploidy)),
               mean          = c(0.223, 0.247)) |>
             dplyr::mutate(theta_nominal = min(dat$theta_nominal, na.rm = TRUE)),
             aes(theta_nominal, mean), inherit.aes = FALSE) +
  facet_grid(stat ~ ploidy,
             scales = "free_y",
             labeller = labeller(
               stat   = as_labeller(stat_labels, label_parsed),
               ploidy = c(`2n` = "Diploid", `8n` = "Octoploid"))) +
  scale_x_log10(breaks = c(0.005, 0.01, 0.025, 0.05, 0.1)) +
  scale_colour_manual(values = c(new = "#0F6F57", new_multi = "#A90264"),
                      labels = estimator_labels_short, name = NULL,
                      guide  = guide_legend(direction = "vertical",
                                            override.aes = list(shape = NA))) +
  scale_linetype_manual(values = c(theoretical = "longdash"), name = NULL) +
  labs(x = expression(theta == 4 * N[e] * mu),
       y = "Mean estimate") +
  theme_pixy() +
  theme(strip.text.x      = element_text(hjust = 0.5, face = "plain", size = (PIXY_BASE_SIZE - 1) * 0.85),
        strip.text.y      = element_text(angle = -90, face = "plain"),
        legend.position   = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.direction  = "vertical",
        legend.spacing    = unit(0, "pt"),
        legend.key.width  = unit(1.4, "lines"),
        legend.margin     = margin(0, 0, 0, 0),
        legend.text       = element_text(size = PIXY_BASE_SIZE - 5),
        legend.background = element_blank(),
        axis.text         = element_text(size = (PIXY_BASE_SIZE - 1) * 0.75))

##########
# compose
##########
fig <- p_a

fig

pixy_save(fig, fig_path, width = 6, height = 6)
