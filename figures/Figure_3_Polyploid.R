# figure 3 — ploidy x missingness sweep, multiallelic estimator (new_multi)
# theoretical lines: finite-sites JC69 per stat, read from the summary TSVs
# (thetaw = E(s*)/a_n, Tajima 1996 eq 15; fst = structured coalescent, alpha = 2)

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
fig_path <- file.path(fig_dir, "Figure3_polyploid.pdf")

##########
# load
##########
miss_pi_dxy_fst <- read_tsv(file.path(POLY_DIR, "missingness_summary.tsv"),
                            show_col_types = FALSE) |>
  filter(stat %in% c("pi", "dxy", "fst_hudson")) |>
  select(stat, ploidy, miss_pct, estimator, mean, ci_lo, ci_hi, theoretical)

miss_tw_td <- read_tsv(file.path(POLY_DIR, "thetaw_tajimaD_summary.tsv"),
                       show_col_types = FALSE) |>
  filter(stat %in% c("thetaw", "tajimaD")) |>
  select(stat, ploidy, miss_pct, estimator, mean, ci_lo, ci_hi, theoretical)

dat <- bind_rows(miss_pi_dxy_fst, miss_tw_td) |>
  filter(estimator == "new_multi") |>
  mutate(
    ploidy = factor(ploidy, levels = c(2, 4, 6, 8),
                    labels = c("2n", "4n", "6n", "8n"))
  )

# per-rep values for the violins (mean-of-windows per rep)
per_rep <- bind_rows(
  read_tsv(file.path(POLY_DIR, "per_rep_missingness.tsv"), show_col_types = FALSE),
  read_tsv(file.path(POLY_DIR, "per_rep_thetaw_tajimaD.tsv"), show_col_types = FALSE)
) |>
  filter(estimator == "new_multi",
         stat %in% c("pi", "dxy", "fst_hudson", "thetaw", "tajimaD")) |>
  mutate(
    ploidy = factor(ploidy, levels = c(2, 4, 6, 8),
                    labels = c("2n", "4n", "6n", "8n")),
    stat   = factor(stat, levels = c("pi", "dxy", "fst_hudson",
                                     "thetaw", "tajimaD"))
  )

stat_labels <- c(
  pi         = "pi",
  dxy        = "d[xy]",
  fst_hudson = "F[ST]",
  thetaw     = "theta[W]",
  tajimaD    = "\"Tajima's D\""
)

dat <- dat |>
  mutate(stat = factor(stat, levels = c("pi", "dxy", "fst_hudson",
                                        "thetaw", "tajimaD")))

# theoretical per (stat, ploidy, miss); flat for pi/dxy/fst/thetaw
# tajimaD is n-normalised so its finite-sites floor slopes with missingness
th <- dat |>
  distinct(stat, ploidy, miss_pct, theoretical) |>
  filter(!is.na(theoretical))

multi_col <- "#0F6F57"

# per-panel y-limits anchored to the theoretical line so every free_y axis
# spans a comparable relative window. multiplicative band for pi/dxy/fst/thetaw;
# tajimaD is near-zero so pad it by an absolute +/-0.25 instead. geom_blank
# carries the limits through free_y
D_PAD <- 0.25
ylim_panel <- dat |>
  group_by(stat, ploidy) |>
  summarise(
    th_lo = min(theoretical, na.rm = TRUE),
    th_hi = max(theoretical, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    ymin = if_else(stat == "tajimaD", th_lo - D_PAD, 0.90 * th_lo),
    ymax = if_else(stat == "tajimaD", th_hi + D_PAD, 1.10 * th_hi)
  ) |>
  select(stat, ploidy, ymin, ymax) |>
  tidyr::pivot_longer(c(ymin, ymax), values_to = "y") |>
  mutate(miss_pct = 0)

fig <- ggplot(dat, aes(miss_pct, mean)) +
  geom_violin(data = per_rep, aes(x = miss_pct, y = value, group = miss_pct),
              width = 0.12, colour = NA, fill = multi_col, alpha = 0.25,
              position = "identity", inherit.aes = FALSE) +
  geom_line(aes(colour = "pixy 2.0.0"), linewidth = 0.7) +
  geom_point(size = 1.6, colour = multi_col) +
  geom_line(data = th,
            aes(miss_pct, theoretical, linetype = "theoretical"),
            colour = "grey20", linewidth = 0.4, inherit.aes = FALSE) +
  geom_blank(data = ylim_panel, aes(miss_pct, y), inherit.aes = FALSE) +
  facet_grid(stat ~ ploidy,
             scales = "free_y",
             labeller = labeller(
               stat   = as_labeller(stat_labels, label_parsed),
               ploidy = c(`2n` = "Diploid", `4n` = "Tetraploid",
                          `6n` = "Hexaploid", `8n` = "Octoploid"))) +
  scale_x_continuous(labels = scales::percent,
                     breaks = c(0, .25, .5, .75)) +
  scale_colour_manual(values = c("pixy 2.0.0" = multi_col), name = NULL) +
  scale_linetype_manual(values = c(theoretical = "dashed"), name = NULL) +
  labs(x = "Missingness", y = "Mean estimate") +
  theme_pixy() +
  theme(strip.text.x = element_text(hjust = 0.5, face = "plain",
                                    size = (PIXY_BASE_SIZE - 1) * 0.85),
        strip.text.y = element_text(angle = -90, face = "plain"),
        legend.position   = c(0.03, 1.00),
        legend.justification = c(0, 1),
        legend.direction  = "horizontal",
        legend.box        = "horizontal",
        legend.spacing    = unit(1, "pt"),
        legend.key.spacing.x = unit(1, "pt"),
        legend.key.width  = unit(0.8, "lines"),
        legend.margin     = margin(0, 0, 0, 0),
        legend.background = element_blank(),
        legend.key        = element_blank(),
        legend.text       = element_text(size = PIXY_BASE_SIZE - 5))

fig

pixy_save(fig, fig_path, width = 9, height = 8)
