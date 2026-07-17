# figure 5 — empirical biallelic vs +multiallelic in anopheles (2n) and arenosa (4n)
# each arm run twice through pixy (biallelic-only, --include_multiallelic_snps),
# per-window pi/dxy/fst joined and compared

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(forcats)
  library(purrr)
  library(ggplot2)
  library(ggh4x)        # facet_grid2: independent scales per panel
  library(patchwork)
})

script_dir <- tryCatch(
  dirname(rstudioapi::getActiveDocumentContext()$path),
  error = function(e) getwd()
)
here <- function(...) file.path(script_dir, ...)

source(here("figure_theme.R"))

EMP_DIR  <- here("..", "analyses", "empirical-tests", "data", "aggregated")
ARMS_TSV <- here("..", "analyses", "empirical-tests", "config", "arms.tsv")

fig_dir  <- here("figs")
fig_path <- file.path(fig_dir, "Figure5_empirical.pdf")

##########
# load
##########
arms <- read_tsv(ARMS_TSV, show_col_types = FALSE) |>
  mutate(
    species_label = {
      parts <- str_split_fixed(species, "_", 2)
      sprintf("%s. %s", str_to_upper(str_sub(parts[, 1], 1, 1)), parts[, 2])
    }
  )

load_arm <- function(arm_id) {
  path <- file.path(EMP_DIR, paste0(arm_id, ".joined.tsv"))
  if (!file.exists(path)) {
    stop("missing aggregated TSV: ", path,
         "\nRun 07_compare.sbatch and sync data/aggregated/ before plotting.")
  }
  read_tsv(path, show_col_types = FALSE) |> mutate(arm_id = arm_id)
}

raw <- map_dfr(arms$arm_id, load_arm) |>
  left_join(arms |> select(arm_id, species_label, ploidy), by = "arm_id") |>
  mutate(species_label = factor(species_label, levels = arms$species_label))

##########
# reshape
##########
# pi: avg_pi_<POP>_<variant>; pivot to (arm, pop, biallelic, multi)
pi_long <- raw |>
  select(species_label, chromosome, window_pos_1, starts_with("avg_pi_")) |>
  pivot_longer(
    cols          = starts_with("avg_pi_"),
    names_to      = c("pop", "variant"),
    names_pattern = "avg_pi_(.+)_(biallelic|multi)",
    values_to     = "value"
  ) |>
  pivot_wider(names_from = variant, values_from = value) |>
  mutate(statistic = "pi", group = pop) |>
  select(species_label, statistic, group, chromosome, window_pos_1, biallelic, multi)

dxy_long <- raw |>
  transmute(species_label,
            statistic = "dxy",
            group     = paste(pop1, pop2, sep = "-"),
            chromosome, window_pos_1,
            biallelic = dxy_biallelic,
            multi     = dxy_multi)

fst_long <- raw |>
  transmute(species_label,
            statistic = "fst",
            group     = paste(pop1, pop2, sep = "-"),
            chromosome, window_pos_1,
            biallelic = fst_biallelic,
            multi     = fst_multi)

long <- bind_rows(pi_long, dxy_long, fst_long) |>
  filter(is.finite(biallelic), is.finite(multi)) |>
  mutate(statistic = factor(statistic, levels = c("pi", "dxy", "fst")))

##########
# per-facet R^2 + slope
##########
fit_stats <- long |>
  group_by(species_label, statistic) |>
  summarise(
    n     = n(),
    r2    = cor(biallelic, multi)^2,
    slope = coef(lm(multi ~ biallelic))[2],
    .groups = "drop"
  ) |>
  mutate(label = sprintf("R^2 == %.3f", r2))

##########
# panel a: scatter, biallelic vs +multiallelic
##########
stat_labeller <- as_labeller(
  c(pi = "pi", dxy = "d[xy]", fst = "F[ST]"),
  default = label_parsed
)

p_a <- ggplot(long, aes(biallelic, multi)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_point(alpha = 0.35, size = 0.9,
             colour = "#3C3489") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              colour = "black", linewidth = 0.5) +
  geom_text(
    data = fit_stats,
    aes(x = -Inf, y = Inf, label = label),
    parse = TRUE, hjust = -0.1, vjust = 1.2, size = 3.4, inherit.aes = FALSE
  ) +
  facet_grid2(species_label ~ statistic,
              scales = "free", independent = "all",
              labeller = labeller(statistic = stat_labeller)) +
  labs(tag = "a",
       x = "biallelic-only estimate",
       y = "+multiallelic estimate") +
  theme_pixy() +
  theme(strip.text.x = element_text(face = "plain"),
        strip.text.y = element_text(face = "italic", angle = 270))

##########
# panel b: percent elevation vs biallelic estimate
##########
# bias <- long |>
#   filter(biallelic > 0) |>
#   mutate(pct_elev = 100 * (multi - biallelic) / biallelic)
# 
# p_b <- ggplot(bias, aes(biallelic, pct_elev, colour = species_label)) +
#   geom_hline(yintercept = 0, linetype = "dashed",
#              colour = "grey40", linewidth = 0.4) +
#   geom_point(alpha = 0.35, size = 0.9) +
#   geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.6) +
#   facet_wrap(~ statistic, scales = "free", ncol = 3,
#              labeller = labeller(statistic = stat_labeller)) +
#   scale_colour_brewer(palette = "Set2", name = NULL) +
#   labs(tag = "b",
#        x = "biallelic-only estimate",
#        y = "% increase vs. biallelic estimate") +
#   theme_pixy() +
#   theme(
#     legend.position      = c(0.28, 0.93),
#     legend.justification = c(1, 1),
#     legend.background    = element_rect(fill = alpha("white", 0.85), colour = NA),
#     legend.text          = element_text(face = "italic", size = PIXY_BASE_SIZE - 3),
#     legend.key.size      = unit(0.5, "lines"),
#     legend.margin        = margin(2, 4, 2, 4)
#   )

##########
# compose
##########

fig <- p_a

# fig <- (p_a / p_b) +
#   plot_layout(heights = c(2, 1)) +
#   plot_annotation(theme = theme(legend.position = "top"))

fig

pixy_save(fig, fig_path, width = 11, height = 10)
