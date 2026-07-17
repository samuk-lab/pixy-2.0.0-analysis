# figure 2 — multicore speedup and peak rss

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
})

script_dir <- tryCatch(
  dirname(rstudioapi::getActiveDocumentContext()$path),
  error = function(e) getwd()
)
here <- function(...) file.path(script_dir, ...)

source(here("figure_theme.R"))

fig_dir  <- here("figs")
fig_path <- file.path(fig_dir, "Figure2_performance.pdf")

##########
# load and clean
##########
dat <- read_tsv(here("..", "analyses", "benchmark_multicore",
                     "data", "results", "all_cells_long.tsv"),
                col_types = cols(.default = "c"), show_col_types = FALSE) |>
  mutate(
    elapsed_s_clean = as.numeric(str_extract(elapsed_s, "^[0-9.]+")),
    n_cores_num     = suppressWarnings(as.integer(n_cores)),
    statistic       = factor(statistic, levels = c("pi", "dxy", "fst")),
    pixy_version    = factor(pixy_version, levels = c("0.93.1", "2.1.2"))
  ) |>
  filter(status == "OK", !is.na(elapsed_s_clean))

# per-cell median
medians <- dat |>
  group_by(pixy_version, statistic, n_cores_num) |>
  summarise(median_s = median(elapsed_s_clean), .groups = "drop")

##########
# speedup vs 0.93.1 single-core
##########
baseline_old <- medians |>
  filter(pixy_version == "0.93.1") |>
  select(statistic, baseline_old = median_s)

speedup_old <- medians |>
  filter(pixy_version == "2.1.2") |>
  inner_join(baseline_old, by = "statistic") |>
  mutate(speedup = baseline_old / median_s)

##########
# panel a: speedup vs 0.93.1
##########
speedup_old_labels <- speedup_old %>%
  arrange(statistic, n_cores_num) %>%
  mutate(median_s = round(median_s, 1)) %>%
  group_by(statistic) %>%
  filter(row_number() %% 2 == 1) %>%  # every other point
  ungroup()

(p_a <- ggplot(speedup_old,
       aes(x = n_cores_num, y = speedup, colour = statistic, group = statistic)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey60",
             linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.4) +
  geom_text_repel(
    data = speedup_old_labels,
    aes(label = paste0(median_s,"s")),
    size = 3,
    segment.size = 0.3,         
    segment.alpha = 0.7,
    segment.linetype = "dotted",
    force_pull = 1,
    show.legend = FALSE,
    min.segment.length = 0,
    nudge_x = -1.0, 
    nudge_y = 0.5,
    box.padding = 0.25
  ) +
  scale_x_continuous(breaks = c(1, 2, 4, 8, 16), trans = "log2") +
  scale_colour_pixy_statistic() +
  scale_y_continuous(labels = function(x) paste0(x, "×")) +
  labs(
    tag = "a",
    x = "pixy 2.0.0 cores",
    y = "Relative speed vs. pixy 0.95.02"
  ) +
  theme_pixy(base_size = 10) +
  theme(
    legend.position = c(0.03, 0.97),
    legend.justification = c(0, 1),
    legend.background = element_rect(colour = "grey80", linewidth = 0.3),
    legend.title = element_blank()
  ))

##########
# panel b: peak rss vs 0.93.1
##########

# reload
dat <- read_tsv(here("..", "analyses", "benchmark_multicore",
                     "data", "results", "all_cells_long.tsv"),
                col_types = cols(.default = "c"), show_col_types = FALSE) %>%
  mutate(
    elapsed_s_clean = as.numeric(str_extract(elapsed_s, "^[0-9.]+")),
    rss_kb_clean = suppressWarnings(as.numeric(rss_kb)),
    n_cores_num  = suppressWarnings(as.integer(n_cores)),
    pixy_version = factor(pixy_version, levels = c("0.93.1", "2.1.2")),
    statistic    = factor(statistic, levels = c("pi", "dxy", "fst"),
                          labels = c("pi", "dxy", "fst")),
    arm_label = case_when(
      pixy_version == "0.93.1" ~ "0.93.1\n(1)",
      n_cores_num == 1         ~ "2.0.0\n(1)",
      TRUE                     ~ paste0("2.0.0\n(", n_cores_num, ")")
    ),
    arm_order = ifelse(pixy_version == "0.93.1", -1L, n_cores_num),
    arm_label = fct_reorder(arm_label, arm_order)
  ) |>
  filter(status == "OK", !is.na(elapsed_s_clean))

rss_summary <- dat %>%
  filter(!is.na(rss_kb_clean), rss_kb_clean > 0) %>%
  mutate(rss_mb = rss_kb_clean / 1024) %>%
  group_by(statistic) %>%
  mutate(
    baseline_rss = mean(rss_mb[arm_label == "0.93.1\n(1)"], na.rm = TRUE),
    rel_rss = rss_mb / baseline_rss
  ) %>%
  ungroup() %>%
  filter(arm_label != "0.93.1\n(1)") %>%
  group_by(statistic, arm_label, n_cores) %>%
  summarise(
    n           = n(),
    mean_rss    = mean(rel_rss),
    mean_rss_mb = mean(rss_mb),
    se_rss      = sd(rel_rss) / sqrt(n),
    ci_rss      = qt(0.975, df = n - 1) * se_rss,
    ymin        = mean_rss - ci_rss,
    ymax        = mean_rss + ci_rss,
    .groups     = "drop"
  ) %>%
  mutate(n_cores = factor(n_cores, levels = c(1, 2, 4, 8, 16)))

rss_labels <- rss_summary %>%
  arrange(statistic, n_cores) %>%
  group_by(statistic) %>%
  filter(row_number() %% 2 == 1) %>%   # every other point
  ungroup()

p_b <- ggplot(rss_summary,
              aes(x = n_cores, y = mean_rss, color = statistic, group = statistic)) +
  geom_point(size = 2.4) +
  geom_line(linewidth = 0.7) +
  geom_text_repel(
    data = rss_labels,
    aes(label = sprintf("%.0f MB", mean_rss_mb)),
    size = 3,
    segment.size = 0.3,
    segment.alpha = 0.7,
    segment.linetype = "dotted",
    force_pull = 1,
    show.legend = FALSE,
    min.segment.length = 0,
    nudge_x = 0.5,
    nudge_y = 0.0005,
    box.padding = 0.25
  ) +
  labs(
    tag = "b",
    x   = "pixy 2.0.0 cores",
    y   = "Percent of pixy 0.95.02 peak RSS"
  ) +
  theme_pixy(base_size = 10) +
  scale_colour_pixy_statistic() +
  scale_y_continuous(labels = scales::label_percent()) +
  theme(
    legend.position = c(0.97, 0.97),
    legend.justification = c(1, 1),
    legend.background = element_rect(colour = "grey80", linewidth = 0.3),
    legend.title = element_blank()
  )

##########
# compose
##########
(fig <- (p_a | p_b))


pixy_save(fig, fig_path, width = 7.1, height = 3.5)
