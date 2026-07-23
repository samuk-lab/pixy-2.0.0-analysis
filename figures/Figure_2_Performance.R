# figure 2 — multicore speedup and peak rss

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
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
    pixy_version    = factor(pixy_version, levels = c("0.95.01", "2.2.3"))
  ) |>
  filter(status == "OK", !is.na(elapsed_s_clean))

# per-cell median
medians <- dat |>
  group_by(pixy_version, statistic, n_cores_num) |>
  summarise(median_s = median(elapsed_s_clean), .groups = "drop")

##########
# speedup vs 0.95.01 single-core
##########
baseline_old <- medians |>
  filter(pixy_version == "0.95.01") |>
  select(statistic, baseline_old = median_s)

speedup_old <- medians |>
  filter(pixy_version == "2.2.3") |>
  inner_join(baseline_old, by = "statistic") |>
  mutate(speedup = baseline_old / median_s)

##########
# panel a: speedup vs 0.95.01
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
    x = "pixy 2.2.3 cores",
    y = "Relative speed vs. pixy 0.95.01"
  ) +
  theme_pixy(base_size = 10) +
  theme(
    legend.position = c(0.03, 0.97),
    legend.justification = c(0, 1),
    legend.background = element_rect(colour = "grey80", linewidth = 0.3),
    legend.title = element_blank()
  ))

##########
# panel b: peak memory of the whole process tree
##########
# Tree_peak_kb (from 05_memory_aggregate.sbatch) is the peak *sum* over all live
# processes, not the largest single process — the %M column in all_cells_long.tsv
# is per-process and undercounts a multicore run. memory arms were run at 1/4/16
# cores only, so panel b has three x positions where panel a has five.
mem_dir <- here("..", "analyses", "benchmark_multicore", "data", "results")

read_mem <- function(path) {
  read_tsv(path, col_types = cols(.default = "c"), show_col_types = FALSE) |>
    mutate(statistic = str_match(basename(path),
                                 "pixy_mem_(?:old_)?(pi|dxy|fst)_")[, 2])
}

mem <- list.files(mem_dir, pattern = "^pixy_mem_.*\\.tsv$", full.names = TRUE) |>
  lapply(read_mem) |>
  bind_rows() |>
  filter(Status == "OK") |>
  mutate(
    tree_peak_mb = as.numeric(Tree_peak_kb) / 1024,
    proc_peak_mb = as.numeric(Proc_peak_kb) / 1024,
    n_cores_num  = as.integer(Cores),
    statistic    = factor(statistic, levels = c("pi", "dxy", "fst")),
    pixy_version = factor(Pixy_version, levels = c("0.95.01", "2.2.3"))
  )

# 0.95.01 is single-process and single-core; its footprint is within ~1% across
# the three statistics, so one pooled baseline line rather than three
baseline_mb <- mem |>
  filter(pixy_version == "0.95.01") |>
  pull(tree_peak_mb) |>
  median()

mem_new <- mem |> filter(pixy_version == "2.2.3")

mem_summary <- mem_new |>
  group_by(statistic, n_cores_num) |>
  summarise(
    n       = n(),
    mean_mb = mean(tree_peak_mb),
    se      = sd(tree_peak_mb) / sqrt(n),
    ci      = qt(0.975, df = n - 1) * se,
    ymin    = mean_mb - ci,
    ymax    = mean_mb + ci,
    .groups = "drop"
  )

# tree peak is linear in cores (per-worker footprint is flat, see proc_peak_mb),
# so solve the fit for the core count at which 2.2.3 reaches the 0.95.01 footprint
crossover <- mem_new |>
  group_by(statistic) |>
  group_modify(~ {
    fit <- lm(tree_peak_mb ~ n_cores_num, data = .x)
    tibble(
      intercept    = coef(fit)[[1]],
      mb_per_core  = coef(fit)[[2]],
      cores_at_baseline = (baseline_mb - coef(fit)[[1]]) / coef(fit)[[2]]
    )
  }) |>
  ungroup()

print(mem_summary, n = Inf)
print(crossover)
message("Pooled 0.95.01 baseline: ", round(baseline_mb), " MB; ",
        "median per-worker peak: ",
        round(median(mem_new$proc_peak_mb)), " MB")

crossover_cores <- mean(crossover$cores_at_baseline)

p_b <- ggplot(mem_summary,
              aes(x = n_cores_num, y = mean_mb, colour = statistic,
                  group = statistic)) +
  geom_hline(yintercept = baseline_mb, linetype = "dashed",
             colour = "grey40", linewidth = 0.4) +
  annotate("segment", x = crossover_cores, xend = crossover_cores,
           y = min(mem_summary$ymin) * 0.9, yend = baseline_mb,
           linetype = "dotted",
           colour = "grey40", linewidth = 0.4) +
  annotate("text", x = 1, y = baseline_mb * 1.12, hjust = 0, vjust = 0,
           size = 2.8, colour = "grey30", family = PIXY_FONT,
           label = sprintf("pixy 0.95.01, 1 core (%.0f MB)", baseline_mb)) +
  annotate("text", x = crossover_cores * 0.93, y = min(mem_summary$ymin) * 0.9,
           hjust = 1, vjust = 0,
           size = 2.8, colour = "grey30", family = PIXY_FONT,
           label = sprintf("parity at %.1f cores", crossover_cores)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.06, linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.4) +
  scale_x_continuous(breaks = c(1, 2, 4, 8, 16), trans = "log2") +
  scale_y_continuous(trans = "log2",
                     breaks = c(125, 250, 500, 1000, 2000),
                     labels = function(x) paste0(x, " MB")) +
  scale_colour_pixy_statistic() +
  labs(
    tag = "b",
    x   = "pixy 2.2.3 cores",
    y   = "Peak memory, whole process tree"
  ) +
  theme_pixy(base_size = 10) +
  theme(
    legend.position = c(0.97, 0.03),
    legend.justification = c(1, 0),
    legend.background = element_rect(colour = "grey80", linewidth = 0.3),
    legend.title = element_blank()
  )

##########
# compose
##########
(fig <- (p_a | p_b))


pixy_save(fig, fig_path, width = 7.1, height = 3.5)
