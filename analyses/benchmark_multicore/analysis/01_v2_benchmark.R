library(tidyverse)

# plot the pixy multicore benchmark:
#   0.93.1 single core (zarr, n_cores == "0.95.01") vs 2.0.0 at 1/2/4/8/16 cores (tabix)
# in:  ../data/results/all_cells_long.tsv (from 04_aggregate_summaries.sh)
# out: figs/
# 2.0.0 per-cell tsvs pack timing into Elapsed_s as "elapsed user sys rss";
# elapsed_s and rss_kb parsed out below

script_dir <- tryCatch(
  dirname(rstudioapi::getActiveDocumentContext()$path),
  error = function(e) getwd()
)
here <- function(...) file.path(script_dir, ...)

fig_dir <- here("figs")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

##########
# load and clean data
##########
dat <- read_tsv(here("../data/results/all_cells_long.tsv"),
                col_types = cols(.default = "c"),
                show_col_types = FALSE)

dat <- dat |>
  mutate(
    # elapsed_s for 2.0.0 is "elapsed user sys rss"
    elapsed_s_clean = as.numeric(str_extract(elapsed_s, "^[0-9.]+")),
    # rss = last token; 0.93.1 rss_kb already clean
    rss_kb_clean = case_when(
      pixy_version == "2.0.0" ~ as.numeric(str_extract(elapsed_s, "[0-9]+$")),
      TRUE                    ~ as.numeric(rss_kb)
    ),
    seed         = as.integer(seed),
    n_cores_num  = suppressWarnings(as.integer(n_cores)),
    pixy_version = factor(pixy_version, levels = c("0.93.1", "2.0.0")),
    arm_label = case_when(
      pixy_version == "0.93.1"              ~ "0.93.1\n(1 core)",
      n_cores_num == 1                       ~ "2.0.0\n(1 core)",
      TRUE                                   ~ paste0("2.0.0\n(", n_cores_num, " cores)")
    ),
    arm_label = fct_reorder(arm_label, ifelse(is.na(n_cores_num), -1L, n_cores_num))
  )

ok <- dat |> filter(status == "OK", !is.na(elapsed_s_clean))

##########
# summary stats
##########
summary_tab <- ok |>
  group_by(pixy_version, statistic, n_cores) |>
  summarise(
    n             = n(),
    median_s      = median(elapsed_s_clean),
    mean_s        = mean(elapsed_s_clean),
    sd_s          = sd(elapsed_s_clean),
    median_rss_mb = median(rss_kb_clean / 1024, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_tab, n = Inf)

##########
# speedup vs 2.0.0 at 1 core
##########
baseline_1core <- ok |>
  filter(pixy_version == "2.0.0", n_cores == "1") |>
  group_by(statistic) |>
  summarise(baseline_s = median(elapsed_s_clean), .groups = "drop")

speedup <- summary_tab |>
  filter(pixy_version == "2.0.0") |>
  left_join(baseline_1core, by = "statistic") |>
  mutate(
    n_cores_num = as.integer(n_cores),
    speedup     = baseline_s / median_s
  )

##########
# fig 1: wall-clock time by arm
##########
p1 <- ok |>
  ggplot(aes(x = arm_label, y = elapsed_s_clean, fill = pixy_version)) +
  geom_boxplot(outlier.size = 0.4, linewidth = 0.4) +
  facet_wrap(~statistic, scales = "free_y", ncol = 3) +
  scale_fill_manual(
    values = c("0.93.1" = "#d95f02", "2.0.0" = "#1b9e77"),
    labels = c("0.93.1" = "pixy 0.93.1", "2.0.0" = "pixy 2.0.0")
  ) +
  labs(
    title = "pixy benchmark: 10 Mb VCF, 25 kb windows, 100 replicates",
    x     = NULL,
    y     = "Wall-clock time (s)",
    fill  = "Version"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

ggsave(file.path(fig_dir, "fig1_runtime_by_arm.pdf"), p1, width = 10, height = 5)
message("Saved fig1_runtime_by_arm.pdf")

##########
# fig 2: multicore speedup (2.0.0 only)
##########
p2 <- speedup |>
  ggplot(aes(x = n_cores_num, y = speedup, colour = statistic, group = statistic)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  scale_x_continuous(breaks = c(1, 2, 4, 8, 16)) +
  scale_y_continuous(breaks = c(1, 2, 4, 8, 16)) +
  scale_colour_brewer(palette = "Dark2") +
  labs(
    title    = "pixy 2.0.0 multicore speedup",
    subtitle = "Relative to 1-core median; dashed line = ideal linear scaling",
    x        = "Cores",
    y        = "Speedup (×)",
    colour   = "Statistic"
  ) +
  theme_bw(base_size = 11)

ggsave(file.path(fig_dir, "fig2_speedup_curve.pdf"), p2, width = 6, height = 4)
message("Saved fig2_speedup_curve.pdf")

##########
# fig 3: speedup vs old pixy baseline
##########
old_baseline <- ok |>
  filter(pixy_version == "0.93.1") |>
  group_by(statistic) |>
  summarise(old_median_s = median(elapsed_s_clean), .groups = "drop")

speedup_vs_old <- summary_tab |>
  filter(pixy_version == "2.0.0") |>
  left_join(old_baseline, by = "statistic") |>
  mutate(
    n_cores_num = as.integer(n_cores),
    speedup     = old_median_s / median_s
  )

p3 <- speedup_vs_old |>
  ggplot(aes(x = n_cores_num, y = speedup, colour = statistic, group = statistic)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey60") +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = c(1, 2, 4, 8, 16)) +
  scale_colour_brewer(palette = "Dark2") +
  labs(
    title    = "pixy 2.0.0 speedup vs. pixy 0.93.1",
    subtitle = "Dashed line = parity with old single-core pixy",
    x        = "pixy 2.0.0 cores",
    y        = "Speedup over pixy 0.93.1 (×)",
    colour   = "Statistic"
  ) +
  theme_bw(base_size = 11)

ggsave(file.path(fig_dir, "fig3_speedup_vs_old.pdf"), p3, width = 6, height = 4)
message("Saved fig3_speedup_vs_old.pdf")

##########
# fig 4: peak RSS by arm
##########
p4 <- ok |>
  filter(!is.na(rss_kb_clean), rss_kb_clean > 0) |>
  mutate(rss_mb = rss_kb_clean / 1024) |>
  ggplot(aes(x = arm_label, y = rss_mb, fill = pixy_version)) +
  geom_boxplot(outlier.size = 0.4, linewidth = 0.4) +
  facet_wrap(~statistic, ncol = 3) +
  scale_fill_manual(
    values = c("0.93.1" = "#d95f02", "2.0.0" = "#1b9e77"),
    labels = c("0.93.1" = "pixy 0.93.1", "2.0.0" = "pixy 2.0.0")
  ) +
  labs(
    title = "Peak resident memory by arm",
    x     = NULL,
    y     = "Peak RSS (MB)",
    fill  = "Version"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

ggsave(file.path(fig_dir, "fig4_peak_rss_by_arm.pdf"), p4, width = 10, height = 5)
message("Saved fig4_peak_rss_by_arm.pdf")

message("\nAll figures written to: ", fig_dir)
