# figure s2 — diploid old vs new agreement per VCF (rows = pi/dxy/WC FST, cols = missingness)
# per-seed value = mean across that replicate's windows; vcftools excluded
# WC-FST "new" from the wcfst_dip_2pop arms, merged into dxy_dip_2pop (merge_wcfst.py)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(ggh4x)
})

script_dir <- tryCatch(
  dirname(rstudioapi::getActiveDocumentContext()$path),
  error = function(e) getwd()
)
here <- function(...) file.path(script_dir, ...)

source(here("figure_theme.R"))

POLY_DIR <- here("..", "analyses", "polyploid_multiallele", "analysis")
fig_path <- here("figs", "FigureS2_diploid_oldnew.pdf")

##########
# load per-seed values
##########
dat <- read_tsv(file.path(POLY_DIR, "per_seed_diploid.tsv"),
                show_col_types = FALSE)

if (nrow(dat) == 0) {
  stop("per_seed_diploid.tsv is empty. Run analysis/per_seed_diploid.R after ",
       "the WC-FST band-aid run and merge_wcfst.py have landed.")
}

##########
# build old-vs-new pairs (drop vcftools)
##########
stat_label <- c(pi = "pi", dxy = "d[xy]", fst_wc = "F[ST]~(WC)")

wide <- dat |>
  filter(estimator %in% c("old", "new")) |>
  select(stat, miss_pct, replicate, estimator, value) |>
  pivot_wider(names_from = estimator, values_from = value) |>
  filter(!is.na(old), !is.na(new))

# row facet: stat (parsed label). col facet: missingness
wide <- wide |>
  mutate(
    stat = factor(stat, levels = c("pi", "dxy", "fst_wc")),
    miss = factor(scales::percent(miss_pct, accuracy = 1),
                  levels = scales::percent(sort(unique(miss_pct)), accuracy = 1))
  )

# per-panel old-vs-new R^2, annotated in each cell so no regression is shown
# cell by cell (pi/dxy are integer ratios and match to machine precision; WC-FST
# is a floating-point ratio, so its R^2 is 1.00 to two decimals, not identical)
r2_lab <- wide |>
  group_by(stat, miss) |>
  summarise(r2 = cor(old, new)^2, .groups = "drop") |>
  # quote the value so plotmath prints "1.00" literally; unquoted it is parsed as
  # the number 1.00 and rendered as "1"
  mutate(label = sprintf('R^2 == "%.2f"', r2))

# shared x limits per stat (old and new are same units, so this fits y too);
# y handled by free_y below
stat_levels <- levels(wide$stat)
row_lims <- wide |>
  group_by(stat) |>
  summarise(lo = min(old, new), hi = max(old, new), .groups = "drop") |>
  arrange(match(stat, stat_levels))
pad <- function(lo, hi) { d <- (hi - lo) * 0.05; c(lo - d, hi + d) }
x_scales <- Map(function(lo, hi) scale_x_continuous(limits = pad(lo, hi)),
                row_lims$lo, row_lims$hi)

##########
# plot
##########
# facet_grid shares x down columns, so use ggh4x to free x per panel and pin
# each stat-row to a shared x range via facetted_pos_scales. y is free_y
fig <- ggplot(wide, aes(old, new)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey30", linewidth = 0.4) +
  geom_point(size = 0.5, alpha = 0.35, colour = "#3C3489") +
  geom_text(data = r2_lab, aes(x = -Inf, y = Inf, label = label),
            parse = TRUE, inherit.aes = FALSE,
            hjust = -0.18, vjust = 1.6, size = 2.5, colour = "grey25") +
  facet_grid2(stat ~ miss, scales = "free", independent = "x",
              labeller = labeller(stat = as_labeller(stat_label, label_parsed))) +
  facetted_pos_scales(x = rep(x_scales, each = nlevels(wide$miss))) +
  # 0.95.01 is the legacy comparator actually installed; it is code-identical to
  # the 0.95.02 of Korunes & Samuk 2021 (that diff touches only README and docs).
  labs(x = "pixy 0.95.01 estimate",
       y = "pixy 2.2.3 estimate") +
  theme_pixy() +
	theme(axis.text = element_text(size = 10))+
  theme(aspect.ratio = 1,
        strip.text.x = element_text(face = "plain"))

fig

pixy_save(fig, fig_path, width = 11, height = 7)
