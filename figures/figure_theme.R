# shared ggplot2 theme + palettes for the manuscript figures
# sourced by each figure script; edit here to restyle all figures at once

suppressPackageStartupMessages({
  library(ggplot2)
  library(RColorBrewer)
})

##########
# global sizing
##########
PIXY_BASE_SIZE <- 14

##########
# global font
##########
# Lato (system-wide); ships the Greek glyphs (pi / theta / mu) for the plotmath
# labels, so no separate math font. resolved by name via cairo_pdf
PIXY_FONT <- "Lato"

##########
# theme
##########
theme_pixy <- function(base_size = PIXY_BASE_SIZE) {
  theme_bw(base_size = base_size, base_family = PIXY_FONT) +
    theme(
      panel.grid.minor = element_line(linewidth = 0.25, colour = "grey90"),
      panel.grid.major = element_line(linewidth = 0.25, colour = "grey90"),
      panel.border     = element_rect(linewidth = 0.4, colour = "grey40"),
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold", size = base_size),
      plot.title       = element_blank(),
      plot.subtitle    = element_blank(),
      plot.tag         = element_text(face = "bold", size = base_size + 2),
      legend.position  = "top",
      legend.title     = element_text(size = base_size - 1),
      legend.text      = element_text(size = base_size - 1),
      legend.key.size  = unit(0.7, "lines"),
      axis.title       = element_text(size = base_size),
      axis.text        = element_text(size = base_size - 1)
    )
}

##########
# font on text geoms
##########
# text-drawing geoms (geom_text/label, ggrepel) default to the device font, not
# the theme, so set their family here for in-panel labels. mutates geom defaults
# for the session on source
update_geom_defaults("text",  list(family = PIXY_FONT))
update_geom_defaults("label", list(family = PIXY_FONT))
if (requireNamespace("ggrepel", quietly = TRUE)) {
  try(update_geom_defaults("text_repel",  list(family = PIXY_FONT)), silent = TRUE)
  try(update_geom_defaults("label_repel", list(family = PIXY_FONT)), silent = TRUE)
}

##########
# palettes (RColorBrewer)
##########
# estimator (pixy variant) — Dark2
.dark2 <- brewer.pal(8, "Set1")
pixy_estimator_colors <- c(
  new        = .dark2[1],
  new_multi  = .dark2[2],
  old        = .dark2[3],
  vcftools   = .dark2[4],
  theory     = "grey20"
)
pixy_estimator_labels <- c(
  new        = "pixy 2.0.0 (biallelic)",
  new_multi  = "pixy 2.0.0 (multiallelic)",
  old        = "pixy 0.95.02",
  vcftools   = "vcftools",
  theory     = "expectation"
)

# statistic (pi / dxy / FST / thetaw / tajimaD) — Set1
.set1 <- brewer.pal(8, "Set1")
pixy_statistic_colors <- c(
  pi      = "#0F6F57",
  dxy     = "#3C3489",
  fst     = "#A90264",
  thetaw  = .set1[3],
  tajimaD = .set1[4]
)
pixy_statistic_labels <- c(
  pi      = "pi",
  dxy     = "d[xy]",
  fst     = "F[ST]",
  thetaw  = "theta[W]",
  tajimaD = "D[Tajima]"
)

# ploidy (2 / 4 / 6 / 8) — sequential YlGnBu, ramps up with ploidy
.ylgnbu <- brewer.pal(9, "YlGnBu")
pixy_ploidy_colors <- c(
  `2` = .ylgnbu[3],
  `4` = .ylgnbu[5],
  `6` = .ylgnbu[7],
  `8` = .ylgnbu[9]
)

# pixy version — Set2
.set2 <- brewer.pal(8, "Set2")
pixy_version_colors <- c(
  `0.93.1` = .set2[2],
  `2.0.0`  = .set2[1]
)

# FST flavour
pixy_fst_flavour_colors <- c(
  hudson = .dark2[1],
  wc     = .dark2[2]
)

##########
# convenience scales
##########
scale_colour_pixy_estimator <- function(...) {
  scale_colour_manual(values = pixy_estimator_colors,
                      labels = pixy_estimator_labels, name = NULL, ...)
}
scale_fill_pixy_estimator <- function(...) {
  scale_fill_manual(values = pixy_estimator_colors,
                    labels = pixy_estimator_labels, name = NULL, ...)
}
scale_colour_pixy_statistic <- function(...) {
  scale_colour_manual(
    values = pixy_statistic_colors,
    labels = function(x) parse(text = pixy_statistic_labels[x]),
    name   = "Statistic",
    ...
  )
}
scale_colour_pixy_ploidy <- function(...) {
  scale_colour_manual(values = pixy_ploidy_colors, name = "Ploidy", ...)
}
scale_fill_pixy_version <- function(...) {
  scale_fill_manual(values = pixy_version_colors,
                    labels = c(`0.93.1` = "pixy 0.95.02", `2.0.0` = "pixy 2.0.0"),
                    name = "Version", ...)
}

##########
# output helper
##########
pixy_save <- function(plot, path, width = 7, height = 5) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  ggsave(path, plot, width = width, height = height, device = cairo_pdf)
  message("Saved ", path)
}
