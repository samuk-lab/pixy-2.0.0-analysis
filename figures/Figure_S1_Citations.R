# figure s1 — pixy citation footprint: citing papers per year + openalex subfields

suppressPackageStartupMessages({
  library(dplyr)
  library(forcats)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

script_dir <- tryCatch(
  dirname(rstudioapi::getActiveDocumentContext()$path),
  error = function(e) getwd()
)
here <- function(...) file.path(script_dir, ...)

source(here("figure_theme.R"))

fig_dir  <- here("figs")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
fig_path <- file.path(fig_dir, "FigureS1_citations.pdf")

##########
# panel a: citing papers per year
##########
per_year <- read_csv(
  here("..", "analyses", "citation_network",
       "pixy_citer_topic_network_output", "pixy_citations_per_year.csv"),
  show_col_types = FALSE
)

cum_scale <- max(per_year$citations_added) / max(per_year$cumulative)

p_year <- ggplot(per_year, aes(x = year)) +
  geom_col(aes(y = citations_added), fill = "grey45", width = 0.75) +
  geom_line(aes(y = cumulative * cum_scale),
            colour = "#A90264", linewidth = 0.8, group = 1) +
  geom_point(aes(y = cumulative * cum_scale),
             colour = "#A90264", size = 1.8) +
  scale_x_continuous(breaks = per_year$year) +
  scale_y_continuous(
    name     = "Citing papers (per year)",
    sec.axis = sec_axis(~ . / cum_scale, name = "Cumulative"),
    expand   = expansion(mult = c(0, 0.05))
  ) +
  labs(x = NULL, tag = "a") +
  theme_pixy() +
  theme(
    axis.title.y.right = element_text(colour = "#A90264"),
    axis.text.y.right  = element_text(colour = "#A90264"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

##########
# panel b: subject distribution
##########
subject_counts <- read_csv(
  here("..", "analyses", "citation_network",
       "pixy_citer_topic_network_output", "subject_counts.csv"),
  show_col_types = FALSE
)

subject_counts_other <- subject_counts |>
  mutate(subject_area = ifelse(n_citing_papers < 5,
                               "Other", as.character(subject_area))) |>
  group_by(subject_area) |>
  summarise(n_citing_papers = sum(n_citing_papers), .groups = "drop") |>
  mutate(subject_area = fct_reorder(subject_area, n_citing_papers),
         subject_area = fct_relevel(subject_area, "Other"))

# label placement: bars >=35% of the longest get white text inside the bar end,
# narrower ones dark text just past it. drops the y-axis entirely
xmax <- max(subject_counts_other$n_citing_papers)
subject_counts_other <- subject_counts_other |>
  mutate(
    label_inside  = n_citing_papers / xmax >= 0.35,
    label_x       = ifelse(label_inside,
                           n_citing_papers - xmax * 0.015,
                           n_citing_papers + xmax * 0.015),
    label_hjust   = ifelse(label_inside, 1, 0),
    label_colour  = ifelse(label_inside, "white", "grey20"),
    label_text    = as.character(subject_area)
  )

p_subj <- ggplot(subject_counts_other,
                 aes(x = n_citing_papers, y = subject_area)) +
  geom_col(fill = "grey45") +
  geom_text(aes(x = label_x, label = label_text,
                hjust = label_hjust, colour = label_colour),
            size = 2.8) +
  scale_colour_identity() +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(
    x = "Citing papers",
    y = NULL,
    tag = "b"
  ) +
  theme_pixy() +
  theme(
    axis.text.y       = element_blank(),
    axis.ticks.y      = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )

##########
# compose
##########
fig_s1 <- p_year | p_subj

fig_s1

pixy_save(fig_s1, fig_path, width = 8, height = 4)
