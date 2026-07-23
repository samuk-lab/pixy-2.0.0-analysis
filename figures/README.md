# figures

The manuscript figure scripts for pixy 2.x. Each script renders one figure into
`figs/`, reading a summary table produced by a pipeline in `../analyses/`. Open
`figures.Rproj` in RStudio (the working directory must be this folder) and run any
script; every script resolves its inputs relative to this folder (`../analyses/...`)
and sources `figure_theme.R` for the shared theme and palettes.

## Scripts

| Script | Output | Reads |
|--------|--------|-------|
| `Figure_1_pixy_overview.sh` | `figs/Figure1_pixy_overview.pdf` | hand-built SVG, no data |
| `Figure_2_Performance.R` | `figs/Figure2_performance.pdf` | `../analyses/benchmark_multicore/data/results/all_cells_long.tsv` (panel a), `.../pixy_mem_*.tsv` (panel b) |
| `Figure_3_Polyploid.R` | `figs/Figure3_polyploid.pdf` | `../analyses/polyploid_multiallele/analysis/{missingness_summary,thetaw_tajimaD_summary,per_rep_*}.tsv` |
| `Figure_4_Multiallelic.R` | `figs/Figure4_multiallelic.pdf` | `../analyses/polyploid_multiallele/analysis/{multiallelic_frac_summary,dxy_fst_theta_summary}.tsv` |
| `Figure_5_Empirical.R` | `figs/Figure5_empirical.pdf` | `../analyses/empirical-tests/data/aggregated/{anopheles,arenosa}.joined.tsv`, `config/arms.tsv` |
| `Figure_S1_Citations.R` | `figs/FigureS1_citations.pdf` | `../analyses/citation_network/pixy_citer_topic_network_output/{pixy_citations_per_year,subject_counts}.csv` |
| `Figure_S2_Diploid_OldNew.R` | `figs/FigureS2_diploid_oldnew.pdf` | `../analyses/polyploid_multiallele/analysis/per_seed_diploid.tsv` |
| `Inline_statistical_tests.R` | `figs/inline_stats.md` | the raw arms in `../analyses/polyploid_multiallele/data/aggregated/`, the summaries above, and both empirical joined TSVs |
| `figure_theme.R` | — | sourced by every figure script (theme + palettes) |

## Reproduce

```sh
cd figures
bash Figure_1_pixy_overview.sh
for f in Figure_2_Performance Figure_3_Polyploid Figure_4_Multiallelic \
         Figure_5_Empirical Figure_S1_Citations Figure_S2_Diploid_OldNew; do
  Rscript $f.R
done
```

The `*_summary.tsv` inputs are produced by `../analyses/polyploid_multiallele`
(`04_run_R_analyses.sh`); `all_cells_long.tsv` by
`../analyses/benchmark_multicore/04_aggregate_summaries.sh`; the citation CSVs by
`../analyses/citation_network/build_pixy_citation_data.R`.

`Inline_statistical_tests.R` needs the raw per-arm data (not included — see the
top-level README) and writes the numbers reported in the Results. Its output is
provided precomputed as `figs/inline_stats.md`.

## Figure 1

`Figure_1_pixy_overview.sh` regenerates the pixy schematic from an SVG heredoc
embedded in the script (the single source of truth), writing `Figure_1_pixy_overview.svg`
and a PDF. Colours live in the `COLORS` block at the top of the script, and special
glyphs use XML numeric entities (`&#960;` for π). The SVG-to-PDF step prefers a native
converter (`rsvg-convert`, `inkscape`, or `cairosvg`) and falls back to Edge via WSL,
so no Linux packages are strictly required. The manuscript embeds an Affinity-designed
`figs/Figure1_pixy_overview.png` built from this schematic.

## Dependencies

`tidyverse`, `patchwork`, `ggh4x` (Figure S2 per-row scales), `RColorBrewer`, and
`scales`. PDF output uses `cairo_pdf`. Tested on R 4.4.2. Edit `figure_theme.R` to
change fonts, palettes, or panel borders, then rerun any script.
