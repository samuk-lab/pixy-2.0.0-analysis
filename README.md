# pixy 2.0.0 — analysis code

Data-generation and figure code for the pixy 2.0.0 manuscript. pixy 2.0.0 extends
[pixy](https://github.com/ksamuk/pixy) — an estimator of population-genetic summary
statistics from VCFs — with multicore support, polyploid genotypes, and
multiallelic-aware estimators of π, d<sub>xy</sub>, Watterson's θ, Tajima's D, and
F<sub>ST</sub>. This repository contains everything needed to reproduce the paper's
five main figures, two supplementary figures, and its tables. It does not contain
the manuscript itself.

> Companion analysis code for the **pixy 2.0.0** manuscript (Samuk et al., full
> citation forthcoming).

---

## Overview

The paper validates the new features against theory, simulation, and real data.
Each figure is drawn by one script in `figures/`, which reads a small summary table
produced by a pipeline in `analyses/`:

- **Multicore performance (Figure 2)** — runtime and peak memory of pixy 2.x versus
  the previous release across 1–16 cores.
- **Polyploid estimators (Figures 3 and S2)** — bias of π, d<sub>xy</sub>, θ<sub>W</sub>,
  and Tajima's D across four ploidy levels and five missingness levels, compared to
  finite-sites expectations, plus a diploid cross-version check of old pixy against 2.0.0.
- **Multiallelic-aware estimators (Figure 4)** — biallelic versus multiallelic-aware
  π, d<sub>xy</sub>, and Hudson's F<sub>ST</sub> across a θ sweep.
- **Empirical demonstration (Figure 5)** — the `--include_multiallelic_snps` flag on
  real data from *Anopheles gambiae* (2n) and *Arabidopsis arenosa* (4n).
- **Citation footprint (Figure S1)** — counts of papers citing pixy, from OpenAlex.

The finite-sites Tajima's D reference line in Figure 3 is calibrated from the theory
in `analyses/tajima_multiallelic_variance/`.

Each figure and its inputs:

| Figure / table | Script (`figures/`) | Input |
|---|---|---|
| Figure 1 (overview) | `Figure_1_pixy_overview.sh` | none — hand-built SVG |
| Figure 2 (performance) | `Figure_2_Performance.R` | `benchmark_multicore/data/results/all_cells_long.tsv` |
| Figure 3 (polyploid) | `Figure_3_Polyploid.R` | `polyploid_multiallele/analysis/{missingness_summary,thetaw_tajimaD_summary,per_rep_*}.tsv` |
| Figure 4 (multiallelic) | `Figure_4_Multiallelic.R` | `polyploid_multiallele/analysis/{multiallelic_frac_summary,dxy_fst_theta_summary}.tsv` |
| Figure 5 (empirical) | `Figure_5_Empirical.R` | `empirical-tests/data/aggregated/{anopheles,arenosa}.joined.tsv`, `config/arms.tsv` |
| Figure S1 (citations) | `Figure_S1_Citations.R` | `citation_network/pixy_citer_topic_network_output/{pixy_citations_per_year,subject_counts}.csv` |
| Figure S2 (diploid old/new) | `Figure_S2_Diploid_OldNew.R` | `polyploid_multiallele/analysis/per_seed_diploid.tsv` |
| In-text stats / Table S2 | `Inline_statistical_tests.R` | the summaries above plus the raw per-arm TSVs |

---

## Repository structure

```
github_analysis/
├── figures/                            # one script per manuscript figure
│   ├── Figure_1_pixy_overview.sh          # Fig 1 schematic (hand-built SVG)
│   ├── Figure_2_Performance.R             # Fig 2  (reads benchmark_multicore)
│   ├── Figure_3_Polyploid.R               # Fig 3  (reads polyploid_multiallele)
│   ├── Figure_4_Multiallelic.R            # Fig 4  (reads polyploid_multiallele)
│   ├── Figure_5_Empirical.R               # Fig 5  (reads empirical-tests)
│   ├── Figure_S1_Citations.R              # Fig S1 (reads citation_network)
│   ├── Figure_S2_Diploid_OldNew.R         # Fig S2 (reads polyploid_multiallele)
│   ├── Inline_statistical_tests.R         # in-text statistics + Table S2
│   ├── figure_theme.R                     # shared ggplot2 theme + palettes
│   └── figs/                              # rendered figure PDFs + inline_stats.md
├── analyses/
│   ├── benchmark_multicore/            # multicore runtime + memory       → Figure 2
│   ├── polyploid_multiallele/          # ploidy × missingness × θ sims     → Figures 3, 4, S2
│   ├── empirical-tests/                # Anopheles + A. arenosa real data  → Figure 5
│   ├── citation_network/              # OpenAlex citation pull            → Figure S1
│   └── tajima_multiallelic_variance/  # finite-sites Tajima's D theory    (Fig 3 reference line)
├── README.md
└── LICENSE
```

Each `analyses/` subfolder has its own README describing the pipeline that
regenerates its summary tables.

---

## Reproducing the figures

The summary tables are included, so all seven figures render without rerunning the
pipelines:

```sh
cd figures
bash Figure_1_pixy_overview.sh
for f in Figure_2_Performance Figure_3_Polyploid Figure_4_Multiallelic \
         Figure_5_Empirical Figure_S1_Citations Figure_S2_Diploid_OldNew; do
  Rscript $f.R
done
```

`Inline_statistical_tests.R` rebuilds the in-text statistics and Table S2. It is the
one script that also needs the large raw per-arm data (not included — see below), so
its output `figures/figs/inline_stats.md` is provided precomputed.

The `analyses/` pipelines were written for the UCR HPCC (SLURM + conda), and their
`.sbatch` scripts carry absolute cluster paths, so they document the analysis rather
than run turnkey elsewhere. Large intermediate data are not included, since they are
bulky and regeneratable: the simulated VCFs, the 88 aggregated per-arm TSVs (~380 MB)
that `Inline_statistical_tests.R` reads, and the per-replicate cells under
`tajima_multiallelic_variance/data/` (~4 GB). Run the matching pipeline to rebuild them.

---

## Dependencies

### R (figures)

R 4.4.2, with `tidyverse`, `patchwork`, `ggh4x`, `RColorBrewer`, and `scales`. PDF
output uses `cairo_pdf`. See `figures/README.md`.

### Pipelines

Each analysis builds its own conda environments from the `envs/*.yml` files in its
folder (pixy, vcfsim, vcftools, GATK, and so on). See the per-module READMEs.

---

## License

MIT — see `LICENSE`.
