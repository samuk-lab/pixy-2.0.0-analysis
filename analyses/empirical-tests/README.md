# empirical-tests → Figure 5

Demonstration of pixy 2.0.0's `--include_multiallelic_snps` flag on real published
data, at two ploidies. For each species, per-window π, d<sub>xy</sub>, and F<sub>ST</sub>
are computed with and without the flag, in two populations. Produces
`data/aggregated/{anopheles,arenosa}.joined.tsv`, read by
`../../figures/Figure_5_Empirical.R` and keyed by `config/arms.tsv`.

## Overview

Two whole autosomes are called from scratch — from public BAMs/FASTQs through GATK
HaplotypeCaller to all-sites VCFs — then run through pixy twice per arm, once with the
multiallelic flag and once without. *Anopheles* is subsampled to 8 individuals per
population (seeded; `config/pops_anopheles.tsv`) to match *arenosa*'s denominator; the
full 25/pop list is kept in `config/pops_anopheles.full.tsv`.

| Arm | Species | Ploidy | Populations | n/pop | Reference | Chromosome |
|---|---|---|---|---|---|---|
| `anopheles` | *A. gambiae* | 2 | BFS, KES | 8 | AgamP4 | 3R |
| `arenosa` | *A. arenosa* | 4 | SPI, TRE | 8 | AARE701a | LR999451.1 |

## Scripts

Numbered in execution order; all run under SLURM on the UCR HPCC.

| Script | Purpose |
|--------|---------|
| `00_create_environments.sh` | Build the conda envs in `envs/` (`fetch`, `gatk`, `pixy`) |
| `01_pick_regions.sh` | Build single-chromosome reference FASTAs and target BEDs |
| `02_fetch_anopheles.sbatch` | Stream chr 3R BAM slices for the *Anopheles* samples from ENA |
| `03_call_anopheles.sbatch` | GATK HaplotypeCaller (`--sample-ploidy 2`) → all-sites VCF |
| `04_fetch_arenosa.sbatch` | FASTQ → BWA-MEM2 → mark duplicates for the *arenosa* samples |
| `05_call_arenosa.sbatch` | GATK HaplotypeCaller (`--sample-ploidy 4`) → all-sites VCF |
| `06_run_pixy.sbatch` | pixy per arm × {biallelic, +multiallelic}, 10 kb windows |
| `07_compare.sbatch` | Join the two pixy runs per arm → `data/aggregated/<arm>.joined.tsv` |

`lib/common.sh` holds shared bash helpers used across the pipeline.

## Dependencies

Conda envs in `envs/`: `fetch`, `gatk`, `pixy`. Sequence data (BAMs, VCFs) are not
tracked; only the joined per-arm TSVs are kept.
