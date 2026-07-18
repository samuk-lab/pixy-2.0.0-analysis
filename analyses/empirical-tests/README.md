# empirical-tests → Figure 5

Demonstration of pixy 2.0.0's `--include_multiallelic_snps` flag on real published
data, at two ploidies. For each species, per-window π, d<sub>xy</sub>, and F<sub>ST</sub>
are computed with and without the flag, in two populations. Produces
`data/aggregated/{anopheles,arenosa}.joined.tsv`, read by
`../../figures/Figure_5_Empirical.R` and keyed by `config/arms.tsv`.

## Overview

Two whole autosomes are called from scratch — from public BAMs/FASTQs through GATK
HaplotypeCaller to raw all-sites VCFs — then filtered (indel + near-indel-SNP removal,
GATK hard filters incl. ExcessHet, a genmap mappability mask, a quantile depth band, and
a per-population max-missing threshold) and run through pixy twice per arm, once with the
multiallelic flag and once without. *Anopheles* is subsampled to 8 individuals per population
(seeded; `config/pops_anopheles.tsv`) to match *arenosa*'s denominator; the full 25/pop list
is kept in `config/pops_anopheles.full.tsv`.

All site filtering is symmetric over variant and invariant sites (indel/mask/depth-band/
max-missing drop whole sites; genotype depth failures become missing), so pixy's per-window
denominators stay correct. Filter thresholds are env-var overridable at the top of
`05c_filter.sbatch` (defaults: genotype `DP>=5`; drop indels + SNPs within 5 bp of an indel,
keep SNP + invariant only; GATK SNP hard filters QD/FS/MQ/MQRankSum/ReadPosRankSum/SOR +
`ExcessHet>54.69`; genmap `K=100,E=1` mappability `==1`; site `INFO/DP` within the 5–95%
quantiles; `>=6/8` genotypes present per population).

| Arm | Species | Ploidy | Populations | n/pop | Reference | Chromosome |
|---|---|---|---|---|---|---|
| `anopheles` | *A. gambiae* | 2 | BFS, KES | 8 | AgamP4 | 3R |
| `arenosa` | *A. arenosa* | 4 | SPI, TRE | 8 | AARE701a | LR999451.1 |

## Scripts

Numbered in execution order; all run under SLURM on the UCR HPCC.

| Script | Purpose |
|--------|---------|
| `00_create_environments.sh` | Build the conda envs in `envs/` (`fetch`, `gatk`, `pixy`, `mappability`) |
| `01_pick_regions.sh` | Build single-chromosome reference FASTAs and target BEDs |
| `02_fetch_anopheles.sbatch` | Stream chr 3R BAM slices for the *Anopheles* samples from ENA |
| `03_call_anopheles.sbatch` | GATK HaplotypeCaller (`--sample-ploidy 2`) → raw all-sites VCF |
| `04_fetch_arenosa.sbatch` | FASTQ → BWA-MEM2 → mark duplicates for the *arenosa* samples |
| `05_call_arenosa.sbatch` | GATK HaplotypeCaller (`--sample-ploidy 4`) → raw all-sites VCF |
| `05b_mappability.sbatch` | genmap `K=100,E=1` mappability mask per arm → `config/regions/<arm>.mappable.bed` |
| `05c_filter.sbatch` | raw VCF → hard filters + mask + depth band + max-missing → canonical `<arm>.all_sites.vcf.gz` |
| `06_run_pixy.sbatch` | pixy per arm × {biallelic, +multiallelic}, 10 kb windows |
| `07_compare.sbatch` | Join the two pixy runs per arm → `data/aggregated/<arm>.joined.tsv` |
| `count_multi.sbatch` | Count multiallelic SNP sites (`bcftools -m3 -v snps`) per 10 kb window → `data/aggregated/<arm>.multi_per_window.tsv` |

`03`/`05` emit only the raw multi-sample VCF (INFO/FORMAT annotations retained);
`05c_filter.sbatch` is the sole producer of the filtered `<arm>.all_sites.vcf.gz`
that `06` consumes. `05b` and `05c` depend on the per-arm reference and raw VCF, so
run them after `01` and `03`/`05`. `count_multi.sbatch` reads the same filtered VCF
and is independent of `06`/`07`; its per-window counts back the uplift-vs-multiallelic-
density Spearman ρ reported in `../../figures/Inline_statistical_tests.R`. `lib/common.sh`
holds shared bash helpers.

## Dependencies

Conda envs in `envs/`: `fetch`, `gatk`, `pixy`, `mappability`. Sequence data (BAMs,
VCFs) are not tracked; only the joined per-arm TSVs and the per-window multiallelic
counts (`<arm>.multi_per_window.tsv`) are kept.
