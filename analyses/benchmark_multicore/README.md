# benchmark_multicore → Figure 2

Multicore runtime and peak-memory benchmark of pixy 2.x against the previous release,
measured on simulated 10 Mb diploid VCFs at 1, 2, 4, 8, and 16 cores. Produces
`data/results/all_cells_long.tsv`, the input to `../../figures/Figure_2_Performance.R`.

## Overview

Each statistic (π, d<sub>xy</sub>, F<sub>ST</sub>) is timed for old pixy (single core)
and for pixy 2.x at each core count, over replicate VCFs. `04_aggregate_summaries.sh`
collates the per-run wall-clock and peak RSS into one long table; the figure then
reports speed relative to the old single-core baseline and peak memory as a fraction
of it.

## Scripts

Scripts are numbered in execution order. The `.sbatch`/`.sh` scripts run under SLURM
on the UCR HPCC; the R script runs natively.

| Script | Purpose |
|--------|---------|
| `00_create_environments.sh` | Build the conda envs in `envs/` (`pixy_old`, `pixy_new`, `vcfsim`) |
| `01_simulate_vcfs.sbatch` | Simulate the 10 Mb diploid VCFs with vcfsim |
| `02_pixy_old_array.sbatch` | Run old pixy, single core — the baseline |
| `03_pixy_new_array.sbatch` | Run pixy 2.x at 1/2/4/8/16 cores (submit with `03_submit_pixy_new_all.sh`) |
| `03b_submit_pixy_new_fst_only.sh` | Resubmit the F<sub>ST</sub> arm alone (wraps `03_submit_pixy_new_all.sh`) |
| `04_aggregate_summaries.sh` | Collate `pixy_old_*` + `pixy_new_*_cores_*` into `data/results/all_cells_long.tsv` |
| `05_memory_aggregate.sbatch` | Peak memory of the whole process tree at one core count (submit with `05_submit_memory_aggregate.sh`) |
| `06_cleanup_vcfs.sbatch` | Remove the simulated VCFs |
| `analysis/01_v2_benchmark.R` | Standalone mirror plot of the same data (diagnostic; the manuscript figure is `figures/Figure_2_Performance.R`) |

`config/pop_map_{1,2}pop.tsv` are the pixy population maps.

## Memory accounting

Two different quantities are recorded, and they are not interchangeable:

- `02`/`03` use `/usr/bin/time -f "%M"`, which reports the **largest single
  process** in the tree, not the sum over concurrently running workers. This is
  the per-process peak RSS.
- `05_memory_aggregate.sbatch` samples the whole process group and records the
  **peak sum across all live processes** (`Tree_peak_kb`) alongside the largest
  single process seen (`Proc_peak_kb`), so the two can be compared directly.

Claims about total job memory at a given core count must come from
`Tree_peak_kb`, not from the `%M` column.

## Software versions

Every arm records the pixy version in the aggregated table. pixy releases before
2.2.3 hardcoded a version string that drifted from the packaged version — for
example the conda 0.95.01 build prints "0.95.0", and 2.2.2 prints "2.2.1". Trust
the conda/pip package version in `envs/*.yml`, not `pixy --version`, for any
release before 2.2.3.

## Dependencies

Conda envs in `envs/`: `pixy_old`, `pixy_new`, `vcfsim`. The simulated VCFs are not
tracked; regenerate them with `01_simulate_vcfs.sbatch`.
