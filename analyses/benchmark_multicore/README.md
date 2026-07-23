# benchmark_multicore → Figure 2

Multicore runtime and peak-memory benchmark of pixy 2.x against the previous release,
measured on simulated 10 Mb diploid VCFs at 1, 2, 4, 8, and 16 cores. Produces
`data/results/all_cells_long.tsv`, the input to `../../figures/Figure_2_Performance.R`.

## Overview

Each statistic (π, d<sub>xy</sub>, F<sub>ST</sub>) is timed for old pixy (single core)
and for pixy 2.x at each core count, over replicate VCFs. `04_aggregate_summaries.sh`
collates the per-run wall-clock into `all_cells_long.tsv` (Figure 2a, speed relative to
the old single-core baseline). Memory comes from the separate tree-sampling runs in
`pixy_mem_*.tsv` (Figure 2b, absolute peak memory of the whole process tree against the
old single-core footprint as a horizontal reference) — see *Memory accounting* below.

## Scripts

Scripts are numbered in execution order. The `.sbatch`/`.sh` scripts run under SLURM
on the UCR HPCC; the R script runs natively.

| Script | Purpose |
|--------|---------|
| `00_create_environments.sh` | Build the conda envs in `envs/` (`pixy_old`, `pixy_new`, `vcfsim`) |
| `00b_reset_previous_run.sh` | Wipe a previous run's summaries, working output and logs before rerunning (dry run unless `--yes`). Keeps the simulated VCFs |
| `01_simulate_vcfs.sbatch` | Simulate the 10 Mb diploid VCFs with vcfsim |
| `02b_pixy_blocked_array.sbatch` | **Timing, current method.** One task = one (statistic, seed) block, running 0.95.01 and pixy 2.x at 1/2/4/8/16 cores back to back on one node — see *Timing design* below |
| `02_pixy_old_array.sbatch` | Run old pixy, single core — the baseline. Superseded by `02b` for Figure 2a; kept for one-off single-cell reruns |
| `03_pixy_new_array.sbatch` | Run pixy 2.x at 1/2/4/8/16 cores (submit with `03_submit_pixy_new_all.sh`). Superseded by `02b` for Figure 2a; kept for one-off single-cell reruns |
| `03b_submit_pixy_new_fst_only.sh` | Resubmit the F<sub>ST</sub> arm alone (wraps `03_submit_pixy_new_all.sh`) |
| `04_aggregate_summaries.sh` | Collate `pixy_old_*` + `pixy_new_*_cores_*` into `data/results/all_cells_long.tsv` |
| `05_memory_aggregate.sbatch` | Peak memory of the whole process tree at one core count (submit with `05_submit_memory_aggregate.sh`) |
| `06_cleanup_vcfs.sbatch` | Remove the simulated VCFs |
| `analysis/01_v2_benchmark.R` | Standalone mirror plot of the same data (diagnostic; the manuscript figure is `figures/Figure_2_Performance.R`) |

`config/pop_map_{1,2}pop.tsv` are the pixy population maps.

## Timing design

Figure 2a is a ratio: 0.95.01 single-core time over pixy 2.x time at *n* cores. Run
one cell per array job, as `02`/`03` do, and the two sides of that ratio are measured
in different jobs at different times, under whatever else happened to be on the
cluster. That cannot be scheduled away — `/bigdata` is shared, and no amount of
self-throttling keeps other users off it.

`02b_pixy_blocked_array.sbatch` removes the problem three ways, in descending order
of how much each is worth:

1. **The shared filesystem leaves the timed path.** The VCF is copied to node-local
   `/scratch/$USER` and every run in the block reads it from there, with zarr
   stores, pixy output and logs all local. Only the summary append touches
   `/bigdata`. `/scratch` is a physical disk on each node (`/dev/sda6` on intel),
   not a shared mount — it is *not* visible from other nodes or from the login
   nodes, which is exactly why it works here and why nothing durable may be left
   in it. Blocks are written under the 0700 per-user directory because `/scratch`
   itself is world-writable with no sticky bit. Nothing purges it automatically,
   so the script clears its block dir on signals as well as on normal exit.
2. **Blocking.** All six measurements for one (statistic, seed) run back to back on
   one node, so each ratio is formed under one ambient load and that load cancels.
3. **Rotated order.** The run order rotates with the seed, so no run type is
   permanently first, where it would pay the cold-cache cost on behalf of the rest.

Because the runs are paired in time, the ambient load largely cancels however the
ratio is summarised. `Figure_2_Performance.R` currently takes the **ratio of per-cell
medians** (`median(0.95.01) / median(2.2.3 @ n)`). The per-seed pairing is preserved
in the summary files (recoverable by joining cells on `Seed`), so a paired
median-of-per-seed-ratios with a bootstrap CI can be substituted if panel a later
needs error bars — but that form is not what is currently plotted.

`02b` writes to the same `pixy_old_*` / `pixy_new_*_cores_*` summary files, with the
same headers, as `02`/`03`. `04_aggregate_summaries.sh` and the figure need no
changes, and the pairing is recovered by joining cells on `Seed`.

Sizing: `--cpus-per-task=16` for all six runs (the 1-core runs leave cores idle —
that is the cost of the design), `--array=1-300%24` for 3 statistics × 100 seeds,
where 24 × 16 = 384 is the `samuklab` group CPU cap. About 3.5 h end to end.

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

The memory arms run 20 seeds per (statistic × cores) cell; the 0.95.01 baseline in
`05b` runs 10, since it is single-process and its footprint is tight across seeds.

`05`/`05b` append to `data/results/pixy_mem_*.tsv` under `flock` and only write the
header when the file is absent (as do `02`/`03` for their own summaries), so a rerun
on top of an existing `data/results/` silently appends to the old rows. Run
`00b_reset_previous_run.sh --yes` first — it clears `data/results/`, the pixy and zarr
working directories, and `logs/`, while keeping the simulated VCFs.

## Simulated input

`01_simulate_vcfs.sbatch` runs each array task in its own `data/vcfs/dm_10Mb/_work/seed_N/`
and matches vcfsim's output by exact name. An earlier version ran every task in the
shared output directory and normalised with `ls "${prefix}"*.vcf | head -n1`; because
vcfsim appends the seed to `--output_file`, seed 11's glob matched seed 1's unrenamed
intermediate, and `dm_sim_vcf_seed_11.vcf.gz` was built from seed 1's simulation. The
two were byte-identical, leaving 99 unique replicates out of 100.

The script **skips seeds whose `.vcf.gz` + `.tbi` already exist**, so affected files
never self-correct. `00b_reset_previous_run.sh` deletes
`dm_sim_vcf_seed_{1,11}.vcf.gz{,.tbi}` and any stray `*.vcf` intermediates, leaving
the other 98 VCFs in place; `01` then backfills only the two deleted seeds.

## Partition

Every measurement job is pinned to `intel`. Timings and memory peaks are not
comparable across heterogeneous hardware, so nothing in this benchmark may be
allowed to land on `short` (or any other partition) — including reruns of single
arms. Only `06_cleanup_vcfs.sbatch`, which measures nothing, is unpinned.

## Software versions

Every arm records the pixy version in the aggregated table. pixy releases before
2.2.3 hardcoded a version string that drifted from the packaged version — for
example the conda 0.95.01 build prints "0.95.0", and 2.2.2 prints "2.2.1". Trust
the conda/pip package version in `envs/*.yml`, not `pixy --version`, for any
release before 2.2.3.

## Dependencies

Conda envs in `envs/`: `pixy_old`, `pixy_new`, `vcfsim`. The simulated VCFs are not
tracked; regenerate them with `01_simulate_vcfs.sbatch`.
