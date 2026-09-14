# Threading pass — results

Date: 2026-09-13  
Grid: 200×120 (default). Machine reported `hardware_concurrency() = 8`.

## Policy

- **Auto** uses **1 worker** on the default 200×120 grid. Red/black pressure barriers cost more than they save at this size.
- Manual **2 / 4 / 6 / 8** still create a persistent pool (capped to the machine). Pressure stays serial unless a red/black list has **≥ 20 000 cells** (larger custom grids).
- `--benchmark`, `--liquid-diag`, and `--rigid-benchmark` force **1 thread**.

## What was parallelized

Pressure red phase → barrier → pressure black phase, using existing `pressureRed` / `pressureBlack` lists and stencil cache. Direct sequential SOR when the pool is empty or the list is below the threshold.

## Left single-threaded

Volume transport, residual consolidation, splash deposit, rigid contacts, surface/advection/chunk metrics. Transport and residual share neighbor writes; conservation wins.

## 200×120 timings (physics ms / pressure ms)

| Scene | 1 thread | 2 thread (pressure on) | 4 thread (pressure on) |
|---|---:|---:|---:|
| A light droplet | 1.5–1.8 / 0.5 | ~1.4 / 0.5 (serial) | ~1.6 / 0.5 (serial) |
| B puddle | 1.3 / 0.45 | 1.3 / 0.45 (serial) | 1.3 / 0.45 (serial) |
| C large pool | **10.3 / 3.1** | **13.4 / 4.6** (parallel, slower) | **17.4 / 6.6** (slower) |
| D dam break | **8.5 / 2.5** | **9.9 / 3.4** (slower) | **20.3 / 8.6** (slower) |
| D2 falling stream | 7.1 / 2.0 | 6.9 / 2.1 (serial) | 6.9 / 2.1 (serial) |

Coarse one-range-per-worker grain did not recover a win. Threshold raised to 20 000 so default scenes stay serial.

## Conservation

Volume error matched across 1/2/4/Auto for each scene (e.g. pool −3.17e-4, dam −3.28e-4, droplet 3.4e-7). No extra leaks/holes from the parallel path.

## Next if a larger grid appears

Re-run `--thread-benchmark` at 400×240+. If pressure lists exceed ~20k cells, try 2 workers first.
