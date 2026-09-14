# Fluid Performance Optimization Pass — Results

Date: 2026-09-12  
Baseline: `build/OPT_BASELINE_before.tsv`  
After: `build/OPT_AFTER.tsv`  
Benchmark mode forces `workerCount = 1` for deterministic comparison.

## What changed

1. **Finer profiling / work counts** — stage timings (surface, viscosity, residual, drain, listsBuild) and counters (pressure cells, surface cells, flux faces, limiter passes, subvisible, splashes).
2. **Nonzero flux face lists** — limiter/transfer only touch faces with flux; dense fallback when >40% of region faces are active.
3. **Pressure red/black lists + stencil cache** — rebuild per projection; sparse when pressure cells ≤70% of solve rectangle; SOR walks lists/stencils.
4. **SurfaceCells list** — normals/curvature/tension use the surface list (face application unchanged).
5. **Deferred paint finalize** — `paintDirty`; flush on mouse-up and at tick start (no per-mousemove boundary rebuild).
6. **Persistent GDI backbuffer + font** — recreate only on resize.
7. **WorkerPool** — optional parallel red/black pressure when ≥2000 cells and `workerCount != 1`. Benches stay single-threaded.
8. **Residue diagnostics** — hover shows prev fill + splash count; `subvisibleCells` / `residualTransfers` counters. **No mass deletion.**

## Before / after (physics ms / tick)

| Scene | Before | After | Speedup | Volume error after |
|---|---:|---:|---:|---:|
| still_pool | 16.39 | 7.69 | **2.1×** | -0.00104 |
| small_puddle_baseline | 4.06 | 0.80 | **5.1×** | ~2e-6 |
| small_puddle_refined | 3.32 | 0.68 | **4.9×** | ~-3e-6 |
| dam_break | 20.10 | 6.85 | **2.9×** | -0.00012 |
| large_moving_body | 22.00 | 6.46 | **3.4×** | -0.00015 |
| single_droplet | 1.96 | 1.07 | **1.8×** | ~1e-6 |
| narrow_nozzle | 7.19 | 5.61 | **1.3×** | -0.00011 |

Pressure and transport stage times dropped sharply on dense scenes (e.g. dam_break pressure 7.88→1.38, transport 9.71→3.95).

## Behavior

- Total volumes match baseline scenes within float noise.
- Conservation errors remain ~1e-3 or better.
- Small-puddle thin-cell behavior with merge on/off preserved in character (refined → 1 thin cell).
- `large_moving_body` shows small differences in thin-cell count / vmax (FP / activity path); volume conserved.

## Risks / invalidation

- Pressure stencils rebuilt every projection (safe; free-surface coeffs still use live fill).
- Sparse flux path must stay equivalent to dense limiter; dense fallback retained for flux-heavy regions.
- Multithreading off by default in benches; enable via `FluidConfig::workerCount` (0=Auto). Expect ulp-level pressure differences when parallel.
- Paint deferral: solids/water edits apply fully on mouse-up or next physics tick.

## Defaults

All sparse lists, stencil cache, deferred paint, and persistent GDI remain **enabled**. Worker count defaults to **1** (serial); set `0` for Auto in interactive builds if desired later.

## Residue / trails note

**Cause (confirmed):** Falling-stream “particles” in Fill view are mostly **Eulerian fractional-fill leftovers** from conservative transport (`0 < fill < MIN_RENDER_FILL`). Normal view hides them. They are **not** SplashParticles (those remain legitimate breakup).

**Amplifier:** `consolidateResidualVolume` used to score empty cells below at +120, so crumbs walked downward into empty air each tick and painted vertical trails in Fill view.

**Fix:** Merge only into neighbors that already hold liquid; require an orthogonal gate for diagonals; never `fill = 0`. Splash spawn thresholds unchanged.
