# Free-Surface Refinement Pass

## Rigid-body contact stability (2026-09-12)

Drawn bodies were hopping because static contacts used a fake 0.65-cell penetration on every overlapped floor pixel and then stacked full COM translations. Contacts now sit on occupied-cell / solid faces with measured depth, positional correction is merged per normal (max depth, not sum), floor velocity constraints are one manifold per supporting plane, sleeping bodies are not shoved, same-cell body overlap is a real contact, gravity is not re-applied while supported, and fluid pressure force is clamped / low-pass filtered so pools do not rocket wood blocks. Collision uses center-sampled occupancy; fluid exclusion stays conservative 5-tap. Shift+F1 now sleeps after landing (`--rigid-benchmark` settle_f1).

## Observed causes

**Rapid spreading:** cells down to 0.2% fill were rendered and participated in pressure/transport like ordinary liquid. Conservative flux could therefore distribute a small puddle over hundreds of barely occupied cells. The older receiver limiter also blocked simultaneous `A -> B -> C` through-flow, leaving incoherent residuals after energetic motion.

**Green-looking water:** fill, depth, speed, and surface terms independently modified RGB channels. At low fill the blue separation was no longer guaranteed, producing dark teal/green-looking pixels.

**Erase lag:** geometry edits woke simulation work, but the numerical stages still swept world-sized arrays. Repeated 64-pass flux limiting and full-domain pressure work amplified the stall.

## Physics fixes

- Preserved the staggered MAC/free-surface architecture.
- Split render, active, and pressure fill thresholds.
- Added conservative residual consolidation; it moves tiny volume into physically plausible lower/substantial/velocity-aligned neighbors and never deletes it.
- Retained previous-frame pressure as the warm start and removed vertical hydrostatic reconstruction.
- Added fractional face occupancy to pressure stencils, pressure gradients, and liquid transport.
- Added deterministic iterative donor/receiver limiting with simultaneous through-flow, a 16-pass safety cap, convergence early-out, and conservative clamp-error redistribution.
- Added a persistent 3x3-smoothed liquid color field, surface mask, normals, curvature, and curvature-normal surface force.
- Replaced random spray eligibility with a breakup-energy score using normal velocity, speed, pressure gradient/impulse, and curvature.
- Removed global water velocity damping. Low-viscosity water skips diffusion; a spatial Laplacian diffusion path is available for future thick liquids.
- Added optional extrema-clamped BFECC velocity advection (`A`); semi-Lagrangian remains the default.
- Defined scale as 4 cells/meter and gravity as 9.81 m/s², yielding 39.24 grid cells/s².
- Playback speed now scales fixed 1/30-second steps instead of changing the physical timestep.

## Rendering and diagnostics

- Normal water uses a blue-constrained palette even at tiny fractional fill.
- Fill debug view maps zero to black and increasing fill through dark blue, blue, cyan, and white.
- Hover data reports fill, pressure, cell velocity, surface flag, and chunk state.
- Added thin-cell count/volume, momentum, kinetic energy, CFL-cap state, adaptive pressure iteration count, and averaged stage timers.
- Added F6 droplet, F7 puddle, F8 incline, F9 nozzle, and F10 thin-channel scenes.

## Performance fixes

- 16x16 chunks sleep from velocity, divergence, fill-change, and pressure-change metrics.
- Painting, erasing, solid changes, fluid transfer, and splashes wake chunks plus a halo.
- The pressure solve covers the connected disturbed liquid region plus halo rather than isolated chunks.
- Advection, gravity, surface work, projection, transport, clamping, and boundary enforcement use active-region bounds.
- Active buffers are cleared/copied only over their working ranges.
- The activity rebuild scans active, fluid-bearing, and solve-halo chunks; this also fixed liquid entering a halo without waking it.
- Pressure uses 8/14/24 iteration budgets plus residual early-out.

## Deterministic before/after checks

All checks use the release build and fixed 1/30-second steps. Run `build\pipace.exe --benchmark` to regenerate the TSV.

| Check | Before | After |
|---|---:|---:|
| F7 thin cells after 240 ticks | 300 | 1 |
| F7 thin volume | 9.73766 | 0.05597 |
| F7 total volume | 36.000002 | 35.999996 |
| F7 volume error | +0.00000165 | -0.00000352 |
| F3 dam-break volume error | -10.3054 (pre clamp-accounting fix) | -0.000386 |
| F1 pool volume error after 180 ticks | n/a | -0.001004 of 6360 |

## Active-region scaling check

The same scene footprint was run on increasingly large empty canvases. This isolates overhead caused by total world size from work caused by the active liquid region.

| Grid | Small puddle ms/tick | Resting pool ms/tick | Dam break ms/tick | Moving body ms/tick |
|---|---:|---:|---:|---:|
| 200x120 | 3.09 | 26.42 | 28.88 | 34.30 |
| 400x240 | 6.52 | 44.70 | 68.23 | 52.25 |
| 800x480 | 7.50 | 50.01 | 57.33 | 40.44 |

The small active puddle is 68-71 active cells in every size. A 16x increase in total cells raises its cost about 2.4x, rather than 16x. Larger connected bodies still dominate pressure and transport work, as expected. Individual timings vary with system load; the stage averages in the toolbar are the best live guide.

## Remaining calibration target

The small-puddle objective is met without deleting volume. The broad F1 pool still carries noticeable low-amplitude numerical motion after six simulated seconds (about 3.35 cells/s maximum in the deterministic run) and keeps a large connected pressure region awake. Improving that case should use a better pressure preconditioner/solver or tighter free-surface reconstruction—not a return to manual vertical pressure, arbitrary friction, or CA movement.
