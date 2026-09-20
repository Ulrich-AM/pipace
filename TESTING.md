# PIPACE testing and headless diagnostics

Run commands from the repository root on Windows. The standard release-style
build is:

```text
build.bat
```

The executable is written to `build/pipace.exe`. CMake is also supported, but
baseline numbers should only be compared when compiler, optimization level,
grid size, and machine load are comparable.

The default grid is 200x120. Set `PIPACE_GRID_WIDTH` and
`PIPACE_GRID_HEIGHT` before building to test another compile-time grid size.
Files whose names contain `<grid>` use the resulting dimensions, for example
`benchmark_200x120.tsv`.

## Command catalog

| Command | Output in `misc/` | Invariant or measurement | Expected pass condition |
|---|---|---|---|
| `build\pipace.exe --liquid-composition-diag` | `liquid_composition_diag.tsv` | Generic fixed-capacity liquid composition storage: SubstanceId slots, sum=fill, overflow reject, take/deposit/splash/advection conservation, Water/Honey gameplay, and composition invariants (no silent Water repair, liquid-capability IDs only, empty dominant is NONE). | Final row is `summary PASS`. Component amounts are nonnegative, sum to fill, overflow does not destroy matter, transferred composition stays with fill, and invalid fill-without-composition is detected without inventing Water. |
| `build\pipace.exe --rigid-contact-diag` | `rigid_contact_diag_<grid>.tsv` | Resting stability, slopes, no-gravity stability, overlap recovery, and off-center torque. | Final row is `summary PASS` and every test row is `PASS`. |
| `build\pipace.exe --solid-diag` | `solid_diag.tsv` | Grab behavior, damage/fatigue/material differences, moisture absorption/spread, fracture, pixel and mass conservation, crack persistence, anchoring, and resting stability. | Final row is `summary PASS`; fracture conserves pixels/mass and absorption conserves free + absorbed liquid within the encoded tolerance. |
| `build\pipace.exe --moisture-diag` | `moisture_diag.tsv` | Wood/stone absorption, nonporous materials, rotated bodies, wet mass, moisture spread, and fracture conservation. | Final row is `summary PASS`. Total free + absorbed + splash + pending liquid remains within 0.12 units; wood absorbs and gains matching mass; glass/metal absorb approximately zero. |
| `build\pipace.exe --moisture-drip-diag` | `moisture_drip_diag.tsv` | Drip onset/decay, pending-volume delivery, saturation threshold, outlet limits, orientation, material differences, fracture, and reabsorption oscillation. | Final row is `summary PASS`. Conserved cases stay within 0.12 liquid units, pending drips are not stranded, damp wood does not drip, and glass/metal emit zero. |
| `build\pipace.exe --thermal-diag` | `thermal_diag.tsv`, `thermal_stats.txt` | Heat mixing/conduction, heat following moving liquid and rigid pixels, fracture inheritance, thermal sleep, closed energy conservation, conduction clamping, and ambient stability. | Every `pass` field is `1`; closed relative energy error is below 0.02; no NaN, infinity, negative-K, or orphan heat storage; ambient seed/idle spread stays within the encoded limits. |
| `build\pipace.exe --thermal-spread-diag` | `thermal_spread_diag.tsv` | Localized hot/cold gas, hot wall, hot liquid, buoyancy rise, thermal sleep, sealed energy, open-boundary heat accounting, and open/walled edge conduction vs an interior control. | Final row is `summary PASS`. A sealed hot gas patch cools at the seed, warms gas above it, and rises. Open `walledBorders=false` edges exchange energy with ambient at `AMBIENT_TEMPERATURE_K`; walled edges do not. External heat (conduction + gas/liquid mass carry) closes `E0 + entered − escaped ≈ E1`. |
| `build\pipace.exe --gas-diag` | `gas_diag.tsv` | Uniform/sealed gas, opening, compression, expansion, rigid displacement, and long sealed conservation. | Every `pass` field is `1`. Conserved cases meet their encoded absolute amount tolerances (0.001 to 0.05 depending on case), and pressure/flow checks pass. |
| `build\pipace.exe --substance-phase-diag` | `substance_phase_diag.tsv` | Phase capability metadata, identity adapters, world queries, mixture reporting, static-wall identity, gas identity, and empty/out-of-bounds queries. | Final row is `summary PASS`; water remains one ID across supported phases and no Ice/Steam IDs exist. |
| `build\pipace.exe --substance-registry-diag` | `substance_registry_diag.tsv` | Registry IDs/names/indexing, invalid lookup fallback, canonical grouped properties, compatibility adapters, selected reference property values, and SACE exact/unknown chemical identities. | Final row is `summary PASS`; all built-in definitions and canonical property checks pass. `SUBSTANCE_COUNT` remains 12. |
| `build\pipace.exe --sace-identity-diag` | `sace_identity_diag.tsv` | Exact elemental composition, canonical signatures, built-in identity lookup, honest unknown materials, and reaction atom-balance (`Balanced` / `Unbalanced` / `Unknown`). | Final row is `summary PASS`. H2+O2 and C+O2 are Balanced; an intentionally unbalanced reaction is Unbalanced; unknown composition is Unknown. Formula text is not unique identity. |
| `build\pipace.exe --phase-transfer-diag` | `phase_transfer_diag.tsv` | Liquid fill ↔ mass ↔ water-vapor amount round trips, latent-heat helpers vs registry, Clausius–Clapeyron `P_sat`/`T_sat` helpers, `canTransition` metadata, and 1000-cycle closed conversion drift. No live boiling/freezing. | Final row is `summary PASS`. Water remains one ID (no Ice/Steam). Round-trip mass error and 1000-cycle drift stay within the encoded tolerances. Invalid transitions fail safely. Saturation at the reference boiling point matches `referencePressurePa`. |
| `build\pipace.exe --water-phase-diag` | `water_phase_diag.tsv` | Live water boiling and condensation: latent plateau, mass transfer, sealed pressure, honey skip, vapor identity, and open-boundary vapor accounting. Temperature plateau is compared to local `T_sat(P)`, not a fixed 373 K. | Final row is `summary PASS`. No Ice/Steam IDs. Honey mixtures do not boil. Closed liquid+vapor mass stays within the encoded tolerance. |
| `build\pipace.exe --water-phase-validation` | `water_phase_validation.tsv` | Tick-rate (20/30/60 Hz), residual superheated water, closed sensible+latent energy, pressure-aware equilibrium, condensation/saturation, hot-box melt, static-wall dryness, and storage integrity. | Final row of each case is `PASS` (no FAIL rows). 20/30/60 Hz phase-mass difference stays under 2%. Closed unforced energy relative error stays under `1e-4`. |
| `build\pipace.exe --water-solid-phase-diag` | `water_solid_phase_diag.tsv` | Live water freezing/melting: ice identity, mass/latent accounting, connected rigid ice, melt splits, tick-rate, and a full solid→liquid→gas→liquid→solid cycle. | Final row is `summary PASS`. No Ice SubstanceId. Honey mixtures do not freeze. Solid water is `SUBSTANCE_WATER` + `MatterPhase::Solid` via `MATERIAL_WATER_SOLID`. |
| `build\pipace.exe --water-phase-stability-diag` | `water_phase_stability_diag.tsv` | Hot-box melt/boil, hot vs cool wall condensation, left/right symmetry, and static-wall moisture exclusion. | Final row is `summary PASS`. Hot walls do not attract condensate. Static `solid[]` walls stay dry; Stone rigid bodies still absorb. |
| `build\pipace.exe --benchmark` | `benchmark_<grid>.tsv` | Main fluid scene timings, stage timings, work counts, liquid volume error, thin cells, and splash counts. | Observational benchmark: all values are finite, volume error remains negligible, and timings/work counts show no unexplained regression against a comparable baseline. |
| `build\pipace.exe --advection-benchmark` | `advection_benchmark_<grid>.tsv` | Six velocity-advection modes across four scenes; timing, volume, velocity, momentum, and kinetic energy. | Observational benchmark: 24 case rows, finite metrics, negligible volume error, and no unexplained mode-specific instability or timing regression. |
| `build\pipace.exe --scale-benchmark` | `scale_benchmark_<grid>.tsv` | Four representative fluid scenes at the compiled grid size. | Observational benchmark: four case rows, finite timings, and negligible volume error. Rebuild at each desired grid size for scaling comparisons. |
| `build\pipace.exe --thread-benchmark` | `thread_benchmark_<grid>.tsv` | Worker modes 1/2/4/Auto across five fluid scenes plus a rigid-pool scene; stage timing and conservation. | Observational benchmark: 24 case rows, finite metrics, negligible volume error, and worker/parallel flags consistent with thresholds. On the default grid Auto is expected to resolve to one worker. |
| `build\pipace.exe --rigid-benchmark` | `rigid_benchmark_<grid>.tsv` | Coupled rigid/fluid timing, liquid conservation/loss, and an F1 settling snapshot. | Observational benchmark: finite timing, negligible liquid error, `lost_rigid == 0`, and settling metrics should not regress from a comparable known-good baseline. |
| `build\pipace.exe --look-bench` | `look_bench.tsv` | Paint cost for each world-rendering look, including the combined effects case. | Observational benchmark: seven finite, nonnegative `ms_per_paint` rows and no unexplained regression on the same machine/configuration. |

Benchmarks do not currently emit PASS/FAIL summaries. Treat conservation and
finite-value checks as correctness gates; treat timings as comparative data,
not fixed universal thresholds.

## Core regression invariants

- Matter is not created or silently discarded to hide topology or rendering
  problems. Closed liquid totals include grid fill, splashes, absorbed
  moisture, and pending drip volume where applicable.
- An open world edge is an explicit sink. Escaped matter must be accounted for
  by the matching expected-total adjustment.
- Rigid fracture conserves occupied pixels, material mass, moisture, pending
  drip volume, and per-pixel heat. Rotation is always rasterized from the
  authoritative local mask.
- Gas amount is conserved in sealed scenes; pressure changes must follow
  compression, expansion, and newly opened connections without inventing gas.
- Closed thermal exchange conserves energy and never produces NaN, infinity,
  negative absolute temperatures, or heat stranded in empty storage.
- `SubstanceId` identifies the substance, `MatterPhase` identifies its current
  represented phase, and engine storage is only the numerical representation.
  Water does not gain separate Ice or Steam IDs.
- Water/honey mixtures remain explicit component fractions and report a
  dominant existing substance; they do not create a fake mixture ID.
- `SubstanceDefinition` and its grouped property tables are canonical.
  Compatibility adapters must agree with the registry.
- Water liquid ⇄ gas is live. Water liquid ⇄ solid (rigid ice pixels) is live.
  Honey/water mixtures do not boil or freeze yet. There is no Ice/Steam SubstanceId.

## Manual water boiling / condensation

Use the default map, Heat/Cool from Energy, and the inspector (hover a cell).

1. **Boiling.** Paint a small pure-water pool. Heat it until the inspector
   temperature sits near the local saturation temperature (373 K at 1 atm;
   lower in vacuum, higher under pressure). Liquid fill should fall gradually;
   gas cells should show Water (or Air/Water composition) rather than a Steam
   id. Temperature should stall near `T_sat(P)` while vapor is produced.
2. **Condensation.** Cool the vapor or nearby walls. Liquid water should
   reappear nearby when the water-vapor partial pressure exceeds saturation;
   undersaturated warm vapor should not rain just because `T < 373 K`.
3. **Sealed chamber.** Draw a closed wall box, put water inside, heat, then
   cool. Pressure should raise the boiling temperature. Ice in a hot box must
   melt from incoming heat, then the liquid follows equilibrium. Static walls
   stay dry.
4. **Mixture.** Paint water+honey in one cell and heat it. It should not boil
   in this pass.

Automated coverage: `build\pipace.exe --water-phase-diag`.


## Pre-phase-change baseline: 2026-09-15

Environment: default 200x120 grid, 8 logical hardware workers reported by the
benchmark, release-style `g++ -O2` build from `build.bat`.

Diagnostic summaries:

| Suite | Result | Conservation note |
|---|---:|---|
| Liquid bug/topology | all topology rows `ok` | Largest ordinary sampled closed-scene error: `0.000324034` units (absolute). |
| Rigid contact | PASS, 7/7 | Stability/response suite; no volume total. |
| Solid | PASS, 18/18 | Wood absorption total delta `-0.000129` units; fragment pixel/mass checks passed. |
| Moisture | PASS, 6/6 | Wood soak total delta `-0.000123` units. |
| Moisture drip | PASS, 9/9 | Reported conserved-case deltas were `0` to `+0.000001` units. |
| Thermal | 11/11 case flags set | Closed relative energy error `3.05157e-07`; zero NaN/Inf/negative-K/orphan-storage counts. |
| Gas | 8/8 case flags set | Largest absolute amount error `1.25019e-05` units. |
| Substance phase | PASS, 37/37 | Identity/phase/world-query invariants passed. |
| Substance registry | PASS, 49/49 | Registry/property/adapter invariants passed. |

Main 200x120 fluid benchmark (`physics_ms` per tick):

| Scene | ms/tick | Volume error |
|---|---:|---:|
| still_pool | 20.3445 | -0.000352543 |
| small_puddle_baseline | 3.11918 | +0.00000277035 |
| small_puddle_refined | 2.34810 | -0.0000220661 |
| dam_break | 16.7845 | -0.000314433 |
| large_moving_body | 15.4407 | -0.000224652 |
| single_droplet | 2.75260 | +0.000000664033 |
| narrow_nozzle | 14.0153 | -0.000716001 |

Additional benchmark records:

- Advection: 24 rows; average wall-clock `per_tick_ms` by mode was None
  21.06, FOU 21.62, NSL 14.24, SL 14.07, MacCormack 19.51, and BFECC
  18.52. Maximum absolute volume error was `0.000561863`.
- Scale benchmark at 200x120: resting pool 20.3105, small puddle 2.23658,
  dam break 18.2565, and large moving body 9.42593 ms/tick. Maximum absolute
  volume error was `0.000200231`.
- Thread benchmark: Auto resolved to one worker; all `parallel_p` values were
  zero at this grid size. Exact per-scene/stage values are in
  `misc/thread_benchmark_200x120.tsv`; maximum absolute volume error was
  `0.000427944`.
- Coupled rigid benchmark: 32.9252 ms/tick overall, 0.112 ms last rigid step,
  liquid error `-0.000455481`, and `lost_rigid = 0`.
- Moisture soak diagnostic: 10.958739 ms/tick.
- Renderer: FlatColor 0.809989, AlphaFlat 0.977010, NoisyFlat 1.21161,
  Detailed 0.994469, Realistic 1.64019, Legacy 2.59822, and
  Realistic+Glow+Outlines 4.07705 ms/paint.

### Baseline warning

The F1 settling snapshot in `rigid_benchmark_200x120.tsv` reports four bodies,
only one sleeping body, and `maxSpeed = 15.0498` after 150 ticks. This benchmark
has no PASS/FAIL gate, but the committed historical output reported one body,
one sleeping body, and zero maximum speed. Treat this as suspicious and
investigate/reproduce it before relying on F1 settling as a phase-change
regression gate. No physics change was made during baseline collection.
