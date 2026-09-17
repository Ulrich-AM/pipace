# PIPACE — Agent Context

Read this before changing the codebase. Prefer this file plus the roadmap/notes over rediscovering the project from scratch.

## What PIPACE is

PIPACE is a **pixel-native physics/chemistry sandbox** (Win32 + C++17, no external deps). The current prototype is a **free-surface liquid engine** (water), not the full game yet.

Long-term architecture (cooperating systems, not one universal solver):

| Responsibility | Answers |
|---|---|
| **Physics engines** (liquid / solid / gas) | How does matter with these properties move? |
| **SACE** (Somewhat Accurate Calculator for Elements) | What is this substance, what can it become, what properties should new substances have? |
| **Chemistry** (`ReactionDefinition`) | What reactions are possible, with what stoichiometry, conditions, and ΔH? Data only today. |
| **World state** | Where is it, how much, temperature/pressure/phase/local conditions? |
| **Rendering** | How should the current state be drawn? |

User-facing UI may call any spawnable material an **"Element"**. Internally prefer **Substance / Material / MaterialDefinition**. Player names are metadata; chemical identity is simulation data.

## Current code layout

```
PIPACE/
  main.cpp                 # Win32 loop, world pixel render, input to engines
  ui/                      # GDI chrome: layout, settings overlay, inspector/properties, UiLanguage
  assets/                  # category tray icons (PNG)
  languages/english.json   # editable UI strings; loadLanguage() in main.cpp
  gas/
    GasEngine.h/.cpp       # Air amount/pressure/flow (temperature-aware, P still isothermal)
  thermal/
    ThermalTypes.h / ThermalConfig.h / ThermalEngine.h/.cpp  # heat storage, conduction, sleep
  fluid/
    FluidTypes.h           # Grid constants, enums, SplashParticle, LiquidProperties alias, TimingAverages
    FluidConfig.h          # FluidConfig (solver knobs; no liquid property copies)
    FluidEngine.h/.cpp     # Owns all fluid arrays + simulation stages (behavior-preserving extract)
  rigid/
    RigidBodyTypes.h       # MaterialId masks + MaterialDefinition adapter (from SubstanceDefinition)
    RigidBodyEngine.h/.cpp # Drawable pixel-native rigid bodies + world/fluid coupling
  substance/
    SubstanceProperties.h  # Grouped intrinsic properties (mechanical/fluid/thermal/phase/porous/…)
    SubstanceTypes.h / SubstanceRegistry.h/.cpp  # SubstanceId, MatterPhase, registry
    PhaseTransfer.h/.cpp  # fill/mass/gas-amount conversion
  world/
    WorldQuery.h/.cpp      # sampleMatterAt (SubstanceId + phase), phase/registry diags
    WaterPhaseChange.h/.cpp # water liquid ⇄ gas / solid; mixture skip via composition API
  chemistry/
    ReactionTypes.h / ReactionRegistry.h/.cpp  # ReactionId + ReactionDefinition data only (no engine)
  docs/
    PHASE_CHANGES.md       # Water solid ⇄ liquid ⇄ gas is live; honey mixtures are not
  CMakeLists.txt / build.bat / run.bat
  README.md
  AGENTS.md
  misc/                   # notes, changelogs, diagnostic TSV dumps
```

Physical properties are canonical on `SubstanceDefinition` grouped structs
(`mechanical` / `fluid` / `thermal` / `phase` / `porous`). Do not add new
duplicate property tables.

**Identity vs phase vs engine:**

- `SubstanceDefinition` = intrinsic identity + properties (WHAT it is)
- `MatterPhase` = current represented phase (WHICH phase this cell/body is)
- engine storage = numerical representation (HOW that phase is simulated today)

These mappings are **implementation, not laws**:

| Current representation | Engine |
|---|---|
| Water + Liquid, Honey + Liquid | FluidEngine (water/honey volume channels) |
| Wood/Glass/Metal + Solid | RigidBodyEngine (MaterialId masks) |
| Water + Solid | RigidBodyEngine (`MATERIAL_WATER_SOLID` mask, still SUBSTANCE_WATER) |
| Stone + Solid | rigid body **or** static `solid[]` walls (`kStaticWallSubstance`) |
| Air + Gas | GasEngine (air component) |
| Water + Gas | GasEngine (`waterVapor` component; same SubstanceId) |

Static `solid[]` walls use Stone thermal/mechanical identity but are
**moisture-inert** (no absorb / drip / wetness). Porous Stone rigid bodies still
absorb. Condensation may form liquid in adjacent free cells, not inside the wall.

`supportsPhase` is capability metadata. Water **liquid ⇄ gas** and **liquid ⇄ solid**
are live. Do **not** add SUBSTANCE_ICE / SUBSTANCE_STEAM.
Density stays phase-specific: mechanical.densityRel (solid), fluid.density (liquid),
gas amount/EoS (gas). Honey/water mixtures do not boil or freeze yet.

Water/honey cells that hold both channels are **mixtures**, not a new SubstanceId.
High-level composition queries use `FluidEngine` SubstanceId APIs
(`liquidComponentAmount` / `liquidComponentFraction` / `dominantLiquidSubstance` /
`liquidComponents`). Storage is a fixed-capacity SoA of up to 4 components per
cell (`liquidCompId` / `liquidCompAmt` / `liquidCompCount`); `fill[]` remains
total occupancy. `MatterSample` reports the dominant component plus fractions.

Solver unit liquid is `sandboxReferenceLiquid()` (currently SUBSTANCE_WATER's fluid
table: relative density 1.0). That is a reference, not “all liquid is water”.

Effective liquid density / Cp / k / surface tension come from
`evaluateLiquidMixture`; viscosity from `evaluateLiquidMixtureViscosity(T)`.
Water phase change requires `liquidCompositionIsPureWater`. Empty or invalid
composition uses reference numbers only (no silent Water identity).

**Compatibility adapters that remain (justified):**

- `MaterialDefinition` / `materialDef(MaterialId)` — rigid masks still store
  `MaterialId`; fracture, moisture, mass, and RGB read this view. Values come
  from `SubstanceDefinition` at first use.
- `LiquidProperties` — name alias for `FluidProperties`.
- `GasSpecies::Air` — single bridge to `SUBSTANCE_AIR` + `MatterPhase::Gas`.
- `LiquidPaint.asHoney` — UI paint flag; simulation uses `paint.substance()`.

Removed copies: `FluidConfig.water` / `.honey`, `kHoneyLiquid()`, `k*Thermal()`
wrappers, `GasConfig.thermal`.

Build: `run.bat` or CMake → `build/pipace.exe`. Headless: `--benchmark`, `--scale-benchmark`, `--rigid-benchmark`, `--thread-benchmark`, `--liquid-diag`, `--liquid-composition-diag`, `--substance-phase-diag`, `--substance-registry-diag`, `--phase-transfer-diag`, `--water-phase-diag`, `--thermal-diag`, `--thermal-spread-diag`. Grid size via `PIPACE_GRID_WIDTH` / `PIPACE_GRID_HEIGHT` (default 200×120). SETTINGS → Simulation threads (Auto / 1 / 2 / 4 / 6 / 8). Auto is 1 worker on the default grid; see `misc/THREAD_PASS_NOTES.md`.

## Fluid engine (what exists)

Hybrid free-surface solver:

- Fractional fill per cell; MAC staggered `u`/`v`
- Semi-Lagrangian velocity advection; optional extrema-clamped **BFECC** (`A`)
- Gravity at 4 cells/m × 9.81 m/s² → 39.24 cells/s²
- Low-viscosity water skips spatial diffusion; path exists for thicker liquids
- Red-black Gauss-Seidel / SOR pressure projection with fractional free-surface weights, warm-start, adaptive 8/14/24 iters + residual early-out
- Conservative face-flux transport with iterative donor/receiver limiter (A→B→C through-flow)
- Residual consolidation; open edges are a void sink
- Surface field (smoothed fill, normals, curvature), surface tension, ballistic splash particles
- Optional vorticity confinement (`O`)
- 16×16 chunk wake/sleep; active solve region + halo
- Volume / momentum / KE diagnostics; F1–F10 test scenes; stage timings
- `FluidConfig.walledBorders` also gates thermal exchange with an infinite ambient
  reservoir at `AMBIENT_TEMPERATURE_K` (no separate thermal-border toggle). Gas
  mass-carried heat stays on `GasEngine::escapedHeat`.

`LiquidProperties` is a compatibility name for `FluidProperties`. Fluid lookup is
`fluidForSubstance` / `mix*` / `sandboxReferenceLiquid()`. Two volume channels
(water/honey) exist; there is no mixture SubstanceId.

## Modularization rules (critical)

This project is mid **behavior-preserving modularization** (roadmap Phase 6).

**Done so far:** shared types/config + `FluidEngine` shell (all fluid arrays and solver stages live as engine methods; `main.cpp` is Win32 UI/render/input).

1. **Do not change physics while extracting** unless the task explicitly asks for a physics fix.
2. Extract one subsystem at a time; compile and run benchmarks after each step.
3. `FluidEngine` owns arrays and scratch buffers — avoid reintroducing file-scope simulation globals.
4. Next peel order: pressure → volume transport → advection → surface/splash → chunks → render → UI → tests (into separate `.cpp` files under `fluid/` / `render/` / `ui/` / `tests/`).
5. Target layout (not all present yet): `fluid/`, `world/`, `rigid/`, `gas/`, `substance/`, `chemistry/`, `render/`, `ui/`, `tests/`.

## Known issues (do not “fix” by redesigning the solver)

1. **Thin-wall / corner leaks** — `depositVolume()` radius search can still bypass solid topology. Residual consolidation now blocks diagonal corner cuts and no longer walks into empty air (trail fix). Full connectivity-aware redeposit remains a follow-up.
2. **Paint/erase stutter** — largely addressed via deferred `paintDirty` flush; further dirty-region rebuilds possible.
3. **Renderer GDI churn** — addressed with persistent backbuffer/font.
4. **Pressure / transport / surface cost** — addressed with sparse lists + stencil cache; further multithreading optional.
5. **Transient one-pixel “holes”** — often render-threshold / residual flicker; do not invent free liquid to fill them.
6. **Falling-stream Fill trails** — were Eulerian sub-visible crumbs extended by residual consolidation preferring empty cells below; consolidation now merges only into cells that already hold liquid (topology-safe).

## Optimization headroom

Performance pass (2026-09-12) landed sparse flux/pressure/surface lists, pressure stencil cache, deferred paint, persistent GDI, and a persistent WorkerPool. Threading pass (2026-09-13) exposes Auto/1/2/4/6/8 in Settings. On 200×120, Auto stays serial: red/black pressure barriers lost to overhead (see `misc/THREAD_PASS_NOTES.md`). Parallel pressure is gated at ≥20k cells/phase.

Remaining headroom:

- Quality presets (Low / Medium / High / Auto) change cost knobs only: pressure iters, CFL substeps, limiter passes, surface tension, spray, gas rate, sim Hz, flat render. Auto drops Low after several overloaded ticks.
- Connectivity-aware residual / splash redeposit (correctness)
- Broader multithreading on **larger grids** (pressure lists ≥20k cells/phase); do not force 200×120 onto many cores
- Avoid guessing; use `--benchmark` and stage timings

Do **not** sacrifice conservation or replace the donor/receiver limiter with scan-ordered movement for speed.

## Design principles for agents

1. Conservation first — never silently create/delete matter to hide bugs.
2. Local/chunk work over full-grid work.
3. Shared engines + material-driven properties; no custom solver per material.
4. Chemistry (SACE) separate from phase motion.
5. Rendering must not change physics.
6. Quality settings change cost/fidelity, not fundamental world laws.
7. Benchmark before optimizing; stop when acceptance criteria are met.
8. Refactor separately from physics changes.

## Near-term roadmap order (condensed)

Completed architecture: SubstanceId, registry, grouped properties, MatterPhase /
MatterIdentity, world query, material migration, moisture stabilization.

**Next major milestone: SACE / mixtures** — water `liquid ⇄ gas` and `liquid ⇄ solid`
are live (`world/WaterPhaseChange.cpp`). Honey mixtures still skip phase change.

Longer sequence (historical): liquid correctness → performance → modularization →
honey composition → rigid coupling → temperature/gas → **phase changes** → SACE.

Detail lives in `misc/PIPACE_Project_Roadmap_TODO.txt` and `misc/PIPACE_SACE_Concept_Notes.txt`.

## Personality / presentation

Homely, utilitarian, slightly garage-lab — not neon sci-fi. Flavor text must never alter simulation behavior.
