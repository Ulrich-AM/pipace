# PIPACE — Agent Context

Read this before changing the codebase. Prefer this file plus the roadmap/notes over rediscovering the project from scratch.

## What PIPACE is

PIPACE is a **pixel-native physics/chemistry sandbox** (Win32 + C++17, no external deps). The current prototype is a **free-surface liquid engine** (water), not the full game yet.

Long-term architecture (cooperating systems, not one universal solver):

| Responsibility | Answers |
|---|---|
| **Physics engines** (liquid / solid / gas) | How does matter with these properties move? |
| **SACE** (Somewhat Accurate Calculator for Elements) | What is this substance, what can it become, what properties should new substances have? |
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
    FluidTypes.h           # Grid constants, enums, SplashParticle, LiquidProperties, TimingAverages
    FluidConfig.h          # FluidConfig (density/viscosity/surface tension, advection, CFL, etc.)
    FluidEngine.h/.cpp     # Owns all fluid arrays + simulation stages (behavior-preserving extract)
  rigid/
    RigidBodyTypes.h       # MaterialId, pixel-mask RigidBody, merged runs, contacts
    RigidBodyEngine.h/.cpp # Drawable pixel-native rigid bodies + world/fluid coupling
  substance/
    SubstanceTypes.h / SubstanceRegistry.h/.cpp  # SubstanceId identity + built-in registry
  CMakeLists.txt / build.bat / run.bat
  README.md
  AGENTS.md
  misc/                   # notes, changelogs, diagnostic TSV dumps
```

Build: `run.bat` or CMake → `build/pipace.exe`. Headless: `--benchmark`, `--scale-benchmark`, `--rigid-benchmark`, `--thread-benchmark`, `--liquid-diag`. Grid size via `PIPACE_GRID_WIDTH` / `PIPACE_GRID_HEIGHT` (default 200×120). SETTINGS → Simulation threads (Auto / 1 / 2 / 4 / 6 / 8). Auto is 1 worker on the default grid; see `misc/THREAD_PASS_NOTES.md`.

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

`LiquidProperties` already exists for future multi-liquid / SACE-driven materials. Behavior is still water-centric.

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

1. Liquid correctness (leaks, void boundary, hole flicker)
2. Performance pass (lists, caches, paint batching, renderer)
3. Finish modularization
4. Dye tracers → honey multi-liquid proof
5. Pixel-native rigid bodies + two-way fluid coupling (first drawable-body milestone is in `rigid/`)
6. Temperature → gas engine → SACE → phase changes

Detail lives in `misc/PIPACE_Project_Roadmap_TODO.txt` and `misc/PIPACE_SACE_Concept_Notes.txt`.

## Personality / presentation

Homely, utilitarian, slightly garage-lab — not neon sci-fi. Flavor text must never alter simulation behavior.
