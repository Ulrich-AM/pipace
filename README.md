# Water Test

A dependency-free C++/Win32 pixel-water playground based on the supplied hybrid fluid-engine specification.

## Run

Double-click `run.bat`. It builds `build\water_test.exe` on the first run with the MSYS2 UCRT64 `g++` compiler, then launches it.

You can also build with CMake:

```powershell
cmake -S . -B build
cmake --build build --config Release
```

## Controls

- Click **SETTINGS** (or press `Esc`) for **Quality** (Low / Medium / High / Auto), air solver rate, 20/30 Hz sim, solver toggles, **simulation threads**, scenes, and the rest of the old keyboard options. Low drops pressure/substep/limiter caps, turns air off, uses flat drawing, and runs at 20 Hz so weaker machines stay smooth.
- Use the top-bar view buttons (**NORM**, **CHNK**, liquid **FILL/LIQP/LVEL/LDIV**, rigid **RGDN/RGDO**, gas **GASP/GASA/GASV**). Hover a button to read what it does in **INSPECTOR**.

Everyday painting:

- Pick **FLUIDS / SOLIDS / TOOLS / MISC** in the bar under the canvas. The right-hand list shows that tab's elements (water, wood/stone, erase, wall). Gases / plasma / energy are empty placeholders.
- **SOLIDS** (wood / stone) fall unless **Anchored** is on in Properties. Anchored bodies stay put and still collide.
- **WALL** is under **MISC** (static world cells).
- Mouse wheel or `[` / `]`: brush size.
- **PAUSE** / **RESUME**, or `Space`: pause. While paused, **Step** in Settings advances one tick.
- Top `-` / `+`, or `-` / `+` keys: playback speed (`0.1x`–`4x`). Physics step stays 1/30 s.
- Right-click a rigid body to delete it. Hover a cell to see it in **INSPECTOR**.
- The search field is a placeholder (no filtering yet). The bottom console bar logs simple `[PHYS]` lines; it does not run commands.

The menu covers:

- Clear / reset world, slosh impulse
- Views: Normal, Fill, Pressure, Velocity, Divergence, Chunks, Rigid
- Glow, rigid overlay, wood/stone draw material
- Vorticity, advection (SL / BFECC), residual merge, pressure iterations
- Fluid test scenes and rigid test scenes

## Implemented model

- Fractional cell-centered liquid fill; rendering remains one crisp simulation pixel per world cell.
- Staggered MAC velocity arrays: horizontal velocity on vertical faces and vertical velocity on horizontal faces.
- Semi-Lagrangian velocity advection with optional extrema-clamped BFECC, gravity, solid no-through boundaries, divergence calculation, and red-black Gauss-Seidel pressure projection.
- Previous-frame pressure is retained as a temporal warm start; free surfaces use approximate fractional face weights and adaptive 8/14/24-iteration pressure budgets with residual early-out.
- Iteratively donor/receiver-limited face fluxes allow simultaneous `A -> B -> C` through-flow while keeping volume bounded and scan-order independent.
- Conservative residual-volume consolidation gathers sub-visible remnants into physically plausible neighboring liquid rather than deleting them or letting them form map-wide films.
- Fixed 1/30-second physical updates with CFL-driven fluid substeps. The toolbar changes playback speed, not physics behavior.
- A persistent 3x3-smoothed color field drives near-surface masks, normals, curvature, coherent surface tension, and physically gated spray breakup.
- Water uses an intentional scale of 4 cells per meter and 9.81 m/s² gravity (39.24 cells/s²). Its low viscosity skips diffusion; a real spatial-diffusion path is available for future thicker liquids.
- Conservative ballistic droplets with volume and momentum transfer on re-entry.
- Open world edges act as a void: outward-flowing liquid and splash particles are removed instead of being reflected or deposited back inside.
- A `LiquidProperties` structure provides density, viscosity, and surface tension for future chemistry materials without implementing them yet.
- 16x16 chunks sleep from measured velocity/divergence/fill/pressure stability. Sleeping chunks stay asleep through leftover MAC speed until fill actually changes or liquid goes airborne. Disturbances wake a halo and the connected fluid region; advection, forces, projection, and transport operate on its bounding region.
- Live diagnostics include thin-cell count/volume, momentum, kinetic energy, hovered-cell state, pressure iterations, CFL cap, and averaged stage timings.
- Run `build\water_test.exe --benchmark` to generate a deterministic TSV report beside the executable.
- Run `build\water_test.exe --thread-benchmark` to compare simulation-step timings at 1 / 2 / 4 / Auto workers.
- Run `build\water_test.exe --rigid-benchmark` for a short rigid-body conservation/timing report.
- Run `build\water_test.exe --gas-diag` for Air conservation / chamber regression rows in `gas_diag.tsv`.
- Run `build\water_test.exe --solid-diag` for grab / fracture / moisture regression rows in `solid_diag.tsv`.

The old falling-sand water mover, per-cell water velocity, lateral target search, directional scan bias, and individual-cell sleeping system have been removed.
