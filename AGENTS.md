# PIPACE — Agent Context

Read this before changing the codebase. Prefer this file plus the roadmap/notes over rediscovering the project from scratch.

## What PIPACE is

PIPACE is a **pixel-native physics/chemistry sandbox** (Win32 + C++17, no external deps). The current prototype is a **free-surface liquid engine** (water), not the full game yet.

Long-term architecture (cooperating systems, not one universal solver):

| Responsibility | Answers |
|---|---|
| **Physics engines** (liquid / solid / gas) | How does matter with these properties move? |
| **SACE** (Somewhat Accurate Calculator for Elements) | What is this substance, what can it become, what properties should new substances have? |
| **Chemistry** (`ReactionDefinition` / `ReactionEngine`) | What reactions are possible, and a conservative local execute step. Player-facing reactions: `2 H2(g) + O2(g) -> 2 H2O(g)` (`HomogeneousCell`, ignition ~850 K) and `C(s) + O2(g) -> CO2(g)` (`SolidGasSurface`, ignition ~900 K). `SUBSTANCE_AIR` is not decomposed into O2/N2; Carbon does not burn in ambient Air. |
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
    GasEngine.h/.cpp       # Generic gas composition + amount/pressure/flow (P still isothermal)
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
    RuntimeSubstance.h     # Phase 18 packed RuntimeSubstanceRef (not SubstanceId)
    GeneratedMaterialRegistry.h/.cpp # Phase 18 generated compiled-profile cache (session-local)
    RuntimeSubstanceProperties.h/.cpp # Phase 19A built-in/generated read-only property bridge
    RuntimeComponentDiagnostics.h/.cpp # Phase 19B1 runtime component payload validation
    PhaseTransfer.h/.cpp  # fill/mass/gas-amount conversion
  world/
    WorldQuery.h/.cpp      # sampleMatterAt (SubstanceId + phase), phase/registry diags
    PhaseChangeEngine.h/.cpp # generic live Liquid ⇄ Gas and Solid ⇄ Liquid by SubstanceId
    WaterPhaseChange.h/.cpp # Water diagnostics / compatibility wrappers
  chemistry/
    ReactionTypes.h / ReactionRegistry.h/.cpp  # ReactionId + ReactionDefinition data only
    ReactionMatterAccess.h/.cpp                # SubstanceId+phase query; liquid commit
    ReactionEngine.h/.cpp                      # local liquid/gas reactions via moles; H2+O2 combustion live
    SaceTypes.h / SaceIdentity.h/.cpp          # Phase 1 exact chemical identity + elemental composition
    SaceCatalog.h/.cpp                         # Phase 2 in-memory generated identity catalog (not live SubstanceId)
    SaceProperties.h/.cpp                      # Phase 3 property provenance + identity-derived molar mass
    SaceMolecule.h/.cpp                        # Phase 4 small-molecule graphs (not canonical identity)
    SaceDescriptors.h/.cpp                     # Phase 5 trusted graph binding helpers + structural descriptors
    SaceFunctional.h/.cpp                      # Phase 6 H/C/O functional motifs (not identity, not properties)
    SaceEstimation.h/.cpp                      # Phase 7–16 Joback thermo + viscosity + Sastri-Rao sigma + Sato-Riedel/Gharagheizi k (StructuralEstimate/Low)
    SaceSimulationReadiness.h/.cpp             # Phase 13–16 generated thermo vs live simulation preflight (not spawnable)
    SaceSimulationCompiler.h/.cpp              # Phase 17 generated engine-profile compiler (not spawnable, not SubstanceId)
  sim/
    SimulationQuality.h/.cpp                   # SimulationQualityProfile + tick scheduler
  docs/
    PHASE_CHANGES.md       # Live Liquid ⇄ Gas and Solid ⇄ Liquid are generic by SubstanceId
  CMakeLists.txt / build.bat / run.bat
  README.md
  AGENTS.md
  misc/                   # notes, changelogs, diagnostic TSV dumps
```

Physical properties are canonical on `SubstanceDefinition` grouped structs
(`mechanical` / `fluid` / `thermal` / `phase` / `porous` / `chemical`). Do not add new
duplicate property tables.

**Identity vs composition vs properties vs phase vs engine:**

- `ChemicalIdentity` = what the substance is / contains (elemental counts, representation kind, optional structure key). Formula text is reference/display only and is not unique identity.
- `ChemicalProperties` = physical/reactive metadata (`molarMass`, flammable, oxidizer). Separate from identity so estimated properties cannot overwrite composition.
- `SubstanceDefinition` = engine built-in record (SubstanceId + grouped properties + chemical identity)
- `MatterPhase` = current represented phase (WHICH phase this cell/body is)
- `SubstanceId` = compact engine identity (built-in domain; no runtime-generated IDs yet)
- `SaceRecordId` = SACE catalog handle for generated identities. Not a `SubstanceId`. Not valid for FluidEngine / GasEngine / RigidBodyEngine / ThermalEngine / ReactionEngine.
- `RuntimeSubstanceRef` = packed 32-bit engine handle (built-in SubstanceId or generated runtime handle). Not canonical identity. Not stored in FluidEngine/GasEngine yet.
- player-facing name (`displayName` / `displayNameKey`) = metadata only; never canonical chemistry
- generated `Element #HEX` = presentation only; not chemical identity
- engine storage = numerical representation (HOW that phase is simulated today)

SACE Phase 1 (`chemistry/SaceIdentity`) authors exact identities for Water, Hydrogen, Oxygen, Carbon, and CO2. Honey, Air, Wood, Stone, Glass, and Metal stay honest unknowns/mixtures/composites. Atom-balance checks are diagnostic/registration-time, not per-cell chemistry.

SACE Phase 2 (`chemistry/SaceCatalog`) is an in-memory canonical cache of generated identities. Built-in exact matches always win. Generated records own their strings, are not spawnable, and are session-local (clear restarts display ordinals). Do not call the catalog from physics ticks. Catalog records live in a `std::deque` so insertion does not invalidate pointers to existing records.

SACE Phase 3 (`chemistry/SaceProperties`) adds property provenance (`SacePropertySource` / `SaceConfidence`) and the first identity-derived property: molar mass from exact elemental composition and a tiny H/C/O atomic-mass table. Identity may be exact while properties remain unknown. Properties never enter the canonical signature. Lower-quality estimates must not silently overwrite higher-quality known data (`Reference` > `IdentityDerived` > structural/empirical > mixture > fallback). Equal source rank uses confidence (`Unknown` < `Low` < `Medium` < `High`); an exact rank+confidence tie keeps the existing value. Do not copy Water as a default for unknown generated matter.

SACE Phase 4 (`chemistry/SaceMolecule`) stores optional small-molecule atom/bond graphs as supporting identity data. The molecular graph does **not** replace `structureKey` canonicalization. Graph atom order is not chemical identity.

SACE Phase 5 (`chemistry/SaceDescriptors`) derives order-independent structural descriptors from a trusted graph. Descriptors are **not** canonical identity. Graph-to-`structureKey` association is trusted provenance (claimed key must equal the record key) until general molecular-graph canonicalization exists. Graph, binding key, descriptors, and functional profile commit atomically. Attachment does not change catalog signatures or display IDs. Optional Joback boiling-point metadata may be assigned only through property precedence.

SACE Phase 6 (`chemistry/SaceFunctional`) derives H/C/O motifs (hydroxyl C–O–H, ether C–O–C, C=O count) and hydrogen-bond donor/acceptor classification from graph connectivity. The functional profile is graph-derived metadata, not canonical identity and not a physical-property estimate. Water is not an organic hydroxyl. O=O is not an ordinary H-bond acceptor. Graphs with elements other than H/C/O store `supported = false` without invented motif counts.

SACE Phase 7 (`chemistry/SaceEstimation`) is the first structure-derived numerical physical property: normal boiling point from a documented Joback–Reid 1987 subset. `Tb [K] = 198.2 + SUM(group contributions)` using published coefficients for `-CH3`, `-CH2-`, `>CH-`, `>C<`, alcohol `-OH`, and non-ring `-O-` only. Applicability is acyclic, net-neutral, fully single-bonded H/C/O molecules completely coverable by those six groups. Source is `StructuralEstimate`, confidence `Low`. Unsupported chemistry (Water, methane/CH4, H2, O2, CO2, rings, charge, N, unsaturation, incomplete coverage) stays `Unknown` — no fallback formula. Joback failure must not fail trusted graph attachment. The estimate is not canonical identity, not live `PhaseProperties`, and is not copied into spawnable generated matter. Same-formula isomers (C2H6O A vs B) may receive different estimates. A later better estimator or reference value should replace Joback through existing precedence.

SACE Phase 8 extends the same six-group Joback subset to critical temperature, critical pressure (stored in Pa), and critical molar volume (stored in m³/mol). Tc uses the Joback-estimated Tb from the same fragmentation, not a later stored/reference Tb. Same-graph reattach may backfill missing estimates through property precedence and must not downgrade higher-quality data. Critical molar volume is not ordinary liquid molar volume; Tc is not Tb.

SACE Phase 9 (`chemistry/SaceEstimation`) adds a Lee–Kesler (1975) acentric factor and saturation vapor pressure Psat(T) from the same self-consistent Joback Tb/Tc/Pc tuple (`StructuralEstimate` / `Low`). Psat is a calculation helper, not live phase physics. Independently replaced generated properties (for example Reference Tc with Structural Pc/omega) must not be treated as a coherent corresponding-states set. There is no ordinary saturation curve above Tc; Psat(Tc) returns Pc. No Antoine coefficients are invented. The public catalog API does not expose mutable generated records.

SACE Phase 10 extends the same Joback subset with enthalpy of vaporization at Tb (`ΔHvap [kJ/mol] = 15.30 + SUM(groups)`, stored as J/mol) and Watson temperature scaling with exponent 0.38. The Watson model is built from one Joback bundle plus identity-derived molar mass, not mixed stored properties. Hvap(Tc) = 0. Values below Tb are extrapolations; generated melting/triple-point stability is unknown. J/kg conversion is a helper, not a stored scalar. Not live `PhaseChangeEngine` latent heat.

SACE Phase 11 adds Joback–Reid ideal-gas Cp(T) = A+BT+CT²+DT³ (298–1000 K, no extrapolation) and Rowlinson–Poling saturated-liquid Cp from the same coherent Joback Cp + Joback Tc + Lee–Kesler omega. Liquid Cp fails closed at Tr ≥ 0.98. Alcohols (C2H6O Graph A) are a difficult associating case; no empirical correction. 298.15 K molar reference scalars are `StructuralEstimate` / `Low`. Mass-specific Cp is J/mol/K divided by kg/mol, not a stored scalar. Not live `ThermalEngine` / `specificHeat` / `gasSpecificHeat`. Unsupported chemistry stays Unknown.

SACE Phase 12 adds COSTALD saturated-liquid molar volume Vs(T) and density rho(T) = M / Vs (`StructuralEstimate` / `Low`). Inputs are one coherent Joback path: Joback Tc, Joback Vc used as characteristic-volume V* fallback (not fitted COSTALD V*), SACE Lee–Kesler omega (not fitted omega_SRK), and identity-derived kg/mol. Strict domain `0.25 < Tr < 0.95`; no clamp or extrapolation. Result is saturation density at Psat(T), not compressed-liquid rho(T,P) and not `FluidProperties::density` (sandbox-relative water=1). Joback Vc is not ambient liquid molar volume. Store only 298.15 K kg/m³ when that T is in range. No alcohol correction. Not live FluidEngine.

SACE Phase 13 (`chemistry/SaceSimulationReadiness`) is a computed preflight for generated records only. It distinguishes gas/liquid thermodynamic readiness and liquid–gas equilibrium readiness from live GasEngine/FluidEngine/phase-change readiness. Built-in `ThermalProperties` currently has one general `conductivity` plus `solidConductivity`; readiness still names liquid vs gas conductivity separately for a later compiler. Physical SACE kg/m³ must not be written into sandbox-relative `FluidProperties`. Readiness is not identity, is not stored on the catalog record, and is not for physics ticks. Generated records remain unspawnable.

SACE Phase 14 adds Joback liquid dynamic viscosity `mu [Pa*s] = MW[g/mol] * exp(A/T + B)` from the same six groups (`StructuralEstimate` / `Low`). A = SUM(mu_a) − 597.82, B = SUM(mu_b) − 11.202. MW stays numerical g/mol in the Joback equation (not kg/mol). Output is physical Pa·s, not `FluidProperties::viscosity` sandbox units; no FluidEngine mapping. The source does not publish a universal validity temperature interval; Phase 14 does not invent one. The stored 298.15 K scalar is gated by a successful COSTALD liquid-density evaluation at that T (domain gate only; density does not enter μ). Graph A/B no longer miss readiness `LiquidViscosity`. Unsupported chemistry stays Unknown. No ethanol correction.

SACE Phase 15 adds Sastri–Rao (1995) liquid surface tension from the same coherent Joback Tb/Tc/Pc plus the SACE functional profile (`StructuralEstimate` / `Low`). General organic: `K=0.158`, `x=0.50`, `y=-1.5`, `z=1.85`, `m=11/9`. Alcohol-like (`supported` profile with `hydroxylCount > 0`): `K=2.28`, `x=0.25`, `y=0.175`, `z=0`, `m=0.8`. No acid branch. `Pc` enters as bar (`Pa × 1e-5`); the correlation returns mN/m and is stored/queried as physical N/m. `σ(Tc)=0`. A mathematical value below an unknown melting/triple point is not proof of liquid stability. Physical N/m is not `FluidProperties::surfaceTension` sandbox units; no FluidEngine mapping. Water is not an alcohol and is not routed through this estimator.

SACE Phase 16 adds Sato–Riedel liquid thermal conductivity and Gharagheizi gas thermal conductivity in W/(m·K) (`StructuralEstimate` / `Low`) from one coherent Joback path. Sato–Riedel uses Joback Tb/Tc plus identity-derived MW[g/mol]; T < Tc only. The stored 298.15 K liquid scalar is gated by COSTALD ρ(298.15) as a liquid-domain check (density does not enter kL). Gharagheizi uses Joback Tb/Pc, Lee–Kesler omega, and MW; internal `P = Pc[Pa] × 1e-4` is the published/corrected scale, not ordinary bar. No associating correction. Not live `ThermalEngine` / `ThermalProperties::conductivity`. Graph A/B become preflight `liveGasReady` / `liveLiquidReady` / `liveLiquidGasPhaseChangeReady` while remaining **unspawnable**.

SACE Phase 17 (`chemistry/SaceSimulationCompiler`) compiles a generated record plus coherent models plus readiness into a `SaceCompiledSimulationProfile`. Physical properties stay separate from identity. Density maps kg/m³ → water-relative (`WATER_DENSITY_KG_M3`). Molar Cp maps to J/(kg·K). Thermal conductivity stays physical W/(m·K) and **liquid k and gas k remain separate fields/models** (not `ThermalProperties::conductivity`). Viscosity (Pa·s) and surface tension (N/m) use explicit Water-relative solver calibration anchors (`kPhysicalWaterViscosityAt298KPaS`, `kPhysicalWaterSurfaceTensionAt298KNPerM`); those anchors are not material fallbacks. Joback viscosity A is the compiled Arrhenius shape. Original T-dependent models are retained. `profile.valid` means at least one coherent subprofile compiled. Generated records remain unspawnable; `SUBSTANCE_COUNT` remains 12; no generated `SubstanceId` or world integration.

SACE Phase 18 (`substance/RuntimeSubstance.h`, `substance/GeneratedMaterialRegistry`) adds a packed 32-bit `RuntimeSubstanceRef` that can name a built-in `SubstanceId` or a generated runtime handle. Generated `SaceRecordId` stays chemistry identity and is **not** packed into the handle (handles are never reused in-process). `SaceCatalog::clear()` resets the generated runtime registry so stale refs cannot alias a later record. Registration compiles via Phase 17 and is idempotent; it is not player spawnability. Built-ins are not copied into the generated registry. Liquid/gas component storage, `MatterIdentity`, and world queries remain `SubstanceId`. `SUBSTANCE_COUNT` remains 12.

Phase 19A (`substance/RuntimeSubstanceProperties`) adds a read-only property query bridge over `RuntimeSubstanceRef`. Built-ins delegate to existing `SubstanceDefinition` tables. Generated liquid/gas refs evaluate the retained Phase-17 SACE models at query temperature (COSTALD density, Rowlinson-Poling/Joback Cp, Sato-Riedel/Gharagheizi conductivity, compiled Joback solver viscosity, Sastri-Rao surface tension calibration). Generated queries fail closed on model-domain failure and never borrow Water/Air material properties. This phase does **not** migrate FluidEngine/GasEngine component storage or make generated matter spawnable.

Phase 19B1 migrates `LiquidComponent` and `GasComponent` payload/view identity to `RuntimeSubstanceRef` while keeping FluidEngine/GasEngine internal SoA arrays built-in `SubstanceId` for one more phase. Built-in compatibility overloads preserve existing callers. Standalone generated payloads can add/find/merge/copy/scale/compact and liquid/gas mixture-property helpers consume the Phase-17/19A generated profiles. Live engine commits explicitly reject generated refs until 19B2 rather than truncating or relabeling identity.

Phase 19B2a migrates FluidEngine's internal liquid composition ID arrays (`liquidCompId`, `nextCompId`) to `RuntimeSubstanceRef`. Public liquid mutation/query APIs remain built-in `SubstanceId`-based in this subphase, and conservative transport explicitly unwraps built-ins before calling those APIs. Generated runtime refs are still not commit/transport-enabled inside FluidEngine until 19B2b. Water/Honey behavior must remain unchanged.

These mappings are **implementation, not laws**:

| Current representation | Engine |
|---|---|
| Water + Liquid, Honey + Liquid | FluidEngine (water/honey volume channels) |
| Wood/Glass/Metal/Carbon + Solid | RigidBodyEngine (MaterialId masks) |
| Water + Solid | RigidBodyEngine (`MATERIAL_WATER_SOLID` mask, still SUBSTANCE_WATER) |
| Water + Solid (sub-pixel pending) | `solidifyPendingKg` / `solidifyPendingHeatJ` until a rigid pixel exists. Carries mass and sensible energy and participates in normal ThermalEngine conduction as an extra reservoir, not the primary cell node. |
| Stone + Solid | rigid body **or** static `solid[]` walls (`kStaticWallSubstance`) |
| Air + Gas | GasEngine (`SUBSTANCE_AIR` component) |
| Hydrogen + Gas | GasEngine (`SUBSTANCE_HYDROGEN` component) |
| Oxygen + Gas | GasEngine (`SUBSTANCE_OXYGEN` component) |
| Carbon Dioxide + Gas | GasEngine (`SUBSTANCE_CARBON_DIOXIDE` component) |
| Water + Gas | GasEngine (`SUBSTANCE_WATER` gas component; same SubstanceId) |

Static `solid[]` walls use Stone thermal/mechanical identity but are
**moisture-inert** (no absorb / drip / wetness). Porous Stone rigid bodies still
absorb. Condensation may form liquid in adjacent free cells, not inside the wall.

`supportsPhase` is capability metadata. Live **Liquid ⇄ Gas** and **Solid ⇄ Liquid**
are generic (`world/PhaseChangeEngine.cpp`) for substances whose metadata supports
both endpoints; Water is currently the only eligible built-in. Do **not** add
SUBSTANCE_ICE / SUBSTANCE_STEAM.
Pending sub-pixel solid matter carries mass and sensible energy and participates
in normal ThermalEngine conduction as an extra reservoir (not the primary cell
node). Density stays phase-specific: mechanical.densityRel (solid), fluid.density (liquid),
gas amount/EoS (gas). Honey/water mixtures do not boil or freeze yet.

Water/honey cells that hold both channels are **mixtures**, not a new SubstanceId.
High-level composition queries use `FluidEngine` SubstanceId APIs
(`liquidComponentAmount` / `liquidComponentFraction` / `dominantLiquidSubstance` /
`liquidComponents`). Storage is a fixed-capacity SoA of up to 4 components per
cell (`liquidCompId` / `liquidCompAmt` / `liquidCompCount`); `fill[]` remains
total occupancy. `MatterSample` reports the dominant component plus fractions.

Gas cells use the same pattern: `amount[]` is total cell-atmospheres; composition
is a 4-slot SoA (`gasCompId` / `gasCompAmt` / `gasCompCount`). Air is an explicit
`SUBSTANCE_AIR` component. Hydrogen, Oxygen, and Carbon Dioxide are player-facing
gases (`SUBSTANCE_HYDROGEN` / `SUBSTANCE_OXYGEN` / `SUBSTANCE_CARBON_DIOXIDE`)
painted from the Gases category. Ambient Air is a mixture identity and is **not**
oxygen; Carbon combustion requires explicit `SUBSTANCE_OXYGEN`.
Water vapor is `SUBSTANCE_WATER` + `MatterPhase::Gas`.
Do not infer missing gas as Air. Pressure is `P = (amount/volume)*(T/T_amb)` with
`T_amb = AMBIENT_TEMPERATURE_K`; `pressure[]` is the cached result.

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
- Volume / momentum / KE diagnostics; internal F1–F10 diagnostic fixtures (not player-facing); stage timings
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

- Quality presets (Performance/Balanced/Accurate/Auto/Custom) flow through `SimulationQualityProfile` (`sim/SimulationQuality.h`). They change solver iterations, substeps, sleep, and thermal/chemistry/phase **update intervals** (accumulated dt, not slower physics). Auto chooses Performance or Balanced only, all subsystems together. Custom composes per-subsystem levels from `profileForQualityLevel` and keeps those selections when leaving Custom. Low/Auto-Low use `GasSimMode::Half`; they never silently set `GasSimMode::Off`. Settings → Gas Off and Thermal Off remain explicit. Quality does not write `physicsHz` or `catchUpTicks`. `PHYSICS_DT` is still 1/30 if the menu Hz is 20 (deferred).
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

**SACE Phase 1 (chemical identity)** is in place. **SACE Phase 2** adds an
in-memory generated identity catalog (`SaceRecordId`); generated records are
not live or spawnable. **SACE Phase 3** adds property provenance and
identity-derived molar mass; identity may be exact while properties remain
unknown. **SACE Phase 4** adds optional small-molecule graphs; they do not
replace `structureKey` canonicalization. **SACE Phase 5** adds trusted
graph-to-structureKey binding and structural descriptors; descriptors are not
canonical identity. **SACE Phase 6** adds H/C/O functional motifs and H-bond
feature classification from graphs; the functional profile is not identity and
not a property estimate. **SACE Phase 7** adds a Joback–Reid 1987 subset
normal-boiling-point estimate (`StructuralEstimate` / `Low`) for acyclic
neutral saturated H/C/O alcohol/ether molecules; unsupported chemistry stays
Unknown and is not live phase data. **SACE Phase 8** extends that subset to
Joback Tc/Pc/Vc (`StructuralEstimate` / `Low`); critical volume is not liquid
density and Tc is not boiling temperature. **SACE Phase 9** adds Lee–Kesler
omega and Psat(T) from the same Joback Tb/Tc/Pc tuple; it is not live phase
integration. **SACE Phase 10** adds Joback ΔHvap(Tb) (stored J/mol) and Watson
Hvap(T) with exponent 0.38; it is not live latent-heat physics. **SACE Phase 11**
adds Joback ideal-gas Cp(T) (298–1000 K) and Rowlinson–Poling liquid Cp; it is
not live ThermalEngine data. **SACE Phase 12** adds COSTALD saturated-liquid
Vs(T)/rho(T) with Joback Vc as V* fallback; it is not live FluidEngine density.
**SACE Phase 13** adds generated thermo vs live simulation readiness. **SACE Phase 14**
adds Joback liquid dynamic viscosity (Pa·s, MW in g/mol). **SACE Phase 15** adds
Sastri–Rao liquid surface tension (N/m) from coherent Joback Tb/Tc/Pc plus
functional hydroxyl classification. **SACE Phase 16** adds Sato–Riedel liquid k and
Gharagheizi gas k (W/(m·K)); Graph A/B become preflight live-ready but remain
unspawnable. **SACE Phase 17** compiles those models into an engine-safe profile
(water-relative density, mass-specific Cp, physical phase-specific k, calibrated
solver viscosity/surface tension) without allocating a SubstanceId or spawning
generated matter. **SACE Phase 18** adds compact `RuntimeSubstanceRef` plus a
session-local generated compiled-profile registry. Registration is not spawnability
and does not occupy FluidEngine/GasEngine cells.
Live `Liquid ⇄ Gas` and `Solid ⇄ Liquid` are generic by
SubstanceId (`world/PhaseChangeEngine.cpp`). Honey mixtures still skip phase
change. Runtime-generated world matter, reaction families, graph isomorphism,
and broader property estimates are not started.

Longer sequence (historical): liquid correctness → performance → modularization →
honey composition → rigid coupling → temperature/gas → **phase changes** → SACE.

Detail lives in `misc/PIPACE_Project_Roadmap_TODO.txt` and `misc/PIPACE_SACE_Concept_Notes.txt`.

## Personality / presentation

Homely, utilitarian, slightly garage-lab — not neon sci-fi. Flavor text must never alter simulation behavior.
