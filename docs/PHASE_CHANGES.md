# Phase changes — design note

This is the next milestone after the Substance / MatterPhase architecture pass.
Live **Liquid ⇄ Gas** is generic (`world/PhaseChangeEngine.cpp`) and driven by
`SubstanceId` + `PhaseProperties`. Water is currently the only built-in that
meets both-endpoint eligibility. Live **Solid ⇄ Liquid** is still Water-specific
(`world/WaterPhaseChange.cpp`). Honey mixtures do not boil or freeze yet.

Target first: **water** `solid ⇄ liquid ⇄ gas` using the existing `SUBSTANCE_WATER`
id. There will be no `SUBSTANCE_ICE` or `SUBSTANCE_STEAM`.

## Current vs future

| Layer | Role today | Role after phase changes |
|---|---|---|
| `SubstanceDefinition` | Intrinsic identity + properties | Unchanged. Identity does not flip. |
| `MatterPhase` | Current represented phase | Same, but cells/bodies can change phase. |
| Engine storage | How that phase is simulated | Matter **moves between** engines when phase changes. |

Today water is only simulated as liquid (`FluidEngine` volume). Capability flags
already say it can be solid and gas. Those flags are metadata, not transfer.

Mass/fill/gas-amount conversion helpers live in `substance/PhaseTransfer.h`.
Live **Liquid ⇄ Gas** transfer is in `world/PhaseChangeEngine.cpp` (Water
uses this path; Honey/CO2/Carbon are ineligible with current metadata).
Live **water liquid ⇄ solid** (rigid `MATERIAL_WATER_SOLID`, still
`SUBSTANCE_WATER`) is in `world/WaterPhaseChange.cpp`. Honey mixtures do not boil or freeze yet.

Vapor placement: the occupancy model is one primary medium per cell, and liquid
with `fill >= MIN_PRESSURE_FILL` has zero gas volume. Boiling deposits vapor
into a neighboring accessible gas cell (fair spatial tie-break among equal
room). Destination storage has a solver safety cap (`kSolverSafetyAtm`,
currently 10000 atm) that is **not** a thermodynamic law; hits are counted as
`blockedBoilSafetyLimit`. Boiling/condensation follow a Clausius–Clapeyron
saturation curve from `PhaseProperties` (`saturationVaporPressurePa` /
`saturationTemperatureK` in `substance/PhaseTransfer.h`): higher total pressure
raises `T_sat`, lower pressure lowers it. Condensation compares each eligible
gas component's **partial** pressure to `P_sat(T)` for that SubstanceId, not
total gas temperature vs 373 K. Condensate prefers existing nearby liquid of
the same SubstanceId, then cells next to surfaces that are actually cooler
than the vapor, then lower neighbors. Hot walls do not get a
near-solid bonus. Equal destinations cycle by a tick-salted spatial hash.

Phase existence uses `kMinFillMove` (1e-6 fill), separate from the fluid-motion
threshold `MIN_ACTIVE_FILL` (0.01). Residual superheated puddles can still
vaporize. Liquid/gas rates are `kMaxFillPerSec * dt` (time-based).

Closed unforced transfers conserve `E_thermal + m_vapor * L_v` by moving the
Cp difference through the energy ledger (liquid Cp ≠ vapor Cp).

Boiling never cools remaining liquid below `T_m + 1 K` to pay latent heat
(`Tpay = max(T_sat, T_m + 1)`), so leftover liquid is not immediately frozen.

Gas `clampSpecies` keeps water vapor when total amount drops below the wipe
threshold (raises amount to match vapor) instead of deleting it. Condensation
latent leftovers go into dest liquid up to `T_sat + 20 K`, then nearby gas/walls,
never as an uncapped dest-liquid dump.

Static world walls (`fluid.solid[]`, `kStaticWallSubstance`) conduct heat and
block flow but never store moisture. Adjacent condensation is allowed; porous
Stone rigid bodies still absorb.

## Current sandbox units

These are the units the helpers use. They match thermal/gas code; they are not SI
lab apparatus except where noted.

| Quantity | Meaning |
|---|---|
| `cellsPerMeter` | Default 4 → cell edge **0.25 m** |
| Implied depth | One cell, so **V = dx³ = 0.015625 m³** |
| Liquid `fill` | Fractional occupancy of that cell volume. `fill = 1` is one full liquid cell. |
| Liquid density | `FluidProperties.density`, relative to water = 1.0 |
| Full water cell mass | `1 × 1000 kg/m³ × 0.015625 m³` = **15.625 kg** |
| Gas `amount` | Conserved cell-atmospheres. **1.0 amount in 1.0 volume = 1 atm** (isothermal `P = amount / volume`). |
| Air mass | `amount × 1.204 kg/m³ × V` (same as `gasMassKg`) |
| Water vapor mass | Component amount × ρ_vapor × V, `ρ_vapor` from water `chemical.molarMass` at reference P and ambient T. Still `SUBSTANCE_WATER` + `MatterPhase::Gas` in generic gas composition slots. Air is an explicit `SUBSTANCE_AIR` component, not `amount − vapor`. |
| Thermal energy | Joules |
| Latent heats | `PhaseProperties.latentHeatFusion` / `latentHeatVaporization` (J/kg) |

## Implementation status

Live Liquid ⇄ Gas is generic (`world/PhaseChangeEngine.cpp`). Only substances
whose registered capabilities support both liquid and gas endpoints, and whose
saturation/latent/molar-mass data are valid, are eligible. Water is eligible.
Honey, CO2, and Carbon are not (do not change their capabilities to force it).
Live Solid ⇄ Liquid remains Water-specific. Remaining:

1. Mixture thermodynamics (honey/water must not boil/freeze until then).
2. Same-cell liquid/gas occupancy if the one-primary-medium model is relaxed.
3. Rigid/fluid buoyancy may not yet make ice float; do not add a special ice force.
4. Clausius–Clapeyron is a two-parameter approximation near the reference
   boiling point, not a steam table / critical-point model.
5. Generic Solid ⇄ Liquid (Phase Change Generalization - Phase 2).

## Invariants

1. **SubstanceId stays constant** across a phase change. Ice and steam are water
   in a different `MatterPhase`, not new substances.

2. **Mass / amount is conserved.** Liquid fill, rigid/static solid mass, and gas
   amount must convert through an explicit mass-equivalent. Do not invent or
   delete matter to hide a transfer bug.

3. **Energy is conserved, including latent heat.** Temperature must not jump
   through the transition energy. Expect plateau behavior at melting/boiling
   while latent heat is absorbed or released.

4. **Phase change transfers representation between engines.** Example paths for
   water:
   - liquid (`FluidEngine` fill) → gas (`GasEngine` amount)
   - liquid → solid (rigid `MATERIAL_WATER_SOLID`, still `SUBSTANCE_WATER`)
   - solid → liquid
   - gas → liquid
   Representation transfer is the hard part; property lookups are already there.

5. **World queries must report the new phase immediately** after a transfer.
   `sampleMatterAt` occupancy order (rigid → wall → liquid → gas) stays the
   source of truth. Do not leave two representations occupying the same cell.

6. **Partial-cell conversion** should be supported eventually (a cell can be
   partly liquid and partly gas/solid during a transition). The current occupancy
   model is one primary medium per cell; that is a known limit, not a law.

7. **Pressure dependence** uses a Clausius–Clapeyron curve from
   `PhaseProperties.referencePressurePa`, `boilingPointK`, `latentHeatVaporization`,
   and `chemical.molarMass`. Condensation is driven by water-vapor partial
   pressure vs `P_sat(T)`. A destination-array safety cap is not a physical law.

8. **Data-driven from `PhaseProperties`.** Melting/boiling points and latent
   heats already live there (copied from thermal authoring at registry init).
   Do not hard-code 273.15 / 373.15 in the transfer solver. Live Liquid ⇄ Gas
   uses the substance's own table; it does not special-case Water identity
   except compatibility counters and Water shading.

9. **Mixtures stay honest.** A water/honey cell is not a new SubstanceId. Phase
   change of mixtures is out of scope for the first water milestone.

10. **SACE is not required** for built-in water. Chemistry stays separate from
    phase motion.

## What not to do

- Do not add Ice/Steam SubstanceIds.
- Do not silently convert honey, wood, or air as a side effect of the water
  path.
- Do not add dead `transferLiquidToGas()` stubs until the first real transfer
  is being written.
- Do not change occupancy order ad hoc to hide two-phase cells.

See `AGENTS.md` for the current identity/engine mapping and
`misc/PIPACE_Project_Roadmap_TODO.txt` for the broader roadmap.
