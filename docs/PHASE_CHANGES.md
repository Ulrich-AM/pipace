# Phase changes — design note

This is the next milestone after the Substance / MatterPhase architecture pass.
Water **liquid ⇄ gas** and **liquid ⇄ solid** are implemented. Honey mixtures
do not boil or freeze yet.

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
Live **water liquid ⇄ gas** transfer is in `world/WaterPhaseChange.cpp`.
Live **water liquid ⇄ solid** (rigid `MATERIAL_WATER_SOLID`, still
`SUBSTANCE_WATER`) is in the same file. Honey mixtures do not boil or freeze yet.

Vapor placement: the occupancy model is one primary medium per cell, and liquid
with `fill >= MIN_PRESSURE_FILL` has zero gas volume. Boiling therefore deposits
vapor into a neighboring accessible gas cell (prefer above, then sides, then BFS
within 64 cells). Condensate prefers existing nearby water, then cells next to
solids, then lower neighbors.

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
| Water vapor mass | `amount × ρ_vapor × V`, `ρ_vapor` from water `chemical.molarMass` at reference P and ambient T. Still `SUBSTANCE_WATER`. Stored as `gas.waterVapor`; air = total − vapor. |
| Thermal energy | Joules |
| Latent heats | `PhaseProperties.latentHeatFusion` / `latentHeatVaporization` (J/kg) |

## Implementation status

Water liquid ⇄ gas is live (`world/WaterPhaseChange.cpp`). Remaining:

1. Mixture thermodynamics (honey/water must not boil/freeze until then).
2. Pressure-dependent boiling curve (currently reference pressure / `PhaseProperties`).
3. Same-cell liquid/gas occupancy if the one-primary-medium model is relaxed.
4. Rigid/fluid buoyancy may not yet make ice float; do not add a special ice force.

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
   - liquid → solid (rigid body and/or static `solid[]`, TBD)
   - solid → liquid
   - gas → liquid
   Representation transfer is the hard part; property lookups are already there.

5. **World queries must report the new phase immediately** after a transfer.
   `sampleMatterAt` occupancy order (rigid → wall → liquid → gas) stays the
   source of truth. Do not leave two representations occupying the same cell.

6. **Partial-cell conversion** should be supported eventually (a cell can be
   partly liquid and partly gas/solid during a transition). The current occupancy
   model is one primary medium per cell; that is a known limit, not a law.

7. **Pressure dependence** may start as reference pressure
   (`PhaseProperties.referencePressurePa`, 1 atm). Clapeyron / vapor-pressure
   curves are later work.

8. **Data-driven from `PhaseProperties`.** Melting/boiling points and latent
   heats already live there (copied from thermal authoring at registry init).
   Do not hard-code 273.15 / 373.15 in the transfer solver.

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
