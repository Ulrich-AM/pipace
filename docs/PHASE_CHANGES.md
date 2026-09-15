# Phase changes — design note (not implemented)

This is the next milestone after the Substance / MatterPhase architecture pass.
**No phase-change solver lives in the code yet.** Do not treat this file as an API.

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

## Suggested first implementation plan (later)

Do **not** implement this in the architecture pass.

1. Document mass conversion: full liquid cell ↔ gas amount ↔ solid pixel mass
   using water's phase-specific densities (fluid.density vs mechanical.densityRel
   vs gas EoS). Water currently has **no solid mechanical table** and **no water
   vapor EoS** — those must be added as water's solid/gas properties, still
   under `SUBSTANCE_WATER`, before transfer can be physical.
2. Thermal plateau: while `T` is at the transition point, incoming heat pays
   latent cost instead of raising temperature.
3. Transfer liquid → gas in a single cell first (boiling at 1 atm), conserving
   mass and energy, then update `sampleMatterAt`.
4. Transfer gas → liquid (condensation) with the same conservation.
5. Solid ⇄ liquid after a water solid representation exists (rigid pixel and/or
   static ice wall). Do not fake ice as stone.
6. Only then consider partial-cell fractions and pressure dependence.

## What not to do

- Do not add Ice/Steam SubstanceIds.
- Do not silently convert honey, wood, or air as a side effect of the water
  path.
- Do not add dead `transferLiquidToGas()` stubs until the first real transfer
  is being written.
- Do not change occupancy order ad hoc to hide two-phase cells.

See `AGENTS.md` for the current identity/engine mapping and
`misc/PIPACE_Project_Roadmap_TODO.txt` for the broader roadmap.
