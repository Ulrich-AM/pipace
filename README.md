# PIPACE

**Purely Inaccurate Physics And Chemical Engine**

PIPACE is a small pixel-based physics sandbox where liquids, gases, heat, rigid bodies, and materials can interact in the same world.

It is inspired by sandboxes where you can simply place things down and see what happens, but the long-term goal is to push the simulation much further with deeper material properties, chemistry, phase changes, and user-created substances.

The name is intentionally honest. PIPACE tries to be believable where it can, approximate where it has to, and avoid pretending every part of the simulation is perfectly accurate.

## What can it do right now?

PIPACE is still in development, but it already includes:

- liquids such as water and honey
- mixing and dyeing liquids
- gases and gas pressure
- temperature and heat transfer
- rigid pixel-built solids
- collisions, fracture, and material strength
- porous materials that can absorb and drip water
- different material properties for wood, stone, glass, metal, and more
- several debug views for seeing what the simulation is doing
- quality and performance settings for weaker computers

Under the hood, materials are treated as **substances** with their own physical properties and supported phases. This is being built so that one substance can eventually exist as a solid, liquid, or gas without being treated as three unrelated materials.

## Where is it going?

Some of the bigger planned systems are:

- melting, freezing, boiling, and condensation
- chemical reactions
- combustion
- more substances and mixtures
- electricity
- more tools for building contraptions
- **SACE**, a system for estimating the properties of substances that are not already built into the game

The idea is to eventually let the sandbox produce things that were never manually added beforehand.

## Running PIPACE

PIPACE currently targets Windows.

Run:

```text
run.bat
```

This builds and launches the program.

You can also build it with CMake if you prefer.

## A note about accuracy

PIPACE is not meant to be a laboratory simulator.

Some systems use real units and physical data, while others use simplified or sandbox-scaled models so the simulation can stay interactive.

The general rule is:

> Exact where known, rule-based where understood, approximate where necessary.

And if something goes catastrophically wrong, PIPACE will probably consider that a learning opportunity.

## Status

PIPACE is actively being developed and a lot can still change.

For now, the best way to understand it is probably to open it, place some things down, and see what happens.
