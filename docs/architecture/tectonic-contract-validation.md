# Tectonics implementation and checks

The only current engine is `src/simulations/geology/plate_tectonics`, built from
the working tree as the library `uw_tectonics`. The app and CLI share
`tectonics::scenario::baseline()` compiled from `configs/tectonics-baseline.json`.
Do not configure a separate checkout or external tectonics library.

`tectonic_api.hpp` exposes the current engine to the app; every target is
C++17. `lithosphere::update()` always advances the force-driven
material columns. Time advances continuously to the scenario duration.
World ingestion is `geology/world_tectonics.cpp`.
It runs at recipe resolution and resamples to the planet grid; it is not a
second engine. Source copyright notices and the module LICENSE must remain.

Use the build and test commands in the root README. `PlateTectonicsTests` checks
current motion, transport, mass budgets, topology, field contracts, serialization,
determinism and world-grid registration. `MineralSystemsTests` checks native
displacement/identity ingestion and resource budgets. The CLI check compares
compiled defaults with the explicit checked-in recipe and rejects invalid options.
Use seed 1. Small deterministic tests do not establish physical realism or
full-size performance; preserve recipe, fields, source and toolchain when making
numerical comparisons. Historical run archives are evidence, not selectable
implementations for new generation.
