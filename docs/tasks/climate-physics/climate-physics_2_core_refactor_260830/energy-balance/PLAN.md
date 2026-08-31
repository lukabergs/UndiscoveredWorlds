# Goal

Replace latitude-band temperature construction with a seasonally forced, energy-conserving surface model whose land/ocean contrast emerges from heat capacity.

# Phase

2 — Surface energy balance and ocean mixed layer.

# Scope

- Compute daily-mean orbital insolation from latitude, obliquity, and eccentricity.
- Integrate land and ocean slab temperatures with temperature-dependent albedo.
- Parameterize conservative poleward atmospheric heat transport.
- Preserve all four simulated seasons instead of reconstructing equinoxes from January and July.
- Export energy-budget diagnostics and compare against monthly WorldClim fields.

# Constraints

- The configured world mean temperature remains the greenhouse-forcing control.
- No observed climate field enters generation.
- Keep the fixed-seed benchmark deterministic and Release-only.
- Retain only changes with a bounded energy residual and interpretable benchmark effect.

# References

- `docs/tasks/climate-physics/PLAN.md`
- `globalclimate.cpp`
- `planet.cpp`
- `climate_ocean.cpp`
- `climate_validation.cpp`

# Steps

1. Add pure orbital-insolation and slab-energy functions with unit tests.
2. Generate four seasonal land/ocean temperature profiles from the energy model.
3. Preserve explicit transition-season temperatures through downstream stages.
4. Export an annual energy budget.
5. Run and inspect a Release Earth benchmark.
6. Retain, revise, or revert based on physical and spatial metrics.

# Done when

Temperature uses orbital energy forcing and mixed-layer heat capacity, reports a closed annual energy budget, and improves the monthly seasonal comparison without empirical latitude bands.

# Verification

- Targeted energy-model unit tests.
- Release build and CTest.
- Fixed-seed Release Earth benchmark.
- Workbook and atlas inspection.

# Open questions

- Whether the current SST current-advection pass should feed back into the surface energy model in this phase or the circulation phase.

# Final state

- Daily orbital insolation, snow/ice albedo, linear outgoing longwave radiation, and conservative heat redistribution replace the active latitude-band temperature generator.
- Land and ocean seasonality emerge from separate thermal reservoirs; all four seasonal states survive downstream processing.
- Run 30 closes the annual energy budget and improves area-weighted annual temperature MAE from `7.00 C` at run 22 to `4.50 C`.
- Precipitation remains the dominant blocker: annual area-weighted precipitation correlation is `0.163`, and dry interiors still dominate the Köppen atlas.
