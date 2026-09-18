# Retained climate alternatives

These implementations are compiled and covered by the existing deterministic
tests. They are comparison tools; their retention does not establish better
Earth-climate accuracy.

- `thermal_response.*` retains empirical local pressure conversion and the
  nonlocal thermal kernel. The kernel has separate meridional smoothing and
  zonal spreading steps, preserves row means, and reverses its asymmetric
  tropical response with rotation. `runs/registry/climate.json` records its
  introduction in run 119, disabling the kernel in run 135, and replacing the
  pressure conversion with hypsometric local forcing in run 136. Its existing
  `climateatmosphere::thermalSurfacePressureAnomalyHpa` and
  `climateatmosphere::nonlocalThermalResponse` entry points remain available
  through `climate_atmosphere.hpp`. An experiment can call them when constructing
  thermal forcing; the production mode-separated circulation uses local
  hypsometric forcing. `climate_atmosphere_tests` checks pressure sign, rotation
  reversal and spatial response.
- `ocean_feedback.*` isolates the empirical SST-to-pressure/wind closure from
  the wind-driven ocean's transport and heat budgets. It is selected by
  `OceanConfig::oneWay = false` in `solveWindDrivenOcean`; the production adapter
  explicitly sets `oneWay = true` and couples SST through the atmospheric solver
  instead. The replacement is documented under “Wind-driven ocean and SST” in
  `docs/tasks/climate-physics/FUTURE_WORK.md`. `climate_ocean_dynamics_tests`
  retains coupled determinism and convergence checks alongside one-way ocean,
  coastal transport, heat-budget and sea-ice tests. It is useful for isolated
  coupling comparisons, not as evidence of atmospheric skill.

Run the retained alternatives with the normal build and
`ctest --test-dir out/build/x64-Debug -C Debug -R "climate_(atmosphere|ocean_dynamics)_tests" --output-on-failure`.

Uncalled latitude-temperature, prescribed-current/SST and directional rain/
monsoon implementations were removed. The all-land rainfall fallback remains
active in `climate_fields.cpp`; runtime weather and stationary/prognostic
circulation remain in their active modules.
