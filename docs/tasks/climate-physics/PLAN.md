# Goal

Replace the climate synthesizer's dominant empirical shortcuts with a reduced, conservative climate model whose circulation, temperature, and precipitation respond to planetary forcing and terrain.

# Why this is a project

The work changes the temperature, circulation, moisture, ocean, cryosphere, and validation boundaries across multiple reviewable branches.

# Phases

1. Earth reference data and conservative water cycle.
2. Surface energy balance and ocean mixed layer.
3. Two-layer atmospheric circulation.
4. Clouds, snow, runoff, and inter-season storage.
5. Calibration, regression hardening, and cleanup.
6. Benchmark pipeline and artifact organization.

# Shared constraints

- Observed Earth data is validation-only, never a hidden input to generated worlds.
- Preserve deterministic fixed-seed generation.
- Maintain explicit energy and water budgets.
- Benchmark every retained physics change in Release.
- Do not tune Köppen classification boundaries to compensate for upstream errors.

# Shared references

- `climate_ocean.cpp`
- `globalclimate.cpp`
- `generation_tuning.hpp`
- `climate_validation.cpp`
- `extra/climate/workbooks/climate.xlsx`
- `data/climate/benchmark_runs.json`
- `docs/tasks/climate-physics/FUTURE_WORK.md`

# Status

Phases 1-4 complete. Phase 5 retained the surface stationary-wave response and documented rejected experiments. Phase 6 made benchmark resolution, map selection, workbook updates, and artifact locations explicit and reproducible. Deferred and temporarily disabled candidates are tracked in `FUTURE_WORK.md`.
