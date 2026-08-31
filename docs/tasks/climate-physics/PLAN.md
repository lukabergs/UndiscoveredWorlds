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
- `extra/climate.xlsx`
- `climate_benchmark_runs.json`

# Status

Phases 1-4 complete. Phase 5 is validating run 118's circulation and moisture transport before retaining further physics changes.
