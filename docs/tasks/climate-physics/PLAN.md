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
7. Cell-centred grid, positive finite-volume transport, and renderer integration.
8. Diagnosed heating and mode-separated surface/upper circulation.
9. Prognostic climate and replayable gameplay weather anomalies.
10. Wind-driven ocean, coupled SST, validation, and cleanup.

# Shared constraints

- Observed Earth data is validation-only, never a hidden input to generated worlds.
- Preserve deterministic fixed-seed generation.
- Maintain explicit energy and water budgets.
- Implement the full roadmap before resuming Release benchmarks (user override for phases 7-10).
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

Phases 1-6 established the retained checkpoint and benchmark pipeline. Phases 7-10 implementation is complete on `climate-physics_7_core_feat_260902/future-work`. Release build, 13 deterministic test suites and generated-world/save smoke checks pass. First post-implementation Earth run 146 and static visual analysis are complete: rainfall realism regressed sharply despite improved upper-air/ascent errors and aggregate Koppen agreement. No tuning was applied. Evidence, reproduction and diagnostic limitations: `extra/climate/benchmarks/analysis/146/REPORT.md`. High-resolution, repeated, multi-world and long-duration validation remain outstanding.
