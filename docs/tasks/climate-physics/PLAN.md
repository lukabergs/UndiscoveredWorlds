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

Phases 1-6 established the retained checkpoint and benchmark pipeline. Phases 7-10 are implemented on `climate-physics_7_core_feat_260902/future-work`, including slab-ice enthalpy and skin/SST separation. Coarse trials 151-158 improved surface wind agreement with rainfall tradeoffs; higher-resolution validation was interrupted. Fresh workbook baseline 159 completed at 64x33 and converged in 13 coupling passes. Verification and remaining limits: `docs/tasks/climate-physics/climate-physics_7_core_feat_260902/future-work/PLAN.md` and `FUTURE_WORK.md`.
