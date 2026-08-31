# Goal

Establish monthly Earth validation inputs and replace non-conserving moisture shortcuts with physical saturation and an auditable water budget.

# Phase

1 — Earth reference data and conservative water cycle.

# Scope

- Select authoritative monthly precipitation and temperature references.
- Add reference-data acquisition and metadata without committing bulky rasters.
- Use physical saturation vapor pressure in the moisture solver.
- Remove silent atmospheric water loss and direct inland moisture creation.
- Add water-budget diagnostics and deterministic checks.
- Record Release benchmark comparisons.

# Constraints

- Preserve fixed-seed determinism.
- Keep observed climate fields validation-only.
- Avoid new production dependencies.
- Retain a change only with an explained physical and benchmark result.

# References

- `docs/tasks/climate-physics/PLAN.md`
- `climate_ocean.cpp`
- `climate_validation.cpp`
- `generation_tuning.hpp`

# Steps

1. Audit existing environment/configuration and reference import paths.
2. Select and document monthly Earth data sources and licensing.
3. Add physical humidity/saturation functions and water-budget accounting.
4. Replace inland maritime injection with conservative transport.
5. Benchmark and inspect spatial/class changes.
6. Commit the verified phase or record the blocker.

# Done when

The solver reports a bounded water-budget residual, uses temperature-dependent physical saturation, and no longer creates or destroys atmospheric water through generic carry or inland-source heuristics.

# Verification

- `cmake --build --preset build-x64-release`
- Fixed-seed Release Earth benchmark.
- Water-budget diagnostic checks.
- Workbook and atlas inspection.

# Open questions

- The conservative solver is physically closed but still under-transports precipitation into continental interiors; circulation and energy forcing are the next upstream constraints.

# Final state

- WorldClim monthly temperature and precipitation are available as validation-only inputs.
- Atmospheric moisture transport conserves water to better than `4e-7` relative residual in the Earth benchmark.
- Run 22 preserves the run 21 simulation and adds area-weighted validation metrics.
