# Final state

Fresh workbook: `outputs/01a062ac-02ef-7710-bf6a-25bd1f12cf87/climate_validation.xlsx`.
Contains only benchmark 159: seed 1, 64x33 output, 64x32 wind/moisture,
selected physics and 24-pass cap. Coupling converged at pass 13; aggregate
heating change 1.49%. No earlier results migrated; original workbook preserved.

28 sheets cover paired reference values, 87 spatial regions, formula-based
scores, distributions/features, budgets and explicit coverage gaps. Koppen is
global only, using common-mask spherical areas rather than resolution-scaled
pixel counts. Numeric exports and `scripts/summarize-climate-validation.py`
provide reusable inputs, including CSVs ordered for Excel pasting.

# Verification

- `cmake --build out/build/x64-Release --config Release --target UndiscoveredWorlds climate_physics_tests climate_ocean_dynamics_tests climate_atmosphere_tests climate_pipeline_tests climate_weather_tests climate_benchmark_outputs_tests --parallel 4`: passed before integration.
- `ctest --test-dir out/build/x64-Release -C Release -R '^climate_(physics|ocean_dynamics|atmosphere|pipeline|weather|benchmark_outputs)_tests$' --output-on-failure`: 6/6 passed before integration.
- Bundled Python: `tests/test_climate_validation_summary.py`: 5/5 deterministic checks.
- One completed Earth benchmark, `--seed 1 --resolution 64`, imported `out/climate-tuning/earth64/` terrain, `--map none`: exit 0.
- Workbook scalar/confusion reconciliation, four-quarter water closure, formula error/cache inspection, all-sheet visual inspection and ZIP/XML validation passed.
- `git -c core.safecrlf=false diff --check`: passed.

# Remaining limits

Higher-resolution, long-duration and additional parameter/seed runs remain
unverified. Reference gaps and temporal/layer proxies are explicit in Coverage.
Ocean summaries cover cells identified by captured currents/ice; the solver's
native wet mask is not exported, and final terrain is unsuitable as a substitute.
Entirely stationary ice-free wet cells may be omitted. Future runs require new
spatial inputs and formula rows; automatic appending to the new workbook is not
implemented. Physics findings and next steps remain in `FUTURE_WORK.md`.
