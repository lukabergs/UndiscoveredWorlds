# Goal

Implement the full climate roadmap before benchmark tuning.

# Phase and final state

Phases 7-10 implementation complete on `climate-physics_7_core_feat_260902/future-work`.
The branch remains local and unmerged. Implementation and the requested first
post-implementation benchmark analysis are complete. This plan is collapsed.

# Scope

Shared cell-centred grid/remapping; spherical MPDATA; diagnosed heating;
independent zonal/surface/upper closures; scaling and solver diagnostics;
prognostic climate; replayable runtime weather/statistics/LOD;
wind-driven ocean and coupled SST; adaptive RK2 trails; both generation pipelines.

# Constraints

Observed Earth data remains validation-only. No production dependencies added.
World save version 17 intentionally breaks version-16 compatibility.
Unrelated `outputs/` content is untouched. No tuning during this evaluation pass.

# Benchmark result

Earth run 146 (512x257, seed 1) completed in 271.744 seconds. Common-mask IMERG
land rainfall RMSE regressed 780.30 -> 3109.27 mm/year; upper-air/ascent errors
improved and temperature was broadly unchanged. The maps reveal excessive
equatorial rain concentration and subtropical drying. Native spatial/thermal
metrics have reference-size gaps; independent aligned Koppen agreement improved
27.49% -> 32.96%, but does not imply improved overall climate realism.

# References

- `docs/tasks/climate-physics/FUTURE_WORK.md`: contracts and remaining validation.
- `tests/climate_*_tests.cpp`: deterministic numerical and pipeline fixtures.
- `out/build/x64-Release/future-work-final.uww`: disposable final smoke save.
- `extra/climate/benchmarks/analysis/146/REPORT.md`: full numerical/visual analysis, exact commands and measurement caveats.
- `extra/climate/benchmarks/analysis/146/analyze.py`: reproducible common-mask calculations and plots.
- `extra/climate/benchmarks/runs/146/`: recorded native diagnostics.

# Verification

- `git diff --check`: passed.
- `cmake --build out/build/x64-Release --config Release --clean-first --parallel 4`: passed.
- Final targeted rebuild: `cmake --build out/build/x64-Release --config Release --target UndiscoveredWorlds climate_weather_tests`: passed.
- `ctest --test-dir out/build/x64-Release -C Release --output-on-failure`: 13/13 passed.
- `out/build/x64-Release/bin/UndiscoveredWorlds.exe --help`: passed.
- Final smoke: `out/build/x64-Release/bin/UndiscoveredWorlds.exe --generate-world --seed 73021 --resolution 32 --plate-cycles 1 --cycle-steps 1 --plates 4 --no-rivers --no-lakes --no-deltas --no-profile-log --save out/build/x64-Release/future-work-final.uww`: passed.
- Final smoke used prognostic samples without fallback; all four final ocean solves converged; maximum absolute seasonal relative water residual was 4.01e-8.
- Earlier 64x33 generation smoke also completed; latest 64x32 prognostic spin-up is covered by a deterministic 30-sample fixture.
- Save header and embedded weather payload inspected; bit-exact replay/invalid-state rejection covered by unit tests.
- Benchmark evaluation: Release target rebuild, Earth run 146, analysis script, workbook row/manifest inspection and static map inspection passed. CTest rerun: 13/13; `git -c core.safecrlf=false diff --check` passed.

# Open questions

First Earth benchmark and static map comparison completed. Repeated runs,
high-resolution/multi-seed, evolving-weather/ocean-current visual validation and
long-duration checks remain unverified. Snow mass is not fully equilibrated;
weather uncertainty lacks autocorrelation adjustment. Native reference alignment
and diagnostic gaps are recorded in the report. Ocean iteration caps can be
reached on other worlds. Full application save/load/runtime-weather round trip is unverified.
Existing MSBuild shared-intermediate-directory warnings remain.
