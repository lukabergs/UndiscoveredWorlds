# Goal

Complete the approved climate-engine integrations before climate-skill tuning.
Gameplay weather and runtime storm systems remain deferred.

# Final checkpoint

Production integrations and verification are complete for this pass. This is
not a claim of converged Earth climate or complete ocean thermodynamics.
Pre-polish code is preserved at `5ae36fa`; run 146 and unrelated `outputs/`
remain untouched. Run 118 is at `b402b8a5f1270545f0278926cbff7a188707fd62`.
No files were deleted, dependencies/save contracts changed, or remote writes made.

# Implemented

Local quadratic drag; separate upper bandwidth and layer heating; grey column
radiation plus actual condensation/re-evaporation/sensible heat; warm-started,
bounded climate coupling; accepted-only ocean states and real atmosphere/SST
feedback; daily finalized-wind statistics; actual MPDATA face-flux diagnostics;
reference alignment, corrected legends/cache, and optional diagnostic map IDs.
Five default maps and mandatory workbook population remain unchanged.
Details: `docs/tasks/climate-physics/FUTURE_WORK.md`.

# Verification

- `cmake --build out/build/x64-Release --config Release --parallel 4`: passed.
- `cmake --build out/build/x64-Debug --config Debug --parallel 4`: passed.
- `ctest --test-dir out/build/x64-Release -C Release --output-on-failure`: 13/13.
- `ctest --test-dir out/build/x64-Debug -C Debug --output-on-failure`: 13/13.
- `git -c core.safecrlf=false diff --check`: passed.
- Earth run 147, seed 1, 512x257: 1281.05 s summed generation stages; exactly five default maps.
- Isolated repeat: identical benchmark metrics, five map hashes and checked budget/process/ocean/statistics CSVs. All 19 selected maps and four ERA5 previews written; PNG dimensions verified.
- `out/climate-engine-verification/verify_maps.ps1`, `verify_outputs.mjs` and `verification.json`: exact commands and assertions.
- Read-only artifact-tool workbook check: run 147 is row 148, HRES 512, counts/WRE agree, no formula errors, inspection preserved file hash. Script/report: `inspect_workbook.mjs`, `workbook-inspection.json` in the same verification folder.
- Final generated-world smoke: `--generate-world --seed 73021 --resolution 32 --plate-cycles 1 --cycle-steps 1 --plates 4 --no-rivers --no-lakes --no-deltas --no-profile-log --save out/climate-engine-verification/final-smoke-32.uww`: passed, coupled convergence at iteration 3.
- A 16-wide terrain smoke stalled before reaching climate and was terminated; that terrain-scale case is unverified and outside this climate pass.
- Visual checks: run-147 LIC, particles, rainfall; optional consistency/current maps; reference LIC and workbook row. ERA5 v2 previews copied into the existing reference folder, preserving prior versions.

# Remaining uncertainty

Run 147's inner solves and snow-inclusive hydrology converged, but its outer
loop hit six passes: relative wind/heating/rainfall changes 0.216/0.231/0.318
against 0.02 tolerance. Seasonal relative water residual is at most 6e-9 in the
rounded CSV; column exchange closure below 4.47e-14 W/m2. Those budgets do not
prove equilibrium. Rainfall remains excessively banded and gyres weak.

Investigate nonlinear heating/rainfall coupling and adaptive relaxation.
Ice-covered ocean cells still need separate ice-skin/liquid-SST and freeze/melt
enthalpy treatment. Upper-layer sampling yields only about three effectively
independent states per quarter. Multi-seed, high-resolution, long-duration,
reference-variability and full runtime/save integration tiers remain unverified.
Existing MSBuild shared-intermediate-directory warnings remain.
