# Goal

Validate and improve the circulation and moisture transport that control run 118's regional rainfall field without losing determinism or conservation.

# Phase

Climate physics phase 5: calibration, regression hardening, and cleanup.

Paused after the first final Earth verification run; repeat-run hashing and benchmark recording remain.

# Scope

- Export seasonal circulation and moisture-flux diagnostics.
- Render particle-trail wind maps with explicit direction glyphs.
- Make Earth wind comparisons populated and like-for-like.
- Characterize transport accuracy with deterministic spherical tracer tests.
- Retain transport and pressure-forcing changes only when diagnostics improve.

# Constraints

- Preserve run 118 at commit `b402b8a5f1270545f0278926cbff7a188707fd62`.
- Preserve fixed-seed determinism, positivity, and area-weighted water conservation.
- Treat Earth reference data as validation-only.
- Keep the pre-existing untracked `outputs/` directory untouched.

# References

- `extra/validation/runs/118/run_manifest.txt`
- `climate_atmosphere.cpp`
- `climate_ocean.cpp`
- `climate_hydrology.cpp`
- `climate_validation.cpp`
- `generation_tuning.hpp`

# Steps

1. Completed: semantic circulation exports and particle-trail rendering.
2. Completed: resampled, like-for-like wind-reference comparisons and the run 118 diagnostic baseline.
3. Completed: deterministic tracer tests cover conservation, positivity, centroid motion, spreading, poles, and Courant numbers above one.
4. Completed: rejected production CFL-controlled transport because it reduced area-weighted column-water correlation from 0.733 to 0.461 and precipitation correlation from 0.392 to 0.313.
5. Completed: retained the nonlocal thermal-pressure response and surface mechanical topographic forcing; rejected the ineffective upper-height extension.
6. In progress: Debug/Release tests and one final Earth run pass; repeat-run hashing and benchmark recording remain.

# Done when

- Seasonal wind structure and moisture-flux convergence are directly inspectable.
- Wind validation reports measured values rather than missing-value zeroes.
- Transport tests cover mass, positivity, centroid motion, spreading, poles, and Courant numbers above one.
- Retained physics changes improve spatial diagnostics without regressing run 118's invariants.

# Verification

- `cmake --build out/build/x64-Debug --config Debug`
- `ctest --test-dir out/build/x64-Debug -C Debug --output-on-failure`
- Release benchmark and repeat run with seed 1 at 512 px.
- Compare wind, precipitation, determinism, timing, and seasonal water budgets with run 118.

# Open questions

- Generalization of the retained dimensionless circulation response beyond the Earth fixture remains unmeasured.
- Upper-level wind agreement remains poor and was not improved by the tested mechanical-height forcing.
