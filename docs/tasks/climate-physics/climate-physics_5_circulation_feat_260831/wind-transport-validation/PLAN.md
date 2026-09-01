# Goal

Validate and improve the circulation and moisture transport controlling regional rainfall while preserving run 118's determinism and conservation.

# Final state

- Run 118 remains reproducible at commit `b402b8a5f1270545f0278926cbff7a188707fd62`.
- Run 119 records the retained nonlocal thermal-pressure response, surface mechanical topographic forcing, circulation exports, and particle-trail maps.
- Run 120 is a visualization-only follow-up to run 119. It adds deterministic LIC surface/upper wind maps, publishes all 28 seasonal circulation previews beside the benchmark maps, and includes a quantity/palette guide; all recorded climate metrics are identical to run 119.
- Runs 121 and 122 add matched ERA5 LIC/particle references and seasonal surface/upper vector-error maps; all climate metrics remain identical to run 120.
- Eligible runs now export both simulated LIC and particle trails. Flow rendering is limited to 1024 horizontal and 525312 total cells; larger runs retain scalar and raw diagnostics.
- The 16 ERA5 maps are cached by renderer version and grid dimensions. Run 122 reused them byte-for-byte while regenerating the run-specific vector errors.
- The 512x257 seasonal mean vector error is 3.836 m/s at the surface and 16.499 m/s aloft. The comparison shows overly meridional, zonally repetitive surface flow and much larger upper-level structural errors.
- The CFL-controlled transport candidate was tested but not activated: it reduced column-water correlation from 0.733 to 0.461 and precipitation correlation from 0.392 to 0.313.
- MPDATA itself was not implemented. Production moisture transport remains the conservative adjacent-row scheme retained for run 118; replacing it needs a separate numerical-transport phase.
- The upper-height mechanical forcing candidate was rejected because it did not improve upper-level wind agreement.
- Relative to the diagnostic run 118 baseline, run 119 changes area-weighted correlation as follows: surface eastward wind 0.471 to 0.487, surface northward wind 0.285 to 0.287, pressure 0.493 to 0.495, column water 0.733 to 0.732, ascent 0.072 to 0.073, and precipitation 0.392 to 0.399.
- Maximum seasonal area-weighted water-budget residual is `2.51442e-9`; wettest 10% of land receives 41.4% of land rainfall.

# Verification

- Debug and Release builds pass.
- All seven Debug and Release test suites pass.
- Two 512x257 seed-1 Earth runs produced byte-identical sets of 85 validation artifacts: aggregate SHA-256 `635A1BFABE433FE199625AB37303262E566E66E11B1E222D18324E13A8FEB065`.
- A second post-fix 512x257 seed-1 Earth run reproduced all eight LIC PNGs byte-for-byte.
- The 28 circulation PNGs in `extra/img/earth/benchmark` match their run-120 diagnostic sources by SHA-256; run 120 contains 40 raw GeoTIFFs and the CSV manifest.
- Runs 121 and 122 produced byte-identical sets of 36 circulation PNGs. All 16 cached ERA5 PNG hashes were unchanged during run 122; its cache status reports reuse.
- Cached run 122 completed in 26.878 seconds at 512x257. The LIC contribution remains approximately 5.3 seconds relative to the pre-LIC diagnostic stage.
- Recorded benchmarks: `extra/validation/runs/119` through `extra/validation/runs/122`.

# Remaining uncertainty

- Generalization beyond the Earth fixture remains unmeasured.
- Upper-level wind correlations remain poor (`u=0.155`, `v=-0.072`) and need a separate circulation-model phase.
- Whether a true positive-definite MPDATA implementation improves spatial rainfall without harming conservation remains untested.
