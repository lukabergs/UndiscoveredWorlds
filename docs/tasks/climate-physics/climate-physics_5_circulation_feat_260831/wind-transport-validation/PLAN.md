# Goal

Validate and improve the circulation and moisture transport controlling regional rainfall while preserving run 118's determinism and conservation.

# Final state

- Run 118 remains reproducible at commit `b402b8a5f1270545f0278926cbff7a188707fd62`.
- Run 119 records the retained nonlocal thermal-pressure response, surface mechanical topographic forcing, circulation exports, and particle-trail maps.
- The CFL-controlled transport candidate was tested but not activated: it reduced column-water correlation from 0.733 to 0.461 and precipitation correlation from 0.392 to 0.313.
- The upper-height mechanical forcing candidate was rejected because it did not improve upper-level wind agreement.
- Relative to the diagnostic run 118 baseline, run 119 changes area-weighted correlation as follows: surface eastward wind 0.471 to 0.487, surface northward wind 0.285 to 0.287, pressure 0.493 to 0.495, column water 0.733 to 0.732, ascent 0.072 to 0.073, and precipitation 0.392 to 0.399.
- Maximum seasonal area-weighted water-budget residual is `2.51442e-9`; wettest 10% of land receives 41.4% of land rainfall.

# Verification

- Debug and Release builds pass.
- All seven Debug and Release test suites pass.
- Two 512x257 seed-1 Earth runs produced byte-identical sets of 85 validation artifacts: aggregate SHA-256 `635A1BFABE433FE199625AB37303262E566E66E11B1E222D18324E13A8FEB065`.
- Recorded benchmark: `extra/validation/runs/119`.

# Remaining uncertainty

- Generalization beyond the Earth fixture remains unmeasured.
- Upper-level wind correlations remain poor (`u=0.155`, `v=-0.072`) and need a separate circulation-model phase.
