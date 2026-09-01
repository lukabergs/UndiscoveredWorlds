# Goal

Validate and improve the circulation and moisture transport controlling regional rainfall while preserving run 118's determinism and conservation.

# Current checkpoint

- Retain runs 139–141: solve the steady damped first-baroclinic pressure equation for the non-zonal surface mode on a 64x33 grid, preserve the zonal overturning mode, and leave the upper solver unchanged.
- Derive the thermal forcing from the 1000-to-500 hPa hypsometric geopotential anomaly. The previous hand-shaped westward thermal convolution is no longer active.
- Use 400 m equivalent depth, nondimensional damping 0.10, quadratic drag linearized at 5 m/s, and a `1e-4` physical-residual target. These are global scalars; there is no basin or latitude mask.
- Runs 127–132 record safe solver fallbacks; run 133 records the over-attenuating six-day/full-column formulation; runs 134–136 isolate equivalent depth, local forcing, and hypsometric amplitude.
- Run 137 solves zonal and stationary modes together. It creates visible gyres and cuts mean vector error to 3.682 m/s, but incorrectly damps the zonal mode: its RMS falls to 2.8–3.1 hPa versus ERA5's 7.5–8.5 hPa and ERA5 precipitation correlation falls to 0.145. Keep the tested full-mode path disabled until the two vertical modes have separate closures.
- Run 123 shows that the existing explicit column-mass pressure update is unstable at its configured timestep; the production response is restored to zero.
- Run 124 confirms that longer surface drag times restore curved circulation but over-accelerate equatorial flow in a linear steady balance.
- Runs 125 and 126 retain bulk quadratic drag only at the surface; the upper layer keeps its linear Rayleigh solver.
- Use a 300 m surface momentum layer and dimensionless drag coefficients of 0.0013 over ocean, 0.0030 over land, and 0.0060 over high relief. The ocean value and approximately 12 h effective damping at 5 m/s agree with the ERA5 pressure-work fit; land and relief remain less certain.
- Run 125 raises the surface rotation/divergence RMS ratio from approximately 0.47 to 0.76 versus ERA5's 1.44, but broad stationary pressure cells remain too weak.
- Treat rejected candidates as checkpoint-specific evidence, not proof that their physical mechanisms are intrinsically harmful.

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
- That result rejects only the tested extension at this checkpoint. It does not establish that mechanical upper forcing is intrinsically harmful; revisit it after correcting upper-layer length scale and amplitude.
- Relative to the diagnostic run 118 baseline, run 119 changes area-weighted correlation as follows: surface eastward wind 0.471 to 0.487, surface northward wind 0.285 to 0.287, pressure 0.493 to 0.495, column water 0.733 to 0.732, ascent 0.072 to 0.073, and precipitation 0.392 to 0.399.
- Maximum seasonal area-weighted water-budget residual is `2.51442e-9`; wettest 10% of land receives 41.4% of land rainfall.
- Run 123's full explicit pressure feedback is not retained: stationary-pressure RMS grows to 8.2–9.6 hPa versus ERA5's 3.1–4.7 hPa, stationary correlation falls to 0.02–0.23, and local surface divergence approaches 900 day-1. The continuity mechanism remains a candidate for implicit or stability-controlled integration.
- Run 125 improves surface `u` correlation from 0.487 to 0.519, surface `v` from 0.287 to 0.293, column water from 0.732 to 0.775, and ERA5 land-precipitation correlation from 0.396 to 0.434. Seasonal mean surface vector error increases from 3.836 to 5.759 m/s and northern 50–70 degree precipitation remains excessive.
- Surface pressure and upper height are separate fields and solvers but share the column-divergence coupling. Run 125 changes upper correlations only from `0.155/-0.072` to `0.157/-0.072`; exact isolation would require a diagnostic coupling switch.
- Runs 139–141 raise mean stationary-pressure RMS from run 125's 1.362 hPa to 2.891 hPa versus ERA5's seasonal 3.1–4.7 hPa. Mean direction error improves from 73.34 to 70.89 degrees, while mean vector error is effectively unchanged at 5.790 m/s.
- Relative to run 125, runs 139–141 change surface `u/v` correlation from `0.519/0.293` to `0.508/0.306`, column-water correlation from `0.775` to `0.773`, and ERA5 precipitation correlation from `0.434` to `0.410`. This checkpoint is retained for its physically derived pressure amplitude and inspectable gyre response, not as a rainfall-score optimum.

# Verification

- Debug and Release builds pass.
- All seven Debug and Release test suites pass.
- Two 512x257 seed-1 Earth runs produced byte-identical sets of 85 validation artifacts: aggregate SHA-256 `635A1BFABE433FE199625AB37303262E566E66E11B1E222D18324E13A8FEB065`.
- A second post-fix 512x257 seed-1 Earth run reproduced all eight LIC PNGs byte-for-byte.
- The 28 circulation PNGs in `extra/img/earth/benchmark` match their run-120 diagnostic sources by SHA-256; run 120 contains 40 raw GeoTIFFs and the CSV manifest.
- Runs 121 and 122 produced byte-identical sets of 36 circulation PNGs. All 16 cached ERA5 PNG hashes were unchanged during run 122; its cache status reports reuse.
- Cached run 122 completed in 26.878 seconds at 512x257. The LIC contribution remains approximately 5.3 seconds relative to the pre-LIC diagnostic stage.
- Runs 125 and 126 match in 103 of 105 recorded artifacts; the only differences are run metadata. All 38 published benchmark PNGs are byte-identical.
- Maximum run-125 seasonal area-weighted water-budget residual is `2.75943e-9`; wettest 10% of land receives 50.2% of land rainfall.
- Runs 139–141 match in 103 of 105 recorded artifacts; only run metadata differs. All 38 published PNGs are byte-identical, and maximum seasonal water-budget residual is `5.09e-10`.
- Recorded benchmarks: `extra/validation/runs/119` through `extra/validation/runs/141`.

# Remaining uncertainty

- Generalization beyond the Earth fixture remains unmeasured.
- Upper-level wind correlations remain poor (`u=0.157`, `v=-0.072`) and need layer-specific length-scale and amplitude work.
- The stationary solve improves pressure-cell amplitude but the zonal bands still dominate. The next circulation phase should represent zonal and stationary/baroclinic modes separately and use diagnosed diabatic heating instead of surface-temperature anomaly as the wave source.
- Whether a true positive-definite MPDATA implementation improves spatial rainfall without harming conservation remains untested.
