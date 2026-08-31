# Goal

Replace run 43's threshold-dominated, spatially sparse precipitation with a conservative continuous closure that improves land rainfall and Köppen geography.

# Phase

Climate physics phase 4: clouds, runoff, and hydrological storage.

# Scope

- Instrument raw precipitation, area-weighted water budgets, and distribution diagnostics.
- Isolate condensation-threshold and soil-initialization effects in Release benchmarks.
- Retain supported precipitation and land-surface changes.
- Follow remaining evidence into circulation or seasonal-cycle corrections when required.
- Measure temperature-locked Köppen area, coarse precipitation skill, and rainfall-storage clipping.
- Export per-run annual-mean temperature previews and unclamped float32 annual-precipitation GeoTIFFs.
- Test land thermal inertia and eddy moisture transport as an interaction rather than isolated effects.
- Acquire independent regional masks and monthly Earth climatologies for spatial and process-level validation.
- Establish an observed-input Köppen control before continuing rainfall-model tuning.

# Constraints

- Earth data remains validation-only.
- Preserve deterministic fixed-seed generation and explicit water conservation.
- Do not change Köppen thresholds to hide upstream errors.
- Benchmark retained physics changes in Release and record their run IDs.

# References

- `climate_ocean.cpp`
- `climate_hydrology.cpp`
- `climate_physics.cpp`
- `generation_tuning.hpp`
- `climate_validation.cpp`
- `extra/validation/seed_1/`
- `extra/climate.xlsx`
- `climate_benchmark_runs.json`
- `extra/reference/ipcc-ar6-wgi-regions-v4/`
- `extra/reference/clara-a3-albedo-1991-2020/`
- `extra/reference/nsidc-g02202-v6-2001-2020/`
- `extra/reference/era5-planette-2001-2020/`

# Steps

1. Add diagnostic metrics that distinguish thresholding, rounding, and land recycling.
2. Run isolated critical-humidity and wet-start experiments.
3. Implement a 365-day prognostic hydrology year over twelve interpolated monthly forcings.
4. Partition precipitation conservatively into stratiform, orographic, and convergence-fed convective components.
5. Remove the geographic condensation threshold and use surface pressure in moisture thermodynamics.
6. Converge atmospheric and soil storage separately and persist process-level diagnostics.
7. Add thermal-reachability, 5-degree precipitation, and storage-clipping diagnostics.
8. Benchmark the completed hydrology architecture against IMERG, ERA5, and Köppen only after Release verification passes.
9. Correct circulation or closure parameters from process diagnostics, one mechanism per run.
10. Verify and commit the retained state.

# Done when

- Zero-rain land area is materially reduced without breaking water conservation.
- Land precipitation spatial correlation and area-weighted Köppen agreement improve over run 43.
- The implementation remains deterministic and generative-world compatible.

# Status

- Run 74 reproduced run 71 and exposed weak zonal pressure structure plus excessive temperature-locked polar climate.
- Runs 75-77 isolated local-ice coupling behavior; annual-temperature coupling snowballed Siberia, while warm-season-only coupling destroyed Antarctica.
- Run 78 retained a converged warm-season ablation plus elevation-seeded ice-sheet proxy: area-weighted kappa 0.1956, WRE 0.5933, northern 50-70 warmest mean 12.67 C, and southern 60-90 EF thermal eligibility 0.9986.
- Antarctica now remains EF without a latitude-clamped band, but eastern Siberia and the Russian Far East remain too polar and Dfc remains deficient.
- IPCC diagnostics now separate ET from EF so mainland-Arctic errors cannot hide inside aggregate E climate.
- Runs 79-81 retained temperature-derived pressure belts and finite-time condensation. Run 79 remains the spatial peak (kappa 0.2035); run 81 improved WRE to 0.5176 and northern 50-70 precipitation ratio to 0.5752 but retained unskilled stationary pressure and wind.
- Run 82 preserved pressure and upper height in floating point through wind generation. Input pressure-gradient noise fell materially, but the circulation solve regenerated stronger unresolved gradients: surface-wind RMS stayed near 10.2 m/s, column-divergence RMS near 52/day, transport-wind correlations near 0.106/-0.010, and kappa reached only 0.1976.
- Run 82 climate-boundary density improved only 0.4% from run 81 and remained 60% above the reference; Sahara-Arabia-Iran BWh fell from 25.0% to 24.2% against 70.9% reference.
- Run 83 projected pressure and upper height to zonal wavenumbers 1-4. Circulation diagnostics improved materially, but independent row projection created horizontal climate striping, increased precipitation concentration, and reduced northern rainfall.
- Run 84 added coherent meridional modes. Pressure and meridional-wind skill improved again, but the four-wave cutoff remained too latitude-dominated and degraded coarse precipitation placement.
- Run 85 retained zonal wavenumbers 1-8, but worsened kappa to 0.1939 and coarse precipitation correlation to 0.4239. Runs 83-85 all show unacceptable horizontal striping.
- Run 86 projected the previously unfiltered zonal mean and removed the narrow striping, confirming its cause, but hard full-state projection over-concentrated rainfall (83.7% in the wettest land decile) and reduced northern rainfall to 15.4% of reference.
- Run 87 left the evolving state intact and filtered only wind-driving gradients. It removed horizontal striping and improved land precipitation correlation to 0.369 and vertical-ascent correlation to 0.059, but kappa fell to 0.1933 and rainfall remained too concentrated.
- Run 88 reproduced run 87 exactly and verified fixed-scale annual-mean temperature and annual precipitation images. The physical rainfall map exposes global comb-like condensation bands that are obscured by the Köppen map.
- The current experiment retains run 87's detailed horizontal transport winds but derives circulation adjustment and vertical ascent from the fully resolved large-scale pressure and upper-height fields.
- Simulator runs are paused after run 89 while an observed-input control isolates the classifier/import pipeline from the atmosphere and hydrology.
- The twelve-month observed-input control reaches area-weighted spatial kappa 0.7685 with WorldClim temperature/precipitation and 0.7013 with WorldClim temperature plus IMERG precipitation. This confirms that the dominant remaining discrepancy is upstream climate-field generation, not the atlas comparison pipeline.
- Benchmark precipitation output now uses the same single-band IEEE float32 raster model as IMERG, with physical mm/year values and WGS84 GeoTIFF metadata; the false-colour PNG is preview-only.
- Run 89's vertical-motion RMS is 98.5-103.5 hPa/day versus ERA5's 30.4-37.4, with 49.3-58.6% of cells limiter-bound; the current experiment corrects that amplitude and its nonzero global mean before changing orographic physics.
- Run 90 reduced vertical-motion RMS to 62.6-68.0 hPa/day, limiter use to 12.4-14.9%, and mean motion below 0.21 hPa/day, but changed precipitation negligibly. The repeated bands therefore come from horizontal moisture convergence, not vertical parcel cooling.
- The current experiment disables the unstable divergent surface-pressure feedback that amplified base surface winds from about 2.1 to 16-17 m/s and produced a 26.8-30.2/day column-divergence field.
- Run 91 retained that removal: pressure correlation improved from 0.169 to 0.481, surface-wind correlations from 0.007/0.069 to 0.460/0.301, and the spurious subtropical rain belts weakened sharply. It under-rains globally and still reverses important South American wet/dry geography.
- The current experiment replaces the dimensionless ridge barrier with parcel displacement from the resolved wind-terrain slope, using moist ascent and dry descent lapse rates to trigger windward condensation and lee warming.
- Run 92 changed only 5,815 precipitation pixels with a 0.07 mm/year global mean absolute difference; one-cell parcel displacement cannot represent the accumulated climb across a resolved mountain range.
- The current experiment accumulates only slope-signed elevation rise/drop along a six-cell wind trajectory and records orographic displacement statistics in the atmospheric diagnostic.
- Run 93 diagnosed maximum terrain-following parcel climbs of 1.1-1.5 km but affected too little rainfall to correct regional geography. Moisture-transport wind skill (0.184/0.024) remains much worse than surface-wind skill (0.460/0.301).
- Run 94 removed the 500 hPa blend. Transport-wind and land-precipitation correlations improved to 0.518/0.186 and 0.432, but 48.0% of land became exactly dry and rainfall collapsed into axisymmetric ocean bands.
- Run 95 showed that removing the ocean-to-continent condensation-threshold transition worsens the dry-continent state and does not remove the bands; the prior threshold was restored.
- Run 96 added row-mass-neutral thermal surface-pressure anomalies while leaving the unstable divergence feedback disabled. Stationary pressure now has 1.28-1.70 hPa RMS with 0.28-0.54 correlation against ERA5, overall pressure correlation rose to 0.492, and area-weighted Köppen kappa reached 0.2364. The stationary amplitude remains below ERA5's 3.08-4.67 hPa and 53.7% of land still has exactly zero raw precipitation.
- A run-96 transect confirms essentially zero simulated rain across the Amazon near 10 S versus roughly 1,700-2,000 mm/year in IMERG, with excess rain offshore west of the Andes.
- The current experiment raises only the thermal stationary-pressure response from 0.12 to 0.30, targeting ERA5's 3.08-4.67 hPa stationary RMS and testing whether stronger summer continental lows draw Atlantic moisture into South America.
- Run 97 reached the ERA5 stationary-pressure amplitude but reduced pressure/wind/rainfall skill, raised raw-zero land to 60.3%, and left the South American monsoon region effectively dry. The 0.30 response was rejected and 0.12 restored.
- The current experiment replaces both surface-only transport and the arbitrary surface/500 hPa blend with an 850 hPa pressure-level wind derived from sea-level pressure and hypsometric thickness over a 150 km stencil. Orographic displacement now uses this same moisture-bearing wind.
- Run 98 reduced raw-zero land to 39.6% and restored northern rainfall totals, but degraded land precipitation correlation to 0.269, coarse correlation to 0.378, transport-v correlation to 0.065, and left the South American monsoon nearly dry. The diagnostic 850 hPa transport was rejected.
- Run 99 restored run 96's surface-directed transport and tested only full-capacity initial soil moisture. After eight cycles the hydrology was still not converged (1.72% storage change/cycle), while only 5.8% of annual-precipitation pixels changed from run 96 with a 3 mm/year mean absolute change.
- Wet initialization reduced area-weighted raw-zero land from 53.7% to 43.7% but did not improve placement: land correlation fell from 0.4168 to 0.4145, 5-degree correlation from 0.5254 to 0.5229, kappa from 0.2364 to 0.2303, and the South American monsoon remained near zero. The experiment was rejected and dry initialization restored.
- Rainfall work is paused after run 99. The retained state is run 96's circulation/hydrology configuration; Köppen remains downstream of an overly zonal, weakly recycled precipitation field.
- The replacement hydrology is now staged but has not been benchmarked: it advances a full 365-day periodic year using twelve monthly forcings interpolated from the four serialized climate snapshots, then stores three-month precipitation means back into the unchanged save format.
- Precipitation is now a conservative sum of uniform-RH stratiform condensation, terrain-attributed condensation, and warm Kuo-style convection limited by smoothed moisture-flux convergence plus surface evaporation. The former land/ocean condensation-threshold interpolation is removed.
- Moisture thermodynamics now use the generated surface-pressure anomaly; soil spin-up begins at field capacity, evapotranspiration stress is referenced to readily available soil water, and atmospheric/soil convergence must each pass tolerance.
- Run 100 benchmarked the full-year closure. The atmosphere stabilized quickly, but soil storage still changed 2.98% in cycle 8 and did not converge.
- Run 100 conserves the area-weighted water budget to numerical precision, and precipitation components close against total precipitation. However, orographic precipitation contributes only 0.009% of the total; stratiform and convective precipitation contribute about 43.0% and 57.0%.
- Run 100 worsened the principal rainfall targets versus runs 96/99: area-weighted land correlation 0.3359, 5-degree correlation 0.4737, raw-zero land 59.3%, wettest-decile share 94.0%, land E/ocean E 1.28%, runoff/land P 88.7%, recycling 1.13x, and area-weighted Koppen kappa 0.2286.
- Visual inspection confirms a narrow, horizontally striped equatorial-ocean rain belt and nearly dry continental interiors. The replacement closure is conservative but not yet a physical solution; its terrestrial recycling, orographic coupling, and convergence geography need correction.
- Benchmark recording now exports images, diagnostics, per-run raw-pixel CSV, and JSON before attempting Excel. Excel failure is non-fatal and leaves a diagnostic error file for later backfill.

# Verification

- `cmake --build --preset build-x64-release`
- `ctest --test-dir out/build/x64-Release -C Release --output-on-failure` (7/7 passed, including hydrology calendar, component conservation, float32 GeoTIFF structure, and sample preservation)
- Fixed-seed Release climate benchmarks 71-100 with logged metrics, images, and persistent diagnostic bundles.
- Visual and class-boundary comparison of runs 71-100 against `extra/img/earth/benchmark/0.png`.
- `python -m py_compile scripts/prepare-imerg-reference.py scripts/benchmark-observed-koppen.py`
- Pillow inspection of the annual IMERG reference and reconstructed run-89 GeoTIFF: mode `F`, 2048x1025, IEEE float32, EPSG:4326.
- Direct observed-input real-fixture classifications at `extra/validation/observed-climate/`.

# Open questions

- How much of run 81's 10.2 m/s surface-wind RMS, 51/day column-divergence RMS, and near-zero meridional-wind skill comes from the whole-hPa/whole-metre intermediate storage round trip?
- After removing that round trip, what stationary pressure scales remain unresolved and require pressure-space regularization?
- Can continental evaporation/recycling increase without amplifying the existing localized precipitation spikes?
- What spin-up convergence criterion and maximum duration are required before soil-state experiments can be compared fairly?
- Does run 96 put tropical rainfall on the windward side of the Andes and southern African escarpment, or is a remaining wind/advection orientation error present?
- Can a generated 850 hPa transport layer preserve run 94's directional skill without its dry interiors and zonal ocean bands?
- Why does the full-year closure drain soil for at least eight years while atmospheric storage stabilizes, and what land-surface flux or storage timescale is missing?
- Why is resolved orographic precipitation only 0.009% of total precipitation despite visually strong coastal mountain maxima?
- Can terrestrial evaporation and retention be corrected without masking the still-axisymmetric moisture-convergence field?
- How should NSIDC sea ice, MODIS snow cover, and CLARA-A3 albedo be aligned and added as validation-only diagnostics?
- How much of the IMERG control's lower boreal score is caused by the 2001-2022 precipitation versus 1970-2000 temperature period mismatch, rather than precipitation-product differences?
