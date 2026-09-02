# Climate roadmap: implemented features and validation queue

Implementation branch: `climate-physics_7_core_feat_260902/future-work`.
Climate-engine mechanisms and their production integrations are implemented;
this is not a claim that every reduced closure is physically complete or that
Earth skill improved. The user requested implementation before tuning and
explicitly deferred gameplay weather. The branch PLAN records verification.

## Numerical foundations

- `climate_grid.*`: shared spherical centres, areas, face lengths and closed
  polar faces. Internal grids are even `W x W/2`; arbitrary pole-inclusive
  input/output rasters use conservative overlap remapping, including `512x257`
  and `2048x1025`. Energy, circulation, hydrology, ocean and validation use this
  contract. Output-only diagnostics retain their explicit pole-inclusive layout.
- `climate_hydrology.*`: donor-cell finite-volume MPDATA with multidimensional
  corrective fluxes, divergent-flow correction, optional monotonicity and local
  CFL substeps. Optional endpoint winds are interpolated at each substep midpoint.
  The former 94% transfer cap and multi-cell centroid jumps are removed.
  Corrective-advector formulation: [Jaruga et al., section 3.1](https://doi.org/10.5194/gmd-8-1005-2015).
- `climate_flow.*`: adaptive midpoint RK2 particle trails, curvature/displacement
  subdivision and explicit elapsed time. Trails are diagnostic, never moisture
  transport; boundary/substep exhaustion is reported rather than slowing flow.

Tests: grid constants/integrals/poles/seam, odd and large raster remaps; tracer
rotation/refinement/diffusion, reversible shear, convergence/divergence,
changing winds, positivity and mass; uniform flow, vortex and polar trails.

## Diagnosed and mode-separated circulation

- `climate_atmosphere.*`: conservative two-layer grey shortwave/longwave
  exchange, bulk sensible heat, actual layer condensation and falling-rain
  re-evaporation diagnose atmospheric heating in W/m2. Surface evaporation
  enters the surface/vapour latent-energy budget, not a duplicate atmospheric
  heating term. Column exchange closure is tested and exported. Optical depths
  and vertical condensation placement remain reduced-model parameters. Framework:
  [GFDL idealized moist model](https://www.gfdl.noaa.gov/idealized-moist-spectral-atmospheric-model-quickstart/).
  Stationary projection removes row means; observations remain validation-only.
- The zonal closure uses thermal contrast and Held-Hou overturning; its upper
  background is hydrostatic thickness between 1000 and 500 hPa. Non-zonal
  surface pressure and upper geopotential have independent equivalent depths,
  restoring times, spectral cutoffs and solvers. Surface quadratic drag uses
  the same local ocean/land/relief coefficients in the pressure solve and final
  wind diagnosis; the upper layer uses linear damping. The old output-grid
  pressure/height relaxation is removed.
- Surface thermal anomalies now add a hydrostatic pressure response, tapering
  to zero at 700 hPa; heating still supplies its separate diabatic response.
  This reduced boundary-layer closure follows the mechanism of
  [Lindzen and Nigam (1987)](https://ntrs.nasa.gov/citations/19880029939).
  The calibrated overturning depth is 16 km and the subpolar low sits 65% of
  the remaining latitude span poleward of the Hadley edge. These are reduced
  closures, not resolved eddy dynamics or observed pressure inputs.
- Upper topography uses damped rotating mountain-wave propagation/evanescence
  of zonal terrain modes, not copied surface windward/lee pressure. Calm
  incident flow gives zero wave forcing. This remains a reduced vertical-mode
  approximation, not resolved three-dimensional topographic flow.
- Simulated dry/moist lapse rate, temperature and gravity diagnose stability;
  mode depth sets wave speed/equivalent depth. Radius, rotation and forcing scale
  constrain adjustment length, damping and spectral bandwidth.
- Weak stationary momentum exchange is equal/opposite and preserves zonal row
  means. Finite-volume layer divergence supplies interface ascent. Energy,
  drag, mass, momentum-exchange and row-transfer diagnostics are exposed.
- GMRES records restart residual histories and stagnation, accepts physical
  residuals and supports configurable restart windows. Each failed stationary
  layer falls back independently to its zonal background. ILU/multigrid remain
  conditional optimization choices, not required dependencies: deterministic
  fixtures converge with the existing preconditioner.

Production: `climate_ocean.cpp`, `globalclimate.cpp`, `generation_workbench.cpp`.
Bounded atmosphere/ocean/hydrology iterations feed diagnosed layer heating back
into circulation without repeating coastal temperature postprocessing. They
retain hydrological stores between passes, under-relax heating, and require
wind, SST, applied-forcing/heating and rainfall changes plus accepted inner
solves before reporting convergence. Snow mass participates in spin-up acceptance;
caps remain explicit nonconvergence. The direct latent-heating ascent feedback is disabled
when that circulation projection is active; thermodynamic moist adjustment stays.
Isolation switches in `generation_tuning.hpp`: zonal, stationary, surface, upper,
lagged coupling and prognostic atmosphere.

Tests: all 16 isolation combinations, zero forcing, hydrostatic units, rotation
reversal/equatorial regularity, heating budgets, stability/bandwidth scaling,
mountain-wave response, residual acceptance and failed-solve fallback.

## Climate sampling; gameplay deferred

`climate_weather.*` implements forced/damped, linear one/two-layer spherical
shallow-water vertical modes. Equivalent depth sets wave speed; signed modal
geopotential is not a positive total water-column thickness. RK2, local wave
viscosity, CFL/source substeps, amplitude bounds and atomic rollback protect the
reduced model. The climate path spins up then averages 30 evolving states; its
blend and independent fallback are explicit controls.
The evolving solve receives the original equilibrium forcing, avoiding a
second attenuation of the already adjusted stationary response. A matching
forced-mode regression covers steady/evolving pressure consistency.

An existing experimental gameplay backend follows interpolated seasonal jets.
`planet::advanceweather(seconds, optionalHorizontalCells)` advances only the
anomaly state. It preserves climate maps, supports conservative LOD changes,
and serializes time, all layer fields and RNG state for replay. World save
version is **17**; version-16 files are not migrated. Climate generation no longer
initializes the runtime gameplay state. Runtime storms, fronts and gameplay
integration are outside the current work, even where reusable backend code exists.

Hydrology continues the climate anomaly state between months and retains at least
90/91/92/92 daily intervals per Jan-Mar/Apr-Jun/Jul-Sep/Oct-Dec season in the final
annual cycle. Samples contain the finalized transport winds for both layers,
with actual elapsed-time weights. Statistics include mean vectors, mean speed,
directional consistency, speed spread, independent-sample error and an estimated
positive-sequence autocorrelation correction. Daily samples are not automatically
independent; the latter is not a validated confidence interval or storm climatology.
Prescribed sinusoidal synoptic forcing is disabled; the explicitly parameterized
coastal day/night cycle remains switchable. Exports include `weather_statistics.csv`.

Tests: forced 64x32 spin-up plus 30 samples, zero state, mass/energy bounds,
jet advection, weighted statistics, malformed-state rejection, LOD and exact
serialized replay. ERA5 comparisons must use matching season/sample definitions.

## Wind-driven ocean and SST

`climate_ocean_dynamics.*` replaces prescribed currents in production with a
coast-bounded Stommel/PV streamfunction using wind stress, Coriolis, drag and
bathymetry. Corner streamfunction is m3/s; closed face transports are exactly
non-divergent, including islands and variable depth. Stored raster currents
are cm/s. Equator-regularized Ekman transport supplies upwelling/downwelling.

Mixed-layer advection/diffusion uses equal/opposite face heat fluxes and explicit
vertical reservoir exchange to close fixed-depth continuity. Surface/deep heat
exchange is included in the budget, with actual diagnosed surface heat flux
when available. Production calls the ocean's one-way numerical solve and returns
SST through the real atmospheric solve in the shared outer iteration. The old
empirical SST-to-pressure/wind surrogate is not used in production. Both basin
and outer fixed-point residuals are reported. Finite-but-unconverged ocean states
are rejected; fallback uses atmospheric SST and zero currents/upwelling.
`oneWayDiagnosticsOnly` prevents diagnosed SST from feeding atmospheric
temperature/evaporation while keeping ocean diagnostics available. Wind stress
shares the atmosphere's ocean drag, air density, radius and rotation constants.
Liquid SST and the deep reservoir cannot cool below the prescribed seawater
freezing point (-1.8 C). A stationary slab-ice enthalpy reservoir receives further
cooling; warming melts it before heating liquid. Zero-capacity conductive ice
skin is diagnosed separately with linearized atmospheric exchange. Production
feeds skin temperature and seasonal ice cover to atmospheric heating/evaporation
while SST/current diagnostics retain liquid temperature. Ocean CSVs include ice
thickness and skin temperature; heat closure includes latent energy. Initial
seasonal ice thickness is prescribed as 1 m in frozen base-climate cells.
This is a fixed seasonal-interval response, not a spun-up three-dimensional
global ocean/sea-ice model; ice transport, brine and snow-on-ice are unresolved.
Atmospheric latent/cloud temperature adjustments are separate from surface
temperature. Surface exchange uses skin saturation and actual vapour storage;
lifted-parcel cooling remains in condensation, not surface sensible exchange.
Albedo and transfer coefficients vary continuously with snow/ice cover.

Both generation pipelines now order winds before ocean/SST before evaporation
and rainfall. Tests cover this order, coastal faces, western intensification,
open polar oceans, variable-depth continuity, uniform SST, heat closure,
deterministic coupling and rejected incomplete basin solves.

## Diagnostics and workflow

- LIC and adaptive RK2 trails remain selectable at run resolution. ERA5 uses the
  same renderer, with versioned reusable caches; its source vectors are unchanged.
- Actual donor-plus-corrector face fluxes are accumulated by layer. Their
  divergence reproduces MPDATA cell changes; exports include mean two-layer
  moisture transport, convergence, ascent, column water, radiation/latent heat,
  ocean currents/SST and statistical wind consistency. Column convergence is
  not total moisture multiplied by the surface wind. These diagnostics describe
  transport separately from any configured post-transport smoothing.
- `--resolution` remains configurable; default outputs remain the five requested
  maps. Additional diagnostic maps are optional `--map` IDs (see README).
  Numerical diagnostics persist per run even when images are disabled.
- Mandatory workbook population, short map names and reference organization are
  retained. Reference categorical masks are aligned to output resolution before
  scoring; northern precipitation sampling uses the reference-grid coordinates.
  Console WRE now uses the workbook's reference counts and normalization.
- Legends give the annual precipitation scale for both PNG and float TIFF,
  distinguish representative-month LICs from quarter-averaged statistics, and
  signed fluxes from speed. Diagnostic budgets are not claims of globally
  equilibrated climate.

## Remaining tuning/validation, not missing engine activation

- Coarse controls 151/158 (seed 1, 64x33, 12 passes): annual ERA5 surface u/v
  correlations improve 0.525/0.338 to 0.735/0.406; quarterly vector RMSE falls
  4.730 to 4.011 m/s. Stationary vector correlation improves only 0.172 to 0.221;
  tropical basin structure and upper stationary circulation remain poor.
  Equatorial rainfall falls 2181 to 2101 mm/year (IMERG 1808), land RMSE 806 to
  790, but mean land rainfall worsens 469 to 404 against 708 mm/year. Better
  wind-belt placement has not solved moisture delivery and rainfall geography.
  Heating remains 2.77%, above the unchanged 2% fixed-point tolerance.
  The earlier 512x257 attempt was interrupted. Fresh workbook baseline 159
  subsequently completed at 64x33, seed 1, selected physics/24-pass cap;
  coupling passed at iteration 13 with 1.49% aggregate heating change.
  Higher-resolution validation remains open.
- Fresh validation workbook: `outputs/01a062ac-02ef-7710-bf6a-25bd1f12cf87/climate_validation.xlsx`.
  Only run 159 is included. Numeric cell exports and
  `scripts/summarize-climate-validation.py` provide paired spatial inputs;
  Excel calculates derived scores. Global Koppen uses spherical areas and
  common masks; regional tables cover other variables. Reference gaps are explicit.
- Export the ocean solver's native wet-cell mask: final terrain changes after
  climate, so it cannot reliably reconstruct that mask. Workbook ocean summaries
  use cells identified by nonzero captured currents or ice and flag partial
  coverage; entirely stationary ice-free cells may be omitted.
- Slab-ice freeze/melt and skin/SST separation are implemented. Validate the
  prescribed initial ice reservoir and linearized ice-surface exchange against
  seasonal ice extent; no observational sea-ice skill is claimed yet.
- Earth moisture/precipitation/Koppen, surface/upper ERA5 vectors, pressure
  decomposition, jets, shear, ascent, ocean currents and SST comparisons.
- Prioritize seasonal basin pressure structure and amplitude next: the thermal
  pressure correction helps winter much more than spring/summer. Both global
  and extratropical-ocean increases in momentum depth failed controlled trials
  (148, 157). Diagnose simulated surface-temperature/heating gradients before
  further drag tuning; then assess moisture delivery and ascent over dry land.
- Multi-seed/planet-parameter and long-duration climate/energy soak suites;
  uncertainty-estimator validation and recovery of parent seasonal climatology.
- Benchmark runtime, forcing-bandwidth/grid-refinement studies and tuning of
  spin-up/iteration budgets. Some coupled worlds may hit the iteration cap.
- If polar CFL subcycling dominates runtime, evaluate a conservative reduced
  longitude grid or flux-form semi-Lagrangian transport; do not recover speed
  by restoring multi-cell centroid jumps or weakening positivity checks.
- Controlled moisture-smoothing, transport-order, mode-isolation and parameter
  ablations; resolution studies. Higher-order advection does not remove diffusion
  introduced by other operators or compensate for inaccurate circulation.
- Visual map/LIC/particle agreement and high-resolution world generation.
- Further observational validation (including optional CFSR cross-checks and
  actual sub-monthly ERA5 variability) needs suitable reference fields; monthly
  mean vectors alone cannot validate directional consistency.
- Full application save/load/runtime-weather round trip across supported hosts;
  binary weather-state replay is covered by deterministic unit tests.

Commands and completed verification are recorded in the branch `PLAN.md`.
