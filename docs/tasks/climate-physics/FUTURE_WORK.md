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
This is a fixed seasonal-interval mixed-layer response with a prescribed deep
reservoir, not a spun-up three-dimensional global ocean/sea-ice model.

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

- Run 147: inner atmosphere/ocean solves and snow-inclusive hydrology converged,
  but the outer fixed point did not (six-pass limit; final relative wind/heating/
  rainfall changes 0.216/0.231/0.318). First numerical follow-up: investigate
  nonlinear heating/rainfall feedback, stability diagnosis and adaptive coupling
  relaxation; do not remove a physical component just to recover a climate score.
- Ocean ice-covered cells still need a freeze/melt enthalpy closure and a clear
  distinction between ice-skin temperature and liquid mixed-layer SST. The
  reduced reservoir presently inherits subfreezing polar temperatures; closed
  heat accounting alone does not make those values physically valid ocean SST.
- Earth moisture/precipitation/Koppen, surface/upper ERA5 vectors, pressure
  decomposition, jets, shear, ascent, ocean currents and SST comparisons.
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
