# Climate roadmap: implemented features and validation queue

Implementation branch: `climate-physics_7_core_feat_260902/future-work`.
The proposed systems are implemented. This is not a claim of improved Earth
skill: the user explicitly requested implementation before benchmark tuning.
Earth comparisons, runtime comparisons and visual acceptance remain pending.

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

- `climate_atmosphere.*`: atmospheric heating in W/m2 from simulated radiation,
  sensible heat and condensation (`mm * latent heat / seconds`). Vertical
  projection removes stationary row means. No observed-world forcing is used.
- The zonal closure uses thermal contrast and Held-Hou overturning; its upper
  background is hydrostatic thickness between 1000 and 500 hPa. Non-zonal
  surface pressure and upper geopotential have independent equivalent depths,
  restoring times and solvers. Surface wind uses quadratic drag; upper wind
  uses linear damping. The old output-grid pressure/height relaxation is removed.
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
One lagged rainfall/heating pass refreshes circulation without repeating coastal
temperature postprocessing. The direct latent-heating ascent feedback is disabled
when that circulation projection is active; thermodynamic moist adjustment stays.
Isolation switches in `generation_tuning.hpp`: zonal, stationary, surface, upper,
lagged coupling and prognostic atmosphere.

Tests: all 16 isolation combinations, zero forcing, hydrostatic units, rotation
reversal/equatorial regularity, heating budgets, stability/bandwidth scaling,
mountain-wave response, residual acceptance and failed-solve fallback.

## Prognostic climate and replayable weather

`climate_weather.*` implements forced/damped, linear one/two-layer spherical
shallow-water vertical modes. Equivalent depth sets wave speed; signed modal
geopotential is not a positive total water-column thickness. RK2, local wave
viscosity, CFL/source substeps, amplitude bounds and atomic rollback protect the
reduced model. The climate path spins up then averages 30 evolving states; its
blend and independent fallback are explicit controls.

Gameplay anomalies follow interpolated seasonal jets and relax toward zero.
`planet::advanceweather(seconds, optionalHorizontalCells)` advances only the
anomaly state. It preserves climate maps, supports conservative LOD changes,
and serializes time, all layer fields and RNG state for replay. World save
version is now **17**; version-16 files are not migrated by this change.

Hydrology continues weather between months and retains 90 daily states per
Jan-Mar/Apr-Jun/Jul-Sep/Oct-Dec season in each final annual cycle. Statistics
use actual elapsed-time weights: mean vectors, mean speed, directional
consistency, speed spread and independent-sample standard error. The last is
not an autocorrelation-corrected confidence interval; storm-track inference
requires longer sampling. Circulation exports include `weather_statistics.csv`.

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
exchange is included in the budget. The bounded outer loop under-relaxes SST,
pressure feedback and winds over the same fixed seasonal interval; it reports
both basin-equation and coupled fixed-point residuals. Reaching the iteration
cap is reported as unconverged, not silently relabelled. One-way mode isolates
ocean response; invalid states fall back to atmospheric SST and zero currents.

Both generation pipelines now order winds before ocean/SST before evaporation
and rainfall. Tests cover this order, coastal faces, western intensification,
open polar oceans, variable-depth continuity, uniform SST, heat closure,
deterministic coupling and rejected incomplete basin solves.

## Remaining validation (intentionally not run between features)

- Earth moisture/precipitation/Koppen, surface/upper ERA5 vectors, pressure
  decomposition, jets, shear, ascent, ocean currents and SST comparisons.
- Multi-seed/planet-parameter and long-duration weather/energy soak suites;
  correlated-sample uncertainty and recovery of parent seasonal climatology.
- Benchmark runtime, forcing-bandwidth/grid-refinement studies and tuning of
  spin-up/iteration budgets. Some coupled worlds may hit the iteration cap.
- Visual map/LIC/particle agreement and high-resolution world generation.
- Full application save/load/runtime-weather round trip across supported hosts;
  binary weather-state replay is covered by deterministic unit tests.

Commands and completed verification are recorded in the branch `PLAN.md`.
