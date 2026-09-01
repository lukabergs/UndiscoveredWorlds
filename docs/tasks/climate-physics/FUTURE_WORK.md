# Deferred climate-physics work

This is the durable register for physically motivated climate work that is not
active in production. Entries are candidates, not assumed improvements. Retain
one only after deterministic tests, conservation checks, an Earth benchmark,
and visual inspection show that its implementation is better than the retained
checkpoint.

## Recommended dependency order

1. Define the cell-centred internal grid contract.
2. Implement and validate spherical finite-volume MPDATA on that grid.
3. Remove the 94% transport limiter only after the replacement proves positive
   and stable without it.
4. Replace surface-temperature stationary-wave forcing with diagnosed diabatic
   heating.
5. Give zonal/barotropic and stationary/baroclinic circulation separate
   closures, then couple them without collapsing either mode.
6. Add a separately closed upper-layer stationary/baroclinic response.
7. Evaluate a prognostic shallow-water climate atmosphere after the steady
   mode-separated system is understood.
8. Replace prescribed ocean-current belts with a wind-driven coupled ocean/SST
   response after the surface winds are credible.
9. Add a separate time-evolving weather anomaly layer for gameplay.

The particle-renderer upgrade below is independent and can be completed at any
point without changing climate physics.

## Cell-centred internal climate grid

**Status:** candidate refactor.

Use `W x W/2` cell-centred latitude-longitude fields for finite-volume climate
state. Keep arbitrary requested raster dimensions as an input/output concern;
do not obtain the new convention by merely deleting a polar row.

Required work:

- centralize coordinate, cell-area, face-length, and pole-boundary definitions;
- remap pole-inclusive inputs and requested outputs explicitly;
- migrate transport, pressure, energy, validation, and diagnostics together;
- preserve support for requested outputs such as `512x257` and `2048x1025`;
- test constant-field remapping, global integrals, pole symmetry, seam wrapping,
  determinism, and water/energy budgets.

## Spherical positive-definite MPDATA moisture transport

**Status:** not implemented. The rejected CFL-subcycled candidate was not a
complete MPDATA implementation.

Implement donor-cell transport followed by MPDATA corrective passes using
spherical face metrics and cell areas. Use local face fluxes, CFL-controlled
substeps, a positivity/monotonicity option, and the same conservative treatment
for zonal and meridional transport. Validate the numerical method independently
before judging it through rainfall.

This is intended to replace, not wrap, both production transport directions:

- replace direct multi-cell zonal centroid jumps with fluxes through every
  intervening longitude face;
- replace the capped adjacent-row meridional pass with the same area-weighted
  finite-volume formulation;
- choose substeps from the multidimensional spherical CFL condition, including
  narrow high-latitude zonal cells;
- evaluate time-centred local face winds on every substep so spatially varying
  flow, convergence, and divergence affect the complete route;
- remove `transportMaxFraction` only after positivity and stability tests pass.

MPDATA reconstructs accurate face fluxes rather than exact moving parcel
polygons. If CFL subcycling proves too expensive, compare it with a conservative
flux-form semi-Lagrangian alternative using RK2/RK3 trajectories and swept-area
remapping; do not return to point-centroid teleportation.

Required tests:

- solid-body rotation and reversible deformational-flow tracer tests;
- constant-field preservation and periodic-seam continuity;
- non-negativity and area-weighted mass conservation;
- convergence with grid refinement and lower diffusion than production;
- analytical centroid and shape error for Courant numbers below and above one;
- spatially varying convergent, divergent, and sheared wind fixtures;
- no multi-row polar-boundary bypass;
- Earth moisture, precipitation, and Koppen comparisons against the retained
  checkpoint.

A successful implementation should reduce numerical blur. It cannot correct a
wrong wind field, condensation closure, or physical mixing coefficient, and it
is not required to improve every aggregate climate score on its first tuning.

## Replace the 94% moisture-transport limiter

**Status:** active production safeguard; removal deferred.

`transportMaxFraction = 0.94` is an empirical stability/positivity limiter, not
an atmospheric parameter. Do not simply change it to 100% in the current
super-CFL centroid remap. Remove it only when a local-flux, CFL-stable,
positive-definite transport method makes it unnecessary.

Acceptance criteria:

- zero negative moisture without post-hoc clipping;
- no empty-cell instability or polar leakage;
- area-weighted water conservation at least as good as the retained checkpoint;
- sharper standard tracer solutions without grid-scale oscillation;
- stable Earth runs over the complete seasonal integration.

## Diagnosed diabatic-heating stationary-wave source

**Status:** candidate circulation feature.

Replace the present surface-temperature proxy for non-zonal stationary-wave
forcing with a vertically projected atmospheric heating budget. Candidate terms
include radiative heating/cooling, surface sensible heat, and latent heat released
by condensation. Avoid counting the same latent energy in both the hydrology
feedback and the stationary-wave source.

The field diagnosis itself should be inexpensive at the current `64x33`
stationary-wave resolution. The main cost risk is an outer circulation-moisture
iteration if winds alter rainfall and rainfall alters latent heating. Begin with
a deterministic lagged seasonal-mean source; add bounded fixed-point iterations
only if diagnostics demonstrate that they are needed.

Acceptance criteria:

- documented units and vertically projected energy budget;
- zero or controlled zonal-mean contribution according to the mode split;
- convergence independent of initial guess if an outer iteration is used;
- improved stationary pressure, Walker/monsoon structure, and surface-vector
  diagnostics without degrading conservation;
- no observed-Earth data used as forcing for generated worlds.

## Mode-separated zonal and stationary circulation

**Status:** partial working solution only. The retained surface solver preserves
the existing zonal mean and solves the non-zonal stationary anomaly; the tested
full-field mode remains disabled.

Replace the shared full-field closure tested in run 137 with separate physical
closures before coupling:

- treat the zonal-mean overturning/barotropic mode with its own equivalent depth,
  damping, meridional energy transport, and angular-momentum constraints;
- treat the non-zonal stationary/baroclinic mode with its own equivalent depth,
  damping, and diabatic/orographic forcing;
- superpose or weakly couple the modes while diagnosing energy, momentum, mass,
  and row-mean transfers explicitly;
- retain isolation switches for zonal-only, stationary-only, surface-only, and
  upper-only experiments;
- require that coupling does not reproduce run 137's zonal-pressure collapse or
  precipitation-correlation regression.

Validate zonal and stationary pressure RMS/correlation separately, alongside
wind vectors, Hadley-cell width/strength, gyres, ascent, precipitation, and
conservation. Preserving the zonal mean is the safe current checkpoint, not the
final coupled solution.

## Upper-layer stationary and baroclinic circulation

**Status:** main unresolved circulation candidate; keep isolated from retained
surface physics until independently validated.

Replace the present upper diagnostic linear balance with a separately closed
upper/free-tropospheric mode. Do not copy surface quadratic drag aloft. Required
work:

- define what pressure/geopotential level or vertical mode the stored upper wind
  represents and use consistent hydrostatic units;
- derive its own equivalent depth, deformation/adjustment length scale,
  spectral bandwidth, and linear wave/eddy damping;
- project diagnosed radiative, sensible, and latent heating onto the upper
  baroclinic mode rather than reuse the surface-temperature pressure anomaly;
- represent vertical momentum and mass coupling to the surface explicitly while
  retaining a switch that isolates each layer during diagnosis;
- derive upper orographic forcing from vertically propagating stationary-wave
  physics and stability, rather than copying the surface windward-high/lee-low
  pressure term;
- preserve distinct surface and upper closures while solving their shared
  pressure, divergence, shear, and ascent relationships consistently;
- reserve freely evolving eddies and transient Rossby waves for the gameplay
  weather model instead of hiding them in seasonal deterministic noise.

Validation must compare seasonal upper `u/v`, geopotential/pressure, jet latitude
and strength, horizontal length scale, vertical shear, ascent, and surface-upper
coupling against ERA5. Require deterministic analytic tests for equatorial
regularity, hemispheric reversal, mode separation, zero-forcing behavior,
solver convergence/fallback, and response to idealized heating and topography.

## Prognostic shallow-water climate atmosphere

**Status:** unimplemented architecture candidate; not required for the retained
steady surface solution.

Evaluate a forced, damped one- or two-layer shallow-water atmosphere that advances
mass and momentum through time and diagnoses seasonal statistics after spin-up.
This would unite transient adjustment, multiple atmospheric states, Rossby-wave
propagation, and vertically coupled circulation more physically than four direct
seasonal equilibria. It would also move the simulator substantially closer to an
intermediate-complexity GCM and require many more timesteps.

Prototype at `64x32` with deterministic forcing and fixed seeds. Compare runtime,
spin-up length, conservation, seasonal climatology, directional consistency, and
Earth reference skill against the steady mode-separated system. Reuse it as the
gameplay anomaly engine only if the climate and runtime contracts are compatible;
otherwise keep the offline climate and runtime-weather solvers separate.

## Stationary-solver parameter diagnosis and scaling

**Status:** future robustness work; current `64x33` surface solve converges and
does not require a solver change.

- diagnose or constrain equivalent depths from the simulated stratification and
  vertical-mode definition instead of treating 400 m as permanently universal;
- verify adjustment length, damping time, and retained spectral bandwidth across
  planet radius, rotation, gravity, atmospheric mass, and resolution changes;
- add bandwidth- and grid-refinement studies so solver resolution increases only
  when the physical forcing contains additional resolved modes;
- record restart-cycle residual histories and detect restarted-GMRES stagnation;
- compare larger restart windows, unrestarted diagnostic solves, ILU, and
  geometric/algebraic multigrid preconditioning when upper/full-field coupling
  makes Jacobi preconditioning inadequate;
- retain physical-residual acceptance and safe fallback rather than applying a
  plausible-looking unconverged field.

Do not optimize the linear algebra below the uncertainty of the forcing and mode
closure. Stronger solvers are justified by failed convergence, scaling cost, or
new physical bandwidth—not by pursuing machine precision for its own sake.

## Time-evolving gameplay weather anomaly layer

**Status:** future gameplay system, separate from the equilibrium climate core.

Use the seasonal climatology as the background state and evolve low-resolution,
zero-mean anomalies over it. Prefer a forced/damped one- or two-layer
shallow-water or quasi-geostrophic solver following the simulated jets. A still
cheaper fallback can evolve a small set of Rossby modes with the physical
dispersion relation, stochastic excitation, and damping.

This layer may generate moving Rossby waves, fronts, cyclones, and temporary
rainfall anomalies while leaving the long-term climate map unchanged. Prescribed
sinusoids or painted storm tracks should be treated as fallback heuristics, not
the preferred design.

Statistics and validation contract:

- retain at least 20-30 statistically distinct, temporally weighted states per
  season for rough mean-vector, mean-speed, and directional-consistency maps;
- use more samples when diagnosing storm-track variability and report sampling
  uncertainty rather than treating 20-30 as a universal sufficient count;
- generate the states from the evolving anomaly model, not duplicated seasonal
  means or independent arbitrary perturbations;
- use fixed seeds and serialized model state so every sequence is replayable;
- compare the simulated statistics with ERA5 statistics calculated using the
  same season definitions and temporal sampling.

Required safeguards:

- fixed-seed replay and save-game serialization;
- anomaly means relax toward zero over climate timescales;
- mass/energy bounds appropriate to the reduced model;
- timestep stability and graceful level-of-detail scaling;
- offline statistical checks that the runtime weather recovers the parent
  seasonal climatology.

## RK2 adaptive particle-trail renderer

**Status:** visualization improvement; climate-neutral.

Replace the current three-hour forward-Euler particle steps and 1.75-cell hard
displacement cap with midpoint RK2 and adaptive substeps selected from local
grid spacing and curvature. Preserve deterministic seeding and existing density
controls. Encode elapsed integration time separately from line length so fast
flow is not silently slowed for display.

Validate straight uniform flow, solid-body rotation, a known vortex, seam
wrapping, polar behavior, and agreement with LIC streamline topology. Keep the
renderer explicitly diagnostic: it must not be reused as the moisture solver.

## Wind-driven ocean circulation and coupled SST

**Status:** candidate replacement for prescribed ocean-current belts.

The current ocean map prescribes latitude bands, equatorial/counter-currents,
and coastal boundary-flow strengths, then smooths and blocks them around land.
It is generated before surface winds. SST receives a one-pass upstream sample
and explicit western-boundary warming/eastern-boundary cooling; there is no
closed wind-current-SST-atmosphere feedback.

After surface winds are credible, replace this with a reduced seasonal ocean:

- derive surface wind stress using the same air-sea drag contract as the
  atmosphere;
- solve a basin-bounded barotropic streamfunction or linear shallow-water
  circulation with Coriolis, drag, coastlines, and bathymetry;
- represent Ekman transport and coastal/equatorial upwelling;
- advect and diffuse mixed-layer heat conservatively with the solved currents;
- feed SST back into atmospheric heating, evaporation, pressure, and wind;
- iterate the seasonal atmosphere-ocean state with bounded under-relaxation to
  a documented residual rather than imposing named gyres or basin masks.

Validate mass continuity, no-normal coastal flow, energy conservation, western
boundary intensification, equatorial/eastern-boundary upwelling, deterministic
coupled convergence, and Earth SST/current patterns. Keep a one-way mode for
diagnosis so atmospheric and oceanic errors can be isolated.
