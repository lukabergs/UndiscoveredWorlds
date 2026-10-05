# World model

Accepted on 2026-10-01. This document describes the target structure for world
state, simulation processes, run recipes and evaluation, and the migration from
the `planet` class. The migration status at the end says what is implemented.

## Why

`planet` (`src/world/planet.hpp`) is a single object that holds every simulated
field, the solver state of climate, tectonic and social entity lists, and more
than 100 presentation settings (colours, shading, marbling). Its 82 grids are
fixed `[2048][1024]` arrays, so a world costs 366 MiB at any resolution and a
workbench snapshot copies all of it. Every module reads and writes single cells
through accessors, so nothing declares what it uses: the climate-only mode can
only check after the run that elevation did not change, and the Earth climate
benchmark runs terrain-editing stages. Stage execution lives in the ImGui
workbench, so headless tools link the desktop stack.

## Concepts

The model separates four things that the current code mixes.

| Concept | What it is |
|---|---|
| Component (sphere) | State only: fields on a grid with units, entity tables, scalars |
| Process | One operation on the world; declares the fields it reads and writes |
| Recipe | Ordered or repeated process invocations with parameters and resolution |
| Presentation and I/O | Palettes, maps, exports, saves and the GUI; never simulation state |

### World and components

```text
World                  one planet at one point in its history
├─ PlanetBody          radius, gravity, rotation, tilt, orbit, tides (constants)
├─ SurfaceGrid         cell-centred W × W/2; models may run on their own grids
├─ Geosphere           elevation, crust type/age/thickness, sediment, landforms;
│  │                   plates, boundaries, deforming regions, mineral deposits
│  └─ Mantle           prescribed flow field, plumes (planned), heat flow
├─ Hydrosphere         ocean state (SST, currents); rivers, lakes, wetlands,
│                      deltas, tides, soil water
├─ Cryosphere          sea ice, land ice, snow
├─ Atmosphere          a column of levels (10 m, 850 hPa, 500 hPa today):
│                      temperature, precipitation, pressure, winds, humidity
├─ Biosphere           biomes, vegetation, habitats
├─ Anthroposphere      settlements, routes, polities, ownership, trade, history
└─ Provenance          which step wrote each field, with which parameters
```

Components hold data only. A thin common base handles field registration,
units, grids, serialization and provenance. Behaviour never lives on a
component, because many processes cross spheres: aeolian landforms read
Atmosphere and write Geosphere, rain-driven erosion reads Atmosphere and
Hydrosphere, isostasy under ice reads Cryosphere. The world component is called
Geosphere, not Lithosphere: it includes the prescribed mantle, and the
tectonics engine's internal solver class is already called `lithosphere`.

In the code today, `planet` holds one data struct per component
(`src/world/planet_components.hpp`: planet body, geosphere, hydrosphere,
cryosphere, atmosphere parameters, biosphere, anthroposphere, presentation
settings) and the climate module's `climateworld`. Grids are `worldgrid<T>`
(`src/core/world/world_grid.hpp`), sized to the world. Solver state that is not
world state (the tectonics engine between its two stages, the climate solve's
diagnostics) is `transientstate<T>` (`src/core/world/transient_state.hpp`): a
copy is empty and the consumer re-creates it.

The mantle is a boundary condition, not a simulation: the engine's seeded
analytic flow field (`tectonic_dynamics.cpp`, `mantle_flow`) and a constant
density for isostasy. Mantle convection would need a 3D Stokes solver over
hundreds of Myr and is out of scope. Plumes fixed in the mantle frame are an
opt-in engine source (`plume_processes`, `--plume-*` flags, off by default): a
steady tail flux and an optional head pulse of basalt, a thermal swell of the
ocean floor above them and optional weakening for progressive rifting. Plates
moving over them carry age-progressive chains. Plume-related deposits and
feeding plume positions to the app's volcano stage are later work.

### Fields

A field has a stable id (`component.name`), element type, layout, units and a
description. Layouts are a surface grid, a seasonal stack (4 × grid) or a
monthly stack (12 × grid). Grids are cell-centred, `W x W/2`, row-major in
exports. Continuous quantities target float32 in physical units. The fields the
climate module owns (temperature, precipitation, pressure, winds, vertical
motion, moisture, evaporation, ocean currents and sea surface temperature) are
stored that way since the Phase 1 re-baseline: the climate processes and their
validation read and write the floats, and the former `int16` encodings (whole
°C, whole m/s, scaled "storage" values) survive only as the integer accessors
of `climateworld` (`climatelegacyscale`), which the terrain, hydrology,
resource and society code, the `.uww` save format and the renderers still read.
The base elevation is float32 metres above the datum (with a separate sea
level): tectonics stores unrounded heights and the climate module reads them
exactly, while the terrain and hydrology passes still read and write whole metres
(`wholemetres`), and the land/sea mask is taken from the whole-metre value so that
every consumer sees the same coast. The elevation added on top of the base
(mountains, canyon uplift, craters, volcanoes, ocean ridges), lake surfaces and
river flow keep their integer types until the code that produces them is
reworked; category codes and masks are discrete and stay integer or boolean. A
process that moves from an integer accessor to the float one changes its results
and needs a new baseline for its module.

The catalogue is `src/core/world/world_fields.hpp` (one X-macro list). Other
code refers to fields by `worldmodel::FieldId`; bundles and reports use the
string id.

### Processes

A process is one operation with a declared access set. Kinds:

| Kind | Meaning | Examples |
|---|---|---|
| initializer | creates fields from a seed and parameters | starting crust and terrain, inherited deposit provinces |
| evolver | advances state through time | tectonics, erosion, climate spin-up |
| diagnoser | derives fields from existing state | flow routing, Köppen classes, biomes, drainage basins |
| detailer | adds plausible small-scale detail from noise | roughness, archipelagos, lake shapes, coastline breakup |
| editor | applies a user-directed, non-physical adjustment | terraforming, inland-sea and basin policies, smoothing |

Rules:

1. Only initializers and detailers create values from seeds. Every other process
   transforms existing state. A process that draws random numbers is marked
   stochastic and gets its own stream derived from the world seed and its step.
2. A coupled solve is one process. Winds, pressure, ocean currents, SST and
   rainfall are iterated to a common fixed point (`convergeclimatecoupling`), so
   Climate is one process with internal sub-models, not five stages.
3. A process that must follow another's time steps observes it rather than being
   merged into it. Mineral formation records deposits during tectonics through
   `TectonicFrame` and stays a separate module: coupled in time, separate in code.
4. Observed Earth data never forces a procedural world. Reference substitution
   (for example the 850 hPa replay in `climate_winds.cpp`) is allowed only in
   explicitly labelled attribution experiments.

The registry is `src/pipeline/world_processes.*`. Each existing generation stage
is one process during the migration; the access sets are declared there and
checked by the scheduler's audit mode (fields whose contents change during a
step must be in its write set).

### Recipes and the scheduler

A recipe lists steps (process id, label, parameters) and may freeze components.
Before running, the scheduler checks that every field a step reads is either
produced by an earlier step or supplied by the recipe's inputs, and that no step
declares a write into a frozen component. With `--audit-access` it also hashes
every field before and after each step and fails on undeclared writes. Each run
records provenance per field.

A run can cover part of a recipe (`world_pipeline run --steps a..b` or `--only`)
and start from the world of a bundle (`--from-bundle`); checkpoints write a bundle
after chosen steps. The same read check then decides whether the start point is
valid: a field or table that an earlier, skipped step provides must have been
loaded from the bundle. Tables, the climate solver's grids, levels and
diagnostics, and `GenerationScratch` state passed between steps cannot be loaded
yet, so a run cannot start between their writer and a reader (migration step 5).

Built-in recipes (`src/pipeline/world_recipes.*`):

| Recipe | Purpose |
|---|---|
| `procedural-world` | Full generation from a seed: tectonics, terrain passes, climate, hydrology, layers, optional society |
| `imported-world` | Imported terrain followed by the remaining terrain passes, climate and hydrology |
| `earth-climate` | Imported Earth terrain, Geosphere frozen, climate processes only; the climate validation recipe |

Coarse and fine climate are two invocations of the Climate process with
different working grids; a later grid-refinement step can raise the surface
resolution between them.

## Atmosphere levels

The governing equations are the same at every height; what differs is which
terms dominate (friction and surface fluxes in the lowest kilometre, a nearly
geostrophic flow above it outside the tropics). The current model represents
heights differently: the 10 m and 500 hPa modes share one stationary solver with
per-level constants (`climate_atmosphere.hpp`, `ModeSeparatedCirculationConfig`)
and two tuned couplings; 850 hPa is a bulk boundary-layer reservoir; radiation
is a two-layer grey column.

Storage is generic: `climateworld::levels()` is a column of `atmospherelevel`
records (name, nominal pressure, seasonal east and north wind, geopotential
height, vertical velocity, availability), filled at 10 m, 850 and 500 hPa after
each wind solve, matching the ERA5 references (`atmosphere.level_*` fields and
the `atmosphere.levels` table in bundles). Temperature and humidity per level
follow when a solver provides them.
A multi-level solver is a research track: either a vertical-mode decomposition
whose modes reuse the existing stationary solver, with couplings and heating
responses projected from friction and heating profiles, or a simplified
time-stepping model (SPEEDY, PlaSim/ExoPlaSim). Either is a new model to be
compared with the current one on the same references.

The first step is an experimental, diagnostic option, `--climate-generic-levels
925,850,700`: levels between the boundary layer and 500 hPa solved by the
850 hPa inversion solver with one rule set instead of per-level constants (a
troposphere-deep mode, the free-tropospheric drag, ideal-gas density at the level
and the surface pressure pattern projected with ln(p/p_r)/ln(p_s/p_r), reversing
at the heat-low depth p_r = 760 hPa). They appear as `generic-<p>hPa` levels; the
model never reads them back, so all other fields are unchanged. Against ERA5 at
850 hPa (seed 1, 512 x 256, 128-wide atmosphere, eight coupling updates) the
generic level's mean loss is 0.620 against 0.603 for the tuned level: one rule set
reaches the tuned level's skill within 3 % but does not improve it. Both 850 hPa
fields are the weakest level (speed correlation 0.13 to 0.25, against about 0.75
at 500 hPa). The current levels do not lie on one smooth vertical profile (drag
about 0.5 days at the surface, 5 days at 850 hPa, 2 days at 500 hPa), and the
500 hPa mode uses another solver form, so a common class for all three reference
levels needs the upper mode recast in the same form.

## Climate module

`climateworld` (`src/simulations/climate/climate_world.hpp`) is the climate
module's world: the climate state (seasonal and legacy climate fields, sea ice,
biomes, subtropical highs, weather state, solver context, levels) plus a bound
`climateinputs` view of terrain and planet parameters. planet holds one and
forwards its climate accessors; `planet::climate()` binds the planet's terrain.
`pipeline/climate_stages` runs every climate stage on a climateworld. The
library set `uw_climate_world` (core, I/O, climate, pipeline and validation, without `uw_world`) needs no
planet, terrain, hydrology or tectonics code; the test binary `world_tests`
links only it. `world_pipeline run --recipe earth-climate` is the one way to run
climate on imported terrain. Tides are hydrology and stay with planet.

## Geosphere initialization and minerals

An initializer creates the starting Geosphere before tectonics: crust type,
thickness and age, elevation, plate seeds, plumes and inherited deposit
provinces. Random, supercontinent, sketch and continuation initializers write
the same fields; tectonics always starts from world state.

The stage "Initial geosphere" (`geology.geosphere_initialization`) is that
initializer today: it creates the tectonics engine, which seeds its own starting
topography (`lithosphere::seedInitialTopography`), begins the mineral history's
inherited provinces (`MineralHistory::begin`) and writes the starting elevation.
"Plate tectonics" continues the same engine. Starting tectonics from arbitrary
world fields (sketch or continuation initializers) needs an engine entry point
that builds plates and columns from fields; the engine has none yet.

Deposits split by how they form. Magmatic and hydrothermal deposits (porphyry,
epithermal and orogenic gold; later plume-related types) need the tectonic
setting at the moment they form and are recorded during tectonics. Surface
deposits (placers, laterite and bauxite, evaporites, coal) depend on climate,
drainage and vegetation and belong to later passes.

## Evaluation

Any run can carry three kinds of evidence:

- **checks**: automatic pass/fail on budgets, finite values, ranges, grid
  registration and determinism;
- **validation**: scores against observed Earth data, cell by cell (climate,
  ocean, runoff, rivers, Köppen) or as distributions where maps cannot match
  (tectonic hypsometry, deposit counts per area);
- **review**: the user's judgement against an approved run, recorded with it.

Modules weight these differently. Tectonics and terrain are judged mainly by
review; climate, atmosphere and hydrosphere by validation, with user visual
overrides. A module can be validated against Earth only if it runs on Earth
inputs alone: everything it reads imported, detailers and editors off, nothing
written outside its own components. Reviewed modules need exact
reproducibility.

The evaluation module (`scripts/evaluation/`) works on run bundles. Each module
supplies a profile (`configs/evaluation/<module>.json`) instead of its own
scripts: the tool and arguments that produce a run, its checks, its validation
comparisons, its review maps and how it weights the three kinds of evidence.

| Part | Script | Output |
|---|---|---|
| experiment | `experiment.py` | run directory: `run.json` (command, seed, code state, timing), log, tool outputs, `bundle/` |
| evaluation | `evaluate.py` (`checks.py`, `validation.py`, `references.py`, `maps.py`) | `evaluation.json`, `report.md`, `index.html`, `maps/` |
| review | `review.py` | `review.json`: verdicts tied to the bundle's manifest hash |

Map comparisons reuse the climate report's per-field metrics, so field losses are
defined as in global-balanced-v1. With every weighted field present the score is
reported as global-balanced-v1 on bundle fields (same formula and weights, the
solvers' native grids); otherwise as a partial score with its coverage. Climate
runs registered in `runs/registry/climate.json` keep their ranking through the
climate experiment archives, which extract inputs on other grids, so the two
values are not interchangeable. Terms:
run, check, validation (skill scores), review; benchmark is kept for fixed
standard experiments such as the Earth climate benchmark.

## Run bundles

A bundle is a directory: `manifest.json` (schema `uw-run-bundle/v1`: seed,
recipe, steps with timings, grid, code revision, per-field id, units, encoding,
dtype, shape, file and provenance) and one NumPy `.npy` file per field under
`fields/<component>/`. Diagnostic fields of the climate solve (column water,
moisture transport and convergence, precipitation rate, ascent, currents,
mixed-layer temperature, storage depth) are on the solvers' own grids and appear
only when the run that recorded them writes the bundle (a run started from a
bundle copies them from it unless it solves climate itself). `world_pipeline --bundle <dir>` writes one;
`scripts/evaluation/compare_bundles.py` compares two field by field and is the
equivalence check for refactors.

## Migration

| Step | Content | Status |
|---|---|---|
| 1 | Field catalogue, run bundles, process registry with declared access and kinds, recipes and scheduler with audit; stage code out of the workbench; headless core without ImGui/SFML | done |
| 2 | Climate inputs and outputs; climate code off `planet`; level-based atmosphere storage; `climate_pipeline` | done |
| 3 | Evaluation module on bundles (climate and world profiles); `uw.py` delegates to it | done |
| 4 | Components behind the `planet` facade, resolution-sized grids, Geosphere initializer; copy-on-write snapshots deferred | done |
| 5 | Climate diagnostics in the world (done); per-step random streams, `GenerationScratch` contents into fields, save format; regional map stays on the facade | in progress |
| A | Mantle plumes in the engine (opt-in, off by default; awaiting visual review) | done |
| B | Multi-level atmosphere: generic levels from one rule set (experimental diagnostic option); full solver open | started |

Equivalence evidence for each step is recorded in
`docs/tasks/main/world-model-migration.md`.
