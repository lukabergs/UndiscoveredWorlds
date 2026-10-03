# World model migration log

Evidence for each step of the migration in
[world model](../../architecture/world-model.md). Equivalence runs use seed 1,
unregistered `world_pipeline generate` runs with every output path redirected to
a scratch directory (nothing under `runs/` changes), and compare saved worlds
(`.uww`, floats at about six significant digits) and, from step 1a on, run
bundles at full precision (`scripts/evaluation/compare_bundles.py`).

Cases:

| Case | Command (all with `--seed 1 --climate-coupling-iterations 2 --no-profile-log`) |
|---|---|
| A procedural | `generate --resolution 512 --climate-resolution 64 --social` |
| B Earth, benchmark stage path | `generate --import-land/--import-sea earth_*_512x256.tif --resolution 512 --climate-resolution 64 --no-rivers --no-lakes --no-deltas` |
| C Earth, climate only | `generate --climate-only` with the same imports and `--map none` |

The baseline is the working tree of 2026-10-01 before the migration (the Release
`world_pipeline` built at 01:25, identical to a fresh build). A repeat of case C
with the baseline binary produced a byte-identical saved world, so these runs
are deterministic.

## Step 1a: field catalogue and run bundles

`src/core/world/world_fields.hpp` lists 148 fields and 10 tables;
`planet_fields.cpp` reads them from `planet` without changing it;
`world_pipeline --bundle <dir>` writes them as NumPy files with a manifest.

- Saved worlds of A, B and C are byte-identical to the baseline.
- Bundles of these runs are the references for later steps.

## Step 1b-1e: processes, recipes, headless core

- Stage code moved from the ImGui workbench to `src/pipeline/world_stages.cpp`;
  the workbench passes `WorldStageParameters` (the hypsometric curve is sampled
  by the editor, as before).
- Recipes (`procedural-world`, `imported-world`, `earth-climate`) and the
  scheduler with `--audit-access` live in `src/pipeline/world_recipes.cpp`.
- `world_pipeline` links the headless `uw_world_core`; it depends on libpng and
  the C runtime only (no SFML, ImGui). PNG I/O outside the desktop uses
  `src/io/raster_image.*` (libpng).
- Climate validation is compiled on its own instead of being included into the
  debug module.

Results with the recipe binary:

- B: 158 of 158 fields and tables identical to step 1a; saved world identical.
- C: 157 of 158 identical. The difference is `geosphere.roughness`: the
  earth-climate recipe no longer runs the noise-only roughness stage, which climate
  never reads. Land and sea cell counts, which that stage also computed, are now
  computed when imported terrain is prepared.
- A: 158 of 158 identical.

Declared access (generated from the analysis build) audited clean: zero undeclared
or frozen writes for A, B, C and the Earth run with rivers, lakes and deltas.
The earth-climate audit found no write into the Geosphere.

## Step 2: climate off planet

- `climateworld` (`src/simulations/climate/climate_world.*`) holds the climate
  state (39 former planet members: seasonal and legacy climate fields, sea ice,
  biomes, subtropical highs, weather state, solver context) and reads terrain and
  planet inputs through a bound `climateinputs` view. planet keeps forwarding
  accessors; `planet::climate()` binds its own terrain on every call, so world
  copies stay correct.
- The climate code (`climate_fields`, `climate_coupling`, `climate_classification`,
  `climate_energy`) and climate validation take `climateworld&`. Shared helpers
  became templates over the world type (`climate_elevation.hpp`,
  `land_distance.hpp`); seeded random numbers and fractal noise moved to
  `core/procedural_noise` (`uw_procedural_noise`).
- `pipeline/climate_stages` runs every climate stage on a climateworld; the planet
  stages delegate to it (tides stay with planet).
- `uw_climate_world` is a static library without planet; `climate_pipeline` links
  only it. Its bundle reader is `climate_world_fields`; planet's bundle reader
  delegates climate-owned fields to it.
- Atmosphere levels: `atmospherelevel` (10 m, 850 hPa, 500 hPa) is filled from the
  model's float wind fields after each wind solve; five `atmosphere.level_*`
  fields and the `atmosphere.levels` table are new catalogue entries.

Equivalence after each sub-step (2a storage move, 2b climate code on climateworld,
noise move, stage runner, bundle reader, import split): B and C identical to
step 1 (158 of 158; from the level step on, 158 of 164 with six new level items);
A identical after 2a. `climate_pipeline run` on the same Earth terrain and
settings as case C: all 82 fields and tables it writes (climate state, levels,
inputs) are identical to `world_pipeline` earth-climate.

Corrections recorded on the way: the land/sea counts, the imported-climate-map
path (`context.importedClimate` is set but no stage reads it; the code that
applied imported maps is only in the unused `generateglobalclimate`), and the
unused `planet::smoothrainmaps`, which was removed.

## Step 3: evaluation module

- `scripts/evaluation/`: `bundle.py` (reader, grid geometry, area weights, sea mask),
  `checks.py`, `references.py` (`.uwclim` reader, month/quarter/annual rules,
  conservative remap), `validation.py` (map comparisons through the climate
  report's `evaluate`, hypsometry through the tectonics Earth reference),
  `maps.py` (climate and tectonics palettes), `evaluate.py`, `review.py`,
  `experiment.py`; profiles `configs/evaluation/climate.json` and `world.json`;
  `bench.py <climate|world> experiment|evaluate` delegate to them.
- Seed-1 checks on scratch bundles: the climate profile on the climate_pipeline
  bundle of step 2 runs in 9 s (checks, 8 comparisons, 13 maps); the world profile
  on case A passes every check, with the repeat check over 148 identical fields.
- Bundles lack column water, moisture transport and convergence and currents in
  m/s, so a bundle score covers 70 % of the global-balanced-v1 weight and is
  reported as a partial score with its coverage, never as global-balanced-v1. The
  numbers from the quick case-C settings (64-wide atmosphere, two coupling
  iterations) say nothing about skill.
- Found and fixed: `climate_pipeline` bundles had no per-field provenance (they
  now record the last declared writer); the access generator derived table
  enumerators from their qualified names and emitted `T::levels` for
  `atmosphere.levels` (it now reads the mapping from the catalogue).
- Declared access regenerated with the level fields (the first regeneration used
  an analysis binary built before them). Audits with the declared access: A, B, C
  and D report no undeclared writes; C and the climate_pipeline bundle are
  identical to step 2 (164 of 164, 82 of 82).

## Step 4: components, resolution-sized grids, Geosphere initializer

- 4a: `worldgrid<T>` (`core/world/world_grid.hpp`) holds exactly the world's cells,
  x-major so `grid[x][y]` reads as before; the 68 planet grids and the 13 climate
  grids use it, so memory follows the resolution (about 160 bytes per cell for the
  planet grids: 21 MiB at 512 x 256 instead of 335 MiB at any size). The
  subtropical-high rows keep the 2048-column capacity of the save format.
- 4b: planet's members are grouped by component in `wip/planet_components.hpp`
  (`planetbodystate`, `geospherestate`, `hydrospherestate`, `cryospherestate`,
  `atmosphereparameters`, `biospherestate`, `anthropospherestate`,
  `presentationsettings`); climate state stays in `climateworld`. Accessors are
  unchanged.
- 4d: stage "Initial geosphere" (`geology.geosphere_initialization`, initializer)
  creates the tectonics engine and the inherited deposit provinces and writes the
  starting elevation; "Plate tectonics" continues that engine. The engine travels
  in the stage scratch as `transientstate` (copies are empty), so a workbench
  commit or a recipe without the initializer re-creates it deterministically.
- Copy-on-write snapshots are deferred: about 60 % of a snapshot is climate solver
  state in plain vectors, and sharing it safely needs equality or change tracking
  per component; declared writes are not a safe basis for the workbench's own edit
  paths. With resolution-sized grids a 1024 x 512 snapshot is roughly 250 MB
  instead of 540 MB.
- Exact sizes exposed reads beyond the grids that the zero-filled capacity had
  hidden. A bounds-checking build (`-DUW_CHECK_GRID_BOUNDS`, stack trace on the
  first violation) found `planet::smoothoverland`/`smooth` reading rows j +- amount
  beyond the poles; those rows now count as zero, which keeps the former averages.
  `planet::deltajul` lacked the bounds check of its siblings (fixed). Resized grids
  now take the values clear() sets (owners -1), because imports size the world
  after clear(); clear() resets every subtropical-high column again.
- Equivalence with step 3 (seed 1): A, B, C, D identical (164 of 164) with
  byte-identical saved worlds; the climate_pipeline bundle identical (82 of 82).
  The bounds-checking build runs A, B, C and D without a violation.

## Step 5a: climate diagnostics in the world

- The per-season process, ocean and weather statistics of the last solve moved from
  statics in the validation module (keyed by the world's address) into
  `climateworld` (`climatediagnostics`, transient: copies do not carry it, as
  before). The climate solver no longer calls the validation module for them.
- Ten diagnostic fields on the solvers' own grids enter the catalogue and bundles:
  column water, moisture transport east/north, moisture convergence, precipitation
  rate and ascent (process grid), currents east/north, mixed-layer temperature and
  storage depth (ocean grid). They are absent until a solve records them.
- Case C and the climate_pipeline case: every former field identical, the ten new
  fields present and identical between `world_pipeline` and `climate_pipeline`;
  saved world byte-identical.
- The climate profile now compares all eleven weighted fields (water, moisture
  transport, convergence and currents from the diagnostics; rain and SST on the
  process and ocean grids, as the climate report does). The score is labelled as
  global-balanced-v1 on bundle fields: same formula and weights, other input grids
  than the registered rankings.

## Step 5b: a random stream per process

- Every stage now seeds the shared generator from the world seed and an FNV-1a
  hash of its process id (`seedprocessstream`, `pipeline/world_processes`), in
  both stage runners. Skipping, repeating or reordering a step no longer moves the
  random draws of the others.
- This changes results wherever a stage drew from the stream left by an earlier
  one. Unchanged: the tectonics engine and the mineral history (own generators),
  every climate stage (each reseeds before drawing), so the earth-climate recipe
  and `climate_pipeline` give identical bundles (case C 174 of 174, climate
  pipeline 92 of 92). Changed: terrain and hydrology passes, so seed-1 procedural
  worlds differ from the shelf rebuild after tectonics on (continents are the
  same), and the imported-world path differs too (case B: base elevation in 87 %
  of cells, mean 81 m; temperature in 1.3 % of cells). `world_pipeline
  climate-benchmark` runs that path, so its next registered run needs a new
  baseline; `climate_pipeline benchmark` does not edit terrain and is unaffected.
- Repeating case B exposed nondeterminism in `river_direction` (1 to 3 cells
  differing between runs of one binary). Access reports now carry the hash of every
  written field after each step (`write_hashes`); comparing four runs put the first
  divergence in "Basin editor", in the parallel drop tracing of
  `accumulaterivercatchments`: drops from different rows that reach the same sea
  cell wrote its river-mouth direction in thread order. The writes are now
  collected per source row and applied in row order after the trace (tracing never
  reads them). Four repeats are identical. The race predates this work; the new
  terrain made it visible.

## Step 5c: continental shelves in the world

- The shelf mask was stage scratch carried from tectonics to the volcano pass,
  ocean-depth refinement and marine resources. It is now
  `geosphere.continental_shelf` (`worldgrid<bool>` in the geosphere component,
  `planet::continentalshelves()`), visible to bundles and access audits; every
  function that took the scratch mask takes the world's grid. Case B is identical
  apart from the new field (and the river-mouth race above).
- Still scratch, still crossing stages: the mountain drainage seeds (mountain bases
  to rivers), and the salt-lake, no-lake and basin-seed grids (rivers to lakes).
  They are hydrology intermediates and need the river and lake stages restructured
  (one process, or hydrosphere fields) before they can move.

## Track A: mantle plumes

Opt-in engine source, logged in `docs/tasks/tectonics/ENGINE_FINDINGS.md` (8).
With plumes off the seed-1 procedural world is identical (175 of 175, saved world
byte-identical) and `tests/scripts/check_tectonics_cli.py` passes (it now also
rejects invalid plume options). With plumes on the column budget closes (residual
below 1e-13) and the plume mass equals the prescribed volume.

## Track B: generic atmosphere levels (experimental)

`--climate-generic-levels` adds diagnostic levels from one rule set (design doc,
"Atmosphere levels"). Off: the climate_pipeline bundle is identical (92 of 92);
on: every field except the level fields and table is identical. At 850 hPa against
ERA5 (seed 1, 512 x 256, 128-wide atmosphere, eight coupling updates) the mean
season loss is 0.620 for the generic level and 0.603 for the tuned one. The
evaluation selects levels by name (`atmosphere.level_wind_east[@generic-850hPa]`).

## Observations

- With the roughness stage removed and no land/sea counts, the climate-only run
  segfaulted in the rainfall stage: a world whose sea count is zero (a desert
  world) skips the ocean, pressure and wind stages and then crashes. Not fixed.
- `generate --resolution 256 --tectonic-duration-myr 100` spent over 16 CPU
  minutes in "FastLEM mountains" (63 ms on the 512-wide case A). Not
  investigated; that case was dropped from the equivalence set.
