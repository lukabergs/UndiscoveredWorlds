# Undiscovered Worlds

Undiscovered Worlds is a world engine for generating a planet's physical environment and simulating the emergence of ecosystems, settlements, and political geography. The intended output is a world snapshot with enough natural and social history to explain its terrain, resources, population, and borders.

The engine combines procedural heuristics with more physically based simulations. Each stage should be independently configurable, inspectable, and repeatable, with computational effort chosen for its purpose.

## Intended simulation lifecycle

1. **Initial conditions.** Start from noisy terrain, another procedural initial state, or imported maps. Record the seed, planetary parameters, and input provenance.
2. **Geology and terrain.** Run plate tectonics, then selected terraforming passes: mountain building, volcanism, coastline refinement, erosion, and deposition. Mineral formation belongs to this geological history; later erosion and deposition can expose or redistribute deposits.
3. **Coarse climate and hydrology.** Run an inexpensive climate pass to estimate precipitation. Use it to shape drainage, erosion, rivers, lakes, and finer terrain detail.
4. **Detailed climate.** Run the same climate system with a larger computational budget against the revised terrain and surface-water state. Explicitly bounded feedback passes may update hydrology again when needed.
5. **Ecology and biological resources.** Simulate coarse vegetation spread and procedurally populate suitable habitats with animals. Evolutionary simulation is outside the intended scope. Combine biological resources with the existing geological resources.
6. **Early settlement.** Seed and develop populations using resource availability and desirability profiles that can differ between fictional species or races. Represent early migration and settlement before introducing more complex institutions.
7. **Settled societies.** After a configurable settlement threshold, simulate trade, political organization, administrative borders, and geopolitical relationships. Stop at the chosen social-development threshold or simulation horizon and export the world snapshot.

These are lifecycle stages, not mandatory one-time function calls. The coarse and detailed climate passes should be two configurations of one implementation. Geological, climatic, ecological, and social time scales should remain explicit.

## Current state

The repository contains terrain generation and import support, an internal plate-tectonics module with world-grid ingestion, a generation workbench and stage registry, a reduced climate simulator, physical-resource heuristics, gold deposits recorded through tectonic history, and initial settlement/trade/polity generation. Climate is the current development focus. Terrain refinement remains work in progress; the full ecological and staged civilizational lifecycle above is a design goal.

The September 2026 cleanup aligns CMake sources, includes, configuration, and preparation/analysis scripts with the new layout. Earth benchmark land/sea imports can be regenerated at arbitrary even resolutions from the retained signed Earth master.

## Repository layout

```text
src/
  core/                    Grids, seasonal indices, deterministic seeds, row workers
    world/                 World field catalogue and field-read tracking
    procedural_noise.*     Seeded random numbers and fractal noise (no world state)
  pipeline/                Stages, process registry, recipes and scheduler, progress
  simulations/
    geology/               Plate-tectonics engine, world adapter, mineral deposit history
    terrain/               Fractals, coasts, ridges, elevation, landforms
    climate/               Atmosphere, ocean, energy, moisture, transport
      legacy/              Tested alternative thermal and ocean feedback closures
    hydrology/             Drainage, basins, lakes, deltas, wetlands, tides
    resources/             Mineral and marine resource potentials, present-day placers
    society/               Settlement, routes, trade, polities, history
  io/                      Map imports, raster/reference readers, export selection
  validation/climate/      Reference comparisons and climate diagnostics
  validation/tectonics/    Tectonic snapshot statistics and records
  app/                     Desktop entry point, windows, controls
    rendering/             Map appearance and CPU/CUDA rendering
    diagnostics/           Profiling and generation debug integration
  wip/                     Mixed responsibilities awaiting extraction
assets/                    Presentation files: textures and appearance presets
configs/                   Application configuration and saved UI layout
  evaluation/              Evaluation profiles: checks, validation and review maps per module
refs/
  source/                  Original downloads, receipts, and licenses
  processed/               Prepared grids grouped by physical quantity
runs/                      Simulation runs, one directory per registered run (see Runs and comparisons)
tests/                     Core, pipeline, I/O, simulation, and script tests
scripts/
  refs/                    Reference download and preparation
  benchmarks/              Benchmark analysis and workbook input preparation
  evaluation/              Run bundles: checks, validation, review verdicts, experiments
docs/                      Architecture, reference guides, and current climate notes
vcpkg-ports/               Dependency overlay ports
out/                       Build products, editor caches, scratch work
```

Root files are limited to the README, license, version, Git ignore rules, and CMake/vcpkg entry points. Git and Codex metadata retain their operational locations.

### Code boundaries

Keep headers beside their implementations. Organize by the responsibility of a module; express execution order in pipeline recipes. Avoid creating a separate climate implementation for each place it appears in a recipe.

- `core/` should contain shared state primitives, coordinates, units, grids, and deterministic randomness. `core/world/world_fields.hpp` catalogues every world field (id, component, units); moving world storage into components is still pending.
- `simulations/` should operate on explicit world inputs and produce defined world layers. Presentation, filesystem operations, and reference-data validation have their own boundaries.
- `pipeline/` selects stages, validates their dependencies, schedules them, and controls checkpoints and stopping conditions. Every stage is a process with a kind and declared field access (`world_processes`); recipes (`world_recipes`) order them and can freeze components. A disabled stage requires either valid existing outputs, an imported replacement, or an explicit fallback.
- `io/` contains code that reads and writes artifacts. `refs/` and `runs/` contain the artifacts themselves.
- `validation/` compares simulated results with observations and checks budgets. Observed climate references are validation inputs, not hidden forcing for procedural worlds.
- `app/` presents and controls the engine. Both interactive and eventual headless workflows should use the same simulation stages.

`globalterrain`, `globalclimate`, and `physical_layers` have been split into simulation modules and pipeline entry points. Society CSV loading belongs to `io/social_definitions`; its stage order belongs to `pipeline/social_generation`. [Module boundaries and retention decisions](docs/architecture/simulation-modules.md) describe the active and legacy paths.

`src/world/` still holds `planet`, `region`, `regionalmap`, the queries on regions, world set-up, the `functions.hpp` umbrella and `generation_tuning`; the workbench UI (`generation_workbench`) is in `src/app/`. World storage/persistence, regional generation and shared utilities remain coupled there; stage execution lives in `pipeline/world_stages`. `functions.hpp` is a compatibility umbrella without SFML; map views and drawing live in `app/rendering/map_views.hpp`. New callers should include the owning module's header. The target structure (data-only sphere components, processes, recipes, evaluation) and its migration status are in [world model](docs/architecture/world-model.md).

Further land-surface hydrology and ecology belong under `src/simulations/hydrology/` and `src/simulations/ecology/`. Existing climate moisture/storage physics remains in `climate/`. Gold deposits form during tectonics in `geology/mineral_systems/`, travel with their host crust and erode as they are exhumed; `resources/` routes liberated gold into placers on the final drainage. Other mineral and marine scores remain derived heuristics in `resources/`. Biological production belongs in ecology, with society consuming both. `world_pipeline run --recipe tectonics` writes the deposit list and gold budget of a native run (mineral history is on by default; `--mineral-history off` disables it).

### Assets, definitions, configurations, and references

| Location | Role | Examples |
| --- | --- | --- |
| `assets/` | Presentation | Sphere texture, map appearance presets |
| `definitions/` | Authored simulation content | Commodity properties, knowledge definitions, species preferences |
| `configs/` | How an application or run is configured | Application paths; future recipes, resolution, duration, stage switches |
| `refs/` | Imported inputs and external reference data | Heightmaps, observed precipitation, prepared Earth grids |

The distinction follows purpose. An imported heightmap is a dataset; a generated heightmap is a run artifact. Current tuning constants remain in `src/world/generation_tuning.hpp`; moving them into configurable recipes is future work.

Keep downloads and their receipts/licenses under `refs/source/<dataset>/`. Prepared grids belong under `refs/processed/<quantity>/`; Earth rasters use `{width}` filenames and metadata under `processed/metadata/{width}/`. See [reference preparation](refs/README.md).

## Runs and comparisons

Each registered run is one directory; the registry supplies its ID:

```text
runs/
  registry/<module>.json      IDs, status, run directory, label, score; `baseline` names the module's baseline run
  ledger/<module>.jsonl       experiment ledger
  <module>/<id>/              run.json  log.txt  bundle/  checkpoints/  outputs/  evaluation/  ledger-entry.json
  reports/<name>/             cross-run comparison pages
  archive/                    sealed layout of runs registered before 2026-10-05 (maps, fields, diagnostics, manifests,
                              parameters, checkpoints, metrics, logs); read-only evidence, see runs/archive/README.md
```

`uv run --offline scripts/uw.py <module> run` makes a registered run (climate and tectonics); `run.json` holds the command, executable hash, code revision and diff hash, inputs with hashes, timings per step and the verdict, and the `evaluation/` page shows checks, validation and review maps. Commands and the directory contents: [benchmark commands](docs/reference-guides/benchmark-commands.md).

Runs registered before 2026-10-05 keep their former paths in their records: a path `runs/<kind>/...` is now `runs/archive/<kind>/...`. That includes the climate map folders (category/season/layer, run ID as file name) and the former Excel workbooks (`runs/archive/metrics/`, with the 304 historical records of `run_history.xlsx`); the climate registry records of those runs keep the C++ benchmark's metrics.

For future runs, record the seed, recipe and stage parameters, code revision, input versions, grid geometry, simulation duration, convergence criteria, and output provenance. Preserve numerical fields as well as previews when they are needed to reproduce a comparison. Changes in palettes, units, grids, or renderer versions should remain visible in the run metadata.

## Tectonics

There is one tectonics implementation: `src/simulations/geology/plate_tectonics`,
built as the internal `uw_tectonics` target from the current working tree.
The desktop generator, workbench, mineral history and the `tectonics` recipe of `world_pipeline`
all use it. No external tectonics checkout or library is needed or supported.
Existing uncommitted work on `main` is part of that implementation.

`configs/tectonics-baseline.json` is the shared production recipe. CMake compiles
it into `tectonics::scenario::baseline()`; changing it triggers reconfiguration.
The current recipe uses seed 1, a 256×128 native grid, 10 initial plates and
600 Myr at 1 Myr per update. Application generation uses the world's seed;
new benchmark runs use seed 1. The app's duration, initial plate count and
timestep overrides use 0 to select the recipe value. `--tectonic-duration-myr`
is the application's duration override; the `tectonics`
recipe uses `--duration`.

Every physical update uses reduced force-balance motion and conservative
isostatic columns. Evolution continues until the positive duration is
reached (rounded up to a whole timestep); local plate splitting, welding and
removal can occur.

Code outside the engine uses only its public headers: `tectonic_api.hpp`,
`tectonic_contract.hpp`, `tectonic_recipe.hpp`, `tectonic_scenario.hpp`,
`tectonic_world_grid.hpp`, `topography_codec.hpp`, plus `material.hpp` (material
codes) and `utils.hpp`. Engine tests may use internal headers;
`scripts/check_layout.py` enforces this.

The solver includes continental/oceanic rock and sediment transport, thermal
subsidence, distributed deformation, root-pressure forces, active margins,
trench bending, arc growth, erosion and sediment redistribution. The recipe
selects parameters within this one solver. Initial ocean depth and land relief
are positive metres relative to sea level, not limits on evolved mountains.
Column budgets track mass sources, sinks and internal transfers. Boundary,
deforming-region and junction fields describe the evolving setting. Read
`tectonic_contract.hpp` and the [native output contract](docs/architecture/tectonics/tectonic-output-contract.md)
for field meanings; thickness proxies and diagnostic scores are not ore grades.

The native domain is a cell-centred planar torus with equal-area square cells.
`world_tectonics.cpp` runs the recipe at its native resolution, then samples it
onto the application's world raster: bilinear terrain, nearest-cell categories
and fields, and world-cell scaling for distances and kinematics. Longitude
wraps; the sampled terrain clamps at polar rows. This ingestion does not turn
the solver into spherical tectonics. Downstream terrain stages can change the
final application terrain, so compare native snapshots when checking the engine.

```powershell
cmake --preset headless
cmake --build out/build/headless --config Release --target world_pipeline tectonics_tests minerals_tests --parallel 4
ctest --test-dir out/build/headless -C Release -L 'tectonics|geology' --output-on-failure
uv run --offline python tests/scripts/check_tectonics_cli.py out/build/headless/bin/Release/world_pipeline.exe
uv run --offline scripts/uw.py tectonics experiment --duration 20 --repeat --explicit-recipe
uv run --offline scripts/uw.py tectonics run
```

`world_pipeline` is the only headless executable: recipe runs and full-world generation (`run`, `generate`). It runs recipes or parts of them (`run --recipe`, `--from-bundle`, `--steps`, `--only`, checkpoints), can write a run bundle (`--bundle`) and audit each step's field access (`--audit-access`), and links no SFML or ImGui. Climate alone on imported terrain is `run --recipe earth-climate`. The desktop `UndiscoveredWorlds.exe` takes no arguments.

Tectonics alone is the `tectonics` recipe: `world_pipeline run --recipe tectonics [--seed N] [--resolution W] [--<scenario-key> value ...] [--duration Myr] [--mineral-history on|off] [--checkpoint-every Myr] --bundle DIR`. It runs the engine on the world grid W x W/2 (default 256) and writes an ordinary run bundle; `world_pipeline --help` lists the scenario options. Two bundles are compared with `scripts/evaluation/compare_bundles.py A B`. See [tectonics runs](docs/architecture/tectonics/pipeline-tooling.md).

Registered runs (`uw.py tectonics run`) take their ID from `runs/registry/tectonics.json`
and never reuse a run directory. They record the git revision and diff hash, the
recipe, the executable's hash and the native fields; executables and sources are
not copied. Archived sources and binaries of earlier runs are evidence, never
development dependencies or defaults for new runs. For normal development, rebuild
the current tree before generating. `uw.py tectonics report --runs <ids>` compares the evaluations of registered runs
(maps and APNG animations in the `configs/tectonics-palettes.json` palette).

Gold history consumes the native column displacement and deterministic topology
identities through `tectonic_frame_reader`. Its point tracers approximate
material movement; they do not reproduce conservative subcell mixing. Formation
and denudation remain provisional resource models. Denudation reads the current
engine elevation and has no private terrain-lowering model; direct per-cell
erosion and climate coupling remain future work. Final drainage routes liberated
gold into placers. Wood, livestock and a full biology module are not implemented
by the tectonics cleanup.

## Development and remaining integration

The project uses C++ with CMake and vcpkg, SFML, Dear ImGui, ImGui-SFML, and ImGuiFileDialog. CUDA rendering is optional in the existing build configuration. Tectonics is built from `src/simulations/geology/plate_tectonics` as the internal `plate_tectonics` target; no standalone checkout is required. Every target is built as C++17. See the [tectonics integration and checks](docs/architecture/tectonic-contract-validation.md) and [native output contract](docs/architecture/tectonics/tectonic-output-contract.md).

The active JSON registry still supplies benchmark IDs; replacing it requires preserving those IDs. Compatible benchmark imports are generated from the retained Earth master. Later refactoring can extract `src/world/` responsibilities.

```powershell
cmake --preset desktop
cmake --build out/build/desktop --config Release --parallel 4
ctest --test-dir out/build/desktop -C Release --output-on-failure -j 4
uv run --offline python -m unittest discover -s tests/scripts
uv run --offline python scripts/check_layout.py
```

Tests mirror the code domains. Prefer deterministic headless checks of state, budgets, and data contracts before visual comparisons. Fixed seeds, explicit inputs, and recorded units/resolutions should accompany benchmark results.

Climate generation accepts `--climate-resolution` (default 128) independently of `--resolution`, the world width that also sets the size of every exported map. The default sets atmosphere/ocean width to 128, moisture/energy width to 256, and weather width to 64, capped by the world grid. Use multiples of four, at least eight. `--climate-temperature-mode calibrated` retains the requested global mean; `radiative` fixes the longwave intercept so the mean emerges from the energy balance. The default uses conservative latitude diffusion, continuous land/soil and ocean/ice heat stores, and monthly Köppen classification. Climate working state belongs to each world; monthly fields and spin-up reservoirs are transient and are not added to the save format.

There are two configure presets: `desktop` (everything, in `out/build/desktop`) and `headless` (`out/build/headless`: no SFML/ImGui or CUDA discovery, with the tests). The configuration is chosen at build time: `--config Release` for every simulation and benchmark, `--config Debug` only to debug. Binaries land in `bin/<Config>/`; scripts take the headless Release paths from `scripts/uwtools/paths.py`.

```powershell
cmake --preset headless
cmake --build out/build/headless --config Release --target world_pipeline --parallel 4
```

`BUILD_TESTING` builds three GoogleTest v1.15.2 binaries (downloaded on first configuration): `world_tests` (core, I/O, pipeline, climate, validation), `tectonics_tests` and `minerals_tests`; ctest lists every `TEST` as `Suite.Name`. The `headless` preset builds `world_pipeline` and the tests. PNG/APNG support uses the PNG/ZLIB packages already supplied with SFML. The public CLI regression command above checks the shared defaults and rejects retired controls.

Git tracks application code, tooling, configuration, small runtime assets and deterministic test fixtures. `docs/`, `runs/registry/` and `runs/ledger/` are tracked; `meta/`, `refs/` and the rest of `runs/` are local-only. Historical `extra/`, research copies, downloads and generated artifacts are also excluded. Local links into these directories require the local workspace data. Keep preparation and reporting code in `scripts/`; the Earth-map entry point is `uv run scripts/refs/prepare_earth_maps.py`.

Enable the repository guards with `git config core.hooksPath .githooks`. They reject local-only paths and files larger than 5 MiB before commits and inspect reachable history before pushes. Run `uv run --offline python scripts/repository_scope.py --staged` or `--history HEAD` manually. After the data-removal history rewrite, use a fresh clone or rebase work onto the cleaned history; merging an old branch would reintroduce the discarded data. The run registry (`runs/registry/`) supplies IDs and is tracked.

## Using the application

Create a seeded world, load a compatible saved world, or import maps through World controls. Click the global map to inspect a location; Zoom opens its region, whose minimap and arrow keys navigate nearby areas. Export options support world, region, and selected-area maps.

Appearance edits affect rendering, not simulated fields. Gradient controls support continuous or discrete colours and editable value/colour anchors. Load or save `.uws` presets through the appearance panel; worlds also save their appearance settings.

Map imports must match the current world dimensions. Numeric imports accept supported grayscale uint16 TIFFs; PNG imports use field-specific channel scales or a one-pixel gradient strip with an explicit minimum and increment. The signed float32 Earth master is converted into the unsigned land/sea import pair during reference processing. See `src/io/map_imports.cpp` for the current decoding rules.

Run `uv run --offline scripts/uw.py climate run` (1024×512 world, 128-wide atmosphere/ocean, 256-wide moisture) for a registered run with its bundle, evaluation and run page in `runs/climate/<id>/`. `uw.py` is the common entry point for tectonics, geology, climate and full-world runs; see `docs/reference-guides/benchmark-commands.md` and the [reference map/palette guide](refs/processed/climate/maps/README.md). Builds and scratch output stay in `out/`; the C++ `climate-benchmark` command writes its maps and run log to `out/climate-benchmark/`.

The [experimental tropical boundary-layer replacement](runs/archive/maps/climate/README.md#experimental-tropical-boundary-layer) is selectable with `--climate-tropical-closure mixed-layer`. The default remains `legacy`: the new wind balance improves tropical wind error but currently regresses rainfall, so it has not passed promotion checks.

The simulator, Earth benchmark imports and exports, and every processed reference use cell-centred `W x W/2` grids. Raw sources keep their native registration until the preparation script performs the ingestion-boundary remap. The default benchmark selection is Koppen, January surface LIC and particles, and precipitation PNG/TIFF. Repeat `--map` (or use a comma-separated list) to select outputs; `all` and `none` are supported. [Map guides](docs/reference-guides/) explain quantity, units, scales, and seasonal meaning.

Generate compatible topography and climate reference maps for one or more widths with `uv run --offline python scripts/refs/prepare_reduced_earth_benchmark.py --width 128 256 512 1024 2048`. PNG previews, float32 GeoTIFF fields, and CSVs go to category-specific `refs/processed/<category>/maps/`, `fields/`, and `csv/` folders. Products include the binary land/ocean mask, rain, temperature, observed Köppen classes, seasonal winds/LIC/particles, and available ERA5 diagnostics. The manifest reports missing reference fields; see [reference preparation](refs/README.md) for coverage, units, and moisture-flux input requirements. Existing uint16 land/sea import TIFFs remain in `maps/`; `configs/app.env` names their base paths, and the executable selects the matching `_WIDTHxHEIGHT` variant.

## Credits and license

Undiscovered Worlds was written by Jonathan Hill, with additional contributions and corrections by Frank Gennari. Sources of adapted code are credited in the code comments. The project is distributed under the [GNU General Public License v3](LICENSE).

The internal plate-tectonics module retains its original copyright notices and [upstream license](src/simulations/geology/plate_tectonics/LICENSE).
