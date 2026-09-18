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

The repository contains terrain generation and import support, an internal plate-tectonics module and world adapter, a generation workbench and stage registry, a reduced climate simulator, physical-resource heuristics, and initial settlement/trade/polity generation. Climate is the current development focus. Terrain refinement remains work in progress; the full ecological and staged civilizational lifecycle above is a design goal.

The September 2026 cleanup aligns CMake sources, includes, configuration, and preparation/analysis scripts with the new layout. Historical and new benchmark maps use category/season/layer paths under `runs/maps/climate/`, with registered run IDs as filenames. Season folders are 1=January, 2=April, 3=July and 4=October. Earth benchmark land/sea imports can be regenerated at arbitrary even resolutions from the retained signed Earth master.

## Repository layout

```text
src/
  core/                    Grids, seasonal indices, deterministic seeds, row workers
  pipeline/                Stage registry, climate and physical/social orchestration
  simulations/
    geology/               Plate-tectonics engine and world adapter
    terrain/               Fractals, coasts, ridges, elevation, FastLEM, landforms
    climate/               Atmosphere, ocean, energy, moisture, transport
      legacy/              Tested alternative thermal and ocean feedback closures
    hydrology/             Drainage, basins, lakes, deltas, wetlands, tides
    resources/             Mineral and marine resource potentials
    society/               Settlement, routes, trade, polities, history
  io/                      Map imports, raster/reference readers, export selection
  validation/climate/      Reference comparisons and climate diagnostics
  app/                     Desktop entry point, windows, controls
    rendering/             Map appearance and CPU/CUDA rendering
    diagnostics/           Profiling and generation debug integration
  wip/                     Mixed responsibilities awaiting extraction
assets/                    Presentation files: textures and appearance presets
definitions/society/       Reserved for authored society content; currently empty
configs/                   Application configuration and saved UI layout
refs/
  source/                  Original downloads, receipts, and licenses
  processed/               Prepared grids grouped by physical quantity
  prepare.py               Reusable Earth-map preparation; configurable width
runs/                      Simulation artifacts, organized by artifact category
tests/                     Core, pipeline, I/O, simulation, and script tests
scripts/
  refs/                    Reference download and preparation
  benchmarks/              Benchmark analysis and workbook input preparation
  archive/                 Previous analysis scripts and tuning experiments
docs/                      Architecture, reference guides, and current climate notes
references/                Local research images, papers, and reference projects
0/                         Files awaiting manual classification
tools/tectonics/           Simulation, snapshot/export tools, heightmap statistics
vcpkg-ports/               Dependency overlay ports
out/                       Build products, editor caches, scratch work
```

Root files are limited to the README, license, version, Git ignore rules, and CMake/vcpkg entry points. Git and Codex metadata retain their operational locations.

### Code boundaries

Keep headers beside their implementations. Organize by the responsibility of a module; express execution order in pipeline recipes. Avoid creating a separate climate implementation for each place it appears in a recipe.

- `core/` should contain shared state primitives, coordinates, units, grids, and deterministic randomness. Extraction of world state into this boundary is still pending.
- `simulations/` should operate on explicit world inputs and produce defined world layers. Presentation, filesystem operations, and reference-data validation have their own boundaries.
- `pipeline/` should select stages, validate their dependencies, schedule them, and control checkpoints and stopping conditions. A disabled stage requires either valid existing outputs, an imported replacement, or an explicit fallback.
- `io/` contains code that reads and writes artifacts. `refs/` and `runs/` contain the artifacts themselves.
- `validation/` compares simulated results with observations and checks budgets. Observed climate references are validation inputs, not hidden forcing for procedural worlds.
- `app/` presents and controls the engine. Both interactive and eventual headless workflows should use the same simulation stages.

`globalterrain`, `globalclimate`, and `physical_layers` have been split into simulation modules and pipeline entry points. Society CSV loading belongs to `io/social_definitions`; its stage order belongs to `pipeline/social_generation`. [Module boundaries and retention decisions](docs/architecture/simulation-modules.md) describe the active and legacy paths.

`src/wip/` still holds `planet`, `region`, `regionalmap`, `classes`, `misc`, `functions`, `generation_tuning`, and `generation_workbench`. World storage/persistence, regional generation, shared utilities, and workbench UI remain coupled there. `functions.hpp` is a compatibility umbrella; new callers should include the owning module's header.

Further land-surface hydrology and ecology belong under `src/simulations/hydrology/` and `src/simulations/ecology/`. Existing climate moisture/storage physics remains in `climate/`. Current mineral and marine scores are derived heuristics in `resources/`; future mineral formation belongs in geology and biological production in ecology, with society consuming both.

### Assets, definitions, configurations, and references

| Location | Role | Examples |
| --- | --- | --- |
| `assets/` | Presentation | Sphere texture, map appearance presets |
| `definitions/` | Authored simulation content | Commodity properties, knowledge definitions, species preferences |
| `configs/` | How an application or run is configured | Application paths; future recipes, resolution, duration, stage switches |
| `refs/` | Imported inputs and external reference data | Heightmaps, observed precipitation, prepared Earth grids |

The distinction follows purpose. An imported heightmap is a dataset; a generated heightmap is a run artifact. Current tuning constants remain in `src/wip/generation_tuning.hpp`; moving them into configurable recipes is future work.

Keep downloads and their receipts/licenses under `refs/source/<dataset>/`. Prepared grids belong under `refs/processed/<quantity>/`; Earth rasters use `{width}` filenames and metadata under `processed/metadata/{width}/`. See [reference preparation](refs/README.md). Research material and reference projects live in `references/`; unresolved files belong in `0/`.

## Runs and comparisons

Artifact category comes before run identity. Map filenames contain only the run ID and extension:

```text
runs/
  maps/
    climate/koppen/
      145.png
      146.png
    climate/wind/lic/1/s/
      145.png
      146.png
    climate/rain/annual/s/
      145.png
      146.png
  fields/climate/<quantity>/<season>/<layer>/<run_id>.tif
  diagnostics/climate/<run_id>/   # new durable diagnostic batches
  diagnostics/<topic>/<run_id>.csv # historical diagnostics
  parameters/<run_id>.txt
  manifests/<manifest_type>/<run_id>.<ext>
  metrics/
    climate.xlsx
    metrics.xlsx
    profiling.xlsx
    run_history.xlsx
  registry/climate.json
  checkpoints/
  logs/
  reports/
```

Open a category to compare successive runs. Seasons `1/2/3/4` mean January/April/July/October; `s/u` mean surface/upper winds. Raw fields mirror map folders and retain their numerical values and units. Annual temperature and rain use `annual/s/`; Koppen omits season/layer folders. `climate/air_temp/annual/s/` contains land surface-air temperature, not land skin temperature.

Unregistered map experiments retain their descriptive identifiers under `runs/work/legacy-maps/`.

**Excel workbooks hold comparable metrics.** `runs/metrics/run_history.xlsx` preserves 304 records from both historical registries and 89 parameter files, including the distinct verification run 147. The main `runs/registry/climate.json` remains required by the executable's run-ID allocator and writer. Numerical CSV diagnostics remain reusable calculation inputs. Reports summarize findings; logs retain verbatim execution output. See the [run 159 report](runs/reports/159.md) for the reporting format.

For future runs, record the seed, recipe and stage parameters, code revision, input versions, grid geometry, simulation duration, convergence criteria, and output provenance. Preserve numerical fields as well as previews when they are needed to reproduce a comparison. Changes in palettes, units, grids, or renderer versions should remain visible in the run metadata.

## Development and remaining integration

The project uses C++ with CMake and vcpkg, SFML, Dear ImGui, ImGui-SFML, and ImGuiFileDialog. CUDA rendering is optional in the existing build configuration. Tectonics is built from `src/simulations/geology/plate_tectonics` as the internal `plate_tectonics` target; no standalone checkout is required. The engine uses C++20 and the app retains C++17. See the [tectonics integration and checks](docs/architecture/tectonic-contract-validation.md) and [native output contract](docs/architecture/tectonics/tectonic-output-contract.md).

The active JSON registry still supplies benchmark IDs; replacing it requires preserving those IDs. Compatible benchmark imports are generated from the retained Earth master. Later refactoring can extract `src/wip/` responsibilities.

```powershell
cmake --preset x64-debug
cmake --build out/build/x64-Debug --config Debug --parallel 4
ctest --test-dir out/build/x64-Debug -C Debug --output-on-failure -j 4
uv run --offline --with numpy --with pillow python -m unittest discover -s tests/scripts
uv run --offline python scripts/check-layout.py
```

Tests mirror the code domains. Prefer deterministic headless checks of state, budgets, and data contracts before visual comparisons. Fixed seeds, explicit inputs, and recorded units/resolutions should accompany benchmark results.

Climate generation accepts `--climate-resolution 64` independently of `--resolution` (the exported world width). This sets atmosphere/ocean width to 64, moisture/energy width to 128, and weather width to 32, capped by the world grid. Use multiples of four, at least eight. `--climate-temperature-mode calibrated` retains the requested global mean; `radiative` fixes the longwave intercept so the mean emerges from the energy balance. The default uses conservative latitude diffusion, continuous land/soil and ocean/ice heat stores, and monthly Köppen classification. Climate working state belongs to each world; monthly fields and spin-up reservoirs are transient and are not added to the save format.

The numerical climate, grid, reference I/O, and stage-registry libraries also build with `-DUW_BUILD_APP=OFF`; this skips SFML/ImGui and CUDA discovery. Tests link these same libraries. For a library-only build using the existing toolchain:

```powershell
cmake --preset x64-debug -B out/build/headless-modules -DUW_BUILD_APP=OFF -DBUILD_TESTING=OFF
cmake --build out/build/headless-modules --config Debug --target uw_climate uw_climate_io uw_generation_stages --parallel 4
```

`BUILD_TESTING` also builds the tectonics regression suite using GoogleTest v1.15.2 (downloaded on first configuration). Enable `UW_BUILD_TECTONICS_TOOLS` to build `tectonic_simulation`, `contract_fixture`, and `tectonic_pipeline`. Tool PNG/APNG support uses the PNG/ZLIB packages already supplied with SFML. The optional `check-tectonics-slow` target runs the full-size deterministic replay and generates comparison images.

Git tracks application code, tooling, configuration, small runtime assets and deterministic test fixtures. `meta/`, `docs/`, `refs/` and `runs/` are entirely local-only, including guides, registries and workbooks. Historical `extra/`, research copies, downloads and generated artifacts are also excluded. Local links into these directories require the local workspace data. Keep preparation and reporting code in `scripts/`; the Earth-map entry point is `uv run scripts/refs/prepare-earth-maps.py`.

Enable the repository guards with `git config core.hooksPath .githooks`. They reject local-only paths and files larger than 5 MiB before commits and inspect reachable history before pushes. Run `uv run --offline python scripts/repository_scope.py --staged` or `--history HEAD` manually. After the data-removal history rewrite, use a fresh clone or rebase work onto the cleaned history; merging an old branch would reintroduce the discarded data. The local run registry still supplies IDs and is preserved outside version control.

## Using the application

Create a seeded world, load a compatible saved world, or import maps through World controls. Click the global map to inspect a location; Zoom opens its region, whose minimap and arrow keys navigate nearby areas. Export options support world, region, and selected-area maps.

Appearance edits affect rendering, not simulated fields. Gradient controls support continuous or discrete colours and editable value/colour anchors. Load or save `.uws` presets through the appearance panel; worlds also save their appearance settings.

Map imports must match the current world dimensions. Numeric imports accept supported grayscale uint16 TIFFs; PNG imports use field-specific channel scales or a one-pixel gradient strip with an explicit minimum and increment. The signed float32 Earth master is converted into the unsigned land/sea import pair during reference processing. See `src/io/map_imports.cpp` for the current decoding rules.

Run `./scripts/benchmarks/run-climate-benchmark.ps1 -Resolution 512 -ClimateResolution 64 -Seed 20260906` for a complete benchmark with logs and a side-by-side map gallery. See [benchmark commands and outputs](runs/maps/climate/README.md) and the [reference map/palette guide](refs/processed/climate/maps/README.md). Builds stay in `out/`; benchmark artifacts default to `runs/`.

The [experimental tropical boundary-layer replacement](runs/maps/climate/README.md#experimental-tropical-boundary-layer) is selectable with `-TropicalClosure mixed-layer`. The default remains `legacy`: the new wind balance improves tropical wind error but currently regresses rainfall, so it has not passed promotion checks.

The simulator, Earth benchmark imports and exports, and every processed reference use cell-centred `W x W/2` grids. Raw sources keep their native registration until the preparation script performs the ingestion-boundary remap. The default benchmark selection is Koppen, January surface LIC and particles, and precipitation PNG/TIFF. Repeat `--map` (or use a comma-separated list) to select outputs; `all` and `none` are supported. [Map guides](docs/reference-guides/) explain quantity, units, scales, and seasonal meaning. The summary command `uv run --offline --with numpy --with pillow python scripts/benchmarks/summarize-climate-validation.py --run-id ID --output DIR` prepares workbook input tables from category-organized diagnostics and carries the latest earlier report as context. It does not append to `metrics.xlsx` automatically.

Generate compatible topography and climate reference maps for one or more widths with `uv run --offline --with numpy --with pillow python scripts/refs/prepare-reduced-earth-benchmark.py --width 128 256 512 1024 2048`. PNG previews, float32 GeoTIFF fields, and CSVs go to category-specific `refs/processed/<category>/maps/`, `fields/`, and `csv/` folders. Products include the binary land/ocean mask, rain, temperature, observed Köppen classes, seasonal winds/LIC/particles, and available ERA5 diagnostics. The manifest reports missing reference fields; see [reference preparation](refs/README.md) for coverage, units, and moisture-flux input requirements. Existing uint16 land/sea import TIFFs remain in `maps/`; `configs/app.env` names their base paths, and the executable selects the matching `_WIDTHxHEIGHT` variant.

## Credits and license

Undiscovered Worlds was written by Jonathan Hill, with additional contributions and corrections by Frank Gennari. Sources of adapted code are credited in the code comments. The project is distributed under the [GNU General Public License v3](LICENSE).

The internal plate-tectonics module retains its original copyright notices and [upstream license](src/simulations/geology/plate_tectonics/LICENSE).
