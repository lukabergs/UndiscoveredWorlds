# Benchmark and run commands

Run everything from the repository root in PowerShell (VS Code terminal or Visual Studio **View → Terminal**).
Every benchmark is headless. The desktop app `UndiscoveredWorlds.exe` takes no arguments and is only for interactive viewing.

## One syntax for every module

```powershell
uv run --offline scripts/benchmarks/bench.py <module> <command> [shared options] [native options]
```

| Module | Commands | Native tool |
| --- | --- | --- |
| `tectonics` | `run` (registered), `inspect` (scratch), `report`, `experiment`, `evaluate` | `world_pipeline.exe run --recipe tectonics` |
| `geology` | `run` | `world_pipeline.exe run --recipe tectonics` + `benchmark-geology.py` |
| `climate` | `run` (registered), `check`, `report`, `experiment`, `evaluate` | `world_pipeline.exe climate-benchmark`; `world_pipeline.exe run --recipe earth-climate` for experiments |
| `world` | `run`, `experiment`, `evaluate` | `world_pipeline.exe generate` |

Shared options, same meaning everywhere:

| Option | Meaning | Default |
| --- | --- | --- |
| `--seed N` | simulation seed | 1 |
| `--resolution W` | world grid width; height is W/2. Every exported map has this size | tectonics 256, geology 128 and 256, climate 1024, world app default |
| `--executable PATH` | native tool to run | see Builds |
| `--output PATH` | output folder or file (scratch commands) | `out/bench/...` |
| `--label TEXT` | change summary stored with a climate run | none |
| `--runs ID [ID...]` | run IDs for `report`/`check` | none |
| `--compare ID` | climate run to compare against | none |

Any other `--flag value` goes to the native tool unchanged. Every module's native tool is `world_pipeline`. Tectonics takes the scenario options of `run --recipe tectonics` (`--plates 14`, `--duration 300`, ...); climate and world take `--climate-resolution 128`, `--climate-tropical-closure mixed-layer`, ... List them with the native tool's help:

```powershell
out/build/headless/bin/Release/world_pipeline.exe --help
```

## What every registered run produces

`bench.py <module> run` and `bench.py <module> report --runs ...` produce the same set of outputs for each module:

| Output | Climate | Tectonics |
| --- | --- | --- |
| Raw log (standard header: module, UTC start, executable hash, git state, command) | `runs/logs/climate/<id>.log` | `runs/logs/tectonics/<id>.log` |
| Maps (world-grid PNGs) | `runs/maps/climate/<quantity>/.../<id>.png` | `runs/maps/tectonics/<field>/<checkpoint>/<id>.png` |
| Projected views (Mollweide, north/south polar) | `runs/maps/climate/projections/<projection>/...` | `runs/maps/tectonics/projections/<projection>/...` |
| Animations | not yet | `runs/maps/tectonics/<field>/animation/<id>.png` (APNG) |
| Numerical fields | `runs/fields/climate/...` (GeoTIFF) | `runs/fields/tectonics/...` (`.npy`) |
| Run page (standard template, map + projection selector, side by side with the previous run) | `runs/reports/<id>.html` | `runs/reports/tectonics/<id>.html` |
| Record / findings (Markdown) | `runs/reports/<id>.md` | `runs/reports/tectonics/<id>.md` |
| Detailed comparison | `runs/reports/<id>-gallery.html` (against observations) | `runs/reports/tectonics/<previous>-<id>.html` |

The run page is built by `scripts/benchmarks/run_report.py` from `run_report_template.html`; projections come from `map_projection.py` and colouring from `map_rendering.py`, shared by all modules. To rebuild a page: `uv run --offline scripts/benchmarks/run_report.py <module> --run <id> [--compare <id>]`.

## Default seed

**Seed 1** is the default for every new run. Use another seed only for a deliberate seed comparison. Tectonics `run` with a seed other than 1 is recorded as a seed comparison and does not replace the tectonics baseline.

Climate runs registered up to run 365 used seed **20260906**. Use it only to reproduce or extend those runs; adaptive batches set up before this change keep 20260906 automatically. The 20260906 inside the report renderers is the texture seed for wind and current images, shared by model and reference, not a simulation seed.

## Builds

Rebuild after every source change; the tools never pick up edits on their own. Configure only the first time, or after changing a `CMakeLists.txt`, a preset, `configs/tectonics-baseline.json`, or the Visual Studio toolset.

```powershell
# Tectonics, geology, climate and world (world_pipeline): headless build (no SFML, ImGui or CUDA); bench.py uses it by default
cmake --preset headless
cmake --build out/build/headless --config Release --target world_pipeline --parallel 4

# Desktop app
cmake --preset desktop
cmake --build out/build/desktop --config Release --target UndiscoveredWorlds --parallel 4
```

Binaries land in `bin/<Config>/` (`bin/Release`, `bin/Debug`). Simulations and benchmarks always use Release; scripts take the paths from `scripts/uwtools/paths.py`.

`world_pipeline` links the headless world core (`uw_world_core`) and the simulation modules; it depends on libpng but not on SFML or ImGui. `bench.py` defaults to `out/build/headless/bin/Release/world_pipeline.exe`; pass `--executable` for another build.

## Recipe runs: `world_pipeline run`

`run` executes a named recipe, or part of one, and registers nothing. A recipe on
imported terrain (`earth-climate`, `imported-world`) imports the Earth terrain of
`configs/app.env` for the width unless `--import-land/--import-sea` or
`--from-bundle` are given.

```powershell
# Climate on Earth terrain, Geosphere frozen (default 512 wide).
world_pipeline.exe run --recipe earth-climate --seed 1 --resolution 1024 --climate-resolution 128 --bundle out/bench/climate-1024
# Experimental: extra diagnostic levels from one rule set, stored as generic-<p>hPa levels in the bundle.
world_pipeline.exe run --recipe earth-climate --seed 1 --climate-generic-levels 925,850,700 --bundle out/bench/climate-levels
# A full world with a checkpoint bundle after the last terrain step -> out/bench/world-checkpoints/17-terrain.roughness
world_pipeline.exe run --recipe procedural-world --seed 1 --checkpoint-at terrain.roughness --bundle out/bench/world
# One part of the recipe on that checkpoint: the climate steps only.
world_pipeline.exe run --recipe procedural-world --from-bundle out/bench/world-checkpoints/17-terrain.roughness --steps climate.temperature..climate.rainfall_and_coupling --bundle out/bench/world-climate
```

`--steps a..b` and `--only <process>` take step indices or process ids (the `steps`
of any bundle manifest list both); `--checkpoint-at`, `--checkpoint-every N` and
`--checkpoint-dir` write one bundle per chosen step. `--from-bundle` loads the
world of a bundle: its seed and grid apply, command-line options (`--climate-*`,
`--no-lakes`, ...) are not inherited and must be repeated. A bundle carries the
world's stored fields. It does not carry tables, the climate solver's own grids,
levels and diagnostics, or state kept between steps outside the world fields
(mountain drainage, the salt-lake masks). The scheduler therefore refuses a start
point whose steps read such state from a skipped step and names the step that
must be in the run; those fields and tables are copied into the new bundle from
the source unless a step of the run writes them.

`world_pipeline climate-benchmark` is the registered form of the same run: the
earth-climate recipe on the Earth terrain, recorded with a run ID. Until run 384
it ran the imported-world stage path, which edited the imported terrain (mountain
bases, smoothing, depression fill, coastline adjustment) before climate, so runs
up to 384 are not comparable with later ones; `--no-rivers`, `--no-lakes` and
`--no-deltas` no longer have any effect on it.

## Recipes, run bundles and access audits

`world_pipeline generate` also runs a recipe: `procedural-world` by default, `imported-world` with `--import-land/--import-sea`, and `earth-climate` with `--climate-only` (climate processes only; Geosphere frozen). `--recipe NAME` selects one explicitly. Recipes and the process registry are described in `docs/architecture/world-model.md`.

```powershell
# Write every world field and table at full precision (NumPy files + manifest.json).
world_pipeline.exe generate --seed 1 --resolution 512 --bundle out/bench/bundle-a
# Check that each step writes only its declared fields (hashes every field around each step).
world_pipeline.exe generate --seed 1 --resolution 512 --audit-access --access-report out/bench/access.json
# Compare two bundles field by field; exit status 1 on any difference.
python scripts/evaluation/compare_bundles.py out/bench/bundle-a out/bench/bundle-b
# Against a bundle from before run 385 (climate fields as int16): compare in the former encoding.
python scripts/evaluation/compare_bundles.py out/bench/old-bundle out/bench/bundle-b --legacy-encoding
```

The bundle comparison is the equivalence check for refactors: unchanged code must give identical fields for the same seed and options. Access reports also hold the hash of every written field after each step (`write_hashes`); comparing the reports of two runs with the same options finds the first step where they diverge. A separate analysis build with `-DCMAKE_CXX_FLAGS_RELWITHDEBINFO="/O2 /Ob1 /Zi /DNDEBUG /DUW_CHECK_GRID_BOUNDS"` (for example in `out/build/grid-check`, `--config RelWithDebInfo`) aborts with a stack trace at the first world-grid index outside the grid. Declared process access is generated from access reports of an analysis build configured with `-DUW_TRACK_FIELD_READS=ON` (`scripts/evaluation/process_access.py`).

If configure fails with "is not a full path to an existing compiler tool", Visual Studio updated its compiler. Add `--fresh` to the configure command once.

## Evaluation: checks, validation and review

`scripts/evaluation/` evaluates run bundles with a per-module profile in `configs/evaluation/<module>.json` (`climate`, `world`). A profile lists checks (grid registration, required and written fields, finite values, ranges, frozen components, repeat determinism), validation comparisons with observed references, review maps in the project palettes, and how the module weights the three kinds of evidence. A new module needs a profile, not new scripts.

```powershell
# Scratch experiment: run the profile's tool (seed 1), write bundle/, then evaluation/ → out/evaluation/<module>/<name>/
uv run --offline scripts/benchmarks/bench.py world experiment --resolution 256 --repeat
uv run --offline scripts/benchmarks/bench.py climate experiment --climate-resolution 64
# Evaluate an existing bundle, optionally against a baseline bundle on common masks.
uv run --offline scripts/benchmarks/bench.py climate evaluate --bundle out/bench/climate-1024 --baseline out/bench/climate-1024-before
# Record the user's verdict with the evaluation (approved, rejected or revise).
uv run --offline scripts/evaluation/review.py out/evaluation/world/<name>/evaluation --verdict approved --notes "coastlines"
```

Each evaluation directory holds `evaluation.json` (inputs with hashes, check results, validation metrics, verdict), `report.md`, `index.html` and `maps/`; it is never overwritten. `experiment.py --archive runs/reports/<experiment>` keeps the runs in an experiment archive routed to the data drive. Tool outputs that would normally go under `runs/` are redirected into the run directory.

Climate comparisons reuse the climate report's per-field metrics, so every field loss is defined as in global-balanced-v1. Bundles do not yet contain column water, moisture transport, moisture convergence or currents in m/s, so a bundle evaluation reports a partial score over the remaining 70 % of the weight with its coverage; it is not a global-balanced-v1 value and does not replace the registered climate ranking. Air temperature against ERA5 2 m temperature is reported with zero weight. The `world` profile reports Earth hypsometry as context for review, not as a target.

## Tectonics

```powershell
# Scratch: final-state run bundle of the baseline recipe (seed 1, 256×128, 10 plates, 600 Myr) → out/bench/tectonics-seed1-256
uv run --offline scripts/benchmarks/bench.py tectonics inspect
uv run --offline scripts/benchmarks/bench.py tectonics inspect --resolution 512 --plates 14 --output out/bench/p14

# Registered run: next ID in runs/registry/tectonics.json; checkpoints in runs/checkpoints/tectonics/<id>/snapshots/, final bundle in runs/checkpoints/tectonics/<id>/final
uv run --offline scripts/benchmarks/bench.py tectonics run
uv run --offline scripts/benchmarks/bench.py tectonics run --plates 14 --duration 300

# Maps, animations, comparison report and run page for registered runs (last ID is the run page)
uv run --offline scripts/benchmarks/bench.py tectonics report --runs 74 75
uv run --offline scripts/benchmarks/bench.py tectonics report --runs 13 --maps-only
```

`run` copies `world_pipeline.exe` into `runs/reports/tectonics/<id>/` and runs it from there with a checkpoint for the starting Geosphere and one every `checkpoint_every_updates` updates. Scenario flags change only that run; the baseline file is not modified. Registered runs up to 75 keep the former snapshot layout (raw arrays plus a manifest per snapshot and a `batch_manifest.json`); the report scripts read both layouts. `bench.py tectonics experiment|evaluate` run and evaluate a bundle with `configs/evaluation/tectonics.json` (grid, required fields, finite values, ranges, repeat; Earth hypsometry as context; review maps).

Mantle plumes are off by default (`--plume-count 0`). They are scenario flags like any other, so `run`, `inspect` and `report` take them unchanged:

```powershell
uv run --offline scripts/benchmarks/bench.py tectonics inspect --plume-count 6 --plume-tail-flux-km3-yr 0.25 --plume-swell-age-myr 10
uv run --offline scripts/benchmarks/bench.py tectonics run --plume-count 6 --plume-head-volume-km3 2000000 --plume-birth-window-myr 400
```

`world_pipeline.exe --help` lists every `--plume-*` flag. Bundles of runs with plumes carry `plume_influence` (palette entry in `configs/tectonics-palettes.json`); the column budget reports `generated_plume_kg` as a subset of `generated_ocean_kg`.

Native tool, for other durations, checkpoints and mineral history (see `docs/architecture/tectonics/pipeline-tooling.md`):

```powershell
$wp = "out/build/headless/bin/Release/world_pipeline.exe"
# Scratch bundle of 300 Myr
& $wp run --recipe tectonics --seed 1 --duration 300 --bundle out/bench/tectonics-300myr
# Checkpoints every 50 Myr -> out/bench/tectonics-ck-checkpoints/01-geology.plate_tectonics-t0050.000, ...
& $wp run --recipe tectonics --seed 1 --checkpoint-every 50 --bundle out/bench/tectonics-ck
# Mineral history off (on by default; tables are written at the end of the run)
& $wp run --recipe tectonics --seed 1 --mineral-history off --bundle out/bench/tectonics-nomin
# Compare two bundles field by field
uv run --offline scripts/evaluation/compare_bundles.py out/bench/a out/bench/b
```

The climate/erosion feedback experiment (run 63) has no runner at present; it returns as a recipe once `terrain.climate_erosion` is an evolver (`docs/architecture/consolidation-plan.md` §6). `bench.py tectonics report` still renders run 63.

## Geology

Tectonics at each width, scored for seafloor age and collisions against Earth references. Output: `out/bench/geology/report.json`.

```powershell
uv run --offline scripts/benchmarks/bench.py geology run
uv run --offline scripts/benchmarks/bench.py geology run --resolution 512 --reference-width 1024
```

## Climate

Earth terrain, simulated climate, compared with observations. `run` takes the next ID from `runs/registry/climate.json`, writes maps and fields under `runs/`, the log to `runs/logs/climate/<id>.log`, and a gallery report.

```powershell
# Standard benchmark (run 373 grids): seed 1, 1024×512 world, 128-wide atmosphere/ocean, 256-wide moisture, 8 coupling updates
uv run --offline scripts/benchmarks/bench.py climate run --label "what changed"

# Quick smoke run on a small world (about 2 minutes)
uv run --offline scripts/benchmarks/bench.py climate run --resolution 128 --climate-resolution 64 --climate-hydrology-resolution 128

# Model switches
uv run --offline scripts/benchmarks/bench.py climate run --climate-tropical-closure mixed-layer
uv run --offline scripts/benchmarks/bench.py climate run --map "koppen,temp,precip"

# Checks and gallery for an existing run, no new simulation
uv run --offline scripts/benchmarks/bench.py climate check --runs 366 --compare 365
uv run --offline scripts/benchmarks/bench.py climate report --runs 366 --compare 365
```

The world grid (`--resolution`) sets the size of every exported map and field. Atmosphere/ocean (`--climate-resolution`, default 128) and moisture (`--climate-hydrology-resolution`, default 256) run on coarser grids and are remapped onto the world grid, as in run 373.

Land surface drag: climate runs use the Earth vegetation drag contrast (`refs/processed/climate/vegetation-drag-20260914/`, prepared from ERA5 vegetation cover, as in runs up to 378). Other world widths get a version area-averaged from the 1024 file, cached beside it. `--surface-drag none` gives uniform land drag; `--surface-drag PATH` uses another file. The log records which input was used. Procedural worlds never use it.

Benchmark defaults added by `bench.py` unless you pass them: `--resolution 1024 --climate-resolution 128 --climate-hydrology-resolution 256 --climate-coupling-iterations 8`, the standard map list, and no rivers, lakes or deltas. Other `world_pipeline` defaults apply as printed by `--help`.

Multi-run comparison report and adaptive experiments keep their own scripts:

```powershell
uv run --offline scripts/benchmarks/climate_report.py --archive runs/reports/<experiment> --baseline 365 --runs 365 366
uv run --offline scripts/benchmarks/climate_experiment.py setup --archive runs/reports/<experiment> --count 10
uv run --offline scripts/benchmarks/climate_experiment.py run --archive runs/reports/<experiment> --run <id> --plan plan.json
uv run --offline scripts/benchmarks/climate_experiment.py analyze --archive runs/reports/<experiment> --run <id>
uv run --offline scripts/benchmarks/climate_experiment.py decision --archive runs/reports/<experiment> --run <id> --plan decision.json
uv run --offline scripts/benchmarks/climate_experiment.py score --archive runs/reports/<experiment>
```

Each experiment run must be analyzed and given a decision before the next one starts. See `docs/reference-guides/climate-reporting.md` for report controls and `runs/maps/climate/README.md` for workflow details.

## Full world

```powershell
# Saved to out/bench/world-seed1-512.uww
uv run --offline scripts/benchmarks/bench.py world run --resolution 512
uv run --offline scripts/benchmarks/bench.py world run --resolution 512 --plates 12 --tectonic-duration-myr 300
uv run --offline scripts/benchmarks/bench.py world run --resolution 256 --social --social-mode historical --history-years 500
```

`--no-rivers`, `--no-lakes`, `--no-deltas` and `--no-profile-log` shorten runs.
