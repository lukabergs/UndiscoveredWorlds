# Benchmark and run commands

Run everything from the repository root in PowerShell (VS Code terminal or Visual Studio **View → Terminal**).
Every benchmark is headless. The desktop app `UndiscoveredWorlds.exe` takes no arguments and is only for interactive viewing.

## One syntax for every module

```powershell
uv run --offline scripts/uw.py <module> <command> [shared options] [native options]
```

| Module | Commands | Native tool |
| --- | --- | --- |
| `tectonics` | `run` (registered), `experiment` (scratch), `evaluate`, `report`, `dashboard` | `world_pipeline.exe run --recipe tectonics` |
| `geology` | `run` | `world_pipeline.exe run --recipe tectonics` + `benchmark_geology.py` |
| `climate` | `run` (registered), `experiment` (scratch), `evaluate`, `report --evaluation`, `dashboard` | `world_pipeline.exe run --recipe earth-climate` |
| `world` | `run`, `report`, `experiment`, `evaluate`, `dashboard` | `world_pipeline.exe generate` |

Shared options, same meaning everywhere:

| Option | Meaning | Default |
| --- | --- | --- |
| `--seed N` | simulation seed | 1 |
| `--resolution W` | world grid width; height is W/2. Every exported map has this size | tectonics 256, geology 128 and 256, climate 1024, world app default |
| `--executable PATH` | native tool to run | see Builds |
| `--output PATH` | output folder or file (scratch commands) | `out/bench/...` |
| `--label TEXT` | change summary stored with a run (run.json and the registry record) | none |
| `--runs ID [ID...]` | run IDs for `tectonics report` | none |
| `--evaluation DIR [DIR...]` | evaluation directories for `report` pages; the first is the baseline | none |

Any other `--flag value` goes to the native tool unchanged. Every module's native tool is `world_pipeline`. Tectonics takes the scenario options of `run --recipe tectonics` (`--plates 14`, `--duration 300`, ...); climate and world take `--climate-resolution 128`, `--climate-tropical-closure mixed-layer`, ... List them with the native tool's help:

```powershell
out/build/headless/bin/Release/world_pipeline.exe --help
```

## What every registered run produces

`uw.py climate run` and `uw.py tectonics run` are one path (`scripts/evaluation/experiment.py --register`). The settings of a registered run are data in the module's profile (`configs/evaluation/<module>.json`, `run` section): `arguments` and declared `inputs` apply to every run of the profile, `registered` adds what a registered run uses (climate: 1024 columns, `--climate-resolution 128 --climate-hydrology-resolution 256 --climate-coupling-iterations 8`; tectonics: the scenario and checkpoint interval of `configs/tectonics-baseline.json`, a checkpoint at `geology.geosphere_initialization` and every `checkpoint_every_updates` × `delta-time-myr` Myr). Native options come last and override them. The run takes the next ID of `runs/registry/<module>.json` and lives in one directory (a junction into the storage drive per module, `scripts/uwtools/storage.py`):

```
runs/<module>/<id>/
  run.json            command, executable path and SHA-256, configuration, code state (git revision, dirty, diff hash),
                      inputs (path, SHA-256), seed, grid, start, finish, seconds, per-step seconds, exit codes, verdict
  log.txt             the tool's output
  bundle/             final run bundle
  checkpoints/        checkpoint bundles (tectonics)
  evaluation/         evaluation.json, review-request.json, index.html (run page), report.md, maps/;
                      analysis/ holds what the tectonics report and review scripts derive from the run
  ledger-entry.json   draft ledger entry: complete it and append it with `uw.py <module> ledger --add`
```

When the run ends a record `{"id", "status", "run": "runs/<module>/<id>", "date", "label", "score"}` is appended to the registry; earlier records stay as they are. Runs registered before 2026-10-05 are sealed under `runs/archive/` (`runs/archive/README.md`): a path `runs/<kind>/...` in their records is now `runs/archive/<kind>/...`. `scripts/uwtools/runs.py` resolves both layouts; nothing writes into the archive (derived outputs of a sealed run go to `runs/reports/<module>/<id>/analysis/`). Cross-run pages live in `runs/reports/<name>/`.

## Default seed

**Seed 1** is the default for every new run. Use another seed only for a deliberate seed comparison. Tectonics `run` with a seed other than 1 is recorded as a seed comparison and does not replace the tectonics baseline.

Climate runs registered up to run 365 used seed **20260906**. Use it only to reproduce or extend those runs; adaptive batches set up before this change keep 20260906 automatically. The 20260906 inside the report renderers is the texture seed for wind and current images, shared by model and reference, not a simulation seed.

## Builds

Rebuild after every source change; the tools never pick up edits on their own. Configure only the first time, or after changing a `CMakeLists.txt`, a preset, `configs/tectonics-baseline.json`, or the Visual Studio toolset.

```powershell
# Tectonics, geology, climate and world (world_pipeline): headless build (no SFML, ImGui or CUDA); uw.py uses it by default
cmake --preset headless
cmake --build out/build/headless --config Release --target world_pipeline --parallel 4

# Desktop app
cmake --preset desktop
cmake --build out/build/desktop --config Release --target UndiscoveredWorlds --parallel 4
```

Binaries land in `bin/<Config>/` (`bin/Release`, `bin/Debug`). Simulations and benchmarks always use Release; scripts take the paths from `scripts/uwtools/paths.py`.

`world_pipeline` links the headless world core (`uw_world_core`) and the simulation modules; it depends on libpng but not on SFML or ImGui. `uw.py` defaults to `out/build/headless/bin/Release/world_pipeline.exe`; pass `--executable` for another build.

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

`world_pipeline climate-benchmark` runs the same earth-climate recipe on the Earth
terrain (a seed-1 bundle at 256 columns is identical to `run --recipe earth-climate`
with the same options) and records it with the C++ validation tables and maps. It is
a scratch tool, as is the desktop app's "Record climate benchmark": their outputs and
their run log go to `out/climate-benchmark/` (`configs/app.env`); registered runs
are made only by `uw.py climate run`. Until run 384 it ran the imported-world stage
path, which edited the imported terrain before climate, so runs up to 384 are not
comparable with later ones.

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

`scripts/evaluation/` evaluates run bundles with a per-module profile in `configs/evaluation/<module>.json` (`climate`, `world`). A profile lists checks (grid registration, required and written fields, finite values, ranges, frozen components, repeat determinism, recipe equality, values of a JSON table, column and mineral budgets), validation comparisons with observed references, review maps in the project palettes, and how the module weights the three kinds of evidence. A new module needs a profile, not new scripts.

```powershell
# Scratch experiment: run the profile's tool (seed 1), write bundle/, then evaluation/ → out/evaluation/<module>/<name>/
uv run --offline scripts/uw.py world experiment --resolution 256 --repeat
uv run --offline scripts/uw.py climate experiment --climate-resolution 64
# Evaluate an existing bundle, optionally against a baseline bundle on common masks.
uv run --offline scripts/uw.py climate evaluate --bundle out/bench/climate-1024 --baseline out/bench/climate-1024-before
# Tectonics: checkpoint bundles give one APNG per map (default <bundle>/../checkpoints).
uv run --offline scripts/uw.py tectonics evaluate --bundle runs/archive/checkpoints/tectonics/76/final --checkpoints runs/archive/checkpoints/tectonics/76/snapshots --output out/evaluation/tectonics/76
# Record the user's verdict with the evaluation (approved, rejected or revise); the run page is rendered again.
uv run --offline scripts/evaluation/review.py out/evaluation/world/<name>/evaluation --verdict approved --notes "coastlines"
# Pages from JSON (scripts/evaluation/report.py and templates/): one directory rewrites its run page,
# several give a comparison page (first = baseline) in --output (default out/evaluation/<module>/comparison).
uv run --offline scripts/uw.py climate report --evaluation <baseline>/evaluation <candidate>/evaluation --output out/evaluation/climate/comparison
# Module dashboard: registry runs (score, status, links to run pages) and ledger entries.
uv run --offline scripts/uw.py tectonics dashboard --output out/evaluation/tectonics/dashboard
```

Every report page is one of three templates in `scripts/evaluation/templates/` (`run.html`, `comparison.html`, `dashboard.html`) with the page's data embedded as JSON; nobody writes per-run HTML. The run page shows the run facts (from `run.json` when present), the verdict, checks, scores and the per-field validation table (loss, baseline loss, delta, RMSE, bias, correlation), the entries of `review.json` and the maps. A profile map with `"validation": "<id>"` (`configs/evaluation/climate.json`) is drawn from that validation entry with its season rule, unit conversion and remap: model, reference in the same palette and grid, and model minus reference on the score's mask, on the diverging `difference` scale of `configs/climate-palettes.json` (vector fields as speed). Quarter seasons are shown as Jan–Mar, Jul–Sep, monthly ones as January, July. `evaluation.json`'s `maps` lists every image with title, path, palette, season and kind (`model`, `reference`, `difference`, `animation`). `report.md` is generated from the same data. Regions come from `configs/regions.json` (the project's latitude/longitude boxes and the 58 IPCC AR6 WGI reference regions, whose rasterised mask is `refs/processed/climate/ipcc_ar6_regions.png`); `scripts/uwtools/regions.py` gives their masks on any grid. `evaluation.json`'s `validation.regions` holds, for every scalar validation entry and season, the area-weighted bias and RMSE of model minus reference per region on the entry's mask; it is reported on the run page (project boxes also in `report.md`) and never enters a score.

Each evaluation directory holds `evaluation.json` (inputs with hashes, check results, validation metrics, verdict), `report.md`, `index.html` and `maps/`; it is never overwritten (pages are re-rendered from it). `experiment.py --archive runs/reports/<experiment>` keeps scratch runs in an experiment archive routed to the data drive. Files the tool writes on its own (C++ climate validation tables, images) are scratch and go to `out/evaluation/outputs/<module>/<name>/`.

Climate comparisons use the per-field metrics of global-balanced-v1 (`scripts/evaluation/score.py`). A bundle that lacks a weighted field gets a partial score over the covered weight, reported with its coverage; a registered climate run carries every weighted field and gets the full score on bundle fields (ledger `climate-0002` compares it with the legacy chain's values). Air temperature against ERA5 2 m temperature is reported with zero weight. The `world` profile reports Earth hypsometry as context for review, not as a target.

## Tectonics

```powershell
# Registered run: next ID in runs/registry/tectonics.json -> runs/tectonics/<id>/ (checkpoints/, bundle/, evaluation/)
uv run --offline scripts/uw.py tectonics run --label "what changed"
uv run --offline scripts/uw.py tectonics run --plates 14 --duration 300 --label "14 plates, 300 Myr"
# Scratch run + evaluation (compiled defaults, final state only) -> out/evaluation/tectonics/<name>/; --repeat checks determinism,
# --explicit-recipe that the compiled defaults equal configs/tectonics-baseline.json given as options
uv run --offline scripts/uw.py tectonics experiment --duration 20 --repeat --explicit-recipe
# Snapshot statistics (counts per regime, crust class and boundary type; segments, deforming regions, junctions) of a bundle
uv run --offline scripts/evaluation/tectonics_statistics.py runs/tectonics/77/bundle
# Maps, animations and comparison page of registered runs, sealed or not -> runs/reports/tectonics/<A>-<B>.html
uv run --offline scripts/uw.py tectonics report --runs 76 77
uv run --offline scripts/uw.py tectonics report --runs 13 --maps-only
# Reviews of a run's snapshots -> its evaluation/analysis/ (sealed runs: runs/reports/tectonics/<id>/analysis/)
uv run --offline scripts/benchmarks/tectonics_review.py 77
uv run --offline scripts/benchmarks/tectonics_batch_analysis.py 77
uv run --offline scripts/benchmarks/tectonics_boundary_review.py 77
uv run --offline scripts/benchmarks/tectonics_coastal_review.py 77
```

Scenario flags change only that run; the baseline file is not modified. A denser sampling of a run is a new run with another `--checkpoint-every`. Registered runs up to 75 keep the former snapshot layout (raw arrays plus a manifest per snapshot); the report scripts read both layouts. The tectonics profile (`configs/evaluation/tectonics.json`) checks the grid, required fields, finite values, ranges, the snapshot record, the column and mineral budgets, repeat and recipe equality, gives Earth hypsometry as context and renders the review maps with one animation per map from the checkpoints.

Mantle plumes are off by default (`--plume-count 0`). They are scenario flags like any other, so `run` and `experiment` take them unchanged:

```powershell
uv run --offline scripts/uw.py tectonics experiment --plume-count 6 --plume-tail-flux-km3-yr 0.25 --plume-swell-age-myr 10
uv run --offline scripts/uw.py tectonics run --plume-count 6 --plume-head-volume-km3 2000000 --plume-birth-window-myr 400
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

The climate/erosion feedback experiment (run 63) has no runner at present; it returns as a recipe once `terrain.climate_erosion` is an evolver (`docs/architecture/consolidation-plan.md` §6). `uw.py tectonics report` still renders run 63.

## Geology

Tectonics at each width, scored for seafloor age and collisions against Earth references. Output: `out/bench/geology/report.json`.

```powershell
uv run --offline scripts/uw.py geology run
uv run --offline scripts/uw.py geology run --resolution 512 --reference-width 1024
```

## Climate

Earth terrain, simulated climate, compared with observations. `run` takes the next ID from `runs/registry/climate.json` and writes `runs/climate/<id>/` (see above); its evaluation gives the global-balanced-v1 score, stored in the registry record. The registry's `baseline` names the run the module's experiments are compared with.

```powershell
# Registered run: seed 1, 1024x512 world, 128-wide atmosphere/ocean, 256-wide moisture, 8 coupling updates (about 8 minutes)
uv run --offline scripts/uw.py climate run --label "what changed"
# Scratch run with the profile's experiment settings (512 wide) or smaller -> out/evaluation/climate/<name>/
uv run --offline scripts/uw.py climate experiment --resolution 128 --climate-resolution 64 --climate-hydrology-resolution 128
# Model switches
uv run --offline scripts/uw.py climate run --climate-tropical-closure mixed-layer --label "mixed-layer closure"
# Comparison page of two evaluations (the first is the baseline)
uv run --offline scripts/uw.py climate report --evaluation runs/climate/<baseline>/evaluation runs/climate/<id>/evaluation --output runs/reports/<name>
```

The world grid (`--resolution`) sets the size of every field. Atmosphere/ocean (`--climate-resolution`, default 128) and moisture (`--climate-hydrology-resolution`, default 256) run on coarser grids and are remapped onto the world grid.

Land surface drag is a declared input of the climate profile, applied to `run` and `experiment` alike: `UW_CLIMATE_SURFACE_DRAG_INPUT` names the Earth vegetation drag contrast (`refs/processed/era5-vegetation-drag/1024/land-drag-contrast-1024.txt`, prepared from ERA5 vegetation cover). Other world widths get a version area-averaged from the 1024 file once, in `out/inputs/`. `run.json` records the file and its SHA-256. Procedural worlds never use it.

The legacy climate chain (`climate_benchmark.py`, `write_benchmark_report.py`, `climate_report.py`, `climate_experiment.py`, `check_climate_baseline.py`, `run_report.py`) reads the sealed run layout and is being retired; its runs and reports are in `runs/archive/` and `runs/reports/`.

## Full world

```powershell
# Saved to out/bench/world-seed1-512.uww
uv run --offline scripts/uw.py world run --resolution 512
uv run --offline scripts/uw.py world run --resolution 512 --plates 12 --tectonic-duration-myr 300
uv run --offline scripts/uw.py world run --resolution 256 --social --social-mode historical --history-years 500
```

`--no-rivers`, `--no-lakes` and `--no-deltas` shorten runs.
