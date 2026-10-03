# Consolidation plan

Written 2026-10-03 as the brief for the sessions that carry out the restructuring. It extends `world-model.md` (components, process kinds, recipes, evaluation) with the target structure of executables, tests, tooling, data layout and the working agreement. Decisions marked **open** need the user's answer before the phase that depends on them.

## 1. Decisions taken on 2026-10-03

- Executables: two programs (`UndiscoveredWorlds.exe` GUI, `world_pipeline.exe` CLI) and three test binaries. `climate_pipeline`, `tectonic_pipeline`, `tectonic_column_probe` and `tectonic_climate_erosion` fold in or go.
- Rule: **C++ produces run bundles; Python compares them.** Comparison, statistics and maps live in `scripts/evaluation/`, never in a C++ subcommand.
- Rule: **a pass is run by the scheduler**, never by a bespoke executable. Per-process benchmarks use `--from-bundle`, step selection and checkpoints.
- Bundles (`uw-run-bundle/v1`: `manifest.json` + one `.npy` per field) replace `.uww` as the save format once the GUI can load them. Reason: open formats readable by any tool, inspectable by the user.
- Reference data migrate to the same open layout (`.npy` + `manifest.json`), replacing `.uwclim`; region masks come from one `configs/regions.json` (seeded from the IPCC AR6 regions already under `refs/source/ipcc-ar6-wgi-regions-v4` plus the project's own regions).
- Builds: the Visual Studio generator stays (needed for debugging in Visual Studio; VS Code's CMake Tools reads the same presets). Two configure presets only, `desktop` and `headless`; **configuration is chosen at build time** (`--config Release` for every simulation and benchmark, `--config Debug` only to debug a crash). Binaries land in `bin/<Config>/`; scripts default to Release and record the configuration in `run.json`.
- Tooling: Python only. `configure-run-storage.ps1` is ported; the Excel/PowerShell automation leaves the simulator (`run.json` carries timings; an optional Python workbook writer if ever wanted).
- Reports: standardized templates rendered from JSON; agents never author per-run HTML or scripts. An append-only experiment ledger per module.
- Review: a fresh-context verifier agent answers a user-written visual checklist per module before the user looks; the user's verdict overrides.
- No cloud or remote sessions are used (subscription via the VS Code extension and Codex); Linux builds are not a goal. Windows-only code is tolerated but should stay out of `src/core` and `src/simulations`.
- The `version.txt` download at GUI start (`src/app/main.cpp`, `getlatestversion`) is unintended and goes.
- Time model: each evolver declares its time step and duration in the recipe (compiled defaults as for `configs/tectonics-baseline.json`); the CLI overrides them (`--duration`, `--checkpoint-every`, `--checkpoint-at`); each checkpoint writes a bundle.

Closed later the same day:

- `climate-benchmark` switches to the frozen-Geosphere `earth-climate` recipe and registers its runs; the first run after the switch is the new baseline (step 5b already requires one).
- Field encodings migrate **all at once** to float32 in physical units (continuous quantities) and to properly typed small integers for categories and codes. Together with the recipe switch this is the single re-baseline of Phase 1. Reason: whole-degree temperatures and integer metres lose real precision and bake rounding into every downstream process; float32 costs twice the memory of int16 per field (8 MB instead of 4 MB at 2048×1024) and no measurable compute time; float64 would double that again for no benefit, so it is not the target.
- **Parameters**: everything that has ever been tuned moves to the recipe. A constant stays in code only when it is a textbook physical constant with units (Stefan–Boltzmann, latent heats, the gas constant) or belongs to the planet (`PlanetBody`: radius, gravity, rotation, tilt). Dimensionless factors, albedos, time scales, amplitudes, widths and switches are tunables: there are 645 `constexpr` values in `src/wip/generation_tuning.hpp`, and even `gravity` there is a scale factor (1.0f), not 9.81. The user's reason: agents have re-tuned "settled" parameters while tuning others, so coupling is high and every value must be visible in the run record. The agent doing the migration classifies by inspection (git history is a single commit and cannot say what was last touched) and lists the few kept constants in this document.
- **Tracking**: `docs/` (23 files, 276 KB) becomes tracked in full, including `docs/tasks/`. From `runs/`, only `runs/registry/` (climate.json 1.3 MB, grows about 3.5 KB per run) and the new `runs/ledger/` are tracked; everything else under `runs/`, and all of `refs/` and `meta/`, stays local. `.gitignore` (lines for `/docs/`, `/docs/tasks/`, `/runs/`) and `LOCAL_ROOTS` in `scripts/repository_scope.py` change accordingly; the 5 MiB blob limit stays.
- **C++17 for every target**, one `CMAKE_CXX_STANDARD` at the root; the engine's C++20 setting goes (a text scan found no C++20-only construct; the headless build confirms).
- **Determinism**: results must be bit-identical across thread counts. Phase 0 tests it (seed 1 at 4 and 8 threads, `compare_bundles.py`); a difference is a bug to fix in the reduction order, not a configuration to pin. Reason: `repeat` checks, refactor equivalence and the ledger's footprint diff are only trustworthy if the machine's core count cannot change a result. If a fix ever proves expensive, the fallback is to record the thread count in `run.json` and compare runs with equal counts, but that is a documented exception, not the policy.

## 2. Target: executables, libraries, tests

```
UndiscoveredWorlds.exe   GUI: a front end over the same recipes; loads bundles; no generation logic of its own
world_pipeline.exe       the only headless entry point
  run        --recipe R [--seed N] [--resolution W] [--from-bundle DIR] [--steps a..b | --only <process>]
             [--duration ...] [--checkpoint-every N] [--audit-access] --bundle DIR
  benchmark  a registered standard run for a module (climate today; tectonics when its profile exists)
world_tests.exe / tectonics_tests.exe / minerals_tests.exe   GoogleTest; one FetchContent at tests/; ctest sees every TEST
```

Libraries, one per directory, STATIC, named `uw_<dir>`: `uw_core`, `uw_io`, `uw_world` (components and the `planet` facade, ex-`wip`), `uw_pipeline`, `uw_validation`, `uw_cli`, `uw_app`, `uw_<module>` under `src/simulations/`, and the engine `plate_tectonics` (rename to `uw_tectonics` when convenient). One `CMAKE_CXX_STANDARD` at the root; one include root (`#include "pipeline/world_processes.hpp"`); `check-layout.py` enforces layering: `core` includes nothing above it; `simulations` and `pipeline` never include `app` or `cli`; only `app` includes SFML/ImGui.

Source tree:

```
src/core          grids, noise, random, field catalogue, world grid, environment (paths), profiler
src/world         planet components + facade, region, regionalmap (GUI-only reachability), tuning until parameters move to recipes
src/simulations/<module>/   processes; one .cpp per process id; tectonics engine inside geology/
src/pipeline      process registry, recipes, scheduler, stages, checkpoints
src/io            bundles (read and write), imports, TIFF/PNG, saves
src/validation/{climate,tectonics}   comparisons with observed data; tectonic snapshot statistics (ex tools/tectonics support)
src/cli           world_pipeline main, command_line_generation, climate_options
src/app           GUI only: main, windows, panels, rendering, workbench
tests/            world_tests (core, io, pipeline, climate, validation), tectonics_tests, minerals_tests
tools/            removed
```

Tectonic snapshots become run bundles (fields registered in the catalogue), so `tectonic_pipeline inspect/batch-export/compare/boundary-stats` reduce to `run --recipe tectonics --checkpoint-every N --bundle` plus Python analysis; `check_tectonics_cli.py` becomes profile checks (`repeat`, recipe equality).

## 3. Tooling

- `scripts/uwtools/` package: bundle reader, `.npy` reference reader, SHA-256, area weights, palettes, build paths, ledger. Replaces the 14/11/5 copies.
- `pyproject.toml` at the root pins the Python dependencies once; scripts drop their inline `# /// script` headers.
- `bench.py <module> <command>` is the single entry: `experiment` (run → evaluate → verifier review → ledger entry, one command), `evaluate`, `review`, `ledger`, `compare`, `report`.
- The legacy climate chain (`climate_benchmark.py`, `write-benchmark-report.py --legacy-layout`, `climate_run_record.py`, `check-climate-baseline.py`, `run_report.py`) retires after one parity check of the global-balanced-v1 ranking against the registry.
- Reports: three templates in `scripts/evaluation/templates/`: run page, comparison page (N runs against baseline and reference; module map lists from the profile; tectonics animations included), module dashboard over the registry. Rendered from `run.json`, `evaluation.json`, `review.json`, ledger. Markdown report generated from the same data.
- Ledger: `runs/ledger/<module>.jsonl`, one line per experiment: `id, parent, hypothesis, change (parameters or files + diff hash), expected, observed, assessment, verdict, decision (keep|revert|park), lessons, do_not_retry`. `bench.py <module> experiment` prints the relevant entries (`--topic`) before running, so reading the ledger does not depend on an instruction. A generated `LEDGER.md` with a "settled / don't retry" section lives beside it.
- Verifier: `.claude/agents/verifier.md` (Opus, high effort) receives the review page (model | reference | difference per field and season, same palette), the metric table with deltas against the parent, and the module checklist (`configs/evaluation/<module>.json`, `review` array); it is blind to the hypothesis and writes a `review.json` entry as `reviewer: agent`. Codex sessions run the same prompt by hand. The user's verdict is `reviewer: user` and overrides.
- Added validation guards (reported, not ranked): variance ratio and centered pattern correlation per field/season; regional table over `configs/regions.json`; feature diagnostics (ITCZ latitude by month, 500 hPa jet core, monsoon reversal index, Köppen area fractions and κ); parent-diff footprint (which fields changed, where).

## 4. Data layout and palettes

Target, run-centric and identical for references:

```
runs/
  registry/<module>.json          ids and status (kept)
  ledger/<module>.jsonl           experiment ledger
  <module>/<id>/                  everything of one run (D: junction per module)
    run.json  log.txt  bundle/  evaluation/{evaluation.json, review.json, index.html, maps/}
  reports/<name>/                 cross-run comparison pages only
  archive/                        sealed pre-consolidation layout (fields/, maps/, diagnostics/, manifests/, parameters/, checkpoints/, metrics/, logs/) — evidence, never migrated
refs/
  source/<dataset>/               downloads and receipts (kept)
  processed/<dataset>/<width>/    fields/ as .npy + manifest.json (bundle schema); maps/ rendered by scripts/evaluation/maps.py
configs/regions.json              masks shared by validation, review and reports
```

Palette rule: every map under `refs/processed/**/maps` and `runs/**` is rendered by the shared renderer from `configs/climate-palettes.json` or `configs/tectonics-palettes.json`; each maps folder carries a `palette.json` sidecar with the contract version; `check-layout.py` fails on a missing or mismatched sidecar. Legacy renders (`refs/processed/climate/maps/climate/`, `refs/processed/climate/previews/`) move to `meta/legacy-maps/` or are deleted.

## 5. Phases

**Phase 0 — hygiene (no behaviour change).** Delete (Recycle Bin) the items in §6; `CMakePresets.json` → `desktop` and `headless`, drop the build presets; remove the per-config `RUNTIME_OUTPUT_DIRECTORY_*` overrides so binaries land in `bin/<Config>/`; scripts default to Release paths from one place (`scripts/uwtools/paths.py`); one `CMAKE_CXX_STANDARD 17` at the root, engine C++20 setting removed, headless build confirms; thread-count determinism test (seed 1 at 4 and 8 threads); remove `getlatestversion` from `main.cpp`; `.gitignore` and `repository_scope.py` updated for tracked `docs/`, `runs/registry/`, `runs/ledger/`; fix README layout entries (`references/`, `definitions/society/`, `0/`), `tests/README.txt`, `.gitignore` dead entries; drop the `check-layout.py` exemption for the probe; `refs/README.md` and `README.md` point to `scripts/refs/prepare-earth-maps.py` only. Everything in this phase is checked by a seed-1 bundle comparison before and after.

**Phase 1 — one CLI and per-pass runs.** Bundle reader in `src/io`; `--from-bundle`, `--steps`, `--only`, checkpoints; fold `climate_pipeline` (its planet-free link guarantee survives as `uw_climate_world` + `world_tests` linking only it); evaluation profile `run.tool` → `world_pipeline`. Then the one re-baseline: `climate-benchmark` runs the frozen-Geosphere recipe and registers, all field encodings move to float32/typed integers, and the first registered run afterwards is the new baseline, compared with the last pre-switch run in the ledger so the shift is documented rather than silent. Then tectonic snapshots as bundles, support code to `src/validation/tectonics` and `src/io`, `tectonic_pipeline` and `tools/` removed, `check_tectonics_cli.py` → profile checks.

**Phase 2 — tests.** Three GoogleTest binaries; one FetchContent at `tests/`; fix `test_plate.cpp`'s own `main`, empty `test_sqrdmd.cpp`, `PLATE_TECTONICS_FIXTURE_DIR`; `gtest_discover_tests`.

**Phase 3 — tooling and data.** `uwtools`, `pyproject.toml`, templates, ledger, verifier agent, `regions.json`, reference migration to `.npy`, new run layout for new runs (old layout sealed under `runs/archive/`), palette sidecars and the check, PowerShell port, Excel removal, legacy climate chain retired after parity.

**Phase 4 — source tree and CMake.** `src/cli`, `src/world`, `src/app` GUI-only, per-directory CMake, one standard, one include root (mechanical rewrite of ~700 includes), layering check. Split the giants: `assetdata.cpp` (data to `assets/`), `climate_coupling.cpp` by sub-model, `climate_validation.cpp` by concern (most of it moves to Python), `main.cpp` by screen; `misc.cpp` and `functions.hpp` dissolve under the world-model migration.

**Phase 5 — GUI as viewer.** Loads bundles; the "Generate elements" buttons call processes; `regionalmap` leaves the headless link line; parameters move from `generation_tuning.hpp` into recipes.

Each phase ends with: seed-1 bundle equivalence (`compare_bundles.py`) where results must not change, `--audit-access`, `check-layout.py`, a compile of both presets, and a ledger entry describing what changed.

## 6. Deletion and cleanup list (Phase 0)

Tracked files (git will show them as deleted; the user commits):
- `tools/tectonics/column_process_probe.cpp`, `tools/tectonics/climate_erosion.cpp` and their targets; `scripts/benchmarks/tectonics_climate_matched_control.py`; `tools/tectonics/heightmap_stats.py`. (`uw_climate_erosion` library stays until a `terrain.climate_erosion` evolver is decided. `scripts/benchmarks/tectonics_climate_feedback.py`, which alternated climate and that erosion tool on a finished tectonics run, went with the tool on 2026-10-03; its one registered run, tectonics 63, stays readable through `tectonics_report.py`. The experiment returns as a recipe when the evolver exists.)
- Unreferenced scripts: `scripts/benchmarks/{climate_history,climate_experiment_review,climate_ocean_depth_replay,tectonics_batch_summary,tectonics_climate_review,tectonics_process_review}.py`, `scripts/benchmarks/rainfall_comparison_dashboard.js`, `scripts/refs/download-climate-reference.ps1`, `tests/scripts/check_tectonics_animations.py`, `tests/scripts/check_map_animation.js`, `check_map_comparison.js`, `check_map_report_selection.js`, `check_climate_report_dashboard.cjs`.
- `refs/prepare.py` (copy of `scripts/refs/prepare-earth-maps.py`), `configs/seed.txt`, the `0/` directory, `docs/licensing/Licence.rtf` (duplicate of `LICENSE`), `docs/archive/TODO.txt`, `refs/TODO.txt`.
- The 4-line "Review (2026-09-03)" notes: `scripts/README.txt`, `tools/README.txt`, `tests/README.txt`, `scripts/archive/README.txt`, `assets/appearance/README.txt`, `assets/palettes/README.txt`, `docs/architecture/README.txt`, `docs/reference-guides/README.txt`, `docs/tasks/climate-physics/README.txt`, `runs/{manifests,parameters,registry,reports}/README.txt`.
- `scripts/archive/` → `meta/archive/scripts/` (historical evidence; one tracked file, `reorganize-climate-assets.ps1`).
- Code nits: `src/validation/climate/climate_circulation_diagnostics.hpp` stale `namespace sf` forward declaration; `getlatestversion` in `src/app/main.cpp`; unused `CMakePresets.json` entries; `.gitignore` entries for paths that no longer exist (`/extra/`, `/data/`, `/references/`, `/tools/ripgrep/`, `saved_worlds/`, the three root files).

Local-only (not in git):
- Stale build directories: `out/build/access-audit` (397 MB), `out/build/grid-check` (544 MB), `out/build/headless` (530 KB), `out/build/state-review` (20 KB). Keep `x64-Debug`, `x64-Release`, `headless-modules` until the presets are renamed, then rebuild into `out/build/desktop` and `out/build/headless` and delete the rest.
- `out/tectonics-report-syntax.js` (28 MB), `out/coastal-report-syntax.js`, `out/coastal-batch-report-syntax.js` (9 MB each): syntax dumps, no consumer. Review the remaining `out/*` experiment leftovers (`scratch` 364 MB, `editor`, `storage-migration`, `tectonics-cleanup`, `geology-fixes`, …) and delete what no report references.
- `D:/dev/undiscovered-worlds/git-recovery/20260918` (203 MB): a working-tree snapshot made on 2026-09-18 by an earlier agent session; nothing references it.
- `refs/processed/climate/vegetation-drag-increase-20260914/` (2.8 MB, cited only by report records): delete. Corrected on 2026-10-03 after a reference check, all deferred to Phase 3: `maps/climate/` is read by `write-benchmark-report.py` and `previews/` is the app's default preview directory, so both go only when their consumers are replaced; `vegetation-drag-20260914/` is `DRAG_DIRECTORY` of `climate_benchmark.py` and moves to a dataset-named folder under the reference layout of §4 (the user: its location and name make no sense); `csv/` (7.9 GB) and `fields/` (3.7 GB) are not empty and are reorganized and renamed with the reference migration.
- `out/` holds only disposable output. Kept until their consumer changes: `out/storage-migration` (receipt directory of `configure-run-storage.ps1`, until the Python port in Phase 3); `out/climate-history-audit` and `out/climate-surface-audit` (cited as evidence by `runs/reports/climate-history-ocean-corrections.md`, `climate-surface-corrections.md` and `climate-historical-comparison-181.json`): Phase 3 moves them under `runs/archive/` with the sealed layout and updates the citing paths. `out/diagnostics` (junction) and its target go, and with them `D:/dev/undiscovered-worlds/out`: nothing under `out/` is routed to D:.
- `D:/undiscovered-worlds-data` (11.2 GB): a copy of `refs/{source,processed}/physical`, every file present in the repository's `refs/` at the same size, unreferenced; delete.
- `runs/diagnostics/<45 per-diagnostic folders>` and `runs/diagnostics/climate/<id>/` are two layouts of the same data; `runs/manifests/{map_manifest,wind_reference_visualization_cache}`, `runs/logs/workbook_update_error`. Sealed: move under `runs/archive/` in Phase 3, never delete.
- Outside the repo, for the user's judgement only: `C:/dev/plate-tectonics` (separate project folder, last touched 2026-09-03, not referenced by this repository).

## 7. Explanations recorded for the user

**Parameters in recipes, not compiled constants.** Today many climate numbers (drag time scales, mixing depths, closure switches) are `constexpr` values in `src/wip/generation_tuning.hpp`. Changing one means recompiling, and a run's bundle cannot say which values produced it. The tectonics side already does it the other way: `configs/tectonics-baseline.json` holds the numbers, CMake compiles them in as defaults, and the CLI can override any of them. Doing the same for every process means an experiment is "edit a JSON value and run", the bundle records the exact values, and the ledger's `change` field is a parameter diff. Decide per constant whether it is physics (a fixed law) or a tunable; only tunables move.

**Field encodings.** Fields are stored as the original program stored them: temperatures as whole degrees in a 16-bit integer, elevation as integer metres above a datum, some categories packed into one array. Every process that reads them gets rounded values, and the rounding is part of today's results. Moving a field to 32-bit floating point in physical units removes that rounding, so downstream numbers change slightly, and every registered score shifts a little: a new baseline is needed. Options: migrate everything at once (one new baseline, one disruption) or one field per experiment (many small baselines). Recommendation: all at once, together with the benchmark-recipe switch, so there is one re-baseline.

**Determinism across threads.** The row-parallel executor splits work over threads. If any reduction depends on the order threads finish, two runs with different thread counts differ in the last bits, and `repeat` checks or refactor-equivalence comparisons fail for the wrong reason. Test: run seed 1 at 4 and 8 threads, `compare_bundles.py`. If identical, require it from now on; if not, either fix the reductions (fixed-order accumulation) or record the thread count in `run.json` and pin it for benchmarks.

**Bundles versus `.uww`.** `.uww` is the GUI's private binary save; only this program reads it. A bundle is a folder: `manifest.json` (readable) plus one `.npy` array per field, which NumPy, Python, QGIS plugins and the evaluation scripts read directly. Maps and reports are rendered from the bundle, not stored in it.

## 8. Working agreement

- One interactive session leads (design, decisions, reviews with the user). One writing session at a time on the working tree; parallel sessions are read-only. Each phase is its own session, started from this document and the previous phase's ledger entry.
- Models and effort (Claude): the lead on Fable 5.1 at `high`; implementation of agreed phases on Opus 5.5 at `high` (or Fable at `medium`); surveys on lower-power subagents; the verifier on Opus at `high` with fresh context. `xhigh` only after a failed attempt at `high`. Codex: `~/.codex/config.toml` already sets its model and effort.
- Instruction files: `AGENTS.md` (repo, provider-neutral) is the only project file; `~/.codex/AGENTS.md` is the global operating contract, imported by `~/.claude/CLAUDE.md`. Do not add a repo `CLAUDE.md` or `CLAUDE.local.md`: by default they would stop `AGENTS.md` from loading.
- A session ends with a ledger entry or a handoff note (`docs/tasks/main/handoff-<date>.md`: what was done, what was verified and what was not, open items), not with a summary only in chat.

## 9. Visual checklist format and sources

Checklist entries live in the module's evaluation profile under `"review"`:

```json
{"id": "itcz-single-band", "question": "Is the equatorial rain maximum a single band that migrates north in JJA and south in DJF?",
 "fields": ["atmosphere.precipitation_rate"], "seasons": "all", "region": "tropics",
 "good": "one band, 5-10 degrees wide, centred near 5N in JJA and near the equator or just south in DJF",
 "severity": "blocker"}
```

`severity` is `blocker` (reject regardless of score), `major` (revise) or `minor` (note). Public verification practice to borrow from: the PCMDI Metrics Package and Gleckler portrait diagrams (per-field, per-season normalized errors), ESMValTool performance-metrics recipes, ILAMB's scoring (bias, RMSE, seasonal phase, spatial distribution combined into 0–1 scores), NCAR AMWG diagnostics plot sets (zonal means, seasonal cycles, model-minus-observation maps), Taylor diagrams (correlation, variance ratio, centered RMSE). For tectonics: hypsometric curve against Earth (present), plate-area distribution (Bird 2003), land fraction, orogens at convergent margins, arcs behind trenches, passive-margin shelf widths, ridge geometry.

## 10. Naming (decided 2026-10-03)

Vocabulary, one word per thing:

- **recipe**: what is run (`procedural-world`, `earth-climate`, `tectonics`, later `geology`, `biology`, `society`). A module is isolated by choosing its recipe and `--from-bundle`/`--only`, never by a separate executable.
- **run**: one execution of a recipe, with an ID from `runs/registry/<module>.json`. Every module's runs are the same kind of thing, whether or not Earth reference fields exist for it; the reference only changes what the evaluation profile can compare.
- **experiment**: a run plus its evaluation, review and ledger entry.
- **baseline**: the registered run a module's experiments are compared with.
- **test**: only a ctest pass/fail check of deterministic correctness. Runs and experiments are never called tests: they are judged by evaluation and by eye, not by pass/fail.
- "benchmark" and "pipeline" stop being kinds of run. `climate-benchmark` becomes `world_pipeline run --recipe earth-climate` and `tectonic_pipeline` becomes `world_pipeline run --recipe tectonics` (Phase 1); the planned `benchmark` subcommand of §2 is dropped in favour of `run`. `world_pipeline` stays as the executable's name. `bench.py` becomes `scripts/uw.py <module> <command>` in Phase 3.

Dashes or underscores, by kind:

- kebab-case: recipe names, CLI subcommands and flags, directories, and data, config and document file names (`earth-climate`, `--from-bundle`, `tectonics-baseline.json`, `consolidation-plan.md`).
- snake_case: everything a compiler or interpreter names: C++ and Python source files, CMake targets and executables, Python packages, JSON keys, field and process ids (`world_pipeline`, `uw_core`, `compare_bundles.py`, `atmosphere.precipitation_rate`).
- Existing violations are renamed in the phase that already touches them, not separately: dashed Python scripts (`check-layout.py`, `write-benchmark-report.py`, `scripts/refs/prepare-*.py`, …) in Phase 3; the dashed keys of `configs/tectonics-baseline.json` when parameters move into recipes (a flag `--trench-depth` maps to the key `trench_depth`); run and reference directory names with the Phase 3 layout.
