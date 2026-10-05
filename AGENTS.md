# Undiscovered Worlds

Procedural planetary simulation (terrain, tectonics, climate, hydrology, early society) in C++17/CMake/vcpkg on Windows. `README.md` describes the intended lifecycle; distinguish intended from implemented behaviour. Design and status live in `docs/architecture/` (`world-model.md`: components, process kinds, recipes; `consolidation-plan.md`: target structure, phases and open decisions). `meta/`, `refs/` and `runs/` (except `runs/registry/` and `runs/ledger/`) are local-only and gitignored.

## Working conventions

- Work directly on `main`; uncommitted changes in the working tree are current work, so build on them and never reset, stash or revert them.
- A compile build is the default check. Run tests, benchmarks, profile experiments, `check_tectonics_cli.py` or the Python suite only when asked or when a bug persists after a fix; say which checks were skipped. Commands are in `docs/reference-guides/verification.md`.
- Seed 1 for every new run; another seed only when the user asks for a seed comparison (20260906 only to reproduce climate runs up to 365). Module checks go through `scripts/uw.py <module> <command>` and the headless CLIs; never launch `UndiscoveredWorlds.exe` for a check.
- Tectonics has exactly one implementation, `src/simulations/geology/plate_tectonics` (namespace `tectonics::`, baseline from `configs/tectonics-baseline.json`). Do not reference or describe any other or former tectonics code. `scripts/check_layout.py` enforces its public-header boundary.
- Delegate read-only surveys to subagents and keep working. Before claiming a run improved anything, have a fresh-context reviewer check the evaluation. Never write per-run HTML or scripts: reports are rendered from JSON by `scripts/evaluation/`.

## Where things go

- Behaviour in `src/simulations/<module>/`, stage order and recipes in `src/pipeline/`, I/O in `src/io/`, observed-data comparisons in `src/validation/`, GUI only in `src/app/`. One translation unit per process, named by its process id. `src/world/functions.hpp` is a compatibility umbrella; include the owning header instead.
- Downloads and receipts in `refs/source/`, prepared inputs in `refs/processed/`, run artefacts in `runs/`, builds and scratch in `out/`. Bulk directories are junctions into `D:/dev/undiscovered-worlds/` (`runs/storage.json`): keep repo-relative paths, never replace a junction with a directory, and create run and report directories through `scripts/uwtools/storage.py` (`run_directory`, `route`). Registered runs live in `runs/<module>/<id>/` (`scripts/uwtools/runs.py` finds a run in either layout); `runs/archive/` is the sealed earlier layout and is never written. `scripts/repository_scope.py` enforces the 5 MiB blob limit and local-only paths.

## Simulation evidence

- Separate deterministic correctness from physical validity: an unchanged seed or a passing test does not show improved skill; compare the affected fields and metrics, and expect a visual review. Observed Earth data are validation inputs, never hidden forcing for procedural worlds; keep assumptions and calibrations explicit.
- A refactor that must not change results is checked with seed-1 `world_pipeline --bundle` runs compared by `scripts/evaluation/compare_bundles.py`; a stage whose writes change keeps its declared access current (`--audit-access`, `scripts/evaluation/process_access.py`).
- New modules get an evaluation profile in `configs/evaluation/<module>.json`, not their own benchmark scripts. Grids are cell-centred `W x W/2`. Record seed, parameters, code version, grid and inputs with every comparison; never overwrite prior run evidence. `runs/registry/<module>.json` supplies run IDs; use numeric IDs in new reports. The fixed global-balanced-v1 score gives the numerical recommendation; the user's visual verdict overrides it. Runs and references use the same palettes, `configs/climate-palettes.json` and `configs/tectonics-palettes.json`.
