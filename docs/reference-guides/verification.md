# Verification commands

Run these only when asked, or when a bug persists after a fix attempt (see `AGENTS.md`). A compile build is the default check for any edit. Two presets: `desktop` (`out/build/desktop`) and `headless` (`out/build/headless`); binaries land in `bin/<Config>/`, and scripts use the headless Release paths from `scripts/uwtools/paths.py`.

## Build

```powershell
cmake --preset desktop                                   # configure (Visual Studio generator, multi-config)
cmake --build out/build/desktop --config Release --parallel 4   # simulations: always Release
cmake --build out/build/desktop --config Debug --target world_pipeline   # debugging only
```

Headless configuration without the desktop stack or CUDA discovery (tectonics tools and tests included):

```powershell
cmake --preset headless
cmake --build out/build/headless --config Release --target world_pipeline --parallel 4
```

## Tests

```powershell
ctest --test-dir out/build/desktop -C Release --output-on-failure -j 4          # all registered C++ tests
ctest --test-dir out/build/desktop -C Release -R "^ClimateHydrology\."         # one suite; -R "^Suite\.Name$" for one test
ctest --test-dir out/build/headless -C Release -L tectonics                      # tectonics_tests (label geology: minerals_tests)
uv run --offline python -m unittest discover -s tests/scripts                    # Python suite
uv run --offline python scripts/check_layout.py                                   # layout and tectonics header boundary
uv run --offline python tests/scripts/check_tectonics_cli.py out/build/headless/bin/Release/world_pipeline.exe
```

Three GoogleTest binaries: `world_tests` (`tests/CMakeLists.txt`: core, I/O, pipeline, climate, validation; it links the planet-free set `uw_climate_world` alone), `tectonics_tests` and `minerals_tests` (`tests/simulations/geology/*/CMakeLists.txt`). ctest discovers every `TEST` and `TEST_F` as `Suite.Name`; since Phase 6 each file of `world_tests` that had one long test is a fixture class named after its suite (`class ClimateOceanDynamics : public ::testing::Test`) with one `TEST_F` per section, so a failure names the section (`ClimateOceanDynamics.PhaseChange`). A new test in such a file is a `TEST_F` on that fixture (GoogleTest refuses `TEST` and `TEST_F` in one suite); the fixture holds the grid and the config states the sections share and rebuilds them per test. The replay and benchmark modes of the climate tests (`--atmosphere-replay`, `--atmosphere-benchmark`, `--momentum-replay`, `--coastal-replay`, `--heat-kernels`, `--accuracy`, `--ocean-heat-benchmark`, `--transport-benchmark`) are disabled tests that ctest never runs: set `UW_TEST_ARGUMENTS` to the mode and its arguments and run `world_tests --gtest_also_run_disabled_tests --gtest_filter=<Suite>.DISABLED_Modes` (`tests/test_modes.hpp`). `check_tectonics_cli.py` checks what only the command line can show: options it refuses, checkpoint times independent of the order of the options, checkpoint labels and update counts, and `--mineral-history off` leaving the tectonic fields unchanged. Values the scenario options or the engine refuse are `TectonicOptions.*` in `tectonics_tests`; recipe equality (a run with no scenario options equals one given every value of `configs/tectonics-baseline.json`), the column and mineral budgets and the snapshot record are checks of the tectonics profile (below). The Python dependencies are declared once in `pyproject.toml` and pinned in `uv.lock`; `uv sync` installs them into `.venv/` (online the first time).

## Equivalence and access checks

```powershell
out/build/headless/bin/Release/world_pipeline.exe generate --seed 1 --bundle out/scratch/before
# ... change code, rebuild ...
out/build/headless/bin/Release/world_pipeline.exe generate --seed 1 --bundle out/scratch/after
uv run --offline python scripts/evaluation/compare_bundles.py out/scratch/before out/scratch/after
out/build/headless/bin/Release/world_pipeline.exe generate --seed 1 --audit-access   # exit status 1 after an undeclared or frozen write (outputs are still written)
# Earth climate (the climate profile's run) and a start from a checkpoint bundle:
out/build/headless/bin/Release/world_pipeline.exe run --recipe earth-climate --seed 1 --bundle out/scratch/climate
out/build/headless/bin/Release/world_pipeline.exe run --recipe procedural-world --from-bundle CHECKPOINT --steps a..b --bundle out/scratch/part
uv run --offline python scripts/evaluation/process_access.py (Get-ChildItem runs/reports/phase4-access/*.json, runs/reports/phase4b-access/*.json, runs/reports/phase4c-access/*.json) --check   # the access file equals what the stored reports give
# Determinism of a profile run: two runs, bundles compared by the profile's repeat check (out/evaluation/<module>/);
# --explicit-recipe adds a run given the recipe file explicitly (the profile's recipe-equality check)
uv run --offline scripts/uw.py tectonics experiment --duration 20 --repeat --explicit-recipe
```

The check of one refactoring step (headless build, the three minimal seed-1 bundles with `--audit-access`, comparison with the previous step's, `parameters.py --check`, `check_layout.py`, `process_access.py --check`):

```
uv run --offline python scripts/evaluation/step_check.py <label> --reference <previous label> [--desktop] [--social]
```

A label is used once: the check refuses a directory that already holds any `<label>-*` entry, a log written there by hand included, so keep your own logs in a subdirectory.

Parameters: edit `configs/parameters/<module>.json`, then `uv run --offline python scripts/parameters.py` writes the headers (`--check` fails on a stale one). One run with another value: `world_pipeline run ... --set hydrology.rivers.minimum_flow=80`, or `--parameters FILE.json`; the bundle's `manifest.json` lists every value under `parameters`.


A registered run (`uw.py <module> run`) reproduces an earlier one when `compare_bundles.py` finds its `bundle/` identical to the earlier run's final bundle (sealed runs: `runs/archive/checkpoints/<module>/<id>/final`).

Headless checks do not establish GUI rendering, historical save compatibility, full Earth benchmark skill or high-resolution performance; verify that boundary when it is part of the requested outcome. Module benchmarks and evaluations: `docs/reference-guides/benchmark-commands.md`.
