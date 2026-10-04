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

Three GoogleTest binaries: `world_tests` (`tests/CMakeLists.txt`: core, I/O, pipeline, climate, validation; it links `uw_climate_world` alone), `tectonics_tests` and `minerals_tests` (`tests/simulations/geology/*/CMakeLists.txt`). ctest discovers every `TEST` as `Suite.Name`. `world_tests` also has explicit replay and benchmark modes that ctest never runs (`--atmosphere-replay`, `--atmosphere-benchmark`, `--momentum-replay`, `--coastal-replay`, `--heat-kernels`, `--accuracy`, `--ocean-heat-benchmark`, `--transport-benchmark`). `check_tectonics_cli.py` checks default/explicit recipe equality, column budgets and rejection of invalid options. The Python dependencies are declared once in `pyproject.toml` and pinned in `uv.lock`; `uv sync` installs them into `.venv/` (online the first time).

## Equivalence and access checks

```powershell
out/build/headless/bin/Release/world_pipeline.exe generate --seed 1 --bundle out/scratch/before
# ... change code, rebuild ...
out/build/headless/bin/Release/world_pipeline.exe generate --seed 1 --bundle out/scratch/after
uv run --offline python scripts/evaluation/compare_bundles.py out/scratch/before out/scratch/after
out/build/headless/bin/Release/world_pipeline.exe generate --seed 1 --audit-access   # exit status 1 after an undeclared or frozen write (outputs are still written)
# Earth climate (the climate profile's run) and a start from a checkpoint bundle:
out/build/headless/bin/Release/world_pipeline.exe run --recipe earth-climate --seed 1 --map none --bundle out/scratch/climate
out/build/headless/bin/Release/world_pipeline.exe run --recipe procedural-world --from-bundle CHECKPOINT --steps a..b --bundle out/scratch/part
uv run --offline python scripts/evaluation/process_access.py (Get-ChildItem runs/reports/phase2-access/*.json) --check   # the access file equals what the stored reports give
```

Headless checks do not establish GUI rendering, historical save compatibility, full Earth benchmark skill or high-resolution performance; verify that boundary when it is part of the requested outcome. Module benchmarks and evaluations: `docs/reference-guides/benchmark-commands.md`.
