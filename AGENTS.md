# Undiscovered Worlds

## Project context

- Procedural planetary simulation with terrain, tectonics, climate, hydrology, and early society systems. The fuller ecology and civilization lifecycle in `README.md` is a design goal; distinguish intended behavior from implemented behavior.
- C++/CMake/vcpkg project on Windows. The desktop app uses C++17 with SFML/ImGui; the internal tectonics module uses C++20. CUDA rendering is optional. `CMakePresets.json` defines the Visual Studio 18 2026 baseline and uses `VCPKG_ROOT`.
- Simulation behavior belongs in `src/simulations/`; stage ordering and dependencies in `src/pipeline/`; artifact I/O in `src/io/`; observed-data comparisons in `src/validation/`; presentation in `src/app/`. Include the owning module's header in new callers; `src/wip/functions.hpp` is a compatibility umbrella.
- Read `README.md` for project scope and commands, `docs/architecture/simulation-modules.md` for module changes, and `docs/architecture/tectonics/tectonic-output-contract.md` for tectonic output changes, as relevant to the task.

## Simulation evidence

- For physics changes, check the affected units, budgets/conservation laws, finite values, grid registration, and stage dependencies. Simulator and processed Earth benchmark grids are cell-centred `W x W/2`; raw sources retain their native registration until ingestion.
- Separate deterministic correctness from physical validity. An unchanged seed or passing unit test does not establish improved climate skill; compare affected numerical fields and metrics on representative inputs when making that claim. A plausible map alone is insufficient.
- Observed climate references are validation inputs, not hidden forcing for procedural worlds. Keep physical assumptions, approximations, and calibration choices explicit.
- For simulation comparisons, record seed, parameters, code/input versions, grid geometry, duration/convergence criteria, and output provenance. Preserve numerical fields needed to reproduce the comparison; avoid overwriting prior run evidence.
- Use `refs/source/` for downloads and receipts, `refs/processed/` for prepared inputs, `runs/` for simulation artifacts, and `out/` for builds and scratch work. The active `runs/registry/climate.json` supplies run IDs; retain its IDs and role unless its replacement is authorized.
- Local bulk-data storage is configured in ignored `runs/storage.json`. Generated output directories and experiment archives use directory junctions into `D:/dev/undiscovered-worlds/`; keep code, references, builds, registry and active working files on C:. Preserve logical repo-relative paths. Before creating a new direct `runs/reports/<experiment>` archive, call `scripts/benchmarks/run_storage.py`'s `ensure_archive_directory`, or `configure-run-storage.ps1 -RelativePaths runs/reports/<experiment> -EnsureArchive`. The shared climate experiment runner does this automatically. Do not replace junctions with ordinary directories or migrate data while a writer is active.
- Climate comparison reports use `scripts/benchmarks/climate_report.py` and the common `configs/climate-palettes.json` palette contract. Use numeric run IDs in new report controls, filenames and narratives; historical archive aliases are input adapters only. Keep visual maps/charts in HTML, hypotheses/findings/decisions in the generated Markdown report, and raw execution output in logs. Use the fixed global-balanced-v1 score for numerical baseline recommendations and honor user visual overrides. No historical wind run is permanently protected. See `docs/reference-guides/climate-reporting.md` for commands and definitions.

## Verification commands

Choose commands for the changed surface. The C++ and Python suite commands below are broader checks, not prerequisites for every edit.

```powershell
cmake --preset x64-debug
cmake --build out/build/x64-Debug --config Debug --parallel 4
ctest --test-dir out/build/x64-Debug -C Debug --output-on-failure -j 4
uv run --offline --with numpy --with pillow python -m unittest discover -s tests/scripts
uv run --offline python scripts/check-layout.py
```

- For targeted C++ work, build the affected test target with `--target` and select its registered test with CTest `-R`. `tests/CMakeLists.txt` and `tests/simulations/geology/plate_tectonics/CMakeLists.txt` define the names. `UW_BUILD_APP=OFF` supports numerical library/test work without desktop or CUDA discovery; see the headless build example in `README.md`.
- For climate benchmark changes or claims about physical skill, use the relevant workflow in `runs/maps/climate/README.md`; preparation and analysis commands are in `refs/README.md` and `scripts/benchmarks/`. Keep resolution, inputs, and comparison settings explicit.
- Default CTest excludes the full-size tectonics replay. Run `cmake --build out/build/x64-Debug --config Debug --target check-tectonics-slow` when that coverage is required. Two historical golden tests skip obsolete baselines; a suite pass does not validate those baselines.
- Headless checks do not establish GUI rendering, historical save compatibility, full Earth benchmark skill, or high-resolution performance. Verify the relevant boundary when it is part of the requested outcome.
