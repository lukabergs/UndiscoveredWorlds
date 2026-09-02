# Goal

Make every climate benchmark produce a predictable, selectable artifact set and a mandatory `climate.xlsx` row from an explicit simulation resolution.

# Phase

Climate physics phase 6: benchmark pipeline and artifact organization.

# Status

Complete on the feature branch.

# Scope

- Add a width-derived `--resolution` CLI option while retaining explicit width and height overrides.
- Add repeatable/selectable benchmark map arguments with a compact naming convention.
- Default to Köppen, January surface LIC, January surface particles, precipitation PNG, and precipitation GeoTIFF.
- Make workbook population part of every recorded climate benchmark.
- Reorganize reference, benchmark, workbook, and legacy generated assets without deleting files.
- Run and inspect one 512-wide Release benchmark.

# Constraints

- Never delete existing artifacts; ask before any later deletion.
- Preserve deterministic simulation behavior and climate physics.
- Keep reference observations validation-only.
- Skip costly flow visualizations automatically at unsupported large resolutions.
- Preserve explicit `--world-width` and `--world-height` compatibility.

# References

- `main.cpp`
- `climate_validation.cpp`
- `climate_circulation_diagnostics.cpp`
- `app.env`
- `docs/tasks/climate-physics/FUTURE_WORK.md`

# Steps

1. Completed: defined and tested benchmark map-selection and resolution contracts.
2. Completed: applied compact names and generated only selected maps.
3. Completed: made workbook updates mandatory for benchmark recording.
4. Completed: moved non-code assets with an idempotent, non-overwriting migration script.
5. Completed: built, ran deterministic tests, ran a 512-wide Release benchmark, and inspected its workbook and maps.

# Done when

- A default Earth benchmark produces exactly the requested five maps.
- Arbitrary supported maps and a 512/2048 width can be selected from the command line.
- The new run has a corresponding workbook row.
- Existing reference and legacy outputs remain available under organized folders.

# Verification

- `cmake --build out/build/x64-Debug --config Debug`
- `ctest --test-dir out/build/x64-Debug -C Debug --output-on-failure`
- `cmake --build out/build/x64-Release --config Release`
- `ctest --test-dir out/build/x64-Release -C Release --output-on-failure`
- Release benchmark run 145 at 512x257.
- Five-file artifact manifest, image dimensions, GeoTIFF metadata, JSON log, and `climate.xlsx` formula/render inspection.

# Open questions

- High-resolution LIC/particle limits remain performance policy rather than climate physics.
