# Tectonic Contract And Validation

## Upstream dependency contract

- `UndiscoveredWorlds` consumes tectonics from `C:\dev\plate-tectonics\src` via `PLATE_TECTONICS_ROOT` in `CMakeLists.txt`.
- The app expects upstream headers/API from that repo, including `platecapi.hpp` and `tectonic_contract.hpp`.
- No local vendored fallback exists. Do not reintroduce `third_party\plate_tectonics`.
- If cleaner CMake packaging is needed (for example, an exported target instead of source-globbing `src/*.cpp`), implement it in `C:\dev\plate-tectonics` and then consume it here.

Ownership boundary:

- Upstream (`plate-tectonics`) owns tectonic simulation, tectonic time semantics, provenance fields, boundary/deforming-region semantics, and tectonic export/checkpoint tools.
- App (`UndiscoveredWorlds`) owns ingestion into `planet`, app persistence, UI/debug surfacing, and downstream interpretation (climate/hydrology/resources/social layers).

## Build commands

App:

```powershell
cmake-vs --preset x64-debug
cmake-vs --build --preset build-x64-debug
```

Release app:

```powershell
cmake-vs --preset x64-release
cmake-vs --build --preset build-x64-release
```

Upstream tooling:

```powershell
cmake-vs --preset x64-debug
cmake-vs --build --preset build-x64-debug --target tectonic_pipeline
```

## Fixed-seed validation matrix

Use these seeds with `--plate-cycles 4`:

- `12345`
- `24680`
- `424242`

For each seed, run app generation:

```powershell
C:\dev\UndiscoveredWorlds\out\build\x64-Debug\bin\UndiscoveredWorlds.exe --generate-world --seed <SEED> --plate-cycles 4 --save C:\dev\UndiscoveredWorlds\out\validation\uw-<SEED>.uww
```

For each seed, run upstream inspection:

```powershell
C:\dev\plate-tectonics\out\build\x64-Debug\bin\tectonic_pipeline.exe inspect --final --seed <SEED> --cycles 4 --bundle-output C:\dev\UndiscoveredWorlds\out\validation\pt-<SEED>\final
C:\dev\plate-tectonics\out\build\x64-Debug\bin\tectonic_pipeline.exe boundary-stats --seed <SEED> --cycles 4 --update 25 --output C:\dev\UndiscoveredWorlds\out\validation\pt-<SEED>\boundary-stats.json
```

Compare for each seed:

- App world metadata: `time_myr`, `delta_time_myr`, `cycle_count` (and any surfaced plate/sea-level metadata).
- App point inspection at representative convergent/divergent/transform/stable points:
  crust age, crust thickness/class, uplift/subsidence tendency, accumulated strain, boundary ids/distances/types, deforming-region ids/types/rates/velocities.
- Upstream bundle summaries (`boundary_segments.json`, `deforming_regions.json`) versus app-inspector ids/types and broad rates/motions.
- Save/reload persistence of tectonic metadata and object collections in `.uww`.

Note on replay fidelity:

- Under the Windows shell harness used in this repo, prefer `cmake-vs` so MSVC include/lib paths are initialized before CMake or MSBuild runs.
- `--time-myr` targets must align exactly with the upstream run timeline; `--final` is the safer default when validating fixed-seed runs.
- The app currently runs tectonics through its adapter flow with app-provided heightmap input. Exact byte-identical replay against upstream tooling is not guaranteed unless adapter inputs are exported or upstream is run in the same seeded-heightmap mode.
