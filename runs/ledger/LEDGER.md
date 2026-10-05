# Experiment ledger

Generated from `runs/ledger/*.jsonl` by `scripts/uwtools/ledger.py`; do not edit. Show entries with `uv run scripts/uw.py <module> ledger [--topic TEXT]`.

## Settled / don't retry

- `world-0003`: Re-enabling terrain.fastlem_mountains as it is: it does not terminate at seed 1, 256 columns.
- `world-0007`: Regenerating the access file from 256-column runs only: it drops writes that occur at the default size.
- `world-0035`: rewriting `accessor(...) > sealevel` without excluding an offset after sealevel or arithmetic before the accessor
- `world-0054`: converting the regional generator to float reads by pattern without launching the desktop application

## climate

| Id | Date | Verdict | Decision | What |
|---|---|---|---|---|
| `climate-0001` | 2026-10-04 | pending user review | keep | consolidation phase 1: the re-baseline |
| `climate-0002` | 2026-10-05 | parity shown | keep | consolidation phase 3: parity of global-balanced-v1 between the legacy chain and the bundle evaluation; repeat of run 385 |
| `climate-0003` | 2026-10-05 | equivalent | keep | consolidation phase 3, step 8: first registered run in runs/climate/<id>/; baseline 388 |

Lessons:

- `climate-0001`: Geosphere and land-hydrosphere fields keep integer types: every producer is integer (uint16 terrain imports, the engine's uint16 datum, legacy integer algorithms), so their type can change later without a new baseline. Category codes stay int16 for the same reason.
- `climate-0001`: atmosphere.precipitation and atmosphere.precipitation_exact are both float now and differ only where the legacy January/July passes adjust rain; the legacy-to-seasonal sync overwrites one and not the other.
- `climate-0001`: Accepted meaning changes of the exact readers: the transition-season fallback of syncseasonalclimatefromlegacy no longer truncates, and a stored SST between -0.5 and 0.5 C no longer counts as unset in pressuresurfacetemperature.
- `climate-0001`: Seen in pre-switch and baseline alike: a precipitation hotspot at the north pole near 87 E, a one-row line along the equator in January column water and convergence, staircase blocks in 850 hPa wind over high terrain, a sharp line along the Antarctic coast in surface wind.
- `climate-0001`: A run started with --from-bundle from a bundle written before this commit fails on the element type of the climate fields.
- `climate-0001`: bench.py climate run swallows --bundle (it is its own option), so a registered run has no bundle; the bundles here come from the same binaries and settings run by hand with the same surface drag input.
- `climate-0002`: A repeat of the climate baseline is bit-identical; any score difference between two runs of one binary comes from the options.
- `climate-0002`: --climate-transport mpdata takes about five times as long as ffsl at the registered settings and changes global-balanced-v1 by +0.025.
- `climate-0002`: weighted_relative_error (registry) and global-balanced-v1 can order a pair of runs differently.
- `climate-0003`: Registered climate runs no longer carry weighted_relative_error or the other C++ benchmark metrics; the registry record holds global-balanced-v1 only.

## tectonics

| Id | Date | Verdict | Decision | What |
|---|---|---|---|---|
| `tectonics-0001` | 2026-10-04 | baseline | keep | consolidation phase 1: module baseline |
| `tectonics-0002` | 2026-10-04 | equivalent | keep | consolidation phase 1: tectonics into world_pipeline |
| `tectonics-0003` | 2026-10-05 | equivalent | keep | consolidation phase 3, step 8: first registered run in runs/tectonics/<id>/ |
| `tectonics-0004` | 2026-10-05 | equivalent except the announced table | keep | consolidation phase 3, step 9: snapshot statistics and comparison in Python, command-line checks as profile checks; the verifier's first review of a tectonics run |

Lessons:

- `tectonics-0002`: The tectonics recipe uses the default planet body; procedural-world draws the body from the seed, so its tectonics step differs from the module run for seeds whose gravity is not 1.
- `tectonics-0002`: A world grid of 16 x 8 never leaves the shelf pass after tectonics (createfractal); world_pipeline run refuses fewer than 32 columns.
- `tectonics-0002`: Editing tectonics_report_data.py invalidates the Earth hypsography cache; refresh with tectonics_earth_reference.py --width 2048 --refresh.
- `tectonics-0002`: Mineral tables are written at the end of the run; checkpoints carry the tectonic fields and the snapshot record only.
- `tectonics-0002`: The snapshot record's summary, and the comparison code used by the engine tests, are still C++ (src/validation/tectonics).
- `tectonics-0004`: On procedural-world bundles the regime counts differ from the engine's own, because later terrain steps rewrite geologic_regime.
- `tectonics-0004`: The phase's reference bundles for later comparisons are out/scratch/phase3/step9-* (the snapshot table of before-* still has the summary).
- `tectonics-0004`: A reviewer cannot tell a plate split from colour reuse: the plate map colours the engine's plate index, which is not a persistent identity.

## world

| Id | Date | Verdict | Decision | What |
|---|---|---|---|---|
| `world-0001` | 2026-10-03 | equivalent | keep | consolidation phase 0 |
| `world-0002` | 2026-10-04 | equivalent | keep | consolidation phase 1, step 1: CLI options and the climate_pipeline fold |
| `world-0003` | 2026-10-04 | changed | keep | consolidation phase 1: FastLEM passes switched off |
| `world-0004` | 2026-10-04 | equivalent with integer producers; changed for procedural worlds | keep | consolidation phase 1: base elevation as float32 |
| `world-0005` | 2026-10-04 | tests repaired | keep | consolidation phase 2, step 1: first test run since before Phase 0 and the fixes |
| `world-0006` | 2026-10-04 | equivalent | keep | consolidation phase 2, step 2: three GoogleTest binaries |
| `world-0007` | 2026-10-04 | equivalent | keep | consolidation phase 2, step 3: process access regenerated |
| `world-0008` | 2026-10-04 | equivalent | keep | consolidation phase 2: review of the test conversion |
| `world-0009` | 2026-10-04 | equivalent | keep | consolidation phase 3, step 1: shared Python helpers in scripts/uwtools |
| `world-0010` | 2026-10-04 | equivalent | keep | consolidation phase 3, step 2: pyproject.toml, scripts/uw.py, snake_case script names, audit exit status |
| `world-0011` | 2026-10-04 | equivalent | keep | consolidation phase 3, step 3: report templates |
| `world-0012` | 2026-10-04 | tooling added | keep | consolidation phase 3, step 4: ledger tooling |
| `world-0013` | 2026-10-04 | tooling added | keep | consolidation phase 3, step 5: verifier agent, review checklists, reported checks |
| `world-0014` | 2026-10-04 | equivalent | keep | consolidation phase 3, step 6: configs/regions.json and the regional table |
| `world-0015` | 2026-10-05 | equivalent | keep | consolidation phase 3, step 7: observed references as .npy with a manifest per dataset |
| `world-0016` | 2026-10-05 | equivalent | keep | consolidation phase 3, step 12 (done before steps 8 to 11): Excel automation removed |
| `world-0017` | 2026-10-05 | equivalent | keep | consolidation phase 3, steps 8 and 11: run layout runs/<module>/<id>/, old layout sealed, storage script ported |
| `world-0018` | 2026-10-05 | retired | keep | consolidation phase 3, step 13: legacy climate chain retired (step 9 is tectonics-0004) |
| `world-0019` | 2026-10-05 | tooling added | keep | consolidation phase 3, step 10: palette sidecars and their check |
| `world-0020` | 2026-10-05 | defects fixed | keep | consolidation phase 3: fresh-context review of the new tooling and its fixes |
| `world-0021` | 2026-10-05 | equivalent | keep | consolidation phase 4, step 1: directory moves (src/world, src/cli, profiler and environment in src/core, workbench in src/app); the step list of the phase |
| `world-0022` | 2026-10-05 | equivalent | keep | consolidation phase 4, step 2: per-directory CMake, one static library per directory, engine target uw_tectonics |
| `world-0023` | 2026-10-05 | equivalent | keep | consolidation phase 4, step 3: one include root (src/) and the layering rules in check_layout.py |
| `world-0024` | 2026-10-05 | tooling added | keep | consolidation phase 4, step 4: review maps a verifier can read (legend, graticule, direction, parent, plate identity); cross-imports of scripts/evaluation |
| `world-0025` | 2026-10-05 | equivalent | keep | consolidation phase 4, step 5: assetdata.cpp split, template data in assets/terrain |
| `world-0026` | 2026-10-05 | equivalent | keep | consolidation phase 4, step 6: climate_validation.cpp split; the C++ climate benchmark, its maps and its reference comparisons removed |
| `world-0027` | 2026-10-05 | equivalent | keep | consolidation phase 4, step 7: climate_coupling.cpp split by process; test modes as disabled tests; callerless field readers removed |
| `world-0028` | 2026-10-05 | equivalent | keep | consolidation phase 4, step 8: main.cpp split by screen |
| `world-0029` | 2026-10-05 | equivalent | keep | consolidation phase 4, step 9: misc.cpp dissolved into owner sources; fractal noise on small grids; planet::smooth() |
| `world-0030` | 2026-10-05 | tooling added | keep | consolidation phase 4, step 9b: validation guards (features, variance ratio, footprint); the tectonics report through the shared templates |
| `world-0031` | 2026-10-05 | pending user review | keep | consolidation phase 4, step 10: the world series; run 1 is the reference before the float steps |
| `world-0032` | 2026-10-05 | equivalent on Earth terrain; procedural worlds shifted as listed | keep | consolidation phase 4, step 11: atmosphere.precipitation and atmosphere.precipitation_exact become one field (module climate) |
| `world-0033` | 2026-10-05 | pending user review | keep | consolidation phase 4, step 11b (added at the user's request): minimal run configuration; climate options kept across a resize; world run 2 is the minimal reference |
| `world-0034` | 2026-10-05 | pending user review | keep | consolidation phase 4, step 12: terrain on the float accessors (module terrain) |
| `world-0035` | 2026-10-05 | pending user review | keep | consolidation phase 4, step 12 corrected: a comparison the float pass had rewritten wrongly (module terrain) |
| `world-0036` | 2026-10-05 | pending user review | keep | consolidation phase 4, step 12c: the six elevation components stored as float (module terrain) |
| `world-0037` | 2026-10-05 | pending user review | keep | consolidation phase 4, step 13: hydrology on the float accessors; lake surfaces stored as float (module hydrology) |
| `world-0038` | 2026-10-05 | pending user review | keep | consolidation phase 4, step 14: resources on the float accessors (module resources) |
| `world-0039` | 2026-10-05 | pending user review | keep | consolidation phase 4, step 15: society on the float accessors (module society) |
| `world-0040` | 2026-10-05 | pending user review | keep | consolidation phase 4, step 16: phase close |
| `world-0041` | 2026-10-05 | equivalent | keep | consolidation phase 5, step 0: step list and the step check |
| `world-0042` | 2026-10-05 | pending user review | keep | consolidation phase 5, step 1: screens as member functions of an application state |
| `world-0043` | 2026-10-05 | equivalent | keep | consolidation phase 5, step 5: the regional map leaves the headless link line |
| `world-0044` | 2026-10-05 | equivalent | keep | consolidation phase 5, steps 6 and 7: parameter mechanism; terrain parameters |
| `world-0045` | 2026-10-05 | equivalent | keep | consolidation phase 5, step 8: hydrology parameters |
| `world-0046` | 2026-10-05 | equivalent | keep | consolidation phase 5, step 9: climate parameters |
| `world-0047` | 2026-10-05 | equivalent | keep | consolidation phase 5, step 10: world defaults and the planet body as parameters; generation_tuning.hpp removed |
| `world-0048` | 2026-10-05 | equivalent | keep | consolidation phase 5, step 10b: snake_case keys in configs/tectonics-baseline.json |
| `world-0049` | 2026-10-05 | equivalent | keep | consolidation phase 5, end of the first session: status note and open steps |
| `world-0050` | 2026-10-05 | pending user review | keep | consolidation phase 5, step 12: procedural worlds keep the default body |
| `world-0051` | 2026-10-05 | pending user review | keep | consolidation phase 5, step 4a: the stage bodies of src/pipeline on the float accessors (module terrain) |
| `world-0052` | 2026-10-05 | pending user review | keep | consolidation phase 5, step 2: the desktop application loads a run bundle (compile-only) |
| `world-0053` | 2026-10-05 | pending user review | keep | consolidation phase 5, step 3: the Generate elements buttons call a process or a terrain function (compile-only) |
| `world-0054` | 2026-10-05 | pending user review | keep | consolidation phase 5, step 11: phase close |
| `world-0055` | 2026-10-06 | equivalent | keep | consolidation phase 6, step 0: baseline test run and the split list |
| `world-0056` | 2026-10-06 | equivalent | keep | consolidation phase 6, step 1: climate_atmosphere_tests.cpp split |
| `world-0057` | 2026-10-06 | equivalent | keep | consolidation phase 6, step 2: climate_ocean_dynamics_tests.cpp split |
| `world-0058` | 2026-10-06 | equivalent | keep | consolidation phase 6, step 3: climate_hydrology_tests.cpp split |
| `world-0059` | 2026-10-06 | equivalent | keep | consolidation phase 6, step 4: climate_ocean_heat_tests.cpp split |
| `world-0060` | 2026-10-06 | equivalent | keep | consolidation phase 6, step 5: climate_transport_tests.cpp split |
| `world-0061` | 2026-10-06 | equivalent | keep | consolidation phase 6, step 6: climate_boundary_layer_tests.cpp split |
| `world-0062` | 2026-10-06 | equivalent | keep | consolidation phase 6, step 7: climate_weather_tests.cpp split |
| `world-0063` | 2026-10-06 | equivalent | keep | consolidation phase 6, step 8: climate_flow_tests.cpp split |

Lessons:

- `world-0001`: A full seed-1 procedural-world run takes 38 min at 20 threads and 57 min at 4 on this machine; background runs need a time limit above 30 min.
- `world-0001`: Plan section 6 was wrong about refs/processed/climate: csv/ and fields/ are not empty; maps/climate/ and previews/ have live consumers (corrected in the plan).
- `world-0001`: scripts/benchmarks/tectonics_climate_feedback.py still copies tectonic_climate_erosion.exe, whose target was removed: it cannot run until a terrain.climate_erosion evolver exists or the script is retired.
- `world-0002`: A procedural world cannot yet be finished from a bundle: layers.finalize reads the mineral tables, rivers need the mountain drainage and lakes the salt-lake masks. These are the GenerationScratch and table items of world-model migration step 5.
- `world-0002`: world_pipeline at seed 1 and 256 columns never leaves the FastLEM mountains step (commit 5cf4bf1 as well): small-resolution world runs are not usable for quick checks until that is fixed.
- `world-0002`: Checkpoints here are per recipe step. Checkpoints inside an evolver (--duration, --checkpoint-every in Myr) come with the tectonics recipe.
- `world-0002`: Timings at 20 threads: seed-1 procedural world 36 min; Earth climate at 512 with default settings about 15 min; climate_pipeline benchmark is gone, so no earth-climate run can be registered until the re-baseline step.
- `world-0003`: At 256 columns the cost is tectonics (6 min for 600 Myr) and the coupled climate step (5.5 min even at a 64-wide atmosphere with two coupling updates); terrain, hydrology and layers are negligible.
- `world-0003`: The FastLEM code stays as legacy (islegacygenerationstage) until it is repaired or replaced by another erosion and drainage model.
- `world-0004`: The land/sea mask is taken from the whole-metre elevation so that integer and float readers agree on the coast.
- `world-0004`: planet::smooth() ignores the grid it is given and smooths the base elevation; kept as it was.
- `world-0004`: Mountain height, canyon uplift, crater and volcano heights, ocean ridge height and lake surfaces are still integer fields.
- `world-0005`: The Python suite needs more packages than verification.md listed; a module that cannot import shows as one error and hides its tests.
- `world-0005`: Ten Python tests are skipped by their own conditions; the reasons were not examined.
- `world-0006`: Most climate files are still one long TEST each (Model): they were one main() and share local fixtures. Splitting them is worthwhile only when a file is next edited.
- `world-0006`: ctest names are Suite.Name; -R '^ClimateHydrology\.' selects a suite and -L tectonics or -L geology a binary.
- `world-0006`: Old test executables stay in bin/<Config> of an existing build directory until they are deleted by hand.
- `world-0007`: A write is measured as a changed hash: a step that writes the value a field already has is not seen, and a small world leaves out whole code paths. Reports need the default size.
- `world-0007`: Timings of the tracking build at 20 threads, seed 1: tectonics under 25 min, earth-climate 16 min, imported-world 16 min, procedural-world 37 min.
- `world-0007`: process_access.py writes the report names into the file, so --check needs the same reports in the same order; the stored directory in name order is the reference.
- `world-0007`: A run with --audit-access exits 0 after reporting an undeclared write in its access report (seen in the default-size procedural run); read the report or the log, not the exit status.
- `world-0009`: The shared .uwclim reader checks the payload size; a file with trailing bytes, which the memmap readers accepted, is now rejected.
- `world-0009`: tectonics_experiment and the palette loaders read UTF-8 instead of the locale encoding.
- `world-0010`: The seed-1 256-column reference bundles of the three recipes for this phase are out/scratch/phase3/before-*; the three audited runs take about 13 minutes together.
- `world-0010`: runs/registry records keep the former script paths in their stored commands; they are history.
- `world-0011`: The climate registry holds no global-balanced-v1 score for runs 379 to 385, so the dashboard's score column is empty for them.
- `world-0011`: The elevation animation of a 61-snapshot tectonics run is 5.3 MB.
- `world-0013`: The review maps have no legend and no graticule, and most are 256 x 128 or smaller (solver grids); the verifier read values by matching colours to the palette file.
- `world-0013`: Wind and transport maps show speed only, so checklist entries about direction (trades, westerlies, monsoon reversal) cannot be answered from them.
- `world-0013`: The request carries the parent's losses but no parent maps; better and worse against the parent is judged from numbers only.
- `world-0013`: A checklist 'good' text must agree with the reference maps: the first ITCZ text contradicted ERA5 in the eastern Pacific.
- `world-0014`: evaluation.json grows from about 0.3 MB to 0.9 MB with the regional table.
- `world-0014`: AR6 land regions are not restricted to land cells in the regional table: each entry's own mask applies.
- `world-0015`: A .npy file carries no variable name: the C++ reader checks that the file name equals the expected variable or starts with it and an underscore.
- `world-0015`: The sealed archive runs/reports/climate-reference-experiments-20260915/prepared-inputs still holds .uwclim inputs that the new binary cannot read; replaying those experiments means preparing the inputs again.
- `world-0015`: The moved receipts still name the former .uwclim files; fields are found by dataset and variable.
- `world-0015`: refs/source/physical and refs/processed/physical are junctions into D:/dev/undiscovered-worlds/refs/ since step 1.
- `world-0016`: Stored commands that pass --no-profile-log (registry records, old run manifests) are rejected by the current binary.
- `world-0017`: A registered tectonics run is 1.4 GB (60 checkpoint bundles 1.1 GB, evaluation with animations and analysis 0.25 GB); a registered climate run is 0.5 GB.
- `world-0017`: The registry files are stored with CRLF in the working tree and LF in git; records are inserted as text.
- `world-0017`: scripts/evaluation still imports map_rendering, climate_palettes and tectonics_earth_reference from scripts/benchmarks and cell_grid from scripts/refs.
- `world-0017`: Old-layout paths are still named in map_projection.py, summarize_climate_validation.py, scripts/refs/organize_climate_maps.py and the map guides under docs/reference-guides.
- `world-0018`: scripts/benchmarks/generate_climate_palette_header.py runs when it is imported (no main guard) and rewrites src/validation/climate/climate_palette.hpp; an import check changed the file's line endings until it was restored.
- `world-0018`: pyshp is declared in pyproject.toml and imported nowhere.
- `world-0018`: scripts/evaluation still imports map_rendering, climate_palettes and tectonics_earth_reference from scripts/benchmarks.
- `world-0019`: A palette change needs a new version string in the palette file: the same version with another hash makes every existing sidecar of that version fail the layout check.
- `world-0020`: A check that compares two bundles must be shown failing once: recipe-equality passed for a week of design because both commands were the same.
- `world-0020`: Sealed directories need a guard in code, not only a README: a default output path beside the bundle pointed into the archive.
- `world-0021`: out/scratch/phase4/step.ps1 <label> <reference-label> builds both presets, makes the three bundles, compares them and runs the layout and access checks; the three audited runs and both builds take about 20 minutes.
- `world-0021`: The reference bundles of the next step are out/scratch/phase4/step01-*.
- `world-0021`: out/scratch/phase4/include_graph.py prints which source directory includes which.
- `world-0022`: The reference bundles of the next step are out/scratch/phase4/step02-*.
- `world-0022`: PNG is now found for every configuration, because uw_io always contains the heightmap I/O of the engine tests.
- `world-0023`: A quoted include is searched beside the including file first, so the compiler accepts a bare sibling name; only check_layout.py keeps the include root.
- `world-0023`: The reference bundles of the next compiled step are out/scratch/phase4/step03-*.
- `world-0024`: Pillow's built-in font has no minus sign or en dash; map_rendering replaces them in drawn text.
- `world-0024`: A legend bar that is linear in value crowds the ticks of the elevation scale near sea level; labels that would overlap are left out, the ticks stay.
- `world-0024`: A Bash heredoc that contains an apostrophe fails before anything runs; scripts are written with the Write tool.
- `world-0025`: The reference bundles of the next step are out/scratch/phase4/step05-*.
- `world-0025`: A single source can be compiled against the tree for a scratch check with ~/bin/msvc-env.cmd cl /std:c++17 /EHsc /I src (out/scratch/phase4/templates-check).
- `world-0026`: The reference bundles of the next step are out/scratch/phase4/step06-*.
- `world-0026`: src/pipeline/world_processes.cpp still describes a process as writing 'the climate validation report'; the text is corrected in step 7, which edits the climate processes.
- `world-0026`: scripts/refs/reference_map_rendering.py says its palettes follow climate_circulation_diagnostics.cpp, which no longer exists.
- `world-0027`: The reference bundles of the next step are out/scratch/phase4/step07-*.
- `world-0027`: src/simulations/climate/climate_context.hpp still names climate_coupling.cpp and a function solvegenericlevel in a comment; no function of that name exists.
- `world-0028`: A second, smaller block of the export area screen (the selection rectangle) stays in main.cpp, beside the drawing of the other screens.
- `world-0029`: The reference bundles of the next compiled step are out/scratch/phase4/step09-*.
- `world-0029`: The grid size macros ARRAYWIDTH and ARRAYHEIGHT are defined identically in six files; grid_shapes.cpp is one more until they become one constant.
- `world-0029`: A function moved upwards in planet.cpp needed a helper template defined further down: the build, not the script, found it.
- `world-0030`: runs/reports/tectonics/76/evaluation now exists (written by the check of this step).
- `world-0030`: The monsoon reversal criterion (120 degrees, 1 m/s in both seasons) gives 24 % of the tropics for ERA5 10 m winds, more than the monsoon regions proper: it also counts cells where weak winds turn.
- `world-0031`: A 512-column procedural world takes 11 minutes; the tectonics step takes the same 6 minutes at 256 and 512 columns.
- `world-0031`: The verifier traced river courses and measured boundaries from the bundle fields with its own scripts (out/scratch/phase4/verifier/world-1); the maps alone would not have answered rivers-reach-water.
- `world-0031`: Lake surfaces are stored in metres above the engine datum (sea level is about 31,000), not above sea level.
- `world-0032`: --climate-coupling-iterations has no effect in a procedural world: climateworld::resize() resets the climate options, and the coupled loop ran 21 iterations with a limit of 2. Every 'climate 64, 2 iterations' procedural run of this phase, including world run 1, ran to convergence. Step 11b fixes it.
- `world-0032`: Of the ten access runs only the default-size procedural world contributes declarations that no other run does.
- `world-0032`: runs/reports/phase4-access holds the report set and the script that produced nine of its runs.
- `world-0033`: The reference bundles of the next step are out/scratch/phase4/step11c-* (minimal configuration).
- `world-0033`: The tectonics engine inside procedural-world runs on its own 256 x 128 grid whatever the world size: 50 Myr take 23 s there and 6 s in the tectonics recipe on 128 x 64.
- `world-0033`: Of the ten access runs of Phase 2, nine repeated recipes at a second size; only the default-size procedural world added declarations.
- `world-0033`: A scratch climate run at the registered size takes 6.5 minutes; on 128 x 64 with the climate on 32 columns it takes 9 seconds and the evaluation remaps the references itself.
- `world-0034`: The reference bundles of the next step are out/scratch/phase4/step12a-* (minimal configuration).
- `world-0034`: The ordinary build prints no conversion warnings; compile a source with /W3 (truncations.ps1) to see every float squeezed into an integer.
- `world-0034`: A regex written into a script through a shell here-document lost its backslash-b: scripts are written with the Write tool.
- `world-0035`: The reference bundles of the next step are out/scratch/phase4/step12b-* (minimal configuration); world run 4 is the parent of the hydrology run.
- `world-0035`: The river-flow palette (0 to 200,000) suits the 512-column world and shows nothing on 128 x 64, where the largest flow is about 2,000.
- `world-0035`: A rewrite by pattern needs a search for what it must never produce, run every time.
- `world-0036`: The minimal world (128 x 64, 50 Myr) has zero mountain height, crater and ocean ridge fields: a parity run at that size says nothing about code that only acts on them.
- `world-0036`: The reference bundles of the next step are out/scratch/phase4/step12c-*; world run 5 is the parent of the hydrology run.
- `world-0037`: A float move in a module that routes water changes the routing wherever whole metres were tied; expect a different world, not a rounding-sized difference.
- `world-0037`: The minimal world has no lakes and, before this step, no deltas; it cannot show a defect in either.
- `world-0037`: The reference bundles of the next step are out/scratch/phase4/step13-*; world run 6 is the parent.
- `world-0038`: The reference bundles of the next step are out/scratch/phase4/step14-*; world run 7 is the parent.
- `world-0039`: Society runs only with --social; no recipe default and no registered run turns it on. A step that touches society needs its own pair of runs with the switch.
- `world-0039`: The reference bundles after this step are out/scratch/phase4/step15-*; world run 8 is the parent.
- `world-0040`: A change of number type did not change any declared read or write: the rule of measuring access only when reads or writes change held.
- `world-0040`: The estimate was wrong by a factor of five because it counted reading every call site by hand; a rename script, a compile with conversion warnings and reading only the flagged lines did the same work.
- `world-0041`: Reference bundles of Phase 5 are out/scratch/phase5/<label>-<recipe>; step_check.py LABEL --reference PREVIOUS is the whole step check.
- `world-0042`: Member functions named after the former lambdas let fifteen fragments compile unchanged except where a lambda was passed as a value (two places).
- `world-0043`: Reference bundles of the next step: out/scratch/phase5/step05-*.
- `world-0044`: Keeping the C++ names and making the values variables moves a module in minutes; only uses that need a compile-time constant (an array size here) have to be edited.
- `world-0044`: Reference bundles of the next step: out/scratch/phase5/step06-*.
- `world-0045`: Reference bundles of the next step: out/scratch/phase5/step08-*.
- `world-0046`: The first registered climate run after this step should compare its solver times with run 388 before any tuning starts.
- `world-0046`: Reference bundles of the next step: out/scratch/phase5/step09-*.
- `world-0047`: Reference bundles of the next step: out/scratch/phase5/step10-*.
- `world-0048`: Run records written before this step hold the dashed scenario keys; a script that reads a stored recipe must accept both.
- `world-0049`: Whole-number reads were still in src/pipeline after Phase 4: its float steps went through src/simulations module by module and the stage bodies written directly in world_stages.cpp and climate_generation.cpp were not a module.
- `world-0049`: A survey by a subagent of who calls what (regional map, tuning constants, GUI buttons) took four minutes and decided the order of the steps; the conversions themselves were scripts over line ranges and names.
- `world-0050`: World run 9 is the parent of the next world-changing step; reference bundles out/scratch/phase5/step12-*.
- `world-0051`: A step-check label is used once: world_pipeline refuses a non-empty bundle directory, and the comparison then reads the old bundles.
- `world-0051`: World run 10 is the parent of the next world-changing step; reference bundles out/scratch/phase5/step4a2-*.
- `world-0054`: Phase 6 reorganises tests that have not been run since Phase 4 step 7: run ctest once before touching them, so that a failure caused by Phases 4 and 5 is not blamed on the split.
- `world-0054`: A claim that a module is on the float accessors needs a count of the remaining whole-number calls by receiver type (planet and climateworld have getters of the same names), not a search for one spelling.
- `world-0055`: A test file whose sections mutate a shared config between them (ocean dynamics: oneWayConfig in three states; weather: the stochastic forcing set midway) needs the config reproduced at the state each section saw, or the split changes what is tested.
- `world-0056`: GoogleTest refuses one suite with both TEST and TEST_F; naming the fixture class after the suite and converting the existing one-liners keeps suite names, ctest -R patterns and the Modes test names stable.
- `world-0056`: A step's duration is measured commit to commit (git log --format='%h %ci').
- `world-0057`: When a split recomputes a config state, write down which of the old body's states each section saw (initial, closed exchange, daily freezing steps) before writing the tests; the replay sections and the bottom-stress sections would otherwise have taken the wrong one.
