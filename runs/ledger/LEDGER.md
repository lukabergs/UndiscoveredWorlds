# Experiment ledger

Generated from `runs/ledger/*.jsonl` by `scripts/uwtools/ledger.py`; do not edit. Show entries with `uv run scripts/uw.py <module> ledger [--topic TEXT]`.

## Settled / don't retry

- `world-0003`: Re-enabling terrain.fastlem_mountains as it is: it does not terminate at seed 1, 256 columns.
- `world-0007`: Regenerating the access file from 256-column runs only: it drops writes that occur at the default size.

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
