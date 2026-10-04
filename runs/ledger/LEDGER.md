# Experiment ledger

Generated from `runs/ledger/*.jsonl` by `scripts/uwtools/ledger.py`; do not edit. Show entries with `uv run scripts/uw.py <module> ledger [--topic TEXT]`.

## Settled / don't retry

- `world-0003`: Re-enabling terrain.fastlem_mountains as it is: it does not terminate at seed 1, 256 columns.
- `world-0007`: Regenerating the access file from 256-column runs only: it drops writes that occur at the default size.

## climate

| Id | Date | Verdict | Decision | What |
|---|---|---|---|---|
| `climate-0001` | 2026-10-04 | pending user review | keep | consolidation phase 1: the re-baseline |

Lessons:

- `climate-0001`: Geosphere and land-hydrosphere fields keep integer types: every producer is integer (uint16 terrain imports, the engine's uint16 datum, legacy integer algorithms), so their type can change later without a new baseline. Category codes stay int16 for the same reason.
- `climate-0001`: atmosphere.precipitation and atmosphere.precipitation_exact are both float now and differ only where the legacy January/July passes adjust rain; the legacy-to-seasonal sync overwrites one and not the other.
- `climate-0001`: Accepted meaning changes of the exact readers: the transition-season fallback of syncseasonalclimatefromlegacy no longer truncates, and a stored SST between -0.5 and 0.5 C no longer counts as unset in pressuresurfacetemperature.
- `climate-0001`: Seen in pre-switch and baseline alike: a precipitation hotspot at the north pole near 87 E, a one-row line along the equator in January column water and convergence, staircase blocks in 850 hPa wind over high terrain, a sharp line along the Antarctic coast in surface wind.
- `climate-0001`: A run started with --from-bundle from a bundle written before this commit fails on the element type of the climate fields.
- `climate-0001`: bench.py climate run swallows --bundle (it is its own option), so a registered run has no bundle; the bundles here come from the same binaries and settings run by hand with the same surface drag input.

## tectonics

| Id | Date | Verdict | Decision | What |
|---|---|---|---|---|
| `tectonics-0001` | 2026-10-04 | baseline | keep | consolidation phase 1: module baseline |
| `tectonics-0002` | 2026-10-04 | equivalent | keep | consolidation phase 1: tectonics into world_pipeline |

Lessons:

- `tectonics-0002`: The tectonics recipe uses the default planet body; procedural-world draws the body from the seed, so its tectonics step differs from the module run for seeds whose gravity is not 1.
- `tectonics-0002`: A world grid of 16 x 8 never leaves the shelf pass after tectonics (createfractal); world_pipeline run refuses fewer than 32 columns.
- `tectonics-0002`: Editing tectonics_report_data.py invalidates the Earth hypsography cache; refresh with tectonics_earth_reference.py --width 2048 --refresh.
- `tectonics-0002`: Mineral tables are written at the end of the run; checkpoints carry the tectonic fields and the snapshot record only.
- `tectonics-0002`: The snapshot record's summary, and the comparison code used by the engine tests, are still C++ (src/validation/tectonics).

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
