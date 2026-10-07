# Experiment ledger

Generated from `runs/ledger/*.jsonl` by `scripts/uwtools/ledger.py`; do not edit. Show entries with `uv run scripts/uw.py <module> ledger [--topic TEXT]`.

## Settled / don't retry

- `erosion-0014`: The rainfall solver's bulk-aerodynamic potential evaporation as the Budyko PET (discharge recall on Earth fell from 0.63 to 0.46).
- `erosion-0016`: Depth from the engine's crust age by the plate-cooling law: the engine's ocean ages are mostly the run's age (median 590 Myr at 600 Myr) and two thirds of the sea lies on continental crust, so 98 percent of the ocean lands at 5 to 6 km.
- `erosion-0021`: Loading fraction 2700/3300 with rebound 2800/3300 on rock volumes (it raised the land by 0.66 m on average in 20 Myr).
- `world-0003`: Re-enabling terrain.fastlem_mountains as it is: it does not terminate at seed 1, 256 columns.
- `world-0007`: Regenerating the access file from 256-column runs only: it drops writes that occur at the default size.
- `world-0035`: rewriting `accessor(...) > sealevel` without excluding an offset after sealevel or arithmetic before the accessor
- `world-0054`: converting the regional generator to float reads by pattern without launching the desktop application
- `world-0086`: Creating a junction from a worktree to the main tree's refs/ (or any data directory).

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

## erosion

| Id | Date | Verdict | Decision | What |
|---|---|---|---|---|
| `erosion-0001` | 2026-10-06 | integration runs; checks pass; look not yet judged | keep | erosion plan, step A4: hydrology.flow_routing and terrain.fluvial_erosion integrated; recipe eroded-world runs at the minimal configuration |
| `erosion-0002` | 2026-10-06 | pending user review | keep | erosion plan, step A6: the judging runs at 256 columns (run 1 at K = 24, run 2 at K = 100); run 1 is the module baseline pending the user's review |
| `erosion-0003` | 2026-10-06 | equivalent | keep | erosion plan, step B0: the Phase B contracts |
| `erosion-0004` | 2026-10-06 | validation in place; observations recorded | keep | erosion plan, step B5: the earth-hydrology profile and the river_mouths validation |
| `erosion-0005` | 2026-10-06 | rivers from the routing; procedural recipes identical | keep | erosion plan, step B1: hydrology.rivers from the routing in eroded-world |
| `erosion-0006` | 2026-10-06 | lake levels by water balance; procedural recipes identical | keep | erosion plan, step B2: lake levels by water balance in eroded-world |
| `erosion-0007` | 2026-10-06 | equivalent at the minimal configuration | keep | erosion plan, step B3: basin classes from the routing and the lake levels |
| `erosion-0008` | 2026-10-06 | consumers on the erosion pass; procedural recipes identical | keep | erosion plan, step B4: the surface potentials, the placers and the mineral exhumation read the erosion pass |
| `erosion-0009` | 2026-10-06 | pending user review | keep | erosion plan, step B6: the Phase B judging run (erosion run 3) and the close |
| `erosion-0010` | 2026-10-06 | accepted by the user | keep | erosion plan, after Phase B: the user's verdict on runs 1 and 3, the river-mouth overlay, the answers and the adversarial review prompt |
| `erosion-0011` | 2026-10-06 | run recorded; estimate wrong (11 min expected, 6.7 min measured) | keep | erosion plan, after Phase B: the 512-column run (erosion run 4), the Earth-reference assessment and the revised Phase C |
| `erosion-0012` | 2026-10-06 | equivalent | keep | erosion plan, step C0: the Phase C contracts |
| `erosion-0013` | 2026-10-06 | validation in place; the Caspian stays open until the runoff is fixed (C2) | keep | erosion plan, step C1: Earth lakes in earth-hydrology and the HydroSHEDS validation kinds |
| `erosion-0014` | 2026-10-06 | law in place; discharge comparison worse than C1 (attributed to the rainfall) | keep | erosion plan, step C2: Budyko runoff, transmission losses and winter baseflow |
| `erosion-0015` | 2026-10-06 | routing unbiased; straightness reduced; look pending (run 5) | keep | erosion plan, step C3: relief roughness, randomised receivers and flat resolution |
| `erosion-0016` | 2026-10-06 | cause fixed; calibration pending the user's acceptance | keep | erosion plan, step C4: the ocean floor's depth distribution |
| `erosion-0017` | 2026-10-06 | contract in place; other recipes identical | keep | erosion plan, step C6: recipe steps with labels and overrides; eroded-world runs a coarse climate before the erosion |
| `erosion-0018` | 2026-10-06 | flux in place; delta chance superseded by C5b | keep | erosion plan, step C5a: the erosion's sediment flux and the delta chance from it |
| `erosion-0019` | 2026-10-07 | equivalent | keep | erosion plan, step C7: the fluvial erosion joins the time model |
| `erosion-0020` | 2026-10-07 | mechanism in place; little land built at the current sediment supply | keep | erosion plan, step C5b: progradation at river mouths and distributaries on the new land |
| `erosion-0021` | 2026-10-07 | law in place; effect to judge at 256 columns | keep | erosion plan, step C8: local isostatic rebound during the erosion run |
| `erosion-0022` | 2026-10-07 | pending user review | keep | erosion plan, step C9: the Phase C judging run (erosion run 5) and the close |
| `erosion-0023` | 2026-10-07 | changed as intended; eroded-world only | keep | erosion plan, step D0: the default erodibility corrected from 24 to the derived 750 (unit error) |

Lessons:

- `erosion-0001`: A dispatch case that calls a helper which writes fields (getlandandseatotals) is an undeclared write of that process even when the audit happens to see no change; keep the stage body to the process and let the measured declaration cover it.
- `erosion-0001`: Inputs gathered once for an evolver that recomputes its own mask must be gathered for every cell (runoff over the sea), or cells that change class during the run get none.
- `erosion-0002`: A land-mean denudation at 20 km cells with kappa = 200 km2/Myr is diffusion-dominated: raising the stream power erodibility 4.2 times changed it by 6 percent. Check which term carries a rate before calibrating by it; the plan's band was written for the fluvial law.
- `erosion-0002`: A band criterion with a single lever needs the lever's sensitivity checked at the minimal configuration first (the A4 bundle showed 0.83 m/Myr at K = 24; one minimal run at K = 100 would have shown the insensitivity in 35 s instead of a 256-column run).
- `erosion-0003`: The plan's step list named one legacy stage per new process, but the legacy stages share state (basin_policy reruns the river map; lakes needs the river stage's scratch masks): a survey of what each legacy stage calls and needs found it before any worker started.
- `erosion-0004`: A fixed runoff coefficient of 0.4 on the climate module's precipitation makes desert rivers like the Nile carry many times their observed flow: the runoff needs the land evaporation the climate module computes but does not store (the evaporation-based runoff named for Phase B in configs/parameters/landscape.json, runoff note).
- `erosion-0005`: A stage that replaces a legacy writer must also write what the legacy writer left on cells outside its own domain (the legacy river map carried directions and flows into the sea, where the delta stage reads them); a consumer survey of the field, not of the stage, finds that.
- `erosion-0006`: A per-depression inflow taken from an accumulated discharge must exclude every cell downstream of the depression that still lies at its level, or the depression's own water returns as inflow; the routing's minimum-depth threshold creates exactly such cells.
- `erosion-0007`: A step whose effect needs a feature the minimal world lacks (closed lakes) shows identical bundles; say so in the expected field before the check, so that identical is a confirmation and not a surprise.
- `erosion-0008`: Cumulative erosion and deposition fields of a scheme with deposition are gross fluxes: any consumer that means 'lowered' or 'buried' needs their difference, and a review question that asks for the net and gross split finds it quickly.
- `erosion-0008`: A quantity defined with a time discount (the liberated gold's placer residence) must keep that discount when a new source adds to it, or the budget closes while the inventory is several times too large.
- `erosion-0009`: A tracked run at the judging size measures writes the minimal world never makes (salt lakes and pans here); one such run per phase close keeps the declarations complete.
- `erosion-0009`: An evaluation records its profile file's hash at evaluation time: changing the profile while a registered run is in progress would record a hash that does not match the maps it rendered, so profile edits wait until the run is evaluated.
- `erosion-0010`: D8 on a latitude-longitude grid prefers east and west away from the equator because those neighbours are closer: count the direction codes of a routing before reading its pattern as terrain.
- `erosion-0011`: A run-time estimate from a registered run of an older code version can be far off: the coupled rainfall fell from 287 s to 34 s at the same settings within a day; take timings from runs of the current code (run.json steps) before estimating.
- `erosion-0012`: A restart from a checkpoint needs every table the later stages read to be loadable from a bundle: the tectonic boundary segments are not, so the earliest useful checkpoint for terrain-stage experiments is after the last stage that reads them, not after the tectonics.
- `erosion-0013`: A whole-tail restart from a checkpoint is limited by the tables and solver fields a bundle cannot load: stopping the restart before the stages that need them (layers.finalize, climate.classification) still gives every earlier field bit for bit in a minute.
- `erosion-0013`: A lake raster that marks water does not give a lake level: the elevation import under it is a bed (ETOPO keeps lake bathymetry); a check that every lake surface lies at or above its terrain finds such a premise in one run.
- `erosion-0014`: A climate model's potential evaporation can be several times Earth's over dry land (a bulk formula with a hot dry surface and no surface resistance): check its land distribution against a temperature law before a Budyko curve reads it.
- `erosion-0014`: A runoff coefficient below one hides an excess of precipitation; replacing it with a physical law can worsen a discharge comparison while improving the model, so report the land water totals beside the comparison.
- `erosion-0015`: A random choice between two grid neighbours that matches the gradient in mean angle still drifts: a river's heading is its mean displacement, so the probabilities must balance the cross products with the gradient, and a test must sum displacements, not angles.
- `erosion-0015`: Noise octaves finer than the grid act as per-cell white noise whose slope grows with the width; band-limit them to twice the cell size so a parameter means the same at every width.
- `erosion-0016`: A random(a, b) that draws 15 bits makes every fractal with a range above 32768 one-sided: check a noise source's mean on its own before reading what a stage does with it.
- `erosion-0016`: A mask whose shapes are drawn in cells (warps, Voronoi spacing) covers a different share of the planet at every width; any field set from such a mask inherits the width dependence.
- `erosion-0017`: A stage that runs twice can read its own earlier output through a world field the first run left non-zero (the climate's SST): a repeated process needs a survey of what its stages read before they write, not only of what they reset.
- `erosion-0017`: A reference built on another step's tree does not show a merged step's effect: rerun the step check on the merged tree before the entry cites numbers.
- `erosion-0018`: A chance law with a reference value must be checked against the model's own distribution of the driver before it is chosen: every outlet here was below a tenth of the reference.
- `erosion-0019`: A checkpoint path that duplicates the end-of-run writer diverges as soon as a later step adds a field to one of them: share one writer before others edit the function.
- `erosion-0020`: A cap measured only by how many cells it turns into land needs an uncapped run of the same code: here every capped prograded cell was land without the cap too, so the cap's effect was elsewhere (heights, spills).
- `erosion-0020`: Access measured on a world where a field stays zero misses its writes: measure a new field's writer on the size where it changes.
- `erosion-0021`: An isostatic response with different fractions for removal and loading creates column whenever material only moves: check a mass-moving kernel's response for volume neutrality before choosing densities.
- `erosion-0022`: A phase whose steps each pass a minimal parity check can still need one tracked run at the judging size before the close: four of Phase C's access gaps showed only at 256 columns.
- `erosion-0023`: A derivation that matches two laws must carry the units of each quantity through: a precipitation of 1 m/yr is 1e-3 km/yr in a discharge counted in km3/yr, and the square root turned that into a factor of 31.6 that four phases of runs read as a calibration question.

## tectonics

| Id | Date | Verdict | Decision | What |
|---|---|---|---|---|
| `tectonics-0001` | 2026-10-04 | baseline | keep | consolidation phase 1: module baseline |
| `tectonics-0002` | 2026-10-04 | equivalent | keep | consolidation phase 1: tectonics into world_pipeline |
| `tectonics-0003` | 2026-10-05 | equivalent | keep | consolidation phase 3, step 8: first registered run in runs/tectonics/<id>/ |
| `tectonics-0004` | 2026-10-05 | equivalent except the announced table | keep | consolidation phase 3, step 9: snapshot statistics and comparison in Python, command-line checks as profile checks; the verifier's first review of a tectonics run |
| `tectonics-0005` | 2026-10-07 | equivalent | keep | erosion plan, step D2a: the landscape kernels on a torus metric (preparation of D2b, the engine's implicit rivers) |
| `tectonics-0006` | 2026-10-07 | equivalent with the option off; option in place | keep | erosion plan, step D2b: an opt-in implicit river solver in the tectonics engine (the landscape kernel on the engine's torus) |
| `tectonics-0007` | 2026-10-07 | changed as intended; contract extended (bundle only) | keep | erosion plan, step D1: the engine's rock uplift rate in the tectonic output contract, read by the erosion as U |

Lessons:

- `tectonics-0002`: The tectonics recipe uses the default planet body; procedural-world draws the body from the seed, so its tectonics step differs from the module run for seeds whose gravity is not 1.
- `tectonics-0002`: A world grid of 16 x 8 never leaves the shelf pass after tectonics (createfractal); world_pipeline run refuses fewer than 32 columns.
- `tectonics-0002`: Editing tectonics_report_data.py invalidates the Earth hypsography cache; refresh with tectonics_earth_reference.py --width 2048 --refresh.
- `tectonics-0002`: Mineral tables are written at the end of the run; checkpoints carry the tectonic fields and the snapshot record only.
- `tectonics-0002`: The snapshot record's summary, and the comparison code used by the engine tests, are still C++ (src/validation/tectonics).
- `tectonics-0004`: On procedural-world bundles the regime counts differ from the engine's own, because later terrain steps rewrite geologic_regime.
- `tectonics-0004`: The phase's reference bundles for later comparisons are out/scratch/phase3/step9-* (the snapshot table of before-* still has the summary).
- `tectonics-0004`: A reviewer cannot tell a plate split from colour reuse: the plate map colours the engine's plate index, which is not a persistent identity.
- `tectonics-0005`: A test that covers a metric change must exercise every code path that reads the changed geometry (here the randomised routing's offsets), and a conservation check cannot detect a face that is missing or doubled on both sides: a symmetry check can.
- `tectonics-0006`: A switch that changes a solver and a boundary condition at once needs a run that changes only the boundary condition before anything is attributed to the solver: here the base level carried almost all of the difference.
- `tectonics-0006`: A budget closed by rescaling is closed by construction: check the solver's own budget before the rescale, or the test of conservation tests nothing.
- `tectonics-0007`: A source term applied only to one class of cell, with the class recomputed every step, moves a boundary one way: apply it to both sides or bound it at the boundary, and count cells gained as well as lost.
- `tectonics-0007`: An exponential mean removes intermittent changes, not a steady advective term: an Eulerian rate over moving plates keeps -v.grad z.

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
| `world-0064` | 2026-10-06 | equivalent | keep | consolidation phase 6, step 9: climate_energy_tests.cpp split |
| `world-0065` | 2026-10-06 | equivalent | keep | consolidation phase 6, step 10: climate_physics_tests.cpp split |
| `world-0066` | 2026-10-06 | equivalent | keep | consolidation phase 6, step 11: climate_koppen_tests.cpp split |
| `world-0067` | 2026-10-06 | equivalent | keep | consolidation phase 6, step 12: core/climate_grid_tests.cpp split |
| `world-0068` | 2026-10-06 | equivalent | keep | consolidation phase 6, step 13: io/climate_reference_tests.cpp split |
| `world-0069` | 2026-10-06 | equivalent | keep | consolidation phase 6, step 14: phase close |
| `world-0070` | 2026-10-06 | equivalent | keep | climate rules on the float accessors (section 11, Later): the count before the move and the reference step check |
| `world-0071` | 2026-10-06 | equivalent | keep | climate rules on the float accessors, 1 of 3: climate_fields.cpp (module climate) |
| `world-0072` | 2026-10-06 | equivalent | keep | climate rules on the float accessors, 2 of 3: climate_classification.cpp (module climate) |
| `world-0073` | 2026-10-06 | equivalent | keep | climate rules on the float accessors, 3 of 3: climate_elevation.hpp (module climate; the whole-number form moves to src/app) |
| `world-0074` | 2026-10-06 | equivalent | keep | climate rules on the float accessors: close (plan section 11, AGENTS.md, verification guide) |
| `world-0075` | 2026-10-06 | equivalent | keep | rounding-getter calls outside src/app (section 11, Later): the count before the move and the reference step check |
| `world-0076` | 2026-10-06 | equivalent | keep | rounding-getter calls outside src/app, 1 of 3: the climate solver files (module climate) |
| `world-0077` | 2026-10-06 | equivalent | keep | rounding-getter calls outside src/app, 2 of 3: the lake queries of hydrology (module hydrology; the regional map rounds at its reads) |
| `world-0078` | 2026-10-06 | equivalent | keep | rounding-getter calls outside src/app, 3 of 3: the marine resources read the float SST (module resources) |
| `world-0079` | 2026-10-06 | equivalent | keep | rounding-getter calls outside src/app: close (plan section 11, AGENTS.md) |
| `world-0080` | 2026-10-06 | equivalent | keep | erosion plan, step 0: the FastLEM stages removed |
| `world-0081` | 2026-10-06 | documents only | keep | erosion plan, step 0.2: the plan, the agent workflow guide and the project agents |
| `world-0082` | 2026-10-06 | documents only | keep | erosion plan, step A0: the landscape kernel contracts |
| `world-0083` | 2026-10-06 | kernel tests pass | keep | erosion plan, step A1: sphere metric, flow graph and depressions kernels with their tests |
| `world-0084` | 2026-10-06 | equivalent | keep | erosion plan, step A3: catalogue fields, stages, the eroded-world recipe and the erosion evaluation module |
| `world-0085` | 2026-10-06 | kernel tests pass | keep | erosion plan, step A2: the stream power kernel with its tests |
| `world-0086` | 2026-10-07 | recovered except one hand-written guide | keep | incident during erosion plan Phase D: refs/processed deleted through a worktree junction, regenerated from refs/source; backup, guard hook and reference palettes |

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
- `world-0069`: Splitting a long test body is fastest when the sections' shared state is listed first (which variables, which config state each section saw) and the assertion literals are compared by script afterwards; reading the body once and writing all tests in one pass took 3 to 9 minutes per file.
- `world-0070`: step_check.py refuses a label whose directory already holds any <label>-* entry, including a log written there by hand; keep logs in a subdirectory.
- `world-0070`: A count of whole-number calls must strip comments: a commented-out block held 16 of the calls counted earlier.
- `world-0071`: climateworld now has basebelowsea; the nom() <= sealevel tests left in climate_winds.cpp, climate_ocean_currents.cpp, climate_sea_surface_temperature.cpp and climate_stages.cpp can use it when those files are next touched.
- `world-0073`: A template helper with int and float overloads hides who uses the rounding one: grep the callers by receiver and argument type before deleting an overload; here the only int caller was in src/app.
- `world-0073`: Reference bundles out/scratch/climate-float/elevation-* are the head state for the next step.
- `world-0074`: Before a float step, check which rules the step check reaches (mountain fields, sea cells, diagnostics switches in the reference bundles); when none is reached, say so before the move instead of predicting a changed world.
- `world-0079`: Counting which recipe reaches a changed line (field nonzero in the reference bundle, step present in the recipe) before the move predicted all three outcomes exactly, including the one field that would differ.
- `world-0080`: generateglobalclimate in src/pipeline/climate_generation.cpp has no caller: the recipes run the climate stages through world_stages.cpp.
- `world-0081`: A plan that will be executed by workers in parallel needs its shared contracts (headers, file ownership per step) written before the workers start; the plan makes that the orchestrator's first step (A0).
- `world-0083`: Steepest descent must be computed on the filled surface, not the terrain, or a rim cell beside a depression drains into it and the receiver chain cycles; the plan's section 3 said terrain.
- `world-0084`: A catalogue addition makes every bundle larger, so the next step check against an older reference reports the new fields as missing and exits 1 even when every existing field is identical; the first step check after A3 (A4's) uses an A3 bundle as its reference.
- `world-0084`: The palette file carries a version name that check_layout compares with the palette sidecars of registered runs: adding a palette entry needs a new version (tectonics-v5-erosion).
- `world-0085`: A budget that closes by construction (the exported flux recomputed from the final erosion and deposition) proves nothing about the deposition coefficients; the closed form needs its own one-cell test against the formula.
- `world-0085`: Float outputs resolve about 6e-8 relative: a plan tolerance of 1e-9 can only be checked on double accumulators.
- `world-0086`: A recursive delete follows a junction inside its target: never leave a link to irreplaceable data inside a directory that will be removed, check that a link removal succeeded (Test-Path) before the delete, and keep a copy of every untracked data directory.
- `world-0086`: When data is regenerated, look for any hash a run recorded (run.json inputs): one surviving hash turns 'probably the same' into 'identical'.
