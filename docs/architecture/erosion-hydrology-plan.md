# Erosion and hydrology plan

Written 2026-10-06 by the lead session as the brief for the sessions that give the simulator an erosion and drainage model. It follows `world-model.md` (components, process kinds, recipes, evaluation) and `consolidation-plan.md` (structure, tooling, working agreement, naming); the roles and checks of the sessions are `docs/reference-guides/agent-workflow.md`. Section 2 records why each candidate algorithm was chosen or rejected, so that later agents do not re-run the comparison. Decisions marked **user** need the user's answer before the step that depends on them; everything else is decided here and undone by a commit.

## 1. Decisions taken on 2026-10-06

- **FastLEM is deleted** (Phase 0). `terrain.fastlem_mountains` solved the Salève steady-state formula on sites about 30 world cells apart with invented uplift and invented rainfall (latitude and coast distance), before any climate stage, and threw the solved elevations away: only ridge lines along basin divides reached the template mountain builder. It was a ridge-mask generator, not an erosion pass; it hung at seed 1, 256 columns (`world-0003`), and had been in no recipe since 2026-10-04. `terrain.post_river_fastlem` raised the base elevation near rivers by distance rules. Neither idea survives; the one reusable idea, the analytical stream-power solution, returns from its papers in Phase E. The code is in git before the removal commit.
- **One law, three forms.** Erosion and drainage rest on the stream power law (erosion rate = erodibility x discharge^m x slope^n, solved along the drainage tree from the sea upwards), in a different numerical form per scale: the **global pass** is an implicit time-stepped solve with sediment deposition and lake routing on the lat-lon sphere (the FastScape and goSPL design, §2); the **regional pass** is the analytical steady solve of Steer 2021 as extended by Tzathas et al. 2024, started from the upsampled global terrain with the global rivers pinned (Phase E); **render detail** is a point-evaluated gully filter in the renderer only (Phase F). Nothing else models erosion.
- **Flow routing is one shared diagnoser.** `hydrology.flow_routing` (receivers, discharge, drainage area, depressions and their spill points) replaces, over time, the four routings that exist today (the engine's river step, `climate_erosion.cpp`, `createrivermap`, the placer routing). Rivers, lakes, basins, placers, exhumation and fertility read its fields.
- **Inputs come from the pipeline, not from new generators.** Runoff is the climate module's precipitation (seasonal `atmosphere.precipitation`, reduced by snow and ice), erodibility is a proxy from the tectonic fields (`crust_class`, `crust_age`), base level is `geosphere.sea_level`. Uplift during the erosion run is zero in Phase A (a post-tectonic denudation of a set duration): the engine's `uplift_tendency` is a 0..1 score, not a rate, and a rate is a change of the tectonic output contract (**user**, Phase D).
- **Coexistence, then promotion.** The new processes run in a new recipe `eroded-world`; `procedural-world` is unchanged until the user's verdict on the first judging run (Phase A close) promotes the eroded recipe. Each phase ends with one 256-column run for the user to look at and nothing larger.
- **No external library.** fastscapelib is the only linkable candidate (GPL-3.0, C++17, header-only) and is early-stage, without deposition, with its spherical grid disabled on Windows; the kernels are written here from the papers (a few hundred lines each) in a planet-free library `uw_landscape` so they can be unit-tested on synthetic grids.
- **Elevation.** The global pass erodes `geosphere.base_elevation` (the surface of procedural worlds now that no recipe adds template mountains) and leaves the six added components alone; an imported world with template mountains is out of scope for `eroded-world`. `float32` storage, `double` accumulation.
- **Planet radius.** The world has no radius (`size_class` 1..3 only; the engine has `world_width_km`). The sphere metric takes `planetbody.radius` (new scalar field, m; default from the parameter `world.body.radius_km` = 6371, scaled by size class later if wanted). Recorded here so no agent asks.
- **Module and names.** Kernels: `src/simulations/landscape/` (library `uw_landscape`, depends on `uw_core` only; parameters `configs/parameters/landscape.json`, `tuning::landscape::*`). Processes: `src/simulations/hydrology/flow_routing.cpp` (`hydrology.flow_routing`, diagnoser) and `src/simulations/terrain/fluvial_erosion.cpp` (`terrain.fluvial_erosion`, evolver), one translation unit per process id. Evaluation module `erosion` (`configs/evaluation/erosion.json`, `runs/registry/erosion.json`, `runs/erosion/<id>/`, ledger `runs/ledger/erosion.jsonl`), because the module is judged in isolation with its own baseline (evaluation philosophy, 2026-10-04). Phase 0 entries are in the `world` ledger because the module does not exist yet.
- **Tests are part of every kernel step** (agent-workflow: new behaviour comes with its own deterministic test). Suite `Landscape` in `world_tests`, synthetic grids of a few hundred cells, exact expected numbers; `world_tests` links `uw_climate_world` and `uw_landscape`, both planet-free, and the comment on its link line says so.

## 2. Algorithms considered

Facts verified against the sources on 2026-10-04 by three read-only research agents (sources at the end). "Role" is the decision for this project.

| Candidate | Verified facts that matter here | Role | Reason |
|---|---|---|---|
| Stream power, implicit, O(n) (Braun & Willett 2013) | one sweep from outlets upstream per time step; unconditionally stable, so 1 Myr steps on a 20 km grid; D8 single receiver | **global kernel** | the cost model fits whole-planet grids (2048 x 1024 in seconds to minutes, inferred, unmeasured); every quantity the layers need falls out of it |
| Sediment deposition (Yuan et al. 2019) | erosion-deposition law with coefficient G; O(n), implicit, Gauss-Seidel sweeps whose count depends on G, not on resolution; G = 0 is FastScape | **global kernel** | without deposition every valley is bedrock to the coast and no floodplain, foreland or delta exists for fertility, placers or settlement |
| Depression handling: Priority-Flood (Barnes et al. 2014), Cordonnier et al. 2019, Fill-Spill-Merge (Barnes et al. 2021) | Priority-Flood: O(n log n), a surface where every cell drains, spill points, basin labels; Cordonnier 2019: linear-time rewiring without changing elevations; Fill-Spill-Merge: finite water volumes, lakes that need not overflow | **global kernel** (Priority-Flood with spill points in Phase A; water-balance lake levels in Phase B) | lakes and closed basins follow from terrain and climate instead of the fill stage and the basin editor |
| goSPL (Salles et al. 2020, 2023) | whole-planet icosahedral mesh at 5-10 km; erosion κ P^d Q^m S^n with discharge from precipitation; Yuan deposition; depressions fill with sediment or spill; flexure; PETSc and MPI, Windows only through WSL or Docker | design reference | the closest published analogue to this pipeline; too heavy to embed |
| Salève (Steer 2021) | analytical solution of the linear (n = 1) stream power law along the drainage tree; steady state and a transient form; converges in tens to about 130 fixed-point iterations; U, K and runoff may vary in space; no deposition; no continuity in time; hillslopes as a threshold slope; MATLAB, MIT | **regional kernel** (Phase E) | ideal where the downstream state is known and a steady relief is wanted at high resolution; wrong as the global model because it cannot deposit and a steady solve under today's uplift flattens every range that is no longer rising |
| Tzathas, Gailleton, Steer, Cordonnier 2024 | extends Salève with a non-flat starting terrain, fixed boundary heights, elapsed time, a prescribed river network, a hillslope term and multigrid acceleration (2048² in 34 s in Python); no deposition; code non-commercial only | **regional kernel** (Phase E), reimplemented from the paper | provides exactly the controls the handover from the global pass needs (pinned rivers, boundary heights) |
| fastlem (TadaTeruki) | steady state only; area exponent fixed at 0.5; no rainfall input; planar, no wrap; archived read-only since June 2025 | none | the formula is three lines and is reimplemented in Phase E; the library offers nothing beyond it |
| Runevision erosion filter (2026) | needs height and gradient; evaluable point by point; gullies follow the local slope and, per the author, can stop "halfway down a mountainside"; no water or sediment; MPL-2.0 | **render detail** (Phase F), below the resolved rivers | the look of gullies at render resolution; it cannot express drainage, so it must never move water or feed a simulation |
| Schott et al. 2023; FastFlow (Jain et al. 2024) | GPU stream power; Schott: explicit, no deposition, no lakes, 1024² in 9 s, MIT; FastFlow: implicit with deposition | none now | the CPU implicit kernel is fast enough at global resolution; revisit only if the regional pass needs the GPU |
| Cordonnier et al. 2016 | uplift map to a stream graph on a Delaunay mesh, lakes by overflow, implicit solve, 100-300 iterations; deposition neglected | none | superseded by the implicit kernel with deposition on the existing grid |
| Badlands | regional to continental triangular mesh, no sphere; Python, C and Fortran | none | neither spherical nor embeddable |
| DLA (setanarut/dla and the terrain variants) | the repository is a generic demo (fixed 854 x 480 grid, PNG out); the terrain variants grow ridges, not valleys; about 5 h for 512²; none found that conditions on an existing heightmap | none | no causal link to drainage, no way to respect the global terrain, too slow |
| Hydraulic particle and shallow-water erosion | time-step limited, metre-scale processes | none | wrong scale for 1-20 km cells; the render filter covers the look |
| Thermal (slope-limited) erosion | mass moved downslope above a critical slope | included as the hillslope term of the global kernel (linear diffusion in Phase A; a critical-slope limiter is a parameter) | cheap, and the Salève hillslope is a threshold slope anyway |
| Glacial erosion in stream-power form (Hergarten 2021) | needs an equilibrium-line altitude and precipitation; 1.7-2.7 times the fluvial cost | optional Phase G | relevant once the cryosphere has land ice; the fjord stage is its consumer |
| Orographic precipitation models (Smith & Barstad 2004; LFPM 2022) | linear theory, planar; LFPM linear-time | none | the climate module already supplies precipitation on the terrain |
| TTLEM, CAESAR-Lisflood | explicit, CFL-limited; catchment to regional scale | none | not plausible on a planet grid |
| Learned terrain synthesis (super-resolution trained on Earth DEMs) | from memory, unverified | excluded | observed Earth data must not force a procedural world (world-model rule 4) |

## 3. Design

### Fields (catalogue additions, `src/core/world/world_fields.hpp`)

| Field | Component | Type, layout, units | Written by |
|---|---|---|---|
| `radius` | planetbody | f32 scalar, m | world setup from `world.body.radius_km` |
| `flow_direction` | hydrosphere | u8 grid, 0 outlet or sea, 1..8 as `river_direction` | flow_routing, fluvial_erosion |
| `drainage_area` | hydrosphere | f32 grid, km² | same |
| `discharge` | hydrosphere | f32 grid, m³/s annual mean | same |
| `filled_elevation` | hydrosphere | f32 grid, m above datum: the Priority-Flood surface; equals the terrain outside depressions | same |
| `depression_id` | hydrosphere | i32 grid, 0 none | same |
| `erodibility` | geosphere | f32 grid, relative to the parameter | fluvial_erosion |
| `erosion_rate`, `deposition_rate` | geosphere | f32 grid, m/Myr, mean over the run | fluvial_erosion |
| `denudation`, `aggradation` | geosphere | f32 grid, m, cumulative | fluvial_erosion |
| `erosion_duration` | geosphere | f32 scalar, Myr | fluvial_erosion |
| `eroded_volume`, `land_deposited_volume`, `marine_deposited_volume` | geosphere | f32 scalars, km³ | fluvial_erosion |

Lake surfaces stay in `hydrosphere.lake_surface` and are written in Phase B. Tables: a `hydrosphere.depressions` table (id, spill cell, spill elevation, volume, area, outlet depression) when a consumer needs it (Phase B).

### Kernels (`src/simulations/landscape/`, namespace `landscape`, planet-free)

- `sphere_metric.{hpp,cpp}`: cell-centred W x W/2 lat-lon grid on a sphere of radius R. Row j spans latitudes φ_n = 90° − 180°·j/H to φ_s = φ_n − 180°/H; area a_j = R² Δλ (sin φ_n − sin φ_s) with Δλ = 2π/W; east-west spacing dx_j = R Δλ cos φ_j at the row centre; north-south spacing dy = π R / H; diagonal = hypot(mean dx of the two rows, dy). Longitude wraps; a neighbour beyond a pole does not exist (`neighbour` returns the cell itself, as `climate_erosion.cpp` does). The whole area sums to 4πR².
- `flow_graph.{hpp,cpp}`: `struct FlowGraph { std::vector<int32_t> receiver; std::vector<uint8_t> direction; std::vector<float> distance_m; std::vector<int32_t> order; std::vector<double> area_m2, discharge_m3_yr; }`. Receiver = the neighbour of steepest descent (drop over distance) among the eight; a sea neighbour counts with elevation = sea level and is an outlet; outlets have `receiver = self`, `direction = 0`. `order` is the Braun-Willett stack: outlets first, each cell after its receiver. Accumulation runs over `order` reversed: area and discharge of a cell include every upstream cell; discharge = Σ runoff · area. Routing inside depressions comes from the depression kernel (below), so every land cell has a receiver.
- `depressions.{hpp,cpp}`: Priority-Flood seeded from the sea cells (nothing else: polar rows are closed). Outputs `filled_elevation` (≥ terrain), a depression id per cell whose filled elevation exceeds the terrain by more than ε, each depression's spill cell and spill elevation, and the receiver of every flooded cell as the neighbour it was reached from in the flood order (Barnes' flow directions), so flow inside a depression leads to its spill. Phase A: lake level = spill elevation (every depression overflows). Phase B: a water balance per depression (inflow − evaporation · surface) gives a level at or below the spill, Fill-Spill-Merge semantics simplified to the depression tree.
- `stream_power.{hpp,cpp}`: one time step of length Δt on a `FlowGraph`, from outlets upstream (`order`). Detachment-limited part (n = 1): with F_i = K_i Q_i^m Δt / L_i, `z_i ← (z_i + U_i Δt + F_i z_r) / (1 + F_i)` where z_r is the receiver's already updated elevation (Braun & Willett 2013). Deposition: the Davy-Lague term of Yuan et al. 2019 with coefficient G, solved by Gauss-Seidel sweeps over the same order until the largest change is below a tolerance (G = 0 gives the detachment-limited scheme exactly). Hillslope: linear diffusion with coefficient κ, explicit with the per-face flux limiter already in `climate_erosion.cpp`. Sea cells: sediment that reaches an outlet is routed down the seabed's steepest descent and deposited with an e-folding length (`sediment_travel_km`, as the engine and `climate_erosion.cpp` do), so the budget closes: eroded = land deposited + marine deposited. Returns a step budget (volumes in km³). Erodibility, uplift and runoff are per-cell inputs; m, n, G, κ, Δt and the travel length are parameters.
- Determinism: the sweeps are serial over `order`; per-cell loops may run in parallel with fixed-order reductions; no random draws (neither process is stochastic). Results must be bit-identical across thread counts (consolidation plan §1).

### Processes

- `hydrology.flow_routing` (diagnoser): reads `base_elevation`, `sea_level`, `precipitation` (seasonal mean → mm/yr), `temperature` (snow and ice reduction as in `climate_erosion.cpp`), `radius`; writes the hydrosphere fields of the table. Runoff = `landscape.runoff.coefficient` x precipitation in Phase A (0.4; evaporation-based runoff is a Phase B refinement). In `earth-hydrology` (Phase B) it runs on Earth terrain with the Geosphere frozen.
- `terrain.fluvial_erosion` (evolver): reads the same plus `crust_class`, `crust_age`; writes `base_elevation`, `max_elevation`, the geosphere fields of the table and, for its final state, the routing fields (so a recipe needs no second routing step after it). Loop: for `duration_myr / step_myr` steps: route (depressions, graph, discharge) → erode and deposit → update. Erodibility K_i = `landscape.stream_power.erodibility` x a per-crust-class multiplier (`crust_class_multipliers`, default 1, 1.3, 1.5 for continental, transitional, oceanic) x an age factor (default off). Default erodibility: the engine's river law is `river_erosion` = 15 m/Myr per km of √(catchment km²) with uniform runoff; with discharge Q = r P A and P = 1 m/yr, r = 0.4, the same incision needs K_Q = 15/√0.4 ≈ 24 m/Myr per √(km³/yr), which is the default of `erodibility` (document the derivation beside the value). m = 0.5, n = 1, G = 1, κ = 200 km²/Myr, `step_myr` = 1, `duration_myr` = 20, `sediment_travel_km` = 300, uplift 0, isostatic rebound fraction 0 (Phase C). The first judging run checks the land-mean denudation rate against 10-100 m/Myr; outside that band the default erodibility is changed before the run is registered, not tuned by look.
- Recipe `eroded-world`: the `procedural-world` steps with `hydrology.flow_routing` and `terrain.fluvial_erosion` inserted after `climate.rainfall_and_coupling` and before `terrain.fjords`. The later legacy stages (`rivers_and_basins` with its depression fill and land noise, lakes, deltas) still run after it in Phase A and partly rework the surface; Phase B replaces them in this recipe. Recipe `earth-hydrology` (Phase B): `earth-climate` followed by `hydrology.flow_routing`, Geosphere frozen.
- Time: `landscape.erosion.duration_myr` and `step_myr` are parameters set with `--set` in Phase A; the evolver joins the time model (`--duration`, checkpoints) in Phase C.

### Evaluation (`configs/evaluation/erosion.json`)

- Run: `world_pipeline run --recipe eroded-world --resolution {width}`; registered settings at the minimal configuration (128 x 64, `--duration 50`, climate 32 columns, 2 coupling iterations); the judging run at 256 columns with the default duration.
- Checks: `grid`; `required` and `written` for the new fields; `finite`; `range` (`erosion_rate` 0..2000 m/Myr, `discharge` ≥ 0, `drainage_area` ≥ 0); three new check kinds in `scripts/evaluation/checks.py`: `budget` (`eroded_volume` = `land_deposited_volume` + `marine_deposited_volume` within 1e-6 relative), `drained` (every land cell has `flow_direction` ≠ 0 or a `depression_id`), `monotone` (`discharge` at the receiver ≥ `discharge` at the cell, following `flow_direction`); `repeat`. The land-mean denudation rate and the lake count are reported numbers, not checks.
- Validation: `hypsometry` as context. Phase B adds the Earth river-mouth comparison (§5).
- Maps: elevation (tectonics palette), `erosion_rate`, `deposition_rate` (log), `discharge` (log), `filled_elevation − base_elevation` (lakes), `depression_id`, drainage basins; review entries drafted for the user to rewrite: rivers reach the sea or a lake; the largest rivers drain the wettest highlands; valleys widen downstream and plains aggrade; no range is erased and no grid pattern appears.

### Tests (`tests/simulations/landscape/landscape_tests.cpp`, suite `Landscape`, one fixture per kernel)

- `SphereMetric`: areas sum to 4πR² (relative 1e-9); the equator spacing is 2πR/W; a neighbour beyond a pole is the cell itself; the wrap neighbour of column 0 is column W−1.
- `FlowGraph` on a 16 x 8 plane tilted towards a sea column: every land cell's receiver is lower or sea; `order` lists each receiver before its cell; the area accumulated at each outlet equals the land area of its band; discharge equals runoff times area.
- `Depressions`: a 3 x 3 bowl in that plane: the filled surface equals the spill elevation; one id; the spill is the lowest rim cell; every bowl cell's receiver chain reaches the spill; no undrained cell.
- `StreamPower`: a 32-cell chain with uniform U and K reaches the analytical steady profile (slope U/(K Q^m) per link) within 1e-3 relative; Δt = 1 Myr and 0.1 Myr give the same steady state within 1e-3; with U = 0 a hill lowers monotonically and never below its receiver; with G > 0 on a domain with one outlet, eroded = deposited + exported within 1e-9; the result is bit-identical at 1 and 8 threads.

### Consumers (Phase B)

Rivers (`river_direction`, the January and July flows in legacy units through `river_factor`) from the routing's discharge split by winter and summer runoff; lakes from the water balance; `basin_class` from whether a depression spills to the sea; mineral exhumation (`mineral-systems.md`, "reference-climate proxy until ... per-cell erosion") from `denudation`; `erosion_potential`, `deposition_potential`, `floodplain_fertility` from the rates; placers routed along `flow_direction`.

## 4. Phases and steps

Durations are commit to commit (`AGENTS.md`); the Estimate column is checked against them later. "Agent" names the project agent of `.claude/agents/` the orchestrator dispatches; "lead" is the planning session, "orch" the orchestrating session itself. Steps in one parallel group have disjoint files and run as background workers at the same time. Every step ends with the checks named in §5 for its kind.

**Phase 0 — decision and cleanup** (lead session, 2026-10-06; `world` ledger).

| Step | Work | Agent | Check | Estimate | Duration |
|---|---|---|---|---|---|
| 0.1 | Delete the FastLEM stages, `broadenfastlemterrainfromrivers`, the `terrain.fastlem` parameters, `islegacygenerationstage` and their references | worker | step check with `--desktop`, every bundle identical | 1 h | done, `3b3395a`, `world-0080`: 55 min (e1254e7 12:13 to 13:08), of which the worker ran about 6 min; the three bundles identical, both presets built |
| 0.2 | This plan; `docs/reference-guides/agent-workflow.md`; the agents `worker`, `mechanic`, `surveyor`, `reviewer`; consolidation plan §8 and §11; `AGENTS.md`; `~/.claude/CLAUDE.md`; memory | lead | none (documents) | 2 h | done, `16ddfb3`, `world-0081`: 6 min of commits (3b3395a 13:08 to 16ddfb3 13:14); the lead session wrote the documents before 3b3395a, so the commit range understates it |

**Phase A — the global kernel and the eroded recipe** (one orchestrating session, Opus 5.5 `high`; `world` ledger until A3 creates the `erosion` module, then `erosion`).

| Step | Work | Agent | Parallel group | Check | Estimate | Duration |
|---|---|---|---|---|---|---|
| A0 | Contracts first: the four kernel headers of §3 with their structs and function signatures, `src/simulations/landscape/CMakeLists.txt` (`uw_landscape`, links `uw_core`), the empty test file registered in `tests/CMakeLists.txt`, `configs/parameters/landscape.json` with every parameter of §3 and its header | orch | — | headless build | 45 min | |
| A1 | `sphere_metric` and `flow_graph` and `depressions` with their tests | worker | 1 | `ctest -R "^Landscape\.(SphereMetric\|FlowGraph\|Depressions)"`, headless build | 3 h | |
| A2 | `stream_power` (implicit SPL, Yuan deposition, hillslope diffusion, marine e-folding deposition, step budget) with its tests, on a hand-built `FlowGraph` | worker | 1 | `ctest -R "^Landscape\.StreamPower"`, headless build | 3 h | |
| A3 | Plumbing: catalogue fields and `radius`; stage ids, process rows (`hydrology.flow_routing` diagnoser, `terrain.fluvial_erosion` evolver), recipe `eroded-world`; `uw.py` module `erosion`, `runs/registry/erosion.json`, `runs/erosion/` through `storage.route`, `configs/evaluation/erosion.json`; check kinds `budget`, `drained`, `monotone` in `checks.py` with a Python unit test on a tiny synthetic bundle; `step_check.py` gains `eroded-world` | mechanic | 1 | headless build; `python -m unittest tests.scripts.test_evaluation_checks`; step check (the three existing recipes identical) | 2 h | |
| A4 | Integration: `flow_routing.cpp` and `fluvial_erosion.cpp` (planet to kernel adapters, declared access), remove `climate_erosion.{hpp,cpp}` and `uw_climate_erosion`; access measured for the two new processes with the tracking build at the minimal configuration and stored as a new report set (`process_access.py`) | worker | after group 1 | step check with `eroded-world` (first bundle, no reference) and the three existing recipes identical; `budget`, `drained`, `monotone` pass on the minimal `eroded-world` bundle (`uw.py erosion evaluate`) | 2.5 h | |
| A5 | Review of A1, A2 and A4 against §3 and their ledger drafts; fixes by the same workers | reviewer, then the workers | — | the step's own checks again | 1 h | |
| A6 | Close: the judging run `uw.py erosion run --resolution 256` (default duration; registered as the module's first run and baseline, "pending user review"); its numbers against world run 10 (land-mean denudation, hypsometry shift, lake count, largest discharges); status note here; §6; memory; the prompt for Phase B; the commit range listed for the user's adversarial review | orch | — | the run's evaluation | 1 h | |

**Phase B — hydrology and the layers on the shared routing** (`erosion` ledger).

| Step | Work | Agent | Parallel group | Check | Estimate | Duration |
|---|---|---|---|---|---|---|
| B1 | `hydrology.rivers` (diagnoser): `river_direction` and the two seasonal flows from the routing's discharge, in `eroded-world` instead of `rivers_and_basins`; the legacy stage stays in `procedural-world` | worker | 1 | step check; river counts and total flow against the A6 run | 3 h | |
| B2 | Lake levels by water balance (inflow − evaporation x surface, level ≤ spill; salt lakes and pans where endorheic and arid), writing `lake_surface` and `special_features`; `hydrology.lakes` out of `eroded-world` | worker | after B1 | test: a closed bowl under a given inflow and evaporation reaches the expected level; step check | 3 h | |
| B3 | `basin_class` from the depression tree; `hydrology.basin_policy` out of `eroded-world` | worker | 1 | step check | 1 h | |
| B4 | Consumers: mineral exhumation from `denudation`, surface potentials from the rates, placers along `flow_direction` (in `eroded-world`; the proxies stay for `procedural-world`) | worker | 1 | step check; deposit and placer totals against the A6 run | 3 h | |
| B5 | `earth-hydrology` recipe (frozen Geosphere, Earth terrain, routing) and a validation entry `river_mouths`: the 20 largest discharge outlets of the routing against a small authored table of Earth's largest river mouths (coordinates), scored by distance; it is observed data used as validation, which the rules allow | worker | 1 | a run of `earth-hydrology` at 512 columns once (the Earth references are 512 wide; one run, recorded) | 2 h | |
| B6 | Review (B1, B2, B4), judging run at 256, close, prompt for Phase C | reviewer, orch | — | | 1.5 h | |

Decided here for B2 because the user said lakes are revisited later and their algorithm will probably change (2026-10-05): B2 is that revision; a 512-column run is not needed because the depression kernel is unit-tested and the 256-column judging run shows lakes at that size (the minimal world has none).

**Phase C — climate feedback and time control** (`erosion` ledger).

| Step | Work | Agent | Check | Estimate | Duration |
|---|---|---|---|---|---|
| C1 | Recipe steps carry parameter overrides and a process may appear twice: `eroded-world` becomes coarse climate (32 columns, 2 coupling iterations as recipe fields) → routing → erosion → the full climate → hydrology, as the README's lifecycle says | worker | step check; the existing recipes identical | 3 h | |
| C2 | `terrain.fluvial_erosion` declares its duration and step and joins `--duration`, `--checkpoint-every` and `--checkpoint-at` (time model of the consolidation plan §1) | worker | step check; a checkpointed minimal run | 2 h | |
| C3 | Local isostatic rebound of eroded columns (`landscape.isostasy.rebound_fraction`, default from the engine's densities) | worker | test: a uniformly denuded plateau rebounds by the fraction; step check | 1.5 h | |
| C4 | Promotion (**user**): `eroded-world` becomes `procedural-world` and the legacy editors it replaced become optional stages or go | mechanic | step check; new world baseline | 1 h | |

**Phase D — engine coupling** (experiment; `tectonics` ledger; `ENGINE_FINDINGS.md`).

| Step | Work | Agent | Check | Estimate | Duration |
|---|---|---|---|---|---|
| D1 | The engine exports a rock uplift rate (m/Myr) per cell in its output contract (**user**: contract change) and `fluvial_erosion` reads it as U | worker | step check; tectonics checks | 2 h | |
| D2 | The engine's `routeSediment` uses the `landscape` kernel on its torus metric with uniform runoff (same law, implicit instead of explicit), then precipitation from a coarse climate at tectonic checkpoints through `erosion_multiplier` | worker | tectonics parity run at the minimal configuration, then one 256-column tectonics run against run 75; hypsometry (engine finding 5) | 4 h | |

**Phase E — regional analytical pass** (parked until the user works with the desktop again; design recorded so it can start cold).

A headless process `terrain.regional_relief` over a window (a longitude-latitude box, factor 16 as the regional map) that writes a regional bundle: upsample the global terrain (bicubic plus seeded roughness), pin the global rivers (cells of `flow_direction` with discharge above a threshold keep their elevations and are the network of Tzathas et al. 2024), fix the window's boundary heights, and solve the analytical steady state of Steer 2021 with the Tzathas hillslope term and multigrid fixed point; discharge from the global routing scaled to the fine cells. Seams between neighbouring windows are the open problem: pinned rivers plus a halo of one global cell is the first thing to try, checked by a test that generates two overlapping windows and compares the overlap. The GUI shows the bundle when it is a viewer (consolidation plan Phase 5).

**Phase F — render detail** (parked, presentation only). A point-evaluated gully filter in `src/app/rendering` after the lineage of the Runevision filter (MPL-2.0 for its code; a reimplementation from the description needs no licence), applied below the resolved channels and faded to zero near rivers, lakes and plains; it never writes a field.

**Phase G — glacial erosion** (optional). Hergarten 2021's ice-flux term in `stream_power` when the cryosphere has land ice and an equilibrium-line altitude; the fjord stage reads it.

## 5. Rules of the task

- Checks per kind of step: a kernel step ends with its suite (`ctest -R "^Landscape\.<Fixture>"`) and the headless build; a plumbing step with the headless build, the step check (existing recipes identical) and any Python unit test it adds; an integration step with the step check including `eroded-world` and `uw.py erosion evaluate` on the minimal bundle; a phase close with one 256-column judging run and nothing larger. No full ctest, no Python suite, no desktop launch unless a step changes what they cover; the desktop preset is built only when `src/app` changes (Phase F).
- Minimal configuration for every "did anything change" question; its world has little relief and no lakes, so a kernel's correctness is shown by its tests and the profile checks, not by that world's look; the judging run is where the user looks.
- Erodibility and the other defaults are changed by a recorded decision with a number (the denudation band of §3), never by eye; every value is a parameter in `configs/parameters/landscape.json` and in the bundle manifest.
- Parallel workers: disjoint files (group 1 of Phase A: `landscape/sphere_metric*`, `flow_graph*`, `depressions*` and their tests; `landscape/stream_power*` and its test; everything outside `src/simulations/landscape` and the test file). The orchestrator writes the shared headers first (A0) and owns `CMakeLists.txt` files and `tests/CMakeLists.txt` after that.
- A worker reports; the orchestrator reads the diff, runs the step check, writes the ledger entry and commits (one commit per step, duration from the previous commit). The reviewer sees the diff, the plan step and the draft entry only.
- The adversarial review (GPT-6 Astra, by the user) gets the phase's commit range and this document; its findings become steps or §6 items.
- Nothing is deleted or migrated under `refs/` or `runs/` by these sessions except the run directories the run tooling creates.

## 6. Questions answered in advance, risks, open items

Answered so that no session asks:

- Where does code go: kernels in `src/simulations/landscape` (planet-free, tested), processes in the module whose fields they write, one translation unit per process id. `world_tests` links `uw_landscape`; the "only project library" comment on its link line becomes "planet-free libraries only".
- `GenerationStageId` gains `flow_routing` and `fluvial_erosion`; the workbench's stage-to-view mapping puts them under relief; the GUI is compile-only (user, 2026-10-05).
- Declared access for new processes is measured (tracking build, minimal `eroded-world` run, report set stored beside the Phase 4 sets, `process_access.py`), as the consolidation plan requires; `step_check.py` keeps the three existing recipes and adds `eroded-world`.
- The second routing after erosion is the evolver writing its final routing fields; a recipe does not list a process twice until C1.
- Discharge in m³/s in the bundle (annual mean), m³/yr inside the kernel; the legacy river flows keep their units through `river_factor` until the GUI draws from bundles.
- Poles: rows 0 and H−1 have no neighbour beyond the pole; D8 directions are biased where cells are elongated (above about 70°); accepted in Phase A, noted for the review.
- If the minimal `eroded-world` run shows a budget residual above tolerance, the cause is in the marine routing or the limiter, not in the tolerance: the tolerance is not raised.
- Procedural worlds never carry template mountains since the FastLEM removal, so `base_elevation` is the land surface; if imports with mountain templates ever enter `eroded-world`, the templates are folded into the base elevation first (a step of its own).
- The `world` ledger holds Phase 0; Phase A's entries move to `erosion` once A3 exists (the orchestrator creates `runs/ledger/erosion.jsonl` through the first `--add`; if `uw.py` needs the module in its `COMMANDS`, A3 adds it).

Found and left (not this plan's work):

- `generateglobalclimate` in `src/pipeline/climate_generation.cpp` has no caller (found in step 0.1 while removing the post-river call from it); `climate_stages.cpp` runs the climate stages. Dead code of about 700 lines, to be deleted by a cleanup step that checks the three recipes stay identical.

Risks:

- The engine's relief at 128 x 64 after 50 Myr may be too low for erosion to do anything visible; the checks still hold (budget, drained, monotone) and the 256-column run is the look.
- A land-mean denudation outside 10-100 m/Myr at the default erodibility: change the default before registering, record the number.
- The legacy stages after erosion in Phase A (`rivers_and_basins` fills depressions and adds land noise) partly undo lakes and plains; Phase B removes them from `eroded-world`. The A6 judging run is read with that in mind.
- Thread determinism of the Gauss-Seidel sweeps is by construction (serial); a parallel per-cell loop with a reduction (volumes) must accumulate in a fixed order.

Open (**user**):

- Promotion of `eroded-world` to the default recipe (C4) after the A6 and B6 verdicts.
- The rock uplift rate in the tectonic output contract (D1).
- The authored river-mouth table of B5 (observed data as validation; the user may prefer another list or none).
- Whether Phase E waits for the desktop or gets a headless image export of its own first.

## 7. Prompt for the Phase A session

Start a fresh session on Opus 5.5 at `high` (`--model opus`) and give it:

> Run Phase A of `docs/architecture/erosion-hydrology-plan.md` (read §1, §3, §4 Phase A, §5 and §6 first, then `docs/reference-guides/agent-workflow.md`, `AGENTS.md`, and the ledger entries of Phase 0 with `uv run scripts/uw.py world ledger --topic "erosion plan"`). Do A0 yourself, then launch A1, A2 and A3 as background agents (`worker`, `worker`, `mechanic`) in one message with briefs as the guide describes, and keep A4's adapter skeleton ready while they run. One commit and one ledger entry per step with its measured duration; the step check at the minimal configuration is the only run until A6; review A1, A2 and A4 with `reviewer` before their entries are final. Stop and record in §6 any decision the plan does not settle; do not stop for ones it does. End with A6: the 256-column judging run registered as the `erosion` baseline pending the user's review, the status note in §4, §6 updated, memory updated, and the prompt for Phase B.

## Sources

Steer 2021 (Salève): https://esurf.copernicus.org/articles/9/1239/2021/ · Tzathas et al. 2024: https://www-sop.inria.fr/reves/Basilic/2024/TGSC24/Analytical_Terrains_EG.pdf · fastlem: https://github.com/TadaTeruki/fastlem · Braun & Willett 2013: https://doi.org/10.1016/j.geomorph.2012.10.008 · Yuan et al. 2019: https://hal.science/hal-02136641/document · Cordonnier, Bovy & Braun 2019: https://esurf.copernicus.org/articles/7/549/2019/ · Priority-Flood: https://arxiv.org/abs/1511.04463 · Fill-Spill-Merge: https://esurf.copernicus.org/articles/9/105/2021/ · goSPL: https://github.com/Geodels/gospl and https://hal.science/hal-04880511/document · fastscapelib: https://github.com/fastscape-lem/fastscapelib · Badlands: https://github.com/badlands-model/badlands · Runevision filter: https://blog.runevision.com/2026/03/fast-and-gorgeous-erosion-filter.html · setanarut/dla: https://github.com/setanarut/dla · ISIM-DLA: https://github.com/Sevenudeloin/ISIM-DLA-rendering · Schott et al. 2023: https://hal.science/hal-04049125/document · FastFlow 2024: https://www-sop.inria.fr/reves/Basilic/2024/JKGFC24/FastFlowPG2024_Author_Version.pdf · Cordonnier et al. 2016: https://www.cs.purdue.edu/cgvlab/www/resources/papers/Cordonnier-Computer_Graphics_Forum-2016-Large_Scale_Terrain_Generation_from_Tectonic_Uplift_and_Fluvial_.pdf · Hergarten 2021 (glacial): https://esurf.copernicus.org/articles/9/937/2021/ · Hergarten 2020 (transport-limited): https://esurf.copernicus.org/articles/8/841/2020/ · LFPM: https://gmd.copernicus.org/articles/15/2063/2022/
