# Master plan: everything left to build, in numbered phases

Written 2026-10-20 (Opus 5.5), after T8 (the cubed-sphere world) closed with the user's verdict: the cube world is adopted. This is the only running plan. It collects every item not yet done or found broken in the earlier plans, which stay as records of what was decided and learned:

- `terrain-roadmap.md` (T0 to T8: sandbox, port, multiscale, engine, the cube);
- `multiscale-plan.md` and `erosion-hydrology-plan.md` (their open items had moved to the roadmap; about fifteen had not, and they are here now);
- `consolidation-plan.md`, `world-model.md`, `mineral-systems.md`, `noise-functions.md`, `simulation-modules.md`, `tectonics/*`;
- `docs/tasks/` (`main/PLAN.md`, `main/world-model-migration.md`, `climate-physics/FUTURE_WORK.md`, `climate-physics/ASTRA_HANDOFF.md`, `tectonics/ENGINE_FINDINGS.md`);
- the README's intended lifecycle, and the resource-chain notes.

Every step names its origin, so its detail is one lookup away. Items found stale are listed in §6 with the reason. Phases are numbered 1 to 17 and steps `phase.step` (3.2 is phase 3, step 2). The old letters (T4.12, S3.17, W2-E, E5, ...) appear only in the Origin column.

## 1. Overview

**Attention**: *now* means the next sessions; *next* starts when a *now* phase on the same lane frees it; *later* waits its turn; *waits* is blocked on the user.

| # | Phase | Lanes | Attention | Waits on | Estimate (agent h) |
|---|---|---|---|---|---|
| 1 | Wrap-up after the cube | tools, world, docs | now | | 4 |
| 2 | Repository hygiene and storage automation | tools | now | | 6 |
| 3 | Speed: engine, climate coupling, cube windows | engine, climate, windows | now | | 8 |
| 4 | Sandbox tools on the cube: tributaries, nested windows | windows, sandbox | next | | 4 |
| 5 | Rivers without the grid's bias (routing item 3) | landscape, sandbox | now | | 8 |
| 6 | Finish the look in the sandbox on base-v10 | sandbox | next | the user's look | 6 + looks |
| 7 | A clean regional base on the cube | windows | next | 7.3: reach verdict | 8 |
| 8 | Port the approved look into the pipeline (the real goal) | terrain, windows, hydrology | next | 6.5 | 25 |
| 9 | Tectonics realism | engine | now (9.1 to 9.3) | | 30 |
| 10 | Climate quality | climate | later (10.1 next) | | 25 |
| 11 | Hydrology and coasts | hydrology, terrain | later | | 10 |
| 12 | The sea floor | sandbox, terrain, engine | waits | the user's research | |
| 13 | Large windows: tiles, GPU, 1 m | windows | later | 7 | 20 |
| 14 | World model and desktop | world, desktop | later (14.2 next) | | 20 |
| 15 | Life, resources and society | life | later | 8, 10 | 40 |
| 16 | Later phenomena (catalogue) | any | later | | |
| 17 | Reviews and the adversarial review | | at the end | | |

### Tracks (the user, 2026-10-21)

The phases are worked in **tracks, one at a time**, each streamlined like T8: parallel workers within the track, the next track when its phases close. An issue seen during a step is triaged (`AGENTS.md`, "Tracks and triage"): fixed in the step if it belongs to it, else one row in the phase that owns it, without interrupting the track.

| Order | Track | Phases and steps |
|---|---|---|
| A (now) | Repository, rules and protocols | 1.5; Phase 2 (2.2 to 2.8); the open items of §5 that tools can close |
| B | Terrain refining (above and below sea) | 3.3; Phases 4 to 9, 11, 13; 6.1b, 6.1c; Phase 12 when the user's research arrives |
| C | Climate rework | 3.2 (the fix); Phase 10 except the ocean steps |
| D | Ocean currents | 10.3 and the ocean items of 10.6 (wet mask, sea ice) |
| E | Mineralogy | 15.3 |
| F | Life, society, world model, desktop | 15.1, 15.2, 15.4; Phase 14 |
| last | Reviews | Phase 17 |

The *Attention* column above predates the tracks; the track order decides.

**Why this order.**
- Phase 1 closes what T8 left open, including the one real cube regression: the hydrology of the lat-lon export.
- Phase 2 is the user's request of 2026-10-18, never started.
- Phase 3 matters because every run waits on it. The full recipe spends most of its 37 to 43 min in a climate coupling that never converges.
- Phase 5 removes the last grid artefact from the rivers.
- Phases 6 to 8 are the goal the sandbox was built for: the approved look, in the pipeline.
- Phase 9 has its own lane, so the engine work runs beside all of them.
- Climate quality, hydrology calibration, the desktop and life follow once the terrain they stand on is final.

## 2. How a session runs: waves of parallel workers

The plan is organised by **lanes** (file ownership), so that independent steps run at once. A session runs one **wave**, not one phase.

**Lanes** (a worker touches only its lane's files; a change outside them goes through the orchestrator on `main`, announced to the running workers):

| Lane | Files |
|---|---|
| engine | `src/simulations/geology/plate_tectonics/**`, `geology/world_tectonics.cpp`, `configs/tectonics-baseline.json`, `src/cli/tectonic_options.*`, engine tests, `docs/architecture/tectonics/*` |
| climate | `src/simulations/climate/**`, `src/core/climate_*`, `src/pipeline/climate_stages.cpp`, `configs/parameters/climate.json`, climate tests |
| landscape | `src/simulations/landscape/*` except `regional_*`, `refinement_*`, `upsampling*` (routing, stream power, the metrics) |
| windows | `src/simulations/landscape/regional_*`, `refinement_*`, `upsampling*`, `src/cli/regional_command.cpp`, `src/io/regional_bundle.*` |
| terrain | `src/simulations/terrain/**`, `src/pipeline/world_stages.cpp`, `src/core/procedural_noise.*`, `configs/parameters/terrain.json` |
| hydrology | `src/simulations/hydrology/**`, `configs/parameters/hydrology.json`, `configs/river-mouths.json` |
| sandbox | `scripts/sandbox/**`, `src/cli/terrain_sandbox_main.cpp`, `configs/sandbox/**` |
| world | `src/world/**`, `src/core/world/**`, `src/core/mesh/**`, `src/io/run_bundle.*`, `src/pipeline/world_processes.*` |
| desktop | `src/app/**`, `src/simulations/regional/**`, `.uww` paths in `src/io` and `src/cli` |
| life | `src/simulations/resources/**`, `geology/mineral_systems/**`, `src/simulations/society/**`, `src/pipeline/social_generation.cpp`, a future `src/simulations/ecology/` |
| tools | `scripts/uwtools/**`, `scripts/evaluation/**`, `scripts/benchmarks/**`, `scripts/*.py`, `configs/evaluation/**`, `.claude/settings.json` hooks |
| docs | `docs/**`, `README.md`, `AGENTS.md` (orchestrator only) |

**A wave.**
1. The orchestrator reads this file's overview and the phases it will draw from, then runs `storage.py status`. Each C++ worktree costs about 1.7 GB on C; skip a C++ worker if C is below 40 GB.
2. It fills up to **four C++ worker slots**, each from a different lane. Each slot takes the lowest-numbered open step on that lane whose dependencies are met. Workers run in `isolation: worktree` and write runs to the main tree's `out/scratch/<phase>/` by absolute path.
3. Python-only steps (sandbox scripts, tools) and surveyors do not count against the slots.
4. **Timing steps (3.1, 3.2, 3.3, any step whose acceptance is a time) run alone on an idle machine**, at the end of the wave. Their profiling can run in parallel.
5. A step's sequence: brief, work, the step's checks, the orchestrator reads the diff, then a `reviewer` for behaviour, contract or kernel changes. Then fixes by the same worker (SendMessage), the merge into `main`, the ledger entry, and one commit (`docs/reference-guides/agent-workflow.md`).
6. Parity references: `step_check.py <label> --reference t8-cube4` (cube) and `--grid latlon --reference t8-ll4`. Each label is used once, and a new label becomes the next reference.
7. Merged worktrees are cleaned up at the wave's end by the orchestrator:
   - `worktrees.py preserve NAME`, then `retire NAME --yes`, which refuses unless all the work is in main;
   - `worktrees.py drop-branches --retired`, which bundles, verifies and checks heads before deleting;
   - scratch no step needs goes to `trash.py move` (AGENTS.md "Autonomy").
8. The wave closes with:
   - the status rows below (with durations);
   - one list of looks for the user, with page paths;
   - what the session learned, written where the next agent reads it;
   - the next wave's prompt, printed in the chat.

**Defaults.**
- Seed 1, the cube grid, and the minimal configuration for parity.
- One judging run at 256 per behaviour phase, where the change acts.
- No 512 runs unless asked; tests only as `AGENTS.md` says.
- The user's look decides anything visual; numbers decide where the user defers.

**Decisions.** Decisions that would go to the user are taken with the recommended option and recorded in the step's row (`AGENTS.md`, "Autonomy"). The exceptions are those listed under *waits* and in §5.

## 3. Phases

Columns: **Agent** (`worker`, `mechanic`, `surveyor`, `reviewer`, `design` = an Opus 5.5 design agent, `orch` = the orchestrator), **Lane**, **Depends**, **Origin** (the old id and document; RM = `terrain-roadmap.md`, EHP = `erosion-hydrology-plan.md`, MSP = `multiscale-plan.md`, CP = `consolidation-plan.md`, WM = `world-model.md`, EF = `ENGINE_FINDINGS.md`, FW = `FUTURE_WORK.md`, AH = `ASTRA_HANDOFF.md`).

### Phase 1: Wrap-up after the cube

Goal: T8's loose ends closed, its one regression fixed, and the documents describing the code as it is.

| Step | What | Agent | Lane | Depends | Origin |
|---|---|---|---|---|---|
| 1.1 | Record the user's verdict (2026-10-20): base-v10 adopted, base-v9 a preset. The borderline land seam rank of 0.96 is accepted with the look. Status rows in RM §3 T8 closed. **Done in this plan's commit.** | orch | docs | | RM "For the user" |
| 1.2 | Register the cube tectonics baseline: the E6 run (`out/scratch/sphere/E6/`, baseline v6, N 64, 400 Myr) as a numbered run in `runs/registry/tectonics.json` (the registry names run 78, torus). Update `tectonics/pipeline-tooling.md`, which still says run 75, W x W/2, recipe v4 and bundle v1. | worker | tools | | RM:352, CP-1, PT-1 |
| 1.3 | **Hydrology of the lat-lon export.** The export reroutes on an averaged surface and raises narrow outlets 0.5 to 16 m above the cube's open lakes. Effects: endorheic land 0.277 on the export against 0.080 native; internal-drainage IoU 0.277 against lat-lon 0.315; river mouths matched 6 of 20 against 8. Carry the native outlets into the export's routing, re-score `uw.py hydrology evaluate`, and refresh the acceptance page's hydrology rows, which still show the values from before the fix. **Wave 1:** these numbers predated W3.3c's lake rule; at HEAD the export is 0.065 against native 0.080 (done-when met), mouths 8 of 20 on every grid, and two thirds of the IoU gap to lat-lon is the cube world's own (11.6). No export change; page `out/evaluation/t8-acceptance-2`. | worker, reviewer | world | | W3.3b, W3.3c |
| 1.4 | Documentation drift (recipe edits). List: the README (torus → cube; `src/world` names that no longer exist; `wip/`; `definitions/`; the `climate-benchmark` command); `simulation-modules.md` (`wip/planet_bundle`; target `uw_tectonics`); `world-model.md` (whole metres, "plumes planned"); `mineral-systems.md` (planar cell areas); `consolidation-plan.md` (`src/app/regional` → `src/simulations/regional`); `tectonic-output-contract.md` (byte 0xB2 at line 93 is not UTF-8; the pre-cube strain figure); `ENGINE_FINDINGS.md` (`tectonic_pipeline.exe`, `check-layout.py`); `FUTURE_WORK.md` (save version 19, archived paths). | mechanic | docs | | surveys of 2026-10-20 |
| 1.5 | Small code faults found and left. `planet::wrapx` adds W - 1, so x = W maps to 1. `planet::latitude` is node-centred (only `src/app` calls it). The sandbox's `make_domains` gives a negative width across the date line. Dead code: `planet::smoothnom`, `prepare`/`finishcustomworldgeneration`. Fix the first one with `step_check.py`; it may change worlds, so compare totals. | worker | world, sandbox | | RM:356, CP-14, CP-16 |
| 1.6 | The scratch files the delete guard refused (RM "For the user": `out/scratch/sphere/W1b/{head,mid,tail}.txt`, `W2-W2f/probe*`, `logs/W3.6-run.txt`, `step-check/ts-base-*` except `ts-base2`) go to the trash (`trash.py move`). | orch | tools | | RM:346-350 |

Done when: the registry names a cube baseline; the export's endorheic share is within 0.02 of native; the documents name only existing paths.

### Phase 2: Repository hygiene and storage automation

Goal: the user's requirements of 2026-10-18: "SCRIPTED, AUTOMATED. No decision up to the agent if it can be avoidable"; only the essential files on C.

| Step | What | Agent | Lane | Depends | Origin |
|---|---|---|---|---|---|
| 2.1 | **Structure enforcer.** A script beside `check_layout.py` and `repository_scope.py`, run by a hook and by `step_check.py`. It refuses folders the layout does not name, per-run HTML outside `scripts/evaluation/` rendering, and files in unnamed places. | worker | tools | | RM §"Repository hygiene" |
| 2.2 | **One naming scheme for experiments.** A design first: the registry's numeric ids plus a readable tag, in one place. Today's `vN`, `tNN`, `base-vN`, `UV`, `W` live in three places (`runs/`, `out/evaluation`, `out/scratch`). Then the tools that create run directories use it. Old labels are not renamed. | design, worker | tools | | same |
| 2.3 | **Archive without regenerable files.** A `scripts/uwtools` command with a dry run and a receipt. It keeps `.npy`, manifests and metrics JSON, and drops PNGs and pages the scripts re-render. | worker | tools | 2.2 | same |
| 2.4 | **Automatic moves.** A new run of a topic sends the previous one to D or the archive by script. Images no check reads are written to D directly. | worker | tools | 2.2 | same |
| 2.5 | `out/scratch.storage-migration-backup` (26 GB on C, verified identical to D): `storage.py finish out/scratch`, which re-verifies every file and marks the receipt complete. Agents may run it since 2026-10-20. | orch | tools | | same |
| 2.7 | **Deletion tools on read-only files.** `worktrees.py retire` and `trash.py move` stop with "access denied" on read-only files (a `.venv` `.pyd`, git pack `.idx` under `out/build/headless/_deps`), and neither resumes a half-done job: `retire` had already removed the worktree's `.git`, so it then refuses ("HEAD moved since the archive"). Clear the read-only attribute before unlinking, and let both resume a partial operation. A reviewer reads the change (AGENTS.md: guards and deleters are changed by reading, never tested with deletes). Then finish the two leftovers of wave 1: `.claude/worktrees/agent-af9bab603c2d9335e` (390 MB left; a full copy is trash entry `20261020T215954Z-0`) and `agent-a8aabce0917bb085c` (1.3 GB); both were preserved to `worktree-archive` and their work is in main (efbeb7f, 1c34a01). | worker, reviewer | tools | | wave 1 |
| 2.8 | **Communication protocols** (the user's "communication protocols", taken as the agents' and the user's exchanges): the brief, the worker's report, the reviewer's verdict, the questions to the user (all in one message; no prompt until answered) and the wave handoff, each as a fixed field list in `docs/reference-guides/agent-workflow.md` and the agent definitions (`.claude/agents/*.md`); a script checks what is machine-readable (report and ledger fields, plan status rows with commit and duration). | design, worker | tools, docs | | the user, 2026-10-21 |
| 2.9 | `uwtools/references.py` resolves `refs/` two ways (`references.ROOT` and the default reference path), which disagree in a worktree: `prepare_reduced_earth_benchmark.py --help` and its full-CLI test fail there with a `relative_to` ValueError. One root for both. | mechanic | tools | | 2.6 (triage) |
| 2.6 | Script inventory. Scripts used only by tests: `climate_convergence.py`, `terrain_wind*.py`, `climate_reference_experiment_prepare.py`. The alias modules `map_rendering` and `cell_grid` in `scripts/benchmarks` and `scripts/refs`. Paths never exercised: `storage.route` into a non-empty directory, reference preparation end to end, an audited run with a violation. Each script is kept with a test or deleted. | surveyor, then mechanic | tools | | CP-4, CP-9, CP-10, PLAN.md P2 |

Decided in 2.1: `scripts/refs/physical_map_index.html` (a renderer's template) stays allowlisted; the sandbox's compare pages (HTML written by `terrain_sandbox.py` under `out/scratch/sandbox/`) move to `out/evaluation` in 2.4; the pipeline's default `out/climate-diagnostics`, `out/images` and `out/inputs` (derived inputs, `scripts/evaluation/inputs.py`) are named places.

Done when: a misplaced file fails a hook; run directories get ids from one registry; an archive command has removed a real run's PNGs with a receipt.

### Phase 3: Speed

Goal: faster runs for every later phase. The full procedural recipe takes 2208 s at 512 and 2572 s at 256; the climate coupling dominates it and runs all 24 iterations without converging on both grids.

| Step | What | Agent | Lane | Depends | Origin |
|---|---|---|---|---|---|
| 3.1 | **Engine transport substeps.** 761 s at N 64 for 400 Myr against the 715 s gate. Two substeps run in 177 of 400 updates with v6's faster plates (223 in E6's v6 profile, `probe-v6-n64-400.json`). Profile per phase, then tune the substep limit or `transport_subdivision`. Acceptance: at most 715 s on an idle machine, with an identical bundle or a measured, reviewed change. | worker, reviewer | engine | | RM:223, RM:353, EF E12b |
| 3.2 | **Climate coupling convergence.** First the residual history per iteration (both grids, 256, seed 1) and the criterion; then relaxation or a stagnation stop. Acceptance: it converges, or stops on a stated criterion, in fewer iterations; `global-balanced-v1` does not get worse (cube 42.961). This includes the climate-run timing never measured after the parameter move. **Wave 1 diagnosis (climate-0004):** the inner gate (periodic ocean, hydrology land drift 0.14 K) is never accepted, so the test cannot pass; heating converges slowly (factor 0.976 per iteration, relaxation at its 0.1 floor); the procedural world has a period-2 tropical-ocean limit cycle in wind and rain; a stop at 11 halves the coupling but costs +0.13 % (SST). Idle timing, earth-climate at 256: cube 1112 s (coupling 1065 s), lat-lon 1039 s (996 s). **Decided for the fix:** (a) the gate requires the inner solves accepted and finite, with ocean periodicity and hydrology drift reported, not required; (b) a stagnation stop when the best heating residual improves by less than 5 % over 4 iterations; (c) a closing ocean spin-up to its tolerance after the stop, with `climate.oceancurrents.periodic_spinup_years` (24) and `periodic_tolerance_k` (0.05) as parameters; (d) the limit cycle (relaxing the ocean skin SST the winds read) is part of 3.2. | worker, reviewer | climate | | RM:353, AH A4, FW F11, CP-27 |
| 3.3 | **Cube window speed.** Cube windows take 1.3 to 2.3 times the lat-lon window's time; the budget is 2.5, the aim 1.3. Profile the level subsets, the margin rings and the neighbour-table indirection. Windows must stay bit-identical. | worker, reviewer | windows | | RM:333, RM:353 |

Done when: all three meet their targets in timing runs made alone. The full recipe's time at 256 is recorded beside 2572 s.

### Phase 4: Sandbox tools on the cube

Goal: the two sandbox tools that the cube left behind work on cube bases.

| Step | What | Agent | Lane | Depends | Origin |
|---|---|---|---|---|---|
| 4.1 | Tributaries on cube windows: the regional tributaries on level subsets, with overlaps bit-identical across a face edge. The default stays off (`tributary_enabled` 0; the user is unsure); a page with tributaries on is made for the user's look. | worker, reviewer | windows | 3.3 | RM:353, MSP M18 |
| 4.2 | `nested_window.py` (the 100 m test tool of T1.5) on cube bases. Nesting as a window mechanism stays off: the user rejected it on sight on 2026-10-17. | worker | sandbox | | RM:353, T1.5 |

### Phase 5: Rivers without the grid's bias (routing item 3)

Goal: rivers whose headings do not prefer the grid's axes. The measure is `scripts/evaluation/d8_axis.py` (mean cos 4φ: cube 0.33, lat-lon 0.43 to 0.82), together with `straight_rivers.py` A8/A4 and the axis-locked share. The design is RM §"Item 3, stage 1".

| Step | What | Agent | Lane | Depends | Origin |
|---|---|---|---|---|---|
| 5.1 | Measures on base-v10 before any code: `d8_axis.py`, A8/A4 against the isotropic calibration (A8 0.0285), the axis-locked share (95 % on base-v2), and river cells in depressions and flats. | worker | sandbox | | RM item 3 step 0 |
| 5.1b | **Measured (erosion-0107):** sandbox rivers on base-v10 are already near isotropic: A8 (lattice) 0.016 on the global and on G (isotropic 0.0285), runs of 8+ cells 6.0 / 8.9 %; the global's off-fill rivers keep A8 0.096; the pipeline's N 64 routing has d8_axis 0.303. **Decided:** 5.2 to 5.4 are judged on the pipeline's d8_axis and the global off-fill A8; the sandbox A8 (lattice) must stay at most 0.057. | | | | |
| 5.2 | **Path-corrected single receiver (D8-LTD, Orlandini et al. 2003)** in `flow_graph`: `RoutingOptions::pathcorrected` (off) and `pathdamping`, exclusive with `randomised`, on the cube's table metric. Off gives identity: `step_check.py` on both grids. | worker, reviewer | landscape | 5.1 | RM item 3 |
| 5.3 | **Euclidean fills.** `sloped_fill` (Python), the exe's `filloffsets` and `limitslopes` measure Euclidean distance instead of 8-neighbour paths, which leave octagonal cones and up to 8.2 % distance error. | worker | sandbox | 5.1 | RM item 3 |
| 5.4 | Runs `path`, `path-fill`, `path-euclid` on base-v10 (global and G) with comparison pages and the measures; the user's look; the default decided. | orch, user | sandbox | 5.2, 5.3 | RM item 3 step 2 |
| 5.5 | Stage 2: child rivers with inflow at a window's edge. It needs single-thread parents, a per-cell discharge source in `accumulateflow` and a sediment source in stream power, `accumulateice` ignoring the inflow (which retires the ELA workaround), and an `--out-directions` raster. | worker, reviewer | landscape | 5.4 | RM item 3 "For stage 2" |

### Phase 6: Finish the look in the sandbox on base-v10

Goal: the target look, frozen as the reference for Phase 8.

| Step | What | Agent | Lane | Depends | Origin |
|---|---|---|---|---|---|
| 6.1 | The details the user saw on t24 / base-v10. Ask for them with the next look list; each becomes a row here. | orch, user | sandbox | | the user, 2026-10-20 |
| 6.1a | The user (2026-10-21): "even though the new baseline looks better, the global render shows square grid-artifacts" in the box (450, 300)-(600, 450) of `t24-base-v10/global-relief.png` (2048 x 1024), and in other areas. First look: rounded-square pockets and stair-stepped shelf edges in the shallow sea west of the peninsula, at about 1 to 5 base cells; land clean at full size. **Diagnosed** (`out/scratch/wave1/6.1a/summary.txt`): the squares are the sandbox's sea floor, not the renderer, land or face edges. M1: base-v10 has one-cell sea-floor steps (-10 to -100 m beside -2 to -5 km); the Catmull-Rom sampling (`cube_domain.py:472-484`) keeps them one base cell wide and `shelf_profile` (`terrain_sandbox.py:1729-1747`) puts a scarp on that contour, so pockets become rounded squares and shelf edges staircases; nothing later breaks it up off oceanic crust. M2: the oceanic blend weight is the per-base-cell `crust_class == 1` flag (line 1671), so the switch between the crust-age floor and the shelf follows the N 64 staircase. M3: the straight N-S cliff is four stratovolcano cells in one column added as whole-cell heights (`world_fields`, lines 908-912). Extent: 148 of 1300 sea tiles (11 %), mostly off oceanic crust, not tied to face edges. | worker | sandbox | | the user, 2026-10-21 |
| 6.1b | **Fix the grid imprint at sea, without redesigning the sea floor** (Phase 12 waits for the user's research; the user agreed on 2026-10-21, since the source is the sea floor): the oceanic blend weight from a smoothed, distance-based field instead of the per-cell flag; the depth fed to `shelf_profile` sampled with a sub-cell warp so its contour leaves the cell boundaries; volcanoes drawn as cones on the fine grid. Each behind a key, compared on a page (global and the worst boxes of 6.1a) for the user's look; the default follows the look. | worker | sandbox | 6.1a | 6.1a |
| 6.1c | Upstream of 6.1b: where base-v10's flat -2.1 km sea cells and its 4-5 km deep cells of crust class 3 come from (pipeline terrain stages or the engine). Survey only. | surveyor | terrain, engine | 6.1a | 6.1a |
| 6.2 | Artefacts left from T1.6, each re-checked on base-v10 first: jagged hard/soft band scarps in the fold-belt foreland; G-coast's micro-texture blotches; the date-line seam (likely gone with the cube); tropical highland glaciation on B; oversized shields (thin engine sediment); crust under 20 km beneath G's land. | worker | sandbox | | RM T1.6, T1 notes |
| 6.3 | Glacial isostasy (optional): the ice load depresses the crust and it rebounds, giving flooded lowlands and raised beaches. | worker | sandbox | 6.2 | RM T1.4 |
| 6.4 | Lakes on base-v10: check the share and sizes against the v17 complaint ("too many and too large") and Earth. | worker | sandbox | 5.4 | RM T1.3, v17 notes |
| 6.5 | **Close.** The user's verdict on the full set (world, B, F, G, O, the k 256 windows) as the target look. `configs/sandbox/terrain.json` is frozen at that commit. | user | sandbox | 6.1 to 6.4 | RM T1.7 |

### Phase 7: A clean regional base on the cube

Goal: the pipeline's regional levels free of grid artefacts before the look is ported onto them. Most of these items were found on lat-lon windows; 7.1 decides which survive.

| Step | What | Agent | Lane | Depends | Origin |
|---|---|---|---|---|---|
| 7.1 | Re-check survey on cube windows of A, B and C:<br>- the level ≥ 2 coast lattice (9.8 km);<br>- A's N–S and C's E–W lattice peaks;<br>- B's worm-like trough texture;<br>- the parent-period ripple;<br>- beds counted twice (S3.2);<br>- routing to the pit bottom;<br>- the lakes-step items (uncapped wet pits, the straight sea–lake contact, the gate over filled pockets, the square pond, the floored sill);<br>- the 51 m step;<br>- the legend cut;<br>- `relief_amplitude` off `localrelief`'s scale;<br>- unaligned `--metres-per-pixel` overlaps;<br>- the share rule's pockets.<br>Output: which still reproduce. | surveyor, worker | windows | | MSP M1, M14, M15, M20, M21; RM T3 open items |
| 7.2 | The coast rule at level ≥ 2 that prints the lattice under the coast warp (the `ShorelineWarp` part is done). | worker, reviewer | windows | 7.1 | RM T3.1 (S3.17) |
| 7.3 | `landscape.refinement.prolong_basis` (0 Catmull-Rom, 1 cubic B-spline) for `prolong` and the level-1 start, with the `--sea` variant. | worker | windows | 7.2 | RM T3.2 (S3.18) |
| 7.4 | Connected seamless solve: descent labels plus a minimax-Dijkstra fill. **Waits:** the user's verdict on reach 0 against 1 at A k 64 (page made in this step). | worker, reviewer, user | windows | 7.3 | RM T3.3 (S4.0) |
| 7.5 | The remaining reproduced items of 7.1, one row each. | worker | windows | 7.1 | |

### Phase 8: Port the approved look into the pipeline (the real goal)

Goal: the pipeline draws what the sandbox draws, inside its contract: the global grid authoritative, every level's block means exact, windows seamless. The design is RM §3 T2 plus §"T1 notes".

| Step | What | Agent | Lane | Depends | Origin |
|---|---|---|---|---|---|
| 8.1 | **Design.** Each sandbox layer's home:<br>- **global:** erodibility from the tectonic fields, the uplands lift, the province and stability fields;<br>- **regional:** the warp, the strike-oriented relief bank, rock bands, the coast term, the micro-texture, glaciers at fine levels, the filter with its strength field, the projection.<br>It also settles what the pipeline already has instead of a sandbox proxy (climate, the regional solve); the fold fabric's fate (a switch, default off: T0.4); and belt zonation (S5.5's share). The sea-floor part waits for Phase 12. It includes the noise calibration of 12.0 (granularity per scale, directionality from the tectonic fields), designed once for above and below sea. Written into this phase. | design, surveyor | docs | 6.5 | RM T2.0, T0.4, MSP M10, M11 |
| 8.2 | Global erodibility from the tectonic fields. World-changing. | worker, reviewer | terrain | 8.1 | RM T2.1 |
| 8.3 | **Hypsometry** (the land's height distribution). First a diagnosis: land hypsometry after each stage of a 256 run. Then the uplands lift and the global erosion's levers, together with the denudation band against the diffusion-dominated rate. Sediment: about 400 times too little reaches the sea. `placerResidenceMyr` becomes a parameter. The user looks at the result. | worker, reviewer, user | terrain, hydrology | 8.2 | RM T2.2, EHP E3, E4, E5 |
| 8.4 | Regional: the warp, the oriented relief bank with the belt zones, the rock bands; the fold fabric behind its switch. | worker, reviewer | windows | 8.1, 7.4 | RM T2.4 |
| 8.5 | Regional: the coast term and the shoreline micro-texture, after the level's solve, against the open sea. | worker, reviewer | windows | 8.4 | RM T2.5 |
| 8.6 | Regional: the erosion filter after each level's solve, with its strength field (noise, rock, precipitation, none under ice), and the projection. This replaces S4.5 and the erosion plan's Phase F. | worker, reviewer | windows | 8.5 | RM T2.6, EHP E2, MSP M8 |
| 8.7 | Regional glaciers at fine levels (fjords, troughs), if the existing regional glacial forms do not give them. | worker, reviewer | windows | 8.6 | RM T2.7 |
| 8.8 | Judging runs: the world at 256; A, B and C at k 64; the k 256 windows. Compared against the sandbox's frozen pages and the Earth atlas. Then the user's verdict. The old regional items (coast ruggedness, the trough cap and wall, cirques against pits, the bed law, shore-floor block edges) are closed or carried forward. | orch, user | | 8.2 to 8.7 | RM T2.8, RM T5 regional items, EHP E26 |

### Phase 9: Tectonics realism

Goal: an engine whose oceans, margins and late history look like Earth's. It has its own lane, so it runs beside Phases 3 to 8. Every behaviour step gets a reviewer; world-changing results go to the user's look.

| Step | What | Agent | Lane | Depends | Origin |
|---|---|---|---|---|---|
| 9.1 | Findings on the cube engine:<br>- a cube corner's junction arm is lost;<br>- no imported-world run since the switch;<br>- the arc distance (re-check after the port);<br>- `crust_class` is read from plate rasters instead of the columns, so Andean margins are labelled collisions;<br>- re-measure H2's rift-edge under-measurement on the cube transport. | worker, reviewer | engine | 3.1 | EF E1, E2, E4, E10a, E10b; RM T4.8 |
| 9.2 | Survey and re-measure on baseline v6: spreading rates, oceanic crust share (19 % on the torus) and ridge ages now that plates rotate at 1.36 cm/yr; whether progressive rifting helps (a run and a look); plume behaviour with `plume_processes.count` > 0. | surveyor, worker | engine | | RM T4.1, T4.3, T4.10 |
| 9.3 | **Late-run drift.** After about 400 Myr the plates slow, the crust thickens and the shelves vanish: the analytic mantle forcing's 600 Myr clock. Design a continent-aware forcing (insulation under supercontinents, slab-driven assembly) so that durations beyond 400 Myr keep Earth's statistics. This is also the open part of the continental growth work (shelf renewal after 500 Myr). | design, then worker, reviewer | engine | 9.2 | RM T4.10, T4.12 (what is left), world duration note |
| 9.4 | **Plumes on.** Hotspot chains; plume positions given to the volcano stage (today `terrain_volcanoes` places hotspots at random); the plumes' visual review. Plume deposits go to 15.3. | worker, reviewer, user | engine, terrain | 9.2 | WM-5, WM-6, EF E8b, RM T6 hotspot chains |
| 9.5 | Back-arc rollback behind old, fast ocean-ocean subduction (back-arc basins, remnant arcs, island arcs, which never emerge today). World-changing. | worker, reviewer, user | engine | 9.2 | RM T4.2, EF E3 |
| 9.6 | **Tectonics preview tool.** Checkpoints every 10 Myr, a time slider over the frames, and the pipeline continued from a chosen checkpoint. It needs 14.2 (a bundle that can resume). | worker | engine, tools | 14.2 | RM T4.9 |
| 9.7 | Material coordinates and the full strain tensor exported; stabilization age and lithosphere keel thickness (cratons, for the mineral work); a persistent plate identity in the bundle. This is a contract change, decided by the orchestrator and recorded. | worker, reviewer | engine | 9.1 | RM T4.4 (S5.4), CP-8, MS-3, MS-7 |
| 9.8 | Drawing in material space, seam bands, belt zonation (what 8.1 leaves to the engine). | worker | engine, windows | 9.7, 8.4 | RM T4.5 (S5.5) |
| 9.9 | Thin viscous sheet: lateral flow, plateaus, curved belts, gradual rifting with rift flanks. Design, then steps; it ends with the user's verdict on a 256 run. | design, worker, user | engine | 9.7 | RM T4.6, T4.7 (S6) |
| 9.10 | Deep mantle dynamics: a coarse convecting mantle, dynamic topography, plume sources. Design only until 9.3 shows whether it is needed. | design | engine | 9.3 | RM T4.11, WM-10 |

Recorded, no step: N 32 has too few mountains because thickening depends on resolution (EF E12). N 32 stays the parity configuration only.

### Phase 10: Climate quality

Goal: Earth-like rain, evaporation and circulation on the cube world.

| Step | What | Agent | Lane | Depends | Origin |
|---|---|---|---|---|---|
| 10.1 | Rain and evaporation: the south-polar hot desert; potential evaporation 2 to 5 times Earth's; land rain 22 % high; residual ice flicker near the ELA; moisture delivery and rainfall geography (the heating change of 2.77 % against the 2 % tolerance). | worker, reviewer | climate | 3.2 | RM T5, EHP E14, E15, FW F5 |
| 10.2 | Köppen: a processed Köppen reference in `refs/processed`; spatial Köppen validation (0 cells compared so far; the C/D boundary convention); the ITCZ by month. | worker | climate, tools | | CP-7, AH A6, FW F8 |
| 10.3 | Tropical closure: the experimental mixed-layer closure (500 hPa return flow, tropical SST gradients) and the handoff's findings 1 to 5. Promote it or retire it. | design, worker | climate | 10.1 | FW F1, AH A1 to A3 |
| 10.4 | Seasonal winds and pressure: the weak Asian summer low, the Southern Ocean amplitude, Antarctic and east-Greenland vectors, basin pressure structure (diagnose heating gradients before drag tuning). | worker | climate | 10.1 | FW F3, F9 |
| 10.5 | Multi-level atmosphere (per-level temperature and humidity; 850 hPa loss 0.620 against 0.603). | worker | climate | 10.4 | WM-7, migration M8 |
| 10.6 | Validation and tests: the ocean solver's wet mask exported; the sea-ice reservoir; multi-seed and long-duration soak suites (only when asked); ablations and resolution studies; four `DISABLED_Modes` tests; `MpdataReversingWinds` without an assertion; five unsplit files. Conditional items: the polar CFL scheme and an ILU or multigrid preconditioner (if 3.2 shows they are needed), CFSR (needs reference data). | worker | climate | | FW F6, F7, F10, F12 to F15, F17; CP-28 |
| 10.7 | Bounded climate–hydrology feedback passes (README lifecycle); coarse-climate runoff at tectonic checkpoints (**deferred by the user**, 2026-10-07/08). | design | climate, hydrology | 10.1 | README R4, EHP E1 (D2c) |

### Phase 11: Hydrology and coasts

Goal: the remaining calibration and the small, known faults. Steps that change every world are batched into one step with one `step_check.py` and totals compared.

| Step | What | Agent | Lane | Depends | Origin |
|---|---|---|---|---|---|
| 11.1 | River mouths: Dai and Trenberth ds551.0 into `refs/source/` with a receipt, replacing the 20 authored mouths; add the Nile. | worker | hydrology, tools | | EHP E13 |
| 11.2 | One world-changing batch: `random(a, b)` draws only 15 bits (it biases `createfractal` low); `normalisecoasts` gets a reach in km instead of one ring; re-check the legacy lowering of coastal land (likely gone). | worker, reviewer | terrain | 8.3 | EHP E16, E17, E18 |
| 11.3 | Lakes and wetlands: three enclosed below-sea cells without a lake surface; delta and wetland thresholds (`river_factor`); `hydrology.deltas_and_wetlands` writes temperature; `terrain.arid_features` (random salt pans and ergs: keep, or leave to the lake balance). | worker | hydrology | | EHP E6, E9, E20; CP-11 |
| 11.4 | Known limits: snow and evaporation read temperatures before `mountain_lapse`; seasonal-flow rounding; volcano-island cells get no direction. Glacial risks: ice on lake cells, checkpoints without present-day ice, a one-step lag in the sea mask, calendar-quarter seasons, fjord heads walking inland, beds under ice caps. | worker | hydrology | | EHP E21, E22 |
| 11.6 | The cube world's own internal-drainage deficit: IoU 0.289 native against lat-lon 0.315, native endorheic share 0.080 against lat-lon 0.099 and Earth 0.177 (world-0088). First a cube-native hydrology evaluation (today the score reads the export). | worker | hydrology, tools | | 1.3 |
| 11.5 | Re-checks: the desert world (sea count 0) segfaulting in rainfall; `earth-climate` and `earth-hydrology` runs on the regenerated Earth data. | worker | hydrology | | migration M4, EHP E27 |

### Phase 12: The sea floor (waits for the user's research)

The user, 2026-10-13: work stops "until I come back with specific proposals". The research input is RM §"The sea floor: paused".

| Step | What | Origin |
|---|---|---|
| 12.0 | The user (2026-10-21): trenches look like paint smudges and are not realistic, being far from continental shelves (seen at (189, 231)-(255, 265) of `t24-base-v10/global-relief.png`, a straight trench along a row of plate-boundary cells). The user suspects much of what is wrong above and below sea is noise calibration: the granularity of each noise at each scale, and directionality from the tectonic fields (noise warped along the geographic features), plus some added noise. That calibration is shared with 8.1 (the warp, the strike-oriented bank) and is designed once for both. | | |
| 12.1 | The sandbox's sea floor follow-up: transform steps, narrow ridges, seamount chains (with 9.4). | RM T1.1 |
| 12.2 | Submarine processes: canyons, gullies, fans, drapes. | RM T1.8 |
| 12.3 | Sea floor and trenches as a global stage from crust age, replacing the calibrated C4 ocean hypsometry and its depth bands. | RM T2.3, EHP E7, E10, E11 |
| 12.4 | The drawing side of spreading rate and ocean area (the engine side is 9.2). | RM T4.3 |

### Phase 13: Large windows: tiles, GPU, 1 m

Goal: a 10 km window at 1 m per pixel. Every item is re-planned on the cube's level subsets first; the old numbers assumed lat-lon levels.

| Step | What | Agent | Lane | Depends | Origin |
|---|---|---|---|---|---|
| 13.1 | Tiles (`world_pipeline region --tiles DIR`). The coast step's reads outside `levelhalo` are a recorded risk. Tiles share their coarse levels. | worker, reviewer | windows | 7.4 | RM T3.4 (S4.1), MSP M17 |
| 13.2 | Kernel inventory test (per-cell functions only). | worker | windows | 13.1 | RM T3.5 (S4.2) |
| 13.3 | OpenGL 4.3 compute backend behind `--gpu`. | worker | windows | 13.2 | RM T3.6 (S4.3) |
| 13.4 | Close: a 10 km window of A at 1 m per pixel as tiles, CPU and GPU timings, a lidar comparison. | orch | windows | 13.3 | RM T3.7 (S4.4), MSP M7 |

### Phase 14: World model and desktop

Goal: runs that resume from a bundle, and a desktop on the cube.

| Step | What | Agent | Lane | Depends | Origin |
|---|---|---|---|---|---|
| 14.1 | Survey: what stops a run from starting at any stage of a bundle: tables, solver state, `GenerationScratch` (mountain drainage seeds), `--from-bundle` ignoring the bundle's parameters, `--parameters` unable to read a manifest. | surveyor | world | | WM-1, WM-2, CP-18, CP-22, CP-23, migration M2 |
| 14.2 | Resume from a bundle: tables and solver state stored, scratch into fields, per-step random streams. 9.6 needs it. | worker, reviewer | world | 14.1 | same |
| 14.3 | Desktop on the cube (today it generates lat-lon worlds, compile-only); progress and step notes on Generate; a parameter panel; the creating-world screen (never entered). The user's hand check, since the desktop has not been launched since Phase 0. | worker, user | desktop | 14.2 | RM:318, EHP E25, CP-15, CP-16, CP-19 |
| 14.4 | Dissolve the `world/functions.hpp` umbrella (31 includers); the regional map in floats (the `wholemetres` and `climatelegacyscale` bridges); the remaining whole-number getters. | worker | desktop, world | 14.3 | CP-13, CP-21, CP-25, SM-1 |
| 14.5 | Sketch and continuation initialisers (plates built from fields); copy-on-write snapshots; the imported-climate context that no stage reads. | design | world, engine | 14.2 | WM-3, WM-4, migration M1, M6 |

### Phase 15: Life, resources and society

Goal: the README's later lifecycle, on the finished terrain and climate.

| Step | What | Agent | Lane | Depends | Origin |
|---|---|---|---|---|---|
| 15.1 | Soils and vegetation (an ecology module), validated against land-cover data; canopy, fetch and snow burial; vegetation feedback to climate later. | design, worker | life | 8.8, 10.1 | resource chain, FW F4, README R5 |
| 15.2 | Fauna and livestock; biological resources. | design, worker | life | 15.1 | resource chain |
| 15.3 | Minerals beyond gold, and surface deposits (laterite, bauxite, evaporites, coal). Per-cell erosion histories instead of the denudation proxy (it removes 2 to 4 times the solver's erosion). Grade, tonnage and arc-rate calibration. Plume deposits (after 9.4). | design, worker | life | 9.7 | MS-1, MS-2, MS-4, MS-5, MS-6, WM-8, EF E5, TC-1 |
| 15.4 | Society as stages: settlement desirability per species, migration, settled societies and geopolitics (today one-shot heuristics, and `definitions/` is empty); the 0 to 100 resource scores. | design, worker | life | 15.2 | README R6, R7, SM-2 |

### Phase 16: Later phenomena (catalogue)

Each item becomes a design question when the user picks it.

- **Water and coasts:** deltas and alluvial fans at fine scale; endorheic basins and salt flats; karst; reefs and atolls; dune fields; meanders (a lateral-migration model); braided rivers and terraces; waterfalls and knickpoints; barrier islands, spits, lagoons, cliffs, marine terraces, tidal flats; ergs, loess, salt pans; wetlands.
- **Glacial:** moraines, drumlins, eskers, hanging valleys, calving, ice shelves, the climate seeing the ice.
- **Volcanic:** stratovolcanoes and calderas on arcs, shield volcanoes, flood basalts.
- **Structural:** great escarpments, half-grabens, salt domes, inselbergs, river capture, landslide scars, mesas, buttes, cuestas, hogbacks; permafrost.
- **Gameplay weather** (storms, fronts, runtime integration; deferred by the user).
- (Origin: RM T6, EHP E23, FW F2.)

### Phase 17: Reviews and the adversarial review

- Behaviour steps not yet reviewed: S5.3b, S5.3c, H4b.
- The adversarial review (GPT-6 Astra through Codex, run by the user) comes after the last phase this plan runs, holistically. Ranges: erosion `16ddfb3..`, multiscale `82fdef0..`, sandbox `0764d05..`, T8 `5edc651..`, and this plan from its commit on.
- (Origin: RM T7, EHP E28, MSP M23.)

## 4. Status

One row per step as it closes. Duration is measured from the previous step's commit to this step's commit (`AGENTS.md`).

| Step | Status | Commit | Duration |
|---|---|---|---|
| 1.1 | done: base-v10 adopted by the user (2026-10-20) | this plan's commit | |
| 1.6 | done: 14 items to the trash (`D:/dev/undiscovered-worlds/trash/20261020T190214Z-0`); the 30 T8 branches dropped after a verified bundle (`worktree-archive/branches-20261020-210153-0.bundle`) | (docs only) | |
| 2.8 | done: field lists for the brief, the report, the verdict, questions to the user, the wave handoff, status rows and ledger entries (`agent-workflow.md`, the four agent definitions, AGENTS.md durations); `scripts/check_records.py` checks strict ledger entries (above climate-0003, erosion-0106, tectonics-0020, world-0087) and §4 rows (outside the eleven older ones) in the pre-commit hook and `step_check.py`, and reports with `--report`; 33 tests; reviewed (pass with notes, all ten fixed) | this commit | worker 6 min; 7 min from 9b03453 |
| 2.6 | done: test-only scripts kept (their tests exist; `climate_reference_experiment_prepare.py` is the only producer for `climate_experiment_input.hpp`); the `map_rendering` and `cell_grid` aliases removed, 16 imports go to `uwtools`; new tests: `storage.route` refuses a non-empty destination, `process_access.py` fails on an unmeasured report and a stale include, `RecipeAudit` (C++) counts an undeclared write; RecipeAudit 2/2, Python 37/38 (the error is 2.9, worktree only). The sweep's seven unreferenced scripts (`scripts/refs` download/verify/prepare CLIs named in the local `refs/README.md`, `prepare_era5_reference.py`, `prepare_seaice_reference.py`, `sandbox/crop_views.py`) are kept: the first four are documented CLIs, the two preparers were maintained on 2026-10-20 | 9b03453 | worker 7 min; 10 min from 55d90db |
| 2.5 | done: the backup was already gone (deleted by the user); `storage.py finish` now handles an absent backup (junction verified, receipt marked complete, nothing removed); `receipts/out__scratch.json` complete | 55d90db | 15 h 47 min from d3c2d61 (session idle until 18:41; orch 3 min) |
| protection (new) | `protect.py lock` applied 2026-10-20 to refs/ (C and D sides), D:/dev/bak, worktree-archive and the trash; every place locked and the permission reaches the existing files (`protect.py status`). The refs writers were converted to copy-based writes (`72a6dd8`, `protected_io.staged`, staging on D) and only compiled: the first reference-data step that runs them (10.2, 11.1) checks them | d333d3d, 72a6dd8, 43debdf | |
| 1.2 | done: run 79 (cube N 64, baseline v6, 400 Myr) registered and the tectonics baseline; `experiment.py` read the cube recipe wrong (KeyError since baseline v5); tectonics-0021 | 5fea750 | worker 26 min (incl. run 79); 72 min from 78272f0 |
| 1.3 | done: done-when already met (export endorheic 0.065, native 0.080); mouths equal; IoU gap two thirds the cube world's (11.6); acceptance page refreshed; world-0088 | efbeb7f | worker 7 min + runs; 57 min from 78272f0 |
| 1.4 | done: eight documents name current paths; the dated sections of `consolidation-plan.md` keep their names (a record); `FUTURE_WORK.md` still names `summarize-climate-validation.py`, which no longer exists | 71ec78e | mechanic 5 min; 55 min from 78272f0 |
| 2.1 | done: `check_structure.py` with `configs/structure.json`, in `step_check.py` and a PreToolUse hook on Write/Edit for new files; 9 tests; 72 MB of legacy `out/` folders trashed | 0c8b3e6, 0217851 | worker 5 min |
| 3.1 | done: the transport's fields in parallel, bit-identical; 455 s idle against the 715 s gate (761 s before); reviewed; tectonics-0022 | 1c34a01 | worker 46 min, review 3 min; 101 min from 78272f0 |
| 3.2 (diagnosis) | done: residual diagnostics (log only, parity identical); the fix is decided in the row; climate-0004 | 2bd4b5d | worker 48 min; 99 min from 78272f0 |
| 5.1 | done: base-v10 baseline; `straight_rivers.py` measures cube domains; erosion-0107 | 74dd111 | worker 59 min; 110 min from 78272f0 |
| safeguards (new) | `guard_catastrophe.py` (global), `delete_guard.py` fails closed, `trash.py`, `worktrees.py drop-branches`; agents retire landed worktrees (the user, 2026-10-20) | 2288474 | |

## 5. For the user

What an agent cannot or may not do alone. The delete guard's refusals land here with path and reason.

- Done 2026-10-20: the user retired the 30 T8 worktrees. Since then agents retire landed worktrees, drop bundled branches and trash scratch themselves (AGENTS.md "Autonomy"; the two guards, `guard_catastrophe.py` and `delete_guard.py`).
- **Emptying the trash** (`D:/dev/undiscovered-worlds/trash`; `trash.py list` shows the entries): the user's only. The trash is delete-protected, so first `uv run python scripts/uwtools/protect.py unlock D:\dev\undiscovered-worlds\trash`, delete, then `protect.py lock`. The same applies to deleting anything by hand in `refs/`, `D:/dev/bak` or `worktree-archive`.
- **Harden the protection (one elevated step):** from a terminal opened with "Run as administrator", in the repository: `uv run python scripts/uwtools/protect.py harden`. That makes Administrators the owner of the protected folders, so removing the deny-delete permission needs an elevated terminal, which agents never get. Before it, only the guards stop an agent from removing the permission.
- `C:/Users/Luka/.claude/hooks/guard_catastrophe.py` is a two-line stub left by the move into the repository; the user may delete it.
- Decisions carried from the old plans (each has a recommendation; nothing blocks on them):
  - the verdict on climate baseline 388 (CP §10);
  - rewriting the draft review checklists (climate, tectonics, world, erosion, hydrology, regional);
  - `refs/processed/climate/{csv,fields,maps}` (12.7 GB): keep or move;
  - the difference maps' half-ranges;
  - the `.uwclim` inputs of `climate-reference-experiments-20260915`;
  - `C:/dev/plate-tectonics`;
  - whether the verifier resumes;
  - tributaries' default (4.1);
  - reach 0 against 1 (7.4);
  - coarse-climate runoff at checkpoints (10.7, deferred).
- Looks: 5.4 (routing), 6.1 and 6.5 (the target look), 8.3 and 8.8 (the port), 9.4, 9.5 and 9.9 (the engine).
- The sea floor research (Phase 12).

## 6. Dropped as stale (with the reason)

- `generateglobalclimate` (RM T5): deleted in `0fd75c5`. The `connectlakes` bug: deleted. The climate maps README: exists and is regenerated.
- FastLEM's 16 CPU minutes (migration M5): FastLEM was removed.
- `docs/tasks/main/PLAN.md`: the reorganisation finished. `ASTRA_HANDOFF.md`'s workflow constraints are overridden by `AGENTS.md`; its scientific findings are in 10.3.
- The torus findings: plates translating without rotation (T4.10's premise: the cube engine rotates plates); the north-south wrap; the 1.57 area excess; `torusmetric`; the engine's planar noise; the pre-cube striping. All are gone with T8.
- The lat-lon-only routing findings (the polar-row D8 band, the east-west routing share; EHP E24): the cube is the default. They apply only to `--grid latlon`.
- The old pending verdicts (regional runs 1 to 15; multiscale runs 37 to 57; S3.13's lake outline; the phasor or Gabor verdict): superseded by the sandbox (T0.1).
- `noise-functions.md`'s diamond-square in terrain and the 4-D torus noise: terrain moved to sphere noise in T8. Climate still uses `createfractal`, which 11.2 covers.
- The water-balance scratch grids (migration M3): overtaken by `lake_levels` (unverified, re-checked in 14.1).

## 7. Mapping from the old ids

| Old | New | Old | New | Old | New |
|---|---|---|---|---|---|
| T1.1, T1.8 | 12.1, 12.2 | T3.1-T3.3 | 7.2-7.4 | T4.9 | 9.6 |
| T1.3 | 6.4 | T3.4-T3.7 | 13.1-13.4 | T4.10 | 9.2, 9.3 |
| T1.4 | 6.3 | T4.1, T4.3 | 9.2, 12.4 | T4.11 | 9.10 |
| T1.6 | 6.2 | T4.2 | 9.5 | T4.12 (left) | 9.3 |
| T1.7 | 6.5 | T4.4, T4.5 | 9.7, 9.8 | T5 | 8.3, 10.1, 11.x |
| T2.0-T2.8 | 8.1-8.8 (T2.3 → 12.3) | T4.6, T4.7 | 9.9 | T6, T7 | 16, 17 |
| T4.8 | 9.1 | item 3 | 5 | T8 open items | 1.2, 1.3, 3, 4 |
