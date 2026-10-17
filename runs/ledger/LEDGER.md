# Experiment ledger

Generated from `runs/ledger/*.jsonl` by `scripts/uwtools/ledger.py`; do not edit. Show entries with `uv run scripts/uw.py <module> ledger [--topic TEXT]`.

## Settled / don't retry

- `erosion-0014`: The rainfall solver's bulk-aerodynamic potential evaporation as the Budyko PET (discharge recall on Earth fell from 0.63 to 0.46).
- `erosion-0016`: Depth from the engine's crust age by the plate-cooling law: the engine's ocean ages are mostly the run's age (median 590 Myr at 600 Myr) and two thirds of the sea lies on continental crust, so 98 percent of the ocean lands at 5 to 6 km.
- `erosion-0021`: Loading fraction 2700/3300 with rebound 2800/3300 on rock volumes (it raised the land by 0.66 m on average in 20 Myr).
- `erosion-0031`: Routing on the ice surface during the run: ice flickers on and off between steps (643 appearances on 305 cells in 20 steps at the minimal size).
- `erosion-0036`: stopping the fixed point when the receivers settle: the decision covers the whole window, so overlapping windows stop at different passes
- `erosion-0037`: ownership of a shared fine cell by the larger drainage area: the bed can rise where a gentle tributary meets a steeper larger link
- `erosion-0043`: coast noise on the elevation of the coastal band (doubled coastal roughness and channel trenches: A coastal channel q1 -740 to -1585 m)
- `erosion-0043`: a shore floor not capped by the global cell's own distance from the sea level (C's block-mean error 29 m)
- `erosion-0044`: springs on every source cell by default (pillars 21 / 166 / 200 -> 24 / 220 / 272)
- `erosion-0045`: the lake zone chosen per global cell (water steps up to 168 m between adjacent wet cells)
- `erosion-0045`: an indicator blind to the terrain (beds carved up to 2.9 km, shores raised up to 744 m)
- `erosion-0049`: the valley floor in fine cells (a floor four times wider on the ground at k 4 than at k 16): kilometres instead
- `erosion-0049`: a side slope of 0.02 at 16 cells of reach (the land around the rivers flattened: block mean -483 m on B)
- `erosion-0050`: the island sink at the legacy margin 0.5 (inert: the ring meets the mainland)
- `erosion-0050`: coast noise 150 m / 0.7 for the default (islands in global sea nearly double, C's block-mean error 15.8 m)
- `erosion-0050`: glacial 100 or denudation 50 (C's block-mean error 15.5 m)
- `erosion-0051`: the sign of the interpolated sea indicator as the gate (narrower than the 3 x 3 rule; cannot open (184, 49))
- `erosion-0051`: the smooth gate with the 3 x 3 lift (the coast stays within one fine cell of the block edge)
- `erosion-0051`: gate reach 2 (158 low-land cells two cells from the sea on A)
- `erosion-0051`: straight_max_cells 6 (runs of 5 and 6 stay)
- `erosion-0054`: trough_width_km 10 or 20 at the 40 km cap (alike; the land's mean falls 63 m on B)
- `erosion-0055`: the candidate's band test on the solved fine elevation (not seam-safe in the innermost halo ring; the review of E5.0)
- `erosion-0058`: exporting accumulated_strain or deformation_rate as shortening: both are unitless 0..1 signals on the global grid with no direction and no advection
- `erosion-0058`: carrying the angle through the remap: angles do not average; carry the tensor components
- `erosion-0060`: reusing the map kind 'reference' for Earth windows: it is the climate validation's and its caption says 'grid of <run>'
- `erosion-0061`: a halo of a constant number of parent cells at every level: the channels and the solve reach 2 global cells whatever the level until S3.1
- `erosion-0066`: sea as below 0 m, or as below 0 m connected to the window's edge: the Caspian is below 0 m and crosses window edges
- `erosion-0068`: weights of 0 on pinned cells with the remainder on the land children: one dry child in a lake parent took 16 lambda (-5.8 km holes)
- `erosion-0068`: assigning a detail to the finest level run (detaillevel clamped to L): levels then depend on L and a k = 16 run averaged to k = 4 differs from a k = 4 run
- `erosion-0068`: the valley band on floor plus side reach: the reach is in level cells, so every valley with sides fell outside [1, factor) at every level
- `erosion-0071`: factor 2: with 4 children per parent the bounds exhaust (237 forced parents at level 4 on C) and steps reach 21.5 km
- `erosion-0072`: the projection's per-parent lift without a smooth correction before it: the solve's per-parent lowering of kilometres prints the parent grid
- `erosion-0073`: a hard clamp of children of land parents to the sea level: their parents (coastal plains about 2 m above the sea) cannot give the volume back, so it makes flats and walls
- `erosion-0073`: unlimited room down for land below the sea: pits to 3.3 km
- `erosion-0075`: treating the level-2 pits as the prolongation's undershoot: the bracketing parents bound them far above the pits
- `erosion-0076`: de-lifting pinned beds (all, or those below the filled cap): the forms-off lambda doubles and level-3 pits pass 150 m
- `erosion-0076`: the C2 B-spline smoothing basis: borders print (step max 481 m)
- `erosion-0077`: a strict rim cap: 46 to 149 dry cells below a lake on C16 L2
- `erosion-0077`: letting the water retreat behind a capped rim: 39 percent of C16's lakes drained
- `erosion-0078`: a valley reach of 10 or 20 km: stepped double lines at the reach's edge
- `erosion-0079`: a tributary bed rising to the surface's running minimum along the whole path: border steps to 1.8 km at A64 level 3
- `erosion-0080`: a halo that lets a global-cell reach shrink with the level (levelhalo before S3.8)
- `erosion-0081`: one multi-seed Priority-Flood with a reach predicate as the catchment labelling: a race whose labels flip along corridors of any length
- `erosion-0081`: labelling by (filled level, hashed seed key) over local floods: up to 13 percent of cells unreached
- `erosion-0082`: transit edge points kept off the shared corner: the peaks rise (C16 6.6 to 8.6)
- `erosion-0082`: an unqualified transit block: trenches, levees and dams (48 m steps, straight lake arms)
- `erosion-0084`: a confluence rule that recomputes ownership to a fixed point (no fixed point in a ring of three crossings; dependency depth 3 to 10 global cells on B at k 16)
- `erosion-0088`: flooring the lake cell at S (b(lake cell) = S) with the lake depth applied to the lake-block channel cells: a cliff at the block edge (98 m level-2 border steps on B)
- `erosion-0088`: holding a coarse wet channel cell at its bed (332 m steps): the mean law at every level
- `erosion-0088`: capping a wet channel cell's mean at its surface at the final level: it would carve B's fjord walls into trenches
- `erosion-0089`: coast_warp_persistence above 0.5 with the finest octaves on (specks); 0.25 cells of displacement (dimension 1.52, block-mean error 13.6 m)
- `erosion-0089`: a pocket core anchored on the interpolated global elevation alone (103 enclosed components) or a core 20 m deep (91): the one-sided rule's base is the anchor
- `erosion-0090`: raising filled pocket cells beside a lake to the shore floor (the lake loses its bed)
- `erosion-0091`: handing a filled pocket's cells to the lake by any per-cell or per-component rule in the coast step while the levels >= 2 flood through the parent's sea gate: level 1's coarse flood is inherited as rectangles
- `erosion-0093`: kernels of width L on a 78 km tile with 8 kernels (28 km apart): the fabric breaks into patches with abrupt fold offsets
- `erosion-0093`: the wavelength from the sediment cover alone: a bare orogen gets the minimum wavelength, a corduroy at two pixels at k 64
- `erosion-0093`: a per-kernel phase origin shared per global cell as the remedy for the along-strike steps: p99 529 -> 520 m on the belt
- `erosion-0095`: a kernel sum with per-splat phases as the fold fabric, whatever the kernel layout or weight: the crests break at 0.17 A on a uniform field
- `erosion-0096`: a vergence sign on k from the strain gradient: a kink drawn as an 80 km fold at every flip
- `erosion-0096`: chasing the along-strike steps with the axis smoothing above 3.5 cells: 2.5 to 3.5 moved p99 by 1 m
- `erosion-0097`: the fine relief, coast terms and filter drawn on unwarped coordinates (v15, with transverse faults, plunging ridges and lower stretches): an isotropic orange peel that broke the directionality; rejected by the user
- `erosion-0097`: the coast term applied again after the physics (after_physics): it put every lowland valley floor into its band and flooded G's coastal plain and F's southern valleys
- `erosion-0097`: the coast term's band chosen by height alone: it drowned G's low inland valleys (it now acts within reach_km of the base's coastline)
- `erosion-0097`: the shoreline micro-texture drawn against all of the sea: flooded glacial troughs became shores and G-coast's troughs turned into blotches (the shore is the open sea less parts narrower than 6 km)
- `erosion-0097`: a vertical micro-texture before the physics: 80 m barely moves a shore, and 20 Myr of physics erases few-km shore detail drawn before it
- `erosion-0097`: the coast term at 60 % of the warp (coast.warp_fraction, kept as a key at 1)
- `erosion-0097`: an unsmoothed strain-driven warp: a 250 km drag folded the grid over itself along the belts
- `erosion-0098`: equal fracture-zone troughs at every segment boundary: they read as a comb
- `erosion-0098`: the smoothed convergent-boundary band drawn as the trench: beaded troughs hundreds of km wide
- `erosion-0099`: starting a procedural-world terrain chain from runs/tectonics/78/bundle: terrain.tectonic_trenches reads a table a bundle cannot load
- `erosion-0100`: a long lowstand at -120 m (3 Myr) as the only sea-level history: it drains inland lakes; use the Earth-shaped history
- `erosion-0100`: flat lake fills: the Quaternary ice cuts comb lakes along their D8 lines
- `erosion-0100`: fill noise at a few km on the world's 19.5 km grid: it fades out and leaves straight radial drainage
- `erosion-0100`: the tectonics sea_level key to change the final land fraction: it only sets the initial ocean coverage
- `erosion-0100`: placing the sea at Earth's land fraction while the continents are thin: the interiors flood
- `erosion-0101`: raising the child's detail amplitude beyond the Earth-fitted law to fake steepness the parent lacks
- `erosion-0101`: a steeper window landslide threshold as the fix for the parent's low slopes (it only moves p99)
- `erosion-0101`: longer child physics (2 Myr) to grow valleys from noise
- `erosion-0102`: damping platform relief to fix the foothills (platform.relief_damp 0.4 barely changes them)
- `erosion-0102`: a steeper or no landslide threshold at k 64 (no effect)
- `erosion-0103`: a higher fluvial K on an unnested window (it drains the land; only the nested projection can hold the mean)
- `erosion-0103`: a uniform per-parent-cell shift as the exactness step after the physics (it draws the grid)
- `erosion-0103`: the lake trap on every lake under ice (it fills the troughs that become fjords)
- `erosion-0104`: showing nested windows to the user before the parent's D8 bias and the block terraces are gone
- `erosion-0104`: a global-scale page for one-pixel features (crop at 1:1)
- `erosion-0105`: comparing a sandbox label's config.json before its render has finished (written at the end)
- `tectonics-0010`: Bilinear semi-Lagrangian sampling of the elevation for the uplift rate (spurious uplift at moving summits).
- `tectonics-0011`: carrying the strain into new crust at divergent boundaries: young ocean floor showed 8 km of shortening
- `tectonics-0012`: measuring continental area with crust_class (plate rasters): use the columns or continental_crust thickness
- `tectonics-0012`: expecting the minimum-thickness cap and collisions alone to regulate thickness: nothing bounds it from above
- `tectonics-0013`: calibrating continental growth at 128 columns and expecting 256 to follow
- `tectonics-0013`: sweeping h_max (32-40 km) or D (10000-30000) to fix the land height at 256 columns: the Airy slope sets it
- `tectonics-0013`: an area-cap threshold below the collapse threshold to make shelves (margins are sharp at 256 columns)
- `tectonics-0014`: a frozen root cap as the way to high mountains: it keeps the root under the orogen
- `tectonics-0014`: placing the land fraction at terraforming: the engine surface has no shelves yet
- `tectonics-0014`: tuning h_max or D to fix the land height: the Airy slope sets it
- `tectonics-0015`: an onset with no collapse at all between the cap and the onset (loses 0.03-0.07 of continental area)
- `tectonics-0015`: judging the Earth-like set at 600 Myr without the 400 Myr checkpoint
- `tectonics-0016`: the margin rule alone (W1, W2: interiors grow past 60 km without bound, 3000+ 3.5-7 %)
- `tectonics-0016`: creep 20000 below a 65 km onset (W4: no mountains)
- `world-0003`: Re-enabling terrain.fastlem_mountains as it is: it does not terminate at seed 1, 256 columns.
- `world-0007`: Regenerating the access file from 256-column runs only: it drops writes that occur at the default size.
- `world-0035`: rewriting `accessor(...) > sealevel` without excluding an offset after sealevel or arithmetic before the accessor
- `world-0054`: converting the regional generator to float reads by pattern without launching the desktop application
- `world-0086`: Creating a junction from a worktree to the main tree's refs/ (or any data directory).
- `world-0087`: the valley-spacing law with the erosion's grid-scale D and K: spacings of hundreds to millions of km

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
| `erosion-0024` | 2026-10-07 | changed as intended | keep | erosion plan, step H3: enclosed below-sea basins are depressions for the lake balance, not sea |
| `erosion-0025` | 2026-10-07 | changed as intended | keep | erosion plan, step H4a (C10, first half): eroded-world becomes procedural-world; imported worlds on the new hydrology |
| `erosion-0026` | 2026-10-07 | equivalent | keep | erosion plan, step H4b (C10, second half): the legacy hydrology deleted |
| `erosion-0027` | 2026-10-07 | pending user review | keep | erosion plan, step H5: the Phase H close (erosion run 8, the promoted procedural-world at 256 columns) |
| `erosion-0028` | 2026-10-08 | pass (contracts; no behaviour change) | keep | erosion plan, step G0: the Phase G survey, design and contracts (glacial erosion) |
| `erosion-0029` | 2026-10-08 | pass | keep | erosion plan, step G2: terrain.fjords deleted |
| `erosion-0030` | 2026-10-08 | pass with notes (review on Opus 5.5) | keep | erosion plan, step G1: the glacial kernel (ELA, ice flux, thickness) and the glacial term of the stream-power step |
| `erosion-0031` | 2026-10-08 | pass with notes (review rework, re-review pass with notes; both on Opus 5.5) | keep | erosion plan, step G3: the glacial term in terrain.fluvial_erosion, its fields, evaluation and access |
| `erosion-0032` | 2026-10-08 | pending user review | keep | erosion plan, step G4: the Phase G close (erosion run 9, glacial erosion at 256 columns) |
| `erosion-0033` | 2026-10-08 | pass | keep | erosion plan, step E1: the legacy regional generator leaves src/app (Phase E, headless regional maps) |
| `erosion-0034` | 2026-10-08 | pass | keep | erosion plan, step E0: survey, the Phase E table and decisions, the contracts of the regional pass |
| `erosion-0035` | 2026-10-08 | pass | keep | erosion plan, step E2a: the window grid, the spherical noise and the upsampling of the regional pass |
| `erosion-0036` | 2026-10-08 | pass | keep | erosion plan, step E2c: the analytical relief of the regional pass |
| `erosion-0037` | 2026-10-08 | pass | keep | erosion plan, step E2b: the pinned channel network of the regional pass |
| `erosion-0038` | 2026-10-08 | pass | keep | erosion plan, step E3: the regional pass (planet-free) and its planet adapter |
| `erosion-0039` | 2026-10-08 | pass | keep | erosion plan, step E2d: the regional bundle, the region subcommand and the evaluation module regional |
| `erosion-0040` | 2026-10-08 | pending user review | keep | erosion plan, step E4: the Phase E close (regional runs 1 to 4 on erosion run 9's bundle) |
| `erosion-0041` | 2026-10-08 | pass | keep | erosion plan, step R0: the contracts of Phase R (regional realism) |
| `erosion-0042` | 2026-10-09 | pass | keep | erosion plan, step R1: channels end at the fine coast, beds follow the start surface, rivers wet in lakes, map outlines |
| `erosion-0043` | 2026-10-09 | pass | keep | erosion plan, step R2: fractal coasts by a warp, a coast noise on the sea decision and faithfulness by land share |
| `erosion-0044` | 2026-10-09 | pass | keep | erosion plan, step R3: rivers after the legacy scheme, inside their own global cells (junctions, edge points, Catmull-Rom legs) |
| `erosion-0045` | 2026-10-09 | pass | keep | erosion plan, step R4: lake shores and ice edges from signed indicators at warped positions |
| `erosion-0046` | 2026-10-09 | pending user review | keep | erosion plan, step R5: the Phase R close (regional runs 5 to 8 on erosion run 9, against runs 1 to 4) |
| `erosion-0047` | 2026-10-09 | pass | keep | erosion plan, after Phase R: the user's review of regional runs 5 to 8; lakes over rivers and a fading-river terrain map; Phase R2 planned |
| `erosion-0048` | 2026-10-09 | pass | keep | erosion plan, step R6: the contracts of Phase R2 (regional realism, second round) |
| `erosion-0049` | 2026-10-09 | pass | keep | erosion plan, step R9: channels incise by the world's uplift or denudation rate, the solve's uplift term, valley floors from the flow |
| `erosion-0050` | 2026-10-09 | pass | keep | erosion plan, step R7: islands and coast ruggedness (one-sided coast noise, the island sink, ruggedness from glaciation and denudation) |
| `erosion-0051` | 2026-10-09 | pass | keep | erosion plan, step R8: grid artefacts on the coast (the smooth sea gate, the straight-run breaker) |
| `erosion-0052` | 2026-10-09 | pending user review | keep | erosion plan, step R10: the Phase R2 close (regional runs 9 to 12 on erosion run 9, against runs 5 to 8) |
| `erosion-0053` | 2026-10-09 | pass | keep | erosion plan, step E5.0: the contracts of Phase E5 (glacial forms on the regional map) |
| `erosion-0054` | 2026-10-09 | pass | keep | erosion plan, step E5.1: glacial troughs (the U-shaped valleys) of the regional pass |
| `erosion-0055` | 2026-10-09 | pass | keep | erosion plan, step E5.2: cirques on the regional map |
| `erosion-0056` | 2026-10-09 | pass | keep | erosion plan, step E5.3: present-day ice above the equilibrium line on the regional map |
| `erosion-0057` | 2026-10-09 | pending user review | keep | erosion plan, step E5.4: the Phase E5 close (regional runs 13 to 15 on erosion run 9, against runs 9 to 11) |
| `erosion-0058` | 2026-10-10 | survey complete | S5.1's row and its user question are written with this design; the question is asked in section 7 at S0.5. | multiscale plan, step S0.3: survey of the tectonics engine's strain (for S5.1's shortening export) |
| `erosion-0059` | 2026-10-10 | survey complete | S0.4 briefed with the projection placed after the cirques and before water and climate, unweighted parent means (as correctblockmeans), and the manifest fields optional in the window block. | multiscale plan, step S0.1: survey of the regional pass's data flow |
| `erosion-0060` | 2026-10-10 | survey complete | S1.2 adds an optional hillshade layer to the regional terrain map and the Earth window alike (recorded in section 6), and Earth windows get their own map kind 'earth' and column, never the climate validation's 'reference'. | multiscale plan, step S0.2: survey of the evaluation renderers and the ETOPO 2022 15-arc-second tiles |
| `erosion-0061` | 2026-10-10 | contracts in place | LevelStep (not LevelDetail: projection and water are steps, not details); projection_weight_channel and detail_cells are floats; level 0 of a halo chain has halo 0 (the global fields are given); the halo of level L is target.halo + 2 (factor^(L-1) - 1) parent cells, so level 1 is today's domain; detail_cells 2 (a feature resolved by two cells of its level), tunable at S2.6 or S3.0; resample at source level 0 with a whole n >= 2 is reported as resampled although its pixels are exact means of whole global cells (safe; intended). | multiscale plan, step S0.4: contracts of the refinement hierarchy (refinement_grid, refinement_projection, refinement_levels, landscape.refinement, test files, optional manifest fields) |
| `erosion-0062` | 2026-10-10 | phase closed | S1.1 and S1.2 in the main tree; S2.1 and S2.2 in worktrees; Earth windows get their own map kind 'earth' and a new 'relief' map with a hillshade beside the unchanged 'terrain'; S5.1 route (a), API and world fields like rock_uplift_rate, recommended to the user. | multiscale plan, step S0.5: Phase S0 close (survey findings into the plan, S1 and S2 rows refined, the S5.1 question) |
| `erosion-0063` | 2026-10-10 | grid implemented | levelwindow also refuses a domain extent above INT_MAX and requires 1 <= columns <= W and rows >= 1 (the worker's addition beyond the contract, so int fields cannot overflow silently). | multiscale plan, step S2.1: refinement_grid implemented and tested (level windows, outward snapping, the halo chain and its report, regionalwindowof) |
| `erosion-0064` | 2026-10-10 | projection implemented | The orchestrator added comment-only clarifications to refinement_projection.hpp from the worker's open questions (border pairs use the weights as applied, pixelcells is n as computed, the two added resample exceptions). | multiscale plan, step S2.2: refinement_projection implemented and tested (weights, the exact projection with its report, the parent-mean error, the resample) |
| `erosion-0065` | 2026-10-10 | reference data in place | The three extra catalogue windows stay (Himalaya in Nepal, the Alps, the Basin and Range). | multiscale plan, step S1.1: the ETOPO 2022 tile module, the tile downloader and the first Earth-window catalogue |
| `erosion-0066` | 2026-10-10 | Earth windows on the pages | The hillshade's vertical exaggeration is 10 (300 m of relief over a 10 km pixel gives about 15 percent brightness); slopes steeper than 1 in 10 saturate, so finer pixels may want a smaller value (a profile option). Sea in an Earth window is the world ocean, not 'below 0 m'. | multiscale plan, step S1.2: Earth reference windows beside the regional maps, and the hillshaded relief map |
| `erosion-0067` | 2026-10-10 | phase closed | Entries: A mariana_arc (nearest catalogue latitude), B norwegian_fjords (glaciated land), C zagros (the plan's choice). | multiscale plan, step S1.3: Phase S1 close (regional runs 13 to 15 evaluated again with the relief map and Earth windows) |
| `erosion-0068` | 2026-10-10 | hierarchy in place | The parent-mean bound is half a float step at the kept box's largest |elevation| (1.95e-3 m on run 9, whose sea level is 31,043 m above datum), not 1e-3 m. The valley band uses max(floor half-width, reach) for its lower end and the floor half-width for its upper end. A ground reach for the valley sides, a limited prolongation in extreme relief and the solve per level against its parent are S3.2/S3.3 items. | multiscale plan, step S2.3: the loop over refinement levels (runrefinement), the per-child projection bounds, detail levels by size band; the regional pass calls the loop for powers of the factor |
| `erosion-0069` | 2026-10-10 | resample in place | Lake surface of a pixel is the weighted mean over its wet cells (zeros would put the surface below the bed); M at or above the global cell size is refused rather than served from level 0. | multiscale plan, step S2.4: world_pipeline region --metres-per-pixel M (the hierarchy to the resample's source level, box-averaged to M) |
| `erosion-0070` | 2026-10-10 | checks in place | level-consistency's allowance is the sum of the finer run's level bounds over the levels between (the coarser run's own level equals the finer's bit for bit, asserted by S2.3's Consistency test). | multiscale plan, step S2.5: evaluation checks parent-mean-exact, level-consistency, resolution-invariance; the per-level table on the run page; pixel-grid bundles in the Python readers |
| `erosion-0071` | 2026-10-10 | contract met; visible parent blocks | factor 4 (recommendation and orchestrator's choice; the user's verdict on the factor is requested in section 7 with this evidence). A step S2.7 is added before Phase S3: a smooth correction before the exact projection (design agent launched), and the k 64 defects (sea holes along channels, square lake edge, staircases). | multiscale plan, step S2.6: Phase S2 close (registered judging runs at k 4, 16 and 64, factors 4 and 2) |
| `erosion-0072` | 2026-10-10 | blocks removed; ripple left | smoothparentmeans's signature and the four bundle keys accepted; 'forced 0' waived to S2.8; registration moved to S2.8. | multiscale plan, step S2.7: smooth before exact (smoothparentmeans before the projection) |
| `erosion-0073` | 2026-10-10 | phase S2 closed | No unlimited sinking room for land below sea; the counted forced fallback stays (forced moves at most 10 m, located); the straight sea-lake contact and the rim raises go to S3.3. | multiscale plan, step S2.8: sea and lake masks inherited from the parent and refined within one parent cell; the judging runs of S2.7 and S2.8 (Phase S2's close) |
| `erosion-0074` | 2026-10-10 | contracts in place | Keys under the regional. prefix; the targets apply from level 2; lakeedgealigned over the whole kept box is S3.4's measure; Génevaux et al. 2013 use no Horton ratios (the plan's section 2 corrected in section 6). | multiscale plan, step S3.0: Phase S3's contracts (the per-step hook interface of the refinement loop, new headers, keys, diagnostics), after the survey of the valley-spacing law and the re-plan of Phase S3 |
| `erosion-0075` | 2026-10-10 | clamp available, off; targets moved to S3.2 | refinement.prolongation_limit 0 until S3.2's measurement; S3.1's pit, rim-raise and forced targets move to S3.2. | multiscale plan, step S3.1: limited prolongation (each child clamped to its 2 x 2 bracketing parents), off by default |
| `erosion-0076` | 2026-10-10 | pits removed; stripes are the forms' | S3.2's lambda and smoothing-move targets apply with the forms off until S3.5; the rim raise is S3.4's; prolongation_limit 1 (a reversible default: shallower pits, slightly larger rim raise); smooth_basis 0 (the B-spline prints borders); the lift taken off unpinned cells only (the bed double count stays: re-rasterising the beds on the de-lifted surface is a later item). | multiscale plan, step S3.2: the solve against its parent, channel cells at their cell mean, the projection's bed floor |
| `erosion-0077` | 2026-10-10 | lake edges and rims improved | The overflow rule kept (join the water, else the old raise); the rim-raise target re-judged at S3.7; regional.noise_relief_field stays 0 (S3.6's relief amplitude is not on localrelief's scale). | multiscale plan, step S3.4: noise amplitude per level from the relief field, the lake-sea contact at a sill, the rim raise capped by the level's noise band |
| `erosion-0078` | 2026-10-10 | stripes and terraces gone | regional.valley_reach_km 40; S3.2's lambda and smoothing-move targets dropped as acceptance criteria (not reachable from the forms' files; forms off already misses them). | multiscale plan, step S3.5: valleys, troughs and cirques by their width on the ground per level |
| `erosion-0079` | 2026-10-10 | tributaries available, off | regional.tributary_enabled stays 0 until S3.7; the inherited bed min(parent link bed, running minimum) kept. | multiscale plan, step S3.3: tributaries grown per level from the level's channels (off by default) |
| `erosion-0080` | 2026-10-10 | level-3 overlaps restored; seamlessness by test, not by bound | The solve's margin of 2 global cells kept at every level (by trial); a design question for the solve's catchment at the ring before S4.1. | multiscale plan, step S3.8: the level-3 overlap defect (the halo's global-cell margin for the solve), the forms' reach counted, a ring-basin diagnostic |
| `erosion-0081` | 2026-10-10 | bounded solve available, off; seamless by test | regional.catchment_reach_global_cells 0 (today's solve) for the judging runs; S4.0 added for a connected local labelling before the tiles; hashed tie keys kept. | multiscale plan, step S3.9: a bounded-catchment solve (catchments labelled by local floods, capped at the outlet's global cell), off by default |
| `erosion-0082` | 2026-10-10 | lattice removed | jitter_edge_margin 0, jitter_junction 0.5, diagonal_transit 1 with tolerance 100 m; the remaining peaks on A and C recorded, not chased. | multiscale plan, step S3.10: channels without a lattice (diagonal links through a qualifying transit block, anchors spread over a full cell period) |
| `erosion-0083` | 2026-10-10 | phase S3 closed | Tributaries off for the judging runs (run 40 for the user's comparison); the judging runs 33 to 40 replace 26 to 32. | multiscale plan, step S3.7: Phase S3 close (the anisotropy check, the striping traced and removed, tributaries decided, judging runs 33 to 40) |
| `erosion-0084` | 2026-10-11 | step done; the user's verdict on run 41 pending | Roughness 0.05 and persistence 0.5 by default until the user's verdict; the island question (S3.13) and the tributaries' default stay open with the user. | multiscale plan, step S3.11 (added from the user's look at run 40's discharge map): confluences of the pinned channel pieces, hashed roughness on the legs, the discharge map masked under lakes, a 1-degree subgrid on the regional graticule |
| `erosion-0085` | 2026-10-11 | step done; box 1 accepted at 1.45 with its cause; the user's verdict on run 42 pending | S4.0's catchment fill keeps the flooded cells, labels and corner rule, the low edges at the flat's filled level and the order rebuilt from the final receivers. | multiscale plan, step S3.12 (added from the user's second look, 2026-10-11): the regional solve's pits and flats resolved as the global routing's (Barnes et al. 2014) with metric distances and hashed ties |
| `erosion-0086` | 2026-10-11 | step done | keep | multiscale plan, step S3.14 (added from the user's verdict of 2026-10-11): the painted channel network closed (no dry tails); the leg roughness 0.08 with a relief term |
| `erosion-0087` | 2026-10-11 | step done; judging pending | keep | multiscale plan, step S3.16 (added from the user's noise question, docs/architecture/noise-functions.md): the OpenSimplex2S 3D kernel behind landscape.noise.kernel, per-octave rotation behind landscape.noise.rotate_octaves, both off by default |
| `erosion-0088` | 2026-10-11 | step done; the user's look at the lake's outline pending | keep | multiscale plan, step S3.13: rivers through lakes (the incision floored at the lake surfaces downstream, the wet rule for channel cells in every block, the drychannelinlake and channelbetweenlakes diagnostics) |
| `erosion-0089` | 2026-10-11 | step done with open items; the user's verdict pending | keep | multiscale plan, step S3.15 (added from the user's verdict of 2026-10-11): the coast's shoreline displaced by a fractal position warp under the level band rule; sea pockets not connected to the core sea filled with their basins; the straight coast at 58 E and the blobby coast |
| `erosion-0090` | 2026-10-11 | fix done; the user's look at runs 49 and 50 pending | keep | multiscale plan, step S3.15's rounds four and five and the lake flood's at-S rule (the orchestrator, from the judging run 47): a filled sea pocket beside a lake belongs to the lake below its surface; the basin raise never touches a global lake cell; a cell at exactly S inside a lake cell's own block is the lake's |
| `erosion-0091` | 2026-10-11 | closed; the user's look at runs 51 and 52 pending | keep | multiscale plan, S3.15 closed: the lake-pocket rule of erosion-0090 reverted (a filled sea pocket beside a lake is land at the floor like every pocket); the basin and pocket exclusion of lake cells and the flood's at-S rule kept |
| `erosion-0092` | 2026-10-11 | done | keep | multiscale plan, the user's notes of 2026-10-11 (evening) on the S3.16 pages: the maps' subgrid at 1 degree under whole-degree graticules; window A read against the Himalaya (a preferred Earth window per window block in the regional profile) |
| `erosion-0093` | 2026-10-11 | done; S5.3 judging runs next (A at k 64 on run 10 beside the Himalaya; the belt beside the Zagros and the Subandean belt; phasor and Gabor) | keep | multiscale plan, S5.2: the fold fabric (regional_folds.{hpp,cpp}): phasor-normalised Gabor stripes from kernels hashed per absolute tile, the wavelength from a 10 km base plus the sediment cover clamped to 5-20 km, the amplitude from the strain and the cell's height above the sea, drawn at the level whose octave band holds the wavelength, channel cells untouched, water gaps within 0.25 L, hook after the solve before the valley floors |
| `erosion-0094` | 2026-10-11 | done | keep | multiscale plan, the kernel default (the user, 2026-10-11 evening: OpenSimplex2S everywhere but the fold fabric): landscape.noise.kernel 1 and rotate_octaves 1 by default; the shoreline warp normalised by its field's own deviation |
| `erosion-0095` | 2026-10-11 | merged as the base of S5.3c; the steps acceptance open | keep | multiscale plan, S5.3b: the fold fabric's second round (Phacelle's jittered splat lattice and compact weight, direction and wavelength interpolated per splat, the relief fraction 0.15, no fabric under the ice or in lake blocks, the lakes flooded on the un-folded surface) |
| `erosion-0096` | 2026-10-12 | merged; the judging round registered next | keep | multiscale plan, S5.3c: the coherent fold phase: the fabric as one phase field (the least-squares potential of the wavevector field on the global grid, oriented by a continuity flood, the axis smoothed on doubled angles, bicubic Hermite at the fine cells, dislocations faded into fold terminations, every fold drawn at the coarsest level with at least 4 cells per wavelength); the kernel sums retired |
| `erosion-0097` | 2026-10-13 | approved by the user: v14 plus the shoreline micro-texture (fine-v16-micro5) is the look; v15 rejected ('much less realistic and aesthetically pleasing than compare-v14') and reverted | keep | terrain sandbox (roadmap T1), the look: sandbox versions v1 to v16 outside the pipeline, 0764d05..8f2a837 (2026-10-12 to 2026-10-13); v15 reverted by 269bba8 |
| `erosion-0098` | 2026-10-13 | awaiting the user's look (roadmap T1.1) | keep | terrain sandbox (roadmap T1.0), the sea floor, first version (522dbdb): step 1 of the user's 2026-10-13 order (sea floor, then sea-level history and lakes, then the 100 m nested test against Copernicus GLO-30, then engine back-arc rollback and S6) |
| `erosion-0099` | 2026-10-14 | reviewed (reviewer on Opus 5.5, Fable off): pass with notes; the user's look at compare-base-v2 pending | keep | terrain roadmap T2.2a: the land branch of terrain.texturing centred at the fractal's measured mean (the missing midlands) |
| `erosion-0100` | 2026-10-16 | the user's look: base-v2, the Earth-shaped history, v18 and v19 approved; the sea floor paused | keep | terrain sandbox, roadmap T1 (T1.1 sea floor to v7d and paused, T1.2 sea level, T1.2b the Earth-shaped history, T1.3 lakes, the v17 diagnosis, lakes and pits, provinces and aggrading basins, the uplands lift) and the world sea-level parameter, 2026-10-13 to 2026-10-16 |
| `erosion-0101` | 2026-10-16 | not reviewed (reviews of evaluations paused since 2026-10-05); not yet seen by the user | keep | terrain roadmap T1.5: the 100 m nested test (sandbox) |
| `erosion-0102` | 2026-10-16 | not reviewed (reviews of evaluations paused since 2026-10-05); awaiting the user's look (pages compare-t16-1 to -5) | keep | terrain roadmap T1.6: sandbox artefacts, diagnosis round |
| `erosion-0103` | 2026-10-16 | not reviewed (reviews of evaluations paused since 2026-10-05); the user has not looked yet | keep | terrain roadmap T1.7 (G's fjords), the window smoothness (T1.6/T1.8b), one nesting mechanism for the window levels |
| `erosion-0104` | 2026-10-17 | the user's visual verdict (the recorded one above); no agent review | keep | terrain roadmap T1.7, T1.8b: the user's verdicts of 2026-10-17 and the decisions taken from them |
| `erosion-0105` | 2026-10-17 | pending the user's look at the pages (asked 2026-10-17) | park | terrain roadmap T1 (step 1 of the 2026-10-17 session): confirming renders of the new defaults, the glacial-ratio pair, the global K question, straight rivers, fill_pits |

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
- `erosion-0024`: A rule that decides which cells are sea must be shared by every consumer (routing, erosion, rivers, basins, checks): a second copy of the old rule in one consumer silently undoes the change there.
- `erosion-0025`: A stage that fills depressions before a water-balance lake stage leaves it nothing to fill: when a new process takes over a feature, look for the older stages that pre-empt it.
- `erosion-0028`: The fjord stage's name promised a carving it never did: read a stage before planning to feed it.
- `erosion-0029`: Stages do not share a random stream across processes (seedprocessstream): removing a stochastic stage changes only what it wrote.
- `erosion-0030`: The flotation floor belongs inside the implicit sweep: a clamp after the solve would have left the erosion, deposition and sediment flux inconsistent with the elevation.
- `erosion-0031`: Feeding a state (the ice thickness) back into the routing one step late makes it flicker; Hergarten's bed directions avoid it.
- `erosion-0031`: A temperature field cooled only above the sea level gives an ELA on max(elevation, sea level), not on the seabed.
- `erosion-0032`: A comparison with a run from before new fields must still report the run's side; the validation dispatcher skipped the whole kind.
- `erosion-0033`: The regional code was GUI-only by place and linkage, not by dependency: no SFML or ImGui in 25,000 lines.
- `erosion-0034`: A halo default must be checked against the sum of the reaches it covers: the first draft set 4 global cells while the block-mean passes alone reach 6.
- `erosion-0034`: Splitting a header between two parallel workers: implement the shared part (here the patch metric) in the contracts step.
- `erosion-0035`: Moving floating-point code into a helper keeps results bit-identical only if the products keep their order: a*(b*c) is not (a*b)*c, so the helper exposes the pieces, not their product.
- `erosion-0035`: An overlap test must make the two windows differ in every origin the code uses (columns, rows, the date line); two windows clipped at the same pole hide row-offset bugs.
- `erosion-0036`: A seam argument must include every rule that moves information sideways: 8-connected pinned channels let a diagonal receiver cut between two pinned cells, which the plan's 'catchments between channels' assumption did not see.
- `erosion-0036`: An early stop of an iteration is a global decision; for window-independent results the pass count must be fixed.
- `erosion-0037`: A rasterised network's guarantees (no capture, no rise, no loop) need a separation bound proved from the geometry; tests over thousands of seeds did not find the capture the reviewer derived.
- `erosion-0038`: A synthetic world of 32 x 16 cells clips every halo at both poles: an overlap test then never varies the row origin. Use a world tall enough for domains clear of the poles.
- `erosion-0038`: A steep synthetic world hides the risk of plains; the seam test needs a low-relief variant where noise and channel gradients are comparable.
- `erosion-0039`: A field id from the catalogue must keep the catalogue's meaning, not only its dtype and units: river_direction is set on every land cell of a world.
- `erosion-0040`: The block-mean error on real terrain (178 m with two passes at k = 16) was seven times the synthetic world's: calibrate a correction on the real bundle before setting its default.
- `erosion-0040`: Evaluations are never overwritten: a map fix after a registration goes into a second evaluation folder of the run.
- `erosion-0041`: A contracts step that lets two workers run in parallel must name the test files too: the shared test world went into its own header so each step can add overlap variants in its own file.
- `erosion-0041`: Inert parameters do not test a seam: the overlap tests build the struct defaults, so a step whose feature is off by default must add its own overlap variants with the feature on.
- `erosion-0042`: A synthetic test world can hide a real-data case: it set the lake surface equal to the filled elevation, while erosion run 9 has closed lakes under their sills, which left rivers dry across them.
- `erosion-0042`: A cell pinned at a global value in a surface that follows the fine terrain stands as a pillar: the junction's bed needs its own rule (R3).
- `erosion-0043`: A noise that sets a mask belongs on the decision surface, not on the elevation the channels and the solve read.
- `erosion-0043`: A per-cell floor in a faithful block-mean scheme must be capped by what the block's mean allows.
- `erosion-0043`: A synthetic world narrower than the box never exercises the box's edge: a seam test needs a world wider than the domain plus the margin.
- `erosion-0043`: A loop that stops when a whole-window condition holds is seam-safe only if every pass after that point would change nothing; check that before trusting the stop.
- `erosion-0044`: A direction histogram mostly measures which way the window's rivers flow; the share of steps that keep their direction, the longest straight run and the sinuosity measure the eight-direction bias the user saw.
- `erosion-0044`: A step built on an earlier version of a parallel step must be re-measured and re-reviewed on the merged tree.
- `erosion-0045`: A mask drawn by an indicator must see the terrain it floods, or the bed and the shore rules reshape the solve's result to fit it.
- `erosion-0045`: A choice between neighbours (which lake a fine cell takes) must be made per fine cell, or adjacent cells of two global cells disagree.
- `erosion-0046`: When the user pauses, stop every background agent (TaskStop), not only the wakeup loop.
- `erosion-0046`: A completion notification can arrive before the agent's report or after a run cut short: check the worktree for an uncommitted diff before waiting.
- `erosion-0047`: A scratch script named after a standard-library module (glob.py) breaks every later Python run from that directory: run scratch scripts with python -I and give probes distinctive names.
- `erosion-0048`: A contracts step can let three workers run at once when every hook has its own translation unit and test file and the composing file only calls them; a hook that must be a no-op keeps the original expression textually (a re-associated sum is identical only by luck), and the parity check on every judging window is the proof.
- `erosion-0048`: A contract that reads the mask beyond a block states its reach against the margin formula of the land-share rule, and a contract that incises or lifts by a rate times a time states the time bound within which its physical argument holds.
- `erosion-0048`: regional_overlap.py compares bundles made with the same parameters only: after a step adds parameters, a reproduction of an earlier run is compared field by field, and the comparison's log is saved.
- `erosion-0049`: A rate times a time must be bounded by what saturates it (here the local relief) and its note must state the time within which the physical argument holds; otherwise a parameter sweep finds a value that lowers whole blocks without bound.
- `erosion-0049`: Widths of fine-scale forms are set in kilometres, never in fine cells, so a window at another scale shows the same ground.
- `erosion-0049`: A step that changes the land's mean states the acceptable block-mean shift before its sweep, and the shift is attributed by turning each term off in turn (the solve's lowering, not the sides, carried it here).
- `erosion-0050`: A hook that reads a neighbouring block adds a whole block to a share pass's reach: derive the margin from what the final mask reads, not from the hook's own radius.
- `erosion-0050`: A selection rule of 'fewest X' cannot tune a term whose purpose is to add X somewhere; the step must name the place where the term should act and measure it there.
- `erosion-0050`: Save attribution probes as scripts or diagnostics, never as transient prints.
- `erosion-0051`: A smooth gate needs its lift to end where the gate ends, or the coast stops where the lift's ramp meets the depth.
- `erosion-0051`: A bicubic indicator's sign is narrower than the hard rule it replaces (its outer weights are negative): a reach around centres, not a sign, replaces a neighbourhood rule.
- `erosion-0051`: A test that asserts a premise about the box's size must compute the box with the composing file's formula, or a later margin change makes it check less than it says.
- `erosion-0052`: A phase whose contracts step gives each step its own translation unit and test file runs its steps at once; the reviews, not the implementation, set the pace (three reviews of 8 to 10 min each, three workers of 24 to 34 min, in about 55 min of wall clock).
- `erosion-0052`: A registered run's code.dirty is false only for the first run registered from a clean tree; the registry write dirties the tree for the rest (as in R5).
- `erosion-0053`: step_check.py writes to out/scratch/phase5 unless --directory names the reference's directory: a reference label in another directory needs the option, or the compare fails on a missing manifest.
- `erosion-0053`: The glacial past's ELA and the present climate's differ by cooling_c over the world's lapse rate (881 m on erosion run 9): the forms cut by ice read the run's ELA, the present-day ice the present one, and a bundle-recorded parameter the regional pass depends on is restored by the region command as the stream power law is.
- `erosion-0053`: A probe that picks a manifest key by its suffix (lapse_rate) can take a climate constant for the world's scalar field: read the field (atmosphere/lapse_rate.npy) the code reads.
- `erosion-0053`: A form centred by a rule that reads the solved elevation at a non-channel cell is not seam-safe at the smallest halo, even when its footprint is within k: only the start surface, the sea and the channel beds agree in the innermost halo ring; a placement or eligibility rule reads a pure global quantity (the upsampled base) there.
- `erosion-0054`: A width law with a cap must be swept against the cap: on an ice-sheet world the catchments of the big glaciated rivers put every setting above a few kilometres at the cap, so the sweep's variable is the cap or the exponent, not the reference width.
- `erosion-0054`: A cross-section metric defined against the lower bank reads a widened trough as a shallower valley (the banks came down): report the lowered zone's width and depth as well.
- `erosion-0055`: A hand-built test surface must set every field the rule reads: when a rule moves from the fine surface to the global base, the tests' flat surfaces must be flat in both, and an exact edge of a band is not asserted on an upsampled constant (bicubic rounding), only values inside and outside it.
- `erosion-0056`: A threshold rule's sweep variable (the fraction) does not move the iced area: the area is the ELA's and the edge noise's; sweep the edge noise when the area is the question.
- `erosion-0057`: When the implementing model's session limit stops parallel workers whose code and tests are complete, their worktree diffs apply cleanly to the main tree when the contracts gave each its own files: the lead finishes the measurements and commits the steps by their files, and the reviews run on those commits.
- `erosion-0057`: A sweep variable that saturates at a cap on the judging world (the ice catchments here) is reported with the measured distribution, not inferred from the sweep's shape; the reviews caught an inferred number twice this phase.
- `erosion-0057`: The registered maps cannot show forms a few hundred metres deep on a 10 km elevation palette: a relief composite (hillshade or on-minus-off) belongs in the evaluation tooling before the next form-cutting step.
- `erosion-0058`: rock_uplift_rate reaches the world outside the tectonic contract's Snapshot and its snapshot files; a new engine field has two possible routes and the choice is the user's contract decision.
- `erosion-0059`: The existing block-mean correction never sees the solve: a correction measured on the start surface cannot show the regional pass's real block-mean drop.
- `erosion-0060`: Check that a renderer feature a plan names exists before planning a step on it (the hillshade).
- `erosion-0061`: The plan's own worked formula (section 2's border step) was wrong by a factor of 20 at factor 4; a contracts review that derives the plan's numbers from the operator catches what a parity check cannot.
- `erosion-0061`: A reviewer cites only logged evidence: save the build and ctest output to files the reviewer can read (out/scratch/multiscale/).
- `erosion-0062`: A plan step that assumes a worktree for a step needing refs/ cannot run there: Python steps that read refs/ run in the main tree on disjoint files.
- `erosion-0063`: A contract change sent to a running worker (the design answer's halo) landed cleanly because the worker had not yet written the tests that assert the old values; send contract changes as soon as they are decided.
- `erosion-0064`: An isolated worktree agent cannot edit the main tree: review fixes come back in its worktree and are copied over (here safe because the main tree's files were the worker's own earlier versions).
- `erosion-0064`: Numbers printed by passing tests reach a log only with ctest -V: ask for the verbose log when an entry will cite them.
- `erosion-0065`: ETOPO 2022's surface carries the bathymetry of some large lakes (the Caspian, Superior) and only the surface of others: below 0 m is not sea; an Earth window's sea is below 0 m and open to the window's edge.
- `erosion-0066`: Read the rendered image before accepting a caption's claim: the Caspian stayed sea after a fix whose caption said lake beds were land.
- `erosion-0067`: A map's note longer than the frame is cut without warning (map_rendering.frame does not wrap): keep legends short or wrap them; left as an open item for the evaluation scripts.
- `erosion-0068`: Measure a new operator on the judging world before trusting synthetic tests: the first version's weights passed every test on the synthetic world and dug 7.5 km holes on run 9.
- `erosion-0068`: A rule converting a parameter to a level must check its unit first: the valley reach is in level cells, not on the ground.
- `erosion-0069`: A bundle on a new grid needs every Python reader checked in the same phase: the window block's scale and fine offsets stayed the source level's while the fields moved to the pixel grid.
- `erosion-0071`: Exactness checks pass on a map that prints its grid: read the relief image of a judging run before calling a phase closed.
- `erosion-0072`: A smoothing that changes a parent's values changes its children's masks at the next level: a class-respecting operator still needs masks inherited from the parent.
- `erosion-0073`: A projection's feasibility argument assumes the parent is on the right side of the sea level; a land parent below the sea breaks it at every level below.
- `erosion-0074`: A contracts step must reserve every field two parallel workers will read or fill (here the relief amplitude, read by S3.4 and filled by S3.6), or the shared struct becomes their conflict.
- `erosion-0075`: Locate a defect's cell and read its parent before designing a fix for it: the S2.8 pit's parent was a channel cell, which named the cause.
- `erosion-0076`: A fix that is right in principle (the bed double count) can be wrong in effect on the judging world: compare variants on the measured targets before keeping one.
- `erosion-0077`: A rule added to a step must carry the step's existing guards (here the start-surface test that keeps the projection feasible): the review found the new branch skipping it.
- `erosion-0078`: A distance in index space at the north-south spacing is not a distance on the ground away from the equator: convert east-west offsets by the row's spacing.
- `erosion-0079`: An overlap test at a new depth (3 levels) found a defect the 2-level tests could not: every guarantee needs a test at the deepest level the runs use.
- `erosion-0080`: Measure the reach a step actually has (here the solve's sub-basins to the domain ring) before writing it into a halo formula: the declared 2 global cells was not a bound.
- `erosion-0081`: A design's 'by construction' claim needs its reviewer to look for a counterexample before the code is written: the race was visible in the design text.
- `erosion-0082`: A new geometric rule needs the same elevation guards as the rule it replaces (here the transit block's filled elevation between the link's ends).
- `erosion-0083`: A user's visual verdict ('a SW-NE striping') named a real artefact whose cause was three steps away (the channel lattice): an azimuth-independent measure turned the look into a number that a stage trace could follow.
- `erosion-0084`: A rule stated as 'iterate until nothing moves' needs a locality argument before a worker builds it: the fixed point can fail to exist, and when it exists its reach is the dependency depth, not one cell.
- `erosion-0084`: The seam bound of a unit must be re-derived whenever a decision reads another piece's bed: the start surface enters through the beds one global cell further than the geometry suggests.
- `erosion-0085`: A breadth-first sweep in a fixed neighbour order prints its spokes on every flat; the global routing had already learned this (resolveflatcells) and the regional solve had not inherited it.
- `erosion-0085`: Barnes' 2:1 weights and an 8-neighbour Dijkstra still favour diagonals in deep closed basins: a depression is a lake question, not a routing one.
- `erosion-0086`: Attribute a visible gap by tracing the chain that should be there, not by the nearest counter: the dry tail was real but elsewhere.
- `erosion-0086`: A rule that reads one more piece per hop than the rule it extends changes the seam bound; measure it in the step and record it as measured, not derived.
- `erosion-0087`: A worker without network access cannot read a reference: give it the constants in the brief, or read them yourself and check its port afterwards (done here: the normaliser to 1e-16, then the file).
- `erosion-0087`: A generated-parameter key cannot be a string in this repository; a plan that names string values is wrong before the worker starts.
- `erosion-0088`: A design that writes b(lake cell) = S for a lake a kilometre deep has not looked at the lake's bed: check the field's range in the fixture before fixing a law.
- `erosion-0088`: A rule that is right at the finest level (a wet cell is at the water's bed) can pin a coarse cell to a wrong mean: every per-cell rule must say what it means for a 39 km cell, and water on a channel cell is a mask over a mean, not a surface.
- `erosion-0088`: Six rounds on one step: when the second round contradicts the design's law, stop and re-derive the law from the fixture's numbers before the third.
- `erosion-0089`: A fill at a floor that the mean-keeping passes cannot lower moves a block's mean for good: every mask rule that raises cells needs a lowering partner or a block-mean check.
- `erosion-0089`: A displacement under an inherited mask is only as free as the parent's share rule; the levels' consistency and the shoreline's freedom pull against each other, and the plan must say which wins before the worker starts.
- `erosion-0089`: A worker's probe and an earlier phase's probe measure 'the dimension' differently (1.156 against 1.02 on the same bundle): record the probe with the number.
- `erosion-0090`: A worktree without the sibling steps' changes cannot judge their interaction: the judging run on the merged tree is part of every step's acceptance, with a look at the maps.
- `erosion-0090`: Three separate rules (the basin raise, the flood's strict inequality, the pocket fill's floor) each left the lake's edge on a block boundary: a mask rule that touches a lake's neighbourhood must say what it does to the lake.
- `erosion-0091`: When two rounds in a row move an artefact without removing it, the cause is in another unit: trace the levels' inheritance before a third rule.
- `erosion-0092`: A subgrid must divide the graticule step, not merely be finer: the next finer catalogue step (2 under 5) prints lines the labels cannot place.
- `erosion-0093`: A design's calibration from one field (the sediment) needs the field's distribution in the target windows before the step: run 10's belts carry almost no sediment.
- `erosion-0093`: A relief map of the scratch run before the review catches what the diagnostics cannot name (the corduroy, the patches).
- `erosion-0093`: A reviewer's causal explanation is a hypothesis: the round that tests it costs 30 min; the test was right to run, the explanation was wrong.
- `erosion-0094`: A per-kernel constant (the one-octave deviation) must be read from the field, never from the kernel's header: the shoreline warp had hard-coded the Perlin one.
- `erosion-0095`: A uniform-field unit test with an along-strike step bound would have caught the kernel sum's dislocations at S5.2; the overlap and bearing tests did not.
- `erosion-0095`: The brief's acceptance figures must come from the same tree as the runs: the belt's 35 forced parents were a Perlin figure, the current tree's no-fold run has 93.
- `erosion-0096`: A fabric at 15 wavelengths per global cell needs a direction field smoother than the engine's per-cell tensor: the field's resolution, not the kernel, sets the coherence.
- `erosion-0096`: Keep the before/after figures on one tree: the no-fold baselines moved with the kernel default.
- `erosion-0097`: A window's outer ring of sea at -1 m grades every river leaving the window to the sea level (canyons, the user's 'rivers jumping up'): the ring is pinned at the lowest terrain within half a world cell.
- `erosion-0097`: A latitude equilibrium line reaching the sea at 62 degrees put an ice sheet on every polar lowland for 20 Myr and cut it below the sea in patches the sea re-flood drowned: the ELA needs a floor and a precipitation term.
- `erosion-0097`: With the delta-top cap lifted, the marine fill-and-spill went quadratic (8 min on B, 21 min on the world against 35 to 40 s) and flooded troughs stacked marine loads on single cells (towers up to 12 km): sea cells are capped at the sea level after each step.
- `erosion-0097`: A slope limiter's averaged cut must stop at the lowest capped height around a cell (it dug closed pits to -15 km beside glacier walls); a cone per channel cell ribs gorge walls (averaged over 3 x 3).
- `erosion-0097`: Switches at the base's zero contour (the relief's sea fraction, the uplands) draw the old coastline as a crease: fade them over a height range.
- `erosion-0097`: Where the warp stretches most, a noise drawn on its coordinates varies only every 150 km along the trend: map the stretch before adding texture.
- `erosion-0097`: Runevision's filter adds only a faint texture at k 64; it belongs to the finer levels, after the physics or before the solve, never as the base.
- `erosion-0098`: The world's crust age is per 156 km cell, so its minimum near a ridge is a broad trough that never reaches 0: the crest needs the age relative to the lowest age nearby.
- `erosion-0098`: Most of this world's sea is flooded continental crust and most oceanic crust is 300 to 600 Myr old: an age-depth law shapes only a third of the sea.
- `erosion-0099`: Measure the land after every stage before tuning physics: one uncentred noise term outweighed every erosion lever discussed for the flatness.
- `erosion-0099`: The hillslope diffusion reaches about sqrt(kappa t) = 4.5 km over the sandbox's run: it cannot remove a band of the hypsometry.
- `erosion-0100`: Measure before tuning: one uncentred noise term (the midlands), one missing sediment path (lakes) and one unfilled pit field (drowning) outweighed every tuning lever discussed.
- `erosion-0100`: Record the per-cell totals: erosion, glacial share and fill fields answered the drowning question directly.
- `erosion-0101`: Measure a level against Earth at its own cell size before nesting below it: the parent at 611 m was already half Earth's slope.
- `erosion-0101`: An fbm with ridging and per-cell mean removal realises about 0.3 of its nominal amplitude: normalise to the measured rms.
- `erosion-0102`: Hillslope diffusivity is a grid parameter: the same km^2/Myr at 19.5 km and 2.44 km cells acts on different shares of a cell.
- `erosion-0102`: A sloped fill with too coarse a noise turns filled basins into cones with straight fall lines.
- `erosion-0103`: Check which lever actually binds: the exe's 5 % diffusion limiter, not kappa, set the window smoothness.
- `erosion-0103`: A physical exception (no trap under ice) beats removing a process: the trap still protects lakes elsewhere.
- `erosion-0103`: Exactness to a coarse parent must be smooth across parent cells, or the parent grid becomes terrain.
- `erosion-0104`: The user judges mechanisms by their artefacts: a grid that shows anywhere fails the whole page.
- `erosion-0104`: Bold physical changes (open troughs, less diffusion) were accepted; numerical exactness that costs a visible grid was not.
- `erosion-0105`: Measure a ratio's cost on the domain it hurt (B above 4 km), not on the land mean: the mean hides an 11 % local loss.

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
| `tectonics-0008` | 2026-10-07 | pending user review | keep | erosion plan, step D3: the Phase D close (tectonics run 78 with implicit engine rivers, erosion run 7 with D0 and D1) and the close documents |
| `tectonics-0009` | 2026-10-07 | changed as intended | keep | erosion plan, step H1: implicit engine rivers become the tectonics baseline (the user's decision after reviewing run 78) |
| `tectonics-0010` | 2026-10-07 | changed as intended; measurement only in the engine | keep | erosion plan, step H2: the engine's rock uplift rate in the plate's frame |
| `tectonics-0011` | 2026-10-10 | fields exported | Route (a) as approved; cumulative (no decay); the tensor turns with the local spin; new ocean floor starts at zero strain. | multiscale plan, step S5.1: the engine accumulates a strain tensor per cell in the plate's frame and exports the cumulative shortening and its direction into the world fields geosphere.shortening and geosphere.shortening_direction (the user approved route (a) on 2026-10-10) |
| `tectonics-0012` | 2026-10-16 | reviewed (reviewer on Opus 5.5, Fable off): pass with notes; the report script's mass check fixed at commit | keep | terrain roadmap T4.12, first round: continental growth as a world parameter (keys off by default) |
| `tectonics-0013` | 2026-10-16 | reviewed (two reviewers on Opus 5.5, Fable off): db9f726 pass with notes, c7cd42a pass with notes (the saturation pulse fixed in d5556a1); the evaluation itself not reviewed (reviews of evaluations paused since 2026-10-05) | keep | terrain roadmap T4.12, second round: seeding fix, column continental_area field, upper thickness regulator (keys off by default) |
| `tectonics-0014` | 2026-10-16 | reviewed (Fable 5.1): 0c64d84 pass with notes (fixed in 5b3f384), 5b3f384 pass with notes (fixed in ae2c114), 0e725bf pass with notes (fixed in 075383e); the evaluation itself not reviewed | keep | terrain roadmap T4.12, third round: continental freeboard (lithospheric root, its cap), area budgets, continent-only controller, resolution check, the sea placed after the coast cleanup |
| `tectonics-0015` | 2026-10-17 | agent measurement against the round's criteria (fails area, median and the 3000+ band); no user look (no candidate) | keep | terrain roadmap T4.12, fourth round: orogenic collapse above an onset thickness (the missing high mountains) |
| `tectonics-0016` | 2026-10-17 | agent measurement against the round's gate (passed by W3 at 400 Myr); the user's look pending | keep | terrain roadmap T4.12, fifth round: collapse into margins and creep below the onset (the area source kept, the interiors thickening) |

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
- `tectonics-0008`: `uw.py <module> run` without arguments uses the profile's registered settings, which for erosion are the minimal configuration: a judging run repeats the previous judging run's command from its run.json.
- `tectonics-0008`: A legacy stage whose writes no measured run ever saw can start writing when an upstream change creates the cells it acts on (here new land near the poles): the close's tracked run at the judging size finds it, a minimal run never would.
- `tectonics-0010`: Sampling cell averages bilinearly is not the same transport as a conservative remap with limited parabolic reconstruction: a 'Lagrangian' measurement must move its fields with the very scheme that moves the material, or curvature turns into spurious uplift.
- `tectonics-0011`: A carried quantity without decay keeps every error forever: test the source of new material (here new ocean floor) explicitly.
- `tectonics-0012`: Name the measure in every table: the same run gave 0.28, 0.336 and 0.545 continental area by three measures.
- `tectonics-0012`: A controller on mass does not hold area when thickness is free.
- `tectonics-0013`: Measure the engine at the resolution the pipeline uses: 128 passed, 256 failed.
- `tectonics-0013`: Elevation per km of crust decides the land height: Earth-like thickness with the engine's Airy slope gives land twice too high.
- `tectonics-0013`: Land fraction placement needs area between 0 and -1 km; cliff margins make it lower the sea by kilometres.
- `tectonics-0014`: Elevation per km of crust is the lever on land height; the root key moves it without touching the crust budget.
- `tectonics-0014`: Measure at 256 and check 128 against it with a script: the core statistics agree, the shelf share does not.
- `tectonics-0014`: The shelf band is transient: stretched margins are oceanised by the lower cap faster than arcs make area; a renewal mechanism is still missing.
- `tectonics-0015`: Before switching off a process, find which budget it feeds: the collapse was the area source, not only a thickness limit.
- `tectonics-0016`: A matched control run with the new exe (W0) turns a guessed cause into a measured one: compare budget terms, not end states.

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
| `world-0087` | 2026-10-10 | fields added | Valley spacing calibrated to Perron et al. 2009's Napa Valley D/K and scaled per cell by erodibility and runoff (orchestrator's decision; the worker had marked it the user's); U = the erosion rate (steady state); procedural-world only; the parameters in terrain.json (S3.0's rules barred landscape.json in group d). | multiscale plan, step S3.6: the global grid carries subgrid statistics (diagnoser terrain.subgrid_statistics writing geosphere.relief_amplitude and geosphere.valley_spacing; the regional adapter passes them to the refinement) |

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
- `world-0087`: A grid-scale closure's coefficients (here the erosion's D and K) are not physical constants: a sub-grid law needs its own calibration.
