# Tectonics engine findings from the gold-deposit work

> **Record since 2026-10-20.** Its open items moved to `docs/architecture/master-plan.md` (the only running plan; §7 maps the old ids to the new steps). Kept for what was decided and learned.

Found while coupling `geology/mineral_systems` to the integrated engine
(force balance, isostatic columns, `configs/tectonics-baseline.json`, seed 1,
256 x 128, 600 Myr, snapshots every 25 Myr). No engine code was changed. The
mineral model now works around items 1 and 2; the others affect it indirectly.
Reproduce with:

```powershell
out/build/headless-modules/bin/tectonic_pipeline.exe batch-export --seed 1 `
  --output <dir> --no-initial --checkpoint-every-updates 25 --mineral-history on
```

## 1. `crust_class` labels deep ocean floor continental

At 600 Myr, 81.5% of cells are `Continental` and none are `Transitional`. Of
the cells below -4,000 m, 55% are `Continental`, yet their column continental
crust has a median of 0 km (90th percentile 4.7 km) over 6.4 km of oceanic crust.
Above -1,000 m the class agrees with the columns: 90% of those cells carry at
least 21 km of continental crust. The class comes from the plate material map
(`classify_crust_class` and `crust_class_from_material`, `lithosphere.cpp:304-327`),
not from the physical columns that drive isostasy and subduction.

## 2. Convergent segment regimes inherit that error

`choose_regime_hint` (`tectonic_boundary_graph.cpp:109`) labels a convergent
contact a collision when both cells are continental or transitional. Classifying
the same segments by median column continental crust within one cell of the
boundary (at least 20 km on both sides = collision) gives, at 600 Myr:

| Engine regime | Column: collision | Ocean-continent | Ocean-ocean |
|---|---|---|---|
| `continent_collision` (49 segments, 502 cells) | 20 / 298 | 10 / 61 | 19 / 143 |
| `convergent_arc` (26 segments, 486 cells) | 0 | 0 | 26 / 486 |

(segments / boundary cells). Across all checkpoints, `continent_collision`
covers 284-783 boundary cells, against 128-550 cells of column collisions.
Andean-type margins therefore appear as collisions, and `convergent_arc` marks
intra-oceanic margins only. Consumers of `geologic_regime` inherit this.

## 3. Intra-oceanic subduction builds no emergent arcs at this resolution

On the overriding (younger) side of ocean-ocean margins at 600 Myr, cells 0-5
cells from the boundary have a median column continental crust of 0-0.1 km and
median elevations of -4,700 to -5,800 m; 1-6% are emergent (7-11% at 225 Myr).
Arc magma is added to the oceanic column (`tectonic_columns.cpp:471`), which
is a reasonable choice, but in seed 1 island arcs never emerge. Ocean-continent
margins do build emergent, thickened overriding edges (32-49 km of continental
crust, 50-90% emergent one to four cells inboard).

## 4. The arc distance rounds up to whole cells

`steps = ceil(arc_distance_km * 1000 / cell_m)` (`tectonic_margins.cpp:78`)
walks two cells inboard at 156 km cells, so the 200 km arc receiver sits about
390 km from the plate edge (receiver centre). The same parameter gives about
273 km at 512 x 256 and 254 km at 1024 x 512. The arc-trench gap therefore
depends on resolution.

## 5. High terrain accumulates beyond Earth's hypsometry

The share of land above 2 km rises from 0.9% at 25 Myr to 15-18% from 300 Myr
on, and above 3 km to 9-11%. The repo's ETOPO 2022 reference
(`refs/processed/metadata/2048/tectonics-hypsography.json`, ice surface) gives
10.7% and 5.0%. The solver erodes 10-16 Gt/yr in total (3.6-5.7 km3/yr), of
which rivers erode 4.4-8.9 Gt/yr, over roughly twice Earth's land area. The
mineral model's reference denudation (5 m/Myr at 1 km, cubic in elevation)
applied to the same terrain removes 2.2-4.1 times the solver's total eroded
volume after 250 Myr. Either the solver erodes high terrain too slowly or the
proxy is too steep. Exposing per-cell erosion would let the mineral model follow
the solver exactly.

## 6. Public API gap (resolved 2026-10-01)

The API now provides `tectonics_get_column_cells`, `tectonics_capture_snapshot`
and `tectonics_get_map_width/height`. The mineral frame reader and
`tectonic_pipeline` use only the API; `scripts/check-layout.py` rejects internal
engine headers outside the engine, its tests and `tectonic_column_probe`.
Heightmap TIFF I/O moved to `src/io/topography_tiff` (`uw_topography_io`).

## 7. Always-zero cycle counter (resolved 2026-10-01)

Removed from the engine, API, snapshot contract (schema 7), pipeline output,
planet state and saves (version 18; version 17 saves still load).

## 8. Mantle plumes (engine change, opt-in, 2026-10-01)

New module `tectonic_plumes.cpp/.hpp` (`tectonics::plumes`), scenario
`plume_processes` (recipe keys `plume-*`, `tectonic_pipeline --plume-*`). The
default `count` 0 changes nothing: the seed-1 procedural world from
`world_pipeline` is bit-identical with the plume-capable engine, and the CLI check
passes.

- Placement: stateless seed-salted hash (salt `0x5bd1e995`), ordered rejection
  for a minimum spacing on the torus, fixed in the native (mantle) frame. Nothing
  draws from `_randsource`.
- Source over each update `[t, t+dt]`: steady tail flux from the birth time and a
  smoothstep head pulse, both integrated exactly over the interval and spread with
  normalized Gaussians, so volumes do not depend on the time step or resolution.
- `columns::Model::advance` adds the basalt after arc emplacement as zero-age
  `ocean` volume, counted in `generated_ocean_kg` and the new subset
  `generated_plume_kg`; the mass residual formula is unchanged. A swell relaxes
  `thermal_rejuvenation_volume` toward `swell_age_myr` with the ridge-heating
  response, scaled by influence. Weakening scales rock strength and is rejected
  unless progressive rifting is on.
- Outputs: optional `plume_influence` (0..1) in snapshots, bundles and the API
  (`tectonics_get_plume_influence_map`, null when off), palette entry, budget
  field in pipeline JSON and the report data checks.

Seed 1, baseline recipe v3 (256 x 128, 600 Myr):

| Configuration | Plume basalt (kg) | Residual | Land fraction | Deep ocean (< -6 km) |
|---|---|---|---|---|
| No plumes | 0 | -1.1e-14 | 40.0 % | 13.8 % |
| 6 plumes, 0.05 km3/yr tails, 2e6 km3 heads, births over 400 Myr, swell 15 Myr | 3.7e20 | -8.6e-15 | 40.7 % | 13.2 % |
| 6 plumes, 0.25 km3/yr tails (120 km), births at 0, swell 10 Myr | 2.6e21 | -7.1e-14 | 42.9 % | 12.8 % |

The plume mass equals the prescribed tails and heads. With 0.05 km3/yr the tracks
are short and low on 156 km cells; with 0.25 km3/yr curved, age-progressive
chains and volcanic islands appear and follow the plate-motion changes. Plume
basalt changes buoyancy and therefore the force balance, so runs with and without
plumes diverge regionally; compare them qualitatively. Maps:
`out/evaluation/plumes-seed1*/index.html` (scratch, not registered).

Not modelled: mantle-supported (dynamic) uplift beyond the oceanic swell, which
would need a signed elevation term through every `elevation_m` caller; plume
heads do not start rifts by themselves; plume-related deposits; the app's
"Tectonic volcanoes" stage still places hotspots at random rather than reading
engine plume positions.

```powershell
out/build/headless-modules/bin/tectonic_pipeline.exe inspect --final --seed 1 `
  --plume-count 6 --plume-tail-flux-km3-yr 0.25 --plume-tail-radius-km 120 `
  --plume-swell-age-myr 10 --plume-swell-radius-km 300 --bundle-output <dir> --output <dir>/summary.json
```

## Mineral-side changes made

- Convergent margins are classified from column continental crust
  (`TectonicFrame::continentalCrustM`, read through
  `lithosphere::getColumnModel()`). Arcs and accretionary orogenic settings
  occur only on the overriding plate.
- Host crust is any cell at or above -1,000 m; `crust_class` is no longer used
  to select hosts.
- Formation rates, maturity (5-20 Myr; 50-75% of convergent length is older
  than 20 Myr) and the erosion proxy are unchanged.

## Engine changes

- 2026-10-04: `TopographyCodec::internal_to_meters_exact` returns the height of
  `internal_to_meters` without rounding (same branches and clamps; the integer
  function is its nearest whole metre and is unchanged). `world_tectonics.cpp`
  stores it as the world's base elevation, which is float since consolidation
  phase 1; code that reads whole metres sees the same values as before.

## 9. Erosion plan Phase D (engine changes, 2026-10-07)

- `columns::Model::advance` measures each update's tectonic elevation change
  (`elevation_m` after `spreadRoots` minus before `transportColumns`, so the
  engine's diffusion, rivers, waves and their isostatic response are excluded)
  and keeps a 10 Myr exponential mean (`uplift_rate_window_myr`), exported as
  `tectonics_get_rock_uplift_rate_map` and `geosphere.rock_uplift_rate`. The
  measurement is Eulerian: plate motion over relief contributes a steady
  -v.grad z that the mean does not remove, and old ocean contributes its
  thermal subsidence (eroded run 7: ocean mean about -11 m/Myr, land p90
  +29 m/Myr). Measurement only: every engine result is unchanged.
- `river_implicit` (off by default) solves river incision with the landscape
  kernel on the engine's torus (Braun-Willett implicit step, Priority-Flood
  depressions, marine e-folding deposition), K = river_erosion / sqrt(1e-3)
  with a uniform 1 m/yr runoff. Most of its effect is the base level, not
  the solver: the baseline's explicit outlets erode down to the seabed
  (`river_outlet_sea_level` off); at 128 x 64, 50 Myr the river erosion is
  4.39e20 kg explicit, 2.19e20 explicit with sea-level outlets, 2.35e20
  implicit. Closed basins pass their load to their spill (the explicit pits
  kept it). At 256 x 128, 600 Myr (run 78 against 77): land above 2 km 18.0
  to 15.8 percent of land, above 3 km 12.2 to 10.1 (finding 5: Earth 10.7
  and 5.0), land fraction 33.6 to 37.1 percent.
- Continental crust (2026-10-16, terrain roadmap T4.12): no engine process
  creates continental crust (arc melt goes to the oceanic reservoir);
  continental area changes only by rift stretching, collision (area divided
  at constant volume) and the start-up seeding; its volume falls only by
  erosion into sediment. Run 78 starts with a surface-mean continental volume
  of 22.7 km, 1.65 times Earth's (CRUST1.0: continental crust 38.3 percent of
  the surface, mean 36.1 km). `crust_class` is drawn from the plate rasters
  (`updateFieldStackMaps`), not the columns: its 81 percent continental and
  14.6 percent ghost cells are a raster artefact; the columns keep 0.62 to
  0.65 of the surface continental. T4.12 added continental_fraction,
  continental_min_thickness_m, arc_continent_fraction,
  subduction_erosion_ratio and a growth controller (all off by default,
  legacy-identical); with the Earth-like set the area shrinks 0.40 to 0.28
  and the thickness grows to 52 km over 600 Myr: nothing bounds the
  thickness from above (the next round: an orogenic-collapse regulator).
  The force-balance solver moves plates as rigid translations (no rotation);
  the mantle forcing is an analytic pattern with a 600 Myr wobble.
- Continental thickness regulation (2026-10-16, T4.12 second round; 128 x 64,
  seed 1, 600 Myr, bundles `out/scratch/t412/R*`, measured with
  `scripts/benchmarks/continental_growth_metrics.py`; area = cells with
  column `continental_area` >= 0.5, the new field `geosphere.continental_area`;
  thickness = sum continent / sum continental_area). The over-thickening sat in
  mixed cells: remap noise against the one-sided lower cap, collisions at
  constant volume, arc continent on outlets with oceanic_area <= 0.5 and root
  flow all add continent without area. With the seeding fix alone (S4b) the
  column area fell 0.400 to 0.20 and the thickness rose to 61 km while crust
  over 20 km per cell held 0.31 (the two measures diverge: name the measure).
  `continental_max_thickness_m` (an upper area cap, oceanic area relaxing into
  continental over 10 Myr where continent > area x h_max, no rock moved) with
  `orogenic_collapse_km2_myr` (linear diffusion of crust above h_max into any
  neighbour) holds the area: h_max 40 km gives 0.33-0.36 and 39-42 km; 36 km
  with D 20000 km^2/Myr (R14, the baseline since this date) 0.364-0.400 and
  37.1-39.5 km, engine land 0.30-0.33, crust over 50 km 0-4.6 % of the
  continental area; D 30000 erases the thick orogens (0-0.7 %). **The 128-column
  calibration does not hold at 256 columns** (base-v5's tectonics and the
  tectonics-recipe runs `out/scratch/t412/T1..T5`, 600 Myr, about 13 min each):
  area 0.40 to 0.345, thickness 37 to 41 km, engine land median 1.2-1.4 km
  (spherical weights; Earth 420 m), and the surface between -200 and 0 m 0.2-0.4 %
  (Earth 5.0): margins are cliffs, so `world.ocean.land_fraction` 0.292 lowered
  the sea by 2.9 km (`out/scratch/sandbox/base-v5-sea`: 21 % of the surface
  above 3 km). An area cap of 32, 33 or 34 km with collapse only above 40 km
  (a split threshold, tried and reverted) changed nothing: at 256 columns the
  margins are sharp (area 0 or 1), so the cap finds few mixed cells, and a
  plateau of uniform thickness has no gradient for the collapse. What sets the
  high land is the Airy slope (below). A cap that
  saturates must convert the last oceanic area fully (d5556a1): otherwise the
  column's basalt sits under a vanishing area and a rounding overlap subducts it
  in one pulse about 300 Myr later. The continental volume drains into sediment
  (surface mean 15.0 to 14.5 km over 600 Myr; the controller counts sediment as
  reservoir), so at a fixed thickness the area slowly falls. The interior land
  stays high (median of land more than one cell from the sea 535-970 m at 128
  columns): the engine's Airy slope is 0.1515 km of height per km of crust,
  twice CRUST1.0's fitted 0.077, and the seeded land relief (2000 m) starts
  the interiors at 660 m.

## 10. The switch to the cube (W2-E, E3, 2026-10-19)

The engine runs on the cube alone (`TectonicGrid::cube`, N = `sphere_face_cells`).
Found while porting its tests; each is fixed in E3:

- A plate's displacement is its rotation's velocity over the update, w x r R. The
  tangent of the rotation's arc (the log map) is not a rigid field: between two
  cells it differs at second order, which read as opening across a slipping
  boundary (ridge heat from pure slip) and as strain.
- The strain increment reads the displacement less the plate's rigid rotation
  (the rotation's vertical spin is added analytically). A least-squares gradient
  of a rotation at a plate edge, one-sided there, gave rigid plates a strain of
  about w sin(latitude) per cell spacing (4e-3 at a divergent edge of a 100 km
  band, N 8).
- Rift partition: each cell to the nearer anchor by great-circle distance, and
  the fault length as the cut edges projected on the cut's normal. A breadth-first
  split and the staircase length made a cut across the grid's axes (at a cube
  corner) about 1.4 times longer, its tension 61 percent of the face-centre case,
  and an island there never failed (now 46 against 40 Myr).
- Orogenic collapse: the single-threshold flux multiplies by the edge weight in
  the same order as the two regimes, so an onset without margin receivers is the
  onset alone bit for bit (the weights are not 1 on the cube).
- Explicit limits on the cube are tighter by the non-orthogonal correction factor:
  surface diffusion stops above 0.2 / factor, and column motion above 0.95 of a
  cell's smallest edge (about 0.71 of a column at a face edge, so 0.6 column per
  update at most where 1 was allowed on a regular grid).
- The per-plate transport's total-restoring scale (|scale - 1| near 1e-7 per
  update) reads as an uplift rate of about 1e-5 of a carried rate away from the
  rock's path, and moves held cells by rounding; tests compare within that.
- Distributed deformation reads the displacement less its area-weighted least-squares
  rotation (the plates' common rotation), so a rotation of the whole sphere changes its
  result by itself alone; the corrections' vector differences, moved between cells, had
  differed by about 2e-3 of a spacing (E3 review).
- One boundary segment id per cell can drop a junction at a cube corner: where three
  sectors meet at a corner whose three cells each border both other plates, a cell keeps
  only one pair's id, and an arm may go unrepresented within the probe radius
  (`TectonicJunctions.FindsAThreePlateJunctionAtACubeCorner` asserts the corner at
  (1, 1, 1) only; the opposite one is lost).
- E3 made no imported-world run: that recipe starts at the mountain bases and has no
  tectonics stage, and E3 changed no code it runs (Earth-climate checked instead).


## 11. Noise and positions on the sphere (W2-E, E4, 2026-10-19)

The position-placed processes read each cell's unit vector; the lat-lon proxy and
the 4-D simplex torus noise are deleted (`tectonic_sphere_noise.hpp`).

- The torus seeding noise had no isotropic equivalent: a field's longitude period
  (the world width) and latitude period (half of it) map onto one noise circle,
  so its features were twice as fine north-south. The sphere fields take the
  geometric mean, R / (sqrt 2 x scale x noise scale) km per noise unit
  (`initial_noise_wavelength_km` 4505 / frequency); the continent macro field is
  29,600 km, larger than the sphere, as on the torus (one or two blobs per world).
- The torus coordinate factors (`sample_x * 1.15` and the like) wrapped with the
  torus' period and so cut a seam at x = 0 in every field read with one; on the
  sphere they only shorten the wavelength.
- The relief warp was 18 cells and the material's 2.4 cells at any resolution, so
  its length in km halved with each doubling of the width; now 2815 km
  (`initial_noise_warp_km`, the 18 cells at 256 columns) and 2.4 / 18 of it.
- One octave's deviation, measured over 4e6 samples: the 4-D simplex 0.297, the
  ocean relief's value noise 0.455, against OpenSimplex2S 0.360 (Perlin 0.268);
  the sphere fields are scaled to the noise they replace and clamped to [-1, 1].
- Unrotated Perlin octaves (the mantle flow's choice) share one lattice with an
  axis on the pole; with it the smoke run (N 64, 100 Myr) had striping ratios of
  0.34 to 0.76 by band (0-80 degrees) and land 0.299, thickness 37.8 km, continental-land
  median 365 m; with rotated OpenSimplex2S octaves (kept) 0.51 to 1.56 by band
  (A4 0.07 over 0-80 degrees, E3's proxy 0.22) and 0.253, 35.5 km, 152 m. One
  realisation each: the band ratios move by 0.3 to 0.6 between them, so the
  striping ratio of one 256-column smoke run does not separate the kernels.
- Seamounts: one per cell of a coarse cube whose cells are about the spacing wide,
  offset within its middle 60 %, the ellipse in its tangent frame by great-circle
  distance (the torus tiles in 3 x 3 neighbourhoods).

## 12. Acceptance on the cube and the one re-tuning round (W2-E, E6, 2026-10-19)

Seed 1, 400 Myr, `tectonics` recipe to `geology.plate_tectonics` with `--grid cube`;
measures in `out/scratch/sphere/E6/acceptance.json` (by true cell areas; W3 the
torus run at the same columns, its cells unweighted).

- The cube engine thickens continents less than the torus did with the same W3
  keys: at N 64 crust over 50 km covers 0.8 % of the surface (torus 2.8 %), the
  continental-land median is 256 m (W3 409 m) and land above 2000 m 0.34 % (W3
  1.8 %, Earth 3.2 %). Area 0.393 and thickness 37.0 km are within W3's
  tolerances. At N 32 nothing stands above 2000 m (W3 at 128 columns 2.2 %).
- The creep below the collapse onset is the lever: 10000 / 5000 / 2500 / 0
  km^2/Myr give 2000-3000 m 0.22 / 1.03 / 0.94 / 2.05 % and 3000+ m 0.12 / 0.33
  / 1.30 / 2.75 % at N 64. 2500 passes every T4.12 gate at N 64 (area 0.374,
  38.9 km, land 0.279, median 375 m, bands within a factor of two of Earth's
  1.7 / 1.5 %) and is the baseline from version 6. N 32 with 2500 still has too
  few mountains (0.24 / 0.37 %) and too much thin area (0.412 at 36.0 km): the
  engine's thickening depends on the resolution more than the torus' did.
- The end-of-run RMS speed is one sample of a quantity that varies by about
  +-20 % over the run (N 64 baseline: p10-p90 of the per-update area-RMS speed
  1.01-1.45 cm/yr, mean 1.29, last 1.13); thicker orogens (creep 2500) end at 1.35.
- The cost per update follows the transport's substep count (0.95 subcells per
  substep over each plate's dilated reference domain, about 0.24 cell at s 4):
  at N 64 this world needs 2 substeps in 177 of 400 updates against 7 of 50 in
  E3's, so the same code runs 22 % slower side by side with E3's executable.
- The process-variation waves now leave the e1-e2 plane (`variation::multiplier`,
  third wave e1 + 2 e3): with three coplanar wave vectors the pattern was
  mirror-symmetric about the great circle normal to e3. The deformation peak's
  two waves stay coplanar; any two plane waves are mirror-symmetric about the
  plane normal to both, and its product with the multiplier is not.
- The seamount reach is 3 r / 0.75, the support of the longest semi-axis
  (it was 3.75 r; seamounts are off in the baseline).
