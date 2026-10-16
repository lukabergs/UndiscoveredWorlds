# Tectonics engine findings from the gold-deposit work

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

Seed 1, baseline recipe (256 x 128, 600 Myr):

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
