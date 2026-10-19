# The spherical world mesh (T8 option (c)): foundation design and frozen interfaces

Written 2026-10-18 night (W0 of `terrain-roadmap.md` §3 T8; Opus 5.5 orchestrator with five surveyors, `W0-S1`..`W0-S5`). The user chose option (c), the spherical mesh through the whole pipeline, and asked for the work to run autonomously: every decision below is the orchestrator's, taken with the recommended option and its reason (§13), for the user to review afterwards. Sections marked **frozen** are the interfaces W2 builds on; after the freeze only the orchestrator changes them, announced to every running worker.

## 1. Principle: one table-driven cell grid, two kinds

Every algorithm that walks cells is written once against `mesh::CellGrid`, a set of tables (centres, areas, neighbours, metric). The tables have two builders:

- **`cube`**: the equiangular cubed sphere (§2). The engine always runs on it; the world runs on it from W3 on (`--grid cube`, the CLI default since W3.2; imported maps are read on lat-lon at the export columns and remapped in, §6).
- **`latlon`**: the cell-centred `W x W/2` grid of today, with tables that reproduce `landscape::SphereMetric` exactly. It is the grid of Earth data, of exports and renders, of climate's process grids, of `.uww`, the desktop and the legacy regional map; and it stays selectable as the world grid (`--grid latlon`, the default until W3.2) for lat-lon parity checks (`step_check.py --grid latlon`) and Earth-raster work; the desktop keeps generating lat-lon worlds until a desktop port.

Why not a cube-only code base: Earth data, renders and the legacy consumers are lat-lon anyway, and a port checked by parity on lat-lon before the grid switch separates "the port is correct" from "the cube changes the world". The lat-lon world grid is not a second simulation: the same code runs on either table.

## 2. The equiangular cubed sphere (frozen)

- Unit-sphere axes: +X at (lat 0, lon 0), +Y at (lat 0, lon 90 E), +Z at the north pole. A point (lat φ, lon λ) is (cos φ cos λ, cos φ sin λ, sin φ).
- Six faces f with centre c, right axis u and up axis v (`uw-cube-1`):

| f | c | u (i grows) | v (j = 0 at the top) | covers |
|---|---|---|---|---|
| 0 | +X | +Y | +Z | lon -45..45 |
| 1 | +Y | -X | +Z | lon 45..135 |
| 2 | -X | -Y | +Z | lon 135..-135 |
| 3 | -Y | +X | +Z | lon -135..-45 |
| 4 | +Z | +Y | -X | north polar cap |
| 5 | -Z | +Y | +X | south polar cap |

  On faces 0-3, u is east and v north; face 4's bottom edge is face 0's top edge and face 5's top edge is face 0's bottom edge, both with the same u (the cross net 4 / 0 1 2 3 / 5).
- Face coordinates a, b in [-π/4, π/4]; the point is normalise(c + tan a u + tan b v). With Δ = π / (2N), cell (f, i, j) spans a in [-π/4 + iΔ, -π/4 + (i+1)Δ] and b in [π/4 - (j+1)Δ, π/4 - jΔ] (row j = 0 is the top, so a face array prints with north up on faces 0-3). Its **centre** is the point at the middle (a, b) (the equiangular centre, not the centroid).
- **Cell id** `mesh::CellId` (uint32) = f N² + j N + i; a field is stored in id order, so as a (6, N, N) array it reads `field[f][j][i]`.
- **Area**, exact: with x = tan a, y = tan b and F(x, y) = atan(x y / sqrt(1 + x² + y²)), the cell's area is R² [F(x1, y1) - F(x0, y1) - F(x1, y0) + F(x0, y0)]. Largest over smallest 1.40 at N 64 (computed, T8 step 1). Amended after W2-W1a: the corner formula loses precision as cells shrink (2e-12 relative at N 128, 5e-10 at 4096), so from N 32 on every cube grid and window (`makecubegrid`, `cubewindowgrid`) uses one rule, `cubekeyunitarea`, a 4 x 4 Gauss-Legendre quadrature of the gnomonic area element (children sum to their parent within 2e-14).
- **Point location** (closed form): the face is the axis of the largest |coordinate| with its sign (ties to the lower face number); then a = atan(p·u / p·c), b = atan(p·v / p·c), i = floor((a + π/4)/Δ), j = floor((π/4 - b)/Δ), clamped to [0, N-1].
- **Resolution mapping**: `--resolution W` keeps its meaning (the lat-lon export's columns, `W x W/2`) and gives N = W / 4 (W must be a multiple of 4; 128 → 32, 256 → 64, 512 → 128). The equatorial cell is then as wide as today's (156 km at 256), and the cube has 0.75 of the lat-lon cells. The engine's own face cells are the scenario key `sphere_face_cells` (default W / 4, the world's N, so the hand-off is an identity on a cube world).
- Radius: `world.body.radius_km` (6371 by default), as the lat-lon metric uses now.

## 3. The neighbour and metric table (frozen)

`mesh::CellGrid` (`src/core/mesh/cell_grid.hpp`) holds, per cell, in id order:

- `centre(c)`: unit vector (double); `latitude_rad(c)`, `longitude_rad(c)`; `east(c)`, `north(c)`: the local unit tangent vectors (at the north pole east is +Y and north -X, at the south pole east +Y and north +X: the limits along lon 0).
- `area_m2(c)`; `totalarea_m2()` (4πR² for both kinds, to rounding).
- **Eight neighbour slots, the D8 codes of today**: 1 N, 2 NE, 3 E, 4 SE, 5 S, 6 SW, 7 W, 8 NW, 0 none (`sphere_metric.hpp:21-27`), in the cell's own frame: geographic on lat-lon; on the cube N is toward v (j - 1), E toward u (i + 1), on every face, also across face edges (the slot names the direction of departure from c).
  - `neighbour(c, k)`: int32 cell id, or -1 where none exists: on lat-lon beyond a pole (where `SphereMetric` returns the cell itself; callers that skip "the cell itself" now skip -1), on the cube the one diagonal that would point at a cube corner (the 24 corner cells have seven neighbours).
  - Across a cube face edge the edge slot goes to the adjacent face's cell that shares the edge; the two diagonals across that edge go to that cell's neighbours along the edge.
  - `reverse(c, k)`: the slot by which `neighbour(c, k)` reaches c. On lat-lon it is `oppositedirection(k)`; on the cube it differs across face edges (the frames rotate), so **no code may compute the reverse as `oppositedirection`**: donor/receiver logic uses `reverse`.
  - `distance_m(c, k)`: on the cube the great-circle distance between the centres; on lat-lon exactly `SphereMetric::distance_m` (dy, the row's dx, or hypot(mean dx, dy)).
  - `facelength_m(c, k)`: on the cube the arc length of the shared edge for slots 1, 3, 5, 7 and 0 for the diagonals (no shared edge); on lat-lon exactly `SphereMetric::facelength_m` (including its mean dx for diagonals). One kernel uses it today (`stream_power.cpp:371`, hillslope flux); a kernel that sums flux over diagonals is reviewed in its port.
  - `bearing_rad(c, k)`: the compass bearing of the slot's direction at c (clockwise from geographic north), for anything with a compass meaning (aspects, "south-facing", winds). On the cube it is the initial great-circle bearing to the neighbour's centre; on lat-lon the planar atan2(dx, dy) of the slot, the diagonal with the mean dx of the two rows (as `distance_m`).
- `kind()`, `facecells()` (cube N) or `columns()`, `rows()` (lat-lon); `cellcount()`; `locate(unit vector) -> CellId`; `cellof(x, y)` and `column(c)`, `row(c)` on lat-lon only (row-major, row 0 north: `SphereMetric::index` today).
- Grids are built once and shared (`std::shared_ptr<const CellGrid>`); `makecubegrid(N, radius_m)`, `makelatlongrid(W, radius_m)`; a **window** is `subgrid(const CellGrid&, global ids)`: a CellGrid of the same kind over a cell subset, local ids in ascending global order, `globalid(local)`, neighbours outside the subset -1 (§8).
- `landscape::SphereMetric` becomes a view of a `CellGrid` with its current per-cell API (`size`, `index`, `area_m2`, `neighbour`, `distance_m`, `facelength_m`); its row API (`rowdx_m`, `rowlatitude_rad`, offsets) stays only on lat-lon and its users are ported (W2 landscape step). `torusmetric` goes with the engine's torus.
- **Geometric helpers** (`src/core/mesh/mesh_geometry.hpp`), shared by every stamp, line and radius loop that today counts pixels: `greatcircle_m(p, q, radius_m)`; `cellswithin(grid, centre, radius_m)` (flood from the located cell, bounded by great-circle distance, deterministic order: ascending id); `cellsalongarc(grid, p, q)` (the great-circle arc stepped at a quarter of the current cell's spacing, located, consecutive duplicates dropped); `tangentframe(centre, radius_m)` and `fromtangent(frame, east_km, north_km)` (gnomonic, so a shape template in km lands on the sphere with its true size); `randomcellbyarea(grid, uniform)` (uniform: a callable returning doubles in [0, 1); a uniform point on the sphere, located: seeds by area).

## 4. Field storage (frozen)

- `worldgrid<T>` (`src/core/world/world_grid.hpp`) stores `cellcount()` values **in the grid's id order** (on lat-lon: row-major, row 0 north, as the kernels, bundles and climate's seasonal vectors already are; today it is x-major). `resize(const mesh::CellGrid&)`; `operator[](CellId)` and `at(CellId)` read a cell; `data()` the flat array. The legacy `grid[x][y]` stays, on lat-lon only, through a strided column proxy (asserts the kind in bounds-checked builds).
- `planet` owns `std::shared_ptr<const mesh::CellGrid> grid_` (`planet::grid()`, `cellcount()`), built by `planet::setgrid(kind, W)`; `setwidth(W - 1)` keeps building a lat-lon grid for its callers. Every per-cell accessor gets a `CellId` overload with the same name (`nomfloat(c)`, `setnomfloat(c, v)`, `basebelowsea(c)`, ...); the `(x, y)` forms forward to `cellof(x, y)` and exist on lat-lon only. Accessors whose bodies read neighbouring x or y stay `(x, y)`-only; the W1b worker lists them in §14.
- **Cube readiness**: `WorldProcessDefinition` (`src/pipeline/world_processes.hpp:37-49`) gains `cubeready` (false for every process at W1); on a cube world `executegenerationstage` refuses a process that is not cube-ready with its id in the message. Each W2 port sets the flag for the processes it ported; W3 starts when every process of the procedural recipe is ready.
- `climateworld`'s grids and seasonal vectors index by cell id the same way; `planet::climateinputsview()` passes id-ordered arrays (its `inputindex` becomes the cell id).
- Records that carry x, y (settlements, polities, routes, crater centres, deposits, segments, regions) carry a `CellId` after their module's port; until then they stay x, y on lat-lon.
- The `width()` / `height()` getters (largest index) stay for lat-lon callers; cube-aware code asks `grid()`.

## 5. Remap and sampling (frozen)

`src/core/mesh/mesh_remap.hpp`:
- `RemapWeights makeremap(const CellGrid& src, const CellGrid& dst, RemapOptions)`: sparse rows per destination cell, entries (source id, overlap area in m²), built deterministically (ascending ids, fixed subdivision) and cacheable.
  - latlon → latlon: exact overlaps (separable in longitude and sin latitude, as `climategrid::makeConservativeRemap`).
  - any pair with a cube: each source cell is subdivided in its native coordinates (face a, b or λ, sin φ) as a quadtree; a piece whose four corners and centre locate in one destination cell is assigned whole with its exact area; others split, to a depth limit `RemapOptions::depth` (default 6: pieces 1/64 of the source edge); a leaf goes to the cell of its centre. Every source cell's area is distributed exactly (column sums exact); a destination row sums to its area within the boundary band of the leaves.
- `remapmean` (intensive: elevations, temperatures, rates): Σ w v / Σ w over the row; a constant field stays exact; a destination row with no weight (lat-lon polar rows at high W finer than the leaves) takes the bilinear sample at its centre (counted in `RemapWeights::emptyrows`). The apply functions take the weights, both grids and the field.
- `remaptotal` (extensive: masses and volumes held by a cell; never accumulated network fluxes such as discharge, drainage area or river flows, which are point values of the drainage network: `nearest` until the export recomputes them with the directions, decision of the W1c review): each source total split by its overlap shares; the global total is exact; a destination row without weight receives zero.
- `remapnearest` (codes, ids, masks): the source cell at the destination centre; `remapmajority` (codes over larger destination cells): the code with the largest overlap, ties to the smaller code.
- `samplebilinear(grid, field, p)` (smooth renders, point sampling, the engine's hand-off of heights): on a cube face bilinear in (a, b) between the four surrounding centres, taking cells across a face edge from the neighbour table; at a cube corner the missing fourth cell's weight is shared by the other three in proportion (exactly their mean at the corner, continuous around it); on lat-lon bilinear with wrapped longitude and clamped polar rows.
- `tangenttoeastnorth(grid, c, vector)`: a 3-D tangent vector to (east, north) components at c; the reverse `eastnorthtotangent`.
- **Direction-coded fields** (`riverdir` and every D8 code) are never remapped: on export they are recomputed on the destination grid by its own routing (§9), and a remap of a direction field is an error.
- **Angles and vectors** (shortening and ridge directions, deformation and plate velocities, winds) are stored relative to the cell's local east and north on every grid kind (never to the grid's axes, which rotate between cube faces), so a remap or export moves them as (east, north) pairs; the engine's hand-off (E5) converts its own (decision of the W1c review).

## 6. The lat-lon export and import of a world (frozen)

`src/world/planet_export.{hpp,cpp}`:
- `exportlatlon(const planet& world, int W) -> planet`: a lat-lon planet of `W x W/2` with every field: floats by `remapmean` (or `remaptotal` where the field catalogue marks a field extensive), codes, masks and ids by `remapnearest`, direction fields recomputed by the lat-lon flow routing on the exported filled elevation; records (settlements, ...) by locating their cells. It feeds `.uww`, the desktop, the legacy regional map, the climate adapter (§10) and renders.
- `importfields(planet& world, const planet& latlon, field ids)`: the reverse for the named fields (climate's outputs, and since W3.2 the imported terrain: on a cube world the CLI reads `--import-land`/`--import-sea` and the imported recipes' Earth maps on the lat-lon grid of the export columns, then moves `base_elevation` (mean) and `imported_lake` (nearest, cleared at sea) onto the cube).
- The field catalogue (`src/core/world/world_fields.cpp`) gains per field a remap rule (`mean`, `total`, `nearest`, `direction`) used by both.

## 7. Bundle v2 (frozen)

v1 (`src/io/run_bundle.cpp:147-210`, W0-S3) stores `grid {columns, rows, registration, layout}`, fields as `.npy` under `fields/<component>/<name>.npy` with shape `[layer,] row, column`, hashes, steps and tables. v2 changes only the grid and adds two optional parts; everything else (manifest keys, field entries, hashes, steps, tables, the catalogue `UW_WORLD_FIELD_LIST`) is as in v1.

- `schema`: `uw-run-bundle/v2`. Readers (C++ and `scripts/uwtools/bundle.py`) accept v1 as `kind latlon` (recorded runs stay readable); writers write v2 only.
- `grid.kind`:
  - `latlon`: `{kind, columns, rows, registration "cell-centred", layout "row-major, north row first", radius_km}`; field shape `[layer,] rows, columns` as v1. No cell table (it is implicit, `uwtools/grid.py`).
  - `cube`: `{kind, face_cells N, orientation "uw-cube-1", layout "face, row (top first), column", radius_km, export_columns W}`; field shape `[layer,] 6, N, N` (id order).
  - `cube-window`: `{kind, face_cells N_L, level L, base_face_cells N, orientation, layout "cell list, ascending id", cell_count n, radius_km, export_columns, export_rows}` (the export box size, added in W2-W1c); export cells the kept cells do not cover are NaN (decision of W2-W1c, no extrapolation); field shape `[layer,] n`; plus the `window` block of today (`src/io/regional_bundle.hpp:62-68`: the box in global pixels or degrees, scale, halo).
  - Regional bundles on a lat-lon world keep their fine lat-lon window grid (`kind latlon` with the `window` block), as now.
- **Cell table** (`cube`, `cube-window`; written by C++, the one implementation of the geometry): `grid/centre.npy` (f8, n x 3), `grid/latitude_deg.npy`, `grid/longitude_deg.npy`, `grid/area_m2.npy` (f8, n), `grid/neighbours.npy` (i4, n x 8, slots 1..8 in columns 0..7, -1 none), `grid/reverse.npy` (i1, n x 8), and for windows `grid/cell_ids.npy` (u8, the level ids). Listed in the manifest under `cells` with hashes. Python never recomputes the topology; `bundle.py` may recompute centres and areas to check them.
- **Lat-lon export** (`cube`, `cube-window`): `latlon/<component>/<name>.npy` at `export_columns` (or the window's box) for every field of layout `grid`, `seasonal`, `monthly`, `levelseasonal`, by the field's remap rule (§6), listed under `exports` with `{id, file, hash, rule}` and `recomputed` where the export was computed on the export grid rather than remapped (directions and accumulated fluxes; decision of the W2-IO review). Fields of layout `zonalprofile` (six values per lat-lon column) live on climate's lat-lon world and are stored in their lat-lon shape in both kinds. Direction fields are exported recomputed (§6).
- `bundle.py`: `Bundle.open` reads v1 and v2; `grid_kind`; `columns`/`rows` are the lat-lon (export) dimensions in every kind; **`load(id)` returns the lat-lon array** (native on lat-lon, the export on a cube), so the ~25 lat-lon consumers (renders, Earth comparisons, `checks.py`, the sandbox's base reading) keep working; `load_native(id)` returns the native array and `cells()` the cell table. `compare_bundles.py` compares native fields (and exports when present).
- C++: `RunBundleWriter`/`Reader` take the grid spec; `writeplanetbundle` writes the planet's kind (and the export for a cube); `loadplanetbundle` builds a planet of the bundle's kind (`--from-bundle` of a cube bundle gives a cube world).

## 8. Windows and the refinement hierarchy on the cube (frozen)

- **Levels**: level L has N_L = N f^L cells per face edge (f the hierarchy's factor, 4 today); the parent of (f, i, j) at L + 1 is (f, i div f, j div f) on the same face: a per-face quadtree with no cross-face parentage. A cell's key at level L is its level-L id `f N_L² + j N_L + i` as uint64 (`mesh::LevelKey`; N_L up to 2^30).
- **Windows** keep the user's syntax: `--pixels x0,y0,x1,y1` in global lat-lon export pixels (W x W/2, rows from the north, x0 > x1 across the date line) or `--window w,s,e,n` in degrees. The box selects the level-L cells whose centres lie inside it (on any faces), plus a halo of `regionalhalo` level-0 cells' width measured in km; the domain is a `subgrid` (§3), so a window across a face edge or a cube corner needs nothing special.
- **Bit-identical overlaps and seamless stitching**: every value of a level-L cell is a pure function of the seed, its key, its centre's unit vector and global fields within the halo (as `regional_pass.hpp:10-13` requires today); two windows that share a cell compute it from the same inputs in the same order. Tests: the overlap suites of today (`RegionalOverlap*`, `RefinementLevelsOverlap*`, `RefinementProjectionBitIdenticalInTwoWindows`) rewritten on the cube, plus a window across a face edge and one around a cube corner.
- A window bundle is a v2 bundle of grid kind `cube-window` (§7) with its cell list; its renders come from its lat-lon export of the box.
- **Window grids (W2-W0 design, 2026-10-19; frozen with this amendment):** `mesh::cubewindowgrid(facecells_L, radius_m, sorted LevelKeys)` builds only a subset's tables in closed form (64-bit `LevelKey` up to N_L 2^30; local ids in ascending key order; neighbours outside -1; a compact mode of about 80 B per cell computing distances on demand bit-identically) with `globalkey(local)` and `localofkey(key)` on `CellGrid`; `mesh::latlonwindowgrid(levelcolumns, radius, domain box)` gives today's lat-lon patches as a CellGrid (today's row-major domain order, no wrap inside a patch, distances equal to `finepatchmetric`) so the regional code has one table-driven implementation with lat-lon parity. A whole-level `subgrid` is not used (6.3 M cells at level 2 for N 64).
- **Halo by rings** (decision): the halo is counted in rings of level-0 cells on the table (8-connected hops), which reproduces today's lat-lon rectangle and matches the reach budget exactly; its km size is reported, not selected (a km radius would have to be rings x the largest diagonal hop and roughly doubles the domain).
- **Cube rules of the regional units** (helpers in `src/core/mesh/window_grid.*`): **G**, a 4 x 4 Catmull-Rom stencil in face coordinates whose nodes beyond one face edge are 1-D Catmull-Rom values along the next face's centre column and beyond two edges a 1-D Catmull-Rom along the third face (amended after W2-W1a: bilinear corner nodes jumped six times more) (decision: ghosts at their true positions, reach 3 rings at face edges, measured in W2-W1a against index-unfolded ghosts); **C**, a cell chart extending a face's index space across one edge by the net rotation; **S**, straight walks (the next slot is the opposite of `reverse`); **R**, rings. Block means are exact over the per-face quadtree; the warp is an (east, north) km displacement placed by `tangentframe`/`fromtangent`; coast and water gates by table neighbours and G footprints; channels hashed by level-0 keys, shared edges named by their key pair, legs within their block, the transit block the common cardinal neighbour (no diagonal link at a corner cell); folds on the level-0 table with axes relative to east/north, Dijkstra weights facelength/distance; global inputs by cell id from the world's grid. Coast rules (decisions after W2-W2c, the worker's recommendations): the warp's level-0 cell size on the cube is π R / (2 N); the smooth gate's footprint is the 2n x 2n chart block around the point without the missing corner quadrant; window boxes widen by dilating their level-0 parents in rings; the land-share offset is bilinear over level 0 with tent weights in face indices; pinch squares use the cell's own E, S, SE neighbours (skipped where SE is missing at a corner); the shoreline displacement is bilinear at the displaced point; the island sink ring is the cells at exactly `margin` rings (missing cells count as land); the straight-run breaker runs along face rows and columns with the side cell from the chart; sea pockets search the level-0 parents within `reach` rings through cardinal neighbours. Rules chosen in W2-W2c-e and accepted: widened boxes keep the lat-lon margin formula in rings; nodes outside a box take the nearest box cell or block by rings; the centre rule takes the child at face offset (k/2, k/2); the straight-run breaker hashes the cell key; forms are discs on the ground by great-circle distance over the level arc π R / (2 N_L); cirque candidates and lake/ice indicators hash the level-0 key; a cardinal channel-edge offset hashes (edge, smaller key, larger key) along the smaller key's block edge; transit candidates are the upstream cell's cardinal neighbours beside the diagonal that are cardinal neighbours of the receiver; links are drawn for level-0 cells within `linkreach` rings. Tributaries have no cube form yet (off by default): the cube driver refuses them. Two further cube rules (decisions after W2-W2b): a stencil node outside the window domain takes the value of the domain cell at the smallest ring distance, ties to the smallest key; prolong's limit is bounded by the 4 nodes of the stencil's central 2 x 2 in face coordinates (true-position ghosts across edges). Lat-lon keeps today's arithmetic in every unit (parity by sha256 of every field against references captured at c481bb1: region-base at level 1, a level-3 refinement, a date-line window, a polar window).

## 9. `.uww`, the legacy regional map and the desktop (compile-only)

They stay lat-lon and are fed by `exportlatlon` (§6) at the world's W: `saveworld` writes the export (the `.uww` layout is unchanged), `loadworld` makes a lat-lon world (a `.uww` has no cube), the legacy `generateregionalmap` and the desktop read the export through the `(x, y)` getters as now. The desktop and `uw_regional` are compiled, never run (`AGENTS.md`).

## 10. Climate

Climate keeps its own lat-lon `SphericalGrid` process grids and its 17,000 lines. Its world-resolution state (`climateworld`) lives on a lat-lon grid of the world's export W. The adapter (W2): before a climate stage, the fields the stage reads (its declared access, `world_process_access.inc`) are exported to climate's lat-lon world (§6); after it, the fields it writes are imported back to the cube (`importfields`). On a lat-lon world the adapter is the identity, so climate results there do not change.

## 11. Noise on the sphere everywhere

- Every noise value is a function of the cell centre's unit vector and the seed: `landscape::FractalNoiseField` (`src/simulations/landscape/spherical_noise.hpp:75`; wavelengths in km, octaves, persistence). The diamond-square `createfractal`, `warp`, `shift`, `flip` (24 + 11 + 4 callers, W0-S5) are replaced at each caller by a field with the wavelength and amplitude the caller's grain stood for at 2048 columns (19.5 km per cell), recorded per caller; `procedural_noise.cpp`'s pixel functions are removed with their last caller (W3).
- Sequence randomness (`random()`, the process streams of `world_processes.cpp:156`) stays; positions are drawn by area (`randomcellbyarea`), never as a uniform row and column.
- The sandbox's Python noise is evaluated at absolute unit vectors for the global and windows alike (T8.8), so overlaps stay bit-identical.

## 12. Module ports and acceptance

Every procedural world changes at the W3 switch: compare totals and counts, not cells. Per module:

| Module (W2 step) | Port | Acceptance |
|---|---|---|
| Engine (W2-E, T8.0-T8.7) | the torus replaced by a cube `CellGrid` (rotations about Euler poles, 3P x 3P torque balance, transport by forward splat on the cube, budgets in true areas, distances in km on the sphere, 3-D noise, seeds by area); the hand-off writes the world's grid through §5 | T8 "Acceptance" (roadmap §3): striping ratio 0.8-1.4 in every band and A4 ≤ 0.08 on the engine's surface; plates against T8.0 and Earth; conservation 1e-6; cube imprint (A4 in face coordinates ≤ 0.028; gradients within two cells of a face edge within 10 % of interiors); 400 Myr at N 64 ≤ 1.5 x 477 s |
| Landscape kernels (W2-L) | `SphereMetric` as a view of `CellGrid`; row-API users and `oppositedirection` donors moved to the table | lat-lon: `step_check.py` parity identical; cube: Tier 0 tests (a cone drains to its apex region, flow accumulation sums to the total area, Priority-Flood leaves no pit, stream power on a face-edge-crossing valley equals the interior case within 1 %) |
| Terrain stages and noise (W2-T1, W2-T2) | (x, y) loops to cell ids; pixel radii to `cellswithin` in km; lines to `cellsalongarc`; templates through `tangentframe`; noise per §11 | lat-lon: totals and counts (land share, hypsometry quantiles, mountain cells, coast length) within the variance of a seed change, recorded; cube: runs at the minimal configuration; no seam (statistics within two cells of face edges within 10 % of interiors) |
| Hydrology (W2-H) | flow_routing, rivers, lake_levels on the table; the planet-getter files (drainage basins, flow queries, wetlands, tides, lakes) to cell ids, hard-coded direction comparisons to `bearing_rad`/`reverse` | lat-lon: parity where mechanical, else totals (river cells, lake count and area, basin count); cube: rivers cross face edges (a test), no polar D8 bias (river direction histogram by band) |
| Society (W2-S) | polities, routes, sites, trade on cell ids; distances in km; A* on the table | `--social`: counts of sites, polities and route length within 10 % on lat-lon; runs on the cube |
| Resources and geology (W2-R) | mineral resources and the mineral frame's areas from `area_m2` (`tectonic_frame.hpp:58,81-82`) | totals of deposits by kind |
| Climate adapter (W2-C) | §10 | lat-lon identical; cube: Earth-data scores (the climate profile's global-balanced-v1 at the minimal configuration) within the seed-change variance of lat-lon |
| io, desktop, legacy regional (W2-IO) | §6, §9, bundle v2 writer/reader (§7) | `.uww` round trip; desktop and `uw_regional` compile |
| Refinement hierarchy and windows (W2-W) | §8 on the table | the overlap tests on the cube, across a face edge and a corner; stitching by `regional_overlap.py` on cell ids |
| Sandbox (W2-SB) | the global on cube fields from the v2 bundle (cell table, neighbour table), noise at unit vectors, `terrain_sandbox.exe` on a cube `CellGrid`, windows as level subsets | the sandbox's overlap test; the T8 striping acceptance on the procedural input and the final |
| Evaluation and Earth references (W2-EV) | `striping_ratio.py` (lat-lon export against ETOPO, and A4 in face coordinates on native fields), `plate_metrics.py` by true area, face-edge seam statistics, profiles read v2 | the scripts run on v1 and v2 bundles; ETOPO's own numbers unchanged |

Timings to record in W3: the full procedural pipeline at the 256 and 512 equivalents (N 64, N 128) against today's lat-lon runs.

## 13. Decisions (the orchestrator's, 2026-10-18 night; for the user's review)

1. **Mesh: equiangular cubed sphere, orientation `uw-cube-1`** (§2). Reason: the T8 design's comparison (four and eight neighbours carry the D8 kernels over, closed-form location, a per-face quadtree for windows); no blocker found by the surveys.
2. **One generic `CellGrid` with a cube and a lat-lon builder; the world grid selectable (`--grid`), lat-lon until W3, cube after.** Reason: one code path for both; lat-lon parity checks each mechanical port before the switch; Earth data, renders and legacy consumers need the lat-lon tables anyway.
3. **Resolution: N = W / 4, `--resolution W` unchanged in meaning (the export width).** Reason: the equatorial cell keeps today's size, so every km-based parameter keeps its meaning; 0.75 of today's cells.
4. **Storage in id order (row-major on lat-lon), `[x][y]` kept only as a lat-lon proxy.** Reason: kernels, bundles and climate's seasonal vectors are already row-major; one layout for both kinds; the legacy consumers compile unchanged.
5. **D8 codes kept as face-local slots with an explicit reverse table and bearings.** Reason: the kernels' eight codes carry over; a stored code is only meaningful with its cell, never as a compass direction.
6. **Remap by exact-area adaptive subdivision; conservative totals exact, constant means exact; directions never remapped.** Reason: exact polygon intersection of great-circle cells with parallels is error-prone; the quadtree is deterministic, simple and accurate to its depth.
7. **Bundle v2: the grid kind, native fields in id order as (6, N, N), a cell table written by C++, a lat-lon export of every gridded field, `load` returning the lat-lon view.** Reason: one format at every scale (memory "maps are maps"); the geometry has one implementation; the lat-lon consumers and Earth comparisons keep working unchanged; v1 stays readable for the recorded evidence.
8. **Windows keep their lat-lon syntax; domains are cell subsets; keys are level ids.** Reason: the user's preference for global pixels (memory "maps are maps"); a subset grid needs no face-edge special cases.
9. **`.uww`, desktop and legacy regional map stay lat-lon, fed by the export, compile-only.** Reason: the T8 decision text; their 24,400 + desktop lines gain nothing from the cube.
10. **Climate keeps its lat-lon process grids; an adapter exports and imports by declared access.** Reason: its physics already handles the sphere; the adapter is the identity on lat-lon.
11. **The W3 re-tuning round (question 4) is allowed: at most one round per module against Earth data.** Reason: the T8 recommendation; true areas move the W3 calibration.

## 14. Status, waves and file ownership

See `terrain-roadmap.md` §3 T8, "W0" and the T8 status table (the single place for state).

**Frozen at 4bdf592 (W1 done, 2026-10-19 00:45).** W1c's departures, part of the frozen text: `exportlatlon` returns `std::unique_ptr<planet>`; export entries carry `dtype` and `shape`; `reverse(c, k)` and `grid/reverse.npy` give slot 0 where there is no neighbour; records in tables are copied as x, y until their module's port relocates them; `cube-window` is implemented by W2-W. `writeplanetbundle` gains the direction callback in W2-IO, and accumulated fluxes (discharge, drainage area, river flows) are recomputed on export by that callback rather than remapped as totals.

### W1b: (x, y)-only accessors

Every other per-cell accessor of `planet` (274 forms) and `climateworld` (164) has its `CellId` overload. These keep only their `(x, y)` form, because their bodies read neighbouring or wrapped x and y, or their y is not a row:
- `planet::outline`, `coast`, `seawrap`, `outlinewrap`, `mountainheightwrap` (`planet.cpp`: the four neighbours, `wrapx`/`clipy`);
- `climateworld::outline`, `seawrap` (the same);
- `planet::horse`/`sethorse`, `climateworld::horse`/`sethorse`: subtropical-high rows, x and one of six bands (`horselats` is 2048 x 6, not a world field on the cell grid);
- not per cell: `planet::longitude(x)`, `latitude(y)`, `reverselatitude`, and the private index helpers `seasonalclimateindex` (both classes).
