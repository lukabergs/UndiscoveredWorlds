# Native tectonic output contract

The implementation is `src/simulations/geology/plate_tectonics` in this working
tree. `tectonic_contract.hpp` defines `tectonics::contract::Snapshot` and its
objects. `capture_snapshot(lithosphere)` copies a state; `tectonic_api.hpp`
exposes borrowed runtime arrays. The app and the `tectonics` recipe of `world_pipeline` use the same
engine and the compiled `configs/tectonics-baseline.json` recipe.

## Geometry and time

Since E3 (W2-E) the engine runs on the equiangular cubed sphere
(`TectonicGrid::cube`, `docs/architecture/sphere-mesh.md`): N = `sphere_face_cells`
cells per face edge, 6 N² cells in uw-cube-1 id order, on a planet of `radius_km`.
Native arrays are that id order; the snapshot and the runtime maps present it as
width N × height 6 N (face f's rows are rows f N to f N + N - 1). Cells have
unequal areas (`area_m2`); the cell spacing is the mean cell side,
sqrt(4 pi R² / 6 N²). Plate motion is a rotation vector per plate; kinematics
and vectors are (east, north) at a cell unless named otherwise, in cell spacings
per update or in physical units where named. `tectonics_get_cell_grid` returns the
engine's `mesh::CellGrid` (shared, so a consumer can cache weights against it).
Array pointers expire when the engine updates or dies.

`time_myr` advances once per physical update by `delta_time_myr` until the finite,
positive `evolution_duration_myr` is reached. A final partial requested timestep
rounds up. Finished updates do nothing.
`iteration_count` and `time_origin_step` encode age timestamps; prefer explicit
elapsed time and `crust_age_myr` for consumers.

## Fields and objects

- `heightmap`: internal topography encoding. Use `TopographyCodec` with the
  snapshot sea-level datum for metres. Initial ocean depth and land relief are
  positive amplitudes relative to sea level, not evolving elevation caps.
- `plate_id`: dense current owner index. Topology `stable_ids` supplies persistent
  deterministic identities across index changes, with split/weld/removal events.
- `crust_age_myr`: chronological material age. Thermal rejuvenation can affect
  cooling without erasing that history. Raw age-step fields are storage metadata.
- `continental_crust_m`, `oceanic_crust_m`, `sediment_m`: physical column material
  volume divided by native cell area. The diagnostic `crust_thickness` field is
  a proxy; do not substitute it for these physical columns.
- `continental_area` (bundle `geosphere.continental_area`): the column's continental share of its area, 0..1; measure continental area with it, not with the plate-raster `crust_class`.
- `crust_class`, `geologic_regime`, boundary and deforming-region type maps:
  categorical classifications defined by the contract enums.
- Convergence, divergence and shear scores, uplift/subsidence tendencies,
  accumulated strain, tension and deformation rate: numerical diagnostics,
  not calibrated stresses or quantities of mineral resources.
- Boundary segment IDs and nearest IDs index the exported boundary objects;
  `boundary_distance` is in cell spacings. Zero is the absent-object sentinel.
- Boundary segments describe plate pairs, centroids, length, persistence, age,
  motion and average scores. Deforming regions describe continental rifting
  and diffuse collision, with centroid and interpolated velocity summaries.
- Tectonic junctions describe three incident plates/arms, orientation, velocity,
  geometry quality and stability, found at the cube's mesh vertices (three cells at
  a cube corner); centroids are longitude and latitude in degrees.
- Plates (`PlateKinematics`, schema 8; `tectonics_get_plate_kinematics`, one per plate
  index below `tectonics_get_plate_kinematics_count`): `rotation_rad_s`, the rotation
  vector omega in rad/s in the planet's axes (uw-cube-1: +X at 0 N 0 E, +Y at 0 N 90 E,
  +Z north), so the plate moves omega x R p at the unit vector p;
  `centroid_longitude_deg` and `centroid_latitude_deg`, the area-weighted mean unit
  vector of its cells; `velocity_east_cm_yr` and `velocity_north_cm_yr`, omega x R c at
  the centroid. The snapshot record writes them as `rotation_rad_s`,
  `centroid_lon_lat_deg` and `velocity_east_north_cm_yr` and reads them back.

Three column measurements are runtime arrays of `tectonic_api.hpp` only, outside
`Snapshot` and the snapshot files (not part of the schema); the world stores them as
`geosphere.rock_uplift_rate`, `geosphere.shortening` and
`geosphere.shortening_direction`, in run bundles only (not in the `.uww`). All three
are measurements in the plate's frame: the column model carries them with the rock
through its own remap of the material (on the transport subcells, with the
material's sub-steps), never feeding back into the columns. Converging rock mixes;
a divergent cell keeps the rock that stayed. New crust has no history of its own:
the uplift rate's mean continues the stayed rock's, the strain of new ocean floor is 0.

- `rock_uplift_rate` (`tectonics_get_rock_uplift_rate_map`, m/Myr, signed): the
  exponential mean (e-folding time `uplift_rate_window_myr`, 10 Myr) of the
  elevation change of the tectonic sub-steps (transport, subduction, rifting,
  trench flexure, arcs, plumes, roots, thermal subsidence), before diffusion,
  rivers and noise- `shortening` (`tectonics_get_shortening_map`, m, >= 0) and `shortening_direction`
  (`tectonics_get_shortening_direction_map`, radians in [0, pi)): derived at export
  from a cumulative strain tensor (xx, yy, xy, dimensionless). Each update first
  rotates the rock's tensor by the local spin W = (duy/dx - dux/dy) / 2 and adds
  the strain E = sym(grad u) of the displacement less the plate's rigid rotation
  (the distributed deformation and rift opening), the gradient by least squares over
  the cardinal neighbours on the same plate, so a velocity jump to another plate
  (subduction, ridges, transforms) is never strain and a rigid rotation turns the
  tensor (by its spin about the vertical) and adds none. The tensor travels with the
  rock as a 3-D tensor, so its components turn with the cells' frames across face
  edges. The ocean floor created in an opening starts at
  0, so a divergent cell keeps the stayed rock's tensor on the rest of its area
  only. Cumulative over the run, no decay; a sum of small strains (the strain
  reaches about 1.4 at 256 columns, seed 1, 600 Myr of the legacy recipe v3), so a small-strain measure.
  `shortening` is minus the least principal component times the engine's cell
  spacing (the mean cell side sqrt(4 pi R� / 6 N�), N = `sphere_face_cells`; 156 km
  at N 64 on Earth's radius), 0 where no axis shortens: strain x the spacing,
  metres that depend on the engine's resolution, so a consumer that needs the
  strain divides by the spacing (the world's `geosphere.shortening` is the same
  value, copied from the nearest engine cell).
  `shortening_direction` is that component's axis from east toward south, an axis
  rather than a vector (pi equals 0), 0 where the principal components are equal
  and meaningful only where `shortening` is positive. Both are taken from the
  nearest native cell into the world, never interpolated.

Column budgets expose initial/current mass, ocean/arc/plume generation,
subduction, erosion, deposition and internal redistribution in kg. Arc and plume
generation are subsets of total generated oceanic material; river erosion is a
subset of total erosion. Avoid double-counting them.

Mantle plumes are an opt-in prescribed source (`plume_processes`, default
`count` 0): plumes fixed in the native (mantle) frame add zero-age basalt at a
steady tail flux and an optional smooth head pulse, reheat oceanic lithosphere
toward a swell age, and can lower rock strength for progressive rifting. Their
volumes integrate exactly over any time step. `plume_influence` (0..1 proximity
to an active plume over the last update) is exported only when plumes are on;
its absence means plumes were off, not zero influence. Plumes are not a mantle
convection solution and add no dynamic topography beyond isostasy and the swell. The relative mass residual checks the
global budget. Motion diagnostics expose force residuals and speeds in cm/year.

## Runtime material sampling

`tectonics_get_material_displacement` exposes the native per-cell displacement
actually supplied to the latest column update ((east, north) in cell spacings per
update); it is null before
the first physical update. The mineral frame reader turns it into a 3-D tangent
vector in metres at the cell's centre. `tectonics_get_plate_uid` returns a topology stable
identity. There are no rigid-raster movement or process-global identity hooks.

`mineral_systems/tectonic_frame_reader` consumes these fields on the engine's cube
(true cell areas, the cell spacing as the unit of `boundary_distance`). Its point
tracers are unit vectors moved along great circles by the displacement field and
located on the cube; deposits leave the run as fractions of longitude and
colatitude; conservative subcell mixing,
root flow and sediment redistribution are not exact marker trajectories. The
mineral denudation proxy reads current native elevation. Per-cell material
erosion histories and physical climate coupling remain unfinished work.

## Application boundary

`geology/world_tectonics.cpp` runs the native recipe and ingests its output into
`planet` (on the world grid N is a cube world's own or W / 4 on a lat-lon world,
otherwise `sphere_face_cells`; the radius is the world body's). The hand-off
(`docs/architecture/sphere-mesh.md` section 5) is built once per engine and cached:
the identity on a cube world of the engine's N and radius; otherwise `makeremap`
weights from the engine's cells to the world's, with `remapmean` for heights and
the rock uplift rate, and per world cell the engine cell at its centre
(`remapnearest`) for codes, ids, scores, the other diagnostic fields, the snapshot's
cell fields and the shortening with its axis (a pair kept one cell's). Vectors are
(east, north) at their cell on both grids, so a nearest copy needs no rotation, and
the shortening axis is measured from the cell's east toward its south on both. In
the world's records centroids stay longitude and latitude in degrees, a segment's
length is in km, cell counts are world cells by true area, and motions, velocities
and `boundary_distance` are in world cells (per update): on lat-lon the equatorial
cell 2 pi R / W, on a cube the mean cell side (both the engine's own spacing on the
identity). Planet save storage keeps its layout; the records' centroid and length
fields now hold degrees and km. Later terrain stages may alter the application's
terrain independently of the native export.

## Serialization and checks

The current numeric field schema is 8 (the plates' rotation vector, centroid in
degrees and centroid velocity in cm/yr; 7 published a planar velocity in cell
spacings per update, a spin and a mass centre; 6 carried an always-zero cycle
counter); optional fields (`rift_damage`, `plume_influence`) do not change it. The existing serialized `phase6` name is
a format identifier, not a selectable implementation. Snapshot and batch JSON use
their existing `tectonic-*/v1` container schemas; snapshot statistics and comparisons
are computed in Python from run bundles. New
scenario metadata is version 5 (the cube engine with the remap hand-off and
contract 8; 4 introduced the cube's keys) and reports `force-balance` / `isostatic-columns`
with explicit parameters: the grid as `sphere_face_cells` and `radius_km` (no width,
height or world width), and the mantle flow's shape as `mantle_flow_wavelength_km`,
`mantle_flow_octaves`, `mantle_flow_wobble`, `mantle_flow_drift` and
`mantle_flow_shear_ratio`, and the noise on the sphere (E4) as
`initial_noise_wavelength_km`, `initial_noise_warp_km` and
`deformation_variation_wavelength_km` (each also a `--key` option of
`world_pipeline`). Deprecated scenario knobs and old snapshot generators
are removed. Previously archived bundles remain evidence and readable comparison
inputs, not a source of defaults for new worlds.

See [pipeline tooling](pipeline-tooling.md) for commands. Run the C++ module
tests, `tests/scripts/check_tectonics_cli.py` and the profile checks from the root README. Native
budget and repeatability checks establish numerical correctness only; climate
skill, geological realism and resolution convergence require separate evidence.
