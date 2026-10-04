# Native tectonic output contract

The implementation is `src/simulations/geology/plate_tectonics` in this working
tree. `tectonic_contract.hpp` defines `tectonics::contract::Snapshot` and its
objects. `capture_snapshot(lithosphere)` copies a state; `tectonic_api.hpp`
exposes borrowed runtime arrays. The app and the `tectonics` recipe of `world_pipeline` use the same
engine and the compiled `configs/tectonics-baseline.json` recipe.

## Geometry and time

Native rasters are cell-centred, row-major width × height, on a planar torus.
Both axes wrap; native cells have equal area. Cell spacing is
`world_width_km / width`, not a latitude-dependent spherical distance.
Coordinates and kinematics are native cells and cells/update unless explicitly
named in physical units. Array pointers expire when the engine updates or dies.

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
- `crust_class`, `geologic_regime`, boundary and deforming-region type maps:
  categorical classifications defined by the contract enums.
- Convergence, divergence and shear scores, uplift/subsidence tendencies,
  accumulated strain, tension and deformation rate: numerical diagnostics,
  not calibrated stresses or quantities of mineral resources.
- Boundary segment IDs and nearest IDs index the exported boundary objects;
  `boundary_distance` is in native cells. Zero is the absent-object sentinel.
- Boundary segments describe plate pairs, centroids, length, persistence, age,
  motion and average scores. Deforming regions describe continental rifting
  and diffuse collision, with centroid and interpolated velocity summaries.
- Tectonic junctions describe three incident plates/arms, orientation, velocity,
  geometry quality and stability. Their diagnostics also inform local processes;
  they do not constitute a spherical plate mesh.

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
actually supplied to the latest column update (cells/update); it is null before
the first physical update. `tectonics_get_plate_uid` returns a topology stable
identity. There are no rigid-raster movement or process-global identity hooks.

`mineral_systems/tectonic_frame_reader` consumes these fields. Its point tracers
approximate movement through the displacement field; conservative subcell mixing,
root flow and sediment redistribution are not exact marker trajectories. The
mineral denudation proxy reads current native elevation. Per-cell material
erosion histories and physical climate coupling remain unfinished work.

## Application boundary

`geology/world_tectonics.cpp` runs the native recipe and ingests its output into
`planet`. `tectonic_world_grid.hpp` defines the cell-centred sampling: bilinear
terrain with longitude wrapping and clamped polar rows; nearest-cell categories
and diagnostic fields. Coordinates, distances, lengths and kinematics convert
to world-grid cells. Object counts scale by cell area. Planet save storage is
unchanged. The world raster is interpreted by spherical climate/hydrology;
this sampling is not spherical tectonic physics. Later terrain stages may alter
the application's terrain independently of the native export.

## Serialization and checks

The current numeric field schema is 7 (6 carried an always-zero cycle counter);
optional fields (`rift_damage`, `plume_influence`) do not change it. The existing serialized `phase6` name is
a format identifier, not a selectable implementation. Snapshot, batch, summary
and comparison JSON use their existing `tectonic-*/v1` container schemas. New
scenario metadata is version 3 and reports `force-balance` / `isostatic-columns`
with explicit parameters. Deprecated scenario knobs and old snapshot generators
are removed. Previously archived bundles remain evidence and readable comparison
inputs, not a source of defaults for new worlds.

See [pipeline tooling](pipeline-tooling.md) for commands. Run the C++ module
tests and `tests/scripts/check_tectonics_cli.py` from the root README. Native
budget and repeatability checks establish numerical correctness only; climate
skill, geological realism and resolution convergence require separate evidence.
