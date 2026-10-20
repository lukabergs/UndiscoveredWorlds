# Pipeline Tooling

Tectonics alone is a recipe of `world_pipeline`, the only headless executable. Behaviour of the engine: [output contract](tectonic-output-contract.md).

## Command

```powershell
world_pipeline.exe run --recipe tectonics [--seed N] [--grid cube|latlon] [--resolution W] [--<scenario-key> value ...] [--duration Myr] [--mineral-history on|off] [--checkpoint-at <step|Myr>[,...]] [--checkpoint-every Myr] [--checkpoint-dir DIR] --bundle DIR
```

- Steps: `geology.geosphere_initialization` (writes the starting Geosphere) and `geology.plate_tectonics`. `--from-bundle`, `--steps` and `--only` work as for other recipes (`docs/reference-guides/benchmark-commands.md`).
- The engine runs on an equiangular cubed sphere of N cells per face edge (`sphere_face_cells`, N = 64 in `configs/tectonics-baseline.json`, so the default width is W = 4 N = 256) on the radius of the default planet body, not a body drawn from the seed. `--resolution W` sets N = W / 4 (W a multiple of 4). On the cube world grid (`--grid cube`, the default) the engine's cube is the world's grid; with `--grid latlon` the engine still runs its cube of N = W / 4 and its result is remapped onto the W x W/2 world (`src/simulations/geology/world_tectonics.cpp`). Seed 1 is the project default.
- Scenario options: every key of `scenario` in `configs/tectonics-baseline.json` and the engine's other scenario members (plumes, noise pass, rifting, ...) as `--key value`; `world_pipeline --help` lists them. A value the engine rejects fails the run with the engine's message and a non-zero exit code; an unknown option is rejected at parsing.
- On the baseline (version 6 of `configs/tectonics-baseline.json`: the Earth-like W3 set, re-tuned once on the cube in E6; compiled into the engine by `tectonic_recipe.cmake`) some keys are switched off only in pairs: `--orogenic-collapse-onset-thickness-m 0` with `--orogenic-collapse-below-onset-km2-myr 0`, `--continental-root-max-thickness-m 0` with `--continental-root-delamination-width-m 0`, and a `--rift-breakup-thickness-m` under the 15 km minimum with `--continental-min-thickness-m 0` (or a lower minimum); test `TectonicOptions.BaselineOverridesComeInPairs`.
- Mineral history runs by default; its tables are written at the end of the run only.

## Checkpoints

- `--checkpoint-at geology.geosphere_initialization` writes the starting Geosphere to `<checkpoint-dir>/00-geology.geosphere_initialization`.
- A number is a time in Myr inside the tectonics step and must be a multiple of `delta-time-myr`; it is written to `<checkpoint-dir>/01-geology.plate_tectonics-t0100.000`. `--checkpoint-every 100` writes one every 100 Myr.
- The end of the run is the `--bundle` directory, not a checkpoint. The default checkpoint directory is `<bundle>-checkpoints`.
- Every checkpoint is an ordinary run bundle.

## Bundle content

A tectonics bundle is a run bundle (schema `uw-run-bundle/v2`, `src/io/run_bundle.hpp`; v1 bundles are still read, as lat-lon): `manifest.json`, `fields/<component>/<name>.npy`, `tables/<component>/<name>.<ext>`. On the cube the manifest's `grid` block is `kind` `cube`, `face_cells` N, `orientation` `uw-cube-1`, `radius_km` and `export_columns` W = 4 N; fields have the shape [layer,] 6, N, N (face, row from the top, column); `grid/<id>.npy` holds the cell table (`centre`, `latitude_deg`, `longitude_deg`, `area_m2`, `neighbours`, `reverse`) and `latlon/<component>/<name>.npy` the cell-centred W x W/2 export of each field, by the remap rule the manifest names. A `--grid latlon` bundle holds W x W/2 fields (`docs/architecture/sphere-mesh.md`, section 7).

Fields under `geosphere.`:

- `base_elevation` (float32 m above the datum, with scalar `sea_level`); `tectonic_height` (float32, engine units: 0.1 at the ocean-floor datum, 1 at sea level on continental crust).
- `plate_id` (int32), `crust_age_steps` (int32), `crust_age` (Myr), `crust_thickness`, `crust_class`, `continental_crust`, `oceanic_crust`.
- `sediment` (float32 m, isostatic columns), `rift_damage`, `plume_influence` (0..1; zero without plumes).
- `uplift_tendency`, `subsidence_tendency`, `accumulated_strain`, `geologic_regime`.
- `boundary_type`, `boundary_distance`, `boundary_segment`, `nearest_boundary`, `deforming_region`, `deforming_region_type`.
- `deformation_rate`, `deformation_velocity_x/y`, `tectonic_convergence`, `tectonic_divergence`, `tectonic_shear`.

Tables under `geosphere.`:

- `tectonic_snapshot` (`tables/geosphere/tectonic_snapshot.json`): the record of the last tectonic update (below).
- `tectonic_boundary_segments` and `tectonic_deforming_regions` (CSV).
- `mineral_deposits` (CSV) and `mineral_budget` (JSON).

## Snapshot record

`tectonic_snapshot` has the keys `schema` (`tectonic-snapshot-bundle/v1`), `label` (`initial`, `updateNNNNNN` or `final`), `source` (`world_pipeline`), `contract_schema_version`, `update_count`, `seed`, `width`, `height` (the engine's grid: N and 6 N on the cube), `iteration_count`, `time_origin_step`, `time_myr`, `delta_time_myr`, `sea_level_m`, `run_scenario`, `motion_diagnostics`, `column_budget`, `topology_diagnostics`, `heightmap_stats`, `crust_age_steps_stats`, `crust_age_myr_stats`, `boundary_segment_count`, `deforming_region_count`, `tectonic_junction_count`, `tectonic_junction_stability_counts`, `plate_cell_counts` and `plates[]` (velocities, mass centres). The code is `src/validation/tectonics/tectonic_snapshots.*`.

## Comparing and analysing

Statistics of a snapshot (counts per geologic regime, crust class and boundary type, active convergence, divergence and shear cells, boundary segment, deforming region and junction statistics) are computed from its bundle by `uv run --offline scripts/evaluation/tectonics_statistics.py BUNDLE`; the snapshot record keeps only what the bundle cannot give, the junction count and junction stability counts. Two bundles are compared with `uv run --offline scripts/evaluation/compare_bundles.py A B`.

## Registered runs

`uw.py tectonics run` (through `scripts/evaluation/experiment.py --register`) runs the current build with the scenario of `configs/tectonics-baseline.json`, a checkpoint for the starting Geosphere and one every `checkpoint_every_updates` updates, in `runs/tectonics/<id>/`: checkpoints in `checkpoints/`, the final bundle in `bundle/`, the evaluation with `configs/evaluation/tectonics.json` in `evaluation/`. `uw.py tectonics experiment` is the unregistered scratch run in `out/evaluation/tectonics/`. Runs up to 76 are sealed in `runs/archive/` (manifests in `runs/archive/manifests/tectonics/`, snapshots in `runs/archive/checkpoints/tectonics/<id>/`). The registry's baseline (`runs/registry/tectonics.json`) is run 79: the cube at N 64 (256-column export), baseline version 6, 400 Myr, seed 1, checkpoints every 10 Myr. Runs 77 and 78 are lat-lon runs of baseline version 3 (600 Myr) from before the cube.

Registered runs up to 75 keep the former snapshot layout (raw arrays plus a manifest per snapshot and a `batch_manifest.json`); the review scripts (`scripts/benchmarks/tectonics_review.py` and its companions) read both layouts, while `uw.py tectonics report` compares evaluations and so needs runs of the bundle format (run 77 holds the fields of run 75).

## Checks

`uv run --offline python tests/scripts/check_tectonics_cli.py out/build/headless/bin/Release/world_pipeline.exe` (options the command line refuses, checkpoint times and labels, the mineral history switch); `ctest ... -L tectonics` (`TectonicOptions.*`: values the scenario options or the engine refuse); `uv run --offline scripts/uw.py tectonics experiment --duration 20 --repeat --explicit-recipe` (the profile's checks, among them repeat, recipe equality and the column budget).

Other code locations: `src/io/tectonic_heightmap_io.*` (heightmap images and metric TIFFs of the engine tests), `src/cli/tectonic_options.*` (scenario options of `world_pipeline`).
