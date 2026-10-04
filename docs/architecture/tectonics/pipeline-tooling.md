# Pipeline Tooling

Tectonics alone is a recipe of `world_pipeline`, the only headless executable. Behaviour of the engine: [output contract](tectonic-output-contract.md).

## Command

```powershell
world_pipeline.exe run --recipe tectonics [--seed N] [--resolution W] [--<scenario-key> value ...] [--duration Myr] [--mineral-history on|off] [--checkpoint-at <step|Myr>[,...]] [--checkpoint-every Myr] [--checkpoint-dir DIR] --bundle DIR
```

- Steps: `geology.geosphere_initialization` (writes the starting Geosphere) and `geology.plate_tectonics`. `--from-bundle`, `--steps` and `--only` work as for other recipes (`docs/reference-guides/benchmark-commands.md`).
- The engine runs on the world grid W x W/2 (default W = 256, the width in `configs/tectonics-baseline.json`) and uses the default planet body, not one drawn from the seed. Seed 1 is the project default.
- Scenario options: every key of `scenario` in `configs/tectonics-baseline.json` and the engine's other scenario members (plumes, noise pass, rifting, ...) as `--key value`; `world_pipeline --help` lists them. A value the engine rejects fails the run with the engine's message and a non-zero exit code; an unknown option is rejected at parsing.
- Mineral history runs by default; its tables are written at the end of the run only.

## Checkpoints

- `--checkpoint-at geology.geosphere_initialization` writes the starting Geosphere to `<checkpoint-dir>/00-geology.geosphere_initialization`.
- A number is a time in Myr inside the tectonics step and must be a multiple of `delta-time-myr`; it is written to `<checkpoint-dir>/01-geology.plate_tectonics-t0100.000`. `--checkpoint-every 100` writes one every 100 Myr.
- The end of the run is the `--bundle` directory, not a checkpoint. The default checkpoint directory is `<bundle>-checkpoints`.
- Every checkpoint is an ordinary run bundle.

## Bundle content

A tectonics bundle is a run bundle (schema `uw-run-bundle/v1`: `manifest.json`, `fields/<component>/<name>.npy`, `tables/<component>/<name>.<ext>`).

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

`tectonic_snapshot` has the keys `schema` (`tectonic-snapshot-bundle/v1`), `label` (`initial`, `updateNNNNNN` or `final`), `source` (`world_pipeline`), `contract_schema_version`, `update_count`, `seed`, `width`, `height`, `iteration_count`, `time_origin_step`, `time_myr`, `delta_time_myr`, `sea_level_m`, `run_scenario`, `motion_diagnostics`, `column_budget`, `topology_diagnostics`, `heightmap_stats`, crust-age statistics, `summary`, `plate_cell_counts` and `plates[]` (velocities, mass centres). `summary` holds the regime and crust-class counts, boundary and deforming-region statistics and junction counts. The code is `src/validation/tectonics/tectonic_snapshots.*`.

## Comparing and analysing

Statistics are in the `summary` block of the snapshot record. Two bundles are compared with `uv run --offline scripts/evaluation/compare_bundles.py A B`.

## Registered runs

`uw.py tectonics run` (through `scripts/benchmarks/tectonics_experiment.py`) copies `world_pipeline.exe` into `runs/reports/tectonics/<id>/` and runs it from there with a checkpoint for the starting Geosphere and one every `checkpoint_every_updates` updates. Checkpoints land in `runs/checkpoints/tectonics/<id>/snapshots/`, the final bundle in `runs/checkpoints/tectonics/<id>/final`. `uw.py tectonics inspect` writes an unregistered final-state bundle to `out/bench/`. `uw.py tectonics experiment|evaluate` use `configs/evaluation/tectonics.json`. The baseline run is registered run 75.

Registered runs up to 75 keep the former snapshot layout (raw arrays plus a manifest per snapshot and a `batch_manifest.json`); `uw.py tectonics report` and the report scripts read both layouts.

## Checks

`uv run --offline python tests/scripts/check_tectonics_cli.py out/build/headless/bin/Release/world_pipeline.exe`

Other code locations: `src/io/tectonic_heightmap_io.*` (heightmap images and metric TIFFs of the engine tests), `src/app/tectonic_options.*` (scenario options of `world_pipeline`).
