# Pipeline Tooling

`tectonic_pipeline` is the sole world-simulation CLI for the in-repo engine. All commands start from the compiled `configs/tectonics-baseline.json` recipe.

## Structured Snapshot Bundles

`tectonic_pipeline batch-export` writes one batch directory with:

- `batch_manifest.json`
- `snapshots/<NNN>_<label>/manifest.json`
- per-snapshot binary field files
- `boundary_segments.json`
- `deforming_regions.json`
- `tectonic_junctions.json`

With `--mineral-history on`, the batch directory also holds `mineral_deposits.csv`
and `mineral_budget.json`; see [mineral systems](../mineral-systems.md).

The per-snapshot bundle format is:

- `schema`: `tectonic-snapshot-bundle/v1`
- `contract_schema_version`: current native contract schema version
- fixed binary filenames for all exported raster fields
- JSON object exports for boundary segments, deforming regions, and tectonic junctions
- per-junction stability diagnostics in `tectonic_junctions.json`
- hashes, scalar summaries, and per-plate motion metadata

The batch manifest format is:

- `schema`: `tectonic-batch-export/v1`
- full requested scenario metadata
- checkpoint request metadata
- exported snapshot directory list
- unresolved requested checkpoints, if any

This stays binary-or-JSON-backed on purpose. The tool does not add georeferenced interchange formats.

## Inspect And Debug Commands

`tectonic_pipeline` exposes four commands:

- `batch-export`: export initial, final, periodic, and explicit time/update checkpoints
- `inspect`: run a seed/scenario to a checkpoint and emit `tectonic-snapshot-summary/v1`
- `boundary-stats`: emit `tectonic-boundary-stats/v1`
- `compare`: emit `tectonic-snapshot-comparison/v1`

Supported checkpoint selectors:

- `--final`
- `--update N`
- `--time-myr X`

`--time-myr` must align exactly to `delta_time_myr`.

## Compare Modes

`tectonic_pipeline compare` supports:

- direct seed-to-seed compare with shared scenario parameters:
  - `--seed-a`
  - `--seed-b`
- bundle-to-bundle compare:
  - `--bundle-a`
  - `--bundle-b`

The bundle compare path is the parameter-set compare path for downstream regression work: export two snapshot bundles from different runs, then compare them.

## Time semantics

`time_myr` is continuous elapsed time; one update advances `delta_time_myr`.
`--evolution-duration-myr` must be finite and positive. The final step may round
the requested duration upward.

```powershell
out/build/headless/bin/Release/tectonic_pipeline.exe batch-export `
  --output out/tectonics-batch --seed 1 `
  --checkpoint-every-updates 25 --checkpoint-at-time-myr 50
```

Use a new output directory. Rebuild after recipe or source changes. See the root
README for build commands and the distinction between native and app grids.
