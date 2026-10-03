# Standard climate reports

`scripts/benchmarks/climate_report.py` is the common batch generator. Use numeric run IDs. The title is `(BASELINE) - FIRST-LAST`, including the baseline when it is in the displayed run list. The input adapter reads already analyzed numerical fields from a batch archive; historical filenames such as `T10-fields.npz` are private input aliases, never case labels. New analysis files should use the run ID directly.

```powershell
uv run --offline scripts/benchmarks/climate_report.py --archive runs/reports/climate-transport-experiments-20260918 --baseline 353 --runs 353 354 355 356 357 358 359 360 361 362 363 --apply-selection
```

The first build renders maps. `--reuse-maps` reuses existing images for a metadata/layout rebuild; omit it after changing palettes, numeric inputs, display scales or rendering code. `--reuse-maps --rebuild-scalars` refreshes scalar/magnitude images while keeping complete ocean textures, for scalar-rendering-only changes. No simulation is executed. Bulk generated reports route to the configured D: storage automatically. The reference data and reference images remain under `refs/processed` on C:.

Each input run has an `*-analysis.json` containing its registered `run` ID and captured circulation paths, a companion `*-fields.npz` and `*-comparison.npz`, and the original native diagnostics. The archive also supplies `references.npz`, `reference-provenance.json` and `mixed-layer-reference.npz`. The generator checks required fields and common geometry rather than silently ranking incomplete cases. Existing native analysis scripts provide this adapter contract; numerical extraction remains separate from report layout.

## Artifact roles

- **HTML comparison:** selectors, baseline/selected/reference maps, global tables, distributions and latitude profiles. Regional crops are visual only and do not change global measurements. Dateline-crossing regions wrap correctly. Download links point at the exact full-size image used by the gallery.
- **Markdown report:** hypotheses, source changes, findings, decisions, numerical recommendation, score definition and limitations. Generated from the run's recorded decision data.
- **Execution log:** raw stdout/stderr, solver output, warnings, timings and failures. It is evidence for debugging, not a second narrative report.
- **JSON receipts:** exact metrics, ranking inputs, source provenance, palette version and hashes. These support later retrieval and rebuilding.

The benchmark wrapper writes a compact per-run record using `climate_run_record.py`; after analysis, the shared batch generator links pending run pages to the appropriate comparison with that run selected. The old layout is available only through `write-benchmark-report.py --legacy-layout` for historical reproduction. Sealed old reports/maps are retained unchanged.

## Palettes

`configs/climate-palettes.json` is the single source for colours and scales. `climate_palettes.py` serves Python exporters and plots. Generate the C++ header after changing the definition:

```powershell
uv run --offline python scripts/benchmarks/generate-climate-palette-header.py
uv run --offline python scripts/benchmarks/generate-climate-palette-header.py --check
```

Positive magnitudes use blue–cyan–green–yellow–orange–red. Signed quantities use white–magenta–blue–cyan–green–yellow–red, with **zero always cyan**, even on asymmetric temperature ranges. Black means missing or masked data. Display endpoints saturate; statistics use unclipped physical values and include overflow/underflow bins. Categorical Köppen colours are unchanged.

The new reference set is `refs/processed/climate/maps/standard-v2/index.html`. Native exporters use these colours on future runs; old sealed images keep their original palettes. Standardized re-exports of the displayed historical runs reside in the new comparison archive. Ocean LIC and particles use the same 45-day, coast-stopping spherical integration for model and reference. LIC/particle texture changes pixel brightness, so exact value reading belongs to magnitude maps.

The display selector supports surface, 850-hPa and 500-hPa wind LIC/particles as well as currents. Wind particles use identical 4.5-day integrations and seed 20260906 for model and reference. **Smooth pixels** switches scalar/magnitude panels between interpolated and native-cell exports; on flow textures it controls display interpolation. **Continent outlines** draws a separate transparent geographic overlay. Neither switch affects metrics.

Build the optional numerical-equivalent rendering accelerator before large report batches:

```powershell
cmake --build out/build/headless --config Release --target climate_flow_texture --parallel 4
```

Python falls back to the original NumPy kernels when the accelerator is absent. `test_climate_standard_report.py` compares the compiled output with those kernels, including masked cells. The helper affects rendering only.

## Ranking and highlighting

The fixed `global-balanced-v1` score is **lower-is-better**. It gives rainfall 25% (annual IMERG 15%, seasonal ERA5 10%), winds 22% (surface 10%, 850 hPa 6%, 500 hPa 6%), pressure 12%, OISST SST 11%, currents 10%, moisture transport 9%, column water 6%, and moisture convergence 5%. Seasonal fields weight all four seasons equally.

These are explicit judgement weights reflecting the project's priorities, not a universal or fitted climate-skill index. Spatial errors prevent a correct histogram in the wrong places from winning; tail and distribution terms expose errors hidden by global means; vector direction distinguishes correct circulation from speed alone. The ocean-current component remains an approximate comparison of model effective-layer velocities with reference near-surface currents.

Scalar loss: 50% normalized spatial RMSE, 20% normalized spatial error P90, 20% area in the wrong physical-value bin, 10% histogram total variation. Normalize errors by reference P95–P05 spread, with a floor of 1% of the fixed display span.

Vector loss: 45% vector RMSE/reference RMS speed, 20% vector error P90/reference speed P90, 20% direction error/180°, 15% speed-distribution total variation. The P90 denominator has a floor of 10% of reference RMS speed. Direction uses a fixed reference-active mask; a calm candidate receives 90° error there. This prevents reducing speed to evade the direction comparison. East/north and speed errors remain visible as separate table rows.

Multiply the weighted mean loss by 100 to obtain the displayed cost. Compare all runs on fixed common valid reference/terrain masks and spherical area weights. Missing required fields stop ranking. There is no Köppen component, regional scoring or claim of equilibrium. Density mixed-layer depth differs from thermal storage depth, and an exact skin-temperature reference is absent; these auxiliary displays have zero weight. ERA5 SST is a secondary check rather than a second SST score.

`--apply-selection` updates `climate-global-current-selection.json` only to the batch winner if its cost is strictly lower than the incoming baseline; ties retain the baseline. The recommendation is numerical, not a substitute for visual review. An explicit user visual choice should replace that selection and be recorded with its reason before the next batch. Updating the baseline record does not silently restore simulation source; the next batch must start from the selected run's archived physics.

To record that choice, add `--override-baseline ID --override-reason "User's reason" --apply-selection --reuse-maps` to the same batch command. A subsequent ordinary rebuild of that batch preserves its recorded visual override. A new batch can generate a new numerical recommendation.

Standalone reference regeneration, independent of any experiment batch:

```powershell
uv run --offline scripts/refs/render-standard-climate-maps.py
uv run --offline scripts/refs/render-standard-flow-references.py
```

In the table, **bold** marks the closer of baseline and selected run to the reference target. **Green** additionally means best across all displayed batch runs and the baseline. Ties share highlighting. Percentiles and histogram proportions compare values with the reference; error rows compare with zero, correlations with one and agreement fractions with 100%.
