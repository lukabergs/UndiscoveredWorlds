# Climate runs and reports

Every step goes through `scripts/uw.py climate <command>`; pages are rendered from JSON by `scripts/evaluation/report.py` and its three templates (`scripts/evaluation/templates/`). Options and the run directory layout are in `benchmark-commands.md`.

```powershell
# Registered run (seed 1, registered settings of configs/evaluation/climate.json) -> runs/climate/<id>/, evaluated, scored, page written
uv run --offline scripts/uw.py climate run --label "what changed"
# Scratch run with the profile's experiment settings -> out/evaluation/climate/<name>/
uv run --offline scripts/uw.py climate experiment --resolution 128 --climate-resolution 64
# Evaluate an existing bundle against the baseline bundle (registry baseline: run 388)
uv run --offline scripts/uw.py climate evaluate --bundle runs/climate/<id>/bundle --baseline runs/climate/388/bundle --output out/evaluation/climate/<id>
# Review: the verifier agent first (reviewer agent), then the user's verdict, which overrides
uv run --offline scripts/uw.py climate review --evaluation runs/climate/<id>/evaluation --verdict approved --notes "..." --reviewer user
# Comparison page of evaluations, the first being the baseline
uv run --offline scripts/uw.py climate report --evaluation runs/climate/388/evaluation runs/climate/<id>/evaluation --output runs/reports/<name>
# Dashboard over the registry and the ledger
uv run --offline scripts/uw.py climate dashboard --output out/evaluation/climate/dashboard
# Ledger: entries on a topic, a draft entry, append a completed entry (LEDGER.md is regenerated)
uv run --offline scripts/uw.py climate ledger --topic transport
uv run --offline scripts/uw.py climate ledger --draft
uv run --offline scripts/uw.py climate ledger --add entry.json
```

An evaluation directory holds `evaluation.json`, `review.json`, `report.md`, `index.html` (run page) and `maps/` (model, reference and difference in the palettes of `configs/climate-palettes.json`). It is never overwritten; pages are rendered again from it.

Runs registered before run 388 were made and reported by a former script chain, retired on 2026-10-05 after the parity check of ledger `climate-0002`. Their data are sealed in `runs/archive/` (`runs/archive/README.md`); their pages (`runs/reports/<id>.html` and `.md`, `runs/reports/climate-comparison-*.html` with their folders) stay as written and are not regenerated.

## Score: global-balanced-v1

`scripts/evaluation/score.py`. The score is **lower-is-better**. It gives rainfall 25% (annual IMERG 15%, seasonal ERA5 10%), winds 22% (surface 10%, 850 hPa 6%, 500 hPa 6%), pressure 12%, OISST SST 11%, currents 10%, moisture transport 9%, column water 6%, and moisture convergence 5%. Seasonal fields weight all four seasons equally.

These are explicit judgement weights reflecting the project's priorities, not a universal or fitted climate-skill index. Spatial errors prevent a correct histogram in the wrong places from winning; tail and distribution terms expose errors hidden by global means; vector direction distinguishes correct circulation from speed alone. The ocean-current component remains an approximate comparison of model effective-layer velocities with reference near-surface currents.

Scalar loss: 50% normalized spatial RMSE, 20% normalized spatial error P90, 20% area in the wrong physical-value bin, 10% histogram total variation. Normalize errors by reference P95–P05 spread, with a floor of 1% of the fixed display span.

Vector loss: 45% vector RMSE/reference RMS speed, 20% vector error P90/reference speed P90, 20% direction error/180°, 15% speed-distribution total variation. The P90 denominator has a floor of 10% of reference RMS speed. Direction uses a fixed reference-active mask; a calm candidate receives 90° error there. This prevents reducing speed to evade the direction comparison.

Multiply the weighted mean loss by 100 to obtain the displayed cost. Compare runs on fixed common valid reference/terrain masks and spherical area weights. There is no Köppen component, regional scoring or claim of equilibrium. Density mixed-layer depth differs from thermal storage depth, and an exact skin-temperature reference is absent; these have zero weight. ERA5 SST is a secondary check rather than a second SST score. The score is a numerical recommendation; the user's visual verdict overrides it.

## Palettes

`configs/climate-palettes.json` is the single source for colours and scales. `scripts/benchmarks/climate_palettes.py` serves Python renderers. Generate the C++ header after changing the definition:

```powershell
uv run --offline python scripts/benchmarks/generate_climate_palette_header.py
uv run --offline python scripts/benchmarks/generate_climate_palette_header.py --check
```

Positive magnitudes use blue–cyan–green–yellow–orange–red. Signed quantities use white–magenta–blue–cyan–green–yellow–red, with **zero always cyan**, even on asymmetric temperature ranges. Black means missing or masked data. Display endpoints saturate; statistics use unclipped physical values and include overflow/underflow bins. Categorical Köppen colours are unchanged.

## Standard reference images

`refs/processed/climate/maps/standard-v2/index.html` is a standalone set of 2048-wide reference images (scalar maps, wind and ocean-current LIC and particles; 45-day ocean and 4.5-day wind particle integrations, texture seed 20260906). The renderers are in `scripts/refs/standard_climate_rendering.py`:

```powershell
uv run --offline scripts/refs/render_standard_climate_maps.py
uv run --offline scripts/refs/render_standard_flow_references.py
```

The optional compiled texture kernel speeds up LIC and particles; Python falls back to the NumPy kernels when it is absent, and `tests/scripts/test_climate_standard_report.py` compares the two:

```powershell
cmake --build out/build/headless --config Release --target climate_flow_texture --parallel 4
```
