# Climate simulation handoff — GPT-6 Astra, Ultra

> **Record since 2026-10-20.** Its open items moved to `docs/architecture/master-plan.md` (the only running plan; §7 maps the old ids to the new steps). Kept for what was decided and learned.

Prepared 2026-09-07. This is a self-contained task prompt for an instance with access to the repository. Verify the current checkout before relying on this snapshot.

**Latest continuation completed 2026-09-07:** [historical comparison and ocean corrections](../../../runs/reports/climate-history-ocean-corrections.md). Runs 179/180 isolate wet-only ocean inputs and spherical ocean stress curl at 128×64; [181 versus 178, reference and historical 145](../../../runs/reports/181.html) is the final 512×256 gallery (78.186 s). Both correctness fixes are retained, but rainfall/moisture/SST regress and the Thai current hotspot worsens; land annual-temperature RMSE improves slightly. Ten strict reference targets fail, and the coupled climate remains unconverged. Checkpoint `10b4aba` contains matching records and source for 142–145, without exact executable provenance. Next isolate shallow South China Sea transport and full-depth mixed-layer heat closure; do not interpret stronger currents or inner-solver convergence as better climate. Historical dry-land percentages mean annual rainfall below **12 mm**, not 1 mm.

**Further continuation completed 2026-09-07:** [coastal winds, highland ice and ocean momentum corrections](../../../runs/reports/climate-surface-corrections.md). Runs 175–177 isolate the three changes at 128×64; [178 versus 174](../../../runs/reports/178.html) is the final 512×256 comparison (79.138 s). Coastal land/ocean wind contrast and mountain temperatures improve; stronger currents worsen cold SST and the Thai hotspot, and Arctic land warms further. Legacy remains default, the strict gate fails seven targets, bands persist, and the climate is not converged. Next investigate tropical ocean heat/pressure coupling and wet forcing/depth reconstruction; the older failed tropical closure is still unvalidated.

**Continuation completed 2026-09-07:** [ocean thermal remap correction and audit](../../../runs/reports/climate-ocean-mask-correction.md). Runs 172/173 isolate the correction at 128×64; run 174 is the corrected 512×256 gallery against 167. Legacy remains default. Rainfall/moisture/SST improve, current tradeoffs and the incorrect tropical SST gradient remain explicit, and none of these runs is converged. The snapshot below describes the earlier 167/171 state.

## Task

Continue improving the climate simulation in `C:\dev\undiscovered-worlds`. Investigate the failed tropical boundary-layer replacement, correct the underlying coupled model where justified, and implement a bounded, tested improvement. Do not stop at another literature review or a plan. You may use sub-agents for independent investigations; coordinate ownership of edits and run benchmarks serially.

## Goal

Build a reasonably accurate, inexpensive seasonal climate baseline for a game: rainfall, surface and upper circulation, ocean currents, air and sea temperatures, and resulting Köppen classification. The user considers historical runs 142–145 the closest visual climate classifications so far. Recover useful skill without copying an Earth reference into the simulated planet or optimizing only Köppen percentages.

Cyclones, Rossby waves, storm tracks and other transient weather belong primarily to a later simplified dynamic simulator. The baseline still needs geographically responsive winds, regional convergence, monsoonal structure and plausible ocean heat transport. Adding arbitrary turbulence to LICs is not a climate fix.

## Constraints

- Work directly in the current branch and preserve existing changes. At handoff `git branch --show-current` returned `main`; the working tree contains a large preexisting staged/unstaged reorganization and many untracked implementation files. Do not reset, clean, broadly revert, switch checkouts or assume all changes belong to this task. No commit was created for the tropical trial. Follow the user's AGENTS.md approval boundaries; do not push, merge or open a PR.
- Keep diffs focused and deterministic. Reuse the Windows/PowerShell/CMake/uv baseline. No new production dependency or save-format/public-contract changes without the required approval.
- The user rejects multi-hour climate runs and eight-minute 128-wide runs. Start with saved diagnostics and small deterministic fixtures; use a bounded 128×64 run when needed. A recent 512×256 **export** benchmark takes about 80 seconds, but it is not converged. Do not promise convergence in that time.
- Default workload: daily-mean `baseline`, conservative `ffsl`, eight outer updates, one carried hydrology year per update plus the initial year: at most nine hydrology years. Avoid returning to nested spin-up or running 24–100 iterations blindly. Do not begin with native 512-wide atmosphere/ocean physics.
- Output must match the requested dimensions. Coarser component grids are acceptable only with explicit continuous field reconstruction before classification/export. Do not undo the accepted seam fix or hide fields with image-only blur.
- Use `runs/` for benchmarks and `out/` for builds. Keep the minimal eleven PNGs: annual rainfall/air temperature/Köppen and Jan/Jul surface LIC, particles, SST and ocean currents. Keep other quantities as numerical diagnostics. Preserve IDs and the current folder schema.
- Do not promote the experimental closure just because its linear solver or conservation tests pass. A major rainfall failure cannot be traded away silently for a smaller wind RMSE. Explain smaller cross-variable tradeoffs rather than hiding them in a composite score.

## Verified current state

**Production/default:** `-TropicalClosure legacy`. **Experimental:** `-TropicalClosure mixed-layer`, or executable `--climate-tropical-closure mixed-layer`.

The previous effort implemented the requested replacement, tested it, and rejected default promotion because coupled rainfall deteriorated severely. The implementation is available for inspection and controlled comparison, not a validated solution to defend.

The seam fix from run 167 is accepted by the user. `makeSmoothRemap` / `remapSmoothField` in `src/core/climate_grid.*` integrate continuous transitions over destination cells in longitude and sin(latitude). They preserve constants, bounds and the global area-weighted integral. Neighboring coarse cells may exchange their reconstructed means; individual coarse-cell means are not preserved during refinement. Same-size grids and restriction retain their original behavior. This is reconstruction, not newly resolved physics.

Runs 167 and 171 used seed **20260906**, **512×256 exports**, **64×32 atmosphere/ocean**, **128×64 hydrology/energy**, eight coupling updates and nine hydrology years. Both were **not converged**. Their controlled comparison is more useful than comparing directly with historical 145.

| Metric | Legacy 167 | Experimental 171 |
| --- | ---: | ---: |
| Wall time including exports | 79.24 s | 81.66 s |
| Area-weighted land rainfall correlation | 0.4369 | 0.0472 |
| ERA5 column-water correlation | 0.5815 | 0.1375 |
| Surface eastward-wind correlation | 0.7081 | 0.6493 |
| Surface northward-wind correlation | 0.4791 | 0.0797 |
| Tropical vector-wind RMSE, Q1 / Q3 | 6.45 / 6.35 m/s | 3.67 / 5.21 m/s |
| Tropical zonal fraction of squared wind, Q1 / Q3 | 94.3% / 93.1% | 81.5% / 80.7% |
| Ocean SST RMSE, Q1 / Q3 | 6.00 / 6.06 °C | 5.08 / 5.45 °C |
| Ocean mean SST bias, Q1 / Q3 | −2.65 / −2.99 °C | +0.29 / −0.06 °C |

Tropical metrics cover 20°S–20°N, land and ocean, with area weighting. ERA5's zonal fractions are **53.7% / 54.1%**. Regional vector RMS remains deficient: 171 gives **1.70 / 1.85 m/s**, versus ERA5 **3.27 / 3.62**. Lower total wind error therefore does not prove better regional circulation. Current-component skill is mixed: three of four checked current-component RMSEs worsened. Near-zero mean SST bias conceals substantial regional temperature errors.

Run 171 hydrology consumed 48.50 s, including 36.11 s of transport and zero transient-weather integration. The maximum water relative residual was 9.06e−9. Its final changes were wind 0.00505, SST 0.00821, heating 0.11205 and rainfall 0.03157; `inner_accepted=0`, `converged=0`. All tropical inversion solves were accepted, which is different from overall coupled convergence.

Other retained trials, not alternative implementations selectable today:

- **168:** heating-only inversion forcing, two updates, 42.9 s. Trades collapsed; rainfall correlation −0.076.
- **169:** full old overturning amplitude as inversion equilibrium forcing, two updates, 42.2 s. Tropical wind RMS became 10.4/9.5 m/s. Too strong.
- **170:** 0.4 amplitude projection, eight updates, 81.9 s. Wind improved; rainfall correlation 0.077. Used the old wave-depth ascent diagnosis.
- **171:** same projection with physical boundary-layer ascent mass, eight updates, 81.7 s. Rainfall correlation fell further to 0.047. Changing ascent depth alone did not solve the coupled problem.

## What the new closure actually does

`src/simulations/climate/climate_boundary_layer.cpp`, declared in `climate_atmosphere.hpp`, implements a steady mixed-layer momentum balance:

`f k × U + grad(Ps)/rho = e(U850 − U) − d U`.

Ocean coefficients are `e=2e−5 /s`, `d=1.5e−5 /s`; land/relief modifies friction. Full strength within 20° latitude fades to legacy by 35°. Pressure gauges are aligned before blending. Linear entrainment/friction blends with the legacy quadratic drag; the equator is handled without dividing by f. The old 500 hPa stationary momentum exchange is suppressed inside the fully active tropical layer.

The mechanism follows Back & Bretherton (2009), https://aos.wisc.edu/~lback/SSTgradients.pdf, DOI 10.1175/2009JCLI2392.1. **Their model uses observed/reanalysis 850 hPa temperature, pressure and winds. Our invented diagnosis of those inputs is not validated by that paper.** Audit it independently.

Our 850 hPa state is a separate forced/damped pressure solve. A temperature anomaly is obtained from net column heating times the existing surface damping time, divided by column heat capacity. Inversion temperature is a global surface-temperature mean minus a lapse-rate height correction plus that anomaly. Inversion equilibrium pressure includes the heating response, mechanical terrain forcing, and **0.4 times the old zonal overturning profile**. The 0.4 is an empirical trial projection, not a derived constant; inversion drag time is five days. Surface thermal thickness assumes a linear anomaly profile from surface to inversion. No observational fields are injected into production.

In `climate_atmosphere.cpp`, ascent in the tropical core is the convergence of a physical **150 hPa** boundary-layer mass flux. Through 20–35°, its flux blends with legacy equivalent-depth fluxes before spherical divergence. The global ascent integral closes. This does **not** close the full column mass budget or make the independently diagnosed 500 hPa winds a consistent return circulation. The pressure/flux blend can itself introduce transition artifacts; inspect it rather than assuming smooth weighting is physically consistent.

`climate_coupling.cpp` supplies surface temperature, heating and terrain, applies the new surface diagnosis, then passes surface/upper winds into hydrology and ocean coupling. Hydrology still transports two moisture stores with the existing upper-flow/departure blend. Its 0.45/0.55 moisture-capacity fractions are not automatically the physical layer pressure masses. Keep these concepts separate.

## Findings that need causal isolation

1. **Upper transport versus new surface flow.** In saved 171 weather statistics near 3°S in Q3, the zonal-mean upper transport is about 15.1 m/s southward, versus 2.1 m/s northward at the surface. Upper moisture export remains strong after boundary inflow weakens. This is evidence of a coupling problem, not proof that a velocity rescale or a particular return-flow correction is sufficient. Diagnose dry-mass continuity, actual moisture-flux divergence and vertical exchange together. Do not infer mass balance from velocities alone.
2. **Wrong thermal gradients.** Ocean-only row profiles show a cold tropical ocean and excessively warm southern subtropics. A thermal-gradient-driven closure can faithfully respond to the wrong SST. Diagnose surface radiative/turbulent fluxes, ocean advection, Ekman/upwelling terms, mixing, ice and forcing/remapping separately. Exclude land zeros when averaging ocean CSVs. A small global bias is not a correct SST pattern.
3. **Unjustified inversion approximation.** The 850 hPa pressure, temperature and winds are weakly constrained and retain projected zonal forcing. Check vertical hydrostatic consistency, pressure/temperature gradients, damping/equivalent-depth choices and the response to topography. Do not fix the 0.4 projection as unquestionable architecture.
4. **Transition effects.** Blending pressure profiles and changing mass-depth weights in the 20–35° belt may generate artificial gradients/divergence and contribute to bands. Global conservation and gauge invariance do not exclude those artifacts.
5. **Rainfall closure compatibility.** Convection depends on moisture convergence, evaporation, upper accession and two-layer saturation/exchange. The new surface flow operates against a hydrology model tuned with the older flow. Inspect the actual budgets before adjusting rain efficiency or convective gains. These hypotheses have not yet been isolated with controlled component experiments.

You may replace or simplify the unsuccessful prototype if a better justified, cheaper closure meets the game's objective. There is no requirement to retain this particular implementation or to restore every historical heuristic.

## Relevant files and evidence

All paths below are relative to `C:\dev\undiscovered-worlds`.

- `runs/maps/climate/README.md`: current commands, reconstruction contract, trial results, limitations and artifact schema.
- `runs/reports/167.html`, `171.html`, `167-checks.json`, `171-checks.json`: galleries and measured baseline/candidate scores. Run 171's gallery compares against 167 and visibly reports the failed reference gate.
- `runs/diagnostics/climate/{167,171}/maps/`: `weather_statistics.csv`, `climate_process_fields.csv`, `ocean_fields.csv`; 171 also has `tropical_boundary_layer.csv`. These are retained native component fields. `runs/fields/climate/rain/annual/s/<ID>.tif` contains unclamped exported rainfall.
- `runs/logs/climate/{167,171}.log`: iteration status and subsystem timing. `runs/registry/climate.json` contains historical metrics and experiment descriptions.
- `src/simulations/climate/climate_boundary_layer.cpp`: inversion diagnosis, thermal thickness, mixed drag solver and application. `climate_atmosphere.cpp`: `diagnoseModeWinds`, momentum exchange, blended mass flux and ascent. `climate_coupling.cpp`: `createvectorwindmap`, ocean forcing/coupling and moisture transport/precipitation loop. Follow the corresponding `climate_ocean_dynamics.*`, `climate_energy.*`, `climate_hydrology.*`, `climate_transport.*`, `climate_physics.*` and tuning values as needed.
- `src/simulations/climate/climate_context.hpp`, `src/app/main.cpp`, `scripts/benchmarks/run-climate-benchmark.ps1`: experimental/default selection. Context and CLI defaults currently have `mixedLayerTropics=false`; wrapper default is `legacy`.
- `src/validation/climate/climate_circulation_diagnostics.cpp`: numerical exports. `scripts/benchmarks/check-climate-baseline.py`: numerical checks, tropical decomposition and reference scoring. `write-benchmark-report.py`: gallery generation.
- `tests/simulations/climate/climate_boundary_layer_tests.cpp`, `climate_context_tests.cpp`, the other climate tests, and `tests/scripts/test_climate_baseline_checks.py`.
- `refs/processed/climate/maps/README.md` and `index.html`: reference semantics and exact palettes. Reference bundles include ERA5 winds/pressure/moisture, OISST and GLORYS. Inspect what exists before assuming a particular level/variable is available.

CSV seasons are **0–3**; map folders are **1=January, 2=April, 3=July, 4=October**. The checker compares day-weighted Jan–Mar and Jul–Sep reference means. Representative-month previews are not identical to quarter averages. CSV meridional wind/current is generally **southward positive**; ERA5/GLORYS northward components require a sign change. Check each schema. Simulator SST and reference SST PNGs use different palettes: compare numerical fields, not equal RGB values. Ocean `sst_c` is an endpoint, while `mean_sst_c` is the quarterly mean used for comparison. GLORYS shallowest currents and the simulated mixed-layer currents are related but not identical quantities.

The registry reports zero independently compared spatial Köppen cells for relevant runs. Do not claim validated spatial classification from aggregate percentages. Historical reference Köppen uses a −3 °C C/D boundary while the simulator uses 0 °C; separate convention mismatch from simulation error.

### Historical baseline

Read `runs/reports/climate-history-140-160.md`. The best recovered source checkpoint for 142–145 is **`10b4aba`**. It establishes a recoverable snapshot, not proof of the exact clean source used by every benchmark executable.

Run 145 exported **512×257** maps from **128×65** moisture physics with bilinear reconstruction. It used 72 forcing updates/year (roughly five-day steps), simpler directionally split transport and one hydrology spin-up of three years. It did not have today's repeated coupled loop. Its old acceptance criterion omitted snow-mass drift; about 10.8% snow-mass drift was still reported as converged. Land rainfall correlation was about 0.476 and land below 1 mm/year about 6.17%, but seeds, inputs and grid conventions differ from 167/171. Exact historical wall time was not recovered.

Source inspection is read-only: for example, `git show 10b4aba:climate_ocean.cpp`. Do not reset the current working tree. The major rewrites are around `5ae36fa` (146) and `23a4bff` (147). Run 160's 85 hydrology years and two-hour runtime came from later nested coupling/spin-up, not the earlier 145 workload. The newer baseline/FFSL work already removed that performance trap.

## Suggested work sequence

1. Read the saved diagnostics and current implementation. Reproduce the reported numerical comparisons without a new climate run. Confirm sign, units, masks, temporal aggregation and native/output grids.
2. If delegating, use independent ownership: atmospheric/vertical-moisture coupling; ocean/SST budgets; numerical/anisotropy/reference or historical comparison. The coordinator should choose and integrate one coherent revision. Agents must not launch concurrent global benchmarks: the registry allocates IDs on completion and is not concurrency-safe.
3. Build targeted component experiments that can distinguish the hypotheses. Consider an offline diagnostic with observed boundary inputs to test the surface momentum mechanism independently of our generated 850 hPa state; label this strictly as validation and never feed those observations into production. Also audit fixed zonal forcing, spectral filtering and spherical metric/transport choices when examining latitude banding; the appearance alone does not identify the culprit.
4. Implement the smallest physically coherent correction supported by the evidence. Do not independently tune two compensating errors until they happen to recover Köppen. Preserve field reconstruction and fast baseline transport unless a targeted check demonstrates a defect.
5. Run focused deterministic tests, then one bounded 128×64 integration if needed. Only use a 512×256 export comparison after the component evidence supports it; its coarse component resolution remains explicit. Use a fixed seed and matching workload. Stop an unproductive run rather than extending it for hours.
6. Compare rainfall, column moisture, wind speed/direction/divergence and regional structure, SST gradients, ocean currents and temperature/classification behavior separately. Keep a working default; promote a replacement only after explaining the evidence and remaining tradeoffs. Finish with one recommended path, exact validation commands, runtime and unresolved uncertainty.

## Commands

PowerShell, repository root. The directory name says `x64-Debug`, but these builds select **Release**.

```powershell
cmake --build out/build/x64-Debug --config Release --target UndiscoveredWorlds climate_boundary_layer_tests climate_atmosphere_tests climate_context_tests climate_grid_tests climate_transport_tests climate_hydrology_tests climate_energy_tests climate_flow_tests climate_ocean_dynamics_tests --parallel 4

ctest --test-dir out/build/x64-Debug -C Release -R '^(climate_boundary_layer_tests|climate_atmosphere_tests|climate_context_tests|climate_grid_tests|climate_transport_tests|climate_hydrology_tests|climate_energy_tests|climate_flow_tests|climate_ocean_dynamics_tests)$' --output-on-failure

uv run --offline --with numpy --with pillow python -m unittest discover -s tests/scripts -p test_climate_baseline_checks.py
```

Read/check the existing trial, without rerunning climate:

```powershell
uv run --offline scripts/benchmarks/check-climate-baseline.py --run-id 171 --compare-run-id 167
uv run --offline python scripts/benchmarks/write-benchmark-report.py --run-id 171 --elapsed-seconds 81.6588620 --compare-run-id 167
```

The comparison checker intentionally exits **2** for 171. It flags any listed metric deterioration beyond 1e−6 and reports eight regressed targets. This strict numerical comparison is not a statistically calibrated rule for judging all future model tradeoffs. Do not weaken it just to conceal the large rainfall regression. Without `--compare-run-id`, the checker checks numerical validity and records reference scores without a regression verdict. It rewrites the check JSON; preserve any additional manual verification annotations you still need.

Bounded diagnostic, after building and only when needed:

```powershell
./scripts/benchmarks/run-climate-benchmark.ps1 -Resolution 128 -ClimateResolution 64 -CouplingIterations 2 -Seed 20260906 -TropicalClosure mixed-layer
```

For a controlled eight-update comparison, change the limit to 8 and compare with an equivalent `-TropicalClosure legacy` control. For a later 512×256 gallery, use:

```powershell
./scripts/benchmarks/run-climate-benchmark.ps1 -Resolution 512 -ClimateResolution 64 -CouplingIterations 8 -Seed 20260906 -TropicalClosure mixed-layer
```

The wrapper prints the new run ID. Check that ID against 167 and regenerate its gallery using its measured elapsed seconds. Keep `-ClimateMode baseline`, `-Transport ffsl` and the default calibrated temperature mode fixed unless deliberately isolating one of them. The reference inputs come from `configs/app.env`; environment overrides may change those inputs.

## Verification already completed / remaining

The retained implementation passed the Release build, **9 targeted climate CTests**, **3 regression-gate unit tests**, Python compilation, PowerShell wrapper parsing, rejection of an invalid CLI closure value and a targeted whitespace check. Run 171 has eleven decoded 512×256 PNGs, finite numerical fields and closed global water budget. Seventy-two gallery links resolve. The rainfall, July SST and July surface LIC were visually inspected.

Unverified: full coupled equilibrium, native 512-wide atmosphere/ocean performance, transient-mode integration with the new closure, and a fresh legacy benchmark after restoring defaults. Numerical tests are not proof that this reduced atmospheric closure or its climate output is correct.

## Done when

Deliver a supported diagnosis and an implemented, tested improvement that addresses a demonstrated coupling/forcing defect within the runtime budget. Report the reference comparisons and material tradeoffs honestly. If a candidate fails, retain the better working default and explain the concrete remaining obstacle; do not present experimental code as a validated upgrade. Keep only useful evidence and give exact reproduction commands.

The human will inspect rainfall belts and dry interiors, ITCZ shape and geographic wind deflection, LIC/particle flow, SST/ocean-current patterns and final Köppen maps in the minimal comparison gallery. Numerical conservation, physical consistency, matching resolution and acceptable runtime remain the agent's responsibility.
