# Goal

Carry atmospheric and soil water through the seasonal cycle and spin the hydrology to a stable periodic state before recording climate outputs.

# Phase

4 — Branch complete.

# Scope

- Persist atmospheric moisture and soil storage across representative months.
- Iterate complete annual cycles to a deterministic convergence criterion.
- Preserve the closed water budget and expose spin-up diagnostics.
- Restore EF in benchmark images and avoid rewriting the fixed reference image.
- Run and compare one retained Release Earth benchmark.

# Constraints

- Earth reference grids remain validation-only.
- Do not change Köppen thresholds or apply a global rainfall multiplier.
- Keep fixed-seed generation deterministic.

# References

- `climate_ocean.cpp`
- `climate_physics.hpp`
- `climate_validation.cpp`
- `generation_tuning.hpp`
- `extra/climate.xlsx`

# Final state

- Atmospheric and soil water persist across the four representative seasons.
- The solver stops after storage changes by at most `0.5%`, bounded to four deterministic cycles.
- Per-season budgets include initial storage and retain relative residuals below `5.1e-9`.
- EF is visible in both reference and simulated benchmark images; `0.png` is regenerated only when missing or outdated.
- Run 43 converged on cycle 4: land precipitation rose from 372.40 to 414.11 mm/year; BWh fell by 1,911 cells and BWk by 3,165.

# Verification

- Release build and CTest.
- Fixed-seed Earth benchmark run 43.
- Workbook formula/error and visual checks.

# Open questions

- Persistence improved land rainfall totals, but IMERG land MAE rose from 862.64 to 882.60 mm/year and correlation fell from 0.2097 to 0.2063.
- Circulation and moisture-convergence geography are the next bottleneck.
