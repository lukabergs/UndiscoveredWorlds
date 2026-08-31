# Goal

Add a deterministic, conservative two-layer rainfall and land-water model using the low-cost physics that transfer well from ExoPlaSim.

# Status

Completed on `climate-physics_4_hydrology_refactor_260831/two-layer-moisture`.

- Boundary/free-troposphere moisture, vertical exchange, moist adjustment, diagnostic cloud fraction, falling-precipitation re-evaporation, rain/snow phase, finite snow storage, melt, infiltration, and runoff are integrated.
- Hydrology retains floating-point winds and rainfall while preserving the existing serialized short rainfall maps.
- Surface exchange accounts for roughness, stability, pressure, and gravity; hot grids are contiguous and reused.
- The benchmark writer now records AG total and AH horizontal resolution while preserving the AI error formula.
- Full-resolution run 101 completed and recorded all images and diagnostics.

# Verification

- `cmake --build --preset build-x64-release`
- `ctest --test-dir out/build/x64-Release --output-on-failure -C Release` — 7/7 passed.
- Factor-4 Earth smoke benchmark, seed 4101, 512×257 — completed in 24.5 s; hydrology 7.84 s over three converged cycles; maximum relative seasonal area-weighted water residual below 3×10⁻⁹ in the final run.
- Full Earth benchmark 101, seed 1, 2048×1025 — hydrology 346.2 s over four converged cycles; maximum relative seasonal area-weighted water residual below 1.3×10⁻⁹.
- `extra/climate.xlsx` inspection and visual render — ID 101, total 737354, HRES 2048, WRE 0.6306901, zero formula errors.

# Remaining uncertainty

- Clouds remain diagnostic and do not feed the radiation solve; moist-adjustment heating does not feed back into the global circulation.
- Run 101 improves land dry-cell and concentration diagnostics but worsens precipitation placement: area-weighted land correlation 0.0815 and 5-degree correlation 0.1378.
