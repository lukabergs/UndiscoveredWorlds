# Goal

Add a deterministic, conservative two-layer rainfall and land-water model using the low-cost physics that transfer well from ExoPlaSim.

# Status

Completed on `climate-physics_4_hydrology_refactor_260831/two-layer-moisture`.

- Boundary/free-troposphere moisture, vertical exchange, moist adjustment, cloud fraction, falling-precipitation re-evaporation, rain/snow phase, finite snow storage, melt, infiltration, and runoff are integrated.
- Hydrology retains floating-point winds and rainfall while preserving the existing serialized short rainfall maps.
- Surface exchange accounts for roughness, stability, pressure, and gravity.
- Hot grids are contiguous and reused; imported terrain accepts explicit reduced dimensions.
- `scripts/prepare-reduced-earth-benchmark.py` creates mutually exclusive reduced Earth land/sea maps and a matching precipitation reference.

# Verification

- `cmake --build --preset build-x64-release`
- `ctest --test-dir out/build/x64-Release --output-on-failure -C Release` — 7/7 passed.
- Factor-4 Earth smoke benchmark, seed 4101, 512×257 — completed in 24.5 s; hydrology 7.84 s over three converged cycles; maximum relative seasonal area-weighted water residual below 3×10⁻⁹ in the final run.

# Remaining uncertainty

- A final 2048×1025 Earth benchmark was not run; runtime and reference skill at full resolution remain unmeasured.
- The reduced benchmark is useful for conservation and failure detection, not tuning: land area-weighted precipitation correlation was -0.146 and tropical precipitation remained too dry.
