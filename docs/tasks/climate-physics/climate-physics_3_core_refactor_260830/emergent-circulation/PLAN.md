# Goal

Replace prescribed pressure belts and arbitrary moisture injection with a coupled two-layer circulation and conservative bulk-flux water cycle.

# Phase

3 — Complete.

# Final state

- Surface pressure, lower/upper winds, divergence, and vertical motion converge together on physical grid spacing.
- Coriolis, Rayleigh drag, hypsometric thickness, mass correction, and bulk-aerodynamic evaporation use explicit units.
- Gaussian pressure belts, latitude-based calm-wind direction, asymmetric rainfall smoothing, and SST-scaled moisture injection are inactive.
- WorldClim validates monthly temperature, precipitation, and 10 m wind; IMERG validates annual global/land/ocean precipitation.
- Run 42 is the retained benchmark: global IMERG precipitation 946.36 versus 994.78 mm/year; water residual below `6e-9` per season.

# Verification

- `cmake --build out/build/x64-Release --config Release --parallel`
- `ctest --test-dir out/build/x64-Release -C Release --output-on-failure`
- Release Earth benchmark run 42.
- `extra/validation/seed_1/climate_atmosphere_budget.csv`
- `extra/validation/seed_1/annual_imerg_precipitation_comparison.csv`

# Open questions

- Land precipitation remains low: 372.40 versus 751.29 mm/year area-weighted.
- WorldClim land wind remains too strong and spatially anticorrelated: 6.41 versus 3.17 m/s.
- Cloud-water conversion and inter-season storage belong in phase 4.
