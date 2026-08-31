# Final state

The ExoPlaSim-inspired rainfall checkpoint is complete and represented by deterministic benchmark runs 117 and 118 at 512 px.

# Included

- 128x65 internal moisture grid with sub-grid land fraction and relief retained during aggregation.
- Two-layer vapor transport, Kuo-style convection, buoyancy closure, elevated moisture accession, shallow/dry mixing, falling-rain re-evaporation, snow storage, and soil-moisture-limited evapotranspiration.
- Deterministic synoptic and coastal day/night phases, seasonal migration, transient eddy mixing, cloud/latent-heating memory, and conservative scale-selective damping.
- Conservative adjacent-row meridional fluxes; multi-row jumps can no longer bypass a zero-transport polar edge.
- Benchmark runs 107-118 recorded in `extra/climate.xlsx` with `HRES=512`; generated images and GeoTIFFs remain under `extra/img/earth/benchmark/`.

# Verification

- `cmake --build out/build/x64-Debug --config Debug`
- `ctest --test-dir out/build/x64-Debug -C Debug --output-on-failure` (7/7 passed)
- `cmake-vs --build --preset build-x64-release`
- Run 118: 16.3 s total, 1.46 s hydrology, maximum seasonal relative water-budget residual `2.49e-9`.
- Runs 117 and 118 produced identical climate counts and `WRE=0.3679873165`.
- Polar diagnostic: cap precipitation fell from 1958/1760 to 0.58/0.23 mm/month north/south after the flux correction.
- Workbook rendered on all four sheets; no formula error values remain.

# Remaining uncertainty

The reference workbook measures class-frequency agreement, not spatial rainfall-field agreement. Further tuning should preserve this checkpoint and evaluate georeferenced precipitation structure directly.
