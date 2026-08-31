# Goal

Add a deterministic, conservative two-layer rainfall model using the low-cost physics that transfer well from ExoPlaSim.

# Final state

- Recorded Earth benchmarks accept explicit imported dimensions; iterations now use 512×257 inputs.
- Tropical background moisture mixing and a two-pass, area-conservative convergence footprint reduce grid-scale filaments.
- Stratiform conversion begins at 88% relative humidity for broader, softer rainfall tails.
- Free-tropospheric vapor uses 65% of the resolved upper-wind departure, retaining more broad-scale structure than run 101.
- Boundary-layer transport slows by 15% over land, providing a low-cost coast-crossing convergence proxy without reducing ocean evaporation.
- Seasonal forcing still migrates the ITCZ; run 106 places it at 5.625°S in January and 9.141°N in July.

# Verification

- `cmake --build --preset build-x64-release`
- `ctest --test-dir out/build/x64-Release --output-on-failure -C Release` — 7/7 passed.
- Runs 102–106, seed 1, 512×257 — recorded with images and diagnostics.
- Run 106 — hydrology 11.37 s over four converged cycles; maximum absolute relative seasonal area-weighted water residual below 2×10⁻⁹.
- `extra/climate.xlsx` — IDs 102–106 present with `HRES=512`; run 106 total 44,872 and cached WRE 0.5707906; formula-error search returned no matches; all four sheets rendered and inspected.

# Remaining uncertainty

- The model does not resolve a diurnal land-breeze/cold-pool cycle or offshore MCS propagation.
- Upper winds transport vapor, not a post-condensation cloud/ice reservoir, so offshore anvil blowback is only represented indirectly.
- Run 106 is a measured compromise: better land coverage and concentration than run 101, but some flow-aligned patchiness remains and the 512px reference correlation is only 0.1188.
