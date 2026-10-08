#pragma once

// The regional pass (docs/architecture/erosion-hydrology-plan.md, Phase E, step E3,
// and section 6, "Decided for Phase E"): global fields and a window in, the window's
// fine fields out. Steps: the start surface (the global base elevation upsampled,
// band-limited noise added with an amplitude from the local relief and tapered near
// the sea level, block means corrected); the pinned channels; the analytical relief
// on the domain; sea, lakes, rivers and climate on the fine grid by local rules;
// the halo discarded (Phase R: the start surface and the sea in regional_surface.hpp,
// lakes and ice in regional_water.hpp). Every step is a pure function of the global fields, the seed
// and absolute position within the halo's reach, so two windows give the same
// values on their overlap (test Landscape.RegionalOverlap). Planet-free: the planet
// adapter is src/simulations/regional/regional_inputs.hpp.

#include <cstdint>
#include <vector>

#include "simulations/landscape/regional_grid.hpp"
#include "simulations/landscape/regional_types.hpp"

namespace landscape
{
// The fine fields of the kept window (finewidth() x fineheight(), row-major, north
// row first); the bundle writer stores them with the window.
struct RegionalFields
{
    RegionalWindow window;
    std::vector<float> elevation_m;
    std::vector<uint8_t> sea;
    std::vector<uint8_t> riverdirection;     // 1..8 on channel cells, 0 elsewhere
    std::vector<float> discharge_m3s;        // pinned channels: the global discharge; elsewhere the fine accumulation
    std::vector<int32_t> januaryflow, julyflow;  // legacy units, on pinned channels
    std::vector<float> lakesurface_m;        // 0 where no lake water
    std::vector<float> januarytemperature_c, julytemperature_c;
    std::vector<float> januaryprecipitation_mm, julyprecipitation_mm;
    std::vector<float> icethickness_m;
    std::vector<int16_t> climate;

    // Diagnostics (plan E3; window-dependent, not fields): the largest |block mean -
    // global value| over the kept window's global cells of the start surface after the
    // block-mean correction, before the channels are pinned and the solve; and the
    // solve's AnalyticalReliefResult::converged (the last two passes' receivers equal)
    // and changedreceivers (over the whole domain); the mean and largest lowering by the
    // solve (start surface minus result) over the kept window's unpinned land cells.
    double blockmeanerror_m = 0.0;
    bool reliefconverged = false;
    int changedreceivers = 0;
    double meanlowering_m = 0.0;
    double maxlowering_m = 0.0;
};

// Throws std::invalid_argument when the window is not on the inputs' grid, a global
// field has the wrong size, the radius is not positive, the sea level not finite, or
// window.halo < 2 x meanpasses + 2 (the block-mean correction reaches 2 global cells
// per pass; two more keep the channels and the solve's edge out of the kept window);
// rasterisechannels and solveanalyticalrelief throw on their own bad inputs.
RegionalFields runregionalpass(const RegionalGlobalInputs& inputs, const RegionalWindow& window,
                               const RegionalParameters& parameters);
}
