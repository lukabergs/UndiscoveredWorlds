#pragma once

// Lakes and present-day ice of the regional pass (docs/architecture/erosion-hydrology-plan.md,
// Phase E, step E3, and Phase R: R0 moved the rules here from regional_pass.cpp
// unchanged; R1 fixes the channel cells inside lakes; R4 owns the shores and the ice
// edges). Local rules only (no flood fill), each a pure function of the global
// fields, the seed and the fine cell's absolute position. Planet-free.

#include <cstdint>
#include <vector>

#include "simulations/landscape/channel_network.hpp"
#include "simulations/landscape/regional_grid.hpp"
#include "simulations/landscape/regional_types.hpp"
#include "simulations/landscape/upsampling.hpp"

namespace landscape
{
// Domain fields (domainfinewidth() x domainfineheight(), row-major).
struct FineWater
{
    std::vector<float> lakesurface_m;   // 0 where no lake water
    std::vector<float> icethickness_m;
};

// Phase E's rules (kept while shores.shorenoise is 0): lake water stands on a fine
// cell that is not sea, whose global cell carries a lake surface and whose
// `elevation_m` lies strictly below it; ice is the upsampled thickness, at least 0,
// on fine cells whose global cell or one of its 8 neighbours has ice. R4: the lake
// and the ice are where a signed per-global-cell indicator (lake or ice cells
// negative, a hashed magnitude), interpolated at the warped positions plus fractal
// noise, lies below zero within those cells and their neighbours; lake beds are
// lowered below the surface there (so `elevation_m` may change), and channel cells
// inside a lake carry its water.
FineWater regionalwater(const RegionalGlobalInputs& inputs, const RegionalWindow& window, const FineWarp& warp,
                        const RegionalParameters& parameters, const std::vector<uint8_t>& sea,
                        const FineChannels& channels, std::vector<float>& elevation_m);
}
