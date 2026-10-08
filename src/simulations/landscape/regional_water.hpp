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
// on fine cells whose global cell or one of its 8 neighbours has ice. R1: a channel
// cell of a lake cell carries the lake where its bed lies at or below the higher of
// the lake surface and the global cell's filled elevation (a closed lake can lie below
// its spill level, where the bed is pinned), and its `elevation_m` is lowered to at
// most the surface, so no lake surface lies below its cell; such a lake bed may lie
// below its outlet's sill, so the elevation along the channel directions rises from
// the lake to the sill (by up to 56 m in window C of erosion run 9), and a lowered
// cell lies below the sea level where its lake's surface does. R4: the lake
// and the ice are where a signed per-global-cell indicator (lake or ice cells
// negative, a hashed magnitude), interpolated at the warped positions plus fractal
// noise, lies below zero within those cells and their neighbours; lake beds are
// lowered below the surface there (so `elevation_m` may change); channel cells
// inside a lake keep carrying its water.
FineWater regionalwater(const RegionalGlobalInputs& inputs, const RegionalWindow& window, const FineWarp& warp,
                        const RegionalParameters& parameters, const std::vector<uint8_t>& sea,
                        const FineChannels& channels, std::vector<float>& elevation_m);
}
