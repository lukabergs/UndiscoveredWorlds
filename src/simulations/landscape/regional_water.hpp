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
// cell lies below the sea level where its lake's surface does.
// R4 (shores.shorenoise > 0; D = shores.terrain_m, z the solve's elevation):
// - The indicator: -m on lake cells (lake surface above zero, not sea), +m elsewhere, m a
//   hash of the cell's indices and the seed in 0.101 .. 0.900, interpolated (bicubic) at
//   the fully warped fine positions.
// - The zone: a fine cell (not sea) with a lake cell among its global cell and the 8
//   neighbours takes the surface S of the one whose term in that interpolation (stencil
//   weight times magnitude) is the strongest at it (ties: the smaller global index).
// - The sum: indicator + shorenoise x RegionalNoise (a salt of its own, persistence
//   coast.noisepersistence) + clamp((z - S) / D, -w, 1), w the lake cells' mask (1 on lake
//   cells, 0 elsewhere) interpolated at the warped position and clamped to 0 .. 1; water
//   where the sum is below zero and z < S + D. So the terrain decides far from the surface
//   (high ground is land; deep ground is water where lake cells dominate, not in the valley
//   beyond an outlet) and the indicator near it; w is continuous across block edges.
// - Steps: a wet cell beside a wet cell of a lower surface is dry (a dam), so adjacent
//   water carries one surface.
// - Channel cells: R1's rule (wet at or below max(surface, filled), lowered to the surface)
//   is the shorenoise-0 rule; R4 replaces it. A channel cell inside the water above is wet
//   like any cell; a channel cell of a lake cell outside it is wet when z < S + D and an
//   8-neighbour that is not a channel is wet, or when z < S, unless a wet 8-neighbour has
//   another surface; otherwise (a bed at or above the surface away from the water, or the
//   outlet between two lakes) it stays a river. Of two adjacent channel cells this rule
//   wets at different surfaces, the higher stays a river too.
// - Beds: a wet cell with z >= S is lowered to S minus a depth growing from 0 where the sum
//   is 0 to min(the global lake depth (S minus the zone's base elevation), D) where the sum
//   reaches minus the zone's magnitude, and to at least the next float below S; a wet cell
//   with z < S keeps z. So the lowering is below 2 D. With `globaldepth` (the refinement
//   loop, multiscale plan S2.3) the depth is the global lake depth S minus the zone's base
//   elevation uncapped, so a deep global lake's fine bed reaches it; false (the single
//   pass) keeps the cap.
// - Shores: a dry fine cell (not sea, not a channel) at or below the level L, the highest
//   surface of its wet 8-neighbours and of its own global cell if that is a lake cell, and
//   at most D below it, is raised to L + (L - z) (mirrored about L), at least the next float
//   above L, so by at most 2 D (a dry cell of a lake cell with no water beside it is raised
//   to the mirror of its own lake's surface the same way); a cell more than D below L is left:
//   the deep bed of a pit lake outside the indicator's outline, left dry rather than dammed by
//   hundreds of metres.
// Ice: the indicator with the ice cells (thickness above zero) negative, a salt of its own
// and no terrain term, on the fine cells of ice cells and their neighbours; the thickness is
// the thickness upsampled at the warped position, at least 1 m; a fine sea cell carries ice
// only where its global cell is sea and iced.
FineWater regionalwater(const RegionalGlobalInputs& inputs, const RegionalWindow& window, const FineWarp& warp,
                        const RegionalParameters& parameters, const std::vector<uint8_t>& sea,
                        const FineChannels& channels, std::vector<float>& elevation_m, bool globaldepth = false);
}
