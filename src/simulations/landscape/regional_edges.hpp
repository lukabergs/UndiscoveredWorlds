#pragma once

// The sea gate and the straight-run breaker of the regional coast (docs/architecture/erosion-hydrology-plan.md,
// Phase R2: the R6 contract, which R8 implements; the diagnosis is in section 6, 'The
// user's review of regional runs 5 to 8', item 2). Every value is a pure function of
// the global fields, the seed and the fine cell's absolute position, so two windows agree
// on their overlap. Called by regional_surface.cpp on the land-share rule's widened box.
// Planet-free.

#include <cstdint>
#include <vector>

#include "simulations/landscape/regional_grid.hpp"
#include "simulations/landscape/regional_types.hpp"
#include "simulations/landscape/upsampling.hpp"

namespace landscape
{
// Phase E's eligibility of a global cell's fine cells for the sea (the 3 x 3 gate): 1
// where a global sea cell lies among the cell and its 8 neighbours and the cell is not a
// lake cell (lake surface above zero) unless the cell itself is sea. W x H.
std::vector<uint8_t> nearseacells(const RegionalGlobalInputs& inputs);

// The global cells that may hold no fine sea (nearseacells 0) as 1.0, the others as 0.0:
// the decision surface is lifted towards them (regional_surface.hpp). W x H.
std::vector<float> inlandcells(const RegionalGlobalInputs& inputs);

// The gate: one value per fine cell of `box`'s domain (domainfinewidth() x
// domainfineheight(), row-major), 1 where the cell may be sea. `warp` is the box's warp
// (regionalwarp; empty without one).
// - edges.smoothgate 0 (Phase E to R): nearseacells of the cell's own global cell, so the
//   gate is a step function of the global grid and cuts the sea on a block edge wherever a
//   global land cell with no sea in its 3 x 3 holds fine cells below the sea level (regional
//   run 6, global cell (184, 49): a straight coast of 29 fine cells on the 183/184 edge).
// - edges.smoothgate 1 (R8): a gate that is not a step function of the block: the global sea
//   indicator (1 on sea cells, 0 elsewhere, lake cells excluded as above) interpolated by
//   `bicubic` at the warped fine position is positive, so the fine sea may reach into a land
//   cell within the interpolation's support (two global cells) wherever the decision surface
//   lies below the sea level, and the land-share rule still keeps a land cell's block at
//   least `landshare` land. R8 states the rule it implements here and in the plan.
std::vector<uint8_t> seagate(const RegionalGlobalInputs& inputs, const RegionalWindow& box, const FineWarp& warp,
                             const RegionalParameters& parameters);

// The straight-run breaker (R8, after the legacy removestraights): no-op while
// edges.straightmax <= 0. Contract: `mask` is the box's mask (1 = sea; domainfinewidth() x
// domainfineheight()); a straight shore run is a row or column of consecutive sea cells
// whose cells all have land on the same side (north, south, east or west); a run longer
// than edges.straightmax fine cells gets chunks of edges.chunkmin to edges.chunkmax cells
// flipped, drawn from a hash of the seed and the run's absolute position (so every window
// draws the same), each chunk inside the global cell it starts in: land in a global land
// cell's block, land or sea in a sea cell's (the legacy: raise where the tile is land).
// Only cells whose `gate` is 1 may become sea, and no cell of a global lake cell. The
// caller repeats the pinch and hole passes after it and counts the land shares on the
// result, so the land-share rule still holds.
void breakstraightruns(const RegionalGlobalInputs& inputs, const RegionalWindow& box, const RegionalParameters& parameters,
                       const std::vector<uint8_t>& gate, std::vector<uint8_t>& mask);
}
