#pragma once

// Glacial troughs: the U-shaped valleys of the regional map (docs/architecture/
// erosion-hydrology-plan.md, Phase E5: the E5.0 contract, which step E5.1 implements;
// section 6, 'Decided for Phase E5'). Called by regional_pass.cpp after the valley floors
// (regional_valleys.hpp) and before the cirques. Planet-free.
//
// Contract (E5.1). No-op while troughs.width_km <= 0 or the glacial fields are inactive. On
// the domain (domainfinewidth() x domainfineheight()) after the solve and the valley
// floors: every land channel cell (channels.channel 1, not sea) whose fine glaciation
// (glacial.glaciation) is at least troughs.glaciationmin and whose global cell's ice
// catchment A_i (glacialcatchment_km2 at channels.globalcell, the upstream global cell of
// the link; 0 or absent: no trough) is positive is the axis of a trough of half-width
// W = min(troughs.widthmax_km, troughs.width_km x (A_i / troughs.reference_km2)^widthexponent)
// km on the ground, in fine cells W / the fine cell's size (finepatchmetric's dy_m, the
// same in every window) and at most k (one global cell, the halo rule as carvevalleys):
// every land cell that is not a channel within Euclidean distance d <= W (fine cells) of
// the axis is lowered to at most b + troughs.wall_m x (d / W)^exponent, b the axis cell's
// bed (its pinned elevation), so the floor is flat near the axis and the walls steepen
// towards the half-width (exponent 2: a parabola, the glacial trough's cross-section);
// a cell within several troughs takes the lowest target; no cell is lowered to or below
// the sea level (the next float above it at least); nothing is raised. A trough that lies
// below the land it would cut leaves it (min), so the trough meets the land without a
// step. Widths in kilometres, never fine cells (R9's lesson). The bed is read before any
// trough is cut (a channel cell is never changed), so the order of the axes does not
// matter. E5.1 keeps the overlap tests bit-identical with the troughs on and adds its
// variants; it measures cross-sections on erosion run 9's window B (k 16, where 48 percent
// of the land was glaciated) and on one scratch window at k 64, where the half-width
// spans several cells, and sets the defaults from them.

#include <cstdint>
#include <vector>

#include "simulations/landscape/channel_network.hpp"
#include "simulations/landscape/regional_glacial.hpp"
#include "simulations/landscape/regional_grid.hpp"
#include "simulations/landscape/regional_types.hpp"

namespace landscape
{
void carvetroughs(const RegionalGlobalInputs& inputs, const RegionalWindow& window, const RegionalParameters& parameters,
                  const FineGlacial& glacial, const FineChannels& channels, const std::vector<uint8_t>& sea,
                  std::vector<float>& elevation_m);
}
