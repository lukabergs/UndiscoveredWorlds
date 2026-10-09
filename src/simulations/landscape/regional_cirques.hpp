#pragma once

// Cirques: the bowls at the heads of glaciated high ground (docs/architecture/
// erosion-hydrology-plan.md, Phase E5: the E5.0 contract, which step E5.2 implements;
// section 6, 'Decided for Phase E5'). Called by regional_pass.cpp after the troughs
// (regional_troughs.hpp) and before the lakes and ice (regional_water.hpp). Planet-free.
//
// Contract (E5.2). No-op while cirques.radius_km <= 0 or the glacial fields are inactive.
// Candidates: for every global land cell whose block lies in the domain and whose
// glaciation (the global field, inputs.glaciation; 0 when empty) is at least
// cirques.glaciationmin, cirques.candidates points, the n-th at the fine offset
// (deterministiccombine(seed, salt, i, j, n) reduced to 0..k-1 in each axis, the two axes
// from different bits) inside the block, a pure function of the global cell, the index
// and the seed, so every window places the same cirques. A candidate is a cirque when its
// fine cell is land (not sea), not a channel cell, and its surface z (the elevation as
// passed in, read before any cirque is cut, so the order of the candidates does not
// matter) lies at or above the run's ELA there (glacial.runela_m, the glacial past's) and
// at most cirques.band_m above it (higher ground was plateau ice, not a cirque's
// headwall). Its floor is f = max(runela at the centre, z - cirques.depth_m), at least the
// next float above the sea level. The bowl: R = cirques.radius_km x 1000 / the fine
// cell's size (finepatchmetric's dy_m) in fine cells, at most k (the halo rule as
// carvevalleys; a radius under one cell carves the centre cell alone); every land cell
// that is not a channel within Euclidean distance d <= R of the centre is lowered to at
// most f + (z_cell - f) x (d / R)^2 (its own surface at the rim, the floor at the centre:
// a parabolic bowl, steepest on the upslope side, which is the headwall); a cell inside
// several bowls takes the lowest target; a cell already below f is left; nothing is
// raised; no cell is lowered to or below the sea level. Radii in kilometres, never fine
// cells (R9's lesson). E5.2 keeps the overlap tests bit-identical with the cirques on
// and adds its variants; it counts cirques and measures their depth and width on erosion
// run 9's window B (k 16) and on one scratch window at k 64, where a bowl spans several
// cells, and sets the defaults from them.

#include <cstdint>
#include <vector>

#include "simulations/landscape/channel_network.hpp"
#include "simulations/landscape/regional_glacial.hpp"
#include "simulations/landscape/regional_grid.hpp"
#include "simulations/landscape/regional_types.hpp"

namespace landscape
{
void carvecirques(const RegionalGlobalInputs& inputs, const RegionalWindow& window, const RegionalParameters& parameters,
                  const FineGlacial& glacial, const FineChannels& channels, const std::vector<uint8_t>& sea,
                  std::vector<float>& elevation_m);
}
