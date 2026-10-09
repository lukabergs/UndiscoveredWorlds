#pragma once

// The coast's ruggedness from the world's fields and the island sink (docs/architecture/erosion-hydrology-plan.md,
// Phase R2: the R6 contract, which R7 implements and tunes; the diagnosis is in section 6,
// 'The user's review of regional runs 5 to 8', item 1). Every value is a pure function of
// the global fields, the seed and the fine cell's absolute position. Called by
// regional_surface.cpp (the Surfacer and the land-share rule's widened box). Planet-free.

#include <cstdint>
#include <vector>

#include "simulations/landscape/regional_grid.hpp"
#include "simulations/landscape/regional_types.hpp"

namespace landscape
{
// Per-global-cell factors of the coast noise (W x H). `active` false means amplitude_m
// 0 and dropfactor 1 everywhere, and the Surfacer reads neither (so the R2 results stay
// bit-identical: an interpolation of a constant field is not exactly the constant).
struct CoastFactors
{
    bool active = false;
    std::vector<float> amplitude_m;   // added to the coast noise's amplitude, interpolated at the sampling position
    std::vector<float> dropfactor;    // multiplies the sea-side push, 0 .. 1, interpolated the same way
};

// R7 (R6 implements the formula, R7 tunes it): amplitude_m = coast.glacialamplitude_m x
// glaciation + coast.denudationamplitude_m x min(1, denudation / coast.denudationreference_m);
// dropfactor = clamp(1 - coast.glacialdropfraction x glaciation, 0, 1); an empty input field
// counts as 0. active when a coefficient is positive and its field is present.
CoastFactors coastfactors(const RegionalGlobalInputs& inputs, const RegionalParameters& parameters);

// The island sink (R7, after the legacy removesmallislands): no-op while coast.islandmargin
// <= 0. Contract: `mask` is the box's mask (1 = sea; domainfinewidth() x domainfineheight()
// of `box`), `gate` the box's sea gate (regional_edges.hpp). For every global sea cell of
// the box whose 3 x 3 neighbourhood holds no one-cell island and no mountain island of the
// world (inputs.island, inputs.mountainisland; empty fields: none), when the ring of fine
// cells at m = max(1, round(coast.islandmargin x k)) cells around its k x k block holds no
// land cell (a ring cell beyond the box or a pole counts as land, so nothing sinks there;
// the box's margin keeps that away from the domain), every land cell of the block becomes
// sea where its gate is 1. The elevation's clamp and shore floor follow the mask. A pure
// function of the mask within m cells of the block and the global fields, so the rule's
// reach stays within a share pass's (regional_surface.cpp, landsharesurface).
void sinksmallislands(const RegionalGlobalInputs& inputs, const RegionalWindow& box, const RegionalParameters& parameters,
                      const std::vector<uint8_t>& gate, std::vector<uint8_t>& mask);
}
