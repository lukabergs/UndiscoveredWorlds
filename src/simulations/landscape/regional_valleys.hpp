#pragma once

// The incision of the pinned channels by the world's uplift and the valley floors
// (docs/architecture/erosion-hydrology-plan.md, Phase R2: the R6 contract, which R9
// implements; the diagnosis is in section 6, 'The user's review of regional runs 5 to 8',
// items 4 and 4b). Every value is a pure function of the global fields, the seed and
// absolute position. Called by regional_pass.cpp. Planet-free.
//
// The frame. The global base elevation is a cell's mean elevation, and the world's own
// erosion run already balanced its uplift against its erosion at the global scale, so the
// regional pass adds no uplift to the mean: in the frame of the mean, uplift is the rivers
// cutting down (the incision) while the land between them keeps the start surface, and the
// characteristic solution's uplift term shapes the flanks within the knickpoint wave's
// reach (analytical_relief.hpp, uplift_m_myr). The block means are not restored: the mean
// drops by about the channel and flank share of the cells times the incision (R9 measures).

#include <cstdint>
#include <vector>

#include "simulations/landscape/channel_network.hpp"
#include "simulations/landscape/regional_grid.hpp"
#include "simulations/landscape/regional_types.hpp"

namespace landscape
{
// The incision per global cell (W x H). `active` false: no incision (the network's filled
// elevations as they are; the vectors stay empty).
struct GlobalIncision
{
    bool active = false;
    std::vector<float> incision_m;          // I' >= 0 per global cell, 0 on sea cells
    std::vector<float> filledelevation_m;   // the network's filled elevation less I'
};

// R9. Per global land cell the rate is the larger of valleys.upliftscale x max(0, the rock
// uplift rate) and valleys.denudationratescale x denudation / erosionduration_myr (0 when
// the duration is 0 or a field is empty), the depth I = rate x relief.elapsed_myr, capped
// at the cell's filled elevation minus the sea level (a channel never incises below the
// sea level) and at 0 below; then made consistent downstream over the global network,
// upstream cells first: I'(r) = max(I(r), max over the donors u of I'(u) - (F(u) - F(r))),
// so F' = F - I' never rises downstream and the bed law's order b(receiver) <= F(receiver)
// <= b(upstream) still holds. A pure function of the global fields alone, so every window
// agrees. active when a scale is positive and its field present and some I' is positive.
GlobalIncision globalincision(const RegionalGlobalInputs& inputs, const RegionalParameters& parameters);

// The valley floors (R9, after the legacy expandrivers and pasterivercircle): no-op while
// valleys.floorcells <= 0. Contract: on the domain (domainfinewidth() x domainfineheight())
// after the solve, for every land channel cell (channels.channel 1, not sea) of discharge Q
// and bed b: the half-width w = min(valleys.floormaxcells, valleys.floorcells x sqrt(Q /
// valleys.reference_m3s)) fine cells; every land cell that is not a channel within
// Chebyshev distance d <= w of it is lowered to at most b, and beyond, within
// valleys.reachcells, to at most b + (d - w) x valleys.sideslope x the fine cell's size
// (sideslope 0: no sides); a cell shaped by several channels takes the lowest target; no
// cell is lowered to or below the sea level (the next float above it at least). reachcells
// is at most k (one global cell), which the halo rule (halo >= 2 x meanpasses + 2) covers:
// the channels agree between windows from 2 x meanpasses + 1 global cells inside the
// domain's edge, the kept window starts one cell further in. R9 keeps the overlap tests
// bit-identical with the valleys on.
void carvevalleys(const RegionalGlobalInputs& inputs, const RegionalWindow& window, const RegionalParameters& parameters,
                  const FineChannels& channels, const std::vector<uint8_t>& sea, std::vector<float>& elevation_m);
}
