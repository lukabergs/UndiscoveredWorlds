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
    std::vector<float> rate_m_myr;          // the incision rate per global cell (m/Myr, >= 0, 0 on sea cells),
                                            // the solve's uplift term after upsampling (regional_pass.cpp)
};

// R9. Per global land cell the rate is the larger of valleys.upliftscale x max(0, the rock
// uplift rate) and valleys.denudationratescale x denudation / erosionduration_myr (0 when
// the duration is 0 or a field is empty), the depth I = rate x relief.elapsed_myr, capped
// at the cell's filled elevation minus s+, the next float above the sea level (a land
// channel never incises to the sea level: the bed law's lower bound), at the cell's local
// global relief (the largest absolute difference between its base elevation and its 8
// neighbours', as regional_surface.cpp's localrelief) and at 0 below. The relief cap: the
// frame of the mean (the true characteristic solution less U t) holds while t stays below
// the plateau's response time; beyond it the whole block would drop without bound, whereas
// the relief saturates, so elapsed_myr must stay below the fine response times (36 to 81
// Myr for a first-order cell on erosion run 9's windows) and the cap bounds the rest. Then
// made consistent downstream over the global network, upstream cells first: I'(r) =
// max(I(r), max over the donors u of I'(u) - (F(u) - F(r))), so F' = F - I' never rises
// downstream and the bed law's order b(receiver) <= F(receiver) <= b(upstream) still holds.
// A pure function of the global fields alone, so every window agrees. active when a scale
// is positive and its field present and some I' is positive. Computed on the incised
// surface in float, F'(r) = max(s+, min(F(r) - I(r), min over the donors of F'(u)))
// (F'(r) = F(r) where F(r) < s+), which is the rule above in exact arithmetic and keeps
// F'(r) <= F'(u) exactly after rounding; I' = F - F'. Receivers as channel_network.cpp
// decodes them (a direction beyond a pole: none); cells on a receiver loop (bad input,
// which rasterisechannels rejects near a window) take their own capped I.
GlobalIncision globalincision(const RegionalGlobalInputs& inputs, const RegionalParameters& parameters);

// The valley floors (R9, after the legacy expandrivers and pasterivercircle): no-op while
// valleys.floor_km <= 0. Contract: on the domain (domainfinewidth() x domainfineheight())
// after the solve, for every land channel cell (channels.channel 1, not sea) of discharge Q
// and bed b: the half-width w = min(valleys.floormax_km, valleys.floor_km x sqrt(Q /
// valleys.reference_m3s)) km on the ground, in fine cells w / the fine cell's size (so the
// floor's ground width does not change with k; at k 4 a few km is under one cell, no floor);
// every land cell that is not a channel within
// Chebyshev distance d <= w of it is lowered to at most b, and beyond, within
// valleys.reachcells, to at most b + (d - w) x valleys.sideslope x the fine cell's size
// (sideslope 0: no sides); a cell shaped by several channels takes the lowest target; no
// cell is lowered to or below the sea level (the next float above it at least). The fine
// cell's size is the north-south spacing pi R / (k H) (finepatchmetric's dy_m), the same in
// every window; b is the elevation at the channel cell (its pinned bed). reachcells and the
// floor's extent are clamped to k (one global cell), which the halo rule (halo >= 2 x meanpasses + 2) covers:
// the channels agree between windows from 2 x meanpasses + 1 global cells inside the
// domain's edge, the kept window starts one cell further in. R9 keeps the overlap tests
// bit-identical with the valleys on.
// Multiscale plan S2.3 (the refinement loop): with bandupper > 0 a form is carved only where
// its full width in fine cells, max(floor half-width, the sides' reach min(k, reachcells)), is
// at least bandlower (detail_cells / 2: resolved) and its floor half-width below bandupper
// (detail_cells x factor / 2: under one parent cell), so a sub-cell floor with sides is carved
// at every level whose floor stays under one parent cell (each a minimum: it only refines);
// bandupper 0 (the single pass) keeps the rule above.
void carvevalleys(const RegionalGlobalInputs& inputs, const RegionalWindow& window, const RegionalParameters& parameters,
                  const FineChannels& channels, const std::vector<uint8_t>& sea, std::vector<float>& elevation_m,
                  double bandlower = 0.0, double bandupper = 0.0);
}
