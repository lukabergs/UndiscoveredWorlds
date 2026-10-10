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
// Multiscale plan S3.5 (the loop, bandupper > 0): the floor and the sides by Euclidean distance
// on the ground (column offsets scaled by eastwestratio, the column extent columnextent), not
// Chebyshev squares; with
// valleys.reach_km > 0 the sides reach reach_km on the ground (in fine cells over the fine
// north-south spacing, at most k) instead of reachcells, and the full width max(floor
// half-width, reach) decides both ends of the band, so the whole form is carved at the one level
// where it lies in [bandlower, bandupper) (bandupper infinite at level 1: the forms above its
// band, clamped to one global cell); reach_km 0 keeps reachcells and S2.3's test. The single pass
// never reads reach_km. In every mode a channel cell keeps its value (its bed in the single pass,
// its cell mean in the loop: channelcellmeans lowers it with the floor around it).
void carvevalleys(const RegionalGlobalInputs& inputs, const RegionalWindow& window, const RegionalParameters& parameters,
                  const FineChannels& channels, const std::vector<uint8_t>& sea, std::vector<float>& elevation_m,
                  double bandlower = 0.0, double bandupper = 0.0);

// Multiscale plan S3.5: a land channel cell (not sea) with a non-channel, non-sea 8-neighbour
// whose elevation differs from `before` (a form changed it) takes min(its elevation, its cell mean
// on the changed field), the mean as regional_solve.cpp's: b + min(h, cap_m), b its bed
// (channels.bed_m), h the smaller of the mean height above b of those neighbours and of the changed
// ones alone (at least 0), so the axis never stands above its carved flanks' mean (a form narrower
// than the diagonal leaves the diagonal neighbours high). So a form lowers a channel cell only
// towards its cell mean, never back to its bed under a sloping side, and never raises it. cap_m <= 0
// (the single pass, refinement.channel_cell_slope 0): nothing.
void channelcellmeans(const RegionalWindow& window, const FineChannels& channels, const std::vector<uint8_t>& sea,
                      double cap_m, const std::vector<float>& before, std::vector<float>& elevation_m);

// Distances on the ground for the forms in the loop (S3.5): the east-west over the north-south
// fine cell size of domain row y (metric.rowdx_m / dy_m, the cosine of the row's latitude: a pure
// function of the absolute row; 1 outside the metric), by which a form scales its column offsets,
// and the column extent of a form of `radius` fine cells (north-south) at that ratio, at most k.
double eastwestratio(const PatchMetric& metric, int y);
int columnextent(double radius, double ratio, int k);

struct LevelContext;  // regional_pass.hpp (the hook interface, multiscale plan S3.0)
struct LevelState;

// The band of a level's forms in its cells (S2.3, S3.5): [detailcells / 2, detailcells x factor
// / 2) in the loop, the upper end infinite at level 1 (the forms above its band, each clamped to
// one global cell); 0, 0 in the single pass (factor < 2).
void formband(const LevelContext& context, double& bandlower, double& bandupper);

// The cap of a channel cell's mean, refinement.channel_cell_slope x the level's north-south cell
// size (regional_solve.hpp); 0 in the single pass.
double channelcellcap_m(const LevelContext& context);

// LevelStep::incision: with context.incision active, state.incisednetwork (the inputs' network
// with the incised filled elevations) and state.bedsurface_m (state.surface_m lowered by the
// incision upsampled, at least 0) for channelslevel; nothing otherwise.
void incisionlevel(const LevelContext& context, LevelState& state);

// LevelStep::valleyfloors: with details.valleys, carvevalleys on state.elevation_m with the
// level's band (formband), then channelcellmeans with channelcellcap_m. The loop's diagnostics
// count the kept cells it changed (LevelDiagnostics::valleycells).
void valleyfloorslevel(const LevelContext& context, LevelState& state);
}
