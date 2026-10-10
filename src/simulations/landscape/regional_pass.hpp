#pragma once

// The regional pass (docs/architecture/erosion-hydrology-plan.md, Phase E, step E3,
// and section 6, "Decided for Phase E"): global fields and a window in, the window's
// fine fields out. Steps: the start surface (the global base elevation upsampled,
// band-limited noise added with an amplitude from the local relief and tapered near
// the sea level, block means corrected); the sea by local rules (Phase R, R1: before
// the channels, which end at the fine coast); the pinned channels; the analytical
// relief on the domain; lakes, rivers and climate on the fine grid by local rules;
// the halo discarded (Phase R: the start surface and the sea in regional_surface.hpp,
// lakes and ice in regional_water.hpp). Every step is a pure function of the global fields, the seed
// and absolute position within the halo's reach, so two windows give the same
// values on their overlap (test Landscape.RegionalOverlap). Planet-free: the planet
// adapter is src/simulations/regional/regional_inputs.hpp.

#include <cstdint>
#include <vector>

#include "simulations/landscape/channel_network.hpp"
#include "simulations/landscape/regional_glacial.hpp"
#include "simulations/landscape/regional_grid.hpp"
#include "simulations/landscape/regional_surface.hpp"
#include "simulations/landscape/regional_types.hpp"
#include "simulations/landscape/regional_valleys.hpp"
#include "simulations/landscape/regional_water.hpp"
#include "simulations/landscape/upsampling.hpp"

namespace landscape
{
// The fine fields of the kept window (finewidth() x fineheight(), row-major, north
// row first); the bundle writer stores them with the window.
struct RegionalFields
{
    RegionalWindow window;
    std::vector<float> elevation_m;
    std::vector<uint8_t> sea;
    std::vector<uint8_t> riverdirection;     // 1..8 on channel cells, 0 elsewhere
    std::vector<float> discharge_m3s;        // pinned channels: the global discharge; elsewhere the fine accumulation
    std::vector<int32_t> januaryflow, julyflow;  // legacy units, on pinned channels
    std::vector<float> lakesurface_m;        // 0 where no lake water
    std::vector<float> januarytemperature_c, julytemperature_c;
    std::vector<float> januaryprecipitation_mm, julyprecipitation_mm;
    std::vector<float> icethickness_m;
    std::vector<int16_t> climate;

    // Diagnostics (plan E3; window-dependent, not fields): the largest |block mean -
    // global value| over the kept window's global cells of the start surface after the
    // block-mean correction, before the channels are pinned and the solve; and the
    // solve's AnalyticalReliefResult::converged (the last two passes' receivers equal)
    // and changedreceivers (over the whole domain); the mean and largest lowering by the
    // solve (start surface minus result) over the kept window's unpinned land cells.
    double blockmeanerror_m = 0.0;
    bool reliefconverged = false;
    int changedreceivers = 0;
    double meanlowering_m = 0.0;
    double maxlowering_m = 0.0;
    // The land-share rule's passes and misses (StartSurface::sharepasses and sharemisses;
    // -1 and 0 without the rule).
    int sharepasses = -1;
    int sharemisses = 0;
};

// Throws std::invalid_argument when the window is not on the inputs' grid, a global
// field has the wrong size, the radius is not positive, the sea level not finite, or
// window.halo < 2 x meanpasses + 2 (the block-mean correction reaches 2 global cells
// per pass; two more keep the channels and the solve's edge out of the kept window);
// rasterisechannels and solveanalyticalrelief throw on their own bad inputs.
// runregionalpass is checkregionalinputs, regionalwarp, regionalstartsurface,
// globalincision, regionallevelsteps and finishregionallevel in that order.
RegionalFields runregionalpass(const RegionalGlobalInputs& inputs, const RegionalWindow& window,
                               const RegionalParameters& parameters);

// --- The pass in parts, for the refinement loop (docs/architecture/multiscale-plan.md,
// section 3, "S2.3's level placement"; refinement_levels.hpp). ---

// runregionalpass's checks of the window and the global fields (without the halo rule).
void checkregionalinputs(const RegionalGlobalInputs& inputs, const RegionalWindow& window);

// Which detail units a pass runs and how (the defaults are the single pass's).
struct RegionalDetails
{
    bool valleys = true;   // carvevalleys (regional_valleys.hpp)
    bool troughs = true;   // carvetroughs (regional_troughs.hpp)
    bool cirques = true;   // carvecirques (regional_cirques.hpp)
    // The refinement loop (multiscale plan S2.3, "S2.3's projection weights and detail
    // levels", item 3): with factor >= 2 a form is carved only where its half-width or radius
    // in the level's cells lies in [detailcells / 2, detailcells x factor / 2); 0 keeps the
    // single pass's clamp to k.
    int factor = 0;
    double detailcells = 0.0;
    // Item 2: the lakes' bed depth is the global lake depth, uncapped (regionalwater).
    bool globallakedepth = false;
};

// The domain fields of one pass after the lakes (domainfinewidth() x domainfineheight(),
// row-major), before the ice above the ELA and the climate.
struct RegionalLevel
{
    std::vector<float> lapsebase_m;          // StartSurface::lapsebase_m: the lapse-rate correction's reference
    std::vector<uint8_t> sea;
    FineChannels channels;
    std::vector<double> solveddischarge_m3_yr;  // the solve's accumulation (AnalyticalReliefResult)
    std::vector<float> elevation_m;          // after the solve, the details and the lakes
    FineWater water;
    FineGlacial glacial;
    // RegionalFields' diagnostics of this pass.
    double blockmeanerror_m = 0.0;
    int sharepasses = -1;
    int sharemisses = 0;
    bool reliefconverged = false;
    int changedreceivers = 0;
    double meanlowering_m = 0.0;
    double maxlowering_m = 0.0;
};

// Steps 2 to 5 of the pass on a start surface (StartSurface's start_m, lapsebase_m and
// diagnostics; start_m is the surface the sea, the channels and the solve start from):
// the sea (regionalsea, with `warp` for its gate), the pinned channels on the incision's
// bed surface, the solve, the details `details` selects (valleys, troughs, cirques), the
// lakes and the present-day ice of regionalwater. `incision` is globalincision's for these
// inputs and parameters (computed once by the caller).
RegionalLevel regionallevelsteps(const RegionalGlobalInputs& inputs, const RegionalWindow& window,
                                 const RegionalParameters& parameters, const FineWarp& warp, StartSurface start,
                                 const GlobalIncision& incision, const RegionalDetails& details);

// Steps 5b to 8: the ice above the ELA (iceabovetheela, adding to level.water), the rivers
// and the lapse-rate climate from level.elevation_m, and the kept window's fields with the
// level's diagnostics; fields.window is `window`.
RegionalFields finishregionallevel(const RegionalGlobalInputs& inputs, const RegionalWindow& window,
                                   const RegionalParameters& parameters, RegionalLevel& level);
}
