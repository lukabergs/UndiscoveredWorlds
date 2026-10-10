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
#include "simulations/landscape/refinement_grid.hpp"
#include "simulations/landscape/refinement_projection.hpp"
#include "simulations/landscape/regional_glacial.hpp"
#include "simulations/landscape/regional_grid.hpp"
#include "simulations/landscape/regional_surface.hpp"
#include "simulations/landscape/regional_tributaries.hpp"
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

// The single pass's smallest halo in global cells: 2 x meanpasses + 2 (the block-mean
// correction reaches 2 global cells per pass; two more keep the channels and the solve's edge
// out of the kept window), plus tributaryreach (regional_tributaries.hpp; the single pass's
// parent cell is the global cell) when the tributaries are on (multiscale plan S3.3).
constexpr int regionalhalo(int meanpasses, bool tributaries)
{
    return 2 * meanpasses + 2 + (tributaries ? tributaryreach : 0);
}

// Throws std::invalid_argument when the window is not on the inputs' grid, a global
// field has the wrong size, the radius is not positive, the sea level not finite, or
// window.halo < regionalhalo(meanpasses, tributaries.enabled != 0);
// rasterisechannels and solveanalyticalrelief throw on their own bad inputs.
// runregionalpass is checkregionalinputs, regionalwarp, globalincision, regionallevelsteps
// (whose noise step is regionalstartsurface) and finishregionallevel in that order.
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
    // S2.8 (masks inherited and refined within one parent cell): at level >= 2 the parent's
    // sea gate (regionalsea's parentgate; null at level 1 and in the single pass) and the lake
    // flood with the parent's lake masks (regionalwater's LakeInheritance; flood off in the
    // single pass). The pointed-to vectors hold one value per domain cell and outlive the call.
    const std::vector<uint8_t>* parentseagate = nullptr;
    LakeInheritance lakes;
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

// --- The hook interface of a level (multiscale plan S3.0, §3 "Phase S3 re-planned"). ---
// Each LevelStep (refinement_levels.hpp) is one function `void <step>level(const LevelContext&,
// LevelState&)` declared in its unit's header: noiselevel and coastlevel (regional_surface.hpp),
// presmoothinglevel (refinement_levels.hpp), incisionlevel and valleyfloorslevel
// (regional_valleys.hpp), channelslevel (channel_network.hpp; it calls tributarieslevel,
// regional_tributaries.hpp), solvelevel (regional_solve.hpp), glacialformslevel
// (regional_glacial.hpp), lakeslevel (regional_water.hpp), smoothinglevel and projectionlevel
// (refinement_projection.hpp). regionallevelsteps calls them in LevelStep order; the ice is
// finishregionallevel's. A step reads the context and the state's fields the earlier steps
// wrote, writes its own, and records itself in the level's diagnostics when it acts
// (recordlevelstep, refinement_levels.hpp). Every rule stays a pure function of absolute position, the seed and the
// parent window within the level's halo.

struct RefinementParameters;  // refinement_levels.hpp
struct LevelDiagnostics;      // refinement_levels.hpp

// The parent level's fields after its projection, on its domain (window.domainsize() values
// each, row-major; a level reads them inside the parent's kept box only, as prolong does).
// Null in the context at level 1 (its parent is the global grid) and in the single pass.
struct ParentFields
{
    LevelWindow window;                  // the parent level's window in the halo chain
    std::vector<float> elevation_m;      // the projected elevation (the next level's prolongation source)
    std::vector<float> lift_m;           // its LevelState::lift_m
    std::vector<float> solved_m;         // its elevation after the solve (LevelState::solved_m)
    std::vector<uint8_t> sea;
    std::vector<float> lakesurface_m;    // after the projection's drying (0: dry)
    FineChannels channels;               // channel, direction, bed, discharge, flows and global cell
    std::vector<float> bed_m;            // its LevelState::bed_m
};

// What a level step reads: the inputs and the level's place in the hierarchy.
struct LevelContext
{
    const RegionalGlobalInputs& inputs;
    const RegionalWindow& window;        // the level's window as the units take it (scale k_l; the single pass's)
    const RegionalParameters& parameters;
    const FineWarp& warp;                // regionalwarp's for `window`
    const GlobalIncision& incision;      // globalincision's, computed once per pass
    const RegionalDetails& details;
    int level = 0;                       // 1 .. L in the refinement loop; 0 in the single pass
    int factor = 0;                      // refinement.factor in the loop; 0 in the single pass
    const RefinementParameters* refinement = nullptr;     // null in the single pass (no smoothing, no projection)
    const LevelWindow* levelwindow = nullptr;            // the level's window in the halo chain (null: single pass)
    const LevelWindow* projectionparent = nullptr;       // the window projected onto (level 1: the global level 0)
    const std::vector<float>* projectionfield = nullptr;  // its elevation on its domain
    const ParentFields* parent = nullptr;                // null at level 1 and in the single pass
};

// What the steps of a level write: RegionalLevel's fields (the domain fields
// finishregionallevel reads) and the steps' intermediate fields, each on the window's domain.
struct LevelState : RegionalLevel
{
    StartSurface start;                  // noiselevel's (level 1, single pass: regionalstartsurface; deeper: the
                                         // prolongation with the octave band added), handed on by presmoothinglevel
    std::vector<float> noiseband_m;      // noiselevel's: the level's share of the noise amplitude per cell (S3.4:
                                         // the band share times noise_amplitude_m + noise_relief_fraction x the
                                         // relief, before the taper); lakeslevel caps its rim raise with it
    std::vector<float> surface_m;        // the working surface: the start surface, its sea centres set by the coast
                                         // rule, then the land channel beds pinned (the lakes' start surface)
    GlobalNetwork incisednetwork;        // incisionlevel's (incision active only): the network's incised filled elevations
    std::vector<float> bedsurface_m;     // incisionlevel's: surface_m lowered by the upsampled incision (empty: inactive)
    std::vector<float> bed_m;            // the channel beds (FineChannels::bed_m, meaningful on channel cells; S3.2
                                         // keeps them here when a channel cell's elevation becomes its cell mean)
    std::vector<uint8_t> pinned;         // the solve's pinned cells (channels, sea, the domain's outer ring)
    std::vector<float> solved_m;         // the elevation after the solve, before the forms and the lakes (S3.2: the
                                         // channel cells' means kept, the parent's lift not added; solvelevel)
    int tributaries = 0;                 // tributarieslevel's count (0: none)
    int tributarycells = 0;              // the domain cells its tributaries pinned
    int tributarybedsabove = 0;          // of them, those whose bed stands above state.surface_m (S3.3: 0 by rule)
    // Refinement loop only (empty in the single pass):
    ProjectionBounds bounds;             // smoothinglevel's, rebuilt after the smoothing; projectionlevel's
    std::vector<float> ownlake_m;        // the own global lake cell's surface per cell (projectionbounds' ownlake_m)
    std::vector<float> lift_m;           // final elevation - the elevation entering the smoothing
    LevelDiagnostics* diagnostics = nullptr;  // the loop's record of the level (null: the single pass)
};

// The steps of a level from the noise to the projection, in LevelStep order (the single pass:
// the noise is regionalstartsurface's, and the smoothing and the projection do nothing).
void regionallevelsteps(const LevelContext& context, LevelState& state);

// Steps 5b to 8: the ice above the ELA (iceabovetheela, adding to level.water), the rivers
// and the lapse-rate climate from level.elevation_m, and the kept window's fields with the
// level's diagnostics; fields.window is `window`.
RegionalFields finishregionallevel(const RegionalGlobalInputs& inputs, const RegionalWindow& window,
                                   const RegionalParameters& parameters, RegionalLevel& level);
}
