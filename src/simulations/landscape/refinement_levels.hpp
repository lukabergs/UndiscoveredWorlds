#pragma once

// The loop over the levels of the refinement hierarchy (docs/architecture/multiscale-plan.md,
// sections 1 and 3, "Kernels"; contracts written in step S0.4, the loop implemented in
// S2.3, which makes regional_pass.cpp its caller). Each level l = 1 .. L of a halo chain
// (refinement_grid.hpp) is the prolongation of level l - 1 (upsampling.hpp's Catmull-Rom
// generalised to one factor step) followed by the steps of LevelStep in their fixed order,
// then the level's diagnostics. S2.3 calls the existing units directly in this order (the
// ones regional_pass.cpp calls today, given a level window); the hook interface through
// which parallel detail steps (S3.x, S5.2) add their units is fixed at S3.0, not here.
// Contracts of the order:
// - Nothing after `projection` writes the elevation: `water` reads it (and the lapse-rate
//   climate after the loop). Today regionalwater (regional_water.hpp, R1 and R4) does
//   write it: it lowers lake beds below the surface and raises dry shores by at most
//   2 shores.terrain_m. S2.3 must move those writes before the projection or bound the
//   parent-mean error they leave (an open item of the plan, not settled here).
// - The lapse-rate climate (temperature, precipitation, climate class) is computed after
//   the loop, at the final level only, from its projected elevation.
// - LevelDiagnostics::blockmeanerror_m is measured after `presmoothing`.
// Each detail runs at the level whose cell size matches its characteristic size on the
// ground (detaillevel). Planet-free: the parameters are filled from
// tuning::landscape::refinement by the planet adapter
// (src/simulations/regional/regional_inputs.hpp, refinementparameters).

#include <vector>

#include "simulations/landscape/refinement_grid.hpp"
#include "simulations/landscape/refinement_projection.hpp"
#include "simulations/landscape/regional_pass.hpp"
#include "simulations/landscape/regional_types.hpp"

namespace landscape
{
// landscape.refinement (configs/parameters/landscape.json).
struct RefinementParameters
{
    int factor = 4;                         // level cells per parent cell edge (2 or 4; S2.6 decides)
    double projectionweightchannel = 0.0;   // the projection's weight on channel-bed cells (0: beds stay)
    int storedscalemax = 256;               // the finest k = factor^L written as a bundle; finer levels are
                                            // render tiles computed on demand
    double detailcells = 2.0;               // a detail is resolved by at least this many cells of its level
    int tilecells = 1024;                   // a render tile's edge in level cells
};

// The steps of a level in their fixed order (section 3, "Kernels"), after the prolongation.
enum class LevelStep : int
{
    noiseoctave = 0,  // the band-limited noise octave of the level (regional_surface.hpp)
    presmoothing,     // the start surface's block-mean correction (upsampling.hpp, correctblockmeans,
                      // regional.mean_passes): keeps the per-parent multipliers slowly varying
    coast,            // the sea by local rules (regional_surface.hpp's regionalsea: regional_coast.hpp,
                      // regional_edges.hpp), before the channels, which end at the fine coast
    channels,         // the pinned channels (channel_network.hpp) and their tributaries (S3.1)
    incision,         // the channels' incision (regional_valleys.hpp, globalincision)
    solve,            // the analytical relief from the pinned channels (analytical_relief.hpp)
    valleyfloors,     // valley floors and sides (regional_valleys.hpp, carvevalleys)
    glacialforms,     // glacial troughs and cirques (regional_troughs.hpp, regional_cirques.hpp)
    foldfabric,       // fold trains (regional_folds.hpp, S5.2)
    projection,       // the exact projection onto the parent level (refinement_projection.hpp)
    water,            // lakes and ice (regional_water.hpp, regionalwater; regional_ice.hpp, iceabovetheela)
};
constexpr int levelstepcount = 11;

inline const char* levelstepname(LevelStep step)
{
    switch (step)
    {
    case LevelStep::noiseoctave: return "noise_octave";
    case LevelStep::presmoothing: return "presmoothing";
    case LevelStep::coast: return "coast";
    case LevelStep::channels: return "channels";
    case LevelStep::incision: return "incision";
    case LevelStep::solve: return "solve";
    case LevelStep::valleyfloors: return "valley_floors";
    case LevelStep::glacialforms: return "glacial_forms";
    case LevelStep::foldfabric: return "fold_fabric";
    case LevelStep::projection: return "projection";
    case LevelStep::water: return "water";
    }
    return "unknown";
}

// A level's cell size on the ground in km: the north-south spacing pi R / (H factor^level)
// with H = worldwidth / 2, which is also the east-west spacing at the equator (PatchMetric's
// dy_m / 1000 at k = factor^level). Throws std::invalid_argument when radius_m is not
// positive and finite, worldwidth < 4 or odd, factor < 2 or level < 0.
double levelcellsize_km(double radius_m, int worldwidth, int factor, int level);

// The level that runs a detail of characteristic size size_km on the ground: the coarsest
// level l in 1 .. levels whose cell size levelcellsize_km(l) is at most size_km /
// detailcells (refinement.detail_cells, 2: the detail is resolved by at least two cells of
// its level, Nyquist); levels when no level is that fine (the detail is finer than the
// finest level run); and 1 when even the global grid's cells are that small (level 0 is
// given, never run). Throws std::invalid_argument when size_km or detailcells is not
// positive and finite, levels < 1, or as levelcellsize_km.
int detaillevel(double size_km, double detailcells, double radius_m, int worldwidth, int factor, int levels);

// The target of a refinement: the global window's kept box and halo in global cells
// (window.scale is not read) and the finest level computed, at least 1; level L's halo is
// levelhalo(window.halo, factor, L) parent cells (refinement_grid.hpp: window.halo at
// level 1, plus the channels' reach of 2 global cells at deeper levels).
struct LevelTarget
{
    RegionalWindow window;
    int level = 1;
};

// What one level of the loop did (window-dependent, not fields).
struct LevelDiagnostics
{
    int level = 0;
    double cellsize_km = 0.0;
    LevelWindow window;                  // the level's window in the halo chain
    double blockmeanerror_m = 0.0;       // measured after LevelStep::presmoothing (kept box's parents)
    ProjectionReport projection;         // the exact projection's report
    bool reliefconverged = false;        // the level's solve (AnalyticalReliefResult)
    int changedreceivers = 0;
    std::vector<LevelStep> steps;        // the steps this level ran, in order
};

struct RefinementResult
{
    std::vector<LevelDiagnostics> levels;  // levels 1 .. target.level, index = level - 1
    // The final level's kept fields as runregionalpass returns them, the lapse-rate climate
    // computed at this level only. fields.window is the target's window with scale
    // factor^level and halo target.window.halo global cells (level 1's halo; the deeper
    // levels' halos are in `levels`), so regional_bundle.hpp writes it unchanged, with
    // Metadata::level and factor set, and the window block's halo is target.window.halo.
    RegionalFields fields;
};

// The loop. Throws std::invalid_argument when the target is not a window on the inputs'
// grid (as runregionalpass), target.level < 1, refinement.factor < 2, a global field has
// the wrong size, the radius is not positive or the sea level not finite. Declared in
// S0.4, implemented in S2.3 (it throws std::logic_error until then).
RefinementResult runrefinement(const RegionalGlobalInputs& inputs, const LevelTarget& target,
                               const RegionalParameters& parameters, const RefinementParameters& refinement);
}
