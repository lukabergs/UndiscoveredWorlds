#pragma once

// The loop over the levels of the refinement hierarchy (docs/architecture/multiscale-plan.md,
// sections 1 and 3, "Kernels" and "S2.3's level placement"; contracts written in step S0.4,
// the loop implemented in S2.3, which makes regional_pass.cpp's parts its steps). Level 1 is
// the single-scale pass at k = factor (regionalstartsurface and regionallevelsteps,
// regional_pass.hpp) on regionalwindowof(chain[1]); each level l = 2 .. L of the halo chain
// (refinement_grid.hpp) is the prolongation of level l - 1 (prolong) plus the level's
// octave band (regional_surface.hpp, addlevelnoise) with its parent means corrected
// (correctparentmeans), then the same units at k_l = factor^l. Every level ends with the
// exact projection onto its parent (level 1 onto the global elevation) and its diagnostics.
// The hook interface through which parallel detail steps (S3.x, S5.2) add their units is
// fixed at S3.0, not here.
// Contracts of the order:
// - Prolonged from the parent: the elevation only. Recomputed at every level from the
//   global fields and the level's elevation: the sea, the channels and their discharge,
//   the solve, the lakes, the glacial fields.
// - The lakes (regionalwater, which writes the elevation: lake beds to the global lake depth,
//   uncapped, and shores) run before the projection. The projection works by bounds
//   (refinement_projection.hpp, projectionbounds: sea and wet cells may sink but not rise
//   past their water, dry land may rise but not sink to the sea level or the water beside
//   it, land channel beds take the channel rate). After it a wet non-channel cell at or above
//   its surface is dried (its lake surface set to 0; LevelDiagnostics::driedcells), and the
//   dry cells at or below a wet 8-neighbour's surface (shoresbelowlake) and at or below the
//   sea level (landbelowsea) are counted. Nothing writes the elevation after the projection.
// - The ice above the ELA (iceabovetheela) and the lapse-rate climate (temperature,
//   precipitation, climate class) run after the loop, at the final level only, from its
//   projected elevation.
// - Valley floors, glacial troughs and cirques run at every level (each whenever its size key
//   is above 0); a form is carved at the level where its half-width or radius in level cells
//   lies in [detail_cells / 2, detail_cells x factor / 2) (RegionalDetails), and finer levels
//   inherit it by the prolongation, so a level's work depends on its own cell size only and
//   level 1 of any L is the one-level result. The noise band, the presmoothing, the coast
//   (deeper levels: regionalsea on the level's surface, the warp only as the position warp of
//   the sea gate and the lakes), the channels, the solve and the lakes run at every level; the
//   incision is globalincision, computed once and upsampled per level.
// - LevelDiagnostics::blockmeanerror_m is measured after `presmoothing`.
// Planet-free: the parameters are filled from tuning::landscape::refinement by the planet
// adapter (src/simulations/regional/regional_inputs.hpp, refinementparameters).

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

// The prolongation's reach in parent cells beyond a child's parent cell (Catmull-Rom:
// the halo chain's `reach`, refinement_grid.hpp).
constexpr int prolongationreach = 2;

// The steps of a level in their fixed order (section 3, "S2.3's level placement", item 5),
// after the prolongation.
enum class LevelStep : int
{
    noiseoctave = 0,  // the band-limited noise octave of the level (regional_surface.hpp)
    presmoothing,     // the parent-mean correction (level 1: correctblockmeans in the start surface,
                      // regional.mean_passes; deeper: correctparentmeans): keeps the per-parent
                      // multipliers slowly varying
    coast,            // the sea by local rules (regional_surface.hpp's regionalsea: regional_coast.hpp,
                      // regional_edges.hpp), before the channels, which end at the fine coast
    incision,         // the channels' incision (regional_valleys.hpp, globalincision upsampled)
    channels,         // the pinned channels (channel_network.hpp) and their tributaries (S3.1)
    solve,            // the analytical relief from the pinned channels (analytical_relief.hpp)
    valleyfloors,     // valley floors and sides (regional_valleys.hpp, carvevalleys)
    glacialforms,     // glacial troughs and cirques (regional_troughs.hpp, regional_cirques.hpp)
    foldfabric,       // fold trains (regional_folds.hpp, S5.2)
    lakes,            // lakes (regional_water.hpp, regionalwater): writes the elevation, gives the lake beds
    projection,       // the exact projection onto the parent level (refinement_projection.hpp)
    ice,              // ice above the ELA (regional_ice.hpp, iceabovetheela): the final level, reads the elevation
};
constexpr int levelstepcount = 12;

inline const char* levelstepname(LevelStep step)
{
    switch (step)
    {
    case LevelStep::noiseoctave: return "noise_octave";
    case LevelStep::presmoothing: return "presmoothing";
    case LevelStep::coast: return "coast";
    case LevelStep::incision: return "incision";
    case LevelStep::channels: return "channels";
    case LevelStep::solve: return "solve";
    case LevelStep::valleyfloors: return "valley_floors";
    case LevelStep::glacialforms: return "glacial_forms";
    case LevelStep::foldfabric: return "fold_fabric";
    case LevelStep::lakes: return "lakes";
    case LevelStep::projection: return "projection";
    case LevelStep::ice: return "ice";
    }
    return "unknown";
}

// A level's cell size on the ground in km: the north-south spacing pi R / (H factor^level)
// with H = worldwidth / 2, which is also the east-west spacing at the equator (PatchMetric's
// dy_m / 1000 at k = factor^level). Throws std::invalid_argument when radius_m is not
// positive and finite, worldwidth < 4 or odd, factor < 2 or level < 0.
double levelcellsize_km(double radius_m, int worldwidth, int factor, int level);

// (S0.4's level assignment; the loop no longer uses it since S2.3's rework, which places
// every form by its own size in level cells, RegionalDetails.)
// The level that runs a detail of characteristic size size_km on the ground: the coarsest
// level l in 1 .. levels whose cell size levelcellsize_km(l) is at most size_km /
// detailcells (refinement.detail_cells, 2: the detail is resolved by at least two cells of
// its level, Nyquist); levels when no level is that fine (the detail is finer than the
// finest level run); and 1 when even the global grid's cells are that small (level 0 is
// given, never run). Throws std::invalid_argument when size_km or detailcells is not
// positive and finite, levels < 1, or as levelcellsize_km.
int detaillevel(double size_km, double detailcells, double radius_m, int worldwidth, int factor, int levels);

// The prolongation of a parent level's field onto a level window's domain (level.level =
// parent.level + 1, the same world grid and factor): Catmull-Rom (upsampling.hpp's
// catmullrom, rows first, as bicubic) through the parent cell centres at every child centre,
// the child at offset o = column mod factor lying (o + 0.5) / factor - 0.5 parent cells from
// its parent's centre (likewise rows); parent columns wrap modulo the parent level's width,
// rows are clamped to the parent level's grid (as bicubic clamps at the poles). The stencil
// reads the parents of the level's domain and `prolongationreach` beyond, which must lie in
// the parent window's kept box (the halo chain guarantees it). A pure function of the
// parent values and the child's absolute position. Throws std::invalid_argument when the
// windows do not fit together, parentfield is not parent.domainsize() long, or a stencil
// cell lies outside the parent's kept box.
std::vector<float> prolong(const LevelWindow& parent, const std::vector<float>& parentfield, const LevelWindow& level);

// The parent-mean correction (upsampling.hpp's correctblockmeans against the parent level
// instead of the global grid), `passes` times on `fine` (level.domainsize() values): per
// pass the residual of every parent cell of the level's domain (the unweighted mean of its
// factor^2 children, summed in row-major order, minus the parent value), stored as float,
// the residuals of the two parent cells beyond the domain's edge taking the nearest domain
// parent's (a clamp; columns wrap when the domain spans the level's whole width), then
// interpolated by Catmull-Rom at the child centres as in prolong and subtracted. The reach
// is 2 x passes parent cells inside the domain's edge (the halo rule of levelhalo). Parent
// values are read from the parent window's kept box. Throws std::invalid_argument as
// prolong, when passes < 0, or when fine is not level.domainsize() long.
void correctparentmeans(const LevelWindow& parent, const std::vector<float>& parentfield, const LevelWindow& level,
                        int passes, std::vector<float>& fine);

// The target of a refinement: the global window's kept box and halo in global cells
// (window.scale is not read) and the finest level computed, at least 1; level l's halo is
// levelhalo(window.halo, factor, l) parent cells (refinement_grid.hpp: window.halo at level 1,
// a whole number of global cells reaching channelreachglobalcells plus the presmoothing's
// share at deeper levels).
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
    // The largest |mean over a global cell - the global elevation| of the projected level
    // field over the global cells of the kept box (the accumulated error of levels 1 .. l).
    double globalmeanerror_m = 0.0;
    // The bound of the mean errors (plan section 3, "S2.3's projection weights and detail
    // levels", item 4): half a float ulp of the largest |F| in the kept box, 2^(floor(log2
    // max |F|) - 24) m; the parent-mean error is at most it, the global-mean error at most
    // the sum of the levels' bounds.
    double meanbound_m = 0.0;
    int driedcells = 0;                  // kept wet non-channel cells at or above their surface after the projection
    int shoresbelowlake = 0;             // kept dry cells (not sea, not channel) at or below a wet 8-neighbour's surface
    int landbelowsea = 0;                // kept dry cells (not sea, not wet) at or below the sea level
    bool reliefconverged = false;        // the level's solve (AnalyticalReliefResult)
    int changedreceivers = 0;
    double seconds = 0.0;                // the level's wall-clock time
    std::vector<LevelStep> steps;        // the steps this level ran, in order
};

struct RefinementResult
{
    std::vector<LevelDiagnostics> levels;  // levels 1 .. target.level, index = level - 1
    // The final level's kept fields as runregionalpass returns them, the lapse-rate climate
    // computed at this level only. fields.window is the target's window with scale
    // factor^level and halo target.window.halo global cells (level 1's; the deeper
    // levels' halos are in `levels`), so regional_bundle.hpp writes it unchanged, with
    // Metadata::level and factor set, and the window block's halo is target.window.halo.
    // The diagnostics of fields are the final level's (blockmeanerror_m after its
    // presmoothing; sharepasses and sharemisses level 1's land-share rule).
    RegionalFields fields;
    // The projected elevation of each coarser level 1 .. target.level - 1 on its chain
    // window's kept box (levels[l - 1].window: columns x rows, row-major, north row first),
    // index = level - 1; the final level's is fields.elevation_m.
    std::vector<std::vector<float>> coarserelevation_m;
};

// The loop. Throws std::invalid_argument when the target is not a window on the inputs'
// grid (as runregionalpass), its halo is below 2 x meanpasses + 2, target.level < 1,
// refinement.factor < 2, factor^level does not fit an int, a global field has the wrong
// size, the radius is not positive or the sea level not finite, and as halochain (a level's
// domain wider than its grid).
RefinementResult runrefinement(const RegionalGlobalInputs& inputs, const LevelTarget& target,
                               const RegionalParameters& parameters, const RefinementParameters& refinement);
}
