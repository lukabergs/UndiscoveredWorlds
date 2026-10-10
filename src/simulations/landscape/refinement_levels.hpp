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
// The hook interface through which the detail steps (S3.x, S5.2) add their units (S3.0):
// LevelContext, LevelState, ParentFields and regionallevelsteps in regional_pass.hpp, one
// function per LevelStep in its unit's header; the loop builds the context per level, puts the
// prolonged parent into state.start at level >= 2, calls regionallevelsteps and keeps the
// level's fields as the next level's ParentFields.
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
// - Between the lakes and the projection, with refinement.smooth_passes > 0, the per-parent
//   correction is smoothed (LevelStep::smoothing, smoothparentmeans with the projection's
//   bounds built on the level's field; S2.7), and the projection's bounds are rebuilt on the
//   smoothed field. Each pass reaches 2 parent cells: the target's halo must be at least
//   2 x meanpasses + 2 + 2 x smooth_passes global cells.
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
// - Masks inherited and refined within one parent cell (S2.8; §3 "S2.7: smooth before exact",
//   item 3, and "S2.8: masks without walls"): at level >= 2 a cell may be sea only where its
//   parent or a parent 8-neighbour is sea (and the global gate allows it); every other cell is
//   land at its own height, below the sea level allowed and counted (LevelDiagnostics::
//   landbelowsea, maxlandbelowsea_m; regionalsea's parentgate); level 1 keeps the single pass's
//   sea. At every level the lakes flood every cell below a lake's surface 4-connected to
//   its core inside the lake cell's 3 x 3 global cells (regionalwater's LakeInheritance), the
//   core at level >= 2 the cells below the surface under the parent's water, the flood within
//   one parent cell of it. The parent's masks are read from its kept box (sea and lake surface
//   after the projection's drying). Reach: the flood is a component inside one lake cell's
//   3 x 3 zone, read on the solved elevation: a kept cell's water depends on the solved field
//   up to 2 global cells beyond it (level 1; 1 parent cell + 3 level cells at level >= 2, also
//   within 2 global cells), on top of the smoothing's reach, so the loop's halo adds
//   lakefloodreach, 2 parent cells at every level (2 global cells at level 1).
// - The halo budget per level (S3.8), in parent cells at level L (s = factor^(L - 1) parent cells
//   per global cell), each reach on top of the one before: the presmoothing 2 x meanpasses; the
//   channels and the solve channelreachglobalcells s plus the solve's margin
//   solvemarginglobalcells s (refinement_grid.hpp; found by trial: the solve's reach is not
//   bounded, LevelDiagnostics::ringbasincells), the forms beside the solve (a form's extent is
//   at most one parent cell: k level cells at level 1, under factor level cells deeper) counted
//   on top as formreach; the channel-cell mean's one level cell (channelmeanreach, one parent
//   cell); the lake flood lakefloodreach; the smoothing 2 x smooth_passes; the tributaries
//   tributaryreach. levelhalo gives s (4 + ceil((halo - 4) / s)): the global-cell reaches kept,
//   the parent-cell ones rounded up to whole global cells. At the defaults (halo 30, factor 4)
//   levels 1, 2, 3 take 30, 44, 96 parent cells (30, 11, 6 global cells; 28, 36, 64 parent
//   cells at halo 28 before S3.8).
// Planet-free: the parameters are filled from tuning::landscape::refinement by the planet
// adapter (src/simulations/regional/regional_inputs.hpp, refinementparameters).

#include <cstdint>
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
    // The smoothing before the projection (S2.7, refinement_projection.hpp, smoothparentmeans):
    // its passes (0: none, the S2.3 loop bit for bit; the parameter file sets 6) and relaxation.
    int smoothpasses = 0;
    double smoothrelaxation = 1.4;
    // Phase S3's keys (contracts of S3.0; the struct keeps the behaviour-keeping values, the
    // parameter file sets the defaults):
    int prolongationlimit = 0;        // S3.1: 1 clamps each prolonged child to the range of its 2 x 2 bracketing
                                      // parents (prolong); 0: the plain Catmull-Rom (the parameter file's
                                      // value until S3.2's measurement)
    double channelcellslope = 0.0;    // S3.2: theta, a channel cell's elevation b + min(h, theta c_l); 0: the bed
    int smoothbasis = 0;              // S3.2's fallback: smoothparentmeans' basis, 0 Catmull-Rom (C1), 1 the cubic
                                      // B-spline (C2)
    int parentlift = 0;               // S3.2: the solve against its parent at level >= 2 (solvelevel): 0 none; 1 the
                                      // parent's lift_m prolonged onto the unpinned cells; 2 the parent's final
                                      // elevation less its solved elevation (ParentFields::solved_m) prolonged
};

// The prolongation's reach in parent cells beyond a child's parent cell (Catmull-Rom:
// the halo chain's `reach`, refinement_grid.hpp).
constexpr int prolongationreach = 2;

// The lake flood's reach in parent cells (S2.8: a component inside one lake cell's 3 x 3 zone: 2
// global cells at level 1, 1 parent cell + 3 level cells at level >= 2).
constexpr int lakefloodreach = 2;

// The forms' reach in parent cells (S3.8): a valley floor, trough or cirque reaches at most k
// level cells (one global cell) at level 1 and under factor level cells (the band) deeper.
constexpr int formreach = 1;

// The channel-cell mean's reach in parent cells (S3.8): one level cell (channelcellmeans and
// solvelevel's cell mean, refinement.channel_cell_slope), counted as one parent cell.
constexpr int channelmeanreach = 1;

// The halo in global cells the loop requires: 2 x meanpasses + channelreachglobalcells (the single
// pass's rule), formreach and channelmeanreach (S3.8), lakefloodreach (S2.8), 2 x smoothpasses
// (S2.7's smoothing); with the tributaries on (S3.3), tributaryreach more (regional_tributaries.hpp:
// parent cells). levelhalo keeps globalreachcells in global cells at every level and turns the rest
// into at least as many parent cells. 30 at the defaults (meanpasses 6, smoothpasses 6).
constexpr int refinementhalo(int meanpasses, int smoothpasses, bool tributaries = false)
{
    return 2 * meanpasses + channelreachglobalcells + formreach + channelmeanreach + lakefloodreach + 2 * smoothpasses +
           (tributaries ? tributaryreach : 0);
}

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
    channels,         // the pinned channels (channel_network.hpp) and their tributaries (S3.3,
                      // regional_tributaries.hpp)
    solve,            // the analytical relief from the pinned channels (regional_solve.hpp, analytical_relief.hpp)
    valleyfloors,     // valley floors and sides (regional_valleys.hpp, carvevalleys)
    glacialforms,     // glacial troughs and cirques (regional_troughs.hpp, regional_cirques.hpp)
    foldfabric,       // fold trains (regional_folds.hpp, S5.2)
    lakes,            // lakes (regional_water.hpp, regionalwater): writes the elevation, gives the lake beds
    smoothing,        // the per-parent correction smoothed with the projection's bounds (refinement_projection.hpp,
                      // smoothparentmeans; S2.7), when refinement.smooth_passes > 0
    projection,       // the exact projection onto the parent level (refinement_projection.hpp)
    ice,              // ice above the ELA (regional_ice.hpp, iceabovetheela): the final level, reads the elevation
};
constexpr int levelstepcount = 13;

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
    case LevelStep::smoothing: return "smoothing";
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
// parent values and the child's absolute position. With `limit` (refinement.prolongation_limit
// 1, step S3.1) each child is clamped to the range [min, max] of its 2 x 2 bracketing parents
// (the stencil's two middle rows and columns, the parents whose centres enclose the child's
// centre), in double before the rounding to float, so the Catmull-Rom undershoot digs no pit
// below a trench's floor and its overshoot raises no child above a ridge; still a pure function
// of the parent values and the child's absolute position (window-independent). A linear field
// is not changed (its interpolant lies inside the bracket). Without it the plain Catmull-Rom of
// S2.3, bit for bit. Throws std::invalid_argument when the windows do not fit together,
// parentfield is not parent.domainsize() long, or a stencil cell lies outside the parent's kept
// box.
struct Prolongation
{
    std::vector<float> values;  // level.domainsize(), row-major
    // The clamp report over the level's kept box: the cells whose float value the limit changed
    // and the largest |change| (0 without the limit).
    int clamped = 0;
    double maxclamp_m = 0.0;
};
Prolongation prolong(const LevelWindow& parent, const std::vector<float>& parentfield, const LevelWindow& level, bool limit);

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
// a whole number of global cells reaching globalreachcells plus the parent-cell reaches' share
// at deeper levels).
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
    SmoothingReport smoothing;           // the smoothing before the projection (lambda before and after it)
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
    double maxlandbelowsea_m = 0.0;      // the largest depth of such a cell below the sea level (S2.8)
    std::int64_t maxlandbelowseacolumn = -1;  // its absolute column and row at the level (-1: none)
    std::int64_t maxlandbelowsearow = -1;
    int rimraised = 0;                   // kept cells the lake flood's rim raise lifted (FineWater::rimraise_m)
    double maxrimraise_m = 0.0;          // the largest such lift
    // S2.8's mask checks (kept box): sea cells whose parent (level 1: global cell) is land; of
    // them, those with no sea among the parent and its 8 neighbours (0 at level >= 2 by the
    // inherited gate); dry cells (not sea, not channel) below the surface of a wet 4-neighbour
    // after the projection.
    int seainlandparents = 0;
    int seabeyondparent = 0;
    int drybelowlake = 0;
    // Phase S3's diagnostics (written by S3.0's hooks; kept box, after the projection unless said):
    // S3.1, the limited prolongation (level >= 2; Prolongation's report): the kept cells it clamped
    // and its largest move (0 at level 1 and with refinement.prolongation_limit 0).
    int prolongationclamped = 0;
    double maxprolongationclamp_m = 0.0;
    // S3.2: the 99th percentile (nearest rank) of the kept cells' |move| by the smoothing, and the
    // largest elevation above its bed (LevelState::bed_m) of a dry land channel cell (not sea, not
    // wet; 0 at least): 0 at S3.0 unless a step after the channels moved a bed (the projection's
    // bedsmoved and forced fallbacks).
    double smoothmove99_m = 0.0;
    double maxchannelrise_m = 0.0;
    // S3.3: the land cells (not sea), the land channel cells and the drainage density, channel
    // cells / (land cells x cellsize_km) in 1/km (each channel cell one cell size of channel, each
    // cell the square of the north-south cell size); the tributaries and the domain cells they pinned.
    int landcells = 0;
    int channelcells = 0;
    double drainagedensity_per_km = 0.0;
    int tributaries = 0;
    int tributarycells = 0;
    int tributarybedsabove = 0;          // domain tributary cells whose bed stands above the start surface (0 by rule)
    // S3.4: the lake edges (pairs of east-west or north-south neighbours in the kept box, one wet,
    // lake surface > 0, the other not) and the share of them on a global cell's boundary
    // (lakeedgealigned; about 1 / k_l for an edge that ignores the global grid, 1 for square lakes);
    // the land roughness, the RMS of z_E + z_W + z_N + z_S - 4 z over the kept dry land cells (not
    // sea, not wet) whose four neighbours are dry land, in m at the level's cell size.
    int lakeedges = 0;
    double lakeedgealigned = 0.0;
    double landroughness_m = 0.0;
    // S3.5: the kept cells the valley floors, the troughs and the cirques changed at this level.
    int valleycells = 0;
    int troughcells = 0;
    int cirquecells = 0;
    // S3.8: the kept cells whose solve sub-basin (the unpinned cells connected to it by orthogonal
    // steps and diagonal steps not between two pinned corners, analytical_relief.hpp) touches the
    // domain's outer ring on a window edge (not at a pole, where every window clips alike). The
    // solve's result there may depend on the window: bit identity with another window is not
    // guaranteed there, only found by the overlap tests (refinement_grid.hpp, solvemarginglobalcells).
    int ringbasincells = 0;
    bool reliefconverged = false;        // the level's solve (AnalyticalReliefResult)
    int changedreceivers = 0;
    double seconds = 0.0;                // the level's wall-clock time
    std::vector<LevelStep> steps;        // the steps this level ran, in order
};

// The hook interface's bookkeeping (regional_pass.hpp): records `step` in
// state.diagnostics->steps (nothing in the single pass, whose state has no diagnostics).
void recordlevelstep(LevelState& state, LevelStep step);

// LevelStep::presmoothing. At level >= 2 (context.parent set): correctparentmeans of
// state.start.start_m against the parent's elevation (regional.mean_passes), its
// levelparentmeanerror as state.start.blockmeanerror_m, and the lapse-rate base (the global
// elevation upsampled at the level); at level 1 and in the single pass the correction is
// regionalstartsurface's (noiselevel). Then the start surface is handed on: the level's
// blockmeanerror_m (also LevelDiagnostics::blockmeanerror_m), sharepasses, sharemisses and
// lapsebase_m, and state.surface_m. Throws std::invalid_argument when the start surface does
// not cover the window's domain.
void presmoothinglevel(const LevelContext& context, LevelState& state);

// The kept cells of `window` (its finewidth() x fineheight() box in the domain) whose value
// differs between two domain fields (the forms' cell counts of LevelDiagnostics).
int keptcellschanged(const RegionalWindow& window, const std::vector<float>& before, const std::vector<float>& after);

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
    // Their sea and lake surface after the projection on the same kept boxes (S2.8).
    std::vector<std::vector<uint8_t>> coarsersea;
    std::vector<std::vector<float>> coarserlakesurface_m;
    // The halo chain of levels 0 .. target.level (refinement_grid.hpp, halochain), index =
    // level; the resample of resampleregionalfields reads it (step S2.4).
    std::vector<LevelWindow> chain;
};

// The loop. Throws std::invalid_argument when the target is not a window on the inputs'
// grid (as runregionalpass), its halo is below refinementhalo(meanpasses, smoothpasses),
// target.level < 1, refinement.factor < 2, refinement.smoothpasses < 0 or
// refinement.smoothrelaxation not positive and finite, factor^level does not fit an int, a global field has the wrong
// size, the radius is not positive or the sea level not finite, and as halochain (a level's
// domain wider than its grid).
RefinementResult runrefinement(const RegionalGlobalInputs& inputs, const LevelTarget& target,
                               const RegionalParameters& parameters, const RefinementParameters& refinement);

// --- Output at a requested resolution (section 1, "Output at any resolution"; step S2.4). ---

// The source level of a request for metresperpixel, as resampleto chooses it: the coarsest
// level l in 0 .. maxlevel whose cell size pi R / (worldheight factor^l) is at most
// metresperpixel (1 + 1e-9); -1 when level maxlevel is still coarser. Throws
// std::invalid_argument when metresperpixel or radius_m is not positive and finite,
// worldwidth < 4 or odd, factor < 2 or maxlevel < 0.
int resamplesourcelevel(double metresperpixel, double radius_m, int worldwidth, int factor, int maxlevel);

// A refinement's final-level fields resampled to a requested metres per pixel.
struct ResampledRegionalFields
{
    // width x height pixels each, row-major, north row first; window and diagnostics are the
    // refinement's (window.scale the source level's factor^level), so window.finewidth() is
    // not the pixel grid.
    RegionalFields fields;
    int width = 0;
    int height = 0;
    ResampleReport report;
};

// result.fields resampled to metresperpixel on resampleto's pixel grid (its source level
// must be result's final level: compute the refinement to resamplesourcelevel's level). The
// pixel's source cells are those resampleto averages (aligned: the round(n) x round(n) block;
// otherwise every cell its box overlaps, weighted by the overlap area). The pixels are
// anchored at the kept box's north-west corner: aligned pixels tile whole global cells, so
// two windows give the same pixels where they overlap; unaligned pixels do not, so two
// windows give different pixel grids (no overlap identity for a resampled output). Per field:
// - elevation, temperatures, precipitations, ice thickness: resampleto's box average;
// - sea: 1 when the sea cells hold at least half the pixel's weight (ties: sea);
// - lake surface: when not sea and the wet cells (lake surface > 0) hold at least half the
//   weight, their weighted mean surface; 0 otherwise;
// - discharge, January and July flow: the largest value of the pixel's cells;
// - river direction: the direction of the cell of the largest discharge (the first in
//   row-major order on a tie);
// - climate class: the class of the cell holding the pixel's centre.
// Throws std::invalid_argument when the chain is empty, a field is not the final level's
// kept size, the source level is not the final level, and as resampleto.
ResampledRegionalFields resampleregionalfields(const RefinementResult& result, double radius_m, double metresperpixel);
}
