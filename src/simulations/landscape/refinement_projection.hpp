#pragma once

// The exact projection of the refinement hierarchy (docs/architecture/multiscale-plan.md,
// section 2, "The exact projection, worked out", and section 3, "Kernels"; contracts
// written in step S0.4, implemented in S2.2). A parent cell c at level L - 1 has the value
// v_c and factor^2 children C(c) at level L holding candidates F* (the prolongation plus
// the level's detail). The projection shifts the children by amounts that sum to
// factor^2 lambda_c, lambda_c = v_c - mean_C(c)(F*) (S2.2 by weights, F = F* + lambda_c w /
// mean_C(c)(w); S2.3's rework by bounds: a common shift capped at each child's headroom,
// projectlevel below), so mean_C(c)(F) = v_c to float rounding. Every mean here is the
// UNWEIGHTED mean over the factor^2 children (no cell areas), as upsampling.hpp's correctblockmeans takes the
// unweighted k x k block mean: unweighted means compose exactly across levels (the mean of
// factor^2 children's means of factor^2 grandchildren each is the mean of the factor^4
// grandchildren), so a level-L field averaged over a global cell's factor^(2L) cells equals
// the global value whatever the levels between. Sums are accumulated in double, the
// children of a parent visited in row-major order of their absolute position (north row
// first, west to east from the parent's west column), the order every window shares, so
// the result is bit-identical on overlapping windows; F is stored as
// float(double(F*) + shift). Derived fields (temperature by the lapse rate,
// climate class) are not projected; lake surfaces are inherited per global cell. The
// projection reads one parent cell per child, so it adds no halo. Planet-free.

#include <cstdint>
#include <vector>

#include "simulations/landscape/refinement_grid.hpp"

namespace landscape
{
// The projection's bounds over a level window's domain (domainsize() values each, row-major;
// S2.3's rework, plan section 3, "S2.3's projection weights and detail levels", item 1):
// per child a rate and the headrooms up and down in metres (+infinity: unbounded).
struct ProjectionBounds
{
    std::vector<float> rate;    // the share of the common shift a child takes (0: it stays)
    std::vector<float> up_m;    // how far the child may rise
    std::vector<float> down_m;  // how far the child may sink
};

// The bounds of a level's children. With s the sea level, s+ and s- the next floats above
// and below it, S a wet cell's lake surface and S- the next float below it:
// - sea (sea[cell] != 0): rate 1, up max(0, s- - F*), down unbounded (it stays sea);
// - wet (lakesurface_m[cell] > 0, not sea, a channel or not): rate 1, up max(0, S- - F*),
//   down unbounded (it stays below its surface);
// - land: rate channelrate on channel beds (channel[cell] != 0;
//   refinement.projection_weight_channel, 0 by default) and 1 elsewhere, up unbounded, down
//   max(0, F* - max(s+, L+)), L the highest lake surface of its wet 8-neighbours in the
//   domain and of `ownlake_m[cell]` (the lake surface of its own global cell when that is a
//   lake cell, 0 otherwise; empty: none), L+ the next float above it (the dry shore stays
//   above the water beside it, regional_water.hpp's step 5).
// F* is `elevation_m`. Throws std::invalid_argument when a field is not domainsize() long
// (ownlake_m may be empty) or channelrate is negative or not finite.
ProjectionBounds projectionbounds(const LevelWindow& window, const std::vector<std::uint8_t>& channel,
                                  const std::vector<std::uint8_t>& sea, const std::vector<float>& lakesurface_m,
                                  const std::vector<float>& elevation_m, float sealevel, double channelrate,
                                  const std::vector<float>& ownlake_m = {});

// What a projection did. The statistics are taken over the parent cells of the level
// window's kept box (whole parent cells: columns / factor by rows / factor), so they do not
// depend on the halo; the projection itself acts on every parent cell of the domain.
struct ProjectionReport
{
    int parents = 0;                    // parent cells of the kept box
    double maxparentmeanerror_m = 0.0;  // largest |mean_C(c)(F) - v_c| after the projection, F read back as float
    int bedsmovedparents = 0;           // parents whose rated children's headroom was exhausted, the rest on the beds
    int forcedparents = 0;              // parents with neither: the rest on every child, headrooms ignored
    double maxforcedmove_m = 0.0;       // the largest such forced shift of a child (kept parents)
    double maxlambda_m = 0.0;           // largest |lambda_c|
    double maxlambdastep_m = 0.0;       // largest |lambda_c - lambda_c'| over kept parents adjacent east-west or north-south
    // The border step a cell sees: d, the shift the projection applied to a child. Over the
    // pairs of edge-adjacent child cells (east-west or north-south) in different parents,
    // both in the kept box and both moved (rate > 0 or a shift applied), |d - d'|: the
    // largest, and the 99th percentile by nearest rank (the sorted value at index
    // ceil(0.99 n) - 1 of n pairs); both 0 when there is no pair. Adjacency is in the kept
    // box's local indices, so a pair across the date line inside the kept box counts; no pair
    // crosses the kept box's own edges. Its size: a one-cell-wide channel line incised by I
    // across a parent of factor^2 children gives lambda = I / factor and, with the line at
    // rate 0, the land beside it lifts by I / (factor - 1) (I / 3 at factor 4, I at factor
    // 2); a single incised cell lifts the others by I / (factor^2 - 1). A neighbouring parent
    // without the line has d near 0, so the step is that lift.
    double maxborderstep_m = 0.0;
    double borderstep99_m = 0.0;
};

// Projects `field` (the level window's domainsize() candidates F*, row-major) onto the
// parent values in place: for every parent cell c of the level window's domain, v_c is
// parentfield at the parent window's domain cell (parent.localcolumn(column of c),
// parent.localrow(row of c)), lambda_c = v_c - mean(F*) and the children's shifts sum to
// factor^2 lambda_c: the same shift delta x rate on every rated child (rate > 0), capped at
// its headroom in lambda's direction, delta found by iteration inside the parent (children
// in row-major absolute order, double, at most factor^2 rounds: delta = remainder / the
// uncapped children's rates, children whose headroom is less take it and leave). When the
// headrooms are exhausted the remainder goes uniformly to the unrated children (the beds;
// bedsmovedparents), or with none to every child ignoring the headrooms (forcedparents,
// maxforcedmove_m): the mean is exact always. F = float(double(F*) + shift). Requires
// level.level >= 1, parent.level == level.level - 1, the same world grid and factor, every
// parent cell of the level's domain inside the parent window's kept box (the halo chain of
// refinement_grid.hpp guarantees it), and field, the bounds and parentfield of their
// windows' domain sizes; throws std::invalid_argument otherwise, and when a rate is
// negative or not finite or a headroom negative or NaN.
ProjectionReport projectlevel(const LevelWindow& level, std::vector<float>& field, const ProjectionBounds& bounds,
                              const LevelWindow& parent, const std::vector<float>& parentfield);

// The largest |mean_C(c)(field) - v_c| over the parent cells of the level window's kept
// box, the mean unweighted and accumulated in double in the order above, v_c read as in
// projectlevel; 0 to float rounding after projectlevel (the S2.2 test: below 1e-3 m on
// 3 km of relief). The same requirements and exceptions as projectlevel, without weights.
// Implemented in S2.2.
double levelparentmeanerror(const LevelWindow& level, const std::vector<float>& field, const LevelWindow& parent,
                            const std::vector<float>& parentfield);

// --- Output at a requested resolution (section 1, "Output at any resolution"). ---

// How a resampled output was made.
struct ResampleReport
{
    int sourcelevel = 0;            // the level whose kept box was averaged
    double metresperpixel = 0.0;    // as requested
    double sourcecell_m = 0.0;      // the source level's north-south cell size, pi R / levelheight()
    double pixelcells = 0.0;        // n = metresperpixel / sourcecell_m, source cells per pixel edge, as computed
                                    // (an aligned output's pixels use round(n))
    bool aligned = false;           // n is whole and divides factor^sourcelevel: pixels tile whole global cells
    bool resampled() const { return !aligned; }
};

struct ResampledField
{
    int width = 0;                  // pixels, west to east from the kept box's west edge
    int height = 0;                 // pixels, north to south from the kept box's north edge
    std::vector<float> values;      // width x height, row-major, north row first
    ResampleReport report;
};

// The box-average resample of a chain's fields to `metresperpixel` (> 0, the north-south
// ground size of a pixel; east-west the pixel spans the same number of source cells, as the
// level grids are square in index space). windows and fields are the levels 0 .. L of a
// halo chain (index = level, fields[l] of windows[l].domainsize()); only the kept boxes
// are read, and every level's kept box covers the target's kept box. The source level is
// the coarsest level l of the chain whose cell size pi R / levelheight(l) is at most
// metresperpixel (1 + 1e-9); n = metresperpixel / that cell size. The output has
// floor(target columns at the source level / n) by floor(target rows / n) pixels, anchored
// at the target's kept box's north-west corner (a partial pixel at the east and south
// edges is dropped), where the target's kept box at level l is windows.back()'s kept box
// scaled by factor^(l - L). Aligned when n is whole (|n - round(n)| <= 1e-9 n) and round(n)
// divides factor^sourcelevel: the pixels then tile whole global cells (no partial pixel at
// any edge, since the kept box is whole global cells), and pixel (p, q) is the unweighted
// mean, accumulated in double in row-major order, of the round(n) x round(n) source cells
// p round(n) .. (p + 1) round(n) - 1 by q round(n) .. (q + 1) round(n) - 1 of the target's
// kept box, so the means of the level grids, and the global values, are preserved exactly.
// Otherwise resampled (a whole n that does not divide factor^sourcelevel included): pixel
// (p, q) covers the source-cell interval [p n, (p + 1) n) x [q n, (q + 1) n) of the
// target's kept box and is the sum of each overlapped source cell times its overlap
// lengths in cell units, divided by n^2 (double; the mean over the pixel's box, not exact against the level grids). Throws std::invalid_argument when
// metresperpixel or radius_m is not positive and finite, the chain is empty, malformed
// (element l not at level l, sizes wrong) or has no level fine enough (metresperpixel
// below the finest level's cell size: no output finer than the data), the target's kept box
// is not whole cells of the source level or the source level's kept box does not cover it,
// or the output would have no pixel. Implemented in S2.2.
ResampledField resampleto(const std::vector<LevelWindow>& windows, const std::vector<std::vector<float>>& fields,
                          double radius_m, double metresperpixel);
}
