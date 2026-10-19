#pragma once

// The coast's ruggedness from the world's fields and the island sink (docs/architecture/erosion-hydrology-plan.md,
// Phase R2: the R6 contract, which R7 implements and tunes; the diagnosis is in section 6,
// 'The user's review of regional runs 5 to 8', item 1). Every value is a pure function of
// the global fields, the seed and the fine cell's absolute position. Called by
// regional_surface.cpp (the Surfacer and the land-share rule's widened box). Planet-free.

#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

#include "simulations/landscape/regional_grid.hpp"
#include "simulations/landscape/regional_types.hpp"
#include "simulations/landscape/spherical_noise.hpp"

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
// function of the mask within m cells of the block and the global fields (every block is
// decided on the mask as it came in); regional_surface.cpp widens the land-share rule's
// margin for the ring (its margin comment).
void sinksmallislands(const RegionalGlobalInputs& inputs, const RegionalWindow& box, const RegionalParameters& parameters,
                      const std::vector<uint8_t>& gate, std::vector<uint8_t>& mask);

// --- Multiscale plan S3.15 (§3 "The user's verdict on runs 41 and 42", S3.15): the coast. ---

// The shoreline's displacement at one level (S3.15 a): a fractal vector noise of two components
// (gradient noise on the unit sphere, spherical_noise.hpp, salts of their own) whose octave o has
// the wavelength of one global north-south spacing divided by 2^o and the standard deviation
// coast.warpcells x coast.warppersistence^o global cells per component (the octave's sum divided
// by octavedeviation), with the octaves of the level's band only: o with 2^(o + 1) <= scale (at
// least two fine cells) and 2^(o + 1) > parentscale (those the parent level does not draw; 1 at
// level 1 and in the single pass, so the levels of a hierarchy add up to the octaves of a single
// pass at the finest scale). The vector u of the band's octaves is shortened to
// b tanh(|u| / b) u / |u| with b = 2 x the band's standard deviation per component
// (coast.warpcells x sqrt(sum of warppersistence^(2 o) over the band); a method constant), so
// |displacement| < bound() = b global cells and a small displacement is kept nearly whole. Inactive
// (0 everywhere) when coast.warpcells <= 0 or the band holds no octave. A pure function of the fine
// cell's absolute position and the seed; dx eastwards, dy southwards, in global cells.
class ShorelineWarp
{
public:
    ShorelineWarp(int worldwidth, int worldheight, int scale, int parentscale, std::uint64_t seed,
                  const CoastParameters& coast);
    bool active() const { return active_; }
    double bound() const { return bound_; }
    // At the fine cell (I, J) of the kW x kH grid (I taken modulo kW).
    void at(std::int64_t finecolumn, std::int64_t finerow, double& dx, double& dy) const;
    // The same at a unit vector (a cube cell's centre; at() is this at the lat-lon centre), in level-0 cells.
    void atpoint(double px, double py, double pz, double& dx, double& dy) const;
    // The fine cells beyond its own that a displaced sample reads: ceil(bound() x scale) + 1
    // (the bilinear stencil); 0 when inactive.
    int reachcells() const;

private:
    bool active_ = false;
    int scale_ = 1;
    int finewidth_ = 0;
    int fineheight_ = 0;
    double amplitude_ = 0.0;
    double bound_ = 0.0;
    std::optional<FractalNoiseField> x_, y_;
};

// A field on `box`'s domain (domainfinewidth() x domainfineheight(), row-major) sampled at every
// cell's displaced position p + warp(p): bilinear between the fine cell centres, the stencil's
// columns wrapped when the domain spans the whole turn and clamped to the domain otherwise, rows
// clamped to the domain (window-independent wherever the stencil lies inside the window-independent
// part of the field: the callers' margins count warp.reachcells()). The field itself where the warp
// is inactive.
std::vector<double> displacedsample(const std::vector<double>& field, const RegionalWindow& box, const ShorelineWarp& warp);
std::vector<double> displacedsample(const std::vector<float>& field, const RegionalWindow& box, const ShorelineWarp& warp);

// The sea pockets (S3.15 b): one value per fine cell of `box`'s domain, 1 on a sea cell of `sea`
// (1 = sea) of any global cell that does not connect through fine sea (4-connected) to a core sea
// cell (`core` 1 and sea) within `reach` global cells of its own global cell (Chebyshev), the path
// inside those (2 reach + 1)^2 global cells; 0 elsewhere. `core` (the callers: the cells of global
// sea cells whose decision without the coast noise and the displacement lies at or below the sea
// level) holds one value per fine cell; empty: the sea cells of global sea cells. Decided per global
// cell on `sea` as it comes, so the result is a pure function of `sea` and `core` within `reach`
// global cells of each block (blocks beyond the box do not exist; columns wrap around a whole turn).
// The path never enters a global lake cell (lake surface above zero, not sea), and a lake cell's own
// sea is never a pocket. All zero when reach <= 0.
std::vector<uint8_t> seapockets(const RegionalGlobalInputs& inputs, const RegionalWindow& box, int reach,
                                const std::vector<uint8_t>& sea, const std::vector<uint8_t>& core = {});

// The basins of filled pockets (S3.15 b): one value per fine cell of `box`'s domain, 1 on the cells
// reached from a global cell's `filled` cells through 4-connected land cells (`sea` 0) that are
// filled or whose `surface` lies below `level` (the callers: the next float above the sea level),
// inside the global cells within `reach` of that global cell (the filled cells included); 0 elsewhere.
// A cell of a global lake cell, or of a global cell beside one (8-neighbours) lying below that lake's
// surface (the highest such), is never a basin cell: it is the lake's bed, which the lake flood fills
// after the coast (the merged run 47 lifted it out of the lake at 58-59 E, 4-5.5 N). The global lake
// mask is absolute, so the rule stays seamless.
// Decided on the inputs as they come (a pure function of them within `reach` global cells).
std::vector<uint8_t> pocketbasins(const RegionalGlobalInputs& inputs, const RegionalWindow& box, int reach,
                                  const std::vector<uint8_t>& sea, const std::vector<uint8_t>& filled,
                                  const std::vector<float>& surface, float level);

// The pockets of a mask in a rectangle of a domain (the kept box): its cells with `pockets` 1 and
// their 4-connected components inside the rectangle.
struct PocketCount
{
    int pockets = 0;
    int pocketcells = 0;
};
PocketCount countpockets(const std::vector<uint8_t>& filled, int domainwidth, int x0, int y0, int width, int height);

// --- The cube forms (T8 W2-W2c, sphere-mesh.md section 8, "Coast rules"; lat-lon boxes take the
// forms above). `box` is a cube window grid refining inputs.grid (a cube level 0; k = its face
// cells over level 0's), whole level-0 blocks (level_blocks.hpp); every per-cell vector holds one
// value per box cell in its local order. ---

// The field sampled bilinearly (mesh::cubebilinearstencil of the box's level) at each cell's
// centre displaced by the warp (warpedpoint, the displacement in level-0 cells), a node outside
// the box taking the box cell at the smallest ring distance (ties to the smallest key).
std::vector<double> displacedsample(const std::vector<double>& field, const mesh::CellGrid& box, const ShorelineWarp& warp,
                                    int levelzerofacecells);
std::vector<double> displacedsample(const std::vector<float>& field, const mesh::CellGrid& box, const ShorelineWarp& warp,
                                    int levelzerofacecells);
// sinksmallislands: a global sea cell's block whose gated land the mask holds is sunk unless a
// one-cell or mountain island of the world lies in its 3 x 3 (the world's table), when every cell
// at exactly `margin` rings from the block (the box level's closed-form neighbours) is sea; a ring
// cell missing from the box counts as land.
void sinksmallislands(const RegionalGlobalInputs& inputs, const mesh::CellGrid& box, const RegionalParameters& parameters,
                      const std::vector<uint8_t>& gate, std::vector<uint8_t>& mask);
// seapockets and pocketbasins: as above with the region of a block the blocks whose level-0 cells
// lie within `reach` rings of its own, the searches spreading through the cardinal neighbours
// (slots 1, 3, 5, 7) on the box's table.
std::vector<uint8_t> seapockets(const RegionalGlobalInputs& inputs, const mesh::CellGrid& box, int reach,
                                const std::vector<uint8_t>& sea, const std::vector<uint8_t>& core = {});
std::vector<uint8_t> pocketbasins(const RegionalGlobalInputs& inputs, const mesh::CellGrid& box, int reach,
                                  const std::vector<uint8_t>& sea, const std::vector<uint8_t>& filled,
                                  const std::vector<float>& surface, float level);
// countpockets: the 4-connected components (cardinal neighbours on the grid's table) of the filled
// cells among those `kept` flags (null: every cell).
PocketCount countpockets(const mesh::CellGrid& grid, const std::vector<uint8_t>& filled, const std::vector<uint8_t>* kept);

}
