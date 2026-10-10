#pragma once

// Lakes and present-day ice of the regional pass (docs/architecture/erosion-hydrology-plan.md,
// Phase E, step E3, and Phase R: R0 moved the rules here from regional_pass.cpp
// unchanged; R1 fixes the channel cells inside lakes; R4 owns the shores and the ice
// edges). Local rules only (no flood fill but the refinement loop's, bounded by a lake cell's
// 3 x 3 global cells, S2.8), each a pure function of the global fields, the seed and the fine
// cell's absolute position. Planet-free.

#include <cstdint>
#include <vector>

#include "simulations/landscape/channel_network.hpp"
#include "simulations/landscape/regional_grid.hpp"
#include "simulations/landscape/regional_types.hpp"
#include "simulations/landscape/upsampling.hpp"

namespace landscape
{
// Domain fields (domainfinewidth() x domainfineheight(), row-major).
struct FineWater
{
    std::vector<float> lakesurface_m;   // 0 where no lake water
    std::vector<float> icethickness_m;
    // S2.8's rim raise (LakeInheritance::flood only, else empty): per domain cell the lift of a
    // dry cell more than D below a wet 8-neighbour's surface raised to just above it, 0 elsewhere.
    std::vector<float> rimraise_m;
};

// Phase E's rules (kept while shores.shorenoise is 0): lake water stands on a fine
// cell that is not sea, whose global cell carries a lake surface and whose
// `elevation_m` lies strictly below it; ice is the upsampled thickness, at least 0,
// on fine cells whose global cell or one of its 8 neighbours has ice. R1: a channel
// cell of a lake cell carries the lake where its bed lies at or below the higher of
// the lake surface and the global cell's filled elevation (a closed lake can lie below
// its spill level, where the bed is pinned), and its `elevation_m` is lowered to at
// most the surface, so no lake surface lies below its cell; such a lake bed may lie
// below its outlet's sill, so the elevation along the channel directions rises from
// the lake to the sill (by up to 56 m in window C of erosion run 9), and a lowered
// cell lies below the sea level where its lake's surface does.
// R4 (shores.shorenoise > 0; D = shores.terrain_m, z the solve's elevation):
// - The indicator: -m on lake cells (lake surface above zero, not sea), +m elsewhere, m a
//   hash of the cell's indices and the seed in 0.101 .. 0.900, interpolated (bicubic) at
//   the fully warped fine positions.
// - The zone: a fine cell (not sea) with a lake cell among its global cell and the 8
//   neighbours takes the surface S of the one whose term in that interpolation (stencil
//   weight times magnitude) is the strongest at it (ties: the smaller global index).
// - The sum: indicator + shorenoise x RegionalNoise (a salt of its own, persistence
//   coast.noisepersistence) + clamp((z - S) / D, -w, 1), w the lake cells' mask (1 on lake
//   cells, 0 elsewhere) interpolated at the warped position and clamped to 0 .. 1; water
//   where the sum is below zero and z < S + D. So the terrain decides far from the surface
//   (high ground is land; deep ground is water where lake cells dominate, not in the valley
//   beyond an outlet) and the indicator near it; w is continuous across block edges.
// - Steps: a wet cell beside a wet cell of a lower surface is dry (a dam), so adjacent
//   water carries one surface.
// - Channel cells: R1's rule (wet at or below max(surface, filled), lowered to the surface)
//   is the shorenoise-0 rule; R4 replaces it. A channel cell inside the water above is wet
//   like any cell; a channel cell of a lake cell outside it is wet when z < S + D and an
//   8-neighbour that is not a channel is wet, or when z < S, unless a wet 8-neighbour has
//   another surface; otherwise (a bed at or above the surface away from the water, or the
//   outlet between two lakes) it stays a river. Of two adjacent channel cells this rule
//   wets at different surfaces, the higher stays a river too.
// - Beds: a wet cell with z >= S is lowered to S minus a depth growing from 0 where the sum
//   is 0 to min(the global lake depth (S minus the zone's base elevation), D) where the sum
//   reaches minus the zone's magnitude, and to at least the next float below S; a wet cell
//   with z < S keeps z. So the lowering is below 2 D. With `globaldepth` (the refinement
//   loop, multiscale plan S2.3) the depth is the global lake depth S minus the zone's base
//   elevation uncapped, so a deep global lake's fine bed reaches it; false (the single
//   pass) keeps the cap.
// - Shores: a dry fine cell (not sea, not a channel) at or below the level L, the highest
//   surface of its wet 8-neighbours and of its own global cell if that is a lake cell, and
//   at most D below it, is raised to L + (L - z) (mirrored about L), at least the next float
//   above L, so by at most 2 D (a dry cell of a lake cell with no water beside it is raised
//   to the mirror of its own lake's surface the same way); a cell more than D below L is left:
//   the deep bed of a pit lake outside the indicator's outline, left dry rather than dammed by
//   hundreds of metres.
// Ice: the indicator with the ice cells (thickness above zero) negative, a salt of its own
// and no terrain term, on the fine cells of ice cells and their neighbours; the thickness is
// the thickness upsampled at the warped position, at least 1 m; a fine sea cell carries ice
// only where its global cell is sea and iced.
//
// The lake flood of the refinement loop (multiscale plan S2.8, §3 "S2.7: smooth before exact",
// item 3; R4's rules only, shorenoise > 0; inactive in the single pass). A lake's water is its
// core and every cell below its surface S 4-connected to the core through cells below S of the
// same zone (so inside the lake cell's 3 x 3 global cells), not sea, where the lake cells
// dominate (their mask interpolated at the warped position, w >= 0.5: the barrier against a lower
// non-lake neighbour, an outlet valley or a low plain, lies on that warped contour, not on the
// cell boundary) or in a global cell above S (the lake's rim), and below S on the start surface
// when one is given. Reach: a component inside one lake cell's 3 x 3 zone (2 global cells beyond
// a kept box at level 1; 1 parent cell + 3 level cells at level >= 2, which is at most 2
// global cells too), read on the solved elevation, so a kept cell depends on it up to 2 global
// cells farther than on the solve: the loop's halo adds them (refinement_levels.hpp). The
// core: at level 1 (no parent masks) the water of the sum above, so the indicator shapes the core
// only; at level >= 2 the cells below S whose parent is wet with the zone's surface, the flood
// kept to the cells whose parent or a parent 8-neighbour is wet (refined within one parent cell).
// The steps from the dam rule on (2 to 5) then act on that water as on the sum's; step 5 raises
// the dry rest, and a dry cell more than D below a wet 8-neighbour's surface (left by step 5
// without the flood) to the next float above that surface. The connected set does not depend on the order of the cells.
// S3.4: the lake-sea contact at a sill: a cell where the sea gate is open (regional_edges.hpp,
// seagate with the lake weight; at level >= 2 also parentseagate) and the elevation lies between
// the sea level and S is passable whatever the lake mask says, so the flood crosses into a
// sea-gated block. The rim cap (rimcap_m, the level's noise band): a dry cell (not sea, not a
// channel) of the lake's zone, below S on the start surface (at level >= 2 within the parent
// reach), beside its lake's water and not beside water of a lower surface, that step 5 would raise
// by more than D and more than its cap joins the water instead, in rounds until none joins; a
// deeper dry cell beyond those bounds keeps the raise to the next float above the surface (rimraise_m).
struct LakeInheritance
{
    bool flood = false;
    // Level >= 2: one value per domain cell, the parent's lake surface (0: the parent is dry)
    // and 1 where the parent or one of its 8 neighbours is wet. Null at level 1.
    const std::vector<float>* parentsurface_m = nullptr;
    const std::vector<uint8_t>* parentreach = nullptr;
    // The level's start surface (before the solve, its parent means kept; regionallevelsteps
    // sets it): the flood enters only cells below S on it too, so a parent above S keeps a dry
    // child (the projection stays feasible where the solve lowered the land below S).
    const std::vector<float>* startsurface_m = nullptr;
    // S3.4: the rim raise's cap per domain cell (lakeslevel sets it to LevelState::noiseband_m, the
    // level's noise band); null: no cap.
    const std::vector<float>* rimcap_m = nullptr;
    // S3.4: at level >= 2 the parent's sea gate per domain cell (RegionalDetails::parentseagate; lakeslevel
    // sets it): the sill rule passes only where it is 1. Null at level 1.
    const std::vector<uint8_t>* parentseagate = nullptr;
};

FineWater regionalwater(const RegionalGlobalInputs& inputs, const RegionalWindow& window, const FineWarp& warp,
                        const RegionalParameters& parameters, const std::vector<uint8_t>& sea,
                        const FineChannels& channels, std::vector<float>& elevation_m, bool globaldepth = false,
                        const LakeInheritance& inheritance = LakeInheritance{});

struct LevelContext;  // regional_pass.hpp (the hook interface, multiscale plan S3.0)
struct LevelState;

// LevelStep::lakes: state.water = regionalwater on state.elevation_m with details.lakes (its
// start surface state.surface_m when the flood is on) and details.globallakedepth. S3.2: every
// land channel cell is read at its bed (state.bed_m) instead of its cell mean; one the lakes leave
// dry takes its mean back, a wet one keeps the lakes' bed.
void lakeslevel(const LevelContext& context, LevelState& state);
}
