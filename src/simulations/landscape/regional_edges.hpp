#pragma once

// The sea gate and the straight-run breaker of the regional coast (docs/architecture/erosion-hydrology-plan.md,
// Phase R2: the R6 contract, which R8 implements; the diagnosis is in section 6, 'The
// user's review of regional runs 5 to 8', item 2). Every value is a pure function of
// the global fields, the seed and the fine cell's absolute position, so two windows agree
// on their overlap. Called by regional_surface.cpp on the land-share rule's widened box.
// Planet-free.

#include <cstdint>
#include <vector>

#include "simulations/landscape/regional_grid.hpp"
#include "simulations/landscape/regional_types.hpp"
#include "simulations/landscape/upsampling.hpp"

namespace landscape
{
// Phase E's eligibility of a global cell's fine cells for the sea (the 3 x 3 gate): 1
// where a global sea cell lies among the cell and its 8 neighbours and the cell is not a
// lake cell (lake surface above zero) unless the cell itself is sea. W x H.
std::vector<uint8_t> nearseacells(const RegionalGlobalInputs& inputs);

// The global cells towards which the decision surface is lifted (regional_surface.hpp) as
// 1.0, the others as 0.0. W x H. edges.smoothgate 0 to 2: the cells that may hold no fine
// sea under the 3 x 3 rule (nearseacells 0). edges.smoothgate n >= 3 (R8): the cells with no
// global sea cell within n - 1 cells (Chebyshev; columns wrap, rows beyond the poles do not
// exist), and lake cells. So the lift is full at the centre of every cell of Chebyshev
// distance n from the sea, which is where the smooth gate of reach n ends (seagate), and
// none within n - 1: the sea does not stop on the gate's edge (that edge runs through global
// cell centres) unless the decision lies more than coast.seadrop_m below the sea level there,
// and inside the reach the start surface, not the lift, sets the coast (R8, window B: with
// the 3 x 3 lift and reach 2 or 3 the coast stayed within one fine cell of the 183/184 block
// edge, the contour where the lift's ramp across global cell (184, 49) met its depth; with
// reach 3 and this lift it moved one to three cells off it).
std::vector<float> inlandcells(const RegionalGlobalInputs& inputs, const RegionalParameters& parameters);

// The lake cells' mask (1 where the lake surface is above zero and the cell is not sea, 0
// elsewhere) interpolated (upsample, bicubic) at the warped fine positions of `box`'s domain:
// regional_water.cpp's lake weight (the flood's half-mask barrier) and seagate's below.
std::vector<float> lakecellweight(const RegionalGlobalInputs& inputs, const RegionalWindow& box, const FineWarp& warp);

// The gate: one value per fine cell of `box`'s domain (domainfinewidth() x
// domainfineheight(), row-major), 1 where the cell may be sea. `warp` is the box's warp
// (regionalwarp; empty without one).
// - edges.smoothgate 0 (Phase E to R): nearseacells of the cell's own global cell, so the
//   gate is a step function of the global grid and cuts the sea on a block edge wherever a
//   global land cell with no sea in its 3 x 3 holds fine cells below the sea level (regional
//   run 6, global cell (184, 49): a straight coast of 29 fine cells on the 183/184 edge).
// - edges.smoothgate n >= 1 (R8, the reach in global cells): 1 where nearseacells of the own
//   global cell is 1, or where a global sea cell's centre lies within n global cells of the
//   cell's warped position in each axis (half-open, -n <= offset < n; rows clamped at the
//   poles as bicubic clamps them), the position gx = (absolute fine column + 0.5) / k + fade x
//   warp.dx, likewise gy, with the Surfacer's fade (from the unwarped base elevation and
//   warp.fade_m). n = 2 is the 4 x 4 stencil that `bicubic` reads (the support of the
//   interpolated sea indicator); the gate's edge runs through global cell centres displaced
//   by the warp, never along a block edge, and the gate contains the 3 x 3 rule (with a warp
//   of up to half a cell, n >= 2 contains it anyway). 0 wherever the own global cell is a
//   lake cell (lake surface above zero and not sea): the sea beside a lake would otherwise
//   reach into it. Not the sign of the interpolated indicator: the Catmull-Rom weights of the
//   stencil's outer rows and columns are negative, so that sign is positive only within about
//   one cell of a sea centre, narrower than the 3 x 3 rule (the west half of a land cell whose
//   only sea neighbour lies east would lose its gate), with stray positive patches beyond
//   corners; it cannot open global cell (184, 49). Where the decision surface lies below the
//   sea level the fine sea may reach into a land cell up to the gate's edge; inlandcells
//   lifts the decision fully at that edge, the land-share rule keeps a land cell's block at
//   least `landshare` land, and every gated land cell is kept above the sea level by the
//   clamp. R8 chose n = 3 (landscape.json's note): with n = 2 the far halves of cells two
//   cells from the sea stayed outside the gate, so their fine cells below the sea level stayed
//   land (window A: 158 such cells, 24 with n = 3).
// - `lakeweight` (multiscale plan S3.4, the refinement loop at level >= 2; null elsewhere): the
//   lake cells' mask (1 on lake cells, 0 elsewhere) interpolated at the box's warped positions,
//   one value per fine cell. A lake cell's block is then not excluded as a whole: a cell is 0
//   where the weight is at least 0.5, and a cell of a lake cell's block below that follows the
//   rule of the other cells (with smoothgate 0 the 3 x 3 rule on its own global cell), so the
//   gate's edge at a lake beside the sea is the warped half-mask contour, not the block edge.
std::vector<uint8_t> seagate(const RegionalGlobalInputs& inputs, const RegionalWindow& box, const FineWarp& warp,
                             const RegionalParameters& parameters, const std::vector<float>* lakeweight = nullptr);

// The cube forms of the two above (T8 W2-W2c, sphere-mesh.md section 8: gates by table neighbours
// and G footprints; lat-lon windows take the forms above), one value per cell of a cube window
// refining inputs.grid (a cube level 0), in its local order; `warp` the window's (regionalwarp's
// cube form, or empty). lakecellweight: the lake mask gathered by G (upsampling.hpp) at each cell's
// fully warped point (warpedpoint, regional_surface.hpp). seagate: the rules above with the cell's
// level-0 ancestor as its global cell (nearseacells and the 3 x 3 on the world's table), and for
// edges.smoothgate n >= 1 a sea cell of the world among the n-footprint (mesh::cubefootprint: the
// C-chart 2n x 2n cells around the point on its face) of the faded warped point (the fade from G of
// the base elevation at the cell's centre and warp.fade_m, as the Surfacer's).
std::vector<float> lakecellweight(const RegionalGlobalInputs& inputs, const mesh::CellGrid& window, const FineWarp& warp);
std::vector<uint8_t> seagate(const RegionalGlobalInputs& inputs, const mesh::CellGrid& window, const FineWarp& warp,
                             const RegionalParameters& parameters, const std::vector<float>* lakeweight = nullptr);

// The straight-run breaker (R8, after the legacy removestraights): no-op while
// edges.straightmax <= 0. `mask` is the box's mask (1 = sea; domainfinewidth() x
// domainfineheight()), `gate` seagate's on the same box. A straight shore run is a row or
// column of consecutive sea cells whose cells all have land on the same side (north, south,
// west or east), measured inside one global cell's block: every run is split at block
// edges (a run along a block edge inside a block is still up to k cells long). The runs
// are found on the mask as it comes; a side beyond the box or a pole is never land. A run
// longer than edges.straightmax gets chunks, walking from its first cell: a gap of 1 to
// straightmax cells (so two chunks never touch; at most up to the run's last cell) drawn from a hash of the seed and
// the piece's first cell's absolute fine position and the side, then a chunk whose size,
// edges.chunkmin to edges.chunkmax cells, and kind are drawn from a hash of the seed and the
// chunk's own start cell's absolute fine position and the side; the walk goes on after the
// chunk until the rest of the run is at most straightmax long. A chunk of n cells is a first
// row of n - n/3 cells along the run (clipped to the run) and a second row of n/3 cells
// centred beyond it, so its own edges stay short: raised (land), the run's cells and the row
// beyond them on the sea side; lowered (sea), the land side's row and the one beyond it.
// In a global land cell's block every chunk is land; in a sea cell's, land or sea by the
// draw, land only while the block keeps at most 1 - landshare land, else sea, and sea only
// where a lowered cell can lie inside the block on a gated cell, else land if the share
// allows it, else no chunk. Chunks are clipped to their block; a cell may become sea only
// where `gate` is 1; a global lake cell's block is never touched. The breaker reads its
// block and one fine cell beyond it (the land-on-one-side test) and writes only its block,
// so it is a pure function of the mask there, the absolute position and the seed (the
// land-share rule's margin counts that reach of one fine cell). The caller repeats the
// pinch and hole passes after it and counts the land shares on the result.
void breakstraightruns(const RegionalGlobalInputs& inputs, const RegionalWindow& box, const RegionalParameters& parameters,
                       const std::vector<uint8_t>& gate, std::vector<uint8_t>& mask);
// The cube form (sphere-mesh.md section 8, "Coast rules"): the runs along each block's face rows and
// columns (blocks never straddle a face), the cell beyond the block from the face's chart
// (mesh::cubechartkey; beyond the box or without a chart cell never land), the hashes of the seed
// and the cell's level key; one value per cell of a cube box refining inputs.grid.
void breakstraightruns(const RegionalGlobalInputs& inputs, const mesh::CellGrid& box, const RegionalParameters& parameters,
                       const std::vector<uint8_t>& gate, std::vector<uint8_t>& mask);
}
