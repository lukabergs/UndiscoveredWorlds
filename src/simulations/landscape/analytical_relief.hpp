#pragma once

// Analytical relief of the regional pass (docs/architecture/erosion-hydrology-plan.md,
// Phase E, step E2c, and section 6, "The solve"): the transient solution of the
// linear stream power law (n = 1, uplift 0 unless uplift_m_myr is given) on a fine patch whose pinned cells keep
// their elevation (Steer 2021, Earth Surface Dynamics 9, 1239; Tzathas, Gailleton,
// Steer and Cordonnier 2024: a non-flat starting terrain, fixed boundary heights, a
// prescribed network, elapsed time and a threshold-slope hillslope). With n = 1 and
// no uplift the solution travels upstream along the flow paths: the response time
// of a cell is tau(x) = sum over the links from x down to its first pinned cell of
// L / (K Q^m); after an elapsed time t a cell takes the start elevation of the point
// at response time tau(x) - t downstream of it (interpolated along the link), or
// the pinned cell's elevation when tau(x) <= t, and never more than its own start
// elevation. Then every slope steeper than the threshold slope is lowered to it,
// from the pinned cells upwards. Receivers come from the result (damped, see
// AnalyticalReliefParameters::damping), so routing and solution are repeated, always
// maxiterations times. Units: elevation m, K in m/Myr per (km^3/yr)^m per unit slope
// (landscape.stream_power.erodibility), discharge m^3/yr inside. Serial in a
// canonical order (each step over the cells in flow order with ties by the absolute
// key of AnalyticalReliefInput), so the result is bit-identical for any thread count.
// Seams: the result at a cell depends only on its sub-basin (the region of unpinned
// cells connected to it by the moves the routing allows: orthogonal steps and
// diagonal steps not between two pinned corner cells), the pinned cells bordering
// it, the inputs and the keys' order there. No receiver, flood step or depression's
// spill path leaves the sub-basin in any pass, so the passes do not widen it, and
// two windows that both hold a sub-basin and its border give it bit-identical
// results. Sub-basins are unbounded (a forest of pinned lines does not cut the plane);
// with AnalyticalReliefParameters::catchmentreach r (S3.9) each cell's result reads only the
// cells within 4r global cells (the bound is stated there). Table-driven (T8 W2-W2a): the
// patch is a mesh::CellGrid (regionalwindowgrid's window, or a hand-built planar patch), its
// neighbour slots the D8 codes, -1 beyond the patch. Planet-free.

#include <cstdint>
#include <vector>

#include "core/mesh/cell_grid.hpp"
#include "simulations/landscape/regional_grid.hpp"

namespace landscape
{
struct AnalyticalReliefParameters
{
    double erodibility = 750.0;     // K, as the global kernel's
    double areaexponent = 0.5;      // m
    double elapsed_myr = 1.0;       // t
    double thresholdslope = 0.6;    // hillslope limit (rise over run)
    int maxiterations = 20;         // routing and solution passes, always all of them (an early
                                    // stop would depend on cells anywhere in the window)
    // Weight of the new solution in the surface the next routing reads, an exponential
    // moving average (Tzathas et al. 2024, section 5.1, which damps the fixed point's
    // oscillations this way without giving a value): z_route(i) = damping z_solved(i)
    // + (1 - damping) z_route(i - 1), z_route(0) = the start. 0 < damping <= 1; 1 is
    // no damping. The solution itself is always solved from the start surface.
    double damping = 0.5;
    // Multiscale plan S3.9 (§3, "Seamless solve: a bounded catchment", with the review's local
    // labelling; regional.catchment_reach_global_cells): r >= 1 bounds every catchment. Labels, once,
    // on the start surface: one Priority-Flood per global cell g (AnalyticalReliefInput::globalcell),
    // seeded at the pinned cells inside g's block and confined to the unpinned cells whose global cell
    // lies within r rings of g (R: 8-connected hops on the level-0 grid's table, AnalyticalReliefInput::
    // globalgrid; on lat-lon the Chebyshev square, columns wrapped, rows stopping at the poles); each unpinned cell
    // takes the seed with the smallest (filled level, hashed absolute key of the seed) over the floods
    // reaching it: its catchment's outlet. Every pass routes inside the label only: steepest descent to
    // a neighbour of the same label (or to the outlet), else the receiver of a flood that steps only
    // into cells of its seed's label; every queue of the bounded solve ties by a hashed key, not the
    // row-major one. A cell no flood reaches, or whose path in its winning flood leaves its label,
    // keeps its start elevation (AnalyticalReliefResult::unreached). Bound: a label reads the cells
    // within 2r global cells (the floods of the global cells within r, each confined to r around its
    // own); a catchment lies within r of its outlet's global cell, so within 2r of any of its cells,
    // and which cells belong to it reads their labels: a cell's result reads the cells within 4r
    // global cells. 0: the solve up to S3.8, bit for bit (its sub-basins are unbounded).
    int catchmentreach = 0;
};

struct AnalyticalReliefInput
{
    std::vector<float> start_m;          // starting surface (the upsampled terrain with noise)
    std::vector<uint8_t> pinned;         // 1: elevation fixed (channels, sea, the patch's outer ring)
    std::vector<float> runoff_m_yr;      // per cell; discharge = sum of runoff x area over the catchment
    // The absolute key of each cell (the window grid's globalkey: on lat-lon the fine index
    // I + J x kW, regional_grid.hpp), supplied by the regional pass, distinct per cell (else a
    // throw); empty: the local index. Every tie of the
    // flood's order (equal levels) is broken by it, never by the local index, so two
    // windows that cross the date line differently order the same cells the same way;
    // the flat resolution's ties (S3.12) go by a hash of it, the other ties (steepest
    // descent, donors) by direction code.
    std::vector<std::uint64_t> key;
    // Erosion plan Phase R2 (the R6 contract, implemented by R9): the rock uplift rate per cell,
    // m/Myr, at least 0; empty: no uplift (the solution above exactly). With it, a cell
    // within the knickpoint wave's reach (tau <= t) takes min(z0(x), the outlet's elevation
    // + U(x) x tau(x)) instead of the outlet's elevation: the characteristic solution's
    // uplift term (Steer 2021; Tzathas et al. 2024, Eqs. 8 and 9 with u > 0), which in the
    // steady state is the profile of slope U / (K Q^m) from the pinned cell up to where it
    // meets the start surface; the cell's own rate stands for the rate along its path (an
    // approximation the plan records). Cells beyond the reach (tau > t) advect the start
    // surface as without uplift: the regional pass works in the frame of the global mean
    // elevation, where uplift appears as the pinned channels' incision (regional_valleys.hpp),
    // not as a lift of the land. The threshold slope applies after it as before.
    std::vector<float> uplift_m_myr;
    // S3.9 (AnalyticalReliefParameters::catchmentreach >= 1 only, read nowhere else): the global
    // cell of each cell, its level-0 ancestor (an id of globalgrid, RegionalWindow::globalcell),
    // and the level-0 grid whose table counts the reach's rings (not null).
    std::vector<mesh::CellId> globalcell;
    const mesh::CellGrid* globalgrid = nullptr;
};

struct AnalyticalReliefResult
{
    std::vector<float> elevation_m;
    // 1..8 towards the receiver of the last pass, 0 on pinned cells. A direction may
    // point uphill on elevation_m: a depression's flood links under the clamp to the
    // start, and the lag of the damped routing surface behind the last solution.
    std::vector<uint8_t> direction;
    // Accumulated over each unpinned catchment; 0 on pinned cells (the regional pass
    // adds the channel network's discharge there).
    std::vector<double> discharge_m3_yr;
    int iterations = 0;                  // routing passes run: maxiterations
    bool converged = false;              // the last two passes' receivers are equal (false with one pass)
    int changedreceivers = 0;            // cells whose receiver differs between the last two passes (0 with one pass)
    // The pinned cell each cell drains to on the last pass's receivers (itself on a pinned cell);
    // with catchmentreach >= 1 an unpinned cell's catchment label.
    std::vector<int32_t> outlet;
    // catchmentreach >= 1: 1 on the unpinned cells no catchment reached or cut from their outlet,
    // kept at their start elevation (direction 0, discharge 0); empty with catchmentreach 0.
    // unreachedcells counts them.
    std::vector<uint8_t> unreached;
    int unreachedcells = 0;
};

// Throws std::invalid_argument when the input sizes differ from the grid's cell count (key and
// uplift_m_myr: empty or one per cell; with catchmentreach >= 1 globalcell one per cell, each an id
// of a non-null globalgrid), two keys are equal, a parameter is out of range (K <= 0,
// m < 0, t < 0 or not finite, threshold slope <= 0, no iteration, damping outside
// 0 < damping <= 1, catchmentreach < 0), or no cell of a non-empty patch is pinned. No cell ends above its start elevation: the
// advected value is min(z0(x), the rule's value), since with no uplift the erosion
// term K A^m |grad z| never raises a cell (it would where the receivers run uphill on
// the start surface: in a depression, or after the routing changed). Receivers: steepest
// descent to a neighbour lower on both the surface and the filled surface of a
// Priority-Flood seeded at the pinned cells (on a surface without depressions: to
// any lower neighbour), otherwise (multiscale plan S3.12) the flat resolution of
// Barnes, Lehman and Mulla 2014 on the flood's pits and flats, with metric (Dijkstra)
// distances and ties by a hash of the key: steepest descent of an increment that falls
// towards the flat's low edges and away from higher terrain; a cell it does not reach
// takes the neighbour from which that flood reached it. Every unpinned cell
// has a path to a pinned cell and no receiver loop exists (steepest descent to any
// lower neighbour would loop between a pit and the rim cell its flood came from);
// the surface is not filled. Neither rule takes a diagonal step whose two corner
// cells (orthogonally adjacent to both ends) are pinned, so no catchment crosses an
// 8-connected channel.
AnalyticalReliefResult solveanalyticalrelief(const mesh::CellGrid& grid, const AnalyticalReliefInput& input,
                                             const AnalyticalReliefParameters& parameters);
}
