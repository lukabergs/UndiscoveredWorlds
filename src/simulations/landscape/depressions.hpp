#pragma once

// Priority-Flood (Barnes et al. 2014) over a terrain on a SphereMetric, seeded
// from the sea cells and nothing else (a missing neighbour, beyond a pole or at a
// cube corner, is skipped; a torus metric wraps rows), or, on a grid without a sea cell, from its terminal sink (terminalsink):
// the filled surface, one id per depression, its spill, and the flood
// receiver of every land cell
// (docs/architecture/erosion-hydrology-plan.md, section 3, Kernels). Phase A: a
// depression's level is its spill elevation (every depression overflows); the
// water-balance level of Phase B lowers it.

#include <cstddef>
#include <cstdint>
#include <vector>

#include "simulations/landscape/sphere_metric.hpp"

namespace landscape
{
struct Depression
{
    int32_t id = 0;                 // 1-based; the depressionid of its cells
    int32_t spillcell = -1;         // the rim cell the flood entered from: the lowest rim cell, not itself flooded
    float spillelevation_m = 0.0f;  // the filled elevation of every cell of the depression
    double area_m2 = 0.0;           // area of the flooded cells
    double volume_m3 = 0.0;         // sum over the flooded cells of (filled - terrain) x area
    // The terminal basin of a grid without a sea cell: no spill (spillcell -1); every
    // other depression drains into it, so a lake in it is closed whatever its inflow.
    bool terminal = false;
};

struct Depressions
{
    std::vector<float> filledelevation;   // >= terrain; equals the terrain outside depressions and on sea cells
    std::vector<int32_t> depressionid;    // 0 none; else the id of the connected flooded region
                                          // whose filled elevation exceeds the terrain by more than minimumdepth
    // The neighbour each land cell was reached from in the flood order (Barnes'
    // flow directions): the chain of flood receivers from any land cell leads
    // through its depression's spill to a sea cell, or to the terminal sink on a grid
    // without a sea cell. Sea cells and the sink: -1.
    std::vector<int32_t> floodreceiver;
    std::vector<Depression> depressions;  // by id - 1
};

// `sea` marks sea cells (1), the seeds of the flood, which enter at their own
// elevation (the seabed): the caller's mask must put every cell at or below the
// sea level that touches a sea cell in the sea (the world ocean of oceanmask, or
// every below-sea cell), so that no land cell lies below a sea neighbour's seabed;
// an enclosed below-sea region is land and is filled from its rim; `elevation` is
// the terrain.
// A flooded cell is raised to the elevation of the cell it was reached from;
// cells raised by less than `minimumdepth_m` keep depressionid 0 (and their
// filledelevation is still the raised value). Deterministic: the priority queue
// breaks ties by cell index.
// Without a sea cell (an all-land grid) the flood starts from the terminal sink, a
// virtual outlet, so every other depression is found and spills, through the chain
// of flood receivers, towards it. The sink's own basin, the 8-connected cells
// below the pass level that hold the sink (the level at which the flood first has
// to raise a cell; the highest terrain when it never does), is raised to that level
// and is depression 1, terminal, whatever its depth; the other depressions follow
// as above. A grid of one elevation has no terminal basin.
Depressions filldepressions(const SphereMetric& metric, const std::vector<float>& elevation,
                            const std::vector<uint8_t>& sea, float minimumdepth_m);

// The terminal sink of a grid without a sea cell: its lowest cell (the lowest
// index on a tie); -1 when `sea` marks a cell or the grid is empty.
int32_t terminalsink(const std::vector<float>& elevation, const std::vector<uint8_t>& sea);
}
