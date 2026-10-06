#pragma once

// Priority-Flood (Barnes et al. 2014) over a terrain on a SphereMetric, seeded
// from the sea cells and nothing else (polar rows are closed): the filled surface,
// one id per depression, its spill, and the flood receiver of every land cell
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
};

struct Depressions
{
    std::vector<float> filledelevation;   // >= terrain; equals the terrain outside depressions and on sea cells
    std::vector<int32_t> depressionid;    // 0 none; else the id of the connected flooded region
                                          // whose filled elevation exceeds the terrain by more than minimumdepth
    // The neighbour each land cell was reached from in the flood order (Barnes'
    // flow directions): the chain of flood receivers from any land cell leads
    // through its depression's spill to a sea cell. Sea cells: -1.
    std::vector<int32_t> floodreceiver;
    std::vector<Depression> depressions;  // by id - 1
};

// `sea` marks sea cells (1), the seeds of the flood, which enter at their own
// elevation (the seabed): the caller's mask must put every cell at or below the
// sea level in the sea (basebelowsea), so that no land cell lies below a sea
// neighbour's seabed; `elevation` is the terrain.
// A flooded cell is raised to the elevation of the cell it was reached from;
// cells raised by less than `minimumdepth_m` keep depressionid 0 (and their
// filledelevation is still the raised value). Deterministic: the priority queue
// breaks ties by cell index.
Depressions filldepressions(const SphereMetric& metric, const std::vector<float>& elevation,
                            const std::vector<uint8_t>& sea, float minimumdepth_m);
}
