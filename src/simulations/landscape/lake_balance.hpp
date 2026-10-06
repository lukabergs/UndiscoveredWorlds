#pragma once

// Lake levels by water balance in the depressions of a routing
// (docs/architecture/erosion-hydrology-plan.md, Phase B, step B2). Planet-free.
//
// For each depression d (the cells with depressionid == d, all at one filled level,
// the spill elevation):
//   inflow Q_d = the discharge of every cell outside d whose receiver lies in d
//              + the runoff of the cells of d (runoff x cell area), in m3/yr;
//   net loss per unit lake area L_d = max(0, E_d - P_d) in m/yr, E_d and P_d the
//              area-weighted means of evaporation and precipitation over d's cells;
//   equilibrium area A*_d = Q_d / L_d (infinite when L_d = 0).
// If A*_d reaches the area of d the lake overflows: every cell of d is lake at the
// spill level (Phase A's rule). Otherwise the lake covers the k lowest cells of d
// (terrain elevation ascending, ties by cell index), k the smallest count whose
// cumulative area reaches A*_d, and its level is midway between the terrain of the
// k-th lowest cell and that of the (k+1)-th (or the spill, if smaller); such a lake
// is closed. A closed depression whose A*_d is below the area of its lowest cell
// holds no lake (k = 0, level 0): it is dry.
// Simplifications of Phase B (plan section 6): the cells of d form one water body
// (two pits below a sill inside one depression share a level); the discharge is
// the routing's, which assumes every upstream depression overflows.

#include <cstddef>
#include <cstdint>
#include <vector>

#include "simulations/landscape/sphere_metric.hpp"

namespace landscape
{
struct LakeBalanceInputs
{
    const SphereMetric* metric = nullptr;
    std::vector<float> elevation;           // terrain, m above datum
    std::vector<float> filledelevation;     // the routing's filled surface (the spill level inside a depression)
    std::vector<int32_t> depressionid;      // 0 none; ids need not be contiguous
    std::vector<int32_t> receiver;          // per cell: the receiver index, -1 for none (sea, outlets into the sea)
    std::vector<double> discharge_m3_yr;    // the routing's accumulated discharge
    std::vector<float> runoff_m_yr;         // runoff depth per cell
    std::vector<float> precipitation_m_yr;  // annual precipitation per cell
    std::vector<float> evaporation_m_yr;    // open-water evaporation per cell
};

struct LakeLevel
{
    int32_t id = 0;                // the depression id
    float spill_m = 0.0f;          // the depression's filled level
    float level_m = 0.0f;          // the lake surface; 0 when dry
    bool overflowing = false;      // level == spill: water leaves through the spill
    int32_t lakecells = 0;         // k
    int32_t lowestcell = -1;       // the lowest cell of the depression (the pan of a dry one)
    double area_m2 = 0.0;          // the lake's area
    double depressionarea_m2 = 0.0;
    double inflow_m3_yr = 0.0;     // Q_d
    double netloss_m_yr = 0.0;     // L_d
};

struct LakeBalance
{
    std::vector<LakeLevel> lakes;   // one per depression id, ascending id
    std::vector<float> surface;     // per cell: the lake level on lake cells, 0 elsewhere
};

// Throws std::invalid_argument when a vector's size differs from metric->size().
LakeBalance balancelakes(const LakeBalanceInputs& inputs);
}
