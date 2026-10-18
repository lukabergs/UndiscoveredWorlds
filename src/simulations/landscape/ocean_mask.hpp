#pragma once

// The sea of a below-sea mask (docs/architecture/erosion-hydrology-plan.md, Phase H,
// step H3): every 8-connected region of cells at or below the sea level whose area
// is at least a given fraction of the planet's, and always the largest region by
// area (the world ocean). Neighbours are those of the metric's grid (on lat-lon
// longitude wraps, no neighbour beyond a pole; a torus wraps rows too); areas are the metric's cell
// areas, so a polar region of many small cells does not outweigh a larger one at
// low latitude. A smaller enclosed below-sea region is not sea: the routing fills it
// from its rim and the lake balance decides whether it holds a lake, a salt lake or
// a salt pan. Planet-free: the caller passes the fraction
// (landscape.routing.ocean_minimum_area_fraction).

#include <cstdint>
#include <vector>

#include "simulations/landscape/sphere_metric.hpp"

namespace landscape
{
// `belowsea` marks the cells at or below the sea level (non-zero). Returns 1 on the
// cells of every 8-connected below-sea region whose area is at least
// `minimumareafraction` x metric.totalarea_m2(), and of the largest region by area
// (of two of equal area, the one holding the lower cell index), 0 elsewhere. A
// fraction of 0 makes every below-sea cell sea; a mask whose below-sea cells form
// one region comes back unchanged (as 0 and 1); a mask without a below-sea cell
// comes back all zero. Throws std::invalid_argument when the mask's size differs
// from metric.size().
std::vector<uint8_t> oceanmask(const SphereMetric& metric, const std::vector<uint8_t>& belowsea,
                               double minimumareafraction);
}
