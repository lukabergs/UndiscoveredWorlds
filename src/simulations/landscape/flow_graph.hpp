#pragma once

// Drainage tree of a terrain on a SphereMetric: D8 receivers, the Braun-Willett
// order and the accumulated area and discharge
// (docs/architecture/erosion-hydrology-plan.md, section 3, Kernels).

#include <cstddef>
#include <cstdint>
#include <vector>

#include "simulations/landscape/sphere_metric.hpp"

namespace landscape
{
struct Depressions;

struct FlowGraph
{
    // Per cell. A land cell's receiver is the neighbour of steepest descent (drop
    // over centre distance) among the eight, a sea neighbour counting with
    // elevation = sea level. A land cell whose steepest descent is a sea neighbour
    // is an outlet: receiver = the cell itself, direction = 0, outletsea = that sea
    // neighbour. A land cell with no lower neighbour, or flooded by a depression
    // (depressionid != 0), takes the depression kernel's flood receiver (towards
    // its spill), so every land cell has a receiver and every chain ends at an
    // outlet. A sea cell's receiver is its lowest sea neighbour (the seabed's
    // steepest descent), itself when none is lower; sea cells have direction 0.
    std::vector<int32_t> receiver;        // index of the receiver; the cell itself for outlets and seabed pits
    std::vector<uint8_t> direction;       // 1..8 towards the receiver; 0 for outlets, sea cells and self-receivers
    std::vector<float> distance_m;        // centre to centre; 0 when the receiver is the cell itself
    std::vector<int32_t> outletsea;       // outlets: the sea cell the water enters; -1 elsewhere
    std::vector<uint8_t> land;            // 1 land, 0 sea (the input mask)
    // The Braun-Willett stack over the land cells: outlets first, every cell after
    // its receiver. Accumulation and the implicit solve run over it; accumulation
    // reversed (upstream cells first).
    std::vector<int32_t> order;
    // Accumulated over the tree: a cell includes every cell upstream of it. Sea cells: 0.
    std::vector<double> area_m2;          // drainage area
    std::vector<double> discharge_m3_yr;  // sum over the upstream cells of runoff x area

    std::size_t size() const { return receiver.size(); }
};

// Receivers, directions, distances, outlets and the order; area and discharge are
// zero until accumulateflow. `sea` marks sea cells (1); `elevation` is the terrain
// of every cell; a sea neighbour of a land cell counts with elevation `sealevel`;
// `depressions` is the Priority-Flood result on the same inputs (depressions.hpp)
// and supplies the receivers of flooded and pit cells. Ties in the steepest descent
// go to the lowest direction code, so the result is deterministic.
FlowGraph buildflowgraph(const SphereMetric& metric, const std::vector<float>& elevation,
                         const std::vector<uint8_t>& sea, float sealevel, const Depressions& depressions);

// Area and discharge accumulated along `order` reversed. `runoff_m_yr` is the
// runoff depth per cell per year (sea cells ignored): discharge = sum over the
// cell and its upstream cells of runoff x area. Serial, so bit-identical for any thread count.
void accumulateflow(const SphereMetric& metric, const std::vector<float>& runoff_m_yr, FlowGraph& graph);
}
