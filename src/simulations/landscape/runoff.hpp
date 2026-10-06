#pragma once

// Runoff and its losses along the rivers (docs/architecture/erosion-hydrology-plan.md,
// Phase C, step C2): Fu's Budyko curve for the runoff depth of a cell, and the
// accumulation of discharge along a FlowGraph with transmission losses, the
// evaporation from a river's channel and floodplain where it crosses land drier
// than open water (the Sudd's role for the Nile). Planet-free; the process that
// calls it (hydrology/flow_routing.cpp) passes the values of landscape.runoff.*
// and landscape.transmission.*.

#include <cstdint>
#include <vector>

#include "simulations/landscape/flow_graph.hpp"
#include "simulations/landscape/sphere_metric.hpp"

namespace landscape
{
// Fu (1981): AET / P = 1 + PET / P - (1 + (PET / P)^w)^(1 / w); runoff = P - AET,
// in the units of the arguments (m/yr). 0 where P <= 0; P where PET <= 0. w > 1.
double budykorunoff(double precipitation, double potentialevaporation, double w);

struct TransmissionLoss
{
    // Per cell, m/yr: the evaporation deficit of open water, max(0, PET - P); 0 turns
    // the loss of the cell off.
    std::vector<float> deficit_m_yr;
    // Per cell, 1 where the cell is open water a river crosses (an imported lake):
    // the deficit acts on the whole cell area instead of the channel.
    std::vector<uint8_t> openwater;
    double widthcoefficient = 0.0;  // channel width (m) = coefficient x (discharge in m3/s)^exponent
    double widthexponent = 0.5;
    double floodplainfactor = 1.0;  // width multiplier where the link's slope is below floodplainslope
    double floodplainslope = 0.0;
};

// As accumulateflow (flow_graph.hpp), but each land cell, before passing its flow to
// its receiver, loses min(flow, width x link length x deficit) (open water: deficit
// x cell area), the width from the flow entering the cell. `elevation` is the
// terrain the graph was built on (for the link slopes). Discharge then no longer
// grows monotonically downstream; drainage area still does. Returns the total loss
// in m3/yr. Serial over the order, so bit-identical for any thread count.
double accumulateflowwithlosses(const SphereMetric& metric, const std::vector<float>& runoff_m_yr,
                                const std::vector<float>& elevation, const TransmissionLoss& loss, FlowGraph& graph);
}
