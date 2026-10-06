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
// Open-water evaporation in m/yr at a temperature (C): perc x max(0, T - zeroc), the
// law of landscape.lakes.evaporation_*. Callers apply it per season and average.
double openwaterevaporation(double temperaturec, double perc, double zeroc);

// Fu (1981): AET / P = 1 + PET / P - (1 + (PET / P)^w)^(1 / w); runoff = P - AET
// = (P^w + PET^w)^(1 / w) - PET, in the units of the arguments (m/yr), within
// [0, P]. 0 where P <= 0; P where PET <= 0. w > 1 (std::invalid_argument otherwise).
double budykorunoff(double precipitation, double potentialevaporation, double w);

struct TransmissionLoss
{
    // Per cell, m/yr: the evaporation deficit of open water, max(0, PET - P); 0 turns
    // the loss of the cell off.
    std::vector<float> deficit_m_yr;
    // Per cell, 1 where the cell is open water a river crosses (an imported lake):
    // the deficit acts on the whole cell area instead of the channel. May be empty.
    std::vector<uint8_t> openwater;
    // Per cell, the routing's depression id (depressions.hpp): a cell with an id other
    // than 0 is flooded by a depression and loses nothing here, because the lake
    // balance owns its evaporation. May be empty (no depressions).
    std::vector<int32_t> depressionid;
    double widthcoefficient = 0.0;  // channel width (m) = coefficient x (discharge in m3/s)^exponent
    double widthexponent = 0.5;
    double floodplainfactor = 1.0;  // width multiplier where the link's slope is below floodplainslope
    double floodplainslope = 0.0;
    // The elevation an outlet's sea cell counts with in its link slope (the sea
    // level the graph was built with, as buildflowgraph).
    float sealevel = 0.0f;
};

// As accumulateflow (flow_graph.hpp), but each land cell outside the depressions
// (depressionid 0), before passing its flow to
// its receiver (an outlet: to the sea), loses min(flow, width x link length x
// deficit) (open water: min(flow, deficit x cell area)). The flow is the cell's
// accumulated discharge, its own runoff plus what its upstream cells passed on;
// the width is widthcoefficient x (flow in m3/s)^widthexponent, times
// floodplainfactor where the link's slope (drop to the receiver, or to the sea
// level for an outlet, over FlowGraph::distance_m) is below floodplainslope; the
// link length is FlowGraph::distance_m (0 for a cell that is its own receiver
// without a sea cell, which then loses nothing unless it is open water).
// discharge_m3_yr of a cell is the flow it passes on, after its loss. `elevation`
// is the terrain the graph was built on. Discharge then no longer grows
// monotonically downstream; drainage area still does. A deficit of 0 everywhere
// reproduces accumulateflow bit for bit. Returns the total loss in m3/yr. Serial
// over the order, so bit-identical for any thread count.
double accumulateflowwithlosses(const SphereMetric& metric, const std::vector<float>& runoff_m_yr,
                                const std::vector<float>& elevation, const TransmissionLoss& loss, FlowGraph& graph);
}
