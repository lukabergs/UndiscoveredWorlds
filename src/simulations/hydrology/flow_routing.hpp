#pragma once

// hydrology.flow_routing (diagnoser): the planet-to-kernel adapter of the
// landscape routing (docs/architecture/erosion-hydrology-plan.md, section 3,
// Processes). The pieces are shared with terrain.fluvial_erosion, which routes
// at every step and writes the final routing itself.

#include <cstdint>
#include <vector>

#include "simulations/landscape/depressions.hpp"
#include "simulations/landscape/flow_graph.hpp"
#include "simulations/landscape/sphere_metric.hpp"

class planet;

namespace flowrouting
{
// The kernel's view of the planet: the metric from planetbody.radius, the base
// elevation per cell (row-major, north row first), the sea mask from
// basebelowsea, the sea level, and the runoff in m/yr from the seasonal
// precipitation reduced by snow and ice (landscape.runoff.*).
struct Inputs
{
    landscape::SphereMetric metric;
    std::vector<float> elevation;
    std::vector<uint8_t> sea;
    std::vector<float> runoff_m_yr;
    float sealevel = 0.0f;
};

// Reads the world. The elevation is the base elevation (nomfloat): procedural
// worlds carry no template mountains since the FastLEM removal (plan, section 6).
Inputs gatherinputs(const planet& world);

// The sea mask of an elevation: at or below the sea level in whole metres, the
// rule of planet::basebelowsea. The evolver recomputes it at every routing, so a
// seabed that marine deposition raises above the sea level becomes land and a
// coastal cell eroded below it becomes sea.
std::vector<uint8_t> seamask(const std::vector<float>& elevation, float sealevel);

// Priority-Flood and the drainage tree with accumulated area and discharge on
// the inputs' elevation (or on `elevation` when given: the evolver's working
// copy, whose sea mask is recomputed).
struct Routing
{
    landscape::Depressions depressions;
    landscape::FlowGraph graph;
    std::vector<uint8_t> sea;
};
Routing route(const Inputs& inputs, const std::vector<float>* elevation = nullptr);

// Writes hydrosphere.flow_direction, drainage_area (km2), discharge (m3/s,
// annual mean), filled_elevation and depression_id.
void writeroutingfields(planet& world, const Inputs& inputs, const Routing& routing);
}

// The stage: gather, route, write.
void applyflowroutingstage(planet& world);
