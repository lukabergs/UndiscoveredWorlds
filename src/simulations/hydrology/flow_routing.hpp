#pragma once

// hydrology.flow_routing (diagnoser): the planet-to-kernel adapter of the
// landscape routing (docs/architecture/erosion-hydrology-plan.md, section 3,
// Processes). The pieces are shared with terrain.fluvial_erosion, which routes
// at every step and writes the final routing itself.

#include <cstdint>
#include <vector>

#include "simulations/landscape/depressions.hpp"
#include "simulations/landscape/flow_graph.hpp"
#include "simulations/landscape/runoff.hpp"
#include "simulations/landscape/sphere_metric.hpp"

class planet;

namespace flowrouting
{
// The kernel's view of the planet: the metric from planetbody.radius, the base
// elevation per cell (row-major, north row first), the sea mask (oceancells: the
// below-sea regions that count as ocean), the sea level, the potential evaporation and the runoff in m/yr
// of every cell (sea included). The potential evaporation is the open-water law
// of landscape.lakes (landscape::openwaterevaporation) applied to each season's
// temperature and averaged over the four; the runoff, the transmission losses and
// the lake balance all use it. With landscape.runoff.budyko on, the runoff is
// Fu's curve of the annual precipitation and that potential evaporation; off, the
// coefficient law reduced by snow and ice (landscape.runoff.*).
struct Inputs
{
    landscape::SphereMetric metric;
    std::vector<float> elevation;
    std::vector<uint8_t> sea;
    std::vector<float> runoff_m_yr;
    std::vector<float> evaporation_m_yr;  // potential (open-water) evaporation
    float sealevel = 0.0f;
    // landscape.routing.*, seeded from the world seed (erosion plan C3).
    landscape::RoutingOptions routingoptions;
    // landscape.transmission.* (erosion plan C2): with transmissionenabled, route()
    // accumulates with the losses; the deficit is max(0, PET - P) per cell (m/yr),
    // open water is hydrosphere.imported_lake. Empty vectors when off. route()
    // adds the depression ids of its own Priority-Flood.
    bool transmissionenabled = false;
    landscape::TransmissionLoss transmission;
};

// Reads the world. The elevation is the base elevation (nomfloat): procedural
// worlds carry no template mountains since the FastLEM removal (plan, section 6).
Inputs gatherinputs(const planet& world);

// Whether an elevation is at or below the sea level in whole metres: the rule of
// planet::basebelowsea.
bool belowsea(float elevation, float sealevel);

// The sea mask of an elevation (erosion plan H3): the below-sea cells (belowsea) of
// the 8-connected below-sea regions that landscape::oceanmask keeps as sea (the
// largest, and every one of at least landscape.routing.ocean_minimum_area_fraction
// of the planet's area); a smaller enclosed below-sea region is land to the
// routing, a depression for the lake balance. oceancells gives the planet's mask by
// the same rule; gatherinputs takes Inputs::sea from it.
// The evolver recomputes it at every routing, so a seabed that marine deposition
// raises above the sea level becomes land and a coastal cell eroded below it
// becomes sea. Without a below-sea cell the mask is all land: Priority-Flood starts
// from the terminal sink, the lowest cell (landscape::filldepressions), every other
// depression drains towards it and its basin is the terminal depression, whose lake
// hydrology.lake_levels keeps closed.
std::vector<uint8_t> seamask(const landscape::SphereMetric& metric, const std::vector<float>& elevation, float sealevel);

// The planet's sea for hydrology (erosion plan H3): seamask of the base elevation on
// the metric of gatherinputs (row-major, north row first, (width + 1) x (height + 1)
// cells) in an eroded world (geosphere.erosion_duration > 0). A world without erosion
// (imported terrain) takes its sea from the import: every cell at or below the sea
// level that hydrosphere.imported_lake does not mark, whatever its region's size. The one rule of the routing, the lake levels, the rivers and the eroded
// drainage basins; planet::sea (the surface without a lake surface) stays the rule
// of the legacy stages and the climate.
std::vector<uint8_t> oceancells(const planet& world);

// Priority-Flood and the drainage tree with accumulated area and discharge on
// the inputs' elevation (or on `elevation` when given: the evolver's working
// copy, whose sea mask is recomputed); the discharge loses the transmission
// losses when the inputs enable them.
struct Routing
{
    landscape::Depressions depressions;
    landscape::FlowGraph graph;
    std::vector<uint8_t> sea;
    double transmissionloss_m3_yr = 0.0;  // total over the land; 0 without losses
};
Routing route(const Inputs& inputs, const std::vector<float>* elevation = nullptr);

// Writes hydrosphere.flow_direction, drainage_area (km2), discharge (m3/s,
// annual mean), filled_elevation, depression_id and runoff (mm/yr, every cell).
void writeroutingfields(planet& world, const Inputs& inputs, const Routing& routing);
}

// The stage: gather, route, write.
void applyflowroutingstage(planet& world);
