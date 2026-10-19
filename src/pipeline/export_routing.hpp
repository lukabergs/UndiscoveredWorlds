#pragma once

// The routing of a lat-lon export and the exports of a world for its consumers
// (docs/architecture/sphere-mesh.md sections 6, 7 and 9). Direction fields are
// never remapped and accumulated fluxes are point values of the drainage network,
// so the export recomputes both on its own grid by hydrology's own law: the routing
// of hydrology.flow_routing (flowrouting::gatherinputs, route and writeroutingfields:
// flow_direction, drainage_area, discharge, filled_elevation, depression_id, runoff)
// on the exported filled surface, so pits the remap left are filled, where the world
// was routed (a drainage area somewhere); then hydrology.rivers (applyriverstage:
// river_direction, river_flow_*, the mouths and, after an erosion run, the
// distributaries delta_*) where the world has river flows. Deterministic from the
// world seed (the routing's draws).

#include <memory>
#include <string>
#include <vector>

#include "core/world/world_fields.hpp"
#include "world/planet_bundle.hpp"

class planet;
class ClimateShadow;

// The fields recomputeexportrouting writes (the export entries it replaces).
const std::vector<worldmodel::FieldId>& exportroutingfields();

// Recomputes the routing and the rivers of a lat-lon planet (an export) in place, as
// described above. The LatlonDirectionRecompute of exportlatlon and writeplanetbundle.
void recomputeexportrouting(planet& latlon);

// The lat-lon export of a world at its export width (4 N on a cube) with the
// recomputed routing; on a cube world with a climate shadow, the climate-owned
// fields the planet stores (the orographic rain directions, the subtropical-high
// profile) are the shadow's. On a lat-lon world, a copy.
std::unique_ptr<planet> exportworld(const planet& world, const ClimateShadow* shadow = nullptr);

// The export of a world loaded from a cube bundle (the desktop): exportworld, then the
// bundle's exports of the climate-owned fields (climateshadowowned) and its stored
// zonal profiles, climate's lat-lon ones. On a lat-lon world, a copy.
std::unique_ptr<planet> exportbundleworld(const planet& world, const PlanetBundleSource& source);

// What a world's bundle takes for its lat-lon export: the recomputed routing and,
// with a shadow that exists, the climate-owned fields of climate's lat-lon world.
PlanetBundleExport bundleexport(const ClimateShadow* shadow = nullptr);

// Saves a world as .uww (sphere-mesh.md section 9): a cube world saves its export
// (exportworld), so the layout is unchanged and loading gives a lat-lon world.
void saveworldfile(planet& world, const std::string& filename, const ClimateShadow* shadow = nullptr);
