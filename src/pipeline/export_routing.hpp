#pragma once

// The routing of a lat-lon export and the exports of a world for its consumers
// (docs/architecture/sphere-mesh.md sections 6, 7 and 9). Direction fields are
// never remapped and accumulated fluxes are point values of the drainage network,
// so the export recomputes both on its own grid with the landscape kernels, from
// the exported filled elevation and runoff:
// - the routing (flow_direction, drainage_area, discharge, filled_elevation): the
//   Priority-Flood of the exported filled surface (pits the remap left are filled),
//   D8 with the routing's flat resolution, and the area and the lossless discharge
//   of the exported runoff accumulated over the tree; only where the world was
//   routed (a drainage area somewhere);
// - the legacy rivers (river_direction, river_flow_*), only where the world has
//   river flows: hydrology.rivers' flow law on the recomputed discharge, January and
//   July split by the native flows' ratio accumulated over the tree (weighted by
//   runoff x area), the river mouths on the sea cells the outlets drain into;
// - the distributaries (delta_direction, delta_flow_*) of a world with an erosion
//   run, by landscape::splitdistributaries on the recomputed fields.
// The legacy units, the minimum flow and the flat resolution are those of
// hydrology.rivers and hydrology.flow_routing (landscape.json).

#include <memory>
#include <string>
#include <vector>

#include "core/world/world_fields.hpp"
#include "world/planet_bundle.hpp"

class planet;
class ClimateShadow;

// The fields recomputeexportrouting writes (the export entries it replaces).
const std::vector<worldmodel::FieldId>& exportroutingfields();

// Recomputes the routing, the rivers and the distributaries of a lat-lon planet
// (an export) in place, as described above. The LatlonDirectionRecompute of
// exportlatlon and writeplanetbundle.
void recomputeexportrouting(planet& latlon);

// The lat-lon export of a world at its export width (4 N on a cube) with the
// recomputed routing; on a cube world with a climate shadow, the climate-owned
// fields the planet stores (the orographic rain directions, the subtropical-high
// profile) are the shadow's. On a lat-lon world, a copy.
std::unique_ptr<planet> exportworld(const planet& world, const ClimateShadow* shadow = nullptr);

// What a world's bundle takes for its lat-lon export: the recomputed routing and,
// with a shadow that exists, the climate-owned fields of climate's lat-lon world.
PlanetBundleExport bundleexport(const ClimateShadow* shadow = nullptr);

// Saves a world as .uww (sphere-mesh.md section 9): a cube world saves its export
// (exportworld), so the layout is unchanged and loading gives a lat-lon world.
void saveworldfile(planet& world, const std::string& filename, const ClimateShadow* shadow = nullptr);
