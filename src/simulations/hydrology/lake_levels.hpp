#pragma once

// hydrology.lake_levels (diagnoser, erosion plan B2): lake levels in the routing's
// depressions (hydrosphere.depression_id, filled_elevation, flow_direction,
// discharge) from the water balance of landscape::balancelakes, written to
// hydrosphere.lake_surface (m above datum; 0 no lake) with special_features 100
// (salt lake) on the cells of a closed lake and 110 (salt pan) on the lowest cell
// of a closed depression too dry to hold a lake. Replaces hydrology.lakes in
// eroded-world. A lake is closed when its level is below its depression's
// filled_elevation (the spill): basin_class (B3) reads exactly that. On imported
// terrain (earth-hydrology, erosion plan C1) the cells of hydrosphere.imported_lake
// not covered by a balance lake are water at their lake's level, the lowest base
// elevation of each 8-connected imported lake.

class planet;

void applylakelevelstage(planet& world);
