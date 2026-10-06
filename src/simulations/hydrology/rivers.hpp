#pragma once

// hydrology.rivers (diagnoser, erosion plan B1): river_direction, river_flow_january,
// river_flow_july and max_river_flow from the routing's flow_direction and discharge
// in recipe eroded-world, where it replaces hydrology.rivers_and_basins and
// hydrology.basin_policy. Flows are in the legacy units of river_factor (a cell
// of mean size with 15 mm/month of precipitation contributes 1 unit before the
// runoff coefficient), split into January and July by the seasonal runoff.

class planet;

void applyriverstage(planet& world);
