#pragma once

// Lake queries for the wetlands, the regional map and the arid features: lake neighbours and
// river start points on a lake's edge.

#include <vector>

#include "core/mesh/cell_grid.hpp"

class planet;

float nexttolake(planet& world, mesh::CellId c);
float nexttolake(planet& world, int x, int y);  // lat-lon worlds (the legacy regional map)
// Lat-lon worlds only: thislake is the caller's (x, y) array, the box its columns and rows.
void makelakestartpoint(planet& world, std::vector<std::vector<int>>& thislake, int lakeno, int leftx, int lefty, int rightx, int righty);
