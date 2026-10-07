#pragma once

// Lake queries for the wetlands, the regional map and the arid features: lake neighbours and
// river start points on a lake's edge.

#include <vector>

class planet;

float nexttolake(planet& world, int x, int y);
void makelakestartpoint(planet& world, std::vector<std::vector<int>>& thislake, int lakeno, int leftx, int lefty, int rightx, int righty);
