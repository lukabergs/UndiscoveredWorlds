#pragma once

// Questions about coasts on the world grid: distance to land and sea, coastal cells, land in one quadrant.

#include <string>
#include <vector>

#include "core/classes.hpp"

class planet;
class region;

int landdistance(planet& world, int x, int y);
twointegers nearestsea(planet& world, int i, int j, bool land, int limit, int grain);
bool vaguelycoastal(planet& world, int x, int y);
bool northwestlandonly(planet& world, int x, int y);
bool northeastlandonly(planet& world, int x, int y);
bool southwestlandonly(planet& world, int x, int y);
bool southeastlandonly(planet& world, int x, int y);
void makeshelvesvoronoi(planet& world, std::vector<std::vector<short>>& voronoi, std::vector<std::vector<bool>>& outline, int pointdist);
