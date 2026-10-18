#pragma once

// Questions about coasts on the world grid: distance to land and sea, coastal cells, land in one quadrant.
// The CellId forms (namespace terrain::coastgrid) work on any cell grid (T8 W2-T2); the (x, y)
// forms are the lat-lon callers' (the legacy regional map) and forward to them.

#include <cstdint>
#include <string>
#include <vector>

#include "core/classes.hpp"
#include "core/mesh/cell_grid.hpp"

class planet;

int landdistance(planet& world, int x, int y);
twointegers nearestsea(planet& world, int i, int j, bool land, int limit, int grain);
bool vaguelycoastal(planet& world, int x, int y);
bool northwestlandonly(planet& world, int x, int y);
bool northeastlandonly(planet& world, int x, int y);
bool southwestlandonly(planet& world, int x, int y);
bool southeastlandonly(planet& world, int x, int y);

namespace terrain::coastgrid
{
// The steps from c toward the equator (toward the poles where a pole is warmer than the
// equator) to the first land with land two steps further and on both sides of the next step;
// -1 when none before the target.
int landdistance(const planet& world, mesh::CellId c);
// The nearest cell (great circle) of sea (land when land) within limit_km of c, -1 when none.
std::int32_t nearestsea(const planet& world, mesh::CellId c, bool land, double limit_km);
// A sea cell next to land or a land cell next to sea (the eight neighbours).
bool vaguelycoastal(const planet& world, mesh::CellId c);
// Land, bordered by sea to the east, south and southeast, with land to the northwest; the other
// three are the same turned. -1 (no cell) is false.
bool northwestlandonly(const planet& world, std::int32_t c);
bool northeastlandonly(const planet& world, std::int32_t c);
bool southwestlandonly(const planet& world, std::int32_t c);
bool southeastlandonly(const planet& world, std::int32_t c);
}

// The Voronoi panels of the continental shelves: per cell the number of its panel (1 upwards),
// panels of seeds about pointdist legacy pixels apart.
void makeshelvesvoronoi(planet& world, std::vector<int>& voronoi, int pointdist);
