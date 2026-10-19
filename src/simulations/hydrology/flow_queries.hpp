#pragma once

// Questions about drainage on the world grid: flatness, flow destinations, inflows, lakes nearby.
// The cell forms run on any cell grid (neighbours from the neighbour table, a direction code a
// slot of its cell); the (x, y) forms are the lat-lon world's (the legacy regional map) and
// forward to them.

#include <cstdint>

#include "core/classes.hpp"
#include "core/mesh/cell_grid.hpp"

class planet;

int getflatness(planet& world, mesh::CellId c);
float getflatelevation(planet& world, mesh::CellId c);
// The cell a direction (0: the cell's river direction) leads to; the cell itself where the slot
// has no neighbour.
mesh::CellId getflowdestination(planet& world, mesh::CellId c, int dir);
// The cell of the 3 x 3 neighbourhood (the cell included) with the largest flow into the cell;
// -1 when none.
std::int32_t getupstreamcell(planet& world, mesh::CellId c);
int checkwaterinflow(planet& world, mesh::CellId c);
twointegers gettotalinflow(planet& world, mesh::CellId c);
// The first lake surface (or, with rift, rift-lake surface) within dist steps of the neighbour
// table (the (2 dist + 1)^2 box of the lat-lon grid); 0 when none.
float nearlake(planet& world, mesh::CellId c, int dist, bool rift);
int getlakeedge(planet& world, mesh::CellId c);

int getflatness(planet& world, int x, int y);
float getflatelevation(planet& world, int x, int y);
twointegers getflowdestination(planet& world, int x, int y, int dir);
twointegers getupstreamcell(planet& world, int x, int y);
int checkwaterinflow(planet& world, int x, int y);
twointegers gettotalinflow(planet& world, int x, int y);
float nearlake(planet& world, int x, int y, int dist, bool rift);
int getlakeedge(planet& world, int x, int y);
