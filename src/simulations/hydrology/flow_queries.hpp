#pragma once

// Questions about drainage on the world grid: slopes, flow destinations, inflows, lakes nearby.

#include <string>
#include <vector>

#include "core/classes.hpp"

class planet;

int findlowestdir(planet& world, int neighbours[8][2], int x, int y);
int findlowestdirriver(planet& world, int neighbours[8][2], int x, int y, std::vector<std::vector<int>>& mountaindrainage);
float getslope(planet& world, int x, int y, int xx, int yy);
int getflatness(planet& world, int x, int y);
float getflatelevation(planet& world, int x, int y);
twointegers findseatile(planet& world, int x, int y, int dir);
twointegers getflowdestination(planet& world, int x, int y, int dir);
twointegers getupstreamcell(planet& world, int x, int y);
int checkwaterinflow(planet& world, int x, int y);
twointegers gettotalinflow(planet& world, int x, int y);
int nearlake(planet& world, int x, int y, int dist, bool rift);
int getlakeedge(planet& world, int x, int y);
