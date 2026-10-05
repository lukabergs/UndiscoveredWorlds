#pragma once

// The same questions on a regional map (coasts, flow, lakes, rivers), and temperature adjusted for regional elevation.

#include <string>
#include <vector>

#include "core/classes.hpp"

class planet;
class region;

int tempelevadd(planet& world, region& region, int temp, int i, int j);
int tempelevremove(planet& world, region& region, int temp, int i, int j);
int getslope(region& region, int x, int y, int xx, int yy);
twointegers nearestsea(region& region, int leftx, int lefty, int rightx, int righty, int i, int j);
bool vaguelycoastal(region& region, int x, int y);
bool northwestlandonly(region& region, int x, int y);
bool lakenorthwestlandonly(region& region, int x, int y);
bool northeastlandonly(region& region, int x, int y);
bool lakenortheastlandonly(region& region, int x, int y);
bool southwestlandonly(region& region, int x, int y);
bool lakesouthwestlandonly(region& region, int x, int y);
bool southeastlandonly(region& region, int x, int y);
bool lakesoutheastlandonly(region& region, int x, int y);
twointegers getregionalflowdestination(region& region, int x, int y, int dir);
int checkregionalwaterinflow(region& region, int x, int y);
twointegers gettotalinflow(region& region, int x, int y);
twointegers findlowesthigher(region& region, int dx, int dy, int x, int y, int janload, int julload, int crount, std::vector<std::vector<bool>>& mainriver);
int getnearestlakelevel(region& region, int x, int y);
int getnearestlakespecial(region& region, int x, int y);
twointegers findclosestriver(region& region, int x, int y, bool delta);
twointegers findclosestriverquickly(region& region, int x, int y);
int countinflows(region& region, int x, int y);
