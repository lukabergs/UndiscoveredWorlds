#pragma once

// River direction planning, downstream catchment accumulation and shared flow edits.

#include <vector>

class planet;

void createrivermap(planet& world, std::vector<std::vector<int>>& mountaindrainage);
void removerivermountains(planet& world);
void checkglobalflows(planet& world);
void markriver(planet& world, int x, int y, std::vector<std::vector<int>>& markedarray, int extra);
int divertriver(planet& world, int x, int y, int destx, int desty, std::vector<std::vector<int>>& removedrivers, int riverno);
int divertlakeriver(planet& world, int x, int y, int destx, int desty, std::vector<std::vector<int>>& removedrivers, int riverno, int outflowx, int outflowy, std::vector<std::vector<int>>& avoidarray);
void removeriver(planet& world, std::vector<std::vector<int>>& removedrivers, int riverno, int x, int y);
void reduceriver(planet& world, int janreduce, int julreduce, std::vector<std::vector<int>>& removedrivers, int riverno, int x, int y);
void addtoriver(planet& world, int x, int y, int janload, int julload);
