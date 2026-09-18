#pragma once

// Lake evaporation rainfall and propagation of the resulting river flow.

#include <vector>

class planet;

void lakerain(planet& world, std::vector<std::vector<int>>& lakewinterrainmap, std::vector<std::vector<int>>& lakesummerrainmap);
void riftlakerain(planet& world, std::vector<std::vector<int>>& lakewinterrainmap, std::vector<std::vector<int>>& lakesummerrainmap);
void lakerivers(planet& world, std::vector<std::vector<int>>& lakewinterrainmap, std::vector<std::vector<int>>& lakesummerrainmap);
