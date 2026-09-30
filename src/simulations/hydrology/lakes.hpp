#pragma once

// Salt, freshwater and rift lake placement, including drainage connections and lake cleanup.

#include <vector>

class planet;
class boolshapetemplate;

void createsaltlakes(planet& world, int& lakesplaced, std::vector<std::vector<std::vector<int>>>& saltlakemap, std::vector<std::vector<int>>& nolake, std::vector<std::vector<int>>& basins, boolshapetemplate smalllake[]);
void placesaltlake(planet& world, int centrex, int centrey, bool large, bool dodepressions, std::vector<std::vector<std::vector<int>>>& saltlakemap, std::vector<std::vector<int>>& basins, std::vector<std::vector<int>>& avoid, std::vector<std::vector<int>>& nolake, boolshapetemplate smalllake[]);
void convertsaltlakes(planet& world, std::vector<std::vector<std::vector<int>>>& saltlakemap);
void createlakemap(planet& world, std::vector<std::vector<int>>& nolake, boolshapetemplate smalllake[], boolshapetemplate largelake[]);
void createriftlakemap(planet& world, std::vector<std::vector<int>>& nolake);
int nexttolake(planet& world, int x, int y);
void removesealakes(planet& world);
void connectlakes(planet& world);
void makelakestartpoint(planet& world, std::vector<std::vector<int>>& thislake, int lakeno, int leftx, int lefty, int rightx, int righty);
