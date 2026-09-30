#pragma once

// Wetland placement and pruning from the completed drainage network.

class planet;
class boolshapetemplate;

void createwetlands(planet& world, boolshapetemplate smalllake[]);
void removeexcesswetlands(planet& world);
