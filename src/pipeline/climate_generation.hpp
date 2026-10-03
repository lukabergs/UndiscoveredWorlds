#pragma once

#include <vector>

#include "world/world_grid.hpp"

class planet;
class boolshapetemplate;
struct ImportedClimateMaps;

// Orders world climate, drainage and biome stages; imported fields retain their seasonal shape.
void generateglobalclimate(planet& world, bool dorivers, bool dolakes, bool dodeltas, boolshapetemplate smalllake[], boolshapetemplate largelake[], boolshapetemplate landshape[], std::vector<std::vector<int>>& mountaindrainage, worldgrid<bool>& shelves, const ImportedClimateMaps* importedClimate = nullptr);
