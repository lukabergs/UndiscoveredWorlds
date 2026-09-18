#pragma once

#include <vector>

class planet;
class boolshapetemplate;
struct ImportedClimateMaps;

// Orders world climate, drainage and biome stages; imported fields retain their seasonal shape.
void generateglobalclimate(planet& world, bool dorivers, bool dolakes, bool dodeltas, boolshapetemplate smalllake[], boolshapetemplate largelake[], boolshapetemplate landshape[], std::vector<std::vector<int>>& mountaindrainage, std::vector<std::vector<bool>>& shelves, const ImportedClimateMaps* importedClimate = nullptr);
