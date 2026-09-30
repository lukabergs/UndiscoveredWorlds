#ifndef physical_layers_hpp
#define physical_layers_hpp

#include <vector>

class planet;

// Derive physical potential layers after terrain, drainage and climate finish.
// Scores preserve the existing 0..100 heuristic interpretation.
void generatephysicalworldlayers(planet& world, const std::vector<std::vector<bool>>& shelves);

#endif
