#ifndef physical_layers_hpp
#define physical_layers_hpp

#include <vector>

#include "core/world/world_grid.hpp"

class planet;

// Derive physical potential layers after terrain, drainage and climate finish.
// Scores preserve the existing 0..100 heuristic interpretation.
void generatephysicalworldlayers(planet& world, const worldgrid<bool>& shelves);

#endif
