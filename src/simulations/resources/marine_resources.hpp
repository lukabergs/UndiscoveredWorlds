#ifndef marine_resources_hpp
#define marine_resources_hpp

#include <vector>

#include "core/world/world_grid.hpp"

class planet;

namespace physical_layers
{
// Compute ocean and coastal fishery potential from finalized seasonal ocean fields.
void generate_marine_resources(planet& world, const worldgrid<bool>& shelves);
}

#endif
