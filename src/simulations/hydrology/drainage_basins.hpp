#ifndef drainage_basins_hpp
#define drainage_basins_hpp

class planet;

namespace physical_layers
{
// Classify sea-connected, coastal and closed basins from the finalized drainage graph.
void generate_drainage_basins(planet& world);
}

#endif
