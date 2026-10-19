#ifndef drainage_basins_hpp
#define drainage_basins_hpp

#include <cstdint>
#include <vector>

class planet;

namespace physical_layers
{
// Per cell: 1 when the cell belongs to a closed depression (erosion plan B3): its depression's
// highest lake surface lies more than 0.01 m below its highest filled elevation.
std::vector<std::uint8_t> closeddepressioncells(const planet& world);

// Classify sea-connected, coastal and closed basins from the finalized drainage graph.
void generate_drainage_basins(planet& world);
}

#endif
