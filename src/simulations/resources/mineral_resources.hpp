#ifndef mineral_resources_hpp
#define mineral_resources_hpp

class planet;

namespace physical_layers
{
// Compute mineral potential scores after basin, erosion and deposition layers.
void generate_mineral_resources(planet& world);
}

#endif
