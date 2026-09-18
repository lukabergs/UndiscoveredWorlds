#ifndef surface_potentials_hpp
#define surface_potentials_hpp

class planet;

namespace physical_layers
{
// Compute erosion, deposition and fertility scores after drainage-basin classification.
void generate_surface_potentials(planet& world);
}

#endif
