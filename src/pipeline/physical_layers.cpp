#include "physical_layers.hpp"
#include "../simulations/hydrology/drainage_basins.hpp"
#include "../simulations/terrain/surface_potentials.hpp"
#include "../simulations/resources/mineral_resources.hpp"
#include "../simulations/resources/marine_resources.hpp"

void generatephysicalworldlayers(planet& world, const std::vector<std::vector<bool>>& shelves)
{
    physical_layers::generate_drainage_basins(world);
    physical_layers::generate_surface_potentials(world);
    physical_layers::generate_mineral_resources(world);
    physical_layers::generate_marine_resources(world, shelves);
}
