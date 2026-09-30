#include "physical_layers.hpp"
#include "generation_progress.hpp"
#include "planet.hpp"
#include "../simulations/hydrology/drainage_basins.hpp"
#include "../simulations/terrain/surface_potentials.hpp"
#include "../simulations/resources/gold_deposits.hpp"
#include "../simulations/resources/mineral_resources.hpp"
#include "../simulations/resources/marine_resources.hpp"

#include <cstdio>

namespace
{
void reportgolddeposits(const planet& world)
{
    using namespace mineralsystems;
    int primary = 0, preserved = 0, placers = 0;
    double lodegold = 0.0, placergold = 0.0;
    for (const Deposit& deposit : world.mineraldeposits())
    {
        if (deposit.system == DepositSystem::placer_gold)
        {
            ++placers;
            placergold += deposit.containedGoldT();
            continue;
        }
        ++primary;
        if (deposit.fate == DepositFate::preserved)
        {
            ++preserved;
            lodegold += deposit.remainingGoldT();
        }
    }
    if (primary == 0)
        return;

    char text[256];
    std::snprintf(text, sizeof(text), "Gold deposits: %d primary (%d preserved, %.0f t), %d placers (%.1f t); %.1f t of liberated gold reached the sea",
        primary, preserved, lodegold, placers, placergold, world.mineralbudget().placerExportedGoldT);
    updatereport(text);
}
}

void generatephysicalworldlayers(planet& world, const std::vector<std::vector<bool>>& shelves)
{
    physical_layers::generate_drainage_basins(world);
    physical_layers::generate_surface_potentials(world);
    physical_layers::generate_mineral_resources(world);
    // Needs the final drainage and deposition potential.
    physical_layers::generate_gold_deposits(world);
    reportgolddeposits(world);
    physical_layers::generate_marine_resources(world, shelves);
}
