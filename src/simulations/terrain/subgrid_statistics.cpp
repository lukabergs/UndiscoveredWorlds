#include "simulations/terrain/subgrid_statistics.hpp"

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <vector>

#include "simulations/landscape/landscape_parameters.hpp"
#include "simulations/landscape/sphere_metric.hpp"
#include "simulations/landscape/subgrid_statistics.hpp"
#include "simulations/terrain/terrain_parameters.hpp"
#include "world/planet.hpp"

void applysubgridstatisticsstage(planet& world)
{
    namespace t = tuning::landscape;
    const int columns = world.width() + 1;
    const int rows = world.height() + 1;
    const landscape::SphereMetric metric(columns, rows, static_cast<double>(world.radius()));
    const std::size_t size = metric.size();

    landscape::SubgridParameters parameters;
    namespace s = tuning::terrain::subgrid;
    parameters.diffusivity_km2_myr = t::hillslope::diffusivityKm2Myr;
    parameters.spacingcoefficient = s::spacingCoefficient;
    parameters.referencedoverk = s::referenceDOverK;
    parameters.referencem = s::referenceM;
    parameters.referencerunoff_mm_yr = s::referenceRunoffMmYr;
    parameters.referenceerodibility = s::referenceErodibility;

    std::vector<float> elevation(size), erodibility(size), runoff(size), erosionrate(size);
    std::vector<std::uint8_t> sea(size);
    // Through a const reference, so that the reads are the tracked getters.
    const planet& reader = world;
    const worldgrid<float>& erodibilityfield = reader.erodibility();
    const worldgrid<float>& runofffield = reader.runoff();
    const worldgrid<float>& erosionratefield = reader.erosionrate();
    for (int y = 0; y < rows; y++)
    {
        for (int x = 0; x < columns; x++)
        {
            const std::size_t cell = metric.index(x, y);
            elevation[cell] = world.nomfloat(x, y);
            sea[cell] = world.basebelowsea(x, y) ? 1 : 0;
            erodibility[cell] = erodibilityfield[x][y];
            runoff[cell] = runofffield[x][y];
            erosionrate[cell] = erosionratefield[x][y];
        }
    }

    std::vector<float> relief, spacing;
    const std::size_t capped =
        landscape::subgridstatistics(metric, parameters, elevation, sea, erodibility, runoff, erosionrate, relief, spacing);
    std::size_t land = 0;
    for (std::uint8_t value : sea)
        land += value == 0 ? 1 : 0;
    std::cout << "Sub-grid statistics: valley spacing capped at the cell size on " << capped << " of " << land
              << " land cells\n";

    worldgrid<float>& relieffield = world.reliefamplitude();
    worldgrid<float>& spacingfield = world.valleyspacing();
    for (int y = 0; y < rows; y++)
    {
        for (int x = 0; x < columns; x++)
        {
            const std::size_t cell = metric.index(x, y);
            relieffield[x][y] = relief[cell];
            spacingfield[x][y] = spacing[cell];
        }
    }
}
