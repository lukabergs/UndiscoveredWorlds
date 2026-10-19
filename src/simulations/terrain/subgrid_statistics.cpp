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
    // The world's own cell grid (a cube or the lat-lon grid), as flowrouting::worldmetric.
    const landscape::SphereMetric metric(world.sharedgrid());
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
    // The metric's cell index is the world's cell id (both in the grid's id order).
    for (std::size_t cell = 0; cell < size; cell++)
    {
        const mesh::CellId c = static_cast<mesh::CellId>(cell);
        elevation[cell] = world.nomfloat(c);
        sea[cell] = world.basebelowsea(c) ? 1 : 0;
        erodibility[cell] = erodibilityfield[c];
        runoff[cell] = runofffield[c];
        erosionrate[cell] = erosionratefield[c];
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
    for (std::size_t cell = 0; cell < size; cell++)
    {
        const mesh::CellId c = static_cast<mesh::CellId>(cell);
        relieffield[c] = relief[cell];
        spacingfield[c] = spacing[cell];
    }
}
