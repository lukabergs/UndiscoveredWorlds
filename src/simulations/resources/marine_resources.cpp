#include "simulations/resources/marine_resources.hpp"
#include "simulations/resources/physical_layers_internal.hpp"

namespace physical_layers
{
using namespace detail;
using namespace cellqueries;

namespace
{
bool isshelfproxy(const planet& world, const worldgrid<bool>& shelves, CellId c)
{
    if (world.sea(c) == 0)
        return false;

    if (c < shelves.size() && shelves[c])
        return true;

    if (world.nomfloat(c) < world.sealevel() - 220)
        return false;

    if (iscoast(world, c))
        return true;

    return anywithinpixels(world, c, 2, [&](CellId n) { return isoutline(world, n); });
}

float averagesst(const planet& world, CellId c)
{
    float total = 0.0f;

    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
        total += world.seasonalsstfloat(season, c);

    return total / static_cast<float>(CLIMATESEASONCOUNT);
}

float averagecurrentspeed(const planet& world, CellId c)
{
    float total = 0.0f;

    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
    {
        const float u = static_cast<float>(world.seasonalcurrentu(season, c));
        const float v = static_cast<float>(world.seasonalcurrentv(season, c));
        total += std::sqrt(u * u + v * v);
    }

    return total / static_cast<float>(CLIMATESEASONCOUNT);
}

float averagemaritimeinfluence(const planet& world, CellId c)
{
    float total = 0.0f;

    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
        total += static_cast<float>(world.seasonalmaritimeinfluence(season, c));

    return total / static_cast<float>(CLIMATESEASONCOUNT);
}

// Mean SST difference to the sea cells among the eight neighbours (slots 1..8).
float sstcontrast(const planet& world, CellId c, float centresst)
{
    const mesh::CellGrid& grid = world.grid();
    float totalcontrast = 0.0f;
    int samples = 0;

    for (int slot = 1; slot <= mesh::slotcount; slot++)
    {
        const std::int32_t n = grid.neighbour(c, slot);

        if (n < 0 || world.sea(static_cast<CellId>(n)) == 0)
            continue;

        totalcontrast += std::abs(centresst - averagesst(world, static_cast<CellId>(n)));
        samples++;
    }

    if (samples == 0)
        return 0.0f;

    return totalcontrast / static_cast<float>(samples);
}

float mean_ocean_sst(const planet& world)
{
    double ssttotal = 0.0;
    int sstcells = 0;

    for (CellId c = 0; c < world.cellcount(); c++)
    {
        if (world.sea(c) == 0)
            continue;

        ssttotal += averagesst(world, c);
        sstcells++;
    }

    return sstcells > 0 ? static_cast<float>(ssttotal / static_cast<double>(sstcells)) : 0.0f;
}

void generate_ocean_fisheries(planet& world, const worldgrid<bool>& shelves)
{
    const float meanoceansst = mean_ocean_sst(world);

    parallelforcells(world, [&](CellId first, CellId last)
    {
        for (CellId c = first; c <= last; c++)
        {
            if (world.sea(c) == 0)
            {
                world.setfisheryreserve(c, 0);
                continue;
            }

            const float shelfbias = isshelfproxy(world, shelves, c) ? 1.0f : 0.0f;
            const float cellsst = averagesst(world, c);
            const float coldwater = clampunit((meanoceansst - cellsst + 6.0f) / 18.0f);
            const float currentstrength = clampunit(averagecurrentspeed(world, c) / 160.0f);
            const float contrast = clampunit(sstcontrast(world, c, cellsst) / 18.0f);
            const float maritime = clampunit(averagemaritimeinfluence(world, c) / 220.0f);
            const float coastal = iscoast(world, c) ? 1.0f : (anywithinpixels(world, c, 2, [&](CellId n) { return isoutline(world, n); }) ? 0.6f : 0.0f);
            const float icepenalty = world.seaice(c) == 2 ? 0.4f : (world.seaice(c) == 1 ? 0.75f : 1.0f);

            float score = 100.0f * (0.38f * shelfbias + 0.18f * coldwater + 0.16f * contrast + 0.14f * currentstrength + 0.08f * maritime + 0.06f * coastal);
            score *= icepenalty;

            world.setfisheryreserve(c, clampscore(score));
        }
    });
}

void generate_coastal_fisheries(planet& world)
{
    const mesh::CellGrid& grid = world.grid();

    parallelforcells(world, [&](CellId first, CellId last)
    {
        for (CellId c = first; c <= last; c++)
        {
            if (!isoutline(world, c))
                continue;

            int nearbyfishery = 0;

            for (int slot = 1; slot <= mesh::slotcount; slot++)
            {
                const std::int32_t n = grid.neighbour(c, slot);

                if (n >= 0 && world.sea(static_cast<CellId>(n)) != 0)
                    nearbyfishery = std::max(nearbyfishery, world.fisheryreserve(static_cast<CellId>(n)));
            }

            world.setfisheryreserve(c, clampscore(static_cast<float>(nearbyfishery) * 0.7f));
        }
    });
}
}

void generate_marine_resources(planet& world, const worldgrid<bool>& shelves)
{
    // Coastal land reads the completed neighbouring ocean fishery scores.
    generate_ocean_fisheries(world, shelves);
    generate_coastal_fisheries(world);
}
}
