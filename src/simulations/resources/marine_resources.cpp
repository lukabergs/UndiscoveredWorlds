#include "marine_resources.hpp"
#include "physical_layers_internal.hpp"

namespace physical_layers
{
using namespace detail;

namespace
{
bool isshelfproxy(const planet& world, const std::vector<std::vector<bool>>& shelves, int x, int y)
{
    if (world.sea(x, y) == 0)
        return false;

    if (x < static_cast<int>(shelves.size()) && y < static_cast<int>(shelves[x].size()) && shelves[x][y])
        return true;

    if (world.nom(x, y) < world.sealevel() - 220)
        return false;

    if (world.coast(x, y))
        return true;

    return hasfeatureinradius(world, x, y, 2, [&](int nx, int ny) { return world.outline(nx, ny); });
}

float averagesst(const planet& world, int x, int y)
{
    float total = 0.0f;

    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
        total += static_cast<float>(world.seasonalsst(season, x, y));

    return total / static_cast<float>(CLIMATESEASONCOUNT);
}

float averagecurrentspeed(const planet& world, int x, int y)
{
    float total = 0.0f;

    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
    {
        const float u = static_cast<float>(world.seasonalcurrentu(season, x, y));
        const float v = static_cast<float>(world.seasonalcurrentv(season, x, y));
        total += std::sqrt(u * u + v * v);
    }

    return total / static_cast<float>(CLIMATESEASONCOUNT);
}

float averagemaritimeinfluence(const planet& world, int x, int y)
{
    float total = 0.0f;

    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
        total += static_cast<float>(world.seasonalmaritimeinfluence(season, x, y));

    return total / static_cast<float>(CLIMATESEASONCOUNT);
}

float sstcontrast(const planet& world, int x, int y, float centresst)
{
    const int width = world.width();
    const int height = world.height();
    float totalcontrast = 0.0f;
    int samples = 0;

    for (const auto& offset : drainageoffsets)
    {
        const int ny = y + offset[1];

        if (ny < 0 || ny > height)
            continue;

        const int nx = wrap(x + offset[0], width);

        if (world.sea(nx, ny) == 0)
            continue;

        totalcontrast += std::abs(centresst - averagesst(world, nx, ny));
        samples++;
    }

    if (samples == 0)
        return 0.0f;

    return totalcontrast / static_cast<float>(samples);
}

float mean_ocean_sst(const planet& world)
{
    const int width = world.width();
    const int height = world.height();

    double ssttotal = 0.0;
    int sstcells = 0;

    for (int y = 0; y <= height; y++)
    {
        for (int x = 0; x <= width; x++)
        {
            if (world.sea(x, y) == 0)
                continue;

            ssttotal += averagesst(world, x, y);
            sstcells++;
        }
    }

    return sstcells > 0 ? static_cast<float>(ssttotal / static_cast<double>(sstcells)) : 0.0f;
}

void generate_ocean_fisheries(planet& world, const std::vector<std::vector<bool>>& shelves)
{
    const int width = world.width();
    const int height = world.height();
    const float meanoceansst = mean_ocean_sst(world);

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int y = startrow; y <= endrow; y++)
        {
            for (int x = 0; x <= width; x++)
            {
                if (world.sea(x, y) == 0)
                {
                    world.setfisheryreserve(x, y, 0);
                    continue;
                }

                const float shelfbias = isshelfproxy(world, shelves, x, y) ? 1.0f : 0.0f;
                const float cellsst = averagesst(world, x, y);
                const float coldwater = clampunit((meanoceansst - cellsst + 6.0f) / 18.0f);
                const float currentstrength = clampunit(averagecurrentspeed(world, x, y) / 160.0f);
                const float contrast = clampunit(sstcontrast(world, x, y, cellsst) / 18.0f);
                const float maritime = clampunit(averagemaritimeinfluence(world, x, y) / 220.0f);
                const float coastal = world.coast(x, y) ? 1.0f : (hasfeatureinradius(world, x, y, 2, [&](int nx, int ny) { return world.outline(nx, ny); }) ? 0.6f : 0.0f);
                const float icepenalty = world.seaice(x, y) == 2 ? 0.4f : (world.seaice(x, y) == 1 ? 0.75f : 1.0f);

                float score = 100.0f * (0.38f * shelfbias + 0.18f * coldwater + 0.16f * contrast + 0.14f * currentstrength + 0.08f * maritime + 0.06f * coastal);
                score *= icepenalty;

                world.setfisheryreserve(x, y, clampscore(score));
            }
        }
    }, 16);
}

void generate_coastal_fisheries(planet& world)
{
    const int width = world.width();
    const int height = world.height();

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int y = startrow; y <= endrow; y++)
        {
            for (int x = 0; x <= width; x++)
            {
                if (world.outline(x, y) == 0)
                    continue;

                int nearbyfishery = 0;

                for (const auto& offset : drainageoffsets)
                {
                    const int ny = y + offset[1];

                    if (ny < 0 || ny > height)
                        continue;

                    const int nx = wrap(x + offset[0], width);

                    if (world.sea(nx, ny) != 0)
                        nearbyfishery = std::max(nearbyfishery, world.fisheryreserve(nx, ny));
                }

                world.setfisheryreserve(x, y, clampscore(static_cast<float>(nearbyfishery) * 0.7f));
            }
        }
    }, 16);
}
}

void generate_marine_resources(planet& world, const std::vector<std::vector<bool>>& shelves)
{
    // Coastal land reads the completed neighbouring ocean fishery scores.
    generate_ocean_fisheries(world, shelves);
    generate_coastal_fisheries(world);
}
}
