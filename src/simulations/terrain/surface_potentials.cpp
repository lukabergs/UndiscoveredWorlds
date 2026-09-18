#include "surface_potentials.hpp"
#include "../resources/physical_layers_internal.hpp"

namespace physical_layers
{
using namespace detail;

namespace
{
bool iswetlandcell(const planet& world, int x, int y)
{
    const int special = world.special(x, y);
    return special >= 130 && special < 140;
}

int averagetemperature(const planet& world, int x, int y)
{
    return (world.mintemp(x, y) + world.maxtemp(x, y)) / 2;
}

int localrelief(const planet& world, int x, int y)
{
    const int width = world.width();
    const int height = world.height();
    const int centre = world.map(x, y);
    int maxdifference = 0;

    for (const auto& offset : drainageoffsets)
    {
        const int ny = y + offset[1];

        if (ny < 0 || ny > height)
            continue;

        const int nx = wrap(x + offset[0], width);
        maxdifference = std::max(maxdifference, std::abs(world.map(nx, ny) - centre));
    }

    return maxdifference;
}

int inflowcount(const planet& world, int x, int y)
{
    const int width = world.width();
    const int height = world.height();
    int inflow = 0;

    for (const auto& offset : drainageoffsets)
    {
        const int ny = y + offset[1];

        if (ny < 0 || ny > height)
            continue;

        const int nx = wrap(x + offset[0], width);

        if (world.riverdir(nx, ny) != 0)
        {
            const twointegers dest = wrappeddestination(world, nx, ny, world.riverdir(nx, ny));

            if (dest.x == x && dest.y == y)
                inflow++;
        }

        if (world.deltadir(nx, ny) > 0)
        {
            const twointegers dest = wrappeddestination(world, nx, ny, world.deltadir(nx, ny));

            if (dest.x == x && dest.y == y)
                inflow++;
        }
    }

    return inflow;
}

float temperaturerangefactor(int temperature)
{
    if (temperature <= -10 || temperature >= 40)
        return 0.0f;

    if (temperature <= 18)
        return clampunit(static_cast<float>(temperature + 10) / 28.0f);

    return clampunit(1.0f - static_cast<float>(temperature - 18) / 22.0f);
}

float rainfertilityfactor(int rainfall)
{
    if (rainfall <= 75)
        return 0.0f;

    if (rainfall <= 700)
        return clampunit(static_cast<float>(rainfall - 75) / 625.0f);

    if (rainfall <= 1800)
        return 1.0f;

    return clampunit(1.0f - static_cast<float>(rainfall - 1800) / 1200.0f);
}
}

void generate_surface_potentials(planet& world)
{
    const int width = world.width();
    const int height = world.height();
    const int maximumflow = std::max(1, world.maxriverflow());
    const int majorriverthreshold = std::max(250, maximumflow / 20);

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int y = startrow; y <= endrow; y++)
        {
            for (int x = 0; x <= width; x++)
            {
                if (world.sea(x, y) != 0)
                {
                    world.seterosionpotential(x, y, 0);
                    world.setdepositionpotential(x, y, 0);
                    world.setfloodplainfertility(x, y, 0);
                    continue;
                }

                const float relief = static_cast<float>(localrelief(world, x, y));
                const float reliefnorm = clampunit(relief / 2200.0f);
                const float lownorm = 1.0f - clampunit(relief / 1400.0f);
                const float rainfallnorm = clampunit(static_cast<float>(annualrainfall(world, x, y)) / 1000.0f);
                const float flownorm = std::log1p(static_cast<float>(averageriverflow(world, x, y))) / std::log1p(static_cast<float>(maximumflow));
                const float inflownorm = clampunit(static_cast<float>(inflowcount(world, x, y)) / 4.0f);
                const bool majorrivernear = hasfeatureinradius(world, x, y, 2, [&](int nx, int ny)
                {
                    return averageriverflow(world, nx, ny) >= majorriverthreshold || isdeltacell(world, nx, ny);
                });
                const bool wetlandnear = hasfeatureinradius(world, x, y, 1, [&](int nx, int ny)
                {
                    return iswetlandcell(world, nx, ny);
                });
                const bool lakemargin = islakecell(world, x, y) || hasfeatureinradius(world, x, y, 1, [&](int nx, int ny)
                {
                    return islakecell(world, nx, ny) || issaltfeature(world, nx, ny);
                });
                const bool terminalbasin = world.basinclass(x, y) == BasinClass::endorheic && (issaltfeature(world, x, y) || lakemargin);
                const float rivercorridor = std::max(flownorm, majorrivernear ? 0.65f : 0.0f);

                float erosionscore = 100.0f * (0.34f * reliefnorm + 0.30f * flownorm + 0.24f * rainfallnorm + 0.12f * inflownorm);

                if (world.mountainheight(x, y) > 0)
                    erosionscore += 8.0f;

                float depositionscore = 100.0f * (0.32f * lownorm + 0.27f * rivercorridor + 0.15f * rainfallnorm
                    + 0.12f * (wetlandnear ? 1.0f : 0.0f)
                    + 0.09f * (lakemargin ? 1.0f : 0.0f)
                    + 0.05f * (terminalbasin ? 1.0f : 0.0f));

                if (isdeltacell(world, x, y))
                    depositionscore += 20.0f;

                if (reliefnorm > 0.7f)
                    depositionscore *= 0.6f;

                const float fertilitytemperature = temperaturerangefactor(averagetemperature(world, x, y));
                const float fertilityrain = rainfertilityfactor(annualrainfall(world, x, y));
                const float riverproximity = std::max(rivercorridor, majorrivernear ? 1.0f : 0.0f);
                const float wetlandfactor = wetlandnear ? 1.0f : (lakemargin ? 0.5f : 0.0f);
                const float gentleterrain = 1.0f - clampunit(relief / 1800.0f);
                const float coldpenalty = clampunit(static_cast<float>(world.glacialtemp() - world.maxtemp(x, y) + 2) / 14.0f);
                const float aridpenalty = clampunit(static_cast<float>(200 - annualrainfall(world, x, y)) / 200.0f);
                const float glacialpenalty = (world.biome(x, y) == biomeice || world.maxtemp(x, y) <= world.glacialtemp()) ? 1.0f : 0.0f;

                float fertilityscore = 0.42f * static_cast<float>(clampscore(depositionscore)) + 22.0f * riverproximity + 14.0f * fertilityrain
                    + 12.0f * fertilitytemperature + 10.0f * wetlandfactor;
                fertilityscore *= 0.55f + 0.45f * gentleterrain;
                fertilityscore -= 24.0f * coldpenalty + 20.0f * aridpenalty + 26.0f * glacialpenalty;

                world.seterosionpotential(x, y, clampscore(erosionscore));
                world.setdepositionpotential(x, y, clampscore(depositionscore));
                world.setfloodplainfertility(x, y, clampscore(fertilityscore));
            }
        }
    }, 16);
}
}
