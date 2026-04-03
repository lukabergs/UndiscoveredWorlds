#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <vector>

#include "physical_layers.hpp"
#include "functions.hpp"

namespace
{
constexpr std::array<std::array<int, 2>, 8> drainageoffsets =
{ {
    { 0, -1 }, { 1, -1 }, { 1, 0 }, { 1, 1 },
    { 0, 1 }, { -1, 1 }, { -1, 0 }, { -1, -1 }
} };

enum class DrainageTerminal : std::uint8_t
{
    unknown = 255,
    exorheic = 1,
    endorheic = 2
};

int physicalindex(int x, int y, int width)
{
    return y * width + x;
}

twointegers wrappeddestination(const planet& world, int x, int y, int dir)
{
    twointegers destination = getdestination(x, y, dir);
    destination.x = wrap(destination.x, world.width());
    destination.y = std::clamp(destination.y, 0, world.height());
    return destination;
}

float clampunit(float value)
{
    return std::clamp(value, 0.0f, 1.0f);
}

int clampscore(float value)
{
    return static_cast<int>(std::round(std::clamp(value, 0.0f, 100.0f)));
}

bool islakecell(const planet& world, int x, int y)
{
    return world.truelake(x, y) != 0 || world.riftlakesurface(x, y) != 0;
}

bool iswatercell(const planet& world, int x, int y)
{
    return world.sea(x, y) != 0 || islakecell(world, x, y);
}

bool issaltfeature(const planet& world, int x, int y)
{
    const int special = world.special(x, y);
    return special == 100 || special == 110;
}

bool iswetlandcell(const planet& world, int x, int y)
{
    const int special = world.special(x, y);
    return special >= 130 && special < 140;
}

bool isdeltacell(const planet& world, int x, int y)
{
    return world.deltadir(x, y) != 0 || world.deltajan(x, y) != 0 || world.deltajul(x, y) != 0;
}

int annualrainfall(const planet& world, int x, int y)
{
    return (world.summerrain(x, y) + world.winterrain(x, y)) / 2;
}

int averagetemperature(const planet& world, int x, int y)
{
    return (world.mintemp(x, y) + world.maxtemp(x, y)) / 2;
}

int averageriverflow(const planet& world, int x, int y)
{
    const int riverflow = world.riveraveflow(x, y);
    const int deltaflow = (std::abs(world.deltajan(x, y)) + std::abs(world.deltajul(x, y))) / 2;
    return std::max(riverflow, deltaflow);
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

template<typename Predicate>
bool hasfeatureinradius(const planet& world, int x, int y, int radius, Predicate predicate)
{
    const int width = world.width();
    const int height = world.height();

    for (int dy = -radius; dy <= radius; dy++)
    {
        const int ny = y + dy;

        if (ny < 0 || ny > height)
            continue;

        for (int dx = -radius; dx <= radius; dx++)
        {
            const int nx = wrap(x + dx, width);

            if (predicate(nx, ny))
                return true;
        }
    }

    return false;
}

int maxscoreinradius(const planet& world, int x, int y, int radius, const std::function<int(int, int)>& getter)
{
    const int width = world.width();
    const int height = world.height();
    int best = 0;

    for (int dy = -radius; dy <= radius; dy++)
    {
        const int ny = y + dy;

        if (ny < 0 || ny > height)
            continue;

        for (int dx = -radius; dx <= radius; dx++)
        {
            const int nx = wrap(x + dx, width);
            best = std::max(best, getter(nx, ny));
        }
    }

    return best;
}

bool findadjacentsea(const planet& world, int x, int y, twointegers& destination)
{
    const int width = world.width();
    const int height = world.height();

    for (const auto& offset : drainageoffsets)
    {
        const int ny = y + offset[1];

        if (ny < 0 || ny > height)
            continue;

        const int nx = wrap(x + offset[0], width);

        if (world.sea(nx, ny) != 0)
        {
            destination.x = nx;
            destination.y = ny;
            return true;
        }
    }

    return false;
}

bool strictdownhilldestination(const planet& world, int x, int y, twointegers& destination)
{
    const int width = world.width();
    const int height = world.height();
    const int currentelevation = world.nom(x, y);
    int bestelevation = currentelevation;
    bool bestiswater = false;
    bool found = false;

    for (const auto& offset : drainageoffsets)
    {
        const int ny = y + offset[1];

        if (ny < 0 || ny > height)
            continue;

        const int nx = wrap(x + offset[0], width);
        const int candidateelevation = world.nom(nx, ny);
        const bool candidatewater = iswatercell(world, nx, ny);

        if (candidateelevation < bestelevation
            || (candidatewater && candidateelevation <= currentelevation && (found == false || candidateelevation < bestelevation || bestiswater == false)))
        {
            bestelevation = candidateelevation;
            bestiswater = candidatewater;
            destination.x = nx;
            destination.y = ny;
            found = true;
        }
    }

    return found && (bestelevation < currentelevation || bestiswater);
}

bool nextdrainagedestination(const planet& world, int x, int y, twointegers& destination)
{
    if (world.sea(x, y) != 0)
        return false;

    if (world.deltadir(x, y) > 0)
    {
        destination = wrappeddestination(world, x, y, world.deltadir(x, y));
        return true;
    }

    if (world.riverdir(x, y) != 0)
    {
        destination = wrappeddestination(world, x, y, world.riverdir(x, y));
        return true;
    }

    if ((world.outline(x, y) || isdeltacell(world, x, y) || islakecell(world, x, y)) && findadjacentsea(world, x, y, destination))
        return true;

    return strictdownhilldestination(world, x, y, destination);
}

bool reachesseaquickly(const planet& world, int x, int y, int maxsteps)
{
    int currentx = x;
    int currenty = y;

    for (int step = 0; step < maxsteps; step++)
    {
        if (world.sea(currentx, currenty) != 0 || world.outline(currentx, currenty))
            return true;

        twointegers destination;

        if (nextdrainagedestination(world, currentx, currenty, destination) == false)
            return false;

        if (destination.x == currentx && destination.y == currenty)
            return false;

        currentx = destination.x;
        currenty = destination.y;
    }

    return world.sea(currentx, currenty) != 0 || world.outline(currentx, currenty);
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

float regimevolcanicbias(GeologicRegime regime)
{
    switch (regime)
    {
    case GeologicRegime::convergent_arc:
        return 1.0f;
    case GeologicRegime::divergent_rift:
        return 0.9f;
    case GeologicRegime::mid_ocean_ridge:
        return 0.85f;
    case GeologicRegime::continent_collision:
        return 0.3f;
    case GeologicRegime::trench_adjacent:
        return 0.25f;
    default:
        return 0.0f;
    }
}

float regimeorebias(GeologicRegime regime)
{
    switch (regime)
    {
    case GeologicRegime::continent_collision:
        return 1.0f;
    case GeologicRegime::convergent_arc:
        return 0.9f;
    case GeologicRegime::transform:
        return 0.55f;
    case GeologicRegime::divergent_rift:
        return 0.4f;
    case GeologicRegime::passive_margin:
        return 0.15f;
    default:
        return 0.0f;
    }
}

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
}

void generatephysicalworldlayers(planet& world, const std::vector<std::vector<bool>>& shelves)
{
    const int width = world.width();
    const int height = world.height();
    const int cellwidth = width + 1;
    const int cellcount = cellwidth * (height + 1);
    const int maximumflow = std::max(1, world.maxriverflow());
    const int majorriverthreshold = std::max(250, maximumflow / 20);

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

    const float meanoceansst = sstcells > 0 ? static_cast<float>(ssttotal / static_cast<double>(sstcells)) : 0.0f;

    std::vector<std::uint8_t> drainage(cellcount, static_cast<std::uint8_t>(DrainageTerminal::unknown));
    std::vector<int> visitstamp(cellcount, 0);
    int currentstamp = 1;

    for (int y = 0; y <= height; y++)
    {
        for (int x = 0; x <= width; x++)
        {
            const int startindex = physicalindex(x, y, cellwidth);

            if (world.sea(x, y) != 0)
            {
                drainage[startindex] = static_cast<std::uint8_t>(DrainageTerminal::exorheic);
                continue;
            }

            if (drainage[startindex] != static_cast<std::uint8_t>(DrainageTerminal::unknown))
                continue;

            std::vector<int> path;
            path.reserve(32);
            int currentx = x;
            int currenty = y;
            DrainageTerminal terminal = DrainageTerminal::endorheic;

            while (true)
            {
                const int index = physicalindex(currentx, currenty, cellwidth);

                if (world.sea(currentx, currenty) != 0)
                {
                    terminal = DrainageTerminal::exorheic;
                    break;
                }

                if (drainage[index] != static_cast<std::uint8_t>(DrainageTerminal::unknown))
                {
                    terminal = static_cast<DrainageTerminal>(drainage[index]);
                    break;
                }

                if (visitstamp[index] == currentstamp)
                {
                    terminal = DrainageTerminal::endorheic;
                    break;
                }

                visitstamp[index] = currentstamp;
                path.push_back(index);

                twointegers destination;

                if (nextdrainagedestination(world, currentx, currenty, destination))
                {
                    if (destination.x == currentx && destination.y == currenty)
                    {
                        terminal = DrainageTerminal::endorheic;
                        break;
                    }

                    currentx = destination.x;
                    currenty = destination.y;
                    continue;
                }

                terminal = world.outline(currentx, currenty) || isdeltacell(world, currentx, currenty) ? DrainageTerminal::exorheic : DrainageTerminal::endorheic;
                break;
            }

            for (int index : path)
                drainage[index] = static_cast<std::uint8_t>(terminal);

            currentstamp++;
        }
    }

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int y = startrow; y <= endrow; y++)
        {
            for (int x = 0; x <= width; x++)
            {
                if (world.sea(x, y) != 0)
                {
                    world.setbasinclass(x, y, BasinClass::none);
                    continue;
                }

                const DrainageTerminal terminal = static_cast<DrainageTerminal>(drainage[physicalindex(x, y, cellwidth)]);
                const bool coastal = terminal == DrainageTerminal::exorheic
                    && world.map(x, y) <= world.sealevel() + 250
                    && (world.outline(x, y) || reachesseaquickly(world, x, y, 3));

                if (coastal)
                    world.setbasinclass(x, y, BasinClass::coastal);
                else if (terminal == DrainageTerminal::exorheic)
                    world.setbasinclass(x, y, BasinClass::exorheic);
                else
                    world.setbasinclass(x, y, BasinClass::endorheic);
            }
        }
    }, 16);

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

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int y = startrow; y <= endrow; y++)
        {
            for (int x = 0; x <= width; x++)
            {
                const GeologicRegime regime = world.geologicregime(x, y);
                const float regimemult = regimevolcanicbias(regime);
                const float divergencenorm = static_cast<float>(world.tectonicdivergence(x, y)) / 100.0f;
                const float convergencenorm = static_cast<float>(world.tectonicconvergence(x, y)) / 100.0f;
                const float ridgebias = (world.oceanridges(x, y) != 0 || world.oceanrifts(x, y) != 0) ? 1.0f : 0.0f;
                const float volcanobias = (world.volcano(x, y) != 0 || world.strato(x, y)) ? 1.0f : 0.0f;
                const float nearbyvolcano = maxscoreinradius(world, x, y, 2, [&](int nx, int ny)
                {
                    return world.volcano(nx, ny) != 0 || world.strato(nx, ny) ? 100 : 0;
                }) / 100.0f;

                float score = 100.0f * (0.34f * regimemult + 0.20f * volcanobias + 0.16f * nearbyvolcano
                    + 0.16f * divergencenorm + 0.14f * convergencenorm);
                score += 12.0f * ridgebias;

                if (world.sea(x, y) != 0 && ridgebias == 0.0f && regimemult < 0.6f)
                    score *= 0.5f;

                world.setvolcanicreserve(x, y, clampscore(score));
            }
        }
    }, 16);

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int y = startrow; y <= endrow; y++)
        {
            for (int x = 0; x <= width; x++)
            {
                const GeologicRegime regime = world.geologicregime(x, y);
                const float regimebias = regimeorebias(regime);
                const float convergencenorm = static_cast<float>(world.tectonicconvergence(x, y)) / 100.0f;
                const float shearnorm = static_cast<float>(world.tectonicshear(x, y)) / 100.0f;
                const float volcanismnear = static_cast<float>(maxscoreinradius(world, x, y, 2, [&](int nx, int ny)
                {
                    return world.volcanicreserve(nx, ny);
                })) / 100.0f;
                const float mountainbias = clampunit(static_cast<float>(world.mountainheight(x, y)) / 5000.0f);

                float score = 100.0f * (0.30f * regimebias + 0.22f * convergencenorm + 0.16f * shearnorm
                    + 0.16f * mountainbias + 0.16f * volcanismnear);

                if (world.sea(x, y) != 0 && world.coast(x, y) == 0)
                    score *= 0.35f;

                world.setmetalorereserve(x, y, clampscore(score));
            }
        }
    }, 16);

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int y = startrow; y <= endrow; y++)
        {
            for (int x = 0; x <= width; x++)
            {
                const float depositionnorm = static_cast<float>(world.depositionpotential(x, y)) / 100.0f;
                const float flownorm = std::log1p(static_cast<float>(averageriverflow(world, x, y))) / std::log1p(static_cast<float>(maximumflow));
                const int upstreamore = maxscoreinradius(world, x, y, 3, [&](int nx, int ny)
                {
                    if (world.map(nx, ny) < world.map(x, y))
                        return 0;

                    return world.metalorereserve(nx, ny);
                });
                const float upstreamorenorm = static_cast<float>(upstreamore) / 100.0f;

                float score = 100.0f * (0.35f * upstreamorenorm + 0.32f * depositionnorm + 0.20f * flownorm + 0.13f * (isdeltacell(world, x, y) ? 1.0f : 0.0f));

                if (world.sea(x, y) != 0 && isdeltacell(world, x, y) == false)
                    score *= 0.3f;

                world.setplacerreserve(x, y, clampscore(score));
            }
        }
    }, 16);

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int y = startrow; y <= endrow; y++)
        {
            for (int x = 0; x <= width; x++)
            {
                if (world.sea(x, y) != 0)
                {
                    world.setevaporitereserve(x, y, 0);
                    continue;
                }

                const float aridity = 1.0f - clampunit(static_cast<float>(annualrainfall(world, x, y)) / 350.0f);
                const float endorheic = world.basinclass(x, y) == BasinClass::endorheic ? 1.0f : 0.0f;
                const float saltbias = issaltfeature(world, x, y) ? 1.0f : 0.0f;
                const float terminalwater = (islakecell(world, x, y) || issaltfeature(world, x, y)) ? 1.0f : 0.0f;

                float score = 100.0f * (0.45f * endorheic + 0.25f * aridity + 0.20f * saltbias + 0.10f * terminalwater);

                world.setevaporitereserve(x, y, clampscore(score));
            }
        }
    }, 16);

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
