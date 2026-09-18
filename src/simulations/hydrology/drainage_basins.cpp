#include "drainage_basins.hpp"
#include "../resources/physical_layers_internal.hpp"
#include <cstdint>

namespace physical_layers
{
using namespace detail;

namespace
{
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

bool iswatercell(const planet& world, int x, int y)
{
    return world.sea(x, y) != 0 || islakecell(world, x, y);
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
}

void generate_drainage_basins(planet& world)
{
    const int width = world.width();
    const int height = world.height();
    const int cellwidth = width + 1;
    const int cellcount = cellwidth * (height + 1);

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
}
}
