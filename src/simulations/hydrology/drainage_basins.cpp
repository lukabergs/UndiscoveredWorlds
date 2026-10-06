#include "simulations/hydrology/drainage_basins.hpp"
#include "simulations/resources/physical_layers_internal.hpp"
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <vector>

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
    const float currentelevation = world.nomfloat(x, y);
    int bestelevation = currentelevation;
    bool bestiswater = false;
    bool found = false;

    for (const auto& offset : drainageoffsets)
    {
        const int ny = y + offset[1];

        if (ny < 0 || ny > height)
            continue;

        const int nx = wrap(x + offset[0], width);
        const float candidateelevation = world.nomfloat(nx, ny);
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

// Eroded worlds (erosion plan, step B3): paths follow flow_direction; the sea is the
// routing's (basebelowsea). A depression is closed when its highest lake surface lies
// more than this below its filled level (a dry closed depression's salt pan carries a lake
// surface at its terrain; a depression without any lake cell is not closed).
constexpr float closedlakemarginm = 0.01f;

// The receiver of (x, y) along flow_direction; false for a cell without a direction
// (a root the flood never reached) or one pointing beyond a pole.
bool flowstep(const worldgrid<std::uint8_t>& direction, int width, int height, int x, int y, int& nextx, int& nexty)
{
    const int code = direction[x][y];

    if (code < 1 || code > static_cast<int>(drainageoffsets.size()))
        return false;

    const auto& offset = drainageoffsets[static_cast<std::size_t>(code - 1)];
    nexty = y + offset[1];

    if (nexty < 0 || nexty > height)
        return false;

    nextx = wrap(x + offset[0], width);
    return true;
}

// Per cell: 1 when the cell belongs to a closed depression.
std::vector<std::uint8_t> closeddepressioncells(const planet& world)
{
    struct DepressionLevel
    {
        float filled = 0.0f;
        float lake = 0.0f;
        bool haslake = false;
    };

    const int width = world.width();
    const int height = world.height();
    const int cellwidth = width + 1;
    const worldgrid<std::int32_t>& depression = world.depressionid();
    const worldgrid<float>& filled = world.filledelevation();
    std::unordered_map<std::int32_t, DepressionLevel> levels;

    for (int y = 0; y <= height; y++)
    {
        for (int x = 0; x <= width; x++)
        {
            const std::int32_t id = depression[x][y];

            if (id == 0)
                continue;

            DepressionLevel& level = levels.try_emplace(id, DepressionLevel{ filled[x][y], 0.0f, false }).first->second;
            level.filled = std::max(level.filled, filled[x][y]);
            const float surface = world.lakesurfacefloat(x, y);

            if (surface != 0.0f)
            {
                level.lake = level.haslake ? std::max(level.lake, surface) : surface;
                level.haslake = true;
            }
        }
    }

    std::vector<std::uint8_t> closed(static_cast<std::size_t>(cellwidth * (height + 1)), 0);

    for (int y = 0; y <= height; y++)
    {
        for (int x = 0; x <= width; x++)
        {
            const std::int32_t id = depression[x][y];

            if (id == 0)
                continue;

            const DepressionLevel& level = levels.find(id)->second;

            // A depression without any lake cell was dropped by hydrology.lake_levels (its routing
            // no longer matched the terrain): it is not counted as closed. Dry closed depressions
            // carry a salt pan whose lake surface is its terrain, so they are found here.
            if (level.haslake && level.lake < level.filled - closedlakemarginm)
                closed[static_cast<std::size_t>(physicalindex(x, y, cellwidth))] = 1;
        }
    }

    return closed;
}

// Per cell: the terminal of its path along flow_direction, memoised (a cell takes its
// receiver's terminal unless it lies in a closed depression). A walk longer than the
// cell count is counted as exorheic.
std::vector<std::uint8_t> classifyerodedbasins(const planet& world)
{
    const int width = world.width();
    const int height = world.height();
    const int cellwidth = width + 1;
    const int cellcount = cellwidth * (height + 1);
    const worldgrid<std::uint8_t>& direction = world.flowdirection();
    const std::vector<std::uint8_t> closed = closeddepressioncells(world);
    std::vector<std::uint8_t> drainage(static_cast<std::size_t>(cellcount), static_cast<std::uint8_t>(DrainageTerminal::unknown));
    std::vector<int> path;

    for (int y = 0; y <= height; y++)
    {
        for (int x = 0; x <= width; x++)
        {
            const int startindex = physicalindex(x, y, cellwidth);

            if (world.basebelowsea(x, y))
            {
                drainage[startindex] = static_cast<std::uint8_t>(DrainageTerminal::exorheic);
                continue;
            }

            if (drainage[startindex] != static_cast<std::uint8_t>(DrainageTerminal::unknown))
                continue;

            path.clear();
            int currentx = x;
            int currenty = y;
            DrainageTerminal terminal = DrainageTerminal::exorheic;

            while (true)
            {
                const int index = physicalindex(currentx, currenty, cellwidth);

                if (world.basebelowsea(currentx, currenty))
                {
                    terminal = DrainageTerminal::exorheic;
                    break;
                }

                if (drainage[index] != static_cast<std::uint8_t>(DrainageTerminal::unknown))
                {
                    terminal = static_cast<DrainageTerminal>(drainage[index]);
                    break;
                }

                if (static_cast<int>(path.size()) >= cellcount)
                {
                    terminal = DrainageTerminal::exorheic;
                    break;
                }

                path.push_back(index);

                if (closed[static_cast<std::size_t>(index)] != 0)
                {
                    terminal = DrainageTerminal::endorheic;
                    break;
                }

                int nextx = 0;
                int nexty = 0;

                if (flowstep(direction, width, height, currentx, currenty, nextx, nexty) == false)
                {
                    terminal = DrainageTerminal::endorheic;
                    break;
                }

                currentx = nextx;
                currenty = nexty;
            }

            for (int index : path)
                drainage[index] = static_cast<std::uint8_t>(terminal);
        }
    }

    return drainage;
}

// The legacy coastal walk along flow_direction: the sea or a coast within maxsteps.
bool reachesseaalongflow(const planet& world, const worldgrid<std::uint8_t>& direction, int x, int y, int maxsteps)
{
    const int width = world.width();
    const int height = world.height();
    int currentx = x;
    int currenty = y;

    for (int step = 0; step < maxsteps; step++)
    {
        if (world.basebelowsea(currentx, currenty) || world.outline(currentx, currenty))
            return true;

        int nextx = 0;
        int nexty = 0;

        if (flowstep(direction, width, height, currentx, currenty, nextx, nexty) == false)
            return false;

        currentx = nextx;
        currenty = nexty;
    }

    return world.basebelowsea(currentx, currenty) || world.outline(currentx, currenty);
}

void generateerodedbasins(planet& world)
{
    const planet& view = world;
    const int width = view.width();
    const int height = view.height();
    const int cellwidth = width + 1;
    const std::vector<std::uint8_t> drainage = classifyerodedbasins(view);
    const worldgrid<std::uint8_t>& direction = view.flowdirection();
    const float coastalceiling = static_cast<float>(view.sealevel()) + 250.0f;

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int y = startrow; y <= endrow; y++)
        {
            for (int x = 0; x <= width; x++)
            {
                if (view.basebelowsea(x, y))
                {
                    world.setbasinclass(x, y, BasinClass::none);
                    continue;
                }

                const DrainageTerminal terminal = static_cast<DrainageTerminal>(drainage[physicalindex(x, y, cellwidth)]);
                const bool coastal = terminal == DrainageTerminal::exorheic
                    && view.mapfloat(x, y) <= coastalceiling
                    && (view.outline(x, y) || reachesseaalongflow(view, direction, x, y, 3));

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

void generate_drainage_basins(planet& world)
{
    if (static_cast<const planet&>(world).erosionduration() > 0.0f)
    {
        generateerodedbasins(world);
        return;
    }

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
                    && world.mapfloat(x, y) <= world.sealevel() + 250
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
