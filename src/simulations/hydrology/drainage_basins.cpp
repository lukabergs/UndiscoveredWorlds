#include "simulations/hydrology/drainage_basins.hpp"
#include "core/parallel_rows.hpp"
#include "simulations/hydrology/flow_routing.hpp"
#include "simulations/hydrology/hydrology_cells.hpp"
#include "world/planet.hpp"
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <vector>

// On the world's cell grid: cells are ids, neighbours come from the neighbour table
// (slots 1..8 in the cell's own frame), a direction code is a slot of its cell.

namespace physical_layers
{
namespace
{
enum class DrainageTerminal : std::uint8_t
{
    unknown = 255,
    exorheic = 1,
    endorheic = 2
};

// Cells per worker of the parallel loops over cells.
constexpr int cellsperworker = 4096;

bool iswatercell(const planet& world, mesh::CellId c)
{
    return world.sea(c) != 0 || hydrologycells::islake(world, c);
}

bool findadjacentsea(const planet& world, mesh::CellId c, mesh::CellId& destination)
{
    const mesh::CellGrid& grid = world.grid();

    for (int slot = 1; slot <= mesh::slotcount; slot++)
    {
        const std::int32_t neighbour = grid.neighbour(c, slot);

        if (neighbour < 0)
            continue;

        if (world.sea(static_cast<mesh::CellId>(neighbour)) != 0)
        {
            destination = static_cast<mesh::CellId>(neighbour);
            return true;
        }
    }

    return false;
}

bool strictdownhilldestination(const planet& world, mesh::CellId c, mesh::CellId& destination)
{
    const mesh::CellGrid& grid = world.grid();
    const float currentelevation = world.nomfloat(c);
    int bestelevation = currentelevation;
    bool bestiswater = false;
    bool found = false;

    for (int slot = 1; slot <= mesh::slotcount; slot++)
    {
        const std::int32_t neighbour = grid.neighbour(c, slot);

        if (neighbour < 0)
            continue;

        const mesh::CellId candidate = static_cast<mesh::CellId>(neighbour);
        const float candidateelevation = world.nomfloat(candidate);
        const bool candidatewater = iswatercell(world, candidate);

        if (candidateelevation < bestelevation
            || (candidatewater && candidateelevation <= currentelevation && (found == false || candidateelevation < bestelevation || bestiswater == false)))
        {
            bestelevation = candidateelevation;
            bestiswater = candidatewater;
            destination = candidate;
            found = true;
        }
    }

    return found && (bestelevation < currentelevation || bestiswater);
}

// The next cell along the legacy drainage: the delta or river direction (the cell
// itself where the slot has no neighbour), the adjacent sea of a coast, delta or lake,
// else the strict downhill neighbour.
bool nextdrainagedestination(const planet& world, mesh::CellId c, mesh::CellId& destination)
{
    if (world.sea(c) != 0)
        return false;

    const mesh::CellGrid& grid = world.grid();

    if (world.deltadir(c) > 0)
    {
        destination = hydrologycells::neighbourorself(grid, c, world.deltadir(c));
        return true;
    }

    if (world.riverdir(c) != 0)
    {
        destination = hydrologycells::neighbourorself(grid, c, world.riverdir(c));
        return true;
    }

    if ((hydrologycells::outline(world, c) || hydrologycells::isdelta(world, c) || hydrologycells::islake(world, c)) &&
        findadjacentsea(world, c, destination))
        return true;

    return strictdownhilldestination(world, c, destination);
}

bool reachesseaquickly(const planet& world, mesh::CellId c, int maxsteps)
{
    mesh::CellId current = c;

    for (int step = 0; step < maxsteps; step++)
    {
        if (world.sea(current) != 0 || hydrologycells::outline(world, current))
            return true;

        mesh::CellId destination = current;

        if (nextdrainagedestination(world, current, destination) == false)
            return false;

        if (destination == current)
            return false;

        current = destination;
    }

    return world.sea(current) != 0 || hydrologycells::outline(world, current);
}

// Eroded worlds (erosion plan, step B3): paths follow flow_direction; the sea is the
// routing's (flowrouting::oceancells, erosion plan H3: an enclosed below-sea basin
// smaller than landscape.routing.ocean_minimum_area_fraction is not sea, so what
// drains into it is endorheic). A depression is closed when its highest lake surface lies
// more than this below its filled level (a dry closed depression's salt pan carries a lake
// surface at its terrain; a depression without any lake cell is not closed).
constexpr float closedlakemarginm = 0.01f;

// The receiver of a cell along flow_direction; false for a cell without a direction
// (a root the flood never reached) or one whose slot has no neighbour.
bool flowstep(const mesh::CellGrid& grid, const worldgrid<std::uint8_t>& direction, mesh::CellId c, mesh::CellId& next)
{
    const int code = direction[c];

    if (code < 1 || code > mesh::slotcount)
        return false;

    const std::int32_t neighbour = grid.neighbour(c, code);

    if (neighbour < 0)
        return false;

    next = static_cast<mesh::CellId>(neighbour);
    return true;
}

}

std::vector<std::uint8_t> closeddepressioncells(const planet& world)
{
    struct DepressionLevel
    {
        float filled = 0.0f;
        float lake = 0.0f;
        bool haslake = false;
    };

    const std::size_t cellcount = world.cellcount();
    const worldgrid<std::int32_t>& depression = world.depressionid();
    const worldgrid<float>& filled = world.filledelevation();
    std::unordered_map<std::int32_t, DepressionLevel> levels;

    for (std::size_t cell = 0; cell < cellcount; cell++)
    {
        const mesh::CellId c = static_cast<mesh::CellId>(cell);
        const std::int32_t id = depression[c];

        if (id == 0)
            continue;

        DepressionLevel& level = levels.try_emplace(id, DepressionLevel{ filled[c], 0.0f, false }).first->second;
        level.filled = std::max(level.filled, filled[c]);
        const float surface = world.lakesurfacefloat(c);

        if (surface != 0.0f)
        {
            level.lake = level.haslake ? std::max(level.lake, surface) : surface;
            level.haslake = true;
        }
    }

    std::vector<std::uint8_t> closed(cellcount, 0);

    for (std::size_t cell = 0; cell < cellcount; cell++)
    {
        const std::int32_t id = depression[static_cast<mesh::CellId>(cell)];

        if (id == 0)
            continue;

        const DepressionLevel& level = levels.find(id)->second;

        // A depression without any lake cell was dropped by hydrology.lake_levels (its routing
        // no longer matched the terrain): it is not counted as closed. Dry closed depressions
        // carry a salt pan whose lake surface is its terrain, so they are found here.
        if (level.haslake && level.lake < level.filled - closedlakemarginm)
            closed[cell] = 1;
    }

    return closed;
}

namespace
{
// Per cell: the terminal of its path along flow_direction, memoised (a cell takes its
// receiver's terminal unless it lies in a closed depression). A walk longer than the
// cell count is counted as exorheic.
std::vector<std::uint8_t> classifyerodedbasins(const planet& world, const std::vector<std::uint8_t>& ocean)
{
    const mesh::CellGrid& grid = world.grid();
    const std::size_t cellcount = world.cellcount();
    const worldgrid<std::uint8_t>& direction = world.flowdirection();
    const std::vector<std::uint8_t> closed = closeddepressioncells(world);
    std::vector<std::uint8_t> drainage(cellcount, static_cast<std::uint8_t>(DrainageTerminal::unknown));
    std::vector<mesh::CellId> path;

    for (std::size_t start = 0; start < cellcount; start++)
    {
        if (ocean[start] != 0)
        {
            drainage[start] = static_cast<std::uint8_t>(DrainageTerminal::exorheic);
            continue;
        }

        if (drainage[start] != static_cast<std::uint8_t>(DrainageTerminal::unknown))
            continue;

        path.clear();
        mesh::CellId current = static_cast<mesh::CellId>(start);
        DrainageTerminal terminal = DrainageTerminal::exorheic;

        while (true)
        {
            if (ocean[current] != 0)
            {
                terminal = DrainageTerminal::exorheic;
                break;
            }

            if (drainage[current] != static_cast<std::uint8_t>(DrainageTerminal::unknown))
            {
                terminal = static_cast<DrainageTerminal>(drainage[current]);
                break;
            }

            if (path.size() >= cellcount)
            {
                terminal = DrainageTerminal::exorheic;
                break;
            }

            path.push_back(current);

            if (closed[current] != 0)
            {
                terminal = DrainageTerminal::endorheic;
                break;
            }

            mesh::CellId next = current;

            if (flowstep(grid, direction, current, next) == false)
            {
                terminal = DrainageTerminal::endorheic;
                break;
            }

            current = next;
        }

        for (const mesh::CellId cell : path)
            drainage[cell] = static_cast<std::uint8_t>(terminal);
    }

    return drainage;
}

// The legacy coastal walk along flow_direction: the sea or a coast within maxsteps.
bool reachesseaalongflow(const planet& world, const std::vector<std::uint8_t>& ocean,
                         const worldgrid<std::uint8_t>& direction, mesh::CellId c, int maxsteps)
{
    const mesh::CellGrid& grid = world.grid();
    mesh::CellId current = c;

    for (int step = 0; step < maxsteps; step++)
    {
        if (ocean[current] != 0 || hydrologycells::outline(world, current))
            return true;

        mesh::CellId next = current;

        if (flowstep(grid, direction, current, next) == false)
            return false;

        current = next;
    }

    return ocean[current] != 0 || hydrologycells::outline(world, current);
}

void generateerodedbasins(planet& world)
{
    const planet& view = world;
    const std::size_t cellcount = view.cellcount();
    const std::vector<std::uint8_t> ocean = flowrouting::oceancells(view);
    const std::vector<std::uint8_t> drainage = classifyerodedbasins(view, ocean);
    const worldgrid<std::uint8_t>& direction = view.flowdirection();
    const float coastalceiling = static_cast<float>(view.sealevel()) + 250.0f;

    parallelforrows(0, static_cast<int>(cellcount) - 1, [&](int first, int last)
    {
        for (int cell = first; cell <= last; cell++)
        {
            const mesh::CellId c = static_cast<mesh::CellId>(cell);

            if (ocean[c] != 0)
            {
                world.setbasinclass(c, BasinClass::none);
                continue;
            }

            const DrainageTerminal terminal = static_cast<DrainageTerminal>(drainage[c]);
            const bool coastal = terminal == DrainageTerminal::exorheic
                && view.mapfloat(c) <= coastalceiling
                && (hydrologycells::outline(view, c) || reachesseaalongflow(view, ocean, direction, c, 3));

            if (coastal)
                world.setbasinclass(c, BasinClass::coastal);
            else if (terminal == DrainageTerminal::exorheic)
                world.setbasinclass(c, BasinClass::exorheic);
            else
                world.setbasinclass(c, BasinClass::endorheic);
        }
    }, cellsperworker);
}
}

void generate_drainage_basins(planet& world)
{
    if (static_cast<const planet&>(world).erosionduration() > 0.0f)
    {
        generateerodedbasins(world);
        return;
    }

    const std::size_t cellcount = world.cellcount();

    std::vector<std::uint8_t> drainage(cellcount, static_cast<std::uint8_t>(DrainageTerminal::unknown));
    std::vector<int> visitstamp(cellcount, 0);
    int currentstamp = 1;

    for (std::size_t start = 0; start < cellcount; start++)
    {
        const mesh::CellId startcell = static_cast<mesh::CellId>(start);

        if (world.sea(startcell) != 0)
        {
            drainage[start] = static_cast<std::uint8_t>(DrainageTerminal::exorheic);
            continue;
        }

        if (drainage[start] != static_cast<std::uint8_t>(DrainageTerminal::unknown))
            continue;

        std::vector<mesh::CellId> path;
        path.reserve(32);
        mesh::CellId current = startcell;
        DrainageTerminal terminal = DrainageTerminal::endorheic;

        while (true)
        {
            if (world.sea(current) != 0)
            {
                terminal = DrainageTerminal::exorheic;
                break;
            }

            if (drainage[current] != static_cast<std::uint8_t>(DrainageTerminal::unknown))
            {
                terminal = static_cast<DrainageTerminal>(drainage[current]);
                break;
            }

            if (visitstamp[current] == currentstamp)
            {
                terminal = DrainageTerminal::endorheic;
                break;
            }

            visitstamp[current] = currentstamp;
            path.push_back(current);

            mesh::CellId destination = current;

            if (nextdrainagedestination(world, current, destination))
            {
                if (destination == current)
                {
                    terminal = DrainageTerminal::endorheic;
                    break;
                }

                current = destination;
                continue;
            }

            terminal = hydrologycells::outline(world, current) || hydrologycells::isdelta(world, current)
                ? DrainageTerminal::exorheic : DrainageTerminal::endorheic;
            break;
        }

        for (const mesh::CellId cell : path)
            drainage[cell] = static_cast<std::uint8_t>(terminal);

        currentstamp++;
    }

    parallelforrows(0, static_cast<int>(cellcount) - 1, [&](int first, int last)
    {
        for (int cell = first; cell <= last; cell++)
        {
            const mesh::CellId c = static_cast<mesh::CellId>(cell);

            if (world.sea(c) != 0)
            {
                world.setbasinclass(c, BasinClass::none);
                continue;
            }

            const DrainageTerminal terminal = static_cast<DrainageTerminal>(drainage[c]);
            const bool coastal = terminal == DrainageTerminal::exorheic
                && world.mapfloat(c) <= world.sealevel() + 250
                && (hydrologycells::outline(world, c) || reachesseaquickly(world, c, 3));

            if (coastal)
                world.setbasinclass(c, BasinClass::coastal);
            else if (terminal == DrainageTerminal::exorheic)
                world.setbasinclass(c, BasinClass::exorheic);
            else
                world.setbasinclass(c, BasinClass::endorheic);
        }
    }, cellsperworker);
}
}
