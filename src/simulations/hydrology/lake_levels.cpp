#include "simulations/hydrology/lake_levels.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <utility>
#include <vector>

#include "core/climate_seasons.hpp"
#include "simulations/hydrology/flow_routing.hpp"
#include "simulations/landscape/depressions.hpp"
#include "simulations/landscape/lake_balance.hpp"
#include "simulations/landscape/landscape_parameters.hpp"
#include "world/planet.hpp"

namespace
{
constexpr double secondsperyear = 365.25 * 86400.0;
constexpr int saltlake = 100;
constexpr int saltpan = 110;
}

void applylakelevelstage(planet& world)
{
    const planet& view = world;  // reads through the const getters, which the access tracking sees
    flowrouting::Inputs routing = flowrouting::gatherinputs(view);
    const landscape::SphereMetric& metric = routing.metric;
    const std::size_t size = metric.size();

    const worldgrid<float>& filled = view.filledelevation();
    const worldgrid<std::int32_t>& depression = view.depressionid();
    const worldgrid<std::uint8_t>& direction = view.flowdirection();
    const worldgrid<float>& discharge = view.discharge();

    landscape::LakeBalanceInputs inputs;
    inputs.metric = &metric;
    inputs.elevation = std::move(routing.elevation);
    inputs.runoff_m_yr = std::move(routing.runoff_m_yr);
    // The routing's potential evaporation: the open-water law per season, averaged (one law
    // for the runoff, the transmission losses and the lakes).
    inputs.evaporation_m_yr = std::move(routing.evaporation_m_yr);
    inputs.filledelevation.resize(size);
    inputs.depressionid.resize(size);
    inputs.receiver.resize(size);
    inputs.discharge_m3_yr.resize(size);
    inputs.precipitation_m_yr.resize(size);
    // The metric's cell index is the world's cell id (both in the grid's id order).
    for (std::size_t cell = 0; cell < size; cell++)
    {
        const mesh::CellId c = static_cast<mesh::CellId>(cell);
        inputs.filledelevation[cell] = filled[c];
        inputs.depressionid[cell] = depression[c];
        inputs.discharge_m3_yr[cell] = static_cast<double>(discharge[c]) * secondsperyear;
        int32_t receiver = -1;
        const int code = direction[c];
        if (code >= 1 && code <= landscape::directioncount)
        {
            const std::size_t neighbour = metric.neighbour(cell, code);
            if (neighbour != cell && routing.sea[neighbour] == 0)
                receiver = static_cast<int32_t>(neighbour);
        }
        inputs.receiver[cell] = receiver;
        // Seasonal precipitation is mm/month; four seasons of three months make the year.
        double rain_mm = 0.0;
        for (int season = 0; season < CLIMATESEASONCOUNT; season++)
            rain_mm += view.seasonalrainfloat(season, c);
        inputs.precipitation_m_yr[cell] = static_cast<float>(3.0 * rain_mm / 1000.0);
    }

    // The routing fields come from the last routing (hydrology.flow_routing, or the final routing of
    // terrain.fluvial_erosion). In every recipe the stage follows hydrology.flow_routing directly, so
    // the guard below drops nothing; it stays against a recipe that edits the terrain in between. The
    // sea here is the routing's (the world ocean, erosion plan H3): an enclosed below-sea basin is a
    // depression, not sea. A depression with a cell that is now sea, borders the sea or is not below
    // its filled level no longer matches the routing (a fresh Priority-Flood never puts a depression
    // next to the sea): it is dropped and reported, so no lake is written next to the sea.
    std::vector<uint8_t> stale(size, 0);
    for (std::size_t cell = 0; cell < size; cell++)
    {
        const int32_t id = inputs.depressionid[cell];
        if (id == 0)
            continue;
        bool consistent = routing.sea[cell] == 0 && inputs.elevation[cell] < inputs.filledelevation[cell];
        for (int code = 1; consistent && code <= landscape::directioncount; code++)
        {
            const std::size_t neighbour = metric.neighbour(cell, code);
            if (neighbour != cell && routing.sea[neighbour] != 0)
                consistent = false;
        }
        if (!consistent)
            stale[cell] = 1;
    }
    std::vector<int32_t> staleids;
    for (std::size_t cell = 0; cell < size; cell++)
        if (stale[cell] != 0)
            staleids.push_back(inputs.depressionid[cell]);
    std::sort(staleids.begin(), staleids.end());
    for (int32_t& id : inputs.depressionid)
        if (id != 0 && std::binary_search(staleids.begin(), staleids.end(), id))
            id = 0;
    staleids.erase(std::unique(staleids.begin(), staleids.end()), staleids.end());
    if (!staleids.empty())
        std::cout << "Lake levels: " << staleids.size() << " depression(s) no longer match the routing and hold no lake\n";
    // A world without sea (landscape::terminalsink): the routing's depression at the lowest cell is
    // its terminal basin, whose lake is closed whatever its inflow (every other depression drains
    // into it).
    const int32_t sink = landscape::terminalsink(inputs.elevation, routing.sea);
    if (sink >= 0)
        inputs.terminalid = inputs.depressionid[static_cast<std::size_t>(sink)];

    // Imported lakes (Earth imports, --import-lakes; erosion plan C1): all zero in a procedural world.
    // Only the cells outside every depression count: a depression's water is the balance's.
    const worldgrid<std::uint8_t>& importedlake = view.importedlake();
    std::vector<uint8_t> imported(size, 0);
    bool anyimported = false;
    for (std::size_t cell = 0; cell < size; cell++)
    {
        if (importedlake[static_cast<mesh::CellId>(cell)] != 0 && inputs.depressionid[cell] == 0)
        {
            imported[cell] = 1;
            anyimported = true;
        }
    }

    // The stage owns the lake surfaces (every cell is rewritten below) and the salt lakes and pans.
    // procedural-world runs it before and after the full climate, so the previous run's features are
    // cleared first; no other stage of a recipe with this one writes special 100 or 110 before it.
    for (std::size_t cell = 0; cell < size; cell++)
    {
        const mesh::CellId c = static_cast<mesh::CellId>(cell);
        const int feature = view.special(c);
        if (feature == saltlake || feature == saltpan)
            world.setspecial(c, 0);
    }

    const landscape::LakeBalance balance = landscape::balancelakes(inputs);
    for (std::size_t cell = 0; cell < size; cell++)
        world.setlakesurface(static_cast<mesh::CellId>(cell), balance.surface[cell]);
    // A dry closed depression is a salt pan on its lowest cell. In every closed depression, dry or
    // holding a lake, each cell at or below the sea level that the lake does not cover is a salt pan
    // too (erosion plan H3: an enclosed below-sea basin is not sea; Death Valley, Qattara, the salt
    // flats around a terminal lake). A pan carries a lake surface at its own terrain, as the legacy pans do, so the arid features
    // keep away from it, basin_class sees a closed depression and planet::sea does not count it.
    // The terrain is above the datum's zero wherever the sea level (at least 1) is above it, so a
    // pan surface is not the 0 that means no lake, except on a cell at the datum floor itself.
    std::vector<uint8_t> belowsea(size);
    for (std::size_t cell = 0; cell < size; cell++)
        belowsea[cell] = flowrouting::belowsea(inputs.elevation[cell], routing.sealevel) ? 1 : 0;
    const std::vector<uint8_t> pan = landscape::saltpans(inputs, balance, belowsea);
    for (std::size_t cell = 0; cell < size; cell++)
    {
        if (pan[cell] == 0)
            continue;
        world.setlakesurface(static_cast<mesh::CellId>(cell), inputs.elevation[cell]);
        world.setspecial(static_cast<mesh::CellId>(cell), saltpan);
    }
    for (std::size_t cell = 0; cell < size; cell++)
    {
        const int32_t id = inputs.depressionid[cell];
        if (id == 0 || balance.surface[cell] == 0.0f)
            continue;
        // balance.lakes is in ascending id.
        const auto lake = std::lower_bound(balance.lakes.begin(), balance.lakes.end(), id,
                                           [](const landscape::LakeLevel& entry, int32_t value) { return entry.id < value; });
        if (lake != balance.lakes.end() && lake->id == id && !lake->overflowing)
            world.setspecial(static_cast<mesh::CellId>(cell), saltlake);
    }

    // Imported lake cells inside a depression (depression_id != 0, after the stale ones are dropped)
    // are left to the balance, whatever it decides: an overflowing lake, a closed one or a dry
    // depression. Imported water is written only on the cells outside depressions, where the balance
    // has no say: each 8-connected group of them is water at one level. The import's terrain is the
    // lake bed (ETOPO keeps the bathymetry of the large lakes), not the lake surface, so the level is
    // the highest filled elevation of the group's cells (outside a depression the filled elevation is
    // the terrain): water fills to the lowest rim, the highest bed cell of a draining bed, and no
    // cell of the group lies above that level. The cells are fresh water with no special feature; a
    // cell the balance already flooded (the flat of a depression's spill) keeps the balance's level.
    if (!anyimported)
        return;
    std::vector<uint8_t> visited(size, 0);
    std::vector<std::size_t> members;
    for (std::size_t start = 0; start < size; start++)
    {
        if (imported[start] == 0 || visited[start] != 0)
            continue;
        members.assign(1, start);
        visited[start] = 1;
        float level = std::max(inputs.filledelevation[start], inputs.elevation[start]);
        for (std::size_t next = 0; next < members.size(); next++)
        {
            const std::size_t cell = members[next];
            level = std::max(level, std::max(inputs.filledelevation[cell], inputs.elevation[cell]));
            for (int code = 1; code <= landscape::directioncount; code++)
            {
                const std::size_t neighbour = metric.neighbour(cell, code);
                if (neighbour != cell && imported[neighbour] != 0 && visited[neighbour] == 0)
                {
                    visited[neighbour] = 1;
                    members.push_back(neighbour);
                }
            }
        }
        for (const std::size_t cell : members)
            if (balance.surface[cell] == 0.0f)
                world.setlakesurface(static_cast<mesh::CellId>(cell), level);
    }
}
