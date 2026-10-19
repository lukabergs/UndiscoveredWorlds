#include "pipeline/export_routing.hpp"

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include "core/mesh/cell_grid.hpp"
#include "pipeline/climate_adapter.hpp"
#include "simulations/hydrology/drainage_basins.hpp"
#include "simulations/hydrology/flow_routing.hpp"
#include "simulations/hydrology/rivers.hpp"
#include "world/planet.hpp"
#include "world/planet_export.hpp"
#include "world/planet_fields.hpp"

using worldmodel::FieldId;

namespace
{
template<typename T>
bool anynonzero(const worldgrid<T>& grid)
{
    const T* values = grid.data();
    return std::any_of(values, values + grid.size(), [](T value) { return value != T{}; });
}

bool anyriverflow(const planet& world)
{
    for (mesh::CellId c = 0; c < world.cellcount(); c++)
        if (world.riverjan(c) != 0 || world.riverjul(c) != 0)
            return true;
    return false;
}

int exportcolumns(const mesh::CellGrid& grid)
{
    return grid.kind() == mesh::GridKind::cube ? 4 * grid.facecells() : grid.columns();
}

// The routing recomputed on an export keeps the world's open or closed answer for its lakes
// (sphere-mesh.md section 6). The export's Priority-Flood can spill above the world's lake
// level (the mean over a narrow outlet), which the closure rule (closeddepressioncells: the
// highest lake surface more than 0.01 m below the highest filled level) would read as a
// closed lake. An export depression none of whose lake cells lies on a closed cell of the
// world (the world cell at the lake cell's centre) is open: its lake cells below its spill
// level (its highest filled elevation) are raised to it. A depression with a lake cell on a
// closed world cell keeps its levels.
void keepopenlakesopen(planet& latlon, const planet& world)
{
    struct Level
    {
        float spill = 0.0f;
        bool lake = false;
        bool closed = false;
    };
    const std::vector<std::uint8_t> closed = physical_layers::closeddepressioncells(world);
    const mesh::CellGrid& grid = latlon.grid();
    const mesh::CellGrid& source = world.grid();
    const worldgrid<std::int32_t>& depression = latlon.depressionid();
    const worldgrid<float>& filled = latlon.filledelevation();
    std::unordered_map<std::int32_t, Level> levels;
    for (mesh::CellId c = 0; c < latlon.cellcount(); c++)
    {
        const std::int32_t id = depression[c];
        if (id == 0)
            continue;
        Level& level = levels.try_emplace(id, Level{ filled[c] }).first->second;
        level.spill = std::max(level.spill, filled[c]);
        if (latlon.lakesurfacefloat(c) == 0.0f)
            continue;
        level.lake = true;
        const mesh::CellId s = source.locate(grid.centre(c));
        if (s != mesh::nocell && closed[s] != 0)
            level.closed = true;
    }
    for (mesh::CellId c = 0; c < latlon.cellcount(); c++)
    {
        const std::int32_t id = depression[c];
        const float surface = latlon.lakesurfacefloat(c);
        if (id == 0 || surface == 0.0f)
            continue;
        const Level& level = levels.find(id)->second;
        if (!level.closed && surface < level.spill)
            latlon.setlakesurface(c, level.spill);
    }
}
}

const std::vector<FieldId>& exportroutingfields()
{
    static const std::vector<FieldId> fields = {
        FieldId::flow_direction, FieldId::drainage_area, FieldId::discharge, FieldId::filled_elevation,
        FieldId::depression_id, FieldId::runoff, FieldId::lake_surface,
        FieldId::river_direction, FieldId::river_flow_january, FieldId::river_flow_july,
        FieldId::delta_direction, FieldId::delta_flow_january, FieldId::delta_flow_july,
    };
    return fields;
}

void recomputeexportrouting(planet& latlon, const planet& world)
{
    const planet& view = latlon;
    const mesh::CellGrid& grid = view.grid();
    if (grid.kind() != mesh::GridKind::latlon)
        throw std::invalid_argument("recomputeexportrouting: the export is a lat-lon planet");
    // A world that was never routed has nothing to recompute.
    if (!anynonzero(view.drainagearea()))
        return;
    const bool rivers = anyriverflow(view);

    // hydrology.flow_routing on the exported filled surface (never below the exported terrain).
    const flowrouting::Inputs inputs = flowrouting::gatherinputs(view);
    std::vector<float> surface(inputs.elevation.size());
    const worldgrid<float>& filled = view.filledelevation();
    for (std::size_t cell = 0; cell < surface.size(); cell++)
        surface[cell] = std::max(filled[static_cast<mesh::CellId>(cell)], inputs.elevation[cell]);
    const flowrouting::Routing routing = flowrouting::route(inputs, &surface, &inputs.sea);
    flowrouting::writeroutingfields(latlon, inputs, routing);
    keepopenlakesopen(latlon, world);

    // hydrology.rivers on the recomputed routing, where the world has rivers.
    if (rivers)
        applyriverstage(latlon);
}

std::unique_ptr<planet> exportworld(const planet& world, const ClimateShadow* shadow)
{
    std::unique_ptr<planet> latlon = exportlatlon(world, exportcolumns(world.grid()), recomputeexportrouting);
    if (world.grid().kind() != mesh::GridKind::cube || shadow == nullptr || shadow->world() == nullptr)
        return latlon;
    const planet& climate = *shadow->world();
    if (climate.grid().kind() != mesh::GridKind::latlon || climate.grid().columns() != latlon->grid().columns())
        throw std::logic_error("exportworld: the climate shadow is not on the export's grid");
    for (std::size_t index = 0; index < worldmodel::kFieldCount; index++)
    {
        const FieldId id = static_cast<FieldId>(index);
        if (!climateshadowowned(id))
            continue;
        const worldmodel::FieldData data = readplanetfield(climate, id);
        if (!data.bytes.empty())
            loadplanetfield(*latlon, data);   // false for what a planet does not store (the solver's grids, levels, diagnostics)
    }
    return latlon;
}

std::unique_ptr<planet> exportbundleworld(const planet& world, const PlanetBundleSource& source)
{
    std::unique_ptr<planet> latlon = exportworld(world);
    if (world.grid().kind() != mesh::GridKind::cube || source.reader == nullptr)
        return latlon;
    const runbundle::Reader& reader = *source.reader;
    const std::vector<std::size_t> profileshape = { 6, static_cast<std::size_t>(latlon->grid().columns()) };
    const auto load = [&](const runbundle::ManifestEntry& entry, FieldId id)
    {
        PlanetFieldData data;
        data.id = id;
        data.shape = entry.shape;
        data.bytes = reader.fielddata(entry);
        loadplanetfield(*latlon, data);   // false for what a planet does not store
    };
    // The climate shadow's fields on the export's grid are the bundle's exports of them.
    for (const runbundle::ManifestEntry& entry : reader.exports())
    {
        FieldId id{};
        if (worldmodel::findfield(entry.id, id) && climateshadowowned(id))
            load(entry, id);
    }
    // Zonal profiles are stored in their lat-lon shape.
    for (const runbundle::ManifestEntry& entry : reader.fields())
    {
        FieldId id{};
        if (worldmodel::findfield(entry.id, id) && worldmodel::fieldspec(id).layout == worldmodel::Layout::zonalprofile &&
            entry.shape == profileshape)
            load(entry, id);
    }
    return latlon;
}

PlanetBundleExport bundleexport(const ClimateShadow* shadow)
{
    PlanetBundleExport sources;
    sources.recompute = recomputeexportrouting;
    sources.recomputed = exportroutingfields();
    if (shadow != nullptr && shadow->world() != nullptr)
    {
        sources.climateshadow = shadow->world();
        sources.shadowowned = climateshadowowned;
    }
    return sources;
}

void saveworldfile(planet& world, const std::string& filename, const ClimateShadow* shadow)
{
    if (world.grid().kind() != mesh::GridKind::cube)
    {
        world.saveworld(filename);
        return;
    }
    exportworld(world, shadow)->saveworld(filename);
}
