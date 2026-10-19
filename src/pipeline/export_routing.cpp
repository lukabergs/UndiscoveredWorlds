#include "pipeline/export_routing.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "core/mesh/cell_grid.hpp"
#include "pipeline/climate_adapter.hpp"
#include "simulations/hydrology/flow_routing.hpp"
#include "simulations/landscape/depressions.hpp"
#include "simulations/landscape/distributaries.hpp"
#include "simulations/landscape/flow_graph.hpp"
#include "simulations/landscape/landscape_parameters.hpp"
#include "simulations/landscape/sphere_metric.hpp"
#include "world/planet.hpp"
#include "world/planet_export.hpp"
#include "world/planet_fields.hpp"

using worldmodel::FieldId;

namespace
{
constexpr double secondsperyear = 365.25 * 86400.0;
constexpr double pi = 3.14159265358979323846;
// hydrology.rivers' legacy flow law (rivers.cpp): 0.012 m/yr is 1 mm/month.
constexpr double metresperyearpermmpermonth = 0.012;

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
}

const std::vector<FieldId>& exportroutingfields()
{
    static const std::vector<FieldId> fields = {
        FieldId::flow_direction, FieldId::drainage_area, FieldId::discharge, FieldId::filled_elevation,
        FieldId::river_direction, FieldId::river_flow_january, FieldId::river_flow_july,
        FieldId::delta_direction, FieldId::delta_flow_january, FieldId::delta_flow_july,
    };
    return fields;
}

void recomputeexportrouting(planet& latlon)
{
    const planet& view = latlon;
    const mesh::CellGrid& grid = view.grid();
    if (grid.kind() != mesh::GridKind::latlon)
        throw std::invalid_argument("recomputeexportrouting: the export is a lat-lon planet");
    // A world that was never routed has nothing to recompute.
    if (!anynonzero(view.drainagearea()))
        return;

    const int columns = grid.columns(), rows = grid.rows();
    const double radius = static_cast<double>(view.radius());
    const landscape::SphereMetric metric(columns, rows, radius);
    const std::size_t size = metric.size();
    const float sealevel = static_cast<float>(view.sealevel());

    // The exported filled surface (never below the exported terrain) and runoff (mm/yr to m/yr).
    std::vector<float> surface(size), runoff(size);
    const worldgrid<float>& filled = view.filledelevation();
    const worldgrid<float>& runoffmm = view.runoff();
    for (std::size_t cell = 0; cell < size; cell++)
    {
        const mesh::CellId c = static_cast<mesh::CellId>(cell);
        surface[cell] = std::max(filled[c], view.nomfloat(c));
        runoff[cell] = std::max(0.0f, runoffmm[c]) / 1000.0f;
    }
    const std::vector<std::uint8_t> sea = flowrouting::oceancells(view);
    const landscape::Depressions depressions =
        landscape::filldepressions(metric, surface, sea, tuning::landscape::depressions::minimumDepthM);
    landscape::RoutingOptions options;
    options.resolveflats = tuning::landscape::routing::resolveFlats;
    landscape::FlowGraph graph = landscape::buildflowgraph(metric, surface, sea, sealevel, depressions, options);
    landscape::accumulateflow(metric, runoff, graph);

    // The native rivers' seasonal ratio, read before the routing fields are replaced.
    const bool rivers = anyriverflow(view);
    std::vector<double> weight(size, 0.0), januaryweight(size, 0.0), julyweight(size, 0.0);
    if (rivers)
    {
        for (std::size_t cell = 0; cell < size; cell++)
        {
            if (sea[cell] != 0)
                continue;
            const mesh::CellId c = static_cast<mesh::CellId>(cell);
            const double january = view.riverjan(c), july = view.riverjul(c);
            const double w = static_cast<double>(runoff[cell]) * metric.area_m2(cell);
            weight[cell] = w;
            januaryweight[cell] = w * (january + july > 0.0 ? 2.0 * january / (january + july) : 1.0);
            julyweight[cell] = w * (january + july > 0.0 ? 2.0 * july / (january + july) : 1.0);
        }
        for (std::size_t index = graph.order.size(); index-- > 0;)
        {
            const std::size_t cell = static_cast<std::size_t>(graph.order[index]);
            const std::size_t next = static_cast<std::size_t>(graph.receiver[cell]);
            if (next == cell || sea[next] != 0)
                continue;
            weight[next] += weight[cell];
            januaryweight[next] += januaryweight[cell];
            julyweight[next] += julyweight[cell];
        }
    }

    worldgrid<std::uint8_t>& direction = latlon.flowdirection();
    worldgrid<float>& area = latlon.drainagearea();
    worldgrid<float>& discharge = latlon.discharge();
    worldgrid<float>& filledout = latlon.filledelevation();
    for (std::size_t cell = 0; cell < size; cell++)
    {
        const mesh::CellId c = static_cast<mesh::CellId>(cell);
        direction[c] = graph.direction[cell];
        area[c] = static_cast<float>(graph.area_m2[cell] / 1.0e6);
        discharge[c] = static_cast<float>(graph.discharge_m3_yr[cell] / secondsperyear);
        filledout[c] = depressions.filledelevation[cell];
    }
    if (!rivers)
        return;

    // hydrology.rivers' legacy units: a mean-sized cell under 1 mm/month at the runoff coefficient, over river_factor.
    const double meanarea = 4.0 * pi * radius * radius / (static_cast<double>(columns) * static_cast<double>(rows));
    const double unitdischarge = tuning::landscape::runoff::coefficient * metresperyearpermmpermonth *
        static_cast<double>(view.riverfactor()) * meanarea;
    const double minimumflow = static_cast<double>(tuning::landscape::rivers::minimumFlow);
    for (std::size_t cell = 0; cell < size; cell++)
    {
        const mesh::CellId c = static_cast<mesh::CellId>(cell);
        int code = 0, january = 0, july = 0;
        if (sea[cell] == 0)
        {
            code = graph.direction[cell];
            const double flow = graph.discharge_m3_yr[cell] / unitdischarge;
            if (flow >= minimumflow)
            {
                const double januaryshare = weight[cell] > 0.0 ? januaryweight[cell] / weight[cell] : 1.0;
                const double julyshare = weight[cell] > 0.0 ? julyweight[cell] / weight[cell] : 1.0;
                january = static_cast<int>(std::lround(flow * januaryshare));
                july = static_cast<int>(std::lround(flow * julyshare));
            }
        }
        latlon.setriverdir(c, code);
        latlon.setriverjan(c, january);
        latlon.setriverjul(c, july);
    }

    // River mouths as hydrology.rivers writes them: the sea cell an outlet drains into takes
    // its direction (the largest outlet's) and the sum of the outlets' flows.
    std::vector<int> mouthbest(size, -1);
    for (std::size_t cell = 0; cell < size; cell++)
    {
        const int code = graph.direction[cell];
        if (sea[cell] != 0 || code < 1 || code > landscape::directioncount)
            continue;
        const std::size_t mouth = metric.neighbour(cell, code);
        if (mouth == cell || sea[mouth] == 0)
            continue;
        const mesh::CellId c = static_cast<mesh::CellId>(cell), m = static_cast<mesh::CellId>(mouth);
        const int january = view.riverjan(c), july = view.riverjul(c);
        if (january == 0 && july == 0)
            continue;
        if (january + july > mouthbest[mouth])
        {
            mouthbest[mouth] = january + july;
            latlon.setriverdir(m, code);
        }
        latlon.setriverjan(m, view.riverjan(m) + january);
        latlon.setriverjul(m, view.riverjul(m) + july);
    }

    if (view.erosionduration() > 0.0f)
    {
        std::vector<std::uint8_t> prograded(size), codes(size);
        std::vector<std::int32_t> january(size), july(size);
        const worldgrid<std::uint8_t>& progradedgrid = view.prograded();
        for (std::size_t cell = 0; cell < size; cell++)
        {
            const mesh::CellId c = static_cast<mesh::CellId>(cell);
            prograded[cell] = progradedgrid[c];
            codes[cell] = graph.direction[cell];
            january[cell] = view.riverjan(c);
            july[cell] = view.riverjul(c);
        }
        const landscape::Distributaries split = landscape::splitdistributaries(metric, sea, prograded,
            depressions.filledelevation, static_cast<double>(view.sealevel()), codes, january, july, minimumflow);
        for (std::size_t cell = 0; cell < size; cell++)
        {
            const mesh::CellId c = static_cast<mesh::CellId>(cell);
            latlon.setdeltadir(c, split.direction[cell]);
            latlon.setdeltajan(c, split.january[cell]);
            latlon.setdeltajul(c, split.july[cell]);
        }
    }
    latlon.setmaxriverflow();
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
