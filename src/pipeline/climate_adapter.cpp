#include "pipeline/climate_adapter.hpp"

#include <chrono>
#include <stdexcept>
#include <string>
#include <string_view>

#include "core/mesh/cell_grid.hpp"
#include "core/world/world_fields.hpp"
#include "pipeline/world_processes.hpp"
#include "simulations/climate/climate_context.hpp"
#include "simulations/climate/climate_world_fields.hpp"
#include "world/planet.hpp"
#include "world/planet_export.hpp"
#include "world/planet_fields.hpp"
#include "world/world_setup.hpp"

using worldmodel::FieldData;
using worldmodel::FieldId;
using worldmodel::Layout;

namespace
{
std::size_t indexof(FieldId id)
{
    return static_cast<std::size_t>(id);
}

bool isclimateprocess(const WorldProcessDefinition& process)
{
    return std::string_view(process.id).rfind("climate.", 0) == 0;
}

// Counts of cells and the grid's shape indices depend on the grid: the shadow
// keeps its own (the counts are recounted on its terrain).
bool shadowderived(FieldId id)
{
    return id == FieldId::land_cells || id == FieldId::sea_cells || id == FieldId::width_index || id == FieldId::height_index;
}

// Climate-owned state a cube world cannot load: the solver's own grids, the
// atmosphere levels and the diagnostics of the last solve (loadclimateworldfield),
// and the zonal profiles, which live on climate's lat-lon world in their lat-lon
// shape (sphere-mesh.md section 7). It stays on the shadow, never exported or imported.
bool shadowonly(FieldId id)
{
    const Layout layout = worldmodel::fieldspec(id).layout;
    if (layout == Layout::levelseasonal || layout == Layout::zonalprofile || isclimatediagnosticfield(id))
        return true;
    switch (id)
    {
    case FieldId::solver_monthly_temperature:
    case FieldId::solver_monthly_precipitation:
    case FieldId::solver_surface_wind_east:
    case FieldId::solver_surface_wind_south:
    case FieldId::solver_upper_wind_east:
    case FieldId::solver_upper_wind_south:
    case FieldId::solver_low_wind_east:
    case FieldId::solver_low_wind_south:
    case FieldId::solver_sea_surface_temperature:
    case FieldId::solver_sea_ice:
        return true;
    default:
        return false;
    }
}

// As hashplanetfield, from data already read.
std::uint64_t hashof(const FieldData& data)
{
    const std::uint64_t hash = fnv1a64(data.shape.data(), data.shape.size() * sizeof(std::size_t));
    return fnv1a64(data.bytes.data(), data.bytes.size(), hash);
}

// The run's climate settings (set on the world by the options and a recipe step's
// overrides); the solver state stays the shadow's own.
void copyclimatesettings(const climatecontext::Context& from, climatecontext::Context& to)
{
    to.resolution = from.resolution;
    to.calibrateMeanTemperature = from.calibrateMeanTemperature;
    to.transientWeather = from.transientWeather;
    to.mixedLayerTropics = from.mixedLayerTropics;
    to.genericLevelPressuresHpa = from.genericLevelPressuresHpa;
    to.genericLevelReversalHpa = from.genericLevelReversalHpa;
    to.transportMethod = from.transportMethod;
    to.oceanHeatStepHours = from.oceanHeatStepHours;
    to.reuseOceanCirculation = from.reuseOceanCirculation;
    to.krylovOceanCirculation = from.krylovOceanCirculation;
    to.gatherOceanHeatFluxes = from.gatherOceanHeatFluxes;
    to.oceanHeatWorkers = from.oceanHeatWorkers;
    to.oceanCoastalScheme = from.oceanCoastalScheme;
    to.oceanCaptureDirectory = from.oceanCaptureDirectory;
    to.atmosphereSolver = from.atmosphereSolver;
    to.atmosphereCaptureDirectory = from.atmosphereCaptureDirectory;
    to.couplingIterationLimit = from.couplingIterationLimit;
}

double secondssince(std::chrono::steady_clock::time_point start)
{
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}
}

ClimateShadow::ClimateShadow() = default;
ClimateShadow::~ClimateShadow() = default;
ClimateShadow::ClimateShadow(ClimateShadow&&) noexcept = default;
ClimateShadow& ClimateShadow::operator=(ClimateShadow&&) noexcept = default;

ClimateShadow::ClimateShadow(const ClimateShadow& other)
    : shadow_(other.shadow_ != nullptr ? std::make_unique<planet>(*other.shadow_) : nullptr),
      source_(other.source_), hashes_(other.hashes_), tracked_(other.tracked_), last_(other.last_)
{
}

ClimateShadow& ClimateShadow::operator=(const ClimateShadow& other)
{
    if (this != &other)
    {
        ClimateShadow copy(other);
        *this = std::move(copy);
    }
    return *this;
}

void ClimateShadow::reset()
{
    *this = ClimateShadow();
}

void runclimateprocess(GenerationStageId stage, planet& world, ClimateShadow& state, const std::function<void(planet&)>& body)
{
    state.last_ = {};
    if (world.grid().kind() != mesh::GridKind::cube)
    {
        body(world);
        return;
    }

    const WorldProcessDefinition& process = getworldprocess(stage);
    ClimateShadowSync& sync = state.last_;
    const auto exportstart = std::chrono::steady_clock::now();
    if (state.shadow_ == nullptr || state.source_ != &world.grid())
    {
        state.remap_.reset();
        state.shadow_ = exportlatlon(world, 4 * world.grid().facecells());
        state.source_ = &world.grid();
        state.hashes_.assign(worldmodel::kFieldCount, 0);
        state.tracked_.assign(worldmodel::kFieldCount, 0);
        for (const WorldProcessDefinition& candidate : getworldprocesses())
            if (isclimateprocess(candidate))
                for (const FieldId id : candidate.reads)
                    state.tracked_[indexof(id)] = 1;
        for (std::size_t index = 0; index < worldmodel::kFieldCount; index++)
            if (state.tracked_[index] != 0)
                state.hashes_[index] = hashof(readplanetfield(world, static_cast<FieldId>(index)));
        getlandandseatotals(*state.shadow_);
        sync.created = true;
    }
    else
    {
        if (state.remap_ == nullptr)
            state.remap_ = std::make_unique<FieldRemap>(world.grid(), state.shadow_->grid());
        for (const FieldId id : process.reads)
        {
            if (shadowderived(id) || shadowonly(id))
                continue;
            const FieldData data = readplanetfield(world, id);
            const std::uint64_t hash = hashof(data);
            if (state.tracked_[indexof(id)] != 0 && hash == state.hashes_[indexof(id)])
                continue;
            state.tracked_[indexof(id)] = 1;
            state.hashes_[indexof(id)] = hash;
            bool loaded = false;
            if (worldmodel::fieldspec(id).layout == Layout::scalar)
                loaded = loadplanetfield(*state.shadow_, data);
            else if (state.remap_->ongrid(data))
                loaded = loadplanetfield(*state.shadow_, state.remap_->remap(data));
            else
                continue;   // not on the world grid: nothing of the cube's to send
            if (!loaded)
                throw std::runtime_error("Climate adapter: the shadow does not load " + worldmodel::qualifiedname(id));
            sync.exported++;
        }
        if (sync.exported > 0)
            getlandandseatotals(*state.shadow_);
    }
    copyclimatesettings(world.climatesimulation(), state.shadow_->climatesimulation());
    sync.exportseconds = secondssince(exportstart);

    body(*state.shadow_);

    const auto importstart = std::chrono::steady_clock::now();
    std::vector<FieldId> fields;
    for (const FieldId id : process.writes)
        if (!shadowderived(id) && !shadowonly(id))
            fields.push_back(id);
    if (!fields.empty())
    {
        importfields(world, *state.shadow_, fields);
        sync.importcalls = 1;
        sync.imported = fields.size();
        // What the cube now holds came from the shadow: not sent back unless changed on the cube.
        for (const FieldId id : fields)
            if (state.tracked_[indexof(id)] != 0)
                state.hashes_[indexof(id)] = hashof(readplanetfield(world, id));
    }
    sync.importseconds = secondssince(importstart);
}
