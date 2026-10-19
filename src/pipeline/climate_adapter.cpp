#include "pipeline/climate_adapter.hpp"

#include <chrono>
#include <cstring>
#include <stdexcept>
#include <string>
#include <string_view>

#include "core/mesh/cell_grid.hpp"
#include "core/mesh/mesh_remap.hpp"
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
using worldmodel::RemapRule;

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

// Climate-owned state that stays on the shadow, never exported or imported: what a
// world cannot load (loadclimateworldfield: the solver's own grids, the atmosphere
// levels, the diagnostics of the last solve), the zonal profiles (they live on
// climate's lat-lon world in their lat-lon shape, sphere-mesh.md section 7) and the
// direction-coded fields (D8 codes are never remapped, section 5; the orographic
// rain directions are climate's own, so no routing on the cube recomputes them).
bool shadowonly(FieldId id)
{
    const worldmodel::FieldSpec& spec = worldmodel::fieldspec(id);
    if (spec.layout == Layout::levelseasonal || spec.layout == Layout::zonalprofile || spec.rule == RemapRule::direction ||
        isclimatediagnosticfield(id))
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

bool synced(FieldId id)
{
    return !shadowderived(id) && !shadowonly(id);
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

// The cube cells whose imported value can depend on a changed shadow cell, by the
// field's rule as the remap applies it (planet_export.cpp, mesh_remap.hpp): mean
// and total from the cell's overlap row (an empty mean row from its bilinear
// stencil), nearest from the shadow cell at the cube cell's centre.
std::vector<std::uint8_t> reachedcells(RemapRule rule, const std::vector<std::uint8_t>& changed, const mesh::RemapWeights& weights,
    const mesh::CellGrid& shadowgrid, const mesh::CellGrid& cubegrid)
{
    std::vector<std::uint8_t> reached(cubegrid.cellcount(), 0);
    for (mesh::CellId c = 0; c < cubegrid.cellcount(); c++)
    {
        if (rule == RemapRule::nearest)
        {
            reached[c] = changed[shadowgrid.locate(cubegrid.centre(c))];
            continue;
        }
        const std::size_t first = weights.rowstart[c], last = weights.rowstart[c + 1];
        for (std::size_t entry = first; entry < last && reached[c] == 0; entry++)
            reached[c] = changed[weights.source[entry]];
        if (first == last && rule == RemapRule::mean)
        {
            const mesh::BilinearStencil stencil = mesh::bilinearstencil(shadowgrid, cubegrid.centre(c));
            for (int n = 0; n < stencil.count && reached[c] == 0; n++)
                reached[c] = changed[stencil.cell[n]];
        }
    }
    return reached;
}
}

// The shadow-only fields, of the direction-coded ones only climate's own (the atmosphere's:
// the orographic rain directions); the routing's and the ridges' are not climate's.
bool climateshadowowned(FieldId id)
{
    if (!shadowonly(id))
        return false;
    const worldmodel::FieldSpec& spec = worldmodel::fieldspec(id);
    return spec.rule != RemapRule::direction || spec.component == worldmodel::Component::atmosphere;
}

ClimateShadow::ClimateShadow() = default;
ClimateShadow::~ClimateShadow() = default;
ClimateShadow::ClimateShadow(ClimateShadow&&) noexcept = default;
ClimateShadow& ClimateShadow::operator=(ClimateShadow&&) noexcept = default;

ClimateShadow::ClimateShadow(const ClimateShadow& other)
    : shadow_(other.shadow_ != nullptr ? std::make_unique<planet>(*other.shadow_) : nullptr),
      facecells_(other.facecells_), hashes_(other.hashes_), tracked_(other.tracked_),
      importweights_(other.importweights_), last_(other.last_)
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
    const mesh::CellGrid& cubegrid = world.grid();
    ClimateShadowSync& sync = state.last_;
    const auto exportstart = std::chrono::steady_clock::now();
    if (state.shadow_ != nullptr && state.facecells_ != cubegrid.facecells())
        throw std::logic_error("Climate adapter: the world's grid changed under its shadow; reset the shadow where the world is cleared");
    if (state.shadow_ == nullptr)
    {
        state.shadow_ = exportlatlon(world, 4 * cubegrid.facecells());
        state.facecells_ = cubegrid.facecells();
        state.importweights_.reset();
        state.hashes_.assign(worldmodel::kFieldCount, 0);
        state.tracked_.assign(worldmodel::kFieldCount, 0);
        for (const WorldProcessDefinition& candidate : getworldprocesses())
        {
            if (!isclimateprocess(candidate))
                continue;
            for (const FieldId id : candidate.reads)
                state.tracked_[indexof(id)] = synced(id) ? 1 : 0;
            for (const FieldId id : candidate.writes)
                state.tracked_[indexof(id)] = synced(id) ? 1 : 0;
        }
        for (std::size_t index = 0; index < worldmodel::kFieldCount; index++)
            if (state.tracked_[index] != 0)
                state.hashes_[index] = hashof(readplanetfield(world, static_cast<FieldId>(index)));
        getlandandseatotals(*state.shadow_);
        sync.created = true;
    }
    else
    {
        // Every field of climate's, not only this process's reads: a field it only
        // writes must start from the cube's value where another process changed it.
        const FieldRemap exportremap(cubegrid, state.shadow_->grid());
        for (std::size_t index = 0; index < worldmodel::kFieldCount; index++)
        {
            if (state.tracked_[index] == 0)
                continue;
            const FieldId id = static_cast<FieldId>(index);
            const FieldData data = readplanetfield(world, id);
            const std::uint64_t hash = hashof(data);
            if (hash == state.hashes_[index])
                continue;
            state.hashes_[index] = hash;
            bool loaded = false;
            if (worldmodel::fieldspec(id).layout == Layout::scalar)
                loaded = loadplanetfield(*state.shadow_, data);
            else if (exportremap.ongrid(data))
                loaded = loadplanetfield(*state.shadow_, exportremap.remap(data));
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

    // The shadow before the process: hashes of every synced field, the data of the declared writes.
    planet& shadow = *state.shadow_;
    std::vector<std::uint64_t> before(worldmodel::kFieldCount, 0);
    std::vector<FieldData> writesbefore(process.writes.size());
    for (std::size_t index = 0; index < worldmodel::kFieldCount; index++)
        if (synced(static_cast<FieldId>(index)))
            before[index] = hashof(readplanetfield(shadow, static_cast<FieldId>(index)));
    for (std::size_t write = 0; write < process.writes.size(); write++)
        if (synced(process.writes[write]))
            writesbefore[write] = readplanetfield(shadow, process.writes[write]);
    sync.exportseconds = secondssince(exportstart);

    body(shadow);

    const auto importstart = std::chrono::steady_clock::now();
    std::vector<bool> declared(worldmodel::kFieldCount, false);
    for (const FieldId id : process.writes)
        declared[indexof(id)] = true;
    std::string undeclared;
    for (std::size_t index = 0; index < worldmodel::kFieldCount; index++)
    {
        const FieldId id = static_cast<FieldId>(index);
        if (!synced(id) || declared[index] || hashof(readplanetfield(shadow, id)) == before[index])
            continue;
        undeclared += (undeclared.empty() ? "" : ", ") + worldmodel::qualifiedname(id);
    }
    if (!undeclared.empty())
        throw std::runtime_error(std::string("Climate adapter: ") + process.id + " changed fields it does not declare "
            "(world_process_access.inc): " + undeclared);

    // The declared writes the process changed, and the cube's values before the import.
    std::vector<FieldId> fields;
    std::vector<std::size_t> writeindex;
    std::vector<FieldData> cubebefore;
    for (std::size_t write = 0; write < process.writes.size(); write++)
    {
        const FieldId id = process.writes[write];
        if (!synced(id) || hashof(readplanetfield(shadow, id)) == hashof(writesbefore[write]))
            continue;
        fields.push_back(id);
        writeindex.push_back(write);
        cubebefore.push_back(readplanetfield(world, id));
    }
    if (!fields.empty())
    {
        importfields(world, shadow, fields);
        sync.importcalls = 1;
        sync.imported = fields.size();

        // Cube cells no changed shadow cell reaches keep their values.
        const mesh::CellGrid& shadowgrid = shadow.grid();
        const std::size_t shadowcells = shadowgrid.cellcount(), cubecells = cubegrid.cellcount();
        for (std::size_t field = 0; field < fields.size(); field++)
        {
            const FieldId id = fields[field];
            const worldmodel::FieldSpec& spec = worldmodel::fieldspec(id);
            if (spec.layout == Layout::scalar)
                continue;
            const FieldData& old = cubebefore[field];
            FieldData updated = readplanetfield(world, id);
            const FieldData& shadowold = writesbefore[writeindex[field]];
            const FieldData shadownew = readplanetfield(shadow, id);
            const std::size_t size = worldmodel::elementsize(spec.element);
            const std::size_t layers = updated.bytes.size() / (cubecells * size);
            if (old.bytes.size() != updated.bytes.size() || shadowold.bytes.size() != shadownew.bytes.size() ||
                shadownew.bytes.size() != layers * shadowcells * size)
                continue;   // the field was not on the cube before: the whole import stands
            if (spec.rule != RemapRule::nearest && state.importweights_ == nullptr)
                state.importweights_ = std::make_shared<const mesh::RemapWeights>(mesh::makeremap(shadowgrid, cubegrid));
            std::size_t kept = 0;
            std::vector<std::uint8_t> changed(shadowcells);
            for (std::size_t layer = 0; layer < layers; layer++)
            {
                const std::uint8_t* a = shadowold.bytes.data() + layer * shadowcells * size;
                const std::uint8_t* b = shadownew.bytes.data() + layer * shadowcells * size;
                for (std::size_t s = 0; s < shadowcells; s++)
                    changed[s] = std::memcmp(a + s * size, b + s * size, size) != 0 ? 1 : 0;
                const std::vector<std::uint8_t> reached = reachedcells(spec.rule, changed,
                    spec.rule == RemapRule::nearest ? mesh::RemapWeights{} : *state.importweights_, shadowgrid, cubegrid);
                for (std::size_t c = 0; c < cubecells; c++)
                {
                    if (reached[c] != 0)
                        continue;
                    const std::size_t offset = (layer * cubecells + c) * size;
                    if (std::memcmp(updated.bytes.data() + offset, old.bytes.data() + offset, size) != 0)
                    {
                        std::memcpy(updated.bytes.data() + offset, old.bytes.data() + offset, size);
                        kept++;
                    }
                }
            }
            if (kept > 0 && !loadplanetfield(world, updated))
                throw std::runtime_error("Climate adapter: the world does not load " + worldmodel::qualifiedname(id));
            sync.keptcells += kept;
        }
        // What the cube now holds came from the shadow: not sent back unless changed on the cube.
        for (const FieldId id : fields)
            if (state.tracked_[indexof(id)] != 0)
                state.hashes_[indexof(id)] = hashof(readplanetfield(world, id));
    }
    sync.importseconds = secondssince(importstart);
}
