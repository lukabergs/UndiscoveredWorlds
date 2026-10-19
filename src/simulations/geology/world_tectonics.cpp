#include "simulations/geology/world_tectonics.hpp"
#include "validation/tectonics/tectonic_snapshots.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <memory>
#include <optional>
#include <utility>
#include <vector>

#include "simulations/geology/mineral_systems/mineral_history.hpp"
#include "simulations/geology/mineral_systems/tectonic_frame_reader.hpp"
#include "simulations/geology/plate_tectonics/tectonic_api.hpp"
#include "simulations/geology/plate_tectonics/tectonic_contract.hpp"
#include "simulations/geology/plate_tectonics/tectonic_recipe.hpp"
#include "simulations/geology/plate_tectonics/topography_codec.hpp"
#include "core/mesh/cell_grid.hpp"
#include "core/mesh/mesh_remap.hpp"
#ifdef max
#undef max
#endif
#ifdef min
#undef min
#endif
#include "world/functions.hpp"
#include "world/planet.hpp"

namespace
{
constexpr auto kVisualizationInterval = std::chrono::milliseconds(50);
struct PlateTectonicsHandle
{
    void* pointer = nullptr;

    PlateTectonicsHandle() = default;
    PlateTectonicsHandle(const PlateTectonicsHandle&) = delete;
    PlateTectonicsHandle& operator=(const PlateTectonicsHandle&) = delete;

    ~PlateTectonicsHandle()
    {
        if (pointer != nullptr)
            tectonics_destroy(pointer);
    }
};

float clamp01(float value)
{
    return std::clamp(value, 0.0f, 1.0f);
}

// Metres on the engine's datum, unrounded: the world stores elevation as float.
float nativeheightmetres(float value, uint16_t sealevelmeters)
{
    return TopographyCodec::internal_to_meters_exact(value, sealevelmeters);
}

// The hand-off of the engine's cube to the world's grid (W2-E, E5; sphere-mesh.md section 5),
// built once per solver state and cached: the identity on a cube world of the engine's N and
// radius; otherwise remapmean weights from the engine's cells to the world's for heights and
// rates (intensive), and per world cell the engine cell at its centre (remapnearest) for
// codes, ids, scores and the values whose pairs must stay one cell's. Vectors stay (east,
// north) at their cell on both grids, so a nearest copy moves them unchanged.
struct HandOff
{
    std::shared_ptr<const mesh::CellGrid> engine;   // the engine's cube
    std::shared_ptr<const mesh::CellGrid> world;    // the world's grid
    bool identity = false;
    mesh::RemapWeights weights;                     // engine to world (empty on the identity)
    std::vector<mesh::CellId> nearest;              // per world cell, the engine cell at its centre
    // Cell-counted lengths and velocities: engine spacings to the world's cells (on lat-lon the
    // equatorial cell 2 pi R / W, on a cube the mean cell side), kilometres per engine spacing,
    // and engine cells to world cells by true area (the mean cell areas).
    double cellscale = 1;
    double spacing_km = 0;
    double areascale = 1;

    HandOff(const planet& planet, std::shared_ptr<const mesh::CellGrid> enginegrid)
        : engine(std::move(enginegrid)), world(planet.sharedgrid())
    {
        const double radius_m = world->radius_m();
        identity = world->kind() == mesh::GridKind::cube && world->facecells() == engine->facecells() &&
            std::fabs(world->radius_m() - engine->radius_m()) <= 1e-9 * engine->radius_m();
        nearest.resize(world->cellcount());
        if (identity)
            for (mesh::CellId c = 0; c < world->cellcount(); c++)
                nearest[c] = c;
        else
        {
            weights = mesh::makeremap(*engine, *world);
            for (mesh::CellId c = 0; c < world->cellcount(); c++)
                nearest[c] = engine->locate(world->centre(c));
        }
        const double enginespacing_m = std::sqrt(engine->totalarea_m2() / static_cast<double>(engine->cellcount()));
        const double worldcell_m = world->kind() == mesh::GridKind::latlon
            ? 2.0 * 3.14159265358979323846 * radius_m / world->columns()
            : std::sqrt(world->totalarea_m2() / static_cast<double>(world->cellcount()));
        cellscale = identity ? 1.0 : enginespacing_m / worldcell_m;
        spacing_km = enginespacing_m / 1000.0;
        areascale = identity ? 1.0 : (engine->totalarea_m2() / static_cast<double>(engine->cellcount())) /
            (world->totalarea_m2() / static_cast<double>(world->cellcount()));
    }
    std::size_t cellcount() const { return world->cellcount(); }
    mesh::CellId cell(mesh::CellId c) const { return nearest[c]; }
    // An intensive field (heights, rates) on the world's cells.
    std::vector<float> mean(const float* field) const
    {
        std::vector<float> values(field, field + engine->cellcount());
        if (identity)
            return values;
        return mesh::remapmean(weights, *engine, *world, values);
    }
};

void applynativeheightmap(planet& world, worldgrid<bool>& shelves, const float* outputheightmap,
    const HandOff& handoff, uint16_t sealevelmeters)
{
    if (outputheightmap == nullptr)
        return;

    const int sealevel = sealevelmeters;
    const int maxelevation = int(TopographyCodec::kMaxHeightMeters) + 1;

    world.setsealevel(sealevel);
    world.setmaxelevation(maxelevation);

    // Metres per engine cell, then the area-weighted mean over each world cell.
    const std::size_t enginecells = handoff.engine->cellcount();
    std::vector<float> metres(enginecells);
    for (std::size_t c = 0; c < enginecells; c++)
        metres[c] = nativeheightmetres(outputheightmap[c], sealevelmeters);
    const std::vector<float> heights = handoff.mean(metres.data());
    for (mesh::CellId c = 0; c < handoff.cellcount(); c++)
    {
        world.setnomfloat(c, heights[c]);
        shelves[c] = false;
    }

    getlandandseatotals(world);
}

GeologicRegime translategeologicregime(std::uint8_t regime)
{
    using NativeRegime = tectonics::contract::GeologicRegime;

    switch (static_cast<NativeRegime>(regime))
    {
    case NativeRegime::Stable:
        return GeologicRegime::stable;
    case NativeRegime::ConvergentArc:
        return GeologicRegime::convergent_arc;
    case NativeRegime::ContinentCollision:
        return GeologicRegime::continent_collision;
    case NativeRegime::DivergentRift:
        return GeologicRegime::divergent_rift;
    case NativeRegime::Transform:
        return GeologicRegime::transform;
    case NativeRegime::PassiveMargin:
        return GeologicRegime::passive_margin;
    case NativeRegime::MidOceanRidge:
        return GeologicRegime::mid_ocean_ridge;
    case NativeRegime::TrenchAdjacent:
        return GeologicRegime::trench_adjacent;
    default:
        return GeologicRegime::stable;
    }
}

CrustClass translatecrustclass(std::uint8_t crustclass)
{
    using NativeCrustClass = tectonics::contract::CrustClass;

    switch (static_cast<NativeCrustClass>(crustclass))
    {
    case NativeCrustClass::Oceanic:
        return CrustClass::oceanic;
    case NativeCrustClass::Transitional:
        return CrustClass::transitional;
    case NativeCrustClass::Continental:
        return CrustClass::continental;
    case NativeCrustClass::None:
    default:
        return CrustClass::none;
    }
}

BoundaryType translateboundarytype(std::uint8_t boundarytype)
{
    using NativeBoundaryType = tectonics::contract::BoundaryType;

    switch (static_cast<NativeBoundaryType>(boundarytype))
    {
    case NativeBoundaryType::Convergent:
        return BoundaryType::convergent;
    case NativeBoundaryType::Divergent:
        return BoundaryType::divergent;
    case NativeBoundaryType::Transform:
        return BoundaryType::transform;
    case NativeBoundaryType::PassiveMargin:
        return BoundaryType::passive_margin;
    case NativeBoundaryType::None:
    default:
        return BoundaryType::none;
    }
}

DeformingRegionType translatedeformingregiontype(std::uint8_t regiontype)
{
    using NativeRegionType = tectonics::contract::DeformingRegionType;

    switch (static_cast<NativeRegionType>(regiontype))
    {
    case NativeRegionType::ContinentalRift:
        return DeformingRegionType::continental_rift;
    case NativeRegionType::DiffuseCollision:
        return DeformingRegionType::diffuse_collision;
    case NativeRegionType::None:
    default:
        return DeformingRegionType::none;
    }
}

float normalizeboundaryhistory(const tectonics::contract::BoundarySegment* segment)
{
    if (segment == nullptr)
        return 0.0f;

    const float agesignal = clamp01(static_cast<float>(segment->age_myr) / 40.0f);
    const float persistencesignal = clamp01(static_cast<float>(segment->persistence_steps) / 24.0f);
    return clamp01(0.55f * agesignal + 0.45f * persistencesignal);
}


}

struct GeologySolverState
{
    explicit GeologySolverState(const tectonics::scenario::Scenario& created) : scenario(created) {}
    GeologySolverState(const GeologySolverState&) = delete;
    GeologySolverState& operator=(const GeologySolverState&) = delete;

    tectonics::scenario::Scenario scenario;
    std::optional<HandOff> handoff;
    PlateTectonicsHandle simulation;
    mineralsystems::TectonicFrameReader mineralframes;
    std::optional<mineralsystems::MineralHistory> minerals;
    // What the state was made for; evolvegeosphere re-creates it on any mismatch.
    long seed = 0;
    float gravity = 0.0f;
    WorldTectonicsOptions options;
    std::uint32_t updates = 0;   // engine updates since the starting state
};

namespace
{
// The world's grid as the hand-off sees it: kind, size and radius.
bool samegrid(const mesh::CellGrid& a, const mesh::CellGrid& b)
{
    return a.kind() == b.kind() && a.facecells() == b.facecells() && a.columns() == b.columns() &&
        a.radius_m() == b.radius_m();
}

std::shared_ptr<GeologySolverState> createsolverstate(const planet& world, const WorldTectonicsOptions& options)
{
    auto scenario = options.scenario ? *options.scenario : tectonics::scenario::baseline();
    scenario.seed = world.seed();
    // The engine's cube on the world's sphere: on the world grid N is the cube world's own, or
    // W / 4 on a lat-lon world (an equatorial cell of the world's size); the body's radius.
    if (options.worldGrid)
    {
        const mesh::CellGrid& grid = world.grid();
        const int facecells = grid.kind() == mesh::GridKind::cube ? grid.facecells() : grid.columns() / 4;
        scenario.sphere_face_cells = static_cast<uint32_t>(std::max(2, facecells));
    }
    scenario.radius_km = world.radius() / 1000.0;
    if (options.durationMyr > 0) scenario.evolution_duration_myr = options.durationMyr;
    if (options.plateCount > 0) scenario.plate_count = static_cast<uint32_t>(options.plateCount);
    if (options.deltaTimeMyr > 0) scenario.delta_time_myr = options.deltaTimeMyr;
    scenario.gravity_mps2 = world.gravity() * tectonics::scenario::kDefaultGravityMps2;

    auto state = std::make_shared<GeologySolverState>(scenario);
    state->simulation.pointer = tectonics_create_from_scenario(&state->scenario);
    if (state->simulation.pointer == nullptr)
        return nullptr;
    state->seed = world.seed();
    state->gravity = world.gravity();
    state->options = options;
    state->handoff.emplace(world, tectonics_get_cell_grid(state->simulation.pointer));

    if (options.mineralHistory)
    {
        state->minerals.emplace(mineralsystems::MineralParameters{}, mineralsystems::mineralseed(world.seed()));
        state->minerals->begin(state->mineralframes.read(state->simulation.pointer));
    }
    return state;
}

bool madefor(const GeologySolverState& state, const planet& world, const WorldTectonicsOptions& options)
{
    return state.seed == world.seed() && state.gravity == world.gravity() &&
        samegrid(*state.handoff->world, world.grid()) && state.options.durationMyr == options.durationMyr &&
        state.options.plateCount == options.plateCount && state.options.deltaTimeMyr == options.deltaTimeMyr &&
        state.options.mineralHistory == options.mineralHistory && state.options.scenario == options.scenario &&
        state.options.worldGrid == options.worldGrid;
}

// The engine's current state as world state: elevation, the tectonic fields and
// tables, the engine's own fields and the snapshot record. Used for the starting
// Geosphere, at every checkpoint and at the end.
void writegeosphere(planet& world, worldgrid<bool>& shelves, GeologySolverState& state, const std::string& label)
{
    PlateTectonicsHandle& simulation = state.simulation;
    const HandOff& handoff = *state.handoff;
    const std::size_t enginecells = handoff.engine->cellcount();
    const auto cellscale = static_cast<float>(handoff.cellscale);
    const auto spacingkm = static_cast<float>(handoff.spacing_km);
    const double areascale = handoff.areascale;

    float* outputheightmap = tectonics_get_heightmap(simulation.pointer);

    if (outputheightmap == nullptr)
        return;

    applynativeheightmap(world, shelves, outputheightmap, handoff, tectonics_get_sea_level_m(simulation.pointer));

    const std::uint8_t* convergencemap = tectonics_get_convergence_map(simulation.pointer);
    const std::uint8_t* divergencemap = tectonics_get_divergence_map(simulation.pointer);
    const std::uint8_t* shearmap = tectonics_get_shear_map(simulation.pointer);
    const std::uint8_t* regimemap = tectonics_get_geologic_regime_map(simulation.pointer);
    const float* crustagemyrmap = tectonics_get_crust_age_myr_map(simulation.pointer);
    const float* crustthicknessmap = tectonics_get_crust_thickness_map(simulation.pointer);
    const std::uint8_t* crustclassmap = tectonics_get_crust_class_map(simulation.pointer);
    const float* uplifttendencymap = tectonics_get_uplift_tendency_map(simulation.pointer);
    const float* rockupliftratemap = tectonics_get_rock_uplift_rate_map(simulation.pointer);
    const float* shorteningmap = tectonics_get_shortening_map(simulation.pointer);
    const float* shorteningdirectionmap = tectonics_get_shortening_direction_map(simulation.pointer);
    const float* subsidencetendencymap = tectonics_get_subsidence_tendency_map(simulation.pointer);
    const float* accumulatedstrainmap = tectonics_get_accumulated_strain_map(simulation.pointer);
    const std::uint8_t* boundarytypemap = tectonics_get_boundary_type_map(simulation.pointer);
    const std::uint16_t* boundarydistancemap = tectonics_get_boundary_distance_map(simulation.pointer);
    const std::uint32_t* boundarysegmentidmap = tectonics_get_boundary_segment_id_map(simulation.pointer);
    const std::uint32_t* nearestboundaryidmap = tectonics_get_nearest_boundary_id_map(simulation.pointer);
    const std::uint32_t* deformingregionidmap = tectonics_get_deforming_region_id_map(simulation.pointer);
    const std::uint8_t* deformingregiontypemap = tectonics_get_deforming_region_type_map(simulation.pointer);
    const float* deformationratemap = tectonics_get_deformation_rate_map(simulation.pointer);
    const float* deformationvelocityxmap = tectonics_get_deformation_velocity_x_map(simulation.pointer);
    const float* deformationvelocityymap = tectonics_get_deformation_velocity_y_map(simulation.pointer);

    if (convergencemap == nullptr || divergencemap == nullptr || shearmap == nullptr || regimemap == nullptr
        || crustagemyrmap == nullptr || crustthicknessmap == nullptr || crustclassmap == nullptr
        || uplifttendencymap == nullptr || rockupliftratemap == nullptr || subsidencetendencymap == nullptr
        || shorteningmap == nullptr || shorteningdirectionmap == nullptr
        || accumulatedstrainmap == nullptr
        || boundarytypemap == nullptr || boundarydistancemap == nullptr || boundarysegmentidmap == nullptr
        || nearestboundaryidmap == nullptr || deformingregionidmap == nullptr || deformingregiontypemap == nullptr
        || deformationratemap == nullptr || deformationvelocityxmap == nullptr || deformationvelocityymap == nullptr)
        return;

    world.settectonictimeoriginstep(static_cast<int>(tectonics_get_time_origin_step(simulation.pointer)));
    world.settectonictimemyr(static_cast<float>(tectonics_get_time_myr(simulation.pointer)));
    world.settectonicdeltatimemyr(static_cast<float>(tectonics_get_delta_time_myr(simulation.pointer)));
    world.settectonicplatecount(static_cast<int>(tectonics_get_plate_count(simulation.pointer)));
    world.settectonicsealevelm(static_cast<int>(tectonics_get_sea_level_m(simulation.pointer)));

    // The tables: centroids stay longitude and latitude in degrees; lengths become km;
    // motions and velocities (east, north) world cells per update; counts world cells by area.
    const std::uint32_t boundarysegmentcount = tectonics_get_boundary_segment_count(simulation.pointer);
    const tectonics::contract::BoundarySegment* boundarysegments = tectonics_get_boundary_segments(simulation.pointer);
    std::vector<tectonics::contract::BoundarySegment> boundarysegmentobjects;
    if (boundarysegments != nullptr && boundarysegmentcount > 0)
        boundarysegmentobjects.assign(boundarysegments, boundarysegments + boundarysegmentcount);
    for (auto& segment : boundarysegmentobjects)
    {
        segment.length_cells *= spacingkm;   // km from here on (TectonicBoundarySegment::lengthKm)
        segment.average_normal_motion *= cellscale;
        segment.average_shear_motion *= cellscale;
        segment.cell_count = static_cast<std::uint32_t>(std::llround(segment.cell_count * areascale));
    }
    world.settectonicboundarysegments(std::move(boundarysegmentobjects));

    const std::uint32_t deformingregioncount = tectonics_get_deforming_region_count(simulation.pointer);
    const tectonics::contract::DeformingRegion* deformingregions = tectonics_get_deforming_regions(simulation.pointer);
    std::vector<tectonics::contract::DeformingRegion> deformingregionobjects;
    if (deformingregions != nullptr && deformingregioncount > 0)
        deformingregionobjects.assign(deformingregions, deformingregions + deformingregioncount);
    for (auto& region : deformingregionobjects)
    {
        region.average_interpolated_velocity_x *= cellscale;
        region.average_interpolated_velocity_y *= cellscale;
        region.average_normal_motion *= cellscale;
        region.average_shear_motion *= cellscale;
        region.cell_count = static_cast<std::uint32_t>(std::llround(region.cell_count * areascale));
    }
    world.settectonicdeformingregions(std::move(deformingregionobjects));

    std::vector<const tectonics::contract::BoundarySegment*> boundarysegmentlookup;
    if (boundarysegments != nullptr && boundarysegmentcount > 0)
    {
        std::size_t maxid = 0;
        for (std::uint32_t i = 0; i < boundarysegmentcount; i++)
            maxid = std::max(maxid, static_cast<std::size_t>(boundarysegments[i].id));

        boundarysegmentlookup.assign(maxid + 1, nullptr);
        for (std::uint32_t i = 0; i < boundarysegmentcount; i++)
            boundarysegmentlookup[boundarysegments[i].id] = &boundarysegments[i];
    }

    // Rates of elevation change: the area-weighted mean, as the heights (no length scale).
    const std::vector<float> rockupliftrate = handoff.mean(rockupliftratemap);
    for (mesh::CellId c = 0; c < handoff.cellcount(); c++)
    {
        const std::size_t index = handoff.cell(c);
        const int convergencescore = static_cast<int>(convergencemap[index]);
        const int divergencescore = static_cast<int>(divergencemap[index]);
        const int shearscore = static_cast<int>(shearmap[index]);
        const CrustClass crustclass = translatecrustclass(crustclassmap[index]);
        const BoundaryType boundarytype = translateboundarytype(boundarytypemap[index]);
        const DeformingRegionType regiontype = translatedeformingregiontype(deformingregiontypemap[index]);
        const std::uint32_t nearestboundaryid = nearestboundaryidmap[index];
        const tectonics::contract::BoundarySegment* boundarysegment =
            nearestboundaryid < boundarysegmentlookup.size() ? boundarysegmentlookup[nearestboundaryid] : nullptr;
        // This history score is app-derived for presentation/stylization, not native contract data.
        const float boundaryhistory = normalizeboundaryhistory(boundarysegment);

        world.settectonicconvergence(c, convergencescore);
        world.settectonicdivergence(c, divergencescore);
        world.settectonicshear(c, shearscore);
        world.setgeologicregime(c, translategeologicregime(regimemap[index]));
        world.settectoniccrustagemyr(c, crustagemyrmap[index]);
        world.settectoniccrustthickness(c, crustthicknessmap[index]);
        world.settectoniccrustclass(c, crustclass);
        world.settectonicuplifttendency(c, uplifttendencymap[index]);
        world.settectonicrockupliftrate(c, rockupliftrate[c]);
        // An axis wraps at pi, so the nearest engine cell, never a mean; the shortening with
        // it, so the pair stays one cell's. Both grids measure the axis from the cell's own east
        // toward its south, so the copy needs no rotation.
        world.settectonicshortening(c, shorteningmap[index]);
        world.settectonicshorteningdirection(c, shorteningdirectionmap[index]);
        world.settectonicsubsidencetendency(c, subsidencetendencymap[index]);
        world.settectonicaccumulatedstrain(c, accumulatedstrainmap[index]);
        world.settectonicboundarytype(c, boundarytype);
        // Engine spacings to world cells (the cell count the consumers' thresholds use).
        world.settectonicboundarydistance(c, static_cast<int>(std::lround(boundarydistancemap[index] * cellscale)));
        world.settectonicboundarysegmentid(c, static_cast<int>(boundarysegmentidmap[index]));
        world.settectonicnearestboundaryid(c, static_cast<int>(nearestboundaryid));
        world.settectonicboundaryhistory(c, boundaryhistory);
        world.settectonicdeformingregionid(c, static_cast<int>(deformingregionidmap[index]));
        world.settectonicdeformingregiontype(c, regiontype);
        world.settectonicdeformationrate(c, deformationratemap[index]);
        // (east, north) at the cell on both grids: only the length unit changes.
        world.settectonicdeformationvelocityx(c, deformationvelocityxmap[index] * cellscale);
        world.settectonicdeformationvelocityy(c, deformationvelocityymap[index] * cellscale);
    }

    const tectonics::contract::Snapshot snapshot = tectonics_capture_snapshot(simulation.pointer);
    // Column and plume fields exist only when the engine keeps them.
    const auto cellvalue = [enginecells](const std::vector<float>& values, std::size_t index)
    {
        return values.size() == enginecells ? values[index] : 0.0f;
    };
    for (mesh::CellId c = 0; c < handoff.cellcount(); c++)
    {
        const std::size_t index = handoff.cell(c);
        world.settectonicsnapshotcell(c, snapshot.heightmap[index], static_cast<int>(snapshot.plate_id[index]),
            static_cast<int>(snapshot.crust_age_steps[index]), cellvalue(snapshot.continental_crust_m, index),
            cellvalue(snapshot.continental_area, index), cellvalue(snapshot.oceanic_crust_m, index), cellvalue(snapshot.sediment_m, index),
            cellvalue(snapshot.rift_damage, index), cellvalue(snapshot.plume_influence, index));
    }

    tectonics::tooling::SnapshotBundleOptions record;
    record.label = label;
    record.source = "world_pipeline";
    record.include_update_count = true;
    record.update_count = state.updates;
    world.settectonicsnapshotrecord(tectonics::tooling::snapshot_record_json(snapshot, record));
}

// A checkpoint time as a count of updates; it must fall on one.
std::uint32_t checkpointupdate(double timemyr, double deltatimemyr)
{
    const double updates = timemyr / deltatimemyr;
    const double whole = std::round(updates);
    if (!(timemyr > 0.0) || std::fabs(updates - whole) > 1e-6)
        throw std::runtime_error("Checkpoint time " + std::to_string(timemyr) +
            " Myr is not a positive multiple of the tectonic time step (" + std::to_string(deltatimemyr) + " Myr)");
    return static_cast<std::uint32_t>(whole);
}

std::string updatelabel(std::uint32_t updates)
{
    std::string digits = std::to_string(updates);
    return "update" + std::string(digits.size() < 6 ? 6 - digits.size() : 0, '0') + digits;
}
}

std::shared_ptr<GeologySolverState> initializegeosphere(planet& world, worldgrid<bool>& shelves,
                                                        const WorldTectonicsOptions& options)
{
    world.cleartectonicprovenance();
    auto state = createsolverstate(world, options);
    if (state == nullptr)
        return nullptr;

    // The starting Geosphere, for review and as the first snapshot; tectonics rewrites it.
    writegeosphere(world, shelves, *state, "initial");
    return state;
}

void generateworldtectonics(planet& world, worldgrid<bool>& shelves, const WorldTectonicsOptions& options)
{
    evolvegeosphere(world, shelves, options, initializegeosphere(world, shelves, options));
}

void evolvegeosphere(planet& world, worldgrid<bool>& shelves, const WorldTectonicsOptions& options,
                     std::shared_ptr<GeologySolverState> state)
{
    if (state == nullptr || !madefor(*state, world, options))
    {
        world.cleartectonicprovenance();
        state = createsolverstate(world, options);
        if (state == nullptr)
            return;
    }
    PlateTectonicsHandle& simulation = state->simulation;
    mineralsystems::TectonicFrameReader& mineralframes = state->mineralframes;
    std::optional<mineralsystems::MineralHistory>& minerals = state->minerals;

    const double deltatimemyr = tectonics_get_delta_time_myr(simulation.pointer);
    const std::uint32_t checkpointevery = options.checkpointEveryMyr > 0.0 ?
        checkpointupdate(options.checkpointEveryMyr, deltatimemyr) : 0;
    std::vector<std::uint32_t> checkpoints;
    for (const double timemyr : options.checkpointsMyr)
        checkpoints.push_back(checkpointupdate(timemyr, deltatimemyr));

    const bool visualizeprogress = hasworldgenvisualizationcallback();
    auto nextvisualization = std::chrono::steady_clock::now();

    while (tectonics_is_finished(simulation.pointer) == 0)
    {
        tectonics_step(simulation.pointer);
        state->updates++;
        if (minerals)
            minerals->step(mineralframes.read(simulation.pointer));

        // The end of the run is the step's result, not a checkpoint.
        if (options.checkpoint && tectonics_is_finished(simulation.pointer) == 0 &&
            ((checkpointevery != 0 && state->updates % checkpointevery == 0) ||
                std::find(checkpoints.begin(), checkpoints.end(), state->updates) != checkpoints.end()))
        {
            writegeosphere(world, shelves, *state, updatelabel(state->updates));
            options.checkpoint(tectonics_get_time_myr(simulation.pointer));
        }

        if (visualizeprogress == false || std::chrono::steady_clock::now() < nextvisualization)
            continue;

        float* previewheightmap = tectonics_get_heightmap(simulation.pointer);

        if (previewheightmap != nullptr)
        {
            applynativeheightmap(world, shelves, previewheightmap, *state->handoff,
                tectonics_get_sea_level_m(simulation.pointer));
            requestworldgenvisualization();
        }

        nextvisualization = std::chrono::steady_clock::now() + kVisualizationInterval;
    }

    writegeosphere(world, shelves, *state, "final");

    if (visualizeprogress)
        requestworldgenvisualization();

    if (minerals)
    {
        minerals->finish(mineralframes.read(simulation.pointer));
        world.setmineralbudget(minerals->budget());
        world.setmineraldeposits(minerals->takedeposits());
    }
}
