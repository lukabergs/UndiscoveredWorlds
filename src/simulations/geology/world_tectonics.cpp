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
#include "simulations/geology/plate_tectonics/tectonic_world_grid.hpp"
#include "simulations/geology/plate_tectonics/topography_codec.hpp"
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

void applynativeheightmap(planet& world, worldgrid<bool>& shelves, const float* outputheightmap, std::size_t cellcount, int simwidth, int simheight, uint16_t sealevelmeters)
{
    if (outputheightmap == nullptr || cellcount == 0)
        return;

    const int width = world.width();
    const int height = world.height();
    const tectonics::WorldGrid grid(simwidth, simheight, width + 1, height + 1);
    const int sealevel = sealevelmeters;
    const int maxelevation = int(TopographyCodec::kMaxHeightMeters) + 1;

    world.setsealevel(sealevel);
    world.setmaxelevation(maxelevation);

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int y = startrow; y <= endrow; y++)
        {
            for (int x = 0; x <= width; x++)
            {
                world.setnomfloat(x, y, nativeheightmetres(grid.height(outputheightmap, x, y), sealevelmeters));
                shelves[x][y] = false;
            }
        }
    });

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
    explicit GeologySolverState(const tectonics::scenario::Scenario& created)
        : scenario(created), mineralframes(created.world_width_km) {}
    GeologySolverState(const GeologySolverState&) = delete;
    GeologySolverState& operator=(const GeologySolverState&) = delete;

    tectonics::scenario::Scenario scenario;
    PlateTectonicsHandle simulation;
    mineralsystems::TectonicFrameReader mineralframes;
    std::optional<mineralsystems::MineralHistory> minerals;
    // What the state was made for; evolvegeosphere re-creates it on any mismatch.
    long seed = 0;
    float gravity = 0.0f;
    int width = 0;
    int height = 0;
    WorldTectonicsOptions options;
    std::uint32_t updates = 0;   // engine updates since the starting state
};

namespace
{
std::shared_ptr<GeologySolverState> createsolverstate(const planet& world, const WorldTectonicsOptions& options)
{
    auto scenario = options.scenario ? *options.scenario : tectonics::scenario::baseline();
    scenario.seed = world.seed();
    if (options.worldGrid)
    {
        scenario.width = static_cast<uint32_t>(world.width() + 1);
        scenario.height = static_cast<uint32_t>(world.height() + 1);
    }
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
    state->width = world.width();
    state->height = world.height();
    state->options = options;

    if (options.mineralHistory)
    {
        state->minerals.emplace(mineralsystems::MineralParameters{}, mineralsystems::mineralseed(world.seed()));
        state->minerals->begin(state->mineralframes.read(state->simulation.pointer));
    }
    return state;
}

bool madefor(const GeologySolverState& state, const planet& world, const WorldTectonicsOptions& options)
{
    return state.seed == world.seed() && state.gravity == world.gravity() && state.width == world.width() &&
        state.height == world.height() && state.options.durationMyr == options.durationMyr &&
        state.options.plateCount == options.plateCount && state.options.deltaTimeMyr == options.deltaTimeMyr &&
        state.options.mineralHistory == options.mineralHistory && state.options.scenario == options.scenario &&
        state.options.worldGrid == options.worldGrid;
}

// The engine's current state as world state: elevation, the tectonic fields and
// tables, the engine's own fields and the snapshot record. Used for the starting
// Geosphere, at every checkpoint and at the end.
void writegeosphere(planet& world, worldgrid<bool>& shelves, GeologySolverState& state, const std::string& label)
{
    const auto& scenario = state.scenario;
    PlateTectonicsHandle& simulation = state.simulation;
    const int width = world.width(), height = world.height();
    const int simwidth = static_cast<int>(scenario.width), simheight = static_cast<int>(scenario.height);
    const std::size_t cellcount = static_cast<std::size_t>(simwidth) * simheight;
    const tectonics::WorldGrid grid(simwidth, simheight, width + 1, height + 1);
    const auto scale = static_cast<float>(grid.scale());

    float* outputheightmap = tectonics_get_heightmap(simulation.pointer);

    if (outputheightmap == nullptr)
        return;

    applynativeheightmap(world, shelves, outputheightmap, cellcount, simwidth, simheight,
        tectonics_get_sea_level_m(simulation.pointer));

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

    const std::uint32_t boundarysegmentcount = tectonics_get_boundary_segment_count(simulation.pointer);
    const tectonics::contract::BoundarySegment* boundarysegments = tectonics_get_boundary_segments(simulation.pointer);
    std::vector<tectonics::contract::BoundarySegment> boundarysegmentobjects;
    if (boundarysegments != nullptr && boundarysegmentcount > 0)
        boundarysegmentobjects.assign(boundarysegments, boundarysegments + boundarysegmentcount);
    for (auto& segment : boundarysegmentobjects)
    {
        segment.centroid_x = grid.worldCoordinate(segment.centroid_x);
        segment.centroid_y = grid.worldCoordinate(segment.centroid_y);
        segment.length_cells *= scale;
        segment.average_normal_motion *= scale;
        segment.average_shear_motion *= scale;
        segment.cell_count = static_cast<std::uint32_t>(std::lround(segment.cell_count * scale * scale));
    }
    world.settectonicboundarysegments(std::move(boundarysegmentobjects));

    const std::uint32_t deformingregioncount = tectonics_get_deforming_region_count(simulation.pointer);
    const tectonics::contract::DeformingRegion* deformingregions = tectonics_get_deforming_regions(simulation.pointer);
    std::vector<tectonics::contract::DeformingRegion> deformingregionobjects;
    if (deformingregions != nullptr && deformingregioncount > 0)
        deformingregionobjects.assign(deformingregions, deformingregions + deformingregioncount);
    for (auto& region : deformingregionobjects)
    {
        region.centroid_x = grid.worldCoordinate(region.centroid_x);
        region.centroid_y = grid.worldCoordinate(region.centroid_y);
        region.average_interpolated_velocity_x *= scale;
        region.average_interpolated_velocity_y *= scale;
        region.average_normal_motion *= scale;
        region.average_shear_motion *= scale;
        region.cell_count = static_cast<std::uint32_t>(std::lround(region.cell_count * scale * scale));
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

    for (int y = 0; y <= height; y++)
    {
        for (int x = 0; x <= width; x++)
        {
            const std::size_t index = grid.nativeIndex(x, y);
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

            world.settectonicconvergence(x, y, convergencescore);
            world.settectonicdivergence(x, y, divergencescore);
            world.settectonicshear(x, y, shearscore);
            world.setgeologicregime(x, y, translategeologicregime(regimemap[index]));
            world.settectoniccrustagemyr(x, y, crustagemyrmap[index]);
            world.settectoniccrustthickness(x, y, crustthicknessmap[index]);
            world.settectoniccrustclass(x, y, crustclass);
            world.settectonicuplifttendency(x, y, uplifttendencymap[index]);
            // A rate of elevation change: sampled bilinearly as the heights are, no length scale.
            world.settectonicrockupliftrate(x, y, grid.height(rockupliftratemap, x, y));
            // An axis wraps at pi, so the nearest engine cell, never a bilinear mean; the
            // shortening with it, so the pair stays one cell's.
            world.settectonicshortening(x, y, shorteningmap[index]);
            world.settectonicshorteningdirection(x, y, shorteningdirectionmap[index]);
            world.settectonicsubsidencetendency(x, y, subsidencetendencymap[index]);
            world.settectonicaccumulatedstrain(x, y, accumulatedstrainmap[index]);
            world.settectonicboundarytype(x, y, boundarytype);
            world.settectonicboundarydistance(x, y, static_cast<int>(std::lround(boundarydistancemap[index] * scale)));
            world.settectonicboundarysegmentid(x, y, static_cast<int>(boundarysegmentidmap[index]));
            world.settectonicnearestboundaryid(x, y, static_cast<int>(nearestboundaryid));
            world.settectonicboundaryhistory(x, y, boundaryhistory);
            world.settectonicdeformingregionid(x, y, static_cast<int>(deformingregionidmap[index]));
            world.settectonicdeformingregiontype(x, y, regiontype);
            world.settectonicdeformationrate(x, y, deformationratemap[index]);
            world.settectonicdeformationvelocityx(x, y, deformationvelocityxmap[index] * scale);
            world.settectonicdeformationvelocityy(x, y, deformationvelocityymap[index] * scale);
        }
    }

    const tectonics::contract::Snapshot snapshot = tectonics_capture_snapshot(simulation.pointer);
    // Column and plume fields exist only when the engine keeps them.
    const auto cellvalue = [cellcount](const std::vector<float>& values, std::size_t index)
    {
        return values.size() == cellcount ? values[index] : 0.0f;
    };
    for (int y = 0; y <= height; y++)
    {
        for (int x = 0; x <= width; x++)
        {
            const std::size_t index = grid.nativeIndex(x, y);
            world.settectonicsnapshotcell(x, y, snapshot.heightmap[index], static_cast<int>(snapshot.plate_id[index]),
                static_cast<int>(snapshot.crust_age_steps[index]), cellvalue(snapshot.continental_crust_m, index),
                cellvalue(snapshot.oceanic_crust_m, index), cellvalue(snapshot.sediment_m, index),
                cellvalue(snapshot.rift_damage, index), cellvalue(snapshot.plume_influence, index));
        }
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
    const auto& scenario = state->scenario;
    PlateTectonicsHandle& simulation = state->simulation;
    mineralsystems::TectonicFrameReader& mineralframes = state->mineralframes;
    std::optional<mineralsystems::MineralHistory>& minerals = state->minerals;
    const int simwidth = static_cast<int>(scenario.width), simheight = static_cast<int>(scenario.height);
    const std::size_t cellcount = static_cast<std::size_t>(simwidth) * simheight;

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
            applynativeheightmap(world, shelves, previewheightmap, cellcount, simwidth, simheight,
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
