#include "pipeline/world_processes.hpp"

#include <algorithm>
#include <stdexcept>

#include "core/deterministic_random.hpp"
#include "core/procedural_noise.hpp"

#include "pipeline/world_process_access.inc"

namespace
{
using K = ProcessKind;
using S = GenerationStageId;

// Kinds and notes describe the current stages. Access lists and the stochastic
// flag come from world_process_access.inc (measured, see world_processes.hpp);
// the flags set here cover generators the analysis cannot see (the tectonics
// engine's own stream, society and mineral streams).
std::vector<WorldProcessDefinition> buildprocesses()
{
    std::vector<WorldProcessDefinition> processes =
    {
        { S::geosphere_initialization, "geology.geosphere_initialization", K::initializer, true, {}, {}, {}, {},
          "Creates the engine (starting topography from its own seeded stream) and the inherited gold provinces." },
        { S::plate_tectonics, "geology.plate_tectonics", K::evolver, true, {}, {}, {}, {},
          "Continues the initializer's engine; re-creates it deterministically when the run carries none.",
          { nullptr, nullptr, "configs/tectonics-baseline.json (evolution_duration_myr, delta_time_myr); the run sets them "
                              "with --duration and --delta-time-myr; checkpoints fall on whole steps" } },
        { S::terraforming, "terrain.terraforming", K::editor, false, {}, {}, {}, {}, "" },
        { S::inland_seas, "terrain.inland_sea_policy", K::editor, false, {}, {}, {}, {}, "" },
        { S::archipelagos, "terrain.archipelagos", K::detailer, false, {}, {}, {}, {}, "" },
        { S::ocean_cleanup, "terrain.ocean_cleanup", K::editor, false, {}, {}, {}, {}, "" },
        { S::coastline_refinement, "terrain.coastline_refinement", K::detailer, false, {}, {}, {}, {}, "" },
        { S::ocean_depth_refinement, "terrain.ocean_depth_refinement", K::detailer, false, {}, {}, {}, {}, "" },
        { S::tectonic_trenches, "terrain.tectonic_trenches", K::detailer, false, {}, {}, {}, {}, "" },
        { S::tectonic_volcanoes, "terrain.tectonic_volcanoes", K::detailer, false, {}, {}, {}, {}, "" },
        { S::terrain_texturing, "terrain.texturing", K::diagnoser, false, {}, {}, {}, {}, "" },
        { S::mountain_bases, "terrain.mountain_bases", K::detailer, false, {}, {}, {}, {}, "" },
        { S::terrain_smoothing, "terrain.smoothing", K::editor, false, {}, {}, {}, {}, "" },
        { S::canyon_uplift, "terrain.canyon_uplift", K::detailer, false, {}, {}, {}, {}, "" },
        { S::depression_fill, "terrain.depression_fill", K::editor, false, {}, {}, {}, {},
          "Adds land noise between two depression fills." },
        { S::coastline_adjustment, "terrain.coastline_adjustment", K::editor, false, {}, {}, {}, {}, "" },
        { S::island_check, "terrain.island_check", K::editor, false, {}, {}, {}, {}, "" },
        { S::terrain_roughness, "terrain.roughness", K::detailer, false, {}, {}, {}, {},
          "Pure seeded noise; climate does not read it." },
        { S::global_temperature, "climate.temperature", K::evolver, false, {}, {}, {}, {},
          "Adds seeded temperature noise (energybalance::fractalTemperatureRangeC)." },
        { S::pressure, "climate.pressure", K::evolver, false, {}, {}, {}, {}, "First guess; the coupled loop in climate.rainfall recomputes it." },
        { S::winds, "climate.winds", K::evolver, false, {}, {}, {}, {}, "First guess; the coupled loop in climate.rainfall recomputes it." },
        { S::ocean_currents, "climate.ocean_currents", K::evolver, false, {}, {}, {}, {}, "First guess; the coupled loop in climate.rainfall recomputes it." },
        { S::sea_surface_temperatures, "climate.sea_surface_temperature", K::evolver, false, {}, {}, {}, {}, "First guess; the coupled loop in climate.rainfall recomputes it." },
        { S::sea_ice_and_tides, "climate.sea_ice_and_tides", K::evolver, false, {}, {}, {}, {}, "Tides are hydrosphere state." },
        { S::rainfall, "climate.rainfall_and_coupling", K::evolver, false, {}, {}, {}, {},
          "Runs the coupled wind/pressure/ocean/SST/rainfall iteration." },
        { S::flow_routing, "hydrology.flow_routing", K::diagnoser, false, {}, {}, {}, {},
          "Runoff from Fu's Budyko curve of the annual precipitation and the open-water evaporation law of the lakes; discharge loses evaporation from channel and floodplain where PET exceeds P (erosion plan C2); every depression overflows at its spill." },
        { S::fluvial_erosion, "terrain.fluvial_erosion", K::evolver, false, {}, {}, {}, {},
          "Writes its final routing fields too, so no second routing step follows it. A checkpoint writes the "
          "same fields as of its time (the routing of the terrain then, rates over the elapsed time). Glaciers "
          "from an equilibrium-line altitude of the coarse climate erode with the rivers; during the run the "
          "routing takes its sea from their ice surface (grounded ice stays land); it writes the ELA, the glaciation, the glacial denudation and, at the end, the "
          "present-day ice thickness on the final terrain (erosion plan Phase G).",
          { "landscape.erosion.duration_myr", "landscape.erosion.step_myr",
            "configs/parameters/landscape.json; the last step is shortened to end at the duration; a checkpoint "
            "follows the step that reaches its time" } },
        { S::rivers, "hydrology.rivers", K::diagnoser, false, {}, {}, {}, {},
          "River directions and the January and July flows in legacy units from the routing (erosion plan B1); in an eroded world also the distributaries on prograded land, in the delta fields, whatever the deltas option (erosion plan C5b)." },
        { S::lake_levels, "hydrology.lake_levels", K::diagnoser, false, {}, {}, {}, {},
          "Lake levels by water balance in the routing's depressions; salt lakes and pans where closed (erosion plan B2)." },
        { S::mountain_temperature_lapse, "climate.mountain_lapse", K::diagnoser, false, {}, {}, {}, {}, "" },
        { S::climates_and_biomes, "climate.classification", K::diagnoser, false, {}, {}, {}, {}, "Writes biomes as well as climate classes." },
        { S::arid_features, "terrain.arid_features", K::detailer, false, {}, {}, {}, {}, "" },
        { S::deltas_wetlands_and_roughness, "hydrology.deltas_and_wetlands", K::detailer, false, {}, {}, {}, {},
          "Also recomputes climate classes and biomes and writes the climate diagnostics. In an eroded world it clears the distributaries of hydrology.rivers when deltas are turned off (erosion plan C5b)." },
        { S::subgrid_statistics, "terrain.subgrid_statistics", K::diagnoser, false, {}, {}, {}, {},
          "The sub-grid statistics the regional levels reproduce (multiscale plan S3.6): the relief amplitude "
          "and Perron's valley spacing from the erosion run's stream-power and hillslope laws, on the final terrain." },
        { S::finalize_layers, "layers.finalize", K::diagnoser, true, {}, {}, {}, {},
          "Physical resource layers, gold placers and optional society. In an eroded world it applies the erosion run's exhumation to the deposits, so it runs once per world." },
    };

    for (WorldProcessDefinition& process : processes)
        applymeasuredaccess(process);
    return processes;
}
}

const char* processkindname(ProcessKind kind)
{
    switch (kind)
    {
    case ProcessKind::initializer: return "initializer";
    case ProcessKind::evolver: return "evolver";
    case ProcessKind::diagnoser: return "diagnoser";
    case ProcessKind::detailer: return "detailer";
    case ProcessKind::editor: return "editor";
    }
    return "";
}

const std::vector<WorldProcessDefinition>& getworldprocesses()
{
    static const std::vector<WorldProcessDefinition> processes = buildprocesses();
    return processes;
}

const WorldProcessDefinition& getworldprocess(GenerationStageId stage)
{
    for (const WorldProcessDefinition& process : getworldprocesses())
        if (process.stage == stage)
            return process;
    throw std::out_of_range("Generation stage has no process definition");
}

const WorldProcessDefinition* findworldprocess(std::string_view id)
{
    for (const WorldProcessDefinition& process : getworldprocesses())
        if (id == process.id)
            return &process;
    return nullptr;
}

std::string timedevolverids()
{
    std::string ids;
    for (const WorldProcessDefinition& process : getworldprocesses())
        if (process.time.source != nullptr)
            ids += (ids.empty() ? "" : ", ") + std::string(process.id);
    return ids;
}

bool declareswrite(const WorldProcessDefinition& process, worldmodel::FieldId field)
{
    return std::find(process.writes.begin(), process.writes.end(), field) != process.writes.end();
}

bool declareswrite(const WorldProcessDefinition& process, worldmodel::TableId table)
{
    return std::find(process.tablewrites.begin(), process.tablewrites.end(), table) != process.tablewrites.end();
}

void seedprocessstream(long worldseed, GenerationStageId stage)
{
    // FNV-1a of the process id: stable across builds and stage reordering.
    std::uint32_t hash = 2166136261u;
    for (const char* character = getworldprocess(stage).id; *character != '\0'; character++)
    {
        hash ^= static_cast<unsigned char>(*character);
        hash *= 16777619u;
    }
    fast_srand(deterministicfastseed(deterministiccontextseed(worldseed, static_cast<int>(hash))));
}
