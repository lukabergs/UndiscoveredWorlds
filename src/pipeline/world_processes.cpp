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
          "Continues the initializer's engine; re-creates it deterministically when the run carries none." },
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
        { S::fjords, "terrain.fjords", K::detailer, false, {}, {}, {}, {}, "" },
        { S::rivers_and_basins, "hydrology.rivers_and_basins", K::diagnoser, false, {}, {}, {}, {},
          "Also places salt lakes, adds land noise and fills depressions (terrain edits)." },
        { S::basin_editor, "hydrology.basin_policy", K::editor, false, {}, {}, {}, {}, "" },
        { S::lakes, "hydrology.lakes", K::detailer, false, {}, {}, {}, {}, "" },
        { S::mountain_temperature_lapse, "climate.mountain_lapse", K::diagnoser, false, {}, {}, {}, {}, "" },
        { S::climates_and_biomes, "climate.classification", K::diagnoser, false, {}, {}, {}, {}, "Writes biomes as well as climate classes." },
        { S::arid_features, "terrain.arid_features", K::detailer, false, {}, {}, {}, {}, "" },
        { S::deltas_wetlands_and_roughness, "hydrology.deltas_and_wetlands", K::detailer, false, {}, {}, {}, {},
          "Also recomputes climate classes and biomes and writes the climate diagnostics." },
        { S::finalize_layers, "layers.finalize", K::diagnoser, true, {}, {}, {}, {},
          "Physical resource layers, gold placers and optional society." },
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
