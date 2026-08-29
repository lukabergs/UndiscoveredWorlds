#include <stdexcept>

#include "generation_stage_registry.hpp"

using namespace std;

const vector<GenerationStageDefinition>& getgenerationstages()
{
    static const vector<GenerationStageDefinition> generationstages =
    {
        { GenerationStageId::plate_tectonics, "Plate tectonics", "Generate the starting global heightmap from the tectonic simulation seed and parameters.", GenerationStageDomain::terrain, true, false },
        { GenerationStageId::terraforming, "Terraforming", "Adjust sea level and remap the world height distribution to reshape continents and ocean basins.", GenerationStageDomain::terrain, true, false },
        { GenerationStageId::inland_seas, "Inland seas", "Review connected inland seas and choose which ones to keep, fill, or drain.", GenerationStageDomain::terrain, true, true },
        { GenerationStageId::archipelagos, "Archipelagos", "Optionally add fragmented coastal islands and archipelago chains.", GenerationStageDomain::terrain, true, true },
        { GenerationStageId::ocean_cleanup, "Ocean cleanup", "Optionally tidy narrow channels and simplify open-ocean structure.", GenerationStageDomain::terrain, true, true },
        { GenerationStageId::coastline_refinement, "Coastline refinement", "Optionally break up straight shorelines and add more natural coastal variation.", GenerationStageDomain::terrain, true, true },
        { GenerationStageId::ocean_depth_refinement, "Ocean depth refinement", "Optionally deepen shelves and open-ocean basins after the coastline passes.", GenerationStageDomain::terrain, true, true },
        { GenerationStageId::tectonic_trenches, "Tectonic trenches", "Optionally imprint convergent-margin trench systems from tectonic metadata.", GenerationStageDomain::terrain, true, true },
        { GenerationStageId::tectonic_volcanoes, "Tectonic volcanoes", "Optionally place volcanic arcs and hotspot-style uplifts from tectonic signals.", GenerationStageDomain::terrain, true, true },
        { GenerationStageId::terrain_texturing, "Terrain texturing", "Refresh derived terrain masks and shelves after major height edits.", GenerationStageDomain::terrain, true, true },
        { GenerationStageId::fastlem_mountains, "FastLEM mountains", "Optionally generate major mountain systems with FastLEM.", GenerationStageDomain::terrain, true, true },
        { GenerationStageId::mountain_bases, "Mountain base uplift", "Optionally broaden mountain roots and connect ranges to surrounding uplands.", GenerationStageDomain::terrain, true, true },
        { GenerationStageId::terrain_smoothing, "Terrain smoothing", "Optionally smooth terrain while preserving the large-scale landform structure.", GenerationStageDomain::terrain, true, true },
        { GenerationStageId::canyon_uplift, "Canyon uplift", "Optionally raise land around incised drainage corridors and canyon systems.", GenerationStageDomain::terrain, true, true },
        { GenerationStageId::depression_fill, "Depression fill", "Optionally fill closed depressions and reduce unwanted sinks.", GenerationStageDomain::terrain, true, true },
        { GenerationStageId::coastline_adjustment, "Coastline adjustment", "Optionally normalize coastal elevations and shallow seas.", GenerationStageDomain::terrain, true, true },
        { GenerationStageId::island_check, "Island check", "Optionally remove tiny isolated artifacts and stabilize island topology.", GenerationStageDomain::terrain, true, true },
        { GenerationStageId::terrain_roughness, "Terrain roughness", "Build the roughness field used by later terrain and climate layers.", GenerationStageDomain::terrain, true, false },
        { GenerationStageId::global_temperature, "Global temperature", "Generate the baseline temperature field from latitude, elevation, and orbital settings.", GenerationStageDomain::climate, true, false },
        { GenerationStageId::ocean_currents, "Ocean currents", "Optionally derive ocean-current structure for later climate passes.", GenerationStageDomain::climate, true, true },
        { GenerationStageId::sea_surface_temperatures, "Sea surface temperatures", "Optionally propagate sea-surface temperatures from the ocean state.", GenerationStageDomain::climate, true, true },
        { GenerationStageId::pressure, "Pressure", "Optionally generate the global pressure field.", GenerationStageDomain::climate, true, true },
        { GenerationStageId::winds, "Winds", "Optionally derive prevailing winds from pressure and planetary settings.", GenerationStageDomain::climate, true, true },
        { GenerationStageId::sea_ice_and_tides, "Sea ice and tides", "Optionally add sea-ice coverage and tidal effects to the climate model.", GenerationStageDomain::climate, true, true },
        { GenerationStageId::rainfall, "Rainfall", "Generate rainfall and moisture transport across the current terrain.", GenerationStageDomain::hydrology, true, false },
        { GenerationStageId::fjords, "Fjords", "Optionally carve fjords along cold coastal mountain fronts.", GenerationStageDomain::hydrology, true, true },
        { GenerationStageId::rivers_and_basins, "Rivers and basins", "Plan rivers and drainage basins from the current terrain and rainfall.", GenerationStageDomain::hydrology, true, true },
        { GenerationStageId::basin_editor, "Basin editor", "Review detected endorheic basins and choose whether to keep, fill, or drain them.", GenerationStageDomain::hydrology, true, true },
        { GenerationStageId::lakes, "Lakes", "Optionally generate lakes and rift lakes from the current drainage network.", GenerationStageDomain::hydrology, true, true },
        { GenerationStageId::post_river_fastlem, "Post-river FastLEM", "Optionally broaden FastLEM terrain around established river systems.", GenerationStageDomain::hydrology, true, true },
        { GenerationStageId::mountain_temperature_lapse, "Mountain temperature lapse", "Optionally apply mountain lapse-rate cooling to seasonal temperatures.", GenerationStageDomain::climate, true, true },
        { GenerationStageId::climates_and_biomes, "Climates and biomes", "Generate climate classes and biome distributions.", GenerationStageDomain::climate, true, false },
        { GenerationStageId::arid_features, "Arid features", "Optionally add dunes, salt pans, and other arid-land features.", GenerationStageDomain::hydrology, true, true },
        { GenerationStageId::deltas_wetlands_and_roughness, "Deltas, wetlands, and roughness", "Optionally add deltas, wetlands, and final roughness refinements.", GenerationStageDomain::hydrology, true, true },
        { GenerationStageId::finalize_layers, "Finalize layers", "Build the remaining physical and social layers for the finished world.", GenerationStageDomain::finalization, false, false }
    };

    return generationstages;
}

const GenerationStageDefinition& getgenerationstage(GenerationStageId id)
{
    for (const GenerationStageDefinition& stage : getgenerationstages())
    {
        if (stage.id == id)
            return stage;
    }

    throw out_of_range("Unknown generation stage");
}

size_t getgenerationstageindex(GenerationStageId id)
{
    const vector<GenerationStageDefinition>& stages = getgenerationstages();

    for (size_t index = 0; index < stages.size(); index++)
    {
        if (stages[index].id == id)
            return index;
    }

    throw out_of_range("Unknown generation stage");
}

bool isterraingenerationstage(GenerationStageId id)
{
    return getgenerationstage(id).domain == GenerationStageDomain::terrain;
}

bool isclimategenerationstage(GenerationStageId id)
{
    const GenerationStageDomain domain = getgenerationstage(id).domain;
    return domain == GenerationStageDomain::climate || domain == GenerationStageDomain::hydrology;
}
