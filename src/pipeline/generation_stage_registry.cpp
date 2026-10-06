#include <stdexcept>

#include "pipeline/generation_stage_registry.hpp"

using namespace std;

bool isclimatefeedbackstage(GenerationStageId id)
{
    switch (id) {
    case GenerationStageId::terrain_roughness:
    case GenerationStageId::global_temperature:
    case GenerationStageId::pressure:
    case GenerationStageId::winds:
    case GenerationStageId::ocean_currents:
    case GenerationStageId::sea_surface_temperatures:
    case GenerationStageId::sea_ice_and_tides:
    case GenerationStageId::rainfall:
    case GenerationStageId::mountain_temperature_lapse:
    case GenerationStageId::climates_and_biomes: return true;
    default: return false;
    }
}

const vector<GenerationStageDefinition>& getgenerationstages()
{
    static const vector<GenerationStageDefinition> generationstages =
    {
        { GenerationStageId::geosphere_initialization, "Initial geosphere", "Create the starting crust, plates and elevation and the inherited deposit provinces from the tectonic seed and parameters.", GenerationStageDomain::terrain, true, false },
        { GenerationStageId::plate_tectonics, "Plate tectonics", "Evolve the starting geosphere through the tectonic simulation to the global heightmap.", GenerationStageDomain::terrain, true, false },
        { GenerationStageId::terraforming, "Terraforming", "Adjust sea level and remap the world height distribution to reshape continents and ocean basins.", GenerationStageDomain::terrain, true, false },
        { GenerationStageId::inland_seas, "Inland seas", "Review connected inland seas and choose which ones to keep, fill, or drain.", GenerationStageDomain::terrain, true, true },
        { GenerationStageId::archipelagos, "Archipelagos", "Optionally add fragmented coastal islands and archipelago chains.", GenerationStageDomain::terrain, true, true },
        { GenerationStageId::ocean_cleanup, "Ocean cleanup", "Optionally tidy narrow channels and simplify open-ocean structure.", GenerationStageDomain::terrain, true, true },
        { GenerationStageId::coastline_refinement, "Coastline refinement", "Optionally break up straight shorelines and add more natural coastal variation.", GenerationStageDomain::terrain, true, true },
        { GenerationStageId::ocean_depth_refinement, "Ocean depth refinement", "Optionally deepen shelves and open-ocean basins after the coastline passes.", GenerationStageDomain::terrain, true, true },
        { GenerationStageId::tectonic_trenches, "Tectonic trenches", "Optionally imprint convergent-margin trench systems from tectonic metadata.", GenerationStageDomain::terrain, true, true },
        { GenerationStageId::tectonic_volcanoes, "Tectonic volcanoes", "Optionally place volcanic arcs and hotspot-style uplifts from tectonic signals.", GenerationStageDomain::terrain, true, true },
        { GenerationStageId::terrain_texturing, "Terrain texturing", "Refresh derived terrain masks and shelves after major height edits.", GenerationStageDomain::terrain, true, true },
        { GenerationStageId::mountain_bases, "Mountain base uplift", "Optionally broaden mountain roots and connect ranges to surrounding uplands.", GenerationStageDomain::terrain, true, true },
        { GenerationStageId::terrain_smoothing, "Terrain smoothing", "Optionally smooth terrain while preserving the large-scale landform structure.", GenerationStageDomain::terrain, true, true },
        { GenerationStageId::canyon_uplift, "Canyon uplift", "Optionally raise land around incised drainage corridors and canyon systems.", GenerationStageDomain::terrain, true, true },
        { GenerationStageId::depression_fill, "Depression fill", "Optionally fill closed depressions and reduce unwanted sinks.", GenerationStageDomain::terrain, true, true },
        { GenerationStageId::coastline_adjustment, "Coastline adjustment", "Optionally normalize coastal elevations and shallow seas.", GenerationStageDomain::terrain, true, true },
        { GenerationStageId::island_check, "Island check", "Optionally remove tiny isolated artifacts and stabilize island topology.", GenerationStageDomain::terrain, true, true },
        { GenerationStageId::terrain_roughness, "Terrain roughness", "Build the roughness field used by later terrain and climate layers.", GenerationStageDomain::terrain, true, false },
        { GenerationStageId::global_temperature, "Global temperature", "Generate the baseline temperature field from latitude, elevation, and orbital settings.", GenerationStageDomain::climate, true, false },
        { GenerationStageId::pressure, "Pressure", "Optionally generate the global pressure field.", GenerationStageDomain::climate, true, true },
        { GenerationStageId::winds, "Winds", "Optionally derive prevailing winds from pressure and planetary settings.", GenerationStageDomain::climate, true, true },
        { GenerationStageId::ocean_currents, "Ocean currents", "Solve wind-driven basin currents and the coupled seasonal SST response.", GenerationStageDomain::climate, true, true },
        { GenerationStageId::sea_surface_temperatures, "Sea surface temperatures", "Diagnose evaporation from the coupled ocean temperature and wind.", GenerationStageDomain::climate, true, true },
        { GenerationStageId::sea_ice_and_tides, "Sea ice and tides", "Optionally add sea-ice coverage and tidal effects to the climate model.", GenerationStageDomain::climate, true, true },
        { GenerationStageId::rainfall, "Rainfall", "Generate rainfall and moisture transport across the current terrain.", GenerationStageDomain::hydrology, true, false },
        { GenerationStageId::flow_routing, "Flow routing", "Route runoff over the terrain: receivers, drainage area, discharge, depressions and their spill points.", GenerationStageDomain::hydrology, true, true },
        { GenerationStageId::fluvial_erosion, "Fluvial erosion", "Erode and deposit along the drainage network with the stream power law over the erosion duration.", GenerationStageDomain::terrain, true, true },
        { GenerationStageId::fjords, "Fjords", "Optionally carve fjords along cold coastal mountain fronts.", GenerationStageDomain::hydrology, true, true },
        { GenerationStageId::rivers_and_basins, "Rivers and basins", "Plan rivers and drainage basins from the current terrain and rainfall.", GenerationStageDomain::hydrology, true, true },
        { GenerationStageId::basin_editor, "Basin editor", "Review detected endorheic basins and choose whether to keep, fill, or drain them.", GenerationStageDomain::hydrology, true, true },
        { GenerationStageId::rivers, "Rivers", "Rivers and their seasonal flows from the routing's discharge and flow directions.", GenerationStageDomain::hydrology, true, true },
        { GenerationStageId::lakes, "Lakes", "Optionally generate lakes and rift lakes from the current drainage network.", GenerationStageDomain::hydrology, true, true },
        { GenerationStageId::lake_levels, "Lake levels", "Lake levels in the routing's depressions from a water balance; salt lakes and pans in closed basins.", GenerationStageDomain::hydrology, true, true },
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
