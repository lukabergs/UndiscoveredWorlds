#pragma once

// Explicit climate feedback subset; excludes all terrain editing and society.

#include <cstddef>
#include <vector>

enum class GenerationStageDomain
{
    terrain,
    climate,
    hydrology,
    finalization
};

enum class GenerationStageId
{
    geosphere_initialization,
    plate_tectonics,
    terraforming,
    inland_seas,
    archipelagos,
    ocean_cleanup,
    coastline_refinement,
    ocean_depth_refinement,
    tectonic_trenches,
    tectonic_volcanoes,
    terrain_texturing,
    mountain_bases,
    terrain_smoothing,
    canyon_uplift,
    depression_fill,
    coastline_adjustment,
    island_check,
    terrain_roughness,
    global_temperature,
    ocean_currents,
    sea_surface_temperatures,
    pressure,
    winds,
    sea_ice_and_tides,
    rainfall,
    flow_routing,
    fluvial_erosion,
    fjords,
    rivers,
    lake_levels,
    mountain_temperature_lapse,
    climates_and_biomes,
    arid_features,
    deltas_wetlands_and_roughness,
    finalize_layers
};

struct GenerationStageDefinition
{
    GenerationStageId id;
    const char* label;
    const char* description;
    GenerationStageDomain domain;
    bool previewEnabledByDefault;
    bool skippable;
};

const std::vector<GenerationStageDefinition>& getgenerationstages();
const GenerationStageDefinition& getgenerationstage(GenerationStageId id);
std::size_t getgenerationstageindex(GenerationStageId id);
bool isterraingenerationstage(GenerationStageId id);
bool isclimategenerationstage(GenerationStageId id);
bool isclimatefeedbackstage(GenerationStageId id);
