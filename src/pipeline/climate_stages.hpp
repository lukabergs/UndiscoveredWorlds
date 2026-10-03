#pragma once

// Climate stages on the climate world alone, without planet. world_stages runs
// them on a planet's climate state; climate_pipeline runs them on imported
// terrain. Seeded noise streams follow the world seed and stage salts.

#include "generation_stage_registry.hpp"

class boolshapetemplate;
class climateworld;

struct ClimateStageContext
{
    boolshapetemplate* smalllake = nullptr;   // lake shapes for worlds without sea
    boolshapetemplate* landshape = nullptr;
};

// True for the stages that run in the earth-climate recipe.
bool isclimateonlystage(GenerationStageId stage);

// Runs one climate stage. Tides are hydrology and stay with the caller.
void executeclimatestage(GenerationStageId stage, climateworld& world, const ClimateStageContext& context);

// Land and sea cell counts of the bound terrain (margins excluded, as in
// getlandandseatotals); sets the climate inputs.
void countclimatelandandsea(climateworld& world);
