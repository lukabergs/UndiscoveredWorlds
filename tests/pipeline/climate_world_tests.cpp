// Links only uw_climate_world: the climate stage runner and the climate
// diagnostics must build without planet, terrain, hydrology or tectonics.

#include <iostream>

#include "pipeline/climate_stages.hpp"
#include <gtest/gtest.h>

#include "simulations/climate/climate_diagnostics.hpp"

TEST(ClimateWorld, BuildsWithoutPlanet)
{
    void (*const runstage)(GenerationStageId, climateworld&, const ClimateStageContext&) = &executeclimatestage;
    void (*const report)(climateworld&) = &writeclimatediagnostics;
    ASSERT_TRUE(runstage != nullptr && report != nullptr);
    EXPECT_TRUE(isclimateonlystage(GenerationStageId::rainfall) && !isclimateonlystage(GenerationStageId::terrain_roughness))
        <<"The climate-only stages are the climate feedback stages without terrain roughness";
}
