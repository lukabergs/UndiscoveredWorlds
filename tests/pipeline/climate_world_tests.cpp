// Links only uw_climate_world: the climate stage runner and the climate
// validation report must build without planet, terrain, hydrology or tectonics.

#include <iostream>

#include "climate_stages.hpp"
#include <gtest/gtest.h>

class climateworld;
void exportclimatevalidationreport(climateworld& world);

TEST(ClimateWorld, BuildsWithoutPlanet)
{
    void (*const runstage)(GenerationStageId, climateworld&, const ClimateStageContext&) = &executeclimatestage;
    void (*const report)(climateworld&) = &exportclimatevalidationreport;
    ASSERT_TRUE(runstage != nullptr && report != nullptr);
    EXPECT_TRUE(isclimateonlystage(GenerationStageId::rainfall) && !isclimateonlystage(GenerationStageId::terrain_roughness))
        <<"The climate-only stages are the climate feedback stages without terrain roughness";
}
