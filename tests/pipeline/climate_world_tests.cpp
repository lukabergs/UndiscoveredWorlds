// Links only uw_climate_world: the climate stage runner and the climate
// validation report must build without planet, terrain, hydrology or tectonics.

#include <iostream>

#include "climate_stages.hpp"

class climateworld;
void exportclimatevalidationreport(climateworld& world);

int main()
{
    void (*const runstage)(GenerationStageId, climateworld&, const ClimateStageContext&) = &executeclimatestage;
    void (*const report)(climateworld&) = &exportclimatevalidationreport;
    if (runstage == nullptr || report == nullptr)
        return 1;
    if (!isclimateonlystage(GenerationStageId::rainfall) || isclimateonlystage(GenerationStageId::terrain_roughness))
    {
        std::cerr << "The climate-only stages are the climate feedback stages without terrain roughness\n";
        return 1;
    }
    return 0;
}
