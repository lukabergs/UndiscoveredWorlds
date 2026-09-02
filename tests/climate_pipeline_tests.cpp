#include "generation_stage_registry.hpp"
#include <iostream>

int main()
{
    const auto temperature = getgenerationstageindex(GenerationStageId::global_temperature);
    const auto pressure = getgenerationstageindex(GenerationStageId::pressure);
    const auto winds = getgenerationstageindex(GenerationStageId::winds);
    const auto ocean = getgenerationstageindex(GenerationStageId::ocean_currents);
    const auto evaporation = getgenerationstageindex(GenerationStageId::sea_surface_temperatures);
    const auto rainfall = getgenerationstageindex(GenerationStageId::rainfall);
    if (!(temperature < pressure && pressure < winds && winds < ocean && ocean < evaporation && evaporation < rainfall))
    {
        std::cerr << "Climate stages must run temperature -> pressure -> winds -> ocean/SST -> evaporation -> rainfall\n";
        return 1;
    }
    return 0;
}
