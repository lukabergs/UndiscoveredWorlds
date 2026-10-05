#include "pipeline/generation_stage_registry.hpp"
#include <iostream>
#include <gtest/gtest.h>

TEST(ClimatePipeline, StageOrder)
{
    const auto temperature = getgenerationstageindex(GenerationStageId::global_temperature);
    const auto pressure = getgenerationstageindex(GenerationStageId::pressure);
    const auto winds = getgenerationstageindex(GenerationStageId::winds);
    const auto ocean = getgenerationstageindex(GenerationStageId::ocean_currents);
    const auto evaporation = getgenerationstageindex(GenerationStageId::sea_surface_temperatures);
    const auto rainfall = getgenerationstageindex(GenerationStageId::rainfall);
    ASSERT_TRUE(temperature < pressure && pressure < winds && winds < ocean && ocean < evaporation && evaporation < rainfall)
        << "Climate stages must run temperature -> pressure -> winds -> ocean/SST -> evaporation -> rainfall";
}
