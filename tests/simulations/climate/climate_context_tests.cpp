#include "climate_context.hpp"
#include <iostream>

int main()
{
    climatecontext::Context first;
    if (first.transientWeather || first.mixedLayerTropics || first.transportMethod != climatecontext::TransportMethod::Ffsl)
    {
        std::cerr << "baseline climate must default to legacy stationary circulation and FFSL\n";
        return 1;
    }
    first.transientWeather = true;
    first.mixedLayerTropics = true;
    first.boundaryLayer[0].thermalPressureHpa = {1.0f, 2.0f};
    first.radiativeZonalTemperature[0] = {10.0f, 20.0f, 20.0f, 10.0f};
    first.transportMethod = climatecontext::TransportMethod::Mpdata;
    first.resolution.hydrologyColumns = 32;
    first.hydrologyYears = 7;
    first.couplingIterationLimit = 8;
    first.storage = {8, 4, {8, 4, 1.0f}, {8, 4, 2.0f}, {8, 4, 3.0f}, {8, 4, 4.0f},
        {8, 4, 10.0}, {8, 4, 8.0}};
    auto& monthly = first.monthly;
    monthly.columns = 8; monthly.rows = 4;
    for (auto& field : monthly.temperatureC) field.assign(32, 12.0f);
    for (auto& field : monthly.rainfallMm) field.assign(32, 40.0f);
    for (auto& field : monthly.referenceTemperatureC) field.assign(32, 12.0f);
    for (auto& field : monthly.referenceRainfallMm) field.assign(32, 40.0f);
    first.winds.reset(8, 4);
    first.winds.drainagewindu[0] = climatecontext::FloatGrid(8, std::vector<float>(4, 2.0f));
    first.winds.hasDrainageProfile[0] = true;
    first.winds.lowlevelwindu[0] = climatecontext::FloatGrid(8, std::vector<float>(4, 7.0f));
    first.winds.hasLowLevelProfile[0] = true;
    auto snapshot = first;
    first.winds.drainagewindu[0][7][3] = 99.0f;
    first.winds.lowlevelwindu[0][7][3] = 99.0f;
    first.winds.reset(8, 4);
    if (first.winds.hasLowLevelProfile[0] || !snapshot.winds.hasLowLevelProfile[0] ||
        snapshot.winds.lowlevelwindu[0][7][3] != 7.0f)
    {
        std::cerr << "low-level wind caches must reset validity and copy independently\n";
        return 1;
    }
    if (first.winds.hasDrainageProfile[0] || !snapshot.winds.hasDrainageProfile[0] ||
        snapshot.winds.drainagewindu[0][7][3] != 2.0f)
    {
        std::cerr << "wind resets must invalidate drainage profiles; snapshots must own their increments\n";
        return 1;
    }
    first.storage.boundary[7][3] = 99.0f;
    monthly.rainfallMm[2].clear();
    if (snapshot.storage.boundary[7][3] != 1.0f || !snapshot.monthly.hasRainfall() ||
        monthly.hasRainfall() || !monthly.hasTemperature())
    {
        std::cerr << "world snapshots must isolate reservoirs and reject incomplete monthly fields\n";
        return 1;
    }
    first.resetDynamics();
    if (first.storage.columns != 0 || first.hasOcean || first.hasProcesses || first.hydrologyYears != 0 ||
        snapshot.hydrologyYears != 7 || first.couplingIterationLimit != 8 ||
        !first.mixedLayerTropics || !snapshot.mixedLayerTropics || !first.boundaryLayer[0].thermalPressureHpa.empty() ||
        snapshot.boundaryLayer[0].thermalPressureHpa.size() != 2 ||
        !first.transientWeather || first.transportMethod != climatecontext::TransportMethod::Mpdata ||
        !first.radiativeZonalTemperature[0].empty() || snapshot.radiativeZonalTemperature[0].size() != 4 ||
        !snapshot.transientWeather || snapshot.transportMethod != climatecontext::TransportMethod::Mpdata ||
        first.resolution.hydrologyColumns != 32 || !first.monthly.hasTemperature())
    {
        std::cerr << "dynamical restart must retain configuration and the base temperature climatology\n";
        return 1;
    }
    return 0;
}
