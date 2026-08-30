#pragma once

#include <vector>

class planet;

namespace climateenergy
{
struct AnnualEnergyBudget
{
    double incomingSolarWm2 = 0.0;
    double absorbedSolarWm2 = 0.0;
    double outgoingLongwaveWm2 = 0.0;
    double atmosphericTransportWm2 = 0.0;
    double storageTendencyWm2 = 0.0;
    double residualWm2 = 0.0;
    double calibratedLongwaveInterceptWm2 = 0.0;
    double areaWeightedMeanTemperatureC = 0.0;
    double areaWeightedLandTemperatureC = 0.0;
    double areaWeightedOceanTemperatureC = 0.0;
    double areaWeightedLandElevationCoolingC = 0.0;
};

float solarDeclinationRadians(float dayOfYear, float obliquityDegrees);
float orbitalDistanceFactor(float dayOfYear, float eccentricity, int perihelionSeason);
float dailyMeanInsolationWm2(
    float latitudeDegrees,
    float declinationRadians,
    float distanceFactor = 1.0f,
    float solarConstantWm2 = 1361.0f);
float implicitSlabTemperatureStep(
    float previousTemperatureC,
    float absorbedSolarWm2,
    float longwaveInterceptWm2,
    float longwaveSlopeWm2K,
    float transportCoefficientWm2K,
    float transportMeanTemperatureC,
    float heatCapacityJm2K,
    float timeStepSeconds);

void createSurfaceEnergyBalanceTemperatureMap(
    planet& world,
    const std::vector<std::vector<int>>& fractal);
const AnnualEnergyBudget& lastAnnualEnergyBudget();
}
