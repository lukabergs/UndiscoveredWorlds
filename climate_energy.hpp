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
    double areaWeightedPermanentLandIceFraction = 0.0;
    double localIceCouplingMaximumRowResidual = 0.0;
    double localIceCouplingTemperatureErrorC = 0.0;
    double northern5070WarmestTemperatureC = 0.0;
    double southern6090WarmestTemperatureC = 0.0;
    double southern6090IcecapThermalFraction = 0.0;
    int localIceCouplingIterations = 0;
    bool localIceCouplingConverged = false;
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
float permanentLandIceFraction(float warmestSeasonTemperatureC);

void createSurfaceEnergyBalanceTemperatureMap(
    planet& world,
    const std::vector<std::vector<int>>& fractal);
const AnnualEnergyBudget& lastAnnualEnergyBudget();
}
