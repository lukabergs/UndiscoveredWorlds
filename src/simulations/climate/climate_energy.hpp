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
struct SolarForcing
{
    float declinationRadians = 0.0f;
    float distanceFactor = 1.0f;
};
SolarForcing solarForcing(float dayOfYear, float obliquityDegrees,
    float eccentricity, int perihelionSeason);
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

struct LandThermalState
{
    double surfaceC = 0.0, deepC = 0.0;
};
// Backward-Euler internal exchange; net surface flux is positive into land.
// Conserves Cs*Ts + Cd*Td, including the supplied heat over the timestep.
LandThermalState stepLandThermalState(LandThermalState state,
    double netSurfaceFluxWm2, double surfaceCapacityJm2K,
    double deepCapacityJm2K, double conductanceWm2K, double timeStepSeconds);

struct MeridionalDiffusionResult
{
    std::vector<double> temperatureC, heatConvergenceWm2;
};
// Solve T = source + response * divergence(D grad(T)) on cell-centred
// spherical latitude bands, with zero polar heat flux.
MeridionalDiffusionResult implicitMeridionalDiffusion(
    const std::vector<double>& sourceTemperatureC,
    const std::vector<double>& responseKPerWm2, double diffusivityWm2K);

void createSurfaceEnergyBalanceTemperatureMap(
    planet& world,
    const std::vector<std::vector<int>>& fractal);
const AnnualEnergyBudget& lastAnnualEnergyBudget();
}
