#include "climate_energy.hpp"

#include "generation_tuning.hpp"
#include "physical_layers.hpp"
#include "planet.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace climateenergy
{
namespace
{
constexpr double pi = 3.14159265358979323846;
constexpr int daysPerYear = 365;
constexpr int spinupYears = 12;
constexpr int calibrationPasses = 20;
constexpr std::array<int, CLIMATESEASONCOUNT> snapshotDays = { 15, 105, 196, 288 };

AnnualEnergyBudget annualEnergyBudget{};

double latitudeDegreesForRow(int row, int height)
{
    return height > 0 ? 90.0 - 180.0 * static_cast<double>(row) / static_cast<double>(height) : 0.0;
}

double gridCellAreaWeight(int row, int height)
{
    return std::max(0.0, std::cos(latitudeDegreesForRow(row, height) * pi / 180.0));
}

double effectiveAlbedo(bool land, double temperatureC)
{
    const double openSurfaceAlbedo = land ?
        tuning::climate::energybalance::landAlbedo :
        tuning::climate::energybalance::oceanAlbedo;
    const double frozenAlbedo = land ?
        tuning::climate::energybalance::snowAlbedo :
        tuning::climate::energybalance::seaIceAlbedo;
    const double freezingPoint = land ? -2.0 : -1.8;
    const double transitionRange = land ? 15.0 : 6.0;
    const double frozenFraction = std::clamp(
        (freezingPoint - temperatureC) / transitionRange,
        0.0,
        1.0);
    return openSurfaceAlbedo + (frozenAlbedo - openSurfaceAlbedo) * frozenFraction;
}

struct ProfileSimulation
{
    std::array<std::array<std::vector<double>, 2>, CLIMATESEASONCOUNT> snapshots;
    AnnualEnergyBudget budget;
    double annualMeanTemperatureC = 0.0;
};

ProfileSimulation simulateProfiles(
    const planet& world,
    const std::vector<int>& landCells,
    const std::vector<int>& oceanCells,
    double longwaveInterceptWm2,
    bool collectOutput)
{
    const int height = world.height();
    const double targetMeanTemperature = static_cast<double>(world.averagetemp());
    const double timeStepSeconds = 86400.0;
    const double longwaveSlope = tuning::climate::energybalance::longwaveSlopeWm2K;
    const double meridionalTransport = tuning::climate::energybalance::meridionalTransportWm2K;
    const double zonalExchange = tuning::climate::energybalance::zonalLandOceanExchangeWm2K;
    const std::array<double, 2> heatCapacity = {
        tuning::climate::energybalance::landHeatCapacityJm2K,
        tuning::climate::energybalance::oceanMixedLayerHeatCapacityJm2K
    };
    std::array<std::vector<double>, 2> temperature = {
        std::vector<double>(height + 1, targetMeanTemperature),
        std::vector<double>(height + 1, targetMeanTemperature)
    };
    std::array<std::vector<double>, 2> nextTemperature = temperature;
    std::array<std::vector<double>, 2> absorbedSolar = temperature;
    std::array<std::vector<double>, 2> baseTemperature = temperature;
    std::array<std::vector<double>, 2> weights = temperature;
    std::vector<double> zonalMeanTemperature(height + 1, targetMeanTemperature);
    std::vector<double> nextZonalMeanTemperature = zonalMeanTemperature;
    std::vector<std::vector<double>> insolation(daysPerYear, std::vector<double>(height + 1, 0.0));
    ProfileSimulation result;

    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
    {
        for (int surface = 0; surface < 2; surface++)
            result.snapshots[season][surface].assign(height + 1, targetMeanTemperature);
    }

    double totalAreaWeight = 0.0;

    for (int y = 0; y <= height; y++)
    {
        const double cellAreaWeight = gridCellAreaWeight(y, height);
        weights[0][y] = cellAreaWeight * static_cast<double>(landCells[y]);
        weights[1][y] = cellAreaWeight * static_cast<double>(oceanCells[y]);
        totalAreaWeight += weights[0][y] + weights[1][y];
        const double latitude = latitudeDegreesForRow(y, height);
        const double polarFraction = std::pow(std::abs(latitude) / 90.0, 4.0);
        const double polarAdjustment = latitude >= 0.0 ?
            static_cast<double>(world.northpolaradjust()) :
            static_cast<double>(world.southpolaradjust());

        temperature[0][y] += polarAdjustment * polarFraction;
        temperature[1][y] += polarAdjustment * polarFraction;

        for (int day = 0; day < daysPerYear; day++)
        {
            const float dayOfYear = static_cast<float>(day) + 0.5f;
            insolation[day][y] = dailyMeanInsolationWm2(
                static_cast<float>(latitude),
                solarDeclinationRadians(dayOfYear, world.tilt()),
                orbitalDistanceFactor(dayOfYear, world.eccentricity(), world.perihelion()));
        }
    }

    double lastYearMeanTemperatureSum = 0.0;
    std::array<double, 2> lastYearSurfaceTemperatureSum{};
    std::array<double, 2> surfaceAreaWeight{};
    double incomingSum = 0.0;
    double absorbedSum = 0.0;
    double outgoingSum = 0.0;
    double transportSum = 0.0;
    double storageSum = 0.0;

    for (int year = 0; year < spinupYears; year++)
    {
        for (int day = 0; day < daysPerYear; day++)
        {
            for (int y = 0; y <= height; y++)
            {
                const double latitude = latitudeDegreesForRow(y, height);
                const double polarFraction = std::pow(std::abs(latitude) / 90.0, 4.0);
                const double polarAdjustment = latitude >= 0.0 ?
                    static_cast<double>(world.northpolaradjust()) :
                    static_cast<double>(world.southpolaradjust());
                const double polarForcing = polarAdjustment * polarFraction *
                    (longwaveSlope + meridionalTransport + zonalExchange);

                for (int surface = 0; surface < 2; surface++)
                {
                    const bool land = surface == 0;
                    const double albedo = effectiveAlbedo(land, temperature[surface][y]);
                    absorbedSolar[surface][y] =
                        insolation[day][y] * (1.0 - albedo) + polarForcing;
                    const double capacityRate = heatCapacity[surface] / timeStepSeconds;
                    const double denominator = capacityRate + longwaveSlope + meridionalTransport + zonalExchange;
                    baseTemperature[surface][y] =
                        (capacityRate * temperature[surface][y] +
                            absorbedSolar[surface][y] - longwaveInterceptWm2) /
                        denominator;
                }
            }

            double globalMeanTemperature = 0.0;

            for (int y = 0; y <= height; y++)
            {
                const double rowWeight = weights[0][y] + weights[1][y];

                if (rowWeight > 0.0)
                {
                    zonalMeanTemperature[y] =
                        (weights[0][y] * temperature[0][y] + weights[1][y] * temperature[1][y]) /
                        rowWeight;
                    globalMeanTemperature += rowWeight * zonalMeanTemperature[y];
                }
            }

            globalMeanTemperature /= totalAreaWeight;

            for (int couplingIteration = 0; couplingIteration < 12; couplingIteration++)
            {
                double nextGlobalMeanTemperature = 0.0;

                for (int y = 0; y <= height; y++)
                {
                    const double rowWeight = weights[0][y] + weights[1][y];

                    for (int surface = 0; surface < 2; surface++)
                    {
                        const double capacityRate = heatCapacity[surface] / timeStepSeconds;
                        const double denominator = capacityRate + longwaveSlope + meridionalTransport + zonalExchange;
                        nextTemperature[surface][y] = baseTemperature[surface][y] +
                            (zonalExchange * zonalMeanTemperature[y] +
                                meridionalTransport * globalMeanTemperature) /
                            denominator;
                    }

                    if (rowWeight > 0.0)
                    {
                        nextZonalMeanTemperature[y] =
                            (weights[0][y] * nextTemperature[0][y] +
                                weights[1][y] * nextTemperature[1][y]) /
                            rowWeight;
                        nextGlobalMeanTemperature += rowWeight * nextZonalMeanTemperature[y];
                    }
                }

                zonalMeanTemperature.swap(nextZonalMeanTemperature);
                globalMeanTemperature = nextGlobalMeanTemperature / totalAreaWeight;
            }

            if (year == spinupYears - 1)
            {
                lastYearMeanTemperatureSum += globalMeanTemperature;

                for (int y = 0; y <= height; y++)
                {
                    for (int surface = 0; surface < 2; surface++)
                        lastYearSurfaceTemperatureSum[surface] += weights[surface][y] * nextTemperature[surface][y];
                }

                if (collectOutput)
                {
                    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
                    {
                        if (day + 1 == snapshotDays[season])
                        {
                            result.snapshots[season][0] = nextTemperature[0];
                            result.snapshots[season][1] = nextTemperature[1];
                        }
                    }

                    for (int y = 0; y <= height; y++)
                    {
                        for (int surface = 0; surface < 2; surface++)
                        {
                            const double weight = weights[surface][y];
                            const double transport =
                                zonalExchange * (zonalMeanTemperature[y] - nextTemperature[surface][y]) +
                                meridionalTransport * (globalMeanTemperature - nextTemperature[surface][y]);
                            const double storage = heatCapacity[surface] *
                                (nextTemperature[surface][y] - temperature[surface][y]) /
                                timeStepSeconds;
                            incomingSum += weight * insolation[day][y];
                            absorbedSum += weight * absorbedSolar[surface][y];
                            outgoingSum += weight *
                                (longwaveInterceptWm2 + longwaveSlope * nextTemperature[surface][y]);
                            transportSum += weight * transport;
                            storageSum += weight * storage;
                        }
                    }
                }
            }

            temperature.swap(nextTemperature);
        }
    }

    result.annualMeanTemperatureC = lastYearMeanTemperatureSum / static_cast<double>(daysPerYear);

    for (int surface = 0; surface < 2; surface++)
    {
        for (int y = 0; y <= height; y++)
            surfaceAreaWeight[surface] += weights[surface][y];
    }

    if (collectOutput)
    {
        const double normalization = totalAreaWeight * static_cast<double>(daysPerYear);
        result.budget.incomingSolarWm2 = incomingSum / normalization;
        result.budget.absorbedSolarWm2 = absorbedSum / normalization;
        result.budget.outgoingLongwaveWm2 = outgoingSum / normalization;
        result.budget.atmosphericTransportWm2 = transportSum / normalization;
        result.budget.storageTendencyWm2 = storageSum / normalization;
        result.budget.residualWm2 =
            result.budget.absorbedSolarWm2 -
            result.budget.outgoingLongwaveWm2 +
            result.budget.atmosphericTransportWm2 -
            result.budget.storageTendencyWm2;
        result.budget.calibratedLongwaveInterceptWm2 = longwaveInterceptWm2;
        result.budget.areaWeightedMeanTemperatureC = result.annualMeanTemperatureC;
        result.budget.areaWeightedLandTemperatureC = lastYearSurfaceTemperatureSum[0] /
            (surfaceAreaWeight[0] * static_cast<double>(daysPerYear));
        result.budget.areaWeightedOceanTemperatureC = lastYearSurfaceTemperatureSum[1] /
            (surfaceAreaWeight[1] * static_cast<double>(daysPerYear));
    }

    return result;
}
}

float solarDeclinationRadians(float dayOfYear, float obliquityDegrees)
{
    const double orbitalLongitude = 2.0 * pi *
        (static_cast<double>(dayOfYear) - 80.0) / static_cast<double>(daysPerYear);
    const double obliquity = static_cast<double>(obliquityDegrees) * pi / 180.0;
    return static_cast<float>(std::asin(std::sin(obliquity) * std::sin(orbitalLongitude)));
}

float orbitalDistanceFactor(float dayOfYear, float eccentricity, int perihelionSeason)
{
    const double clampedEccentricity = std::clamp(static_cast<double>(eccentricity), 0.0, 0.95);
    const double perihelionDay = perihelionSeason == 1 ? 186.0 : 3.0;
    const double anomaly = 2.0 * pi *
        (static_cast<double>(dayOfYear) - perihelionDay) / static_cast<double>(daysPerYear);
    const double numerator = 1.0 + clampedEccentricity * std::cos(anomaly);
    const double denominator = 1.0 - clampedEccentricity * clampedEccentricity;
    return static_cast<float>((numerator * numerator) / (denominator * denominator));
}

float dailyMeanInsolationWm2(
    float latitudeDegrees,
    float declinationRadians,
    float distanceFactor,
    float solarConstantWm2)
{
    const double latitude = std::clamp(static_cast<double>(latitudeDegrees), -90.0, 90.0) * pi / 180.0;
    const double declination = static_cast<double>(declinationRadians);
    const double sunsetArgument = -std::tan(latitude) * std::tan(declination);
    double sunsetHourAngle = 0.0;

    if (sunsetArgument <= -1.0)
        sunsetHourAngle = pi;
    else if (sunsetArgument < 1.0)
        sunsetHourAngle = std::acos(sunsetArgument);

    const double dailyMean = static_cast<double>(solarConstantWm2) *
        static_cast<double>(distanceFactor) / pi *
        (sunsetHourAngle * std::sin(latitude) * std::sin(declination) +
            std::cos(latitude) * std::cos(declination) * std::sin(sunsetHourAngle));
    return static_cast<float>(std::max(0.0, dailyMean));
}

float implicitSlabTemperatureStep(
    float previousTemperatureC,
    float absorbedSolarWm2,
    float longwaveInterceptWm2,
    float longwaveSlopeWm2K,
    float transportCoefficientWm2K,
    float transportMeanTemperatureC,
    float heatCapacityJm2K,
    float timeStepSeconds)
{
    const double capacityRate = static_cast<double>(heatCapacityJm2K) /
        std::max(1.0, static_cast<double>(timeStepSeconds));
    const double denominator = capacityRate + static_cast<double>(longwaveSlopeWm2K) +
        static_cast<double>(transportCoefficientWm2K);
    return static_cast<float>(
        (capacityRate * static_cast<double>(previousTemperatureC) +
            static_cast<double>(absorbedSolarWm2) -
            static_cast<double>(longwaveInterceptWm2) +
            static_cast<double>(transportCoefficientWm2K) *
                static_cast<double>(transportMeanTemperatureC)) /
        denominator);
}

void createSurfaceEnergyBalanceTemperatureMap(
    planet& world,
    const std::vector<std::vector<int>>& fractal)
{
    const int width = world.width();
    const int height = world.height();
    const int seaLevel = world.sealevel();
    std::vector<int> landCells(height + 1, 0);
    std::vector<int> oceanCells(height + 1, 0);

    for (int y = 0; y <= height; y++)
    {
        for (int x = 0; x <= width; x++)
        {
            if (world.sea(x, y) == 1)
                oceanCells[y]++;
            else
                landCells[y]++;
        }
    }

    double longwaveIntercept = tuning::climate::energybalance::earthLongwaveInterceptWm2;

    for (int pass = 0; pass < calibrationPasses; pass++)
    {
        const ProfileSimulation calibration = simulateProfiles(
            world, landCells, oceanCells, longwaveIntercept, false);
        const double temperatureError = calibration.annualMeanTemperatureC -
            static_cast<double>(world.averagetemp());
        longwaveIntercept += temperatureError * tuning::climate::energybalance::longwaveSlopeWm2K;
    }

    const ProfileSimulation simulation = simulateProfiles(
        world, landCells, oceanCells, longwaveIntercept, true);
    annualEnergyBudget = simulation.budget;
    const double fractalScale = static_cast<double>(std::max(1, world.maxelevation()));
    double landElevationCoolingSum = 0.0;
    double landAreaWeight = 0.0;

    for (int y = 0; y <= height; y++)
    {
        for (int x = 0; x <= width; x++)
        {
            const int surface = world.sea(x, y) == 1 ? 1 : 0;
            const double fractalAnomaly =
                (static_cast<double>(fractal[x][y]) / fractalScale - 0.5) *
                tuning::climate::energybalance::fractalTemperatureRangeC;
            const double elevationKm = static_cast<double>(std::max(0, world.map(x, y) - seaLevel)) / 1000.0;
            const double elevationCooling = surface == 0 ?
                elevationKm * static_cast<double>(world.tempdecrease()) : 0.0;

            if (surface == 0)
            {
                const double areaWeight = gridCellAreaWeight(y, height);
                landElevationCoolingSum += areaWeight * elevationCooling;
                landAreaWeight += areaWeight;
            }

            for (int season = 0; season < CLIMATESEASONCOUNT; season++)
            {
                const double temperature = simulation.snapshots[season][surface][y] +
                    fractalAnomaly - elevationCooling;
                world.setseasonaltemp(season, x, y, static_cast<int>(std::round(temperature)));
            }

            world.setjantemp(x, y, world.seasonaltemp(seasonjanuary, x, y));
            world.setjultemp(x, y, world.seasonaltemp(seasonjuly, x, y));
        }
    }

    auto rowAnnualMean = [&](int y)
    {
        double total = 0.0;

        for (int season = 0; season < CLIMATESEASONCOUNT; season++)
        {
            const double cells = static_cast<double>(landCells[y] + oceanCells[y]);
            total += cells > 0.0 ?
                (simulation.snapshots[season][0][y] * static_cast<double>(landCells[y]) +
                    simulation.snapshots[season][1][y] * static_cast<double>(oceanCells[y])) / cells :
                static_cast<double>(world.averagetemp());
        }

        return total / static_cast<double>(CLIMATESEASONCOUNT);
    };

    world.setnorthpolartemp(static_cast<int>(std::round(rowAnnualMean(0))));
    world.setsouthpolartemp(static_cast<int>(std::round(rowAnnualMean(height))));
    world.seteqtemp(static_cast<int>(std::round(rowAnnualMean(height / 2))));
    annualEnergyBudget.areaWeightedLandElevationCoolingC =
        landAreaWeight > 0.0 ? landElevationCoolingSum / landAreaWeight : 0.0;
}

const AnnualEnergyBudget& lastAnnualEnergyBudget()
{
    return annualEnergyBudget;
}
}
