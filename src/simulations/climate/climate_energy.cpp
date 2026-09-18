#include "climate_energy.hpp"

#include "climate_grid.hpp"
#include "detail/land_ice.hpp"
#include "detail/tropical_cloud_albedo.hpp"
#include "generation_tuning.hpp"
#include "planet.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

namespace climateenergy
{
LandThermalState stepLandThermalState(LandThermalState state,
    double flux, double surfaceCapacity, double deepCapacity,
    double conductance, double dt)
{
    if (!(surfaceCapacity > 0.0) || !(deepCapacity > 0.0) ||
        !(dt > 0.0) || !(conductance >= 0.0) ||
        !std::isfinite(surfaceCapacity) || !std::isfinite(deepCapacity) ||
        !std::isfinite(conductance) || !std::isfinite(dt) ||
        !std::isfinite(flux) || !std::isfinite(state.surfaceC) || !std::isfinite(state.deepC))
        return state;
    const double surfaceRate = surfaceCapacity / dt;
    const double deepRate = deepCapacity / dt;
    const double exchange = conductance * deepRate / (deepRate + conductance);
    const double surface = (surfaceRate * state.surfaceC + flux + exchange * state.deepC) /
        (surfaceRate + exchange);
    return {surface, (deepRate * state.deepC + conductance * surface) / (deepRate + conductance)};
}

namespace
{
constexpr double pi = 3.14159265358979323846;
constexpr int daysPerYear = 365;
constexpr int spinupYears = 12;
constexpr std::array<int, CLIMATESEASONCOUNT> snapshotDays = { 15, 105, 196, 288 };

AnnualEnergyBudget annualEnergyBudget{};

double latitudeDegreesForRow(int row, int height)
{
    return climategrid::latitudeCentreRadians(row, height + 1) * 180.0 / pi;
}

double gridCellAreaWeight(int row, int height)
{
    return climategrid::latitudeBandMeasure(row, height + 1);
}

double smoothstep01(double value)
{
    const double clamped = std::clamp(value, 0.0, 1.0);
    return clamped * clamped * (3.0 - 2.0 * clamped);
}

double effectiveAlbedo(
    bool land,
    double temperatureC,
    double permanentIceFraction,
    double latitudeDegrees)
{
    const double openSurfaceAlbedo = land ?
        tuning::climate::energybalance::landAlbedo :
        tuning::climate::energybalance::oceanAlbedo;
    double frozenAlbedo = land ?
        tuning::climate::energybalance::snowAlbedo :
        tuning::climate::energybalance::seaIceAlbedo;
    const double freezingPoint = land ? -2.0 : -1.8;
    const double transitionRange = land ? 15.0 : 6.0;
    double frozenFraction = std::clamp(
        (freezingPoint - temperatureC) / transitionRange,
        0.0,
        1.0);

    if (land)
    {
        permanentIceFraction = std::clamp(permanentIceFraction, 0.0, 1.0);
        frozenFraction = std::max(frozenFraction, permanentIceFraction);
        frozenAlbedo +=
            (tuning::climate::energybalance::permanentIceAlbedo - frozenAlbedo) *
            permanentIceFraction;
    }

    const double cloud = detail::tropicalCloudAlbedo(latitudeDegrees,
        tuning::climate::energybalance::tropicalCloudAlbedoAmplitude,
        tuning::climate::energybalance::tropicalCloudAlbedoWidthDegrees);
    return std::clamp(openSurfaceAlbedo + (frozenAlbedo - openSurfaceAlbedo) * frozenFraction + cloud, 0.0, 1.0);
}

struct ProfileSimulation
{
    std::array<std::array<std::vector<double>, 2>, CLIMATESEASONCOUNT> snapshots;
    std::array<std::array<std::vector<double>, 2>, 12> monthly;
    AnnualEnergyBudget budget;
    double annualMeanTemperatureC = 0.0;
};

struct CoupledProfileSimulation
{
    ProfileSimulation profile;
    std::vector<double> permanentIceFraction;
    std::array<std::vector<double>, CLIMATESEASONCOUNT> localTemperatureCorrection;
    int couplingIterations = 0;
    double maximumRowIceFractionResidual = 0.0;
    double temperatureErrorC = 0.0;
    bool converged = false;
};

ProfileSimulation simulateProfiles(
    const planet& world,
    const std::vector<int>& outputLandCells,
    const std::vector<int>& outputOceanCells,
    const std::vector<double>& outputPermanentIceFraction,
    double longwaveInterceptWm2,
    bool collectOutput)
{
    // Prognostic energy profiles and world fields share cell-centred bands.
    const int rows = std::max(1, std::min(world.width() + 1,
        world.climatesimulation().resolution.hydrologyColumns) / 2);
    const int height = rows - 1;
    const auto floats = [](const auto& field)
    {
        std::vector<float> result(field.size());
        std::transform(field.begin(), field.end(), result.begin(), [](auto value) { return static_cast<float>(value); });
        return result;
    };
    const auto remapInput = [&](const auto& field)
    {
        return climategrid::remapField(
            1, world.height() + 1, floats(field), 1, rows);
    };
    const auto landCells = remapInput(outputLandCells);
    const auto oceanCells = remapInput(outputOceanCells);
    std::vector<double> outputIceArea(world.height() + 1, 0.0);
    std::vector<double> outputElevationCooling(world.height() + 1, 0.0);
    for (int y = 0; y <= world.height(); y++)
    {
        outputIceArea[y] = outputPermanentIceFraction[y] * outputLandCells[y];
        for (int x = 0; x <= world.width(); x++)
        {
            if (world.sea(x, y) == 0)
                outputElevationCooling[y] += std::max(0, world.map(x, y) - world.sealevel()) *
                    static_cast<double>(world.tempdecrease()) / 1000.0;
        }
    }
    auto rowPermanentIceFraction = remapInput(outputIceArea);
    auto meanLandElevationCooling = remapInput(outputElevationCooling);
    for (int y = 0; y < rows; y++)
    {
        const float denominator = std::max(1.0e-10f, landCells[y]);
        rowPermanentIceFraction[y] /= denominator;
        meanLandElevationCooling[y] /= denominator;
    }
    const auto latitudeForCell = [rows](int y)
    {
        return climategrid::latitudeCentreRadians(y, rows) * 180.0 / pi;
    };
    const double targetMeanTemperature = static_cast<double>(world.averagetemp());
    const double timeStepSeconds = 86400.0;
    const double longwaveSlope = tuning::climate::energybalance::longwaveSlopeWm2K;
    const double meridionalTransport = tuning::climate::energybalance::meridionalTransportWm2K;
    const double zonalExchange = tuning::climate::energybalance::zonalLandOceanExchangeWm2K;
    const double deepLandCapacity = tuning::climate::energybalance::landDeepLayerHeatCapacityJm2K;
    const double deepLandCoupling = tuning::climate::energybalance::landDeepLayerCouplingWm2K;
    const std::array<double, 2> heatCapacity = {
        tuning::climate::energybalance::landHeatCapacityJm2K,
        tuning::climate::energybalance::oceanMixedLayerHeatCapacityJm2K
    };
    std::array<std::vector<double>, 2> temperature = {
        std::vector<double>(height + 1, targetMeanTemperature),
        std::vector<double>(height + 1, targetMeanTemperature)
    };
    std::array<std::vector<double>, 2> nextTemperature = temperature;
    std::vector<double> deepLandTemperature(height + 1, targetMeanTemperature);
    std::vector<double> nextDeepLandTemperature = deepLandTemperature;
    std::array<std::vector<double>, 2> absorbedSolar = temperature;
    std::array<std::vector<double>, 2> weights = temperature;
    std::vector<double> zonalMeanTemperature(height + 1, targetMeanTemperature);
    std::vector<std::vector<double>> insolation(daysPerYear, std::vector<double>(height + 1, 0.0));
    ProfileSimulation result;
    for (auto& month : result.monthly)
        for (auto& surface : month) surface.assign(rows, 0.0);

    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
    {
        for (int surface = 0; surface < 2; surface++)
            result.snapshots[season][surface].assign(height + 1, targetMeanTemperature);
    }

    double totalAreaWeight = 0.0;

    for (int y = 0; y <= height; y++)
    {
        const double cellAreaWeight = climategrid::latitudeBandMeasure(y, rows);
        weights[0][y] = cellAreaWeight * static_cast<double>(landCells[y]);
        weights[1][y] = cellAreaWeight * static_cast<double>(oceanCells[y]);
        totalAreaWeight += weights[0][y] + weights[1][y];

        const double latitude = latitudeForCell(y);
        const double polarFraction = std::pow(std::abs(latitude) / 90.0, 4.0);
        const double polarAdjustment = latitude >= 0.0 ?
            static_cast<double>(world.northpolaradjust()) :
            static_cast<double>(world.southpolaradjust());

        temperature[0][y] += polarAdjustment * polarFraction;
        temperature[1][y] += polarAdjustment * polarFraction;
        deepLandTemperature[y] = temperature[0][y];
        nextDeepLandTemperature[y] = deepLandTemperature[y];

        for (int day = 0; day < daysPerYear; day++)
        {
            const float dayOfYear = static_cast<float>(day) + 0.5f;
            const auto solar = solarForcing(dayOfYear, world.tilt(), world.eccentricity(), world.perihelion());
            insolation[day][y] = dailyMeanInsolationWm2(
                static_cast<float>(latitude),
                solar.declinationRadians,
                solar.distanceFactor);
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
        int month = 0, monthEnd = climatehydrology::calendarMonth(0).days;
        for (int day = 0; day < daysPerYear; day++)
        {
            if (day >= monthEnd) monthEnd += climatehydrology::calendarMonth(++month).days;
            for (int y = 0; y <= height; y++)
            {
                const double latitude = latitudeForCell(y);
                const double polarFraction = std::pow(std::abs(latitude) / 90.0, 4.0);
                const double polarAdjustment = latitude >= 0.0 ?
                    static_cast<double>(world.northpolaradjust()) :
                    static_cast<double>(world.southpolaradjust());
                const double polarForcing = polarAdjustment * polarFraction *
                    (longwaveSlope + meridionalTransport + zonalExchange);

                for (int surface = 0; surface < 2; surface++)
                {
                    const bool land = surface == 0;
                    const double albedoTemperature = temperature[surface][y] -
                        (land ? meanLandElevationCooling[y] : 0.0);
                    const double albedo = effectiveAlbedo(
                        land,
                        albedoTemperature,
                        land ? rowPermanentIceFraction[y] : 0.0,
                        latitude);
                    absorbedSolar[surface][y] =
                        insolation[day][y] * (1.0 - albedo) + polarForcing;
                }
            }

            // Eliminate local surface/deep-layer exchange first, then solve one
            // tridiagonal system for the area-weighted row temperature. The
            // same conservative face heat flux heats both surface fractions.
            std::vector<double> rowSource(rows), rowResponse(rows);
            std::array<std::vector<double>, 2> denominator = {std::vector<double>(rows), std::vector<double>(rows)};
            std::array<std::vector<double>, 2> rhs = denominator;
            const double deepRate = deepLandCapacity / timeStepSeconds;
            for (int y = 0; y < rows; ++y)
            {
                double response = 0.0, source = 0.0;
                const double rowWeight = weights[0][y] + weights[1][y];
                for (int surface = 0; surface < 2; ++surface)
                {
                    const double rate = heatCapacity[surface] / timeStepSeconds;
                    denominator[surface][y] = rate + longwaveSlope + zonalExchange;
                    rhs[surface][y] = rate * temperature[surface][y] + absorbedSolar[surface][y] - longwaveInterceptWm2;
                    if (surface == 0 && deepLandCoupling > 0.0 && deepLandCapacity > 0.0)
                    {
                        denominator[surface][y] += deepLandCoupling * deepRate / (deepRate + deepLandCoupling);
                        rhs[surface][y] += deepLandCoupling * deepRate * deepLandTemperature[y] / (deepRate + deepLandCoupling);
                    }
                    const double fraction = weights[surface][y] / rowWeight;
                    response += fraction / denominator[surface][y];
                    source += fraction * rhs[surface][y] / denominator[surface][y];
                }
                const double local = 1.0 - zonalExchange * response;
                rowSource[y] = source / local;
                rowResponse[y] = response / local;
            }
            const auto diffusion = implicitMeridionalDiffusion(rowSource, rowResponse,
                tuning::climate::energybalance::meridionalDiffusivityWm2K);
            zonalMeanTemperature = diffusion.temperatureC;
            double globalMeanTemperature = 0.0;
            for (int y = 0; y < rows; ++y)
            {
                for (int surface = 0; surface < 2; ++surface)
                    nextTemperature[surface][y] = (rhs[surface][y] + zonalExchange * zonalMeanTemperature[y] +
                        diffusion.heatConvergenceWm2[y]) / denominator[surface][y];
                nextDeepLandTemperature[y] = (deepRate * deepLandTemperature[y] + deepLandCoupling * nextTemperature[0][y]) /
                    (deepRate + deepLandCoupling);
                globalMeanTemperature += (weights[0][y] + weights[1][y]) * zonalMeanTemperature[y] / totalAreaWeight;
            }

            if (year == spinupYears - 1)
            {
                for (int surface = 0; surface < 2; ++surface)
                    for (int y = 0; y < rows; ++y)
                        result.monthly[month][surface][y] += nextTemperature[surface][y] /
                            climatehydrology::calendarMonth(month).days;
                lastYearMeanTemperatureSum += globalMeanTemperature;

                for (int y = 0; y <= height; y++)
                {
                    for (int surface = 0; surface < 2; surface++)
                        lastYearSurfaceTemperatureSum[surface] += weights[surface][y] * nextTemperature[surface][y];
                }

                for (int season = 0; season < CLIMATESEASONCOUNT; season++)
                {
                    if (day + 1 == snapshotDays[season])
                    {
                        result.snapshots[season][0] = nextTemperature[0];
                        result.snapshots[season][1] = nextTemperature[1];
                    }
                }

                if (collectOutput)
                {
                    for (int y = 0; y <= height; y++)
                    {
                        for (int surface = 0; surface < 2; surface++)
                        {
                            const double weight = weights[surface][y];
                            const double transport =
                                zonalExchange * (zonalMeanTemperature[y] - nextTemperature[surface][y]) +
                                diffusion.heatConvergenceWm2[y];
                            const double storage = heatCapacity[surface] *
                                (nextTemperature[surface][y] - temperature[surface][y]) /
                                timeStepSeconds;
                            const double deepStorage = surface == 0 ?
                                deepLandCapacity *
                                    (nextDeepLandTemperature[y] - deepLandTemperature[y]) /
                                    timeStepSeconds :
                                0.0;
                            incomingSum += weight * insolation[day][y];
                            absorbedSum += weight * absorbedSolar[surface][y];
                            outgoingSum += weight *
                                (longwaveInterceptWm2 + longwaveSlope * nextTemperature[surface][y]);
                            transportSum += weight * transport;
                            storageSum += weight * (storage + deepStorage);
                        }
                    }
                }
            }

            temperature.swap(nextTemperature);
            deepLandTemperature.swap(nextDeepLandTemperature);
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
        result.budget.areaWeightedLandTemperatureC = surfaceAreaWeight[0] > 0.0
            ? lastYearSurfaceTemperatureSum[0] / (surfaceAreaWeight[0] * daysPerYear) : 0.0;
        result.budget.areaWeightedOceanTemperatureC = surfaceAreaWeight[1] > 0.0
            ? lastYearSurfaceTemperatureSum[1] / (surfaceAreaWeight[1] * daysPerYear) : 0.0;
    }
    const auto remapOutput = [&](std::vector<double>& field)
    {
        const auto raster = climategrid::remapSmoothField(
            1, rows, floats(field), 1, world.height() + 1);
        field.assign(raster.begin(), raster.end());
    };
    for (auto& season : result.snapshots)
        for (auto& surface : season)
            remapOutput(surface);
    for (auto& month : result.monthly)
        for (auto& surface : month) remapOutput(surface);
    return result;
}

CoupledProfileSimulation simulateCoupledProfiles(
    const planet& world,
    const std::vector<int>& landCells,
    const std::vector<int>& oceanCells,
    const std::vector<double>& localTemperatureOffset,
    const std::vector<double>& rowAlbedoTemperatureOffset,
    bool collectOutput)
{
    const int width = world.width();
    const int height = world.height();
    const int rowWidth = width + 1;
    const std::size_t cellCount =
        static_cast<std::size_t>(rowWidth) * static_cast<std::size_t>(height + 1);
    std::vector<double> rowPermanentIceFraction(height + 1, 0.0);
    std::vector<double> permanentIceFraction(cellCount, 0.0);
    std::array<std::vector<double>, CLIMATESEASONCOUNT> localTemperatureCorrection;

    for (auto& correction : localTemperatureCorrection)
        correction.assign(cellCount, 0.0);

    const double responseFeedback =
        tuning::climate::energybalance::longwaveSlopeWm2K +
        tuning::climate::energybalance::meridionalTransportWm2K +
        tuning::climate::energybalance::zonalLandOceanExchangeWm2K +
        tuning::climate::energybalance::landDeepLayerCouplingWm2K;
    const double annualAngularFrequency =
        2.0 * pi / (static_cast<double>(daysPerYear) * 86400.0);
    const double seasonalThermalResponse = std::hypot(
        responseFeedback,
        annualAngularFrequency * tuning::climate::energybalance::landHeatCapacityJm2K);
    const double targetMeanTemperatureC = static_cast<double>(world.averagetemp());
    double longwaveInterceptWm2 = tuning::climate::energybalance::earthLongwaveInterceptWm2;

    auto calculateLocalCorrections = [&](const ProfileSimulation& profile)
    {
        for (auto& correction : localTemperatureCorrection)
            std::fill(correction.begin(), correction.end(), 0.0);

        for (int y = 0; y <= height; y++)
        {
            if (landCells[y] == 0)
                continue;

            std::array<double, CLIMATESEASONCOUNT> rowCorrectionSum{};

            for (int x = 0; x <= width; x++)
            {
                if (world.sea(x, y) == 1)
                    continue;

                const std::size_t index =
                    static_cast<std::size_t>(y) * static_cast<std::size_t>(rowWidth) +
                    static_cast<std::size_t>(x);

                for (int season = 0; season < CLIMATESEASONCOUNT; season++)
                {
                    const float dayOfYear = static_cast<float>(snapshotDays[season]) - 0.5f;
                    const auto solar = solarForcing(dayOfYear, world.tilt(), world.eccentricity(), world.perihelion());
                    const double insolation = dailyMeanInsolationWm2(
                        static_cast<float>(latitudeDegreesForRow(y, height)),
                        solar.declinationRadians,
                        solar.distanceFactor);
                    const double rowTemperature =
                        profile.snapshots[season][0][y] + rowAlbedoTemperatureOffset[y];
                    const double localTemperature =
                        profile.snapshots[season][0][y] + localTemperatureOffset[index];
                    const double rowAlbedo = effectiveAlbedo(
                        true, rowTemperature, rowPermanentIceFraction[y], latitudeDegreesForRow(y, height));
                    const double localAlbedo = effectiveAlbedo(
                        true, localTemperature, permanentIceFraction[index], latitudeDegreesForRow(y, height));
                    const double correction =
                        insolation * (rowAlbedo - localAlbedo) / seasonalThermalResponse;
                    localTemperatureCorrection[season][index] = correction;
                    rowCorrectionSum[season] += correction;
                }
            }

            for (int season = 0; season < CLIMATESEASONCOUNT; season++)
            {
                const double rowMeanCorrection =
                    rowCorrectionSum[season] / static_cast<double>(landCells[y]);

                for (int x = 0; x <= width; x++)
                {
                    if (world.sea(x, y) == 1)
                        continue;

                    const std::size_t index =
                        static_cast<std::size_t>(y) * static_cast<std::size_t>(rowWidth) +
                        static_cast<std::size_t>(x);
                    localTemperatureCorrection[season][index] -= rowMeanCorrection;
                }
            }
        }
    };

    auto targetIceFractionForCell = [&](const ProfileSimulation& profile, int y, std::size_t index)
    {
        std::array<double, CLIMATESEASONCOUNT> temperatures;
        for (int season = 0; season < CLIMATESEASONCOUNT; season++)
            temperatures[season] = profile.snapshots[season][0][y];
        return detail::seasonalPermanentLandIceFraction(temperatures, localTemperatureOffset[index]);
    };

    ProfileSimulation profile;
    double maximumRowIceFractionResidual = 0.0;
    double temperatureErrorC = 0.0;
    int couplingIterations = 0;
    int stableIterations = 0;
    bool converged = false;

    for (int iteration = 0;
        iteration < tuning::climate::energybalance::localIceCouplingMaximumIterations;
        iteration++)
    {
        profile = simulateProfiles(
            world,
            landCells,
            oceanCells,
            rowPermanentIceFraction,
            longwaveInterceptWm2,
            false);
        maximumRowIceFractionResidual = 0.0;
        temperatureErrorC = profile.annualMeanTemperatureC - targetMeanTemperatureC;

        for (int y = 0; y <= height; y++)
        {
            if (landCells[y] == 0)
                continue;

            double targetRowIceFraction = 0.0;

            for (int x = 0; x <= width; x++)
            {
                if (world.sea(x, y) == 1)
                    continue;

                const std::size_t index =
                    static_cast<std::size_t>(y) * static_cast<std::size_t>(rowWidth) +
                    static_cast<std::size_t>(x);
                targetRowIceFraction += targetIceFractionForCell(profile, y, index);
            }

            targetRowIceFraction /= static_cast<double>(landCells[y]);
            const double iceFractionChange =
                std::abs(targetRowIceFraction - rowPermanentIceFraction[y]);

            if (iceFractionChange > maximumRowIceFractionResidual)
                maximumRowIceFractionResidual = iceFractionChange;

            rowPermanentIceFraction[y] +=
                tuning::climate::energybalance::localIceCouplingRelaxation *
                (targetRowIceFraction - rowPermanentIceFraction[y]);
        }

        if (world.climatesimulation().calibrateMeanTemperature) longwaveInterceptWm2 +=
            tuning::climate::energybalance::temperatureCalibrationRelaxation *
            tuning::climate::energybalance::longwaveSlopeWm2K *
            temperatureErrorC;
        couplingIterations = iteration + 1;

        if (maximumRowIceFractionResidual <=
                tuning::climate::energybalance::localIceCouplingTolerance &&
            (!world.climatesimulation().calibrateMeanTemperature || std::abs(temperatureErrorC) <=
                tuning::climate::energybalance::temperatureCalibrationToleranceC))
        {
            stableIterations++;

            if (stableIterations >=
                tuning::climate::energybalance::localIceCouplingStableIterations)
            {
                converged = true;
                break;
            }
        }
        else
            stableIterations = 0;
    }

    profile = simulateProfiles(
        world,
        landCells,
        oceanCells,
        rowPermanentIceFraction,
        longwaveInterceptWm2,
        collectOutput);
    temperatureErrorC = profile.annualMeanTemperatureC - targetMeanTemperatureC;
    maximumRowIceFractionResidual = 0.0;

    for (int y = 0; y <= height; y++)
    {
        if (landCells[y] == 0)
            continue;

        double targetRowIceFraction = 0.0;

        for (int x = 0; x <= width; x++)
        {
            if (world.sea(x, y) == 1)
                continue;

            const std::size_t index =
                static_cast<std::size_t>(y) * static_cast<std::size_t>(rowWidth) +
                static_cast<std::size_t>(x);
            permanentIceFraction[index] = targetIceFractionForCell(profile, y, index);
            targetRowIceFraction += permanentIceFraction[index];
        }

        targetRowIceFraction /= static_cast<double>(landCells[y]);
        maximumRowIceFractionResidual = std::max(
            maximumRowIceFractionResidual,
            std::abs(targetRowIceFraction - rowPermanentIceFraction[y]));
    }

    calculateLocalCorrections(profile);

    converged =
        maximumRowIceFractionResidual <= tuning::climate::energybalance::localIceCouplingTolerance &&
        (!world.climatesimulation().calibrateMeanTemperature || std::abs(temperatureErrorC) <=
            tuning::climate::energybalance::temperatureCalibrationToleranceC);

    CoupledProfileSimulation result;
    result.profile = std::move(profile);
    result.permanentIceFraction = std::move(permanentIceFraction);
    result.localTemperatureCorrection = std::move(localTemperatureCorrection);
    result.couplingIterations = couplingIterations;
    result.maximumRowIceFractionResidual = maximumRowIceFractionResidual;
    result.temperatureErrorC = temperatureErrorC;
    result.converged = converged;
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
    return solarForcing(dayOfYear, 0.0f, eccentricity, perihelionSeason).distanceFactor;
}

SolarForcing solarForcing(float dayOfYear, float obliquityDegrees,
    float eccentricity, int perihelionSeason)
{
    const double e = std::clamp(static_cast<double>(eccentricity), 0.0, 0.95);
    const double perihelionDay = perihelionSeason == 1 ? 186.0 : 3.0;
    const auto orbit = [&](double day)
    {
        const double mean = std::remainder(2.0 * pi * (day - perihelionDay) / daysPerYear, 2.0 * pi);
        double anomaly = mean;
        // Safeguarded Newton solve of E - e sin(E) = M.
        double lower = -pi, upper = pi;
        for (int iteration = 0; iteration < 48; ++iteration)
        {
            const double residual = anomaly - e * std::sin(anomaly) - mean;
            if (std::abs(residual) < 1.0e-13) break;
            if (residual > 0.0) upper = anomaly; else lower = anomaly;
            const double next = anomaly - residual / (1.0 - e * std::cos(anomaly));
            anomaly = next > lower && next < upper ? next : 0.5 * (lower + upper);
        }
        const double trueAnomaly = 2.0 * std::atan2(std::sqrt(1.0 + e) * std::sin(anomaly / 2.0),
            std::sqrt(1.0 - e) * std::cos(anomaly / 2.0));
        return std::pair<double, double>{trueAnomaly, 1.0 - e * std::cos(anomaly)};
    };
    const auto position = orbit(dayOfYear);
    // Keep the existing calendar convention: vernal equinox is day 80.
    const double longitude = position.first - orbit(80.0).first;
    return {static_cast<float>(std::asin(std::sin(obliquityDegrees * pi / 180.0) * std::sin(longitude))),
        static_cast<float>(1.0 / (position.second * position.second))};
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

float permanentLandIceFraction(float warmestSeasonTemperatureC)
{
    const double iceTransitionRange =
        tuning::climate::energybalance::permanentIceWarmTransitionC -
        tuning::climate::energybalance::permanentIceColdTransitionC;
    return static_cast<float>(smoothstep01(
        (tuning::climate::energybalance::permanentIceWarmTransitionC -
            static_cast<double>(warmestSeasonTemperatureC)) /
        iceTransitionRange));
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

MeridionalDiffusionResult implicitMeridionalDiffusion(const std::vector<double>& source,
    const std::vector<double>& response, double diffusivity)
{
    const int rows = static_cast<int>(source.size());
    if (rows == 0 || response.size() != source.size() || !std::isfinite(diffusivity) || diffusivity < 0.0) return {};
    for (int y = 0; y < rows; ++y)
        if (!std::isfinite(source[y]) || !std::isfinite(response[y]) || response[y] < 0.0) return {};
    std::vector<double> north(rows), south(rows), diagonal(rows), upper(rows), rhs = source;
    const double spacing = pi / rows;
    for (int y = 0; y < rows; ++y)
    {
        const double area = climategrid::latitudeBandMeasure(y, rows);
        north[y] = y > 0 ? diffusivity * std::sin(y * spacing) / (spacing * area) : 0.0;
        south[y] = y + 1 < rows ? diffusivity * std::sin((y + 1) * spacing) / (spacing * area) : 0.0;
        diagonal[y] = 1.0 + response[y] * (north[y] + south[y]);
        upper[y] = -response[y] * south[y];
        if (y > 0)
        {
            const double factor = -response[y] * north[y] / diagonal[y - 1];
            diagonal[y] -= factor * upper[y - 1];
            rhs[y] -= factor * rhs[y - 1];
        }
    }
    MeridionalDiffusionResult result;
    result.temperatureC.resize(rows); result.heatConvergenceWm2.resize(rows);
    for (int y = rows - 1; y >= 0; --y)
        result.temperatureC[y] = (rhs[y] - (y + 1 < rows ? upper[y] * result.temperatureC[y + 1] : 0.0)) / diagonal[y];
    for (int y = 0; y < rows; ++y)
        result.heatConvergenceWm2[y] =
            (y > 0 ? north[y] * (result.temperatureC[y - 1] - result.temperatureC[y]) : 0.0) +
            (y + 1 < rows ? south[y] * (result.temperatureC[y + 1] - result.temperatureC[y]) : 0.0);
    return result;
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

    const int rowWidth = width + 1;
    const double fractalScale = static_cast<double>(std::max(1, world.maxelevation()));
    std::vector<double> localTemperatureOffset(
        static_cast<std::size_t>(rowWidth) * static_cast<std::size_t>(height + 1),
        0.0);
    std::vector<double> rowAlbedoTemperatureOffset(height + 1, 0.0);

    for (int y = 0; y <= height; y++)
    {
        double landElevationCoolingSum = 0.0;

        for (int x = 0; x <= width; x++)
        {
            const std::size_t index =
                static_cast<std::size_t>(y) * static_cast<std::size_t>(rowWidth) +
                static_cast<std::size_t>(x);
            const double fractalAnomaly =
                (static_cast<double>(fractal[x][y]) / fractalScale - 0.5) *
                tuning::climate::energybalance::fractalTemperatureRangeC;
            const double elevationKm = static_cast<double>(
                std::max(0, world.map(x, y) - seaLevel)) / 1000.0;
            const double elevationCooling = world.sea(x, y) == 0 ?
                elevationKm * static_cast<double>(world.tempdecrease()) : 0.0;
            localTemperatureOffset[index] = fractalAnomaly - elevationCooling;

            if (world.sea(x, y) == 0)
                landElevationCoolingSum += elevationCooling;
        }

        if (landCells[y] > 0)
        {
            rowAlbedoTemperatureOffset[y] =
                -landElevationCoolingSum / static_cast<double>(landCells[y]);
        }
    }

    const CoupledProfileSimulation simulation = simulateCoupledProfiles(
        world,
        landCells,
        oceanCells,
        localTemperatureOffset,
        rowAlbedoTemperatureOffset,
        true);
    annualEnergyBudget = simulation.profile.budget;
    auto& monthly = world.climatesimulation().monthly;
    monthly = {};
    monthly.columns = width + 1; monthly.rows = height + 1;
    for (auto& field : monthly.temperatureC) field.resize(localTemperatureOffset.size());
    for (auto& field : monthly.referenceTemperatureC) field.resize(localTemperatureOffset.size());
    double landElevationCoolingSum = 0.0;
    double landAreaWeight = 0.0;
    double permanentLandIceAreaWeight = 0.0;
    double northernWarmestTemperatureSum = 0.0;
    double northernAreaWeight = 0.0;
    double southernWarmestTemperatureSum = 0.0;
    double southernIcecapAreaWeight = 0.0;
    double southernAreaWeight = 0.0;

    for (int y = 0; y <= height; y++)
    {
        for (int x = 0; x <= width; x++)
        {
            const std::size_t index =
                static_cast<std::size_t>(y) * static_cast<std::size_t>(rowWidth) +
                static_cast<std::size_t>(x);
            const int surface = world.sea(x, y) == 1 ? 1 : 0;
            const double elevationKm = static_cast<double>(std::max(0, world.map(x, y) - seaLevel)) / 1000.0;
            const double elevationCooling = surface == 0 ?
                elevationKm * static_cast<double>(world.tempdecrease()) : 0.0;

            if (surface == 0)
            {
                const double areaWeight = gridCellAreaWeight(y, height);
                landElevationCoolingSum += areaWeight * elevationCooling;
                landAreaWeight += areaWeight;
                permanentLandIceAreaWeight +=
                    areaWeight * simulation.permanentIceFraction[index];
            }

            for (int season = 0; season < CLIMATESEASONCOUNT; season++)
            {
                const double localCorrection = surface == 0 ?
                    simulation.localTemperatureCorrection[season][index] : 0.0;
                const double temperature = simulation.profile.snapshots[season][surface][y] +
                    localTemperatureOffset[index] + localCorrection;
                world.setseasonaltemp(season, x, y, static_cast<int>(std::round(temperature)));
                monthly.referenceTemperatureC[season][index] = static_cast<float>(world.seasonaltemp(season, x, y));
            }

            for (int month = 0; month < 12; ++month)
            {
                const auto calendar = climatehydrology::calendarMonth(month);
                const double correction = surface == 0 ? climatehydrology::interpolateSeasonal(
                    static_cast<float>(simulation.localTemperatureCorrection[calendar.firstSeason][index]),
                    static_cast<float>(simulation.localTemperatureCorrection[calendar.secondSeason][index]), calendar.interpolation) : 0.0;
                monthly.temperatureC[month][index] = static_cast<float>(simulation.profile.monthly[month][surface][y] +
                    localTemperatureOffset[index] + correction);
            }

            world.setjantemp(x, y, world.seasonaltemp(seasonjanuary, x, y));
            world.setjultemp(x, y, world.seasonaltemp(seasonjuly, x, y));

            if (surface == 0)
            {
                int warmestTemperature = world.seasonaltemp(0, x, y);

                for (int season = 1; season < CLIMATESEASONCOUNT; season++)
                {
                    warmestTemperature = std::max(
                        warmestTemperature,
                        world.seasonaltemp(season, x, y));
                }

                const double latitude = latitudeDegreesForRow(y, height);
                const double areaWeight = gridCellAreaWeight(y, height);

                if (latitude >= 50.0 && latitude <= 70.0)
                {
                    northernWarmestTemperatureSum +=
                        areaWeight * static_cast<double>(warmestTemperature);
                    northernAreaWeight += areaWeight;
                }

                if (latitude <= -60.0)
                {
                    southernWarmestTemperatureSum +=
                        areaWeight * static_cast<double>(warmestTemperature);
                    southernIcecapAreaWeight += warmestTemperature < 0 ? areaWeight : 0.0;
                    southernAreaWeight += areaWeight;
                }
            }
        }
    }

    auto rowAnnualMean = [&](int y)
    {
        double total = 0.0;

        for (int season = 0; season < CLIMATESEASONCOUNT; season++)
        {
            const double cells = static_cast<double>(landCells[y] + oceanCells[y]);
            total += cells > 0.0 ?
                (simulation.profile.snapshots[season][0][y] * static_cast<double>(landCells[y]) +
                    simulation.profile.snapshots[season][1][y] * static_cast<double>(oceanCells[y])) / cells :
                static_cast<double>(world.averagetemp());
        }

        return total / static_cast<double>(CLIMATESEASONCOUNT);
    };

    world.setnorthpolartemp(static_cast<int>(std::round(rowAnnualMean(0))));
    world.setsouthpolartemp(static_cast<int>(std::round(rowAnnualMean(height))));
    world.seteqtemp(static_cast<int>(std::round(rowAnnualMean(height / 2))));
    annualEnergyBudget.areaWeightedLandElevationCoolingC =
        landAreaWeight > 0.0 ? landElevationCoolingSum / landAreaWeight : 0.0;
    annualEnergyBudget.areaWeightedPermanentLandIceFraction =
        landAreaWeight > 0.0 ? permanentLandIceAreaWeight / landAreaWeight : 0.0;
    annualEnergyBudget.localIceCouplingIterations = simulation.couplingIterations;
    annualEnergyBudget.localIceCouplingMaximumRowResidual =
        simulation.maximumRowIceFractionResidual;
    annualEnergyBudget.localIceCouplingTemperatureErrorC = simulation.temperatureErrorC;
    annualEnergyBudget.localIceCouplingConverged = simulation.converged;
    annualEnergyBudget.northern5070WarmestTemperatureC =
        northernAreaWeight > 0.0 ?
            northernWarmestTemperatureSum / northernAreaWeight : 0.0;
    annualEnergyBudget.southern6090WarmestTemperatureC =
        southernAreaWeight > 0.0 ?
            southernWarmestTemperatureSum / southernAreaWeight : 0.0;
    annualEnergyBudget.southern6090IcecapThermalFraction =
        southernAreaWeight > 0.0 ?
            southernIcecapAreaWeight / southernAreaWeight : 0.0;
    std::cout
        << "Climate energy coupling iterations=" << simulation.couplingIterations
        << " maximum_row_ice_fraction_residual="
        << simulation.maximumRowIceFractionResidual
        << " temperature_error_c=" << simulation.temperatureErrorC
        << " converged=" << (simulation.converged ? 1 : 0)
        << '\n';
    std::cout
        << "Climate thermal gates northern_50_70_warmest_temperature_c="
        << annualEnergyBudget.northern5070WarmestTemperatureC
        << " southern_60_90_warmest_temperature_c="
        << annualEnergyBudget.southern6090WarmestTemperatureC
        << " southern_60_90_icecap_thermal_fraction="
        << annualEnergyBudget.southern6090IcecapThermalFraction
        << '\n';
}

const AnnualEnergyBudget& lastAnnualEnergyBudget()
{
    return annualEnergyBudget;
}
}
