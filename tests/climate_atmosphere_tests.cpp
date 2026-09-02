#include "climate_atmosphere.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace
{
int failures = 0;

void expect(bool condition, const char* message)
{
    if (condition)
        return;

    std::cerr << "FAIL: " << message << '\n';
    failures++;
}
}

int main()
{
    constexpr float pi = 3.14159265358979323846f;
    constexpr float earthRotation = 7.2921159e-5f;
    expect(
        std::abs(climateatmosphere::coriolisParameterPerSecond(0.0f, earthRotation)) < 1.0e-9f,
        "Coriolis acceleration must vanish at the equator");
    expect(
        climateatmosphere::coriolisParameterPerSecond(45.0f, earthRotation) > 0.0f &&
            climateatmosphere::coriolisParameterPerSecond(-45.0f, earthRotation) < 0.0f,
        "Coriolis sign must reverse across the equator");

    const float heightResponse = climateatmosphere::hypsometricHeightResponseMetresPerKelvin(
        100000.0f, 50000.0f);
    expect(
        std::abs(heightResponse - 20.27f) < 0.05f,
        "1000-to-500 hPa thickness must change by about 20.27 metres per kelvin");

    const float earthHadleyEdge = climateatmosphere::heldHouHadleyEdgeLatitudeDegrees(
        60.0f, 10000.0f, 288.0f, 9.80665f, earthRotation, 6371000.0f);
    expect(
        earthHadleyEdge > 20.0f && earthHadleyEdge < 26.0f,
        "Held-Hou scaling must place the Earth-like Hadley edge in the subtropics");
    expect(
        climateatmosphere::heldHouHadleyEdgeLatitudeDegrees(
            60.0f, 10000.0f, 288.0f, 9.80665f, 0.0f, 6371000.0f) == 90.0f,
        "a non-rotating atmosphere must permit pole-to-pole overturning");

    const float warmSurfacePressure = climateatmosphere::thermalSurfacePressureAnomalyHpa(
        10.0f, 1000.0f, 288.0f, 0.12f);
    const float coldSurfacePressure = climateatmosphere::thermalSurfacePressureAnomalyHpa(
        -10.0f, 1000.0f, 288.0f, 0.12f);
    expect(
        warmSurfacePressure < 0.0f && coldSurfacePressure > 0.0f &&
            std::abs(warmSurfacePressure + coldSurfacePressure) < 0.0001f,
        "warm columns must form thermal lows and cold columns thermal highs");

    const float warmModePressure = climateatmosphere::thermalModePressureAnomalyHpa(
        10.0f, 1.225f, 100000.0f, 50000.0f);
    expect(
        std::abs(warmModePressure + 24.37f) < 0.05f,
        "the thermal-mode pressure must equal the hypsometric geopotential anomaly");

    const float equatorialPressure =
        climateatmosphere::axisymmetricOverturningPressureAnomalyHpa(0.0f, 0.0f, 25.0f, 10.0f);
    const float subtropicalPressure =
        climateatmosphere::axisymmetricOverturningPressureAnomalyHpa(25.0f, 0.0f, 25.0f, 10.0f);
    const float subpolarPressure =
        climateatmosphere::axisymmetricOverturningPressureAnomalyHpa(57.5f, 0.0f, 25.0f, 10.0f);
    const float polarPressure =
        climateatmosphere::axisymmetricOverturningPressureAnomalyHpa(90.0f, 0.0f, 25.0f, 10.0f);
    expect(
        equatorialPressure < 0.0f && subtropicalPressure > 0.0f &&
            subpolarPressure < 0.0f && polarPressure > 0.0f,
        "axisymmetric overturning must produce equatorial and subpolar lows with subtropical and polar highs");

    constexpr int responseColumns = 72;
    constexpr int responseRows = 37;
    const auto responseIndex = [=](int x, int y)
    {
        return static_cast<size_t>(y * responseColumns + x);
    };
    std::vector<float> thermalImpulse(
        static_cast<size_t>(responseColumns * responseRows),
        0.0f);
    const int impulseColumn = responseColumns / 2;
    const int equatorRow = responseRows / 2;
    thermalImpulse[responseIndex(impulseColumn, equatorRow)] = 1.0f;
    const auto progradeResponse = climateatmosphere::nonlocalThermalResponse(
        responseColumns,
        responseRows,
        thermalImpulse,
        20.0f,
        20.0f,
        8.0f,
        1.0f);
    const auto retrogradeResponse = climateatmosphere::nonlocalThermalResponse(
        responseColumns,
        responseRows,
        thermalImpulse,
        20.0f,
        20.0f,
        8.0f,
        -1.0f);
    expect(
        progradeResponse[responseIndex(impulseColumn - 1, equatorRow)] >
            progradeResponse[responseIndex(impulseColumn + 1, equatorRow)] &&
        retrogradeResponse[responseIndex(impulseColumn + 1, equatorRow)] >
            retrogradeResponse[responseIndex(impulseColumn - 1, equatorRow)],
        "the tropical response must extend westward relative to planetary rotation");
    const int extratropicalRow = 10;
    thermalImpulse.assign(static_cast<size_t>(responseColumns * responseRows), 0.0f);
    thermalImpulse[responseIndex(impulseColumn, extratropicalRow)] = 1.0f;
    const auto extratropicalResponse = climateatmosphere::nonlocalThermalResponse(
        responseColumns,
        responseRows,
        thermalImpulse,
        20.0f,
        20.0f,
        8.0f,
        1.0f);
    expect(
        std::abs(
            extratropicalResponse[responseIndex(impulseColumn - 1, extratropicalRow)] -
            extratropicalResponse[responseIndex(impulseColumn + 1, extratropicalRow)]) < 1.0e-6f,
        "the extratropical thermal response must be zonally symmetric");
    double equatorialResponseSum = 0.0;
    for (int x = 0; x < responseColumns; x++)
        equatorialResponseSum += progradeResponse[responseIndex(x, equatorRow)];
    expect(
        equatorialResponseSum > 0.0 && equatorialResponseSum < 1.0,
        "the nonlocal response must redistribute an impulse without creating a new extremum");

    constexpr int topographicColumns = 21;
    constexpr int topographicRows = 9;
    const auto topographicIndex = [=](int x, int y)
    {
        return static_cast<size_t>(y * topographicColumns + x);
    };
    std::vector<float> ridgeTerrain(
        static_cast<size_t>(topographicColumns * topographicRows),
        0.0f);
    std::vector<float> ridgeEastWind(ridgeTerrain.size(), 10.0f);
    std::vector<float> ridgeSouthWind(ridgeTerrain.size(), 0.0f);
    for (int y = 0; y < topographicRows; y++)
        ridgeTerrain[topographicIndex(10, y)] = 2000.0f;
    const auto ridgeForcing = climateatmosphere::mechanicalTopographicPressureForcingHpa(
        topographicColumns,
        topographicRows,
        ridgeTerrain,
        ridgeEastWind,
        ridgeSouthWind,
        3.0f,
        2000.0f,
        2.0f,
        1.0f,
        8.0f,
        10.0f,
        30.0f);
    expect(
        ridgeForcing[topographicIndex(7, 2)] > 0.0f &&
            ridgeForcing[topographicIndex(13, 2)] < 0.0f,
        "westerly flow across a ridge must produce a windward high and lee trough");
    ridgeEastWind.assign(ridgeTerrain.size(), -10.0f);
    const auto reversedRidgeForcing =
        climateatmosphere::mechanicalTopographicPressureForcingHpa(
            topographicColumns,
            topographicRows,
            ridgeTerrain,
            ridgeEastWind,
            ridgeSouthWind,
            3.0f,
            2000.0f,
            2.0f,
            1.0f,
            8.0f,
            10.0f,
            30.0f);
    expect(
        reversedRidgeForcing[topographicIndex(13, 2)] > 0.0f &&
            reversedRidgeForcing[topographicIndex(7, 2)] < 0.0f,
        "topographic pressure forcing must reverse when the background wind reverses");

    const auto spacing = climateatmosphere::cellSpacingMetres(0.0f, 2048, 1025, 6371000.0f);
    expect(
        spacing.zonalMetres > 19000.0f && spacing.zonalMetres < 20000.0f &&
            spacing.meridionalMetres > 19000.0f && spacing.meridionalMetres < 20000.0f,
        "Earth benchmark cells must be about 19.5 kilometres at the equator");

    const auto unrotated = climateatmosphere::steadyRayleighCoriolisWind(
        0.001f, 0.002f, 0.0f, 1000.0f, earthRotation);
    expect(
        std::abs(unrotated.eastMetresPerSecond - 1.0f) < 0.001f &&
            std::abs(unrotated.southMetresPerSecond + 2.0f) < 0.001f,
        "Rayleigh flow at the equator must follow the pressure-gradient acceleration");

    const auto northern = climateatmosphere::steadyRayleighCoriolisWind(
        0.0f, 0.001f, 45.0f, 86400.0f, earthRotation);
    const auto southern = climateatmosphere::steadyRayleighCoriolisWind(
        0.0f, 0.001f, -45.0f, 86400.0f, earthRotation);
    expect(
        northern.eastMetresPerSecond > 0.0f && southern.eastMetresPerSecond < 0.0f,
        "the same meridional height force must produce opposite zonal flow across the equator");

    const auto quadraticEquatorial = climateatmosphere::steadyQuadraticDragCoriolisWind(
        0.001f, 0.0f, 0.0f, 0.001f, 100.0f, earthRotation);
    const auto quadraticEquatorialStronger = climateatmosphere::steadyQuadraticDragCoriolisWind(
        0.004f, 0.0f, 0.0f, 0.001f, 100.0f, earthRotation);
    expect(
        std::abs(quadraticEquatorial.eastMetresPerSecond - 10.0f) < 0.001f &&
            std::abs(quadraticEquatorial.southMetresPerSecond) < 0.001f &&
            std::abs(quadraticEquatorialStronger.eastMetresPerSecond - 20.0f) < 0.001f,
        "quadratic drag must keep equatorial flow finite and scale speed with the square root of force");

    const auto quadraticNorthern = climateatmosphere::steadyQuadraticDragCoriolisWind(
        0.0f, 0.001f, 45.0f, 0.0013f, 300.0f, earthRotation);
    const auto quadraticSouthern = climateatmosphere::steadyQuadraticDragCoriolisWind(
        0.0f, 0.001f, -45.0f, 0.0013f, 300.0f, earthRotation);
    expect(
        quadraticNorthern.eastMetresPerSecond > 0.0f &&
            quadraticSouthern.eastMetresPerSecond < 0.0f,
        "quadratic surface drag must preserve the Coriolis reversal across the equator");

    constexpr int waveColumns = 64;
    constexpr int waveRows = waveColumns / 2;
    constexpr int waveCentreX = 16;
    constexpr int waveCentreY = waveRows / 2;
    const auto waveIndex = [=](int x, int y)
    {
        return static_cast<size_t>(y) * waveColumns + x;
    };
    std::vector<float> waveForcing(waveColumns * waveRows, 0.0f);
    std::vector<float> waveDragTime(waveForcing.size(), 43200.0f);
    for (int y = 0; y < waveRows; y++)
    {
        for (int x = 0; x < waveColumns; x++)
        {
            int distanceX = x - waveCentreX;
            if (distanceX > waveColumns / 2)
                distanceX -= waveColumns;
            if (distanceX < -waveColumns / 2)
                distanceX += waveColumns;
            const int distanceY = y - waveCentreY;
            waveForcing[waveIndex(x, y)] = std::exp(
                -static_cast<float>(distanceX * distanceX) / 25.0f -
                static_cast<float>(distanceY * distanceY) / 9.0f);
        }
    }
    const auto progradeWave = climateatmosphere::solveSteadyStationaryWavePressure(
        waveColumns,
        waveRows,
        waveForcing,
        waveDragTime,
        48.0f,
        2.2f * 86400.0f,
        1.225f,
        6371000.0f,
        earthRotation,
        1.0f,
        true,
        200,
        1.0e-4f);
    const auto retrogradeWave = climateatmosphere::solveSteadyStationaryWavePressure(
        waveColumns,
        waveRows,
        waveForcing,
        waveDragTime,
        48.0f,
        2.2f * 86400.0f,
        1.225f,
        6371000.0f,
        earthRotation,
        -1.0f,
        true,
        200,
        1.0e-4f);
    const auto fullPressureWave = climateatmosphere::solveSteadyStationaryWavePressure(
        waveColumns,
        waveRows,
        waveForcing,
        waveDragTime,
        48.0f,
        2.2f * 86400.0f,
        1.225f,
        6371000.0f,
        earthRotation,
        1.0f,
        false,
        300,
        1.0e-4f);
    double waveMagnitude = 0.0;
    double waveAsymmetry = 0.0;
    double reversedWaveDifference = 0.0;
    double maximumRowMean = 0.0;
    double fullPressureAreaTotal = 0.0;
    double fullPressureAreaWeight = 0.0;
    double fullPressureRowMeanMagnitude = 0.0;
    for (int y = 0; y < waveRows; y++)
    {
        double rowMean = 0.0;
        double fullPressureRowMean = 0.0;
        for (int x = 0; x < waveColumns; x++)
        {
            const int mirroredX = (2 * waveCentreX - x + waveColumns) % waveColumns;
            const float value = progradeWave.pressureAnomalyHpa[waveIndex(x, y)];
            waveMagnitude += std::abs(value);
            waveAsymmetry += std::abs(
                value - progradeWave.pressureAnomalyHpa[waveIndex(mirroredX, y)]);
            reversedWaveDifference += std::abs(
                value - retrogradeWave.pressureAnomalyHpa[waveIndex(mirroredX, y)]);
            rowMean += value;
            fullPressureRowMean +=
                fullPressureWave.pressureAnomalyHpa[waveIndex(x, y)];
        }
        maximumRowMean = std::max(
            maximumRowMean,
            std::abs(rowMean / static_cast<double>(waveColumns)));
        fullPressureRowMean /= static_cast<double>(waveColumns);
        fullPressureRowMeanMagnitude += std::abs(fullPressureRowMean);
        const double latitudeRadians =
            (90.0 - 180.0 * (static_cast<double>(y) + 0.5) /
                    static_cast<double>(waveRows)) *
            3.14159265358979323846 / 180.0;
        const double areaWeight = std::max(0.0, std::cos(latitudeRadians));
        fullPressureAreaTotal += fullPressureRowMean * areaWeight;
        fullPressureAreaWeight += areaWeight;
    }
    if (!progradeWave.converged || !retrogradeWave.converged)
    {
        std::cerr
            << "stationary-wave diagnostics prograde_iterations=" << progradeWave.iterations
            << " prograde_residual=" << progradeWave.relativeResidual
            << " retrograde_iterations=" << retrogradeWave.iterations
            << " retrograde_residual=" << retrogradeWave.relativeResidual << '\n';
    }
    expect(
        progradeWave.converged && retrogradeWave.converged &&
            progradeWave.relativeResidual < 1.0e-4f &&
            retrogradeWave.relativeResidual < 1.0e-4f,
        "stationary-wave pressure solves must converge to their requested residual");
    expect(
        waveMagnitude > 0.0 && waveAsymmetry / waveMagnitude > 0.01 &&
            reversedWaveDifference / waveMagnitude < 0.001,
        "rotation must create a longitudinally asymmetric response that mirrors when rotation reverses");
    expect(
        maximumRowMean < 1.0e-5,
        "the stationary-wave response must preserve each row's zonal mean");
    expect(
        fullPressureWave.converged && fullPressureRowMeanMagnitude > 0.01 &&
            std::abs(fullPressureAreaTotal / fullPressureAreaWeight) < 1.0e-5,
        "the full pressure solve must retain zonal structure while conserving global mean pressure");
    expect(!progradeWave.residualHistory.empty() && progradeWave.restartCycles > 0 &&
            progradeWave.residualHistory.back() <= progradeWave.residualHistory.front(),
        "stationary solves must retain restart-cycle physical residual histories");

    constexpr int modeColumns = 16;
    constexpr int modeRows = 8;
    constexpr std::size_t modeCellCount = modeColumns * modeRows;
    std::vector<float> absorbed(modeCellCount, 220.0f);
    std::vector<float> outgoing(modeCellCount, 220.0f);
    std::vector<float> sensible(modeCellCount, 0.0f);
    std::vector<float> condensation(modeCellCount, 0.0f);
    for (int y = 0; y < modeRows; y++)
    {
        for (int x = 0; x < modeColumns; x++)
        {
            const std::size_t cell = static_cast<std::size_t>(y) * modeColumns + x;
            absorbed[cell] += 20.0f * std::cos(2.0f * pi * x / modeColumns);
            sensible[cell] = 4.0f * std::sin(2.0f * pi * x / modeColumns);
            condensation[cell] = x < modeColumns / 2 ? 10.0f : 0.0f;
        }
    }
    const auto heating = climateatmosphere::diagnoseDiabaticHeating(
        modeColumns,
        modeRows,
        absorbed,
        outgoing,
        sensible,
        condensation,
        86400.0f,
        0.35f,
        1.0f);
    expect(heating.areaWeightedLatentHeatingWm2 > 140.0 &&
            heating.areaWeightedLatentHeatingWm2 < 150.0,
        "latent heating must convert millimetres per day to watts per square metre");
    expect(heating.maximumAbsoluteRowMeanWm2 < 1.0e-5f,
        "stationary heating projection must have a controlled zero zonal mean");
    const auto noLatentHeating = climateatmosphere::diagnoseDiabaticHeating(
        modeColumns,
        modeRows,
        absorbed,
        outgoing,
        sensible,
        condensation,
        86400.0f,
        0.35f,
        0.0f);
    expect(std::abs(noLatentHeating.areaWeightedLatentHeatingWm2) < 1.0e-9,
        "the latent projection switch must prevent double-counting hydrology energy");

    const auto earthParameters = climateatmosphere::diagnoseStationaryParameters(
        0.01f, 10000.0f, 9.80665f, 6371000.0f, earthRotation,
        64, 32, 0.10f, 1200000.0f);
    const auto slowRotationParameters = climateatmosphere::diagnoseStationaryParameters(
        0.01f, 10000.0f, 9.80665f, 6371000.0f, earthRotation * 0.5f,
        64, 32, 0.10f, 1200000.0f);
    expect(earthParameters.equivalentDepthMetres > 0.0f &&
            earthParameters.maximumZonalWavenumber <= modeColumns * 2 &&
            slowRotationParameters.adjustmentLengthMetres >
                earthParameters.adjustmentLengthMetres,
        "stationary parameters must derive from stratification, rotation, and resolution");

    std::vector<float> zonalPressure(modeRows, 0.0f);
    std::vector<float> orographic(modeCellCount, 0.0f);
    for (int y = 0; y < modeRows; y++)
        zonalPressure[y] = 3.0f * std::cos(pi * (y + 0.5f) / modeRows);
    climateatmosphere::ModeSeparatedCirculationConfig modeConfig;
    modeConfig.maximumIterations = 500;
    const auto separated = climateatmosphere::solveModeSeparatedCirculation(
        modeColumns,
        modeRows,
        zonalPressure,
        heating.stationaryProjectedHeatingWm2,
        orographic,
        modeConfig);
    double surfaceMagnitude = 0.0;
    double upperMagnitude = 0.0;
    double maximumZonalTransfer = 0.0;
    for (int y = 0; y < modeRows; y++)
    {
        double rowMean = 0.0;
        for (int x = 0; x < modeColumns; x++)
        {
            const std::size_t cell = static_cast<std::size_t>(y) * modeColumns + x;
            rowMean += separated.surfacePressureAnomalyHpa[cell];
            surfaceMagnitude += std::abs(separated.surfaceEastWindMps[cell]) +
                std::abs(separated.surfaceSouthWindMps[cell]);
            upperMagnitude += std::abs(separated.upperEastWindMps[cell]) +
                std::abs(separated.upperSouthWindMps[cell]);
        }
        maximumZonalTransfer = std::max(
            maximumZonalTransfer,
            std::abs(rowMean / modeColumns - zonalPressure[y]));
    }
    expect(separated.surfaceStationarySolver.converged &&
            separated.upperStationarySolver.converged &&
            maximumZonalTransfer < 1.0e-4,
        "surface and upper stationary modes must converge without collapsing the zonal mode");
    expect(surfaceMagnitude > 0.0 && upperMagnitude > 0.0,
        "separately closed surface and upper modes must both respond to heating");
    modeConfig.enabled.stationary = false;
    modeConfig.enabled.upper = false;
    const auto zonalOnly = climateatmosphere::solveModeSeparatedCirculation(
        modeColumns,
        modeRows,
        zonalPressure,
        heating.stationaryProjectedHeatingWm2,
        orographic,
        modeConfig);
    expect(std::all_of(
            zonalOnly.upperHeightAnomalyMetres.begin(),
            zonalOnly.upperHeightAnomalyMetres.end(),
            [](float value) { return value == 0.0f; }),
        "upper-only and stationary-only responses must remain isolatable");

    const std::vector<float> zeroHeating(modeCellCount, 0.0f);
    const std::vector<float> zeroZonal(modeRows, 0.0f);
    modeConfig = {};
    const auto zeroResponse = climateatmosphere::solveModeSeparatedCirculation(modeColumns, modeRows,
        zeroZonal, zeroHeating, zeroHeating, modeConfig);
    expect(zeroResponse.surfaceStationarySolver.converged && zeroResponse.upperStationarySolver.converged &&
        zeroResponse.areaWeightedKineticEnergyJm2 == 0.0, "zero forcing must produce zero flow without noise");
    for (int mask = 0; mask < 16; ++mask)
    {
        modeConfig.enabled = {(mask & 1) != 0, (mask & 2) != 0, (mask & 4) != 0, (mask & 8) != 0};
        modeConfig.zonalUpperHeightMetres = zonalPressure;
        const auto isolated = climateatmosphere::solveModeSeparatedCirculation(modeColumns, modeRows,
            zonalPressure, heating.stationaryProjectedHeatingWm2, orographic, modeConfig);
        const auto zero = [](const auto& values) { return std::all_of(values.begin(), values.end(), [](float v) { return v == 0.0f; }); };
        expect((modeConfig.enabled.surface || (zero(isolated.surfaceEastWindMps) && zero(isolated.surfacePressureAnomalyHpa))) &&
            (modeConfig.enabled.upper || (zero(isolated.upperEastWindMps) && zero(isolated.upperHeightAnomalyMetres))),
            "all 16 mode-isolation combinations must respect disabled layers");
    }
    modeConfig = {};
    modeConfig.maximumIterations = 1;
    modeConfig.relativeTolerance = 1.0e-12f;
    const auto failed = climateatmosphere::solveModeSeparatedCirculation(modeColumns, modeRows,
        zonalPressure, heating.stationaryProjectedHeatingWm2, orographic, modeConfig);
    expect(!failed.surfaceStationarySolver.converged && failed.surfacePressureAnomalyHpa[0] == zonalPressure[0],
        "a failed stationary solve must fall back to the independently closed zonal mode");
    std::vector<float> mountain(modeCellCount), jets(modeCellCount, 20.0f), calm(modeCellCount, 0.0f);
    for (int y = 0; y < modeRows; ++y)
        for (int x = 0; x < modeColumns; ++x)
            mountain[y * modeColumns + x] = 200.0f * std::cos(2.0f * pi * x / modeColumns);
    modeConfig = {};
    const auto lowerMountain = climateatmosphere::upperOrographicHeightForcing(modeColumns, modeRows,
        mountain, jets, 0.01f, 3000.0f, 86400.0f, modeConfig);
    const auto upperMountain = climateatmosphere::upperOrographicHeightForcing(modeColumns, modeRows,
        mountain, jets, 0.01f, 8000.0f, 86400.0f, modeConfig);
    const auto calmMountain = climateatmosphere::upperOrographicHeightForcing(modeColumns, modeRows,
        mountain, calm, 0.01f, 5000.0f, 86400.0f, modeConfig);
    expect(lowerMountain != upperMountain && std::all_of(upperMountain.begin(), upperMountain.end(),
        [](float v) { return std::isfinite(v); }) && std::all_of(calmMountain.begin(), calmMountain.end(),
        [](float v) { return v == 0.0f; }), "mountain waves must propagate/damp vertically and vanish without incident wind");
    expect(climateatmosphere::diagnoseBruntVaisalaFrequency(280.0f, 0.006f, 9.80665f) > 0.0f &&
        climateatmosphere::diagnoseBruntVaisalaFrequency(280.0f, 0.012f, 9.80665f) == 0.0f,
        "stratification must distinguish stable and convectively unstable lapse rates");
    const auto refined = climateatmosphere::diagnoseStationaryParameters(0.01f, 10000.0f, 9.80665f,
        6371000.0f, earthRotation, 128, 64, 0.1f, 1200000.0f);
    const auto twiceRefined = climateatmosphere::diagnoseStationaryParameters(0.01f, 10000.0f, 9.80665f,
        6371000.0f, earthRotation, 256, 128, 0.1f, 1200000.0f);
    expect(refined.maximumZonalWavenumber == twiceRefined.maximumZonalWavenumber,
        "refinement must not invent physical forcing bandwidth");
    const auto largeRestart = climateatmosphere::solveSteadyStationaryWavePressure(
        modeColumns, modeRows, heating.stationaryProjectedHeatingWm2, std::vector<float>(modeCellCount, 43200.0f),
        48.0f, 190080.0f, 1.225f, 6371000.0f, earthRotation, 1.0f, true, 500, 1.0e-4f, 128);
    expect(largeRestart.converged && largeRestart.relativeResidual <= 1.0e-4f,
        "configurable larger GMRES restart windows must retain physical residual acceptance");
    modeConfig = {};
    modeConfig.enabled = {true, false, true, false};
    climateatmosphere::ColumnHeatingInput column;
    column.incomingSolarWm2 = 340.0;
    column.sensibleHeatingWm2 = 20.0;
    column.condensationMm = {1.0, 4.0};
    column.reevaporationMm = 2.0;
    column.surfaceEvaporationMm = 3.0;
    const auto heatColumn = climateatmosphere::diagnoseColumnHeating(column);
    expect(std::abs(heatColumn.closureResidualWm2) < 1.0e-10 &&
        heatColumn.latentWm2[0] < 0.0 && heatColumn.latentWm2[1] > 0.0,
        "grey radiation and phase changes must close energy and locate re-evaporative cooling below condensation");
    column.longwaveOpticalDepth = {0.0, 0.0};
    column.shortwaveOpticalDepth = {0.0, 0.0};
    const auto transparent = climateatmosphere::diagnoseColumnHeating(column);
    expect(transparent.radiativeWm2[0] == 0.0 && transparent.radiativeWm2[1] == 0.0,
        "transparent air must not receive surface radiative heating");
    const auto uniformDrag = climateatmosphere::solveModeSeparatedCirculation(modeColumns, modeRows,
        zonalPressure, zeroHeating, zeroHeating, modeConfig);
    modeConfig.surfaceDragCoefficients.assign(modeCellCount, modeConfig.surfaceDragCoefficient);
    for (int y = 0; y < modeRows; ++y)
        for (int x = modeColumns / 2; x < modeColumns; ++x)
            modeConfig.surfaceDragCoefficients[y * modeColumns + x] *= 8.0f;
    const auto roughLand = climateatmosphere::solveModeSeparatedCirculation(modeColumns, modeRows,
        zonalPressure, zeroHeating, zeroHeating, modeConfig);
    const int equatorialRow = modeRows / 2;
    const int oceanCell = equatorialRow * modeColumns;
    const int landCell = oceanCell + modeColumns / 2;
    const auto speed = [](const auto& flow, int cell) {
        return std::hypot(flow.surfaceEastWindMps[cell], flow.surfaceSouthWindMps[cell]); };
    expect(speed(roughLand, landCell) < speed(roughLand, oceanCell) &&
        speed(roughLand, oceanCell) == speed(uniformDrag, oceanCell),
        "the final quadratic wind solver must use local drag without changing ocean cells");
    modeConfig = {};
    modeConfig.interlayerMomentumCoupling = 0.0f;
    modeConfig.upperMaximumZonalWavenumber = 0;
    const auto upperFiltered = climateatmosphere::solveModeSeparatedCirculation(modeColumns, modeRows,
        zonalPressure, heating.stationaryProjectedHeatingWm2, orographic, modeConfig);
    expect(upperFiltered.surfacePressureAnomalyHpa == separated.surfacePressureAnomalyHpa &&
        std::all_of(upperFiltered.upperHeightAnomalyMetres.begin(), upperFiltered.upperHeightAnomalyMetres.end(),
            [](float v) { return std::abs(v) < 1.0e-4f; }),
        "upper bandwidth must filter upper forcing independently of surface pressure");
    return failures == 0 ? 0 : 1;
}
