#include "climate_weather.hpp"
#include "climate_grid.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <vector>

namespace
{
int failures = 0;

void expect(bool condition, const char* message)
{
    if (!condition)
    {
        std::cerr << "FAIL: " << message << '\n';
        failures++;
    }
}
}

int main()
{
    constexpr int columns = 16;
    constexpr int rows = 8;
    constexpr std::size_t cellCount = columns * rows;
    climateweather::ShallowWaterConfig config;
    config.layerCount = 2;
    auto resting = climateweather::makeState(columns, rows, 2, 42);
    climateweather::ShallowWaterForcing zeroForcing;
    auto diagnostics = climateweather::advance(resting, config, zeroForcing, 3600.0f);
    expect(diagnostics.finite && diagnostics.bounded &&
            std::abs(diagnostics.areaWeightedMassAnomaly) < 1.0e-9 &&
            diagnostics.areaWeightedEnergy == 0.0,
        "zero forcing must preserve the exact resting shallow-water state");

    auto perturbed = climateweather::makeState(columns, rows, 2, 91);
    const std::size_t centre = static_cast<std::size_t>(rows / 2) * columns + columns / 2;
    perturbed.layers[0].heightAnomalyMetres[centre] = 12.0f;
    perturbed.layers[1].heightAnomalyMetres[centre] = -6.0f;
    diagnostics = climateweather::advance(perturbed, config, zeroForcing, 21600.0f);
    const double lowerSpeed = [&]
    {
        double total = 0.0;
        for (std::size_t cell = 0; cell < cellCount; cell++)
            total += std::hypot(
                perturbed.layers[0].eastWindMps[cell],
                perturbed.layers[0].southWindMps[cell]);
        return total;
    }();
    expect(diagnostics.substeps > 1 && diagnostics.finite && diagnostics.bounded &&
            std::abs(diagnostics.areaWeightedMassAnomaly) < 1.0e-5,
        "the prognostic atmosphere must subcycle stably and conserve layer mass");
    expect(lowerSpeed > 0.0 && diagnostics.areaWeightedEnergy > 0.0,
        "an idealized height anomaly must launch a resolved dynamical response");

    const auto bytes = climateweather::serializeState(perturbed);
    climateweather::ShallowWaterState restored;
    expect(climateweather::deserializeState(bytes, restored) &&
            climateweather::serializeState(restored) == bytes,
        "weather state serialization must round-trip bit-exactly");
    auto replayFirst = restored;
    auto replaySecond = restored;
    config.stochasticHeightForcingMetresPerSecond = 2.0e-5f;
    climateweather::advance(replayFirst, config, zeroForcing, 10800.0f);
    climateweather::advance(replaySecond, config, zeroForcing, 10800.0f);
    expect(climateweather::serializeState(replayFirst) ==
            climateweather::serializeState(replaySecond),
        "fixed seed and serialized RNG state must replay the same weather step");

    auto weatherInitial = climateweather::makeState(columns, rows, 2, 123456);
    const auto sequence = climateweather::generateWeatherSequence(
        weatherInitial, config, zeroForcing, 24, 10800.0f);
    const auto statistics = climateweather::calculateStatistics(sequence, 0);
    expect(sequence.size() == 24 && statistics.sampleCount == 24,
        "a season must retain at least 20 statistically distinct weather states");
    bool distinct = false;
    if (sequence.size() >= 2)
    {
        distinct = sequence.front().layers[0].heightAnomalyMetres !=
            sequence.back().layers[0].heightAnomalyMetres;
    }
    expect(distinct,
        "weather samples must come from one evolving state rather than duplicated means");
    expect(std::all_of(
            statistics.directionalConsistency.begin(),
            statistics.directionalConsistency.end(),
            [](float value) { return value >= 0.0f && value <= 1.0f; }),
        "directional consistency must remain a bounded temporal statistic");

    std::vector<climateweather::ShallowWaterState> irregular(3,
        climateweather::makeState(columns, rows, 1, 7));
    irregular[0].elapsedSeconds = 0.0;
    irregular[1].elapsedSeconds = 1.0;
    irregular[2].elapsedSeconds = 4.0;
    std::fill(irregular[1].layers[0].eastWindMps.begin(),
        irregular[1].layers[0].eastWindMps.end(), 10.0f);
    const auto weighted = climateweather::calculateStatistics(irregular);
    expect(weighted.sampleCount == 3 && weighted.durationSeconds == 4.0 &&
            std::abs(weighted.meanEastWindMps[0] - 5.0f) < 1.0e-6f &&
            weighted.speedStandardErrorMps[0] > 0.0f && weighted.effectiveSampleCount < 3.0,
        "irregular weather sampling must be temporally weighted with reported sampling spread");
    irregular[2].layers[0].southWindMps.clear();
    expect(climateweather::calculateStatistics(irregular).sampleCount == 0,
        "malformed sampled states must not produce misleading statistics");
    std::vector<climateweather::ShallowWaterState> daily(90, climateweather::makeState(4, 2, 1, 9));
    for (int day = 0; day < 90; ++day)
    {
        daily[day].elapsedSeconds = (day + 0.5) * 86400.0;
        std::fill(daily[day].layers[0].eastWindMps.begin(), daily[day].layers[0].eastWindMps.end(),
            5.0f + static_cast<float>(day) / 20.0f);
    }
    const auto dailyStatistics = climateweather::calculateStatistics(daily, 0, std::vector<double>(90, 86400.0));
    expect(dailyStatistics.durationSeconds == 90.0 * 86400.0 && dailyStatistics.sampleCount == 90 &&
        dailyStatistics.decorrelatedSampleCount[0] < dailyStatistics.effectiveSampleCount &&
        dailyStatistics.correlatedSpeedStandardErrorMps[0] >= dailyStatistics.speedStandardErrorMps[0],
        "daily midpoint samples must cover the full quarter without claiming correlated states are independent");

    const auto coarse = climateweather::resampleState(replayFirst, columns / 2);
    const double originalMass = climategrid::areaWeightedIntegral(columns, rows,
        climategrid::LatitudeLayout::cellCentred, replayFirst.layers[0].heightAnomalyMetres);
    const double coarseMass = climategrid::areaWeightedIntegral(columns / 2, rows / 2,
        climategrid::LatitudeLayout::cellCentred, coarse.layers[0].heightAnomalyMetres);
    expect(coarse.columns == columns / 2 && coarse.rows == rows / 2 &&
            coarse.randomState == replayFirst.randomState &&
            coarse.elapsedSeconds == replayFirst.elapsedSeconds &&
            std::abs(coarseMass - originalMass) < 1.0e-5,
        "weather LOD changes must conservatively remap state without resetting time or randomness");

    auto advected = restored;
    auto unadvected = restored;
    climateweather::ShallowWaterForcing jets;
    jets.backgroundEastWindMps.assign(2, std::vector<float>(cellCount, 15.0f));
    climateweather::advance(advected, config, jets, 10800.0f);
    climateweather::advance(unadvected, config, zeroForcing, 10800.0f);
    expect(advected.layers[0].heightAnomalyMetres != unadvected.layers[0].heightAnomalyMetres,
        "the evolving weather must propagate along its parent seasonal jets");

    auto truncated = bytes;
    truncated.resize(48);
    const auto savedRestored = climateweather::serializeState(restored);
    expect(!climateweather::deserializeState(truncated, restored) &&
            climateweather::serializeState(restored) == savedRestored,
        "truncated replay data must fail atomically without mutating the previous state");
    config.lowerMeanDepthMetres = -1.0f;
    diagnostics = climateweather::advance(restored, config, zeroForcing, 3600.0f);
    expect(!diagnostics.bounded && climateweather::serializeState(restored) == savedRestored,
        "invalid weather configurations must fail without mutating replay state");

    auto momentumState = climateweather::makeState(columns, rows, 2, 8);
    climateweather::ShallowWaterConfig momentumConfig;
    momentumConfig.layerCount = 2;
    momentumConfig.rotationRatePerSecond = 0.0f;
    momentumConfig.lowerMeanDepthMetres = 200.0f;
    momentumConfig.upperMeanDepthMetres = 100.0f;
    momentumConfig.lowerDragTimeSeconds = momentumConfig.upperDragTimeSeconds = 1.0e20f;
    std::fill(momentumState.layers[0].eastWindMps.begin(), momentumState.layers[0].eastWindMps.end(), 4.0f);
    std::fill(momentumState.layers[1].eastWindMps.begin(), momentumState.layers[1].eastWindMps.end(), -2.0f);
    climateweather::advance(momentumState, momentumConfig, {}, 3600.0f);
    expect(std::abs(200.0f * momentumState.layers[0].eastWindMps[0] + 100.0f * momentumState.layers[1].eastWindMps[0] - 600.0f) < 0.001f,
        "interlayer exchange must conserve depth-weighted momentum for unequal equivalent depths");

    climateweather::ShallowWaterConfig seasonalConfig;
    seasonalConfig.layerCount = 2;
    seasonalConfig.lowerMeanDepthMetres = 120.0f;
    seasonalConfig.upperMeanDepthMetres = 60.0f;
    seasonalConfig.lowerDragTimeSeconds = 43200.0f;
    seasonalConfig.upperDragTimeSeconds = 86400.0f;
    seasonalConfig.heightRelaxationTimeSeconds = 68000.0f;
    auto seasonal = climateweather::makeState(64, 32, 2, 73021);
    climateweather::ShallowWaterForcing seasonalForcing;
    seasonalForcing.equilibriumHeightMetres.assign(2, std::vector<float>(64 * 32));
    seasonalForcing.backgroundEastWindMps.assign(2, std::vector<float>(64 * 32, 20.0f));
    for (int y = 0; y < 32; ++y)
        for (int x = 0; x < 64; ++x)
        {
            const float shape = std::cos(6.283185307f * x / 64.0f) * std::sin(3.141592654f * (y + 0.5f) / 32.0f);
            seasonalForcing.equilibriumHeightMetres[0][y * 64 + x] = 10.0f * shape;
            seasonalForcing.equilibriumHeightMetres[1][y * 64 + x] = 230.0f * shape;
        }
    for (int step = 0; step < 12; ++step)
    {
        const auto result = climateweather::advance(seasonal, seasonalConfig, seasonalForcing, 21600.0f);
        expect(result.finite && result.bounded, "64x32 forced climate spin-up must remain inside modal amplitude bounds");
    }
    const auto climateSamples = climateweather::generateWeatherSequence(seasonal, seasonalConfig, seasonalForcing, 30, 21600.0f);
    expect(climateSamples.size() == 30 && climateweather::calculateStatistics(climateSamples, 1).sampleCount == 30,
        "prognostic climatology must retain thirty evolving states after spin-up");
    return failures == 0 ? 0 : 1;
}
