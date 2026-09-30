#pragma once

#include <cstdint>
#include <vector>

namespace climateweather
{
struct ShallowWaterLayer
{
    std::vector<float> heightAnomalyMetres;
    std::vector<float> eastWindMps;
    std::vector<float> southWindMps;
};

struct ShallowWaterState
{
    int columns = 0;
    int rows = 0;
    double elapsedSeconds = 0.0;
    std::uint64_t randomState = 1;
    std::vector<ShallowWaterLayer> layers;
};

struct ShallowWaterForcing
{
    std::vector<std::vector<float>> equilibriumHeightMetres;
    std::vector<std::vector<float>> heightTendencyMetresPerSecond;
    // Optional seasonal jets advect anomalies; the background itself is never
    // modified. These fields use the state grid and one vector per layer.
    std::vector<std::vector<float>> backgroundEastWindMps;
    std::vector<std::vector<float>> backgroundSouthWindMps;
};

struct ShallowWaterConfig
{
    // Linear vertical-mode equations. Equivalent depths set wave speed;
    // signed geopotential anomalies are NOT the total thickness of a column.
    int layerCount = 1;
    float planetRadiusMetres = 6371000.0f;
    float rotationRatePerSecond = 7.2921159e-5f;
    float rotationDirection = 1.0f;
    float gravityMetresPerSecondSquared = 9.80665f;
    float lowerMeanDepthMetres = 250.0f;
    float upperMeanDepthMetres = 100.0f;
    float lowerDragTimeSeconds = 432000.0f;
    float upperDragTimeSeconds = 864000.0f;
    float heightRelaxationTimeSeconds = 864000.0f;
    float interlayerMomentumTimeSeconds = 432000.0f;
    float baroclinicCoupling = 0.25f;
    float maximumCourant = 0.45f;
    float maximumHeightAnomalyMetres = 2000.0f;
    float maximumAnomalyWindMps = 250.0f;
    float stochasticHeightForcingMetresPerSecond = 0.0f;
};

struct ShallowWaterDiagnostics
{
    int substeps = 1;
    float maximumCourant = 0.0f;
    double areaWeightedMassAnomaly = 0.0;
    double areaWeightedEnergy = 0.0;
    bool finite = true;
    bool bounded = true;
};

struct WeatherStatistics
{
    int sampleCount = 0;
    double durationSeconds = 0.0;
    double effectiveSampleCount = 0.0;
    std::vector<float> meanHeightAnomalyMetres;
    std::vector<float> meanEastWindMps;
    std::vector<float> meanSouthWindMps;
    std::vector<float> meanSpeedMps;
    std::vector<float> directionalConsistency;
    std::vector<float> speedStandardDeviationMps;
    // Independent-sample estimate; correlated storm-track samples need more
    // sampling, not a claim that this is a confidence bound.
    std::vector<float> speedStandardErrorMps;
    // Positive-sequence autocorrelation estimate, per cell. Descriptive only:
    // finite seasonal samples do not establish climate-estimation confidence.
    std::vector<float> decorrelatedSampleCount;
    std::vector<float> correlatedSpeedStandardErrorMps;
};

ShallowWaterState makeState(
    int columns,
    int rows,
    int layerCount,
    std::uint64_t seed = 1);
ShallowWaterState resampleState(const ShallowWaterState& state, int columns);
ShallowWaterDiagnostics advance(
    ShallowWaterState& state,
    const ShallowWaterConfig& config,
    const ShallowWaterForcing& forcing,
    float elapsedSeconds);
std::vector<ShallowWaterState> generateWeatherSequence(
    ShallowWaterState initialState,
    const ShallowWaterConfig& config,
    const ShallowWaterForcing& forcing,
    int sampleCount,
    float sampleIntervalSeconds);
WeatherStatistics calculateStatistics(
    const std::vector<ShallowWaterState>& states,
    int layer = 0,
    const std::vector<double>& sampleDurationsSeconds = {});
std::vector<std::uint8_t> serializeState(const ShallowWaterState& state);
bool deserializeState(
    const std::vector<std::uint8_t>& bytes,
    ShallowWaterState& state);
}
