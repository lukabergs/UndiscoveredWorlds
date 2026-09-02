#include "climate_weather.hpp"

#include "climate_grid.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <type_traits>

namespace climateweather
{
namespace
{
constexpr double pi = 3.1415926535897932384626433832795;

std::uint64_t nextRandom(std::uint64_t& state)
{
    if (state == 0)
        state = 0x9e3779b97f4a7c15ULL;
    state ^= state >> 12;
    state ^= state << 25;
    state ^= state >> 27;
    return state * 0x2545f4914f6cdd1dULL;
}

float randomUnit(std::uint64_t& state)
{
    return static_cast<float>((nextRandom(state) >> 40) * (1.0 / 16777216.0));
}

bool validForcing(
    const ShallowWaterForcing& forcing,
    int layerCount,
    std::size_t cellCount)
{
    const auto validFields = [=](const std::vector<std::vector<float>>& fields)
    {
        if (fields.empty())
            return true;
        if (fields.size() != static_cast<std::size_t>(layerCount))
            return false;
        return std::all_of(fields.begin(), fields.end(), [=](const auto& field)
            { return field.size() == cellCount && std::all_of(field.begin(), field.end(),
                [](float value) { return std::isfinite(value); }); });
    };
    return validFields(forcing.equilibriumHeightMetres) &&
        validFields(forcing.heightTendencyMetresPerSecond) &&
        validFields(forcing.backgroundEastWindMps) && validFields(forcing.backgroundSouthWindMps);
}

float configuredDepth(const ShallowWaterConfig& config, int layer)
{
    return layer == 0 ? config.lowerMeanDepthMetres : config.upperMeanDepthMetres;
}

float configuredDragTime(const ShallowWaterConfig& config, int layer)
{
    return layer == 0 ? config.lowerDragTimeSeconds : config.upperDragTimeSeconds;
}

template<typename T>
void appendScalar(std::vector<std::uint8_t>& bytes, T value)
{
    static_assert(std::is_trivially_copyable<T>::value, "serialized values must be POD");
    const auto* first = reinterpret_cast<const std::uint8_t*>(&value);
    bytes.insert(bytes.end(), first, first + sizeof(T));
}

template<typename T>
bool readScalar(const std::vector<std::uint8_t>& bytes, std::size_t& offset, T& value)
{
    if (offset + sizeof(T) > bytes.size())
        return false;
    std::memcpy(&value, bytes.data() + offset, sizeof(T));
    offset += sizeof(T);
    return true;
}
}

ShallowWaterState makeState(
    int columns,
    int rows,
    int layerCount,
    std::uint64_t seed)
{
    ShallowWaterState state;
    if (columns <= 0 || rows <= 0 || layerCount < 1 || layerCount > 2)
        return state;
    state.columns = columns;
    state.rows = rows;
    state.randomState = seed == 0 ? 1 : seed;
    const std::size_t cellCount = static_cast<std::size_t>(columns) * rows;
    state.layers.resize(layerCount);
    for (ShallowWaterLayer& layer : state.layers)
    {
        layer.heightAnomalyMetres.assign(cellCount, 0.0f);
        layer.eastWindMps.assign(cellCount, 0.0f);
        layer.southWindMps.assign(cellCount, 0.0f);
    }
    return state;
}

ShallowWaterState resampleState(const ShallowWaterState& state, int columns)
{
    auto result = makeState(columns, columns / 2, static_cast<int>(state.layers.size()), state.randomState);
    if (result.layers.empty())
        return result;
    result.elapsedSeconds = state.elapsedSeconds;
    const auto remap = [&](const auto& field)
    {
        return climategrid::remapField(state.columns, state.rows,
            climategrid::LatitudeLayout::cellCentred, field, result.columns, result.rows,
            climategrid::LatitudeLayout::cellCentred);
    };
    for (std::size_t layer = 0; layer < state.layers.size(); layer++)
    {
        result.layers[layer].heightAnomalyMetres = remap(state.layers[layer].heightAnomalyMetres);
        result.layers[layer].eastWindMps = remap(state.layers[layer].eastWindMps);
        result.layers[layer].southWindMps = remap(state.layers[layer].southWindMps);
    }
    return result;
}

ShallowWaterDiagnostics advance(
    ShallowWaterState& state,
    const ShallowWaterConfig& config,
    const ShallowWaterForcing& forcing,
    float elapsedSeconds)
{
    ShallowWaterDiagnostics diagnostics;
    const int columns = state.columns;
    const int rows = state.rows;
    const int layerCount = static_cast<int>(state.layers.size());
    const std::size_t cellCount = static_cast<std::size_t>(std::max(0, columns)) *
        static_cast<std::size_t>(std::max(0, rows));
    if (columns < 3 || rows < 2 || !std::isfinite(elapsedSeconds) || elapsedSeconds <= 0.0f ||
        layerCount != std::clamp(config.layerCount, 1, 2) ||
        !validForcing(forcing, layerCount, cellCount) ||
        !(config.planetRadiusMetres > 0.0f) || !(config.gravityMetresPerSecondSquared > 0.0f) ||
        !(config.lowerMeanDepthMetres > 0.0f) || !(config.upperMeanDepthMetres > 0.0f))
    {
        diagnostics.finite = false;
        diagnostics.bounded = false;
        return diagnostics;
    }
    for (const ShallowWaterLayer& layer : state.layers)
    {
        if (layer.heightAnomalyMetres.size() != cellCount ||
            layer.eastWindMps.size() != cellCount || layer.southWindMps.size() != cellCount ||
            !std::all_of(layer.heightAnomalyMetres.begin(), layer.heightAnomalyMetres.end(),
                [](float value) { return std::isfinite(value); }) ||
            !std::all_of(layer.eastWindMps.begin(), layer.eastWindMps.end(),
                [](float value) { return std::isfinite(value); }) ||
            !std::all_of(layer.southWindMps.begin(), layer.southWindMps.end(),
                [](float value) { return std::isfinite(value); }))
        {
            diagnostics.finite = false;
            diagnostics.bounded = false;
            return diagnostics;
        }
    }

    const climategrid::SphericalGrid grid = climategrid::makeSphericalGrid(
        columns, rows, config.planetRadiusMetres);
    const auto originalState = state;
    const auto background = [](const auto& fields, int layer, std::size_t cell)
    {
        return fields.empty() ? 0.0f : fields[layer][cell];
    };
    double minimumSpacing = std::numeric_limits<double>::infinity();
    double maximumSignalSpeed = 0.0;
    for (int y = 0; y < rows; y++)
    {
        const double zonalSpacing = grid.cellAreasSquareMetres[y] /
            std::max(1.0, grid.zonalFaceLengthsMetres[y]);
        minimumSpacing = std::min(
            minimumSpacing,
            std::min(zonalSpacing, config.planetRadiusMetres * grid.latitudeSpacingRadians));
    }
    for (int layer = 0; layer < layerCount; layer++)
    {
        const double waveSpeed = std::sqrt(
            config.gravityMetresPerSecondSquared * configuredDepth(config, layer));
        for (std::size_t cell = 0; cell < cellCount; cell++)
        {
            maximumSignalSpeed = std::max(
                maximumSignalSpeed,
                waveSpeed * std::sqrt(1.0 + 2.0 * std::abs(config.baroclinicCoupling)) +
                    std::hypot(state.layers[layer].eastWindMps[cell] +
                        background(forcing.backgroundEastWindMps, layer, cell),
                        state.layers[layer].southWindMps[cell] +
                        background(forcing.backgroundSouthWindMps, layer, cell)));
        }
    }
    diagnostics.maximumCourant = static_cast<float>(
        2.0 * maximumSignalSpeed * elapsedSeconds / minimumSpacing);
    diagnostics.substeps = std::max(
        1,
        static_cast<int>(std::ceil(
            diagnostics.maximumCourant / std::clamp(config.maximumCourant, 0.05f, 0.9f))));
    const double sourceRate = 2.0 * std::abs(config.rotationRatePerSecond) +
        1.0 / std::max(1.0f, config.lowerDragTimeSeconds) +
        1.0 / std::max(1.0f, config.upperDragTimeSeconds) +
        (config.interlayerMomentumTimeSeconds > 0.0f ? 2.0 / config.interlayerMomentumTimeSeconds : 0.0) +
        (config.heightRelaxationTimeSeconds > 0.0f ? 1.0 / config.heightRelaxationTimeSeconds : 0.0);
    diagnostics.substeps = std::max(diagnostics.substeps,
        static_cast<int>(std::ceil(elapsedSeconds * sourceRate / 0.25)));
    const float timeStep = elapsedSeconds / static_cast<float>(diagnostics.substeps);

    std::vector<std::vector<float>> stochastic(
        layerCount, std::vector<float>(cellCount, 0.0f));
    if (config.stochasticHeightForcingMetresPerSecond != 0.0f)
    {
        for (int layer = 0; layer < layerCount; layer++)
        {
            const int wavenumber = 2 + static_cast<int>(nextRandom(state.randomState) % 4);
            const float phase = 2.0f * static_cast<float>(pi) * randomUnit(state.randomState);
            const float sign = randomUnit(state.randomState) < 0.5f ? -1.0f : 1.0f;
            for (int y = 0; y < rows; y++)
            {
                const float latitude = static_cast<float>(grid.latitudeCentresRadians[y]);
                const float envelope = std::cos(latitude) *
                    std::sin(2.0f * latitude) * sign;
                for (int x = 0; x < columns; x++)
                {
                    const float longitude = 2.0f * static_cast<float>(pi) *
                        (static_cast<float>(x) + 0.5f) / static_cast<float>(columns);
                    stochastic[layer][grid.index(x, y)] =
                        config.stochasticHeightForcingMetresPerSecond * envelope *
                        std::sin(wavenumber * longitude + phase);
                }
            }
        }
    }

    struct Tendency
    {
        std::vector<float> height;
        std::vector<float> east;
        std::vector<float> south;
    };
    std::vector<Tendency> first(layerCount);
    std::vector<Tendency> midpoint(layerCount);
    std::vector<ShallowWaterLayer> middle = state.layers;
    for (int layer = 0; layer < layerCount; layer++)
    {
        for (Tendency* tendency : { &first[layer], &midpoint[layer] })
        {
            tendency->height.assign(cellCount, 0.0f);
            tendency->east.assign(cellCount, 0.0f);
            tendency->south.assign(cellCount, 0.0f);
        }
    }

    const auto calculateTendencies = [&](const std::vector<ShallowWaterLayer>& layers,
                                         std::vector<Tendency>& tendencies)
    {
        for (int layer = 0; layer < layerCount; layer++)
        {
            const float depth = configuredDepth(config, layer);
            const float dragTime = std::max(1.0f, configuredDragTime(config, layer));
            for (int y = 0; y < rows; y++)
            {
                const int north = std::max(0, y - 1);
                const int south = std::min(rows - 1, y + 1);
                const double northSouthSpacing = config.planetRadiusMetres *
                    grid.latitudeSpacingRadians * static_cast<double>(south - north);
                const float latitude = static_cast<float>(grid.latitudeCentresRadians[y]);
                const float coriolis = 2.0f * config.rotationRatePerSecond *
                    config.rotationDirection * std::sin(latitude);
                const double zonalSpacing = grid.cellAreasSquareMetres[y] /
                    std::max(1.0, grid.zonalFaceLengthsMetres[y]);
                for (int x = 0; x < columns; x++)
                {
                    const std::size_t cell = grid.index(x, y);
                    const std::size_t west = grid.index(x - 1, y);
                    const std::size_t east = grid.index(x + 1, y);
                    const std::size_t northCell = grid.index(x, north);
                    const std::size_t southCell = grid.index(x, south);
                    const auto effectiveHeight = [&](std::size_t at)
                    {
                        float height = layers[layer].heightAnomalyMetres[at];
                        if (layerCount == 2)
                        {
                            const int other = 1 - layer;
                            height += config.baroclinicCoupling *
                                (layers[layer].heightAnomalyMetres[at] -
                                    layers[other].heightAnomalyMetres[at]);
                        }
                        return height;
                    };
                    const float gradientEast = static_cast<float>(
                        (effectiveHeight(east) - effectiveHeight(west)) /
                        (2.0 * zonalSpacing));
                    const float gradientSouth = north == south ? 0.0f : static_cast<float>(
                        (effectiveHeight(southCell) - effectiveHeight(northCell)) /
                        northSouthSpacing);
                    float eastTendency = -config.gravityMetresPerSecondSquared * gradientEast -
                        coriolis * layers[layer].southWindMps[cell] -
                        layers[layer].eastWindMps[cell] / dragTime;
                    float southTendency = -config.gravityMetresPerSecondSquared * gradientSouth +
                        coriolis * layers[layer].eastWindMps[cell] -
                        layers[layer].southWindMps[cell] / dragTime;
                    const float backgroundEast = background(forcing.backgroundEastWindMps, layer, cell);
                    const float backgroundSouth = background(forcing.backgroundSouthWindMps, layer, cell);
                    const auto materialAdvection = [&](const std::vector<float>& field)
                    {
                        const double dx = backgroundEast >= 0.0f
                            ? (field[cell] - field[west]) / zonalSpacing
                            : (field[east] - field[cell]) / zonalSpacing;
                        const double dy = backgroundSouth >= 0.0f
                            ? (field[cell] - field[northCell]) /
                                (config.planetRadiusMetres * grid.latitudeSpacingRadians)
                            : (field[southCell] - field[cell]) /
                                (config.planetRadiusMetres * grid.latitudeSpacingRadians);
                        return static_cast<float>(backgroundEast * dx + backgroundSouth * dy);
                    };
                    eastTendency -= materialAdvection(layers[layer].eastWindMps);
                    southTendency -= materialAdvection(layers[layer].southWindMps);
                    // Local Lax-Friedrichs wave viscosity closes the otherwise
                    // unstable centred-gravity-wave/RK2 combination. Pairwise
                    // face fluxes conserve the layer mass on the sphere.
                    const double waveSpeed = std::sqrt(config.gravityMetresPerSecondSquared * depth *
                        (1.0 + 2.0 * std::abs(config.baroclinicCoupling)));
                    const auto waveDiffusion = [&](const std::vector<float>& field)
                    {
                        double flux = grid.zonalFaceLengthsMetres[y] * (field[east] + field[west] - 2.0 * field[cell]);
                        if (y > 0) flux += grid.northFaceLengthsMetres[y] * (field[northCell] - field[cell]);
                        if (y + 1 < rows) flux += grid.southFaceLengthsMetres[y] * (field[southCell] - field[cell]);
                        return static_cast<float>(0.5 * waveSpeed * flux / grid.cellAreasSquareMetres[y]);
                    };
                    eastTendency += waveDiffusion(layers[layer].eastWindMps);
                    southTendency += waveDiffusion(layers[layer].southWindMps);
                    if (layerCount == 2 && config.interlayerMomentumTimeSeconds > 0.0f)
                    {
                        const int other = 1 - layer;
                        const float exchangeRate = 2.0f * configuredDepth(config, other) /
                            (depth + configuredDepth(config, other)) / config.interlayerMomentumTimeSeconds;
                        eastTendency += (layers[other].eastWindMps[cell] -
                            layers[layer].eastWindMps[cell]) * exchangeRate;
                        southTendency += (layers[other].southWindMps[cell] -
                            layers[layer].southWindMps[cell]) * exchangeRate;
                    }

                    const auto backgroundHeightFlux = [&](std::size_t left, std::size_t right,
                                                           const auto& winds, double faceLength)
                    {
                        const float velocity = 0.5f * (background(winds, layer, left) +
                            background(winds, layer, right));
                        return velocity * faceLength * layers[layer].heightAnomalyMetres[
                            velocity >= 0.0f ? left : right];
                    };
                    const double eastFlux = backgroundHeightFlux(cell, east,
                        forcing.backgroundEastWindMps, grid.zonalFaceLengthsMetres[y]) + 0.5 * depth *
                        (layers[layer].eastWindMps[cell] + layers[layer].eastWindMps[east]) *
                        grid.zonalFaceLengthsMetres[y];
                    const double westFlux = backgroundHeightFlux(west, cell,
                        forcing.backgroundEastWindMps, grid.zonalFaceLengthsMetres[y]) + 0.5 * depth *
                        (layers[layer].eastWindMps[west] + layers[layer].eastWindMps[cell]) *
                        grid.zonalFaceLengthsMetres[y];
                    const double southFlux = y == rows - 1 ? 0.0 :
                        backgroundHeightFlux(cell, southCell, forcing.backgroundSouthWindMps,
                            grid.southFaceLengthsMetres[y]) + 0.5 * depth *
                        (layers[layer].southWindMps[cell] +
                            layers[layer].southWindMps[southCell]) *
                        grid.southFaceLengthsMetres[y];
                    const double northFlux = y == 0 ? 0.0 :
                        backgroundHeightFlux(northCell, cell, forcing.backgroundSouthWindMps,
                            grid.northFaceLengthsMetres[y]) + 0.5 * depth *
                        (layers[layer].southWindMps[northCell] +
                            layers[layer].southWindMps[cell]) *
                        grid.northFaceLengthsMetres[y];
                    float heightTendency = static_cast<float>(
                        -(eastFlux - westFlux + southFlux - northFlux) /
                        grid.cellAreasSquareMetres[y]);
                    heightTendency += waveDiffusion(layers[layer].heightAnomalyMetres);
                    if (config.heightRelaxationTimeSeconds > 0.0f)
                    {
                        const float target = forcing.equilibriumHeightMetres.empty()
                            ? 0.0f : forcing.equilibriumHeightMetres[layer][cell];
                        heightTendency += (target -
                            layers[layer].heightAnomalyMetres[cell]) /
                            config.heightRelaxationTimeSeconds;
                    }
                    if (!forcing.heightTendencyMetresPerSecond.empty())
                        heightTendency += forcing.heightTendencyMetresPerSecond[layer][cell];
                    heightTendency += stochastic[layer][cell];
                    tendencies[layer].east[cell] = eastTendency;
                    tendencies[layer].south[cell] = southTendency;
                    tendencies[layer].height[cell] = heightTendency;
                }
            }
        }
    };

    for (int substep = 0; substep < diagnostics.substeps; substep++)
    {
        calculateTendencies(state.layers, first);
        for (int layer = 0; layer < layerCount; layer++)
        {
            for (std::size_t cell = 0; cell < cellCount; cell++)
            {
                middle[layer].heightAnomalyMetres[cell] =
                    state.layers[layer].heightAnomalyMetres[cell] +
                    0.5f * timeStep * first[layer].height[cell];
                middle[layer].eastWindMps[cell] = state.layers[layer].eastWindMps[cell] +
                    0.5f * timeStep * first[layer].east[cell];
                middle[layer].southWindMps[cell] = state.layers[layer].southWindMps[cell] +
                    0.5f * timeStep * first[layer].south[cell];
            }
        }
        calculateTendencies(middle, midpoint);
        for (int layer = 0; layer < layerCount; layer++)
        {
            double heightTotal = 0.0;
            double areaTotal = 0.0;
            for (int y = 0; y < rows; y++)
            {
                const double area = grid.cellAreasSquareMetres[y];
                for (int x = 0; x < columns; x++)
                {
                    const std::size_t cell = grid.index(x, y);
                    state.layers[layer].heightAnomalyMetres[cell] +=
                        timeStep * midpoint[layer].height[cell];
                    state.layers[layer].eastWindMps[cell] +=
                        timeStep * midpoint[layer].east[cell];
                    state.layers[layer].southWindMps[cell] +=
                        timeStep * midpoint[layer].south[cell];
                    heightTotal += area * state.layers[layer].heightAnomalyMetres[cell];
                    areaTotal += area;
                }
            }
            const float meanHeight = areaTotal > 0.0
                ? static_cast<float>(heightTotal / areaTotal)
                : 0.0f;
            for (std::size_t cell = 0; cell < cellCount; cell++)
            {
                state.layers[layer].heightAnomalyMetres[cell] -= meanHeight;
                const float height = state.layers[layer].heightAnomalyMetres[cell];
                const float east = state.layers[layer].eastWindMps[cell];
                const float south = state.layers[layer].southWindMps[cell];
                diagnostics.finite = diagnostics.finite && std::isfinite(height) &&
                    std::isfinite(east) && std::isfinite(south);
                diagnostics.bounded = diagnostics.bounded &&
                    std::abs(height) < config.maximumHeightAnomalyMetres &&
                    std::hypot(east, south) < config.maximumAnomalyWindMps;
            }
        }
        if (!diagnostics.finite || !diagnostics.bounded)
        {
            state = originalState;
            return diagnostics;
        }
    }
    state.elapsedSeconds += elapsedSeconds;

    double areaTotal = 0.0;
    for (int layer = 0; layer < layerCount; layer++)
    {
        const float depth = configuredDepth(config, layer);
        for (int y = 0; y < rows; y++)
        {
            const double area = grid.cellAreasSquareMetres[y];
            for (int x = 0; x < columns; x++)
            {
                const std::size_t cell = grid.index(x, y);
                const double height = state.layers[layer].heightAnomalyMetres[cell];
                const double speedSquared =
                    state.layers[layer].eastWindMps[cell] *
                        state.layers[layer].eastWindMps[cell] +
                    state.layers[layer].southWindMps[cell] *
                        state.layers[layer].southWindMps[cell];
                diagnostics.areaWeightedMassAnomaly += area * height;
                diagnostics.areaWeightedEnergy += area *
                    (0.5 * depth * speedSquared +
                        0.5 * config.gravityMetresPerSecondSquared * height * height);
                areaTotal += area;
            }
        }
    }
    if (areaTotal > 0.0)
    {
        diagnostics.areaWeightedMassAnomaly /= areaTotal;
        diagnostics.areaWeightedEnergy /= areaTotal;
    }
    return diagnostics;
}

std::vector<ShallowWaterState> generateWeatherSequence(
    ShallowWaterState initialState,
    const ShallowWaterConfig& config,
    const ShallowWaterForcing& forcing,
    int sampleCount,
    float sampleIntervalSeconds)
{
    std::vector<ShallowWaterState> sequence;
    if (sampleCount <= 0 || sampleIntervalSeconds <= 0.0f)
        return sequence;
    sequence.reserve(sampleCount);
    for (int sample = 0; sample < sampleCount; sample++)
    {
        const ShallowWaterDiagnostics diagnostics = advance(
            initialState, config, forcing, sampleIntervalSeconds);
        if (!diagnostics.finite || !diagnostics.bounded)
            break;
        sequence.push_back(initialState);
    }
    return sequence;
}

WeatherStatistics calculateStatistics(
    const std::vector<ShallowWaterState>& states,
    int layer)
{
    WeatherStatistics statistics;
    if (states.empty() || layer < 0 ||
        layer >= static_cast<int>(states.front().layers.size()))
    {
        return statistics;
    }
    const std::size_t cellCount = states.front().layers[layer].eastWindMps.size();
    if (cellCount != static_cast<std::size_t>(states.front().columns) * states.front().rows)
        return {};
    std::vector<double> weights(states.size(), 1.0);
    if (states.size() > 1)
    {
        std::fill(weights.begin(), weights.end(), 0.0);
        for (std::size_t sample = 1; sample < states.size(); sample++)
        {
            const double duration = states[sample].elapsedSeconds - states[sample - 1].elapsedSeconds;
            if (!(duration > 0.0) || !std::isfinite(duration))
                return {};
            weights[sample - 1] += duration * 0.5;
            weights[sample] += duration * 0.5;
            statistics.durationSeconds += duration;
        }
    }
    double totalWeight = 0.0, squaredWeight = 0.0;
    for (double weight : weights)
    {
        totalWeight += weight;
        squaredWeight += weight * weight;
    }
    statistics.effectiveSampleCount = totalWeight * totalWeight / squaredWeight;
    std::vector<double> heightSum(cellCount), eastSum(cellCount), southSum(cellCount),
        speedSum(cellCount), squaredSpeedSum(cellCount);
    statistics.meanHeightAnomalyMetres.assign(cellCount, 0.0f);
    statistics.meanEastWindMps.assign(cellCount, 0.0f);
    statistics.meanSouthWindMps.assign(cellCount, 0.0f);
    statistics.meanSpeedMps.assign(cellCount, 0.0f);
    statistics.directionalConsistency.assign(cellCount, 0.0f);
    statistics.speedStandardDeviationMps.assign(cellCount, 0.0f);
    statistics.speedStandardErrorMps.assign(cellCount, 0.0f);
    for (std::size_t sample = 0; sample < states.size(); sample++)
    {
        const auto& state = states[sample];
        if (layer >= static_cast<int>(state.layers.size()) ||
            state.columns != states.front().columns || state.rows != states.front().rows ||
            state.layers[layer].eastWindMps.size() != cellCount ||
            state.layers[layer].southWindMps.size() != cellCount ||
            state.layers[layer].heightAnomalyMetres.size() != cellCount)
        {
            return {};
        }
        for (std::size_t cell = 0; cell < cellCount; cell++)
        {
            const float east = state.layers[layer].eastWindMps[cell];
            const float south = state.layers[layer].southWindMps[cell];
            const double height = state.layers[layer].heightAnomalyMetres[cell];
            const double speed = std::hypot(east, south);
            if (!std::isfinite(speed) || !std::isfinite(height))
                return {};
            heightSum[cell] += weights[sample] * height;
            eastSum[cell] += weights[sample] * east;
            southSum[cell] += weights[sample] * south;
            speedSum[cell] += weights[sample] * speed;
            squaredSpeedSum[cell] += weights[sample] * speed * speed;
        }
    }
    statistics.sampleCount = static_cast<int>(states.size());
    for (std::size_t cell = 0; cell < cellCount; cell++)
    {
        statistics.meanHeightAnomalyMetres[cell] = static_cast<float>(heightSum[cell] / totalWeight);
        statistics.meanEastWindMps[cell] = static_cast<float>(eastSum[cell] / totalWeight);
        statistics.meanSouthWindMps[cell] = static_cast<float>(southSum[cell] / totalWeight);
        statistics.meanSpeedMps[cell] = static_cast<float>(speedSum[cell] / totalWeight);
        const double variance = std::max(0.0, squaredSpeedSum[cell] / totalWeight -
            std::pow(speedSum[cell] / totalWeight, 2));
        statistics.speedStandardDeviationMps[cell] = static_cast<float>(std::sqrt(variance));
        statistics.speedStandardErrorMps[cell] = static_cast<float>(
            std::sqrt(variance / statistics.effectiveSampleCount));
        statistics.directionalConsistency[cell] = statistics.meanSpeedMps[cell] > 0.0f
            ? std::clamp(
                std::hypot(
                    statistics.meanEastWindMps[cell], statistics.meanSouthWindMps[cell]) /
                    statistics.meanSpeedMps[cell],
                0.0f,
                1.0f)
            : 0.0f;
    }
    return statistics;
}

std::vector<std::uint8_t> serializeState(const ShallowWaterState& state)
{
    std::vector<std::uint8_t> bytes;
    bytes.reserve(32);
    appendScalar<std::uint32_t>(bytes, 0x57535755U);
    appendScalar<std::uint32_t>(bytes, 1U);
    appendScalar<std::int32_t>(bytes, state.columns);
    appendScalar<std::int32_t>(bytes, state.rows);
    appendScalar<std::uint32_t>(bytes, static_cast<std::uint32_t>(state.layers.size()));
    appendScalar<double>(bytes, state.elapsedSeconds);
    appendScalar<std::uint64_t>(bytes, state.randomState);
    for (const ShallowWaterLayer& layer : state.layers)
    {
        for (const std::vector<float>* field : {
                 &layer.heightAnomalyMetres, &layer.eastWindMps, &layer.southWindMps })
        {
            appendScalar<std::uint64_t>(bytes, static_cast<std::uint64_t>(field->size()));
            for (float value : *field)
                appendScalar<float>(bytes, value);
        }
    }
    return bytes;
}

bool deserializeState(
    const std::vector<std::uint8_t>& bytes,
    ShallowWaterState& state)
{
    std::size_t offset = 0;
    std::uint32_t magic = 0;
    std::uint32_t version = 0;
    std::int32_t columns = 0;
    std::int32_t rows = 0;
    std::uint32_t layerCount = 0;
    ShallowWaterState decoded;
    if (!readScalar(bytes, offset, magic) || !readScalar(bytes, offset, version) ||
        !readScalar(bytes, offset, columns) || !readScalar(bytes, offset, rows) ||
        !readScalar(bytes, offset, layerCount) || magic != 0x57535755U || version != 1U ||
        columns < 3 || rows < 2 || columns > 4096 || rows > 2048 || layerCount < 1 || layerCount > 2 ||
        !readScalar(bytes, offset, decoded.elapsedSeconds) ||
        !readScalar(bytes, offset, decoded.randomState) ||
        !std::isfinite(decoded.elapsedSeconds) || decoded.elapsedSeconds < 0.0)
    {
        return false;
    }
    decoded.columns = columns;
    decoded.rows = rows;
    decoded.layers.resize(layerCount);
    const std::uint64_t expectedCount = static_cast<std::uint64_t>(columns) * rows;
    for (ShallowWaterLayer& layer : decoded.layers)
    {
        for (std::vector<float>* field : {
                 &layer.heightAnomalyMetres, &layer.eastWindMps, &layer.southWindMps })
        {
            std::uint64_t count = 0;
            if (!readScalar(bytes, offset, count) || count != expectedCount ||
                count > (bytes.size() - offset) / sizeof(float))
                return false;
            field->resize(static_cast<std::size_t>(count));
            for (float& value : *field)
            {
                if (!readScalar(bytes, offset, value) || !std::isfinite(value))
                    return false;
            }
        }
    }
    if (offset != bytes.size())
        return false;
    state = std::move(decoded);
    return true;
}
}
