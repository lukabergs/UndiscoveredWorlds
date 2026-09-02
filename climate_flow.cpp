#include "climate_flow.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

namespace climateflow
{
namespace
{
constexpr float pi = 3.14159265358979323846f;

float wrap(float x, int columns)
{
    return std::fmod(std::fmod(x, static_cast<float>(columns)) +
        static_cast<float>(columns), static_cast<float>(columns));
}

std::pair<float, float> sample(
    const std::vector<float>& east,
    const std::vector<float>& south,
    int columns,
    int rows,
    float x,
    float y)
{
    const float wrappedX = wrap(x, columns);
    const float boundedY = std::clamp(y, 0.0f, static_cast<float>(rows - 1));
    const int x0 = static_cast<int>(std::floor(wrappedX));
    const int x1 = (x0 + 1) % columns;
    const int y0 = static_cast<int>(std::floor(boundedY));
    const int y1 = std::min(rows - 1, y0 + 1);
    const float fractionX = wrappedX - std::floor(wrappedX);
    const float fractionY = boundedY - std::floor(boundedY);
    const auto interpolate = [&](const std::vector<float>& field)
    {
        const auto at = [&](int column, int row)
        {
            return field[static_cast<std::size_t>(row) * columns + column];
        };
        const float north = at(x0, y0) * (1.0f - fractionX) +
            at(x1, y0) * fractionX;
        const float southValue = at(x0, y1) * (1.0f - fractionX) +
            at(x1, y1) * fractionX;
        return north * (1.0f - fractionY) + southValue * fractionY;
    };
    return { interpolate(east), interpolate(south) };
}

std::pair<float, float> cellsPerSecond(
    std::pair<float, float> wind,
    int columns,
    int rows,
    float y,
    const ParticleTraceConfig& config)
{
    const float latitudeDegrees = config.poleInclusiveRows
        ? 90.0f - 180.0f * y / static_cast<float>(rows - 1)
        : 90.0f - 180.0f * (y + 0.5f) / static_cast<float>(rows);
    const float cosine = std::max(
        0.02f, std::abs(std::cos(latitudeDegrees * pi / 180.0f)));
    const float zonalSpacing = 2.0f * pi * config.planetRadiusMetres * cosine /
        static_cast<float>(columns);
    const float meridionalSpacing = pi * config.planetRadiusMetres /
        static_cast<float>(config.poleInclusiveRows ? rows - 1 : rows);
    return { wind.first / zonalSpacing, wind.second / meridionalSpacing };
}
}

ParticleStep advectParticleRk2(
    int columns,
    int rows,
    const std::vector<float>& eastWindMps,
    const std::vector<float>& southWindMps,
    float startX,
    float startY,
    float elapsedSeconds,
    const ParticleTraceConfig& config)
{
    ParticleStep result;
    const std::size_t cellCount = static_cast<std::size_t>(std::max(0, columns)) *
        static_cast<std::size_t>(std::max(0, rows));
    if (columns < 2 || rows < 2 || eastWindMps.size() != cellCount ||
        southWindMps.size() != cellCount || elapsedSeconds <= 0.0f ||
        config.planetRadiusMetres <= 0.0f || config.maximumSubstepCells <= 0.0f ||
        config.maximumSubsteps <= 0)
    {
        result.remainedInDomain = false;
        return result;
    }

    ParticlePoint current{ wrap(startX, columns), startY, 0.0f };
    result.points.push_back(current);
    if (startY < 0.0f || startY > static_cast<float>(rows - 1))
    {
        result.remainedInDomain = false;
        return result;
    }

    float remaining = elapsedSeconds;
    for (int substep = 0; substep < config.maximumSubsteps && remaining > 0.0f; substep++)
    {
        const auto initialWind = sample(
            eastWindMps, southWindMps, columns, rows, current.x, current.y);
        const auto initialRate = cellsPerSecond(
            initialWind, columns, rows, current.y, config);
        const float speedCellsPerSecond = std::max(
            std::abs(initialRate.first), std::abs(initialRate.second));
        if (speedCellsPerSecond < 1.0e-10f)
        {
            current.elapsedSeconds += remaining;
            result.points.push_back(current);
            remaining = 0.0f;
            break;
        }
        float timeStep = std::min(
            remaining,
            config.maximumSubstepCells / speedCellsPerSecond);

        for (int curvatureAttempt = 0; curvatureAttempt < 8; curvatureAttempt++)
        {
            const float midpointX = wrap(
                current.x + 0.5f * timeStep * initialRate.first, columns);
            const float midpointY = current.y + 0.5f * timeStep * initialRate.second;
            if (midpointY < 0.0f || midpointY > static_cast<float>(rows - 1))
                break;
            const auto midpointWind = sample(
                eastWindMps, southWindMps, columns, rows, midpointX, midpointY);
            const float initialMagnitude = std::hypot(initialWind.first, initialWind.second);
            const float midpointMagnitude = std::hypot(midpointWind.first, midpointWind.second);
            if (initialMagnitude <= 1.0e-6f || midpointMagnitude <= 1.0e-6f)
                break;
            const float cosine = std::clamp(
                (initialWind.first * midpointWind.first +
                    initialWind.second * midpointWind.second) /
                    (initialMagnitude * midpointMagnitude),
                -1.0f,
                1.0f);
            if (std::acos(cosine) <= config.maximumDirectionChangeRadians)
                break;
            timeStep *= 0.5f;
        }

        const float midpointX = wrap(
            current.x + 0.5f * timeStep * initialRate.first, columns);
        const float midpointY = current.y + 0.5f * timeStep * initialRate.second;
        if (midpointY < 0.0f || midpointY > static_cast<float>(rows - 1))
        {
            result.remainedInDomain = false;
            break;
        }
        const auto midpointWind = sample(
            eastWindMps, southWindMps, columns, rows, midpointX, midpointY);
        const auto midpointRate = cellsPerSecond(
            midpointWind, columns, rows, midpointY, config);
        const float nextX = wrap(current.x + timeStep * midpointRate.first, columns);
        const float nextY = current.y + timeStep * midpointRate.second;
        if (nextY < 0.0f || nextY > static_cast<float>(rows - 1))
        {
            result.remainedInDomain = false;
            break;
        }
        current = { nextX, nextY, current.elapsedSeconds + timeStep };
        result.points.push_back(current);
        remaining -= timeStep;
    }
    if (remaining > std::max(1.0e-3f, elapsedSeconds * 1.0e-6f))
        result.remainedInDomain = false;
    return result;
}
}
