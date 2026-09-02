#include "climate_flow.hpp"

#include <cmath>
#include <iostream>
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
    constexpr int columns = 32;
    constexpr int rows = 17;
    constexpr float radius = 6371000.0f;
    constexpr float pi = 3.14159265358979323846f;
    constexpr std::size_t cellCount = columns * rows;
    const auto index = [=](int x, int y)
    {
        return static_cast<std::size_t>(y) * columns + x;
    };
    std::vector<float> east(cellCount, 0.0f);
    std::vector<float> south(cellCount, 0.0f);
    const int equator = rows / 2;
    const float equatorialCellWidth = 2.0f * pi * radius / columns;
    std::fill(east.begin(), east.end(), equatorialCellWidth / 3600.0f);
    climateflow::ParticleTraceConfig config;
    config.planetRadiusMetres = radius;
    config.maximumSubstepCells = 0.2f;
    const auto straight = climateflow::advectParticleRk2(
        columns, rows, east, south, 3.0f, static_cast<float>(equator), 3600.0f, config);
    expect(straight.remainedInDomain && straight.points.size() >= 6 &&
            std::abs(straight.points.back().x - 4.0f) < 1.0e-4f &&
            std::abs(straight.points.back().elapsedSeconds - 3600.0f) < 1.0e-3f,
        "adaptive midpoint integration must preserve straight uniform flow and elapsed time");

    std::fill(east.begin(), east.end(), 5.0f * equatorialCellWidth / 3600.0f);
    const auto fast = climateflow::advectParticleRk2(
        columns, rows, east, south, 3.0f, static_cast<float>(equator), 3600.0f, config);
    expect(fast.remainedInDomain && fast.points.size() > straight.points.size() &&
            std::abs(fast.points.back().x - 8.0f) < 1.0e-3f &&
            std::abs(fast.points.back().elapsedSeconds - 3600.0f) < 1.0e-3f,
        "fast flow must use more geometry without silently reducing physical elapsed time");

    std::fill(east.begin(), east.end(), 0.0f);
    const auto stationary = climateflow::advectParticleRk2(
        columns, rows, east, south, 7.0f, static_cast<float>(equator), 3600.0f, config);
    expect(stationary.remainedInDomain && stationary.points.size() == 2 &&
            stationary.points.back().x == 7.0f &&
            stationary.points.back().elapsedSeconds == 3600.0f,
        "a valid zero-wind trace must advance time without inventing line length");

    const float angularVelocity = 2.0f * pi / 86400.0f;
    for (int y = 0; y < rows; y++)
    {
        const float latitude = (90.0f - 180.0f * y / (rows - 1)) * pi / 180.0f;
        for (int x = 0; x < columns; x++)
            east[index(x, y)] = angularVelocity * radius * std::cos(latitude);
    }
    const auto rotation = climateflow::advectParticleRk2(
        columns, rows, east, south, 30.0f, 5.0f, 10800.0f, config);
    expect(rotation.remainedInDomain &&
            std::abs(rotation.points.back().x - 2.0f) < 0.02f &&
            std::abs(rotation.points.back().y - 5.0f) < 1.0e-5f,
        "solid-body rotation must wrap the periodic seam with the analytical phase speed");

    const float vortexCentreX = columns * 0.5f;
    const float vortexCentreY = rows * 0.5f;
    const float vortexAngularVelocity = 2.0f * pi / 86400.0f;
    const float meridionalCellHeight = pi * radius / (rows - 1);
    for (int y = 0; y < rows; y++)
    {
        const float latitude = (90.0f - 180.0f * y / (rows - 1)) * pi / 180.0f;
        const float zonalCellWidth = 2.0f * pi * radius *
            std::max(0.02f, std::abs(std::cos(latitude))) / columns;
        for (int x = 0; x < columns; x++)
        {
            east[index(x, y)] = -vortexAngularVelocity *
                (static_cast<float>(y) - vortexCentreY) * zonalCellWidth;
            south[index(x, y)] = vortexAngularVelocity *
                (static_cast<float>(x) - vortexCentreX) * meridionalCellHeight;
        }
    }
    const auto vortex = climateflow::advectParticleRk2(
        columns, rows, east, south,
        vortexCentreX + 3.0f, vortexCentreY, 21600.0f, config);
    expect(vortex.remainedInDomain &&
            std::abs(vortex.points.back().x - vortexCentreX) < 0.2f &&
            std::abs(vortex.points.back().y - (vortexCentreY + 3.0f)) < 0.2f,
        "adaptive midpoint integration must follow a known vortex through a quarter turn");

    std::fill(east.begin(), east.end(), 0.0f);
    std::fill(south.begin(), south.end(), -100.0f);
    const auto polar = climateflow::advectParticleRk2(
        columns, rows, east, south, 5.0f, 0.0f, 10800.0f, config);
    expect(!polar.remainedInDomain && polar.points.front().y == 0.0f,
        "particle integration must stop cleanly at a polar boundary");

    return failures == 0 ? 0 : 1;
}
