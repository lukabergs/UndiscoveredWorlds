#include "simulations/climate/climate_flow.hpp"

#include <cmath>
#include <iostream>
#include <vector>
#include <gtest/gtest.h>

namespace
{

void expect(bool condition, const char* message)
{
    EXPECT_TRUE(condition) << message;
}
}

namespace
{

// The 32 x 16 particle grid of the flow tests with calm winds and the trace config.
class ClimateFlow : public ::testing::Test
{
protected:
    static constexpr int columns = 32;
    static constexpr int rows = 16;
    static constexpr float radius = 6371000.0f;
    static constexpr float pi = 3.14159265358979323846f;
    static constexpr std::size_t cellCount = columns * rows;
    static std::size_t index(int x, int y)
    {
        return static_cast<std::size_t>(y) * columns + x;
    }
    std::vector<float> east = std::vector<float>(cellCount, 0.0f);
    std::vector<float> south = std::vector<float>(cellCount, 0.0f);
    static constexpr int equator = rows / 2;
    const float equatorLatitude =
        (90.0f - 180.0f * (equator + 0.5f) / rows) * pi / 180.0f;
    const float equatorialCellWidth =
        2.0f * pi * radius * std::cos(equatorLatitude) / columns;
    climateflow::ParticleTraceConfig config;

    ClimateFlow()
    {
        config.planetRadiusMetres = radius;
        config.maximumSubstepCells = 0.2f;
    }
};

}

TEST_F(ClimateFlow, UniformFlow)
{
    std::fill(east.begin(), east.end(), equatorialCellWidth / 3600.0f);
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
}

TEST_F(ClimateFlow, ZeroWind)
{
    const auto stationary = climateflow::advectParticleRk2(
        columns, rows, east, south, 7.0f, static_cast<float>(equator), 3600.0f, config);
    expect(stationary.remainedInDomain && stationary.points.size() == 2 &&
            stationary.points.back().x == 7.0f &&
            stationary.points.back().elapsedSeconds == 3600.0f,
        "a valid zero-wind trace must advance time without inventing line length");
}

TEST_F(ClimateFlow, SolidBodyRotation)
{
    const float angularVelocity = 2.0f * pi / 86400.0f;
    for (int y = 0; y < rows; y++)
    {
        const float latitude =
            (90.0f - 180.0f * (y + 0.5f) / rows) * pi / 180.0f;
        for (int x = 0; x < columns; x++)
            east[index(x, y)] = angularVelocity * radius * std::cos(latitude);
    }
    const auto rotation = climateflow::advectParticleRk2(
        columns, rows, east, south, 30.0f, 5.0f, 10800.0f, config);
    expect(rotation.remainedInDomain &&
            std::abs(rotation.points.back().x - 2.0f) < 0.02f &&
            std::abs(rotation.points.back().y - 5.0f) < 1.0e-5f,
        "solid-body rotation must wrap the periodic seam with the analytical phase speed");
}

TEST_F(ClimateFlow, Vortex)
{
    const float vortexCentreX = columns * 0.5f;
    const float vortexCentreY = rows * 0.5f;
    const float vortexAngularVelocity = 2.0f * pi / 86400.0f;
    const float meridionalCellHeight = pi * radius / rows;
    for (int y = 0; y < rows; y++)
    {
        const float latitude =
            (90.0f - 180.0f * (y + 0.5f) / rows) * pi / 180.0f;
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
}

TEST_F(ClimateFlow, PolarBoundary)
{
    std::fill(south.begin(), south.end(), -100.0f);
    const auto polar = climateflow::advectParticleRk2(
        columns, rows, east, south, 5.0f, 0.0f, 10800.0f, config);
    expect(!polar.remainedInDomain && polar.points.front().y == 0.0f,
        "particle integration must stop cleanly at a polar boundary");
}
