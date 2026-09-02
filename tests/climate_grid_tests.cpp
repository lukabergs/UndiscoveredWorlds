#include "climate_grid.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
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
    constexpr double pi = 3.14159265358979323846;
    const auto grid = climategrid::makeSphericalGrid(128, 64, 6371000.0);
    expect(grid.columns == 128 && grid.rows == 64,
        "the internal climate grid must be W by W/2");
    expect(grid.latitudeCentresRadians.front() < pi * 0.5 &&
            grid.latitudeCentresRadians.back() > -pi * 0.5,
        "finite-volume state must be cell-centred rather than pole-inclusive");
    expect(std::abs(grid.latitudeCentresRadians.front() +
            grid.latitudeCentresRadians.back()) < 1.0e-12,
        "cell centres must be equator-symmetric");
    expect(grid.wrapColumn(-1) == 127 && grid.wrapColumn(128) == 0 &&
            grid.index(-1, 4) == grid.index(127, 4),
        "spherical indexing must wrap continuously across the periodic seam");

    double sphereArea = 0.0;
    for (double rowArea : grid.cellAreasSquareMetres)
        sphereArea += rowArea * grid.columns;
    const double analyticalArea = 4.0 * pi * grid.radiusMetres * grid.radiusMetres;
    expect(std::abs(sphereArea - analyticalArea) / analyticalArea < 1.0e-12,
        "spherical cell areas must integrate to the sphere area");
    expect(grid.northFaceLengthsMetres.front() < 1.0e-8 &&
            grid.southFaceLengthsMetres.back() < 1.0e-8,
        "polar boundary faces must have zero length");
    expect(std::abs(climategrid::latitudeCentreRadians(
            0, 257, climategrid::LatitudeLayout::poleInclusive) - pi * 0.5) < 1.0e-12 &&
            std::abs(climategrid::latitudeCentreRadians(
                256, 257, climategrid::LatitudeLayout::poleInclusive) + pi * 0.5) < 1.0e-12,
        "pole-inclusive raster coordinates must remain an explicit I/O convention");
    double bandMeasure = 0.0;
    for (int y = 0; y < 257; y++)
        bandMeasure += climategrid::latitudeBandMeasure(
            y, 257, climategrid::LatitudeLayout::poleInclusive);
    expect(std::abs(bandMeasure - 2.0) < 1.0e-12,
        "pole-inclusive raster bands must still cover the sphere exactly");

    std::vector<float> poleInclusive(512 * 257, 7.25f);
    const auto internal = climategrid::remapField(
        512, 257, climategrid::LatitudeLayout::poleInclusive, poleInclusive,
        128, 64, climategrid::LatitudeLayout::cellCentred);
    const auto restored = climategrid::remapField(
        128, 64, climategrid::LatitudeLayout::cellCentred, internal,
        2048, 1025, climategrid::LatitudeLayout::poleInclusive);
    expect(internal.size() == 128 * 64 && restored.size() == 2048 * 1025,
        "explicit remapping must support requested pole-inclusive outputs");
    expect(std::all_of(internal.begin(), internal.end(), [](float value)
        { return std::abs(value - 7.25f) < 1.0e-6f; }) &&
        std::all_of(restored.begin(), restored.end(), [](float value)
        { return std::abs(value - 7.25f) < 1.0e-6f; }),
        "constant fields must survive both remap directions");
    const auto repeatedInternal = climategrid::remapField(
        512, 257, climategrid::LatitudeLayout::poleInclusive, poleInclusive,
        128, 64, climategrid::LatitudeLayout::cellCentred);
    expect(repeatedInternal == internal,
        "explicit grid remapping must be deterministic");

    std::vector<float> bands(512 * 257, 0.0f);
    for (int y = 0; y < 257; y++)
    {
        const float value = y < 128 ? 1.0f : -1.0f;
        std::fill_n(bands.begin() + static_cast<std::size_t>(y) * 512, 512, value);
    }
    const double sourceIntegral = climategrid::areaWeightedIntegral(
        512, 257, climategrid::LatitudeLayout::poleInclusive, bands);
    const auto remappedBands = climategrid::remapField(
        512, 257, climategrid::LatitudeLayout::poleInclusive, bands,
        128, 64, climategrid::LatitudeLayout::cellCentred);
    const double destinationIntegral = climategrid::areaWeightedIntegral(
        128, 64, climategrid::LatitudeLayout::cellCentred, remappedBands);
    expect(std::abs(sourceIntegral - destinationIntegral) < 1.0e-5,
        "finite-volume remapping must preserve global integrals");

    return failures == 0 ? 0 : 1;
}
