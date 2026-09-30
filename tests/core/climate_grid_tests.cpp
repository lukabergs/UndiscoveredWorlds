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
        "finite-volume state must exclude the coordinate singularities");
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
    constexpr int invalidRows = 64 + 1;
    expect(climategrid::validGlobalGridDimensions(128, 64) &&
            !climategrid::validGlobalGridDimensions(128, invalidRows) &&
            !climategrid::validGlobalGridDimensions(127, 63),
        "global climate grids must have an even width and exactly half as many rows");
    expect(climategrid::makeSphericalGrid(128, invalidRows, 6371000.0).columns == 0,
        "the spherical-grid constructor must reject an extra latitude row");
    double bandMeasure = 0.0;
    for (int y = 0; y < 64; y++)
        bandMeasure += climategrid::latitudeBandMeasure(y, 64);
    expect(std::abs(bandMeasure - 2.0) < 1.0e-12,
        "cell-centred latitude bands must cover the sphere exactly");

    std::vector<float> fine(512 * 256, 7.25f);
    const auto internal = climategrid::remapField(
        512, 256, fine, 128, 64);
    const auto restored = climategrid::remapField(
        128, 64, internal, 2048, 1024);
    expect(internal.size() == 128 * 64 && restored.size() == 2048 * 1024,
        "explicit remapping must preserve the global grid contract");
    expect(std::all_of(internal.begin(), internal.end(), [](float value)
        { return std::abs(value - 7.25f) < 1.0e-6f; }) &&
        std::all_of(restored.begin(), restored.end(), [](float value)
        { return std::abs(value - 7.25f) < 1.0e-6f; }),
        "constant fields must survive both remap directions");
    const auto repeatedInternal = climategrid::remapField(
        512, 256, fine, 128, 64);
    expect(repeatedInternal == internal,
        "explicit grid remapping must be deterministic");

    std::vector<float> bands(512 * 256, 0.0f);
    for (int y = 0; y < 256; y++)
    {
        const float value = y < 128 ? 1.0f : -1.0f;
        std::fill_n(bands.begin() + static_cast<std::size_t>(y) * 512, 512, value);
    }
    const double sourceIntegral = climategrid::areaWeightedIntegral(
        512, 256, bands);
    const auto remappedBands = climategrid::remapField(
        512, 256, bands, 128, 64);
    const double destinationIntegral = climategrid::areaWeightedIntegral(
        128, 64, remappedBands);
    expect(std::abs(sourceIntegral - destinationIntegral) < 1.0e-5,
        "finite-volume remapping must preserve global integrals");

    // Noninteger refinement used by the 128-cell hydrology / 600-cell world.
    std::vector<float> polarRain(128 * 64, 0.0f);
    std::fill_n(polarRain.begin() + 128, 128, 100.0f);
    const auto rainPlan = climategrid::makeConservativeRemap(128, 64, 600, 300);
    std::vector<float> outputRain(600 * 300);
    for (int y = 0; y < 300; ++y)
        for (int x = 0; x < 600; ++x)
            outputRain[y * 600 + x] = climategrid::sampleConservative(rainPlan, x, y,
                [&](int xx, int yy) { return polarRain[yy * 128 + xx]; });
    expect(std::abs(climategrid::areaWeightedIntegral(600, 300, outputRain) /
        climategrid::areaWeightedIntegral(128, 64, polarRain) - 1.0) < 1.0e-6,
        "hydrology output remapping must preserve polar rainfall at noninteger refinement");
    expect(outputRain == climategrid::remapField(128, 64, polarRain, 600, 300),
        "production accessor remapping must agree with flat-field remapping");
    std::vector<float> smooth(32 * 16);
    for (int y = 0; y < 16; ++y)
        for (int x = 0; x < 32; ++x)
            smooth[y * 32 + x] = static_cast<float>(5.0 + 2.0 * std::cos(2.0 * pi * (x + 0.5) / 32) *
                std::cos(climategrid::latitudeCentreRadians(y, 16)));
    for (int outputColumns : {128, 150})
    {
        const auto reconstruction = climategrid::remapField(32, 16, smooth, outputColumns, outputColumns / 2);
        expect(std::abs(climategrid::areaWeightedIntegral(outputColumns, outputColumns / 2, reconstruction) /
            climategrid::areaWeightedIntegral(32, 16, smooth) - 1.0) < 1.0e-7,
            "linear reconstruction conserves spherical integrals at integer and noninteger refinement");
        const auto bounds = std::minmax_element(smooth.begin(), smooth.end());
        expect(std::all_of(reconstruction.begin(), reconstruction.end(), [&](float value) {
            return value >= *bounds.first && value <= *bounds.second; }),
            "limited reconstruction introduces no new extrema or negative rainfall");
        if (outputColumns == 128)
        {
            expect(reconstruction[32 * 128 + 32] != reconstruction[32 * 128 + 33],
                "smooth subcell variation must replace replicated 4x4 pixel blocks");
            const auto restoredMeans = climategrid::remapField(128, 64, reconstruction, 32, 16);
            for (std::size_t i = 0; i < smooth.size(); ++i)
                expect(std::abs(restoredMeans[i] - smooth[i]) < 1.0e-6,
                    "each reconstructed coarse cell retains its original area-weighted mean");
        }
    }
    // Output reconstruction may exchange neighbouring cell means, but must be
    // continuous, bounded and conservative, including the seam and polar caps.
    std::vector<float> peak(16 * 8, 0.0f);
    peak[3 * 16] = 100.0f;
    peak[0] = 40.0f;
    double previousJump = 0.0;
    for (int width : {128, 256, 300})
    {
        const auto image = climategrid::remapSmoothField(16, 8, peak, width, width / 2);
        expect(std::abs(climategrid::areaWeightedIntegral(width, width / 2, image) /
            climategrid::areaWeightedIntegral(16, 8, peak) - 1.0) < 1.0e-7,
            "continuous reconstruction conserves global mass for arbitrary refinement");
        double jump = 0.0;
        for (int y = 0; y < width / 2; ++y)
            for (int x = 0; x < width; ++x)
            {
                const float value = image[y * width + x];
                expect(value >= 0.0f && value <= 100.0f, "continuous reconstruction stays in donor bounds");
                jump = std::max(jump, static_cast<double>(std::abs(value - image[y * width + (x + 1) % width])));
                if (y > 0) jump = std::max(jump, static_cast<double>(std::abs(value - image[(y - 1) * width + x])));
            }
        if (width == 128) expect(jump < 25.0, "isolated wet cells have no discontinuous faces");
        if (width == 256) expect(jump < previousJump * 0.6,
            "adjacent-pixel jumps shrink under refinement rather than retaining coarse face discontinuities");
        previousJump = jump;
        const auto constant = climategrid::remapSmoothField(16, 8, std::vector<float>(128, 7.25f), width, width / 2);
        expect(std::all_of(constant.begin(), constant.end(), [](float v) { return std::abs(v - 7.25f) < 1.0e-6f; }),
            "continuous reconstruction preserves constants, including polar boundaries");
    }
    expect(climategrid::remapSmoothField(16, 8, peak, 16, 8) == peak,
        "native-resolution fields are not smoothed");
    expect(climategrid::remapSmoothField(16, 8, peak, 8, 4) == climategrid::remapField(16, 8, peak, 8, 4),
        "downsampling retains exact finite-volume restriction");
    std::vector<float> latitudeProfile(8, 0.0f);
    latitudeProfile[0] = 1.0f;
    const auto fineProfile = climategrid::remapSmoothField(1, 8, latitudeProfile, 1, 128);
    expect(std::abs(climategrid::areaWeightedIntegral(1, 128, fineProfile) /
        climategrid::areaWeightedIntegral(1, 8, latitudeProfile) - 1.0) < 1.0e-7,
        "one-dimensional energy profiles conserve polar integrals during reconstruction");
    return failures == 0 ? 0 : 1;
}
