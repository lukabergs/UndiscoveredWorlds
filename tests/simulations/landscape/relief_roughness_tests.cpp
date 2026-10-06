// Landscape.ReliefRoughness*: the seeded relief perturbation on synthetic grids
// (docs/architecture/erosion-hydrology-plan.md, Phase C, step C3).

#include "landscape_fixture.hpp"
#include "simulations/landscape/relief_roughness.hpp"
#include "simulations/landscape/sphere_metric.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace
{
constexpr double radius = 6371000.0;

landscape::ReliefRoughness defaults(std::uint64_t seed)
{
    landscape::ReliefRoughness roughness;
    roughness.amplitude_m = 30.0;
    roughness.wavelength_km = 2000.0;
    roughness.octaves = 6;
    roughness.persistence = 0.6;
    roughness.seed = seed;
    return roughness;
}

// The perturbation of an all-land grid at 0 m (sea level far below).
std::vector<float> perturbation(const landscape::SphereMetric& metric, const landscape::ReliefRoughness& roughness)
{
    const std::vector<uint8_t> sea(metric.size(), 0);
    std::vector<float> elevation(metric.size(), 0.0f);
    landscape::addreliefroughness(metric, roughness, sea, -10000.0f, elevation);
    return elevation;
}
}

TEST_F(Landscape, ReliefRoughnessZeroAmplitudeChangesNothing)
{
    const landscape::SphereMetric metric(64, 32, radius);
    landscape::ReliefRoughness roughness = defaults(1);
    roughness.amplitude_m = 0.0;
    std::vector<uint8_t> sea(metric.size(), 0);
    std::vector<float> elevation(metric.size());
    for (std::size_t cell = 0; cell < metric.size(); cell++)
        elevation[cell] = static_cast<float>(cell % 97) * 3.5f;
    const std::vector<float> before = elevation;
    landscape::addreliefroughness(metric, roughness, sea, 0.0f, elevation);
    EXPECT_EQ(elevation, before);
}

// Sea cells keep their elevation. A lowering is the untapered perturbation (the
// same noise on land far above the sea) times min(1, (old - floor) / amplitude),
// floor = sealevel + 1 m, so low land is not flattened onto the floor: no cell ends
// exactly at it here (the backstop acts only beyond the amplitude). A land cell
// that starts at or below the floor is not lowered.
TEST_F(Landscape, ReliefRoughnessSeaUnchangedCoastTapered)
{
    const landscape::SphereMetric metric(64, 32, radius);
    const float sealevel = 0.0f;
    const float floor = sealevel + 1.0f;
    std::vector<uint8_t> sea(metric.size());
    std::vector<float> elevation(metric.size());
    for (std::size_t cell = 0; cell < metric.size(); cell++)
    {
        const int x = metric.column(cell);
        sea[cell] = x < 16 ? 1 : 0;
        // Land from 0.5 m to 40 m, most of it within the amplitude of the sea level.
        elevation[cell] = sea[cell] != 0 ? -50.0f : 0.5f + static_cast<float>((x * 7 + metric.row(cell) * 3) % 40);
    }
    const std::vector<float> before = elevation;
    landscape::ReliefRoughness roughness = defaults(1);
    roughness.wavelength_km = 4000.0;  // two octaves resolved at 64 x 32: 4000 and 2000 km
    roughness.octaves = 2;
    landscape::addreliefroughness(metric, roughness, sea, sealevel, elevation);
    const std::vector<float> raw = perturbation(metric, roughness);
    int lowered = 0, atfloor = 0;
    for (std::size_t cell = 0; cell < metric.size(); cell++)
    {
        if (sea[cell] != 0)
        {
            EXPECT_EQ(elevation[cell], before[cell]) << "sea cell " << cell;
            continue;
        }
        const double old = before[cell];
        double delta = raw[cell];
        if (delta < 0.0)
            delta *= std::clamp((old - floor) / roughness.amplitude_m, 0.0, 1.0);
        EXPECT_NEAR(elevation[cell], old + delta, 1.0e-3) << "land cell " << cell;
        if (old > floor)
            EXPECT_GT(elevation[cell], floor) << "land cell " << cell;
        else
            EXPECT_GE(elevation[cell], before[cell]) << "land cell " << cell;
        lowered += elevation[cell] < before[cell] ? 1 : 0;
        atfloor += elevation[cell] == floor ? 1 : 0;
    }
    EXPECT_GT(lowered, 100);
    EXPECT_EQ(atfloor, 0);
}

// On 256 x 128 (all land), with the four octaves of 4000 to 500 km all resolved,
// the area-weighted mean is within 10 percent of the amplitude of zero (about 30
// features of the first octave on the sphere) and the largest magnitude within 30
// percent of the amplitude (measured at 30 m: means 0.17, -0.36 and 0.79 m, peaks
// 30.9, 25.4 and 24.9 m).
TEST_F(Landscape, ReliefRoughnessMeanAndPeak)
{
    const landscape::SphereMetric metric(256, 128, radius);
    for (std::uint64_t seed : { 1u, 2u, 3u })
    {
        landscape::ReliefRoughness roughness = defaults(seed);
        roughness.wavelength_km = 4000.0;
        roughness.octaves = 4;
        const std::vector<float> delta = perturbation(metric, roughness);
        double sum = 0.0, area = 0.0, peak = 0.0;
        for (std::size_t cell = 0; cell < metric.size(); cell++)
        {
            sum += delta[cell] * metric.area_m2(cell);
            area += metric.area_m2(cell);
            peak = std::max(peak, std::fabs(static_cast<double>(delta[cell])));
        }
        EXPECT_LT(std::fabs(sum / area), 0.1 * 30.0) << "seed " << seed;
        EXPECT_GT(peak, 0.7 * 30.0) << "seed " << seed;
        EXPECT_LT(peak, 1.3 * 30.0) << "seed " << seed;
    }
}

// Octaves below twice the north-south cell size get weight 0 and the scale stays
// that of the full set: on 64 x 32 (dy 625 km) only the 2000 km octave of the six
// is resolved, so the field is the one-octave field divided by
// sqrt(sum of persistence^2k over the six).
TEST_F(Landscape, ReliefRoughnessBandLimitedToGrid)
{
    const landscape::SphereMetric metric(64, 32, radius);
    const std::vector<float> six = perturbation(metric, defaults(4));
    landscape::ReliefRoughness single = defaults(4);
    single.octaves = 1;
    const std::vector<float> one = perturbation(metric, single);
    double squared = 0.0;
    for (int k = 0; k < 6; k++)
        squared += std::pow(0.36, k);
    for (std::size_t cell = 0; cell < metric.size(); cell++)
        EXPECT_NEAR(six[cell], one[cell] / std::sqrt(squared), 1.0e-4) << "cell " << cell;
}

// The value is a function of the position on the sphere up to the octaves the
// grid resolves: the centres of the 64 x 32 grid are centres of the 192 x 96 grid
// (column 3x + 1, row 3y + 1); both resolve the octaves of 4000 and 2000 km.
// Across the date line columns 0 and W - 1 differ no more than other neighbours.
TEST_F(Landscape, ReliefRoughnessIndependentOfGridWidth)
{
    const landscape::SphereMetric coarse(64, 32, radius), fine(192, 96, radius);
    landscape::ReliefRoughness roughness = defaults(5);
    roughness.wavelength_km = 4000.0;
    roughness.octaves = 2;
    const std::vector<float> a = perturbation(coarse, roughness), b = perturbation(fine, roughness);
    for (int y = 0; y < coarse.height; y++)
        for (int x = 0; x < coarse.width; x++)
            EXPECT_NEAR(a[coarse.index(x, y)], b[fine.index(3 * x + 1, 3 * y + 1)], 1.0e-3)
                << "column " << x << " row " << y;

    double seam = 0.0, interior = 0.0;
    for (int y = 0; y < fine.height; y++)
    {
        seam = std::max(seam, std::fabs(static_cast<double>(b[fine.index(0, y)] - b[fine.index(fine.width - 1, y)])));
        for (int x = 0; x + 1 < fine.width; x++)
            interior = std::max(interior, std::fabs(static_cast<double>(b[fine.index(x, y)] - b[fine.index(x + 1, y)])));
    }
    EXPECT_LE(seam, interior);
}

TEST_F(Landscape, ReliefRoughnessDeterministic)
{
    const landscape::SphereMetric metric(64, 32, radius);
    const std::vector<float> first = perturbation(metric, defaults(9));
    EXPECT_EQ(first, perturbation(metric, defaults(9)));
    const std::vector<float> other = perturbation(metric, defaults(10));
    std::size_t differing = 0;
    for (std::size_t cell = 0; cell < metric.size(); cell++)
        differing += first[cell] != other[cell] ? 1 : 0;
    EXPECT_GT(differing, metric.size() / 2);
}
