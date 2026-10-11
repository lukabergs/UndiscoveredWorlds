// Tests of landscape/spherical_noise.hpp (erosion plan Phase E, step E2a), suite Landscape.

#include <gtest/gtest.h>

#include "landscape_fixture.hpp"
#include "core/deterministic_random.hpp"
#include "simulations/landscape/landscape_parameters.hpp"
#include "simulations/landscape/opensimplex2s.hpp"
#include "simulations/landscape/spherical_noise.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <stdexcept>

namespace
{
constexpr double pi = 3.14159265358979323846;

landscape::FractalNoise noise(double shortest_km, double longest_km)
{
    landscape::FractalNoise result;
    result.seed = 20261008;
    result.radius_km = 6371.0;
    result.wavelength_km = 2000.0;  // octaves of 2000, 1000, 500, 250, 125 and 62.5 km
    result.octaves = 6;
    result.persistence = 0.6;
    result.shortest_km = shortest_km;
    result.longest_km = longest_km;
    return result;
}

double valueat(const landscape::FractalNoiseField& field, double latitude, double longitude)
{
    double px = 0.0, py = 0.0, pz = 0.0;
    landscape::unitvector(latitude, longitude, px, py, pz);
    return field.sum(px, py, pz);
}

// The kernels and rotations of S3.16.
constexpr landscape::NoiseKernel kernels[] = { landscape::NoiseKernel::perlin, landscape::NoiseKernel::opensimplex2s };

landscape::FractalNoise noisewith(landscape::NoiseKernel kernel, int rotate, double shortest_km, double longest_km)
{
    landscape::FractalNoise result = noise(shortest_km, longest_km);
    result.kernel = kernel;
    result.rotateoctaves = rotate;
    return result;
}

double kernelvalue(landscape::NoiseKernel kernel, std::uint64_t seed, double x, double y, double z)
{
    return kernel == landscape::NoiseKernel::opensimplex2s ? landscape::opensimplex2s(seed, x, y, z)
                                                           : landscape::gradientnoise(seed, x, y, z);
}

double unitnumber(std::uint64_t hash) { return static_cast<double>(hash >> 11) * (1.0 / 9007199254740992.0); }

// The i-th point of a hashed sequence uniform on the unit sphere.
void spherepoint(std::uint64_t i, double& px, double& py, double& pz)
{
    const double z = 2.0 * unitnumber(deterministiccombine(20261011u, i, 0)) - 1.0;
    const double phi = 2.0 * pi * unitnumber(deterministiccombine(20261011u, i, 1));
    const double r = std::sqrt(std::max(0.0, 1.0 - z * z));
    px = r * std::cos(phi);
    py = r * std::sin(phi);
    pz = z;
}
}

TEST_F(Landscape, SphericalNoiseIsIndependentOfTheCallOrder)
{
    const landscape::FractalNoiseField field(noise(0.0, 0.0));
    const double first = valueat(field, 0.4, 1.3);
    double others = 0.0;
    for (int i = 0; i < 50; i++)
        others += valueat(field, -1.2 + 0.05 * i, 0.13 * i);
    const double again = valueat(field, 0.4, 1.3);
    EXPECT_EQ(first, again);
    // A second field built from the same parameters, evaluated in another order.
    const landscape::FractalNoiseField other(noise(0.0, 0.0));
    EXPECT_NE(valueat(other, -0.7, 4.0), first);
    EXPECT_EQ(valueat(other, 0.4, 1.3), first);
    EXPECT_NE(first, 0.0);
    EXPECT_TRUE(std::isfinite(others));
    // The lattice noise itself has no state either.
    const double lattice = landscape::gradientnoise(7, 1.25, -3.5, 0.75);
    landscape::gradientnoise(7, 9.0, 9.5, 10.25);
    EXPECT_EQ(landscape::gradientnoise(7, 1.25, -3.5, 0.75), lattice);
    EXPECT_EQ(landscape::gradientnoise(7, 2.0, -3.0, 5.0), 0.0);  // zero at a lattice point
}

TEST_F(Landscape, SphericalNoiseResolvesTheOctavesInsideItsBand)
{
    EXPECT_EQ(landscape::FractalNoiseField(noise(0.0, 0.0)).resolved(), 6);
    EXPECT_EQ(landscape::FractalNoiseField(noise(100.0, 0.0)).resolved(), 5);   // 62.5 km is below
    EXPECT_EQ(landscape::FractalNoiseField(noise(100.0, 600.0)).resolved(), 3); // 500, 250, 125 km
    EXPECT_EQ(landscape::FractalNoiseField(noise(125.0, 500.0)).resolved(), 2); // 250, 125: the upper end is outside
    EXPECT_EQ(landscape::FractalNoiseField(noise(3000.0, 0.0)).resolved(), 0);
    // Bands that meet at an octave's wavelength count it once: [0, 125) holds 62.5,
    // [125, no limit) the other five.
    EXPECT_EQ(landscape::FractalNoiseField(noise(0.0, 125.0)).resolved(), 1);
    EXPECT_EQ(landscape::FractalNoiseField(noise(125.0, 0.0)).resolved(), 5);

    // The normalisation is that of the full octave set, whatever the band.
    double squared = 0.0;
    for (int o = 0; o < 6; o++)
        squared += std::pow(0.6, o) * std::pow(0.6, o);
    const landscape::FractalNoiseField banded(noise(100.0, 600.0));
    EXPECT_DOUBLE_EQ(banded.squaredweights(), squared);
    EXPECT_DOUBLE_EQ(banded.normalisation(), landscape::octavedeviation * std::sqrt(squared));
    EXPECT_EQ(banded.normalisation(), landscape::FractalNoiseField(noise(0.0, 0.0)).normalisation());

    // An octave outside the band contributes nothing: the banded sum equals the sum
    // of the full set minus the octaves outside it, built as bands of their own.
    const landscape::FractalNoiseField low(noise(0.0, 0.0));
    const landscape::FractalNoiseField outsidelong(noise(600.0, 0.0));
    const landscape::FractalNoiseField outsideshort(noise(0.0, 100.0));
    const double whole = valueat(low, 0.3, 2.0);
    const double parts = valueat(banded, 0.3, 2.0) + valueat(outsidelong, 0.3, 2.0) + valueat(outsideshort, 0.3, 2.0);
    EXPECT_NEAR(parts, whole, 1e-12);
    // The same with bands that meet exactly at an octave (125 km): no octave twice.
    const double meeting = valueat(landscape::FractalNoiseField(noise(0.0, 125.0)), 0.3, 2.0) +
                           valueat(landscape::FractalNoiseField(noise(125.0, 0.0)), 0.3, 2.0);
    EXPECT_NEAR(meeting, whole, 1e-12);
}

TEST_F(Landscape, SphericalNoiseHasNoSeamAtTheDateLine)
{
    const landscape::FractalNoiseField field(noise(0.0, 0.0));
    for (double latitude : { -1.4, -0.6, 0.0, 0.5, 1.2 })
        EXPECT_NEAR(valueat(field, latitude, 0.0), valueat(field, latitude, 2.0 * pi), 1e-12) << "latitude " << latitude;
}

TEST_F(Landscape, SphericalNoiseThrowsOnBadInput)
{
    landscape::FractalNoise bad = noise(0.0, 0.0);
    bad.wavelength_km = 0.0;
    EXPECT_THROW(landscape::FractalNoiseField{ bad }, std::invalid_argument);
    bad = noise(0.0, 0.0);
    bad.octaves = 0;
    EXPECT_THROW(landscape::FractalNoiseField{ bad }, std::invalid_argument);
    bad = noisewith(static_cast<landscape::NoiseKernel>(2), 0, 0.0, 0.0);
    EXPECT_THROW(landscape::FractalNoiseField{ bad }, std::invalid_argument);
    bad = noisewith(landscape::NoiseKernel::perlin, 2, 0.0, 0.0);
    EXPECT_THROW(landscape::FractalNoiseField{ bad }, std::invalid_argument);
}

// S3.16: both kernels lie in [-1, 1] with zero mean, sampled at 1e5 points of a
// sphere 64 lattice units in radius (about 5e4 lattice cells on its surface).
TEST_F(Landscape, SphericalNoiseKernelsStayInRangeWithZeroMean)
{
    for (const landscape::NoiseKernel kernel : kernels)
    {
        double lowest = 0.0, highest = 0.0, total = 0.0, squares = 0.0;
        const int samples = 100000;
        for (int i = 0; i < samples; i++)
        {
            double px = 0.0, py = 0.0, pz = 0.0;
            spherepoint(static_cast<std::uint64_t>(i), px, py, pz);
            const double value = kernelvalue(kernel, 7, 64.0 * px, 64.0 * py, 64.0 * pz);
            lowest = std::min(lowest, value);
            highest = std::max(highest, value);
            total += value;
            squares += value * value;
        }
        const double mean = total / samples;
        const double deviation = std::sqrt(squares / samples - mean * mean);
        EXPECT_GE(lowest, -1.0) << "kernel " << static_cast<int>(kernel);
        EXPECT_LE(highest, 1.0) << "kernel " << static_cast<int>(kernel);
        EXPECT_LT(std::abs(mean), 0.01) << "kernel " << static_cast<int>(kernel);
        EXPECT_GT(deviation, 0.1) << "kernel " << static_cast<int>(kernel);
    }
}

// S3.16: both kernels are continuous: a step of 1e-6 lattice units along an axis
// changes the value by at most slopebound x the step, at 1e4 points (measured
// steepest: Perlin 2.72, OpenSimplex2S 4.64 per lattice unit).
TEST_F(Landscape, SphericalNoiseKernelsAreContinuous)
{
    const double step = 1e-6;
    const double slopebound = 8.0;
    for (const landscape::NoiseKernel kernel : kernels)
    {
        double steepest = 0.0;
        for (int i = 0; i < 10000; i++)
        {
            const std::uint64_t index = static_cast<std::uint64_t>(i);
            const double x = 100.0 * unitnumber(deterministiccombine(5u, index, 0)) - 50.0;
            const double y = 100.0 * unitnumber(deterministiccombine(5u, index, 1)) - 50.0;
            const double z = 100.0 * unitnumber(deterministiccombine(5u, index, 2)) - 50.0;
            const double here = kernelvalue(kernel, 11, x, y, z);
            steepest = std::max(steepest, std::abs(kernelvalue(kernel, 11, x + step, y, z) - here) / step);
            steepest = std::max(steepest, std::abs(kernelvalue(kernel, 11, x, y + step, z) - here) / step);
            steepest = std::max(steepest, std::abs(kernelvalue(kernel, 11, x, y, z + step) - here) / step);
        }
        EXPECT_LE(steepest, slopebound) << "kernel " << static_cast<int>(kernel);
    }
}

// S3.16: the deviation of one octave on the sphere is the stored constant of its
// kernel, measured as octavedeviation was (1024 x 512 cell centres, area-weighted,
// three seeds, two wavelengths).
TEST_F(Landscape, SphericalNoiseDeviationsAreTheMeasuredOnes)
{
    const int width = 1024, height = 512;
    for (const landscape::NoiseKernel kernel : kernels)
    {
        const double expected = kernel == landscape::NoiseKernel::opensimplex2s ? landscape::opensimplexdeviation
                                                                                : landscape::octavedeviation;
        for (const std::uint64_t seed : { 1u, 2u, 3u })
            for (const double wavelength : { 2000.0, 500.0 })
            {
                landscape::FractalNoise one = noisewith(kernel, 0, 0.0, 0.0);
                one.seed = seed;
                one.wavelength_km = wavelength;
                one.octaves = 1;
                const landscape::FractalNoiseField field(one);
                double area = 0.0, total = 0.0, squares = 0.0;
                for (int row = 0; row < height; row++)
                {
                    const double latitude = pi / 2.0 - pi * (row + 0.5) / height;
                    const double weight = std::cos(latitude);
                    for (int column = 0; column < width; column++)
                    {
                        const double value = valueat(field, latitude, 2.0 * pi * (column + 0.5) / width);
                        area += weight;
                        total += weight * value;
                        squares += weight * value * value;
                    }
                }
                const double mean = total / area;
                const double deviation = std::sqrt(squares / area - mean * mean);
                EXPECT_NEAR(deviation, expected, 0.006)
                    << "kernel " << static_cast<int>(kernel) << " seed " << seed << " wavelength " << wavelength;
            }
    }
}

// S3.16: an octave's rotation is a pure function of its seed and a rotation
// (orthonormal, determinant +1); rotating keeps the band and the normalisation.
TEST_F(Landscape, SphericalNoiseOctaveRotationIsOrthonormal)
{
    for (int k = 0; k < 200; k++)
    {
        const std::uint64_t seed = deterministiccombine(5u, k);
        const std::array<double, 9> m = landscape::octaverotation(seed);
        EXPECT_EQ(m, landscape::octaverotation(seed));
        EXPECT_NE(m, landscape::octaverotation(seed + 1));
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                double dot = 0.0;
                for (int c = 0; c < 3; c++)
                    dot += m[3 * i + c] * m[3 * j + c];
                EXPECT_NEAR(dot, i == j ? 1.0 : 0.0, 1e-12) << "seed " << seed << " rows " << i << ", " << j;
            }
        const double determinant = m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6]) +
                                   m[2] * (m[3] * m[7] - m[4] * m[6]);
        EXPECT_NEAR(determinant, 1.0, 1e-12) << "seed " << seed;
    }

    const landscape::FractalNoiseField plain(noisewith(landscape::NoiseKernel::perlin, 0, 100.0, 600.0));
    const landscape::FractalNoiseField rotated(noisewith(landscape::NoiseKernel::perlin, 1, 100.0, 600.0));
    EXPECT_EQ(rotated.resolved(), plain.resolved());
    EXPECT_EQ(rotated.normalisation(), plain.normalisation());
    EXPECT_NE(valueat(rotated, 0.3, 2.0), valueat(plain, 0.3, 2.0));
    // The defaults (landscape.noise.kernel 0, rotate_octaves 0) are Perlin without rotation.
    const landscape::FractalNoiseField configured(noise(100.0, 600.0));
    EXPECT_EQ(configured.normalisation(), plain.normalisation());
    for (int i = 0; i < 20; i++)
        EXPECT_EQ(valueat(configured, -1.3 + 0.13 * i, 0.31 * i), valueat(plain, -1.3 + 0.13 * i, 0.31 * i));
}

// S3.16: a FractalNoise left at "configured" takes the kernel and the rotation from
// landscape.noise.kernel and landscape.noise.rotate_octaves when its field is built.
TEST_F(Landscape, SphericalNoiseConfiguredReadsTheParameters)
{
    const int kernel = tuning::landscape::noise::kernel;
    const int rotate = tuning::landscape::noise::rotateOctaves;
    tuning::landscape::noise::kernel = 1;
    tuning::landscape::noise::rotateOctaves = 1;
    const landscape::FractalNoiseField configured(noise(0.0, 0.0));
    tuning::landscape::noise::kernel = kernel;
    tuning::landscape::noise::rotateOctaves = rotate;
    const landscape::FractalNoiseField explicit_(noisewith(landscape::NoiseKernel::opensimplex2s, 1, 0.0, 0.0));
    const landscape::FractalNoiseField perlin(noisewith(landscape::NoiseKernel::perlin, 0, 0.0, 0.0));
    EXPECT_EQ(configured.normalisation(), explicit_.normalisation());
    EXPECT_EQ(configured.deviation(), landscape::opensimplexdeviation);
    for (int i = 0; i < 10; i++)
    {
        const double latitude = -1.3 + 0.27 * i, longitude = 0.61 * i;
        EXPECT_EQ(valueat(configured, latitude, longitude), valueat(explicit_, latitude, longitude)) << "point " << i;
        EXPECT_NE(valueat(configured, latitude, longitude), valueat(perlin, latitude, longitude)) << "point " << i;
    }
}

// S3.16: every kernel and rotation is seamless (the same unit vector reached from
// another longitude, or by another field built from the same parameters after other
// points, gives the same value) and keeps the dropped-octave rule.
TEST_F(Landscape, SphericalNoiseIsSeamlessForEveryKernel)
{
    for (const landscape::NoiseKernel kernel : kernels)
        for (const int rotate : { 0, 1 })
        {
            const landscape::FractalNoiseField field(noisewith(kernel, rotate, 0.0, 0.0));
            for (double latitude : { -1.4, -0.6, 0.0, 0.5, 1.2 })
                EXPECT_NEAR(valueat(field, latitude, 0.0), valueat(field, latitude, 2.0 * pi), 1e-12)
                    << "kernel " << static_cast<int>(kernel) << " rotate " << rotate << " latitude " << latitude;
            const double first = valueat(field, 0.4, 1.3);
            const landscape::FractalNoiseField other(noisewith(kernel, rotate, 0.0, 0.0));
            EXPECT_NE(valueat(other, -0.7, 4.0), first);
            EXPECT_EQ(valueat(other, 0.4, 1.3), first);

            const double whole = valueat(field, 0.3, 2.0);
            const double parts = valueat(landscape::FractalNoiseField(noisewith(kernel, rotate, 100.0, 600.0)), 0.3, 2.0) +
                                 valueat(landscape::FractalNoiseField(noisewith(kernel, rotate, 600.0, 0.0)), 0.3, 2.0) +
                                 valueat(landscape::FractalNoiseField(noisewith(kernel, rotate, 0.0, 100.0)), 0.3, 2.0);
            EXPECT_NEAR(parts, whole, 1e-12) << "kernel " << static_cast<int>(kernel) << " rotate " << rotate;
        }
}
