// Tests of landscape/spherical_noise.hpp (erosion plan Phase E, step E2a), suite Landscape.

#include <gtest/gtest.h>

#include "landscape_fixture.hpp"
#include "simulations/landscape/spherical_noise.hpp"

#include <cmath>
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
}
