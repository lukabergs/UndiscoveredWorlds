// Landscape.Subgrid*: the sub-grid statistics of the global grid
// (docs/architecture/multiscale-plan.md, step S3.6): Perron's valley spacing calibrated to a
// site of Perron et al. 2009, Table 1, and the grid-scale steady-state relief held to the
// local relief.

#include "landscape_fixture.hpp"
#include "simulations/landscape/sphere_metric.hpp"
#include "simulations/landscape/subgrid_statistics.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace
{
constexpr double radius = 6371000.0;

// The erosion's grid-scale D (landscape.json, hillslope), c = 8 and the Napa Valley reference
// (terrain.json, subgrid).
landscape::SubgridParameters defaults()
{
    landscape::SubgridParameters parameters;
    parameters.diffusivity_km2_myr = 200.0;
    parameters.spacingcoefficient = 8.0;
    parameters.referencedoverk = 86.0;
    parameters.referencem = 0.35;
    parameters.referencerunoff_mm_yr = 400.0;
    parameters.referenceerodibility = 1.0;
    return parameters;
}

// The kernel's valley spacing (km) on a flat 64 x 32 all-land grid with uniform erodibility
// and runoff, at one mid-latitude cell (far below its cell size).
float kernelspacing_km(const landscape::SubgridParameters& parameters, float erodibility, float runoff_mm_yr)
{
    const landscape::SphereMetric metric(64, 32, radius);
    const std::size_t size = metric.size();
    const std::vector<float> elevation(size, 100.0f), erodibilityfield(size, erodibility),
        runoff(size, runoff_mm_yr), rate(size, 1.0f);
    const std::vector<std::uint8_t> sea(size, 0);
    std::vector<float> relief, spacing;
    landscape::subgridstatistics(metric, parameters, elevation, sea, erodibilityfield, runoff, rate, relief, spacing);
    return spacing[metric.index(20, 10)];
}
}

// Perron et al. 2009, Table 1, Napa Valley: D/K = 86 m^(2m+1) at m = 0.35 gives Lc = 13.7 m;
// at the reference runoff and erodibility the cell's D/K is the site's.
TEST_F(Landscape, SubgridCharacteristicLengthReproducesPerronNapaValley)
{
    const landscape::SubgridParameters parameters = defaults();
    const double length_m = 1000.0 * landscape::characteristiclength_km(parameters, 1.0, 400.0);
    EXPECT_NEAR(length_m, 13.7, 0.05);
    // lambda = 8 Lc = 110 m against the site's measured 128 m (lambda / Lc 9.3 there).
    EXPECT_NEAR(1000.0 * kernelspacing_km(parameters, 1.0f, 400.0f), 8.0 * length_m, 1.0e-3);
}

// Perron et al. 2009, Table 1, Dragon's Back: D/K = 12 m^(2m+1) at m = 0.42; the table gives
// Lc = 4.0 m and lambda = 30 m. The kernel, given that site as its reference, reproduces both
// (lambda with c = 8 against the table's lambda / Lc of 7.5 there).
TEST_F(Landscape, SubgridSpacingReproducesPerronDragonsBack)
{
    landscape::SubgridParameters parameters = defaults();
    parameters.referencedoverk = 12.0;
    parameters.referencem = 0.42;
    EXPECT_NEAR(1000.0 * landscape::characteristiclength_km(parameters, 1.0, 400.0), 4.0, 0.2);
    EXPECT_NEAR(1000.0 * kernelspacing_km(parameters, 1.0f, 400.0f), 30.0, 3.0);
}

// Per cell, K_A = K x e x r^m_ref: doubling the runoff multiplies lambda by
// 2^(-m/(2m+1)), doubling the erodibility by 2^(-1/(2m+1)); no runoff or no erodibility is
// no spacing.
TEST_F(Landscape, SubgridSpacingScalesWithRunoffAndErodibility)
{
    const landscape::SubgridParameters parameters = defaults();
    const double m = parameters.referencem;
    const double reference = kernelspacing_km(parameters, 1.0f, 400.0f);
    EXPECT_NEAR(kernelspacing_km(parameters, 1.0f, 800.0f) / reference, std::pow(2.0, -m / (2.0 * m + 1.0)), 1.0e-5);
    EXPECT_NEAR(kernelspacing_km(parameters, 2.0f, 400.0f) / reference, std::pow(2.0, -1.0 / (2.0 * m + 1.0)), 1.0e-5);
    EXPECT_EQ(kernelspacing_km(parameters, 1.0f, 0.0f), 0.0f);
    EXPECT_EQ(kernelspacing_km(parameters, 0.0f, 400.0f), 0.0f);
}

// E l^2 / (2 D): 100 m/Myr over 10 km at 200 km2/Myr is 25 m; a negative rate is no relief.
TEST_F(Landscape, SubgridHillslopeReliefIsTheParabola)
{
    const landscape::SubgridParameters parameters = defaults();
    EXPECT_NEAR(landscape::hillsloperelief_m(parameters, 100.0, 10.0), 25.0, 1.0e-9);
    EXPECT_NEAR(landscape::hillsloperelief_m(parameters, 1.0, 10.0), 0.25, 1.0e-9);
    EXPECT_EQ(landscape::hillsloperelief_m(parameters, -5.0, 10.0), 0.0);
}

// A 64 x 32 grid: a checkerboard of 0 and 500 m (local relief 500 m on every cell), sea in
// column 0, no runoff in column 1, a runoff of 1e-20 mm/yr in column 2 (a spacing far beyond
// the cell: stored as the cell size), 400 mm/yr elsewhere (8 x 13.7 m), 3 m/Myr everywhere.
// The relief's half-length is half the cell size (sqrt of its area: about 625 km at the
// equator), and the relief is the parabola's, at most 500 m (the cap binds near the equator,
// not near the poles).
TEST_F(Landscape, SubgridStatisticsOnASyntheticGrid)
{
    const landscape::SphereMetric metric(64, 32, radius);
    const landscape::SubgridParameters parameters = defaults();
    const std::size_t size = metric.size();
    std::vector<float> elevation(size), erodibility(size, 1.0f), runoff(size, 400.0f), rate(size, 3.0f);
    std::vector<std::uint8_t> sea(size, 0);
    for (int y = 0; y < metric.height; y++)
        for (int x = 0; x < metric.width; x++)
        {
            const std::size_t cell = metric.index(x, y);
            elevation[cell] = (x + y) % 2 == 0 ? 500.0f : 0.0f;
            sea[cell] = x == 0 ? 1 : 0;
            runoff[cell] = x == 1 ? 0.0f : x == 2 ? 1.0e-20f : 400.0f;
        }
    std::vector<float> relief, spacing;
    const std::size_t capped =
        landscape::subgridstatistics(metric, parameters, elevation, sea, erodibility, runoff, rate, relief, spacing);
    ASSERT_EQ(relief.size(), size);
    ASSERT_EQ(spacing.size(), size);
    EXPECT_EQ(capped, static_cast<std::size_t>(metric.height));

    const double expectedspacing = 8.0 * std::pow(86.0, 1.0 / 1.7) / 1000.0;
    int reliefcapped = 0, belowcap = 0;
    for (int y = 0; y < metric.height; y++)
        for (int x = 0; x < metric.width; x++)
        {
            const std::size_t cell = metric.index(x, y);
            if (x == 0)
            {
                EXPECT_EQ(relief[cell], 0.0f);
                EXPECT_EQ(spacing[cell], 0.0f);
                continue;
            }
            const double cellsize = std::sqrt(metric.area_m2(cell)) / 1000.0;
            const double cellspacing = x == 1 ? 0.0 : x == 2 ? cellsize : expectedspacing;
            EXPECT_NEAR(spacing[cell], cellspacing, 1.0e-6 * std::max(1.0, cellspacing));
            const double half = 0.5 * cellsize;
            const double expected = std::min(3.0 * half * half / 400.0, 500.0);
            EXPECT_NEAR(relief[cell], expected, 1.0e-3 * std::max(1.0, expected));
            (expected >= 500.0 ? reliefcapped : belowcap)++;
        }
    // Both regimes of the relief occur on this grid.
    EXPECT_GT(reliefcapped, 0);
    EXPECT_GT(belowcap, 0);

    // A flat land cell has no local relief, so no grid-scale relief either.
    std::fill(elevation.begin(), elevation.end(), 200.0f);
    landscape::subgridstatistics(metric, parameters, elevation, sea, erodibility, runoff, rate, relief, spacing);
    EXPECT_EQ(*std::max_element(relief.begin(), relief.end()), 0.0f);
    EXPECT_NEAR(spacing[metric.index(10, 16)], expectedspacing, 1.0e-6);

    std::vector<float> shortgrid(size - 1);
    EXPECT_THROW(landscape::subgridstatistics(metric, parameters, shortgrid, sea, erodibility, runoff, rate, relief, spacing),
                 std::invalid_argument);
}
