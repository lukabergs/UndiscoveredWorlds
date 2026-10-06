// Landscape.Isostasy*: the local isostatic response of the erosion run
// (docs/architecture/erosion-hydrology-plan.md, Phase C, step C8).

#include "landscape_fixture.hpp"
#include "simulations/landscape/isostasy.hpp"
#include "simulations/landscape/sphere_metric.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace
{
constexpr double radius = 6371000.0;
// Deliberately unequal fractions, so that the tests tell the two parts apart: the
// engine's continental crust and sediment over its mantle (tectonics::columns). The
// default loading fraction equals the rebound fraction (landscape.isostasy.* note).
constexpr double rebound = 2800.0 / 3300.0;
constexpr double loading = 2700.0 / 3300.0;

landscape::IsostasyParameters defaults()
{
    landscape::IsostasyParameters parameters;
    parameters.reboundfraction = rebound;
    parameters.loadingfraction = loading;
    return parameters;
}
}

// A plateau at 2000 m denuded uniformly by 100 m (the erosion step left it at
// 1900 m) rises back by the fraction of the removed thickness on every cell.
TEST_F(Landscape, IsostasyDenudedPlateauReboundsByFraction)
{
    const std::size_t size = 200;
    const std::vector<float> erosion(size, 100.0f), deposition(size, 0.0f);
    std::vector<float> elevation(size, 1900.0f);
    landscape::applyisostasy(defaults(), erosion, deposition, elevation);
    const float expected = static_cast<float>(1900.0 + rebound * 100.0);
    for (float height : elevation)
        EXPECT_FLOAT_EQ(height, expected);
    // The net lowering is the removed thickness times 1 - fraction (about 15 m).
    EXPECT_NEAR(2000.0 - elevation[0], 100.0 * (1.0 - rebound), 1.0e-3);
}

// A deposit of 50 m (the step raised a basin from -200 m to -150 m) sinks by the
// loading fraction of its thickness, on sea cells as on land.
TEST_F(Landscape, IsostasyDepositSinksByLoadingFraction)
{
    const std::size_t size = 50;
    const std::vector<float> erosion(size, 0.0f), deposition(size, 50.0f);
    std::vector<float> elevation(size, -150.0f);
    landscape::applyisostasy(defaults(), erosion, deposition, elevation);
    const float expected = static_cast<float>(-150.0 - loading * 50.0);
    for (float height : elevation)
        EXPECT_FLOAT_EQ(height, expected);
}

TEST_F(Landscape, IsostasyZeroFractionsChangeNothing)
{
    const std::size_t size = 97;
    std::vector<float> erosion(size), deposition(size), elevation(size);
    for (std::size_t cell = 0; cell < size; cell++)
    {
        erosion[cell] = static_cast<float>(cell % 7) * 3.0f;
        deposition[cell] = static_cast<float>(cell % 5) * 2.0f;
        elevation[cell] = static_cast<float>(cell) * 11.5f - 300.0f;
    }
    const std::vector<float> before = elevation;
    landscape::applyisostasy(landscape::IsostasyParameters{}, erosion, deposition, elevation);
    EXPECT_EQ(elevation, before);
}

// Each cell responds to its own thicknesses only (local compensation, no flexure):
// a cell without erosion or deposition keeps its elevation next to cells that move,
// and a cell that both lost and gained moves by the difference of the two responses.
TEST_F(Landscape, IsostasyCellsAreIndependent)
{
    const std::vector<float> erosion = {0.0f, 40.0f, 0.0f, 10.0f, 0.0f};
    const std::vector<float> deposition = {0.0f, 0.0f, 0.0f, 20.0f, 30.0f};
    const std::vector<float> start = {500.0f, 460.0f, 300.0f, 120.0f, -20.0f};
    std::vector<float> elevation = start;
    landscape::applyisostasy(defaults(), erosion, deposition, elevation);
    for (std::size_t cell = 0; cell < start.size(); cell++)
    {
        const float expected = static_cast<float>(start[cell] + rebound * erosion[cell] - loading * deposition[cell]);
        EXPECT_FLOAT_EQ(elevation[cell], expected) << "cell " << cell;
    }
    EXPECT_EQ(elevation[0], start[0]);
    EXPECT_EQ(elevation[2], start[2]);
}

// With equal fractions (the default: the kernel's volumes are rock volumes) a step
// that moves rock without losing any leaves the area-weighted response zero: the
// column is neither created nor destroyed. On a sphere grid (cell areas differ by
// row) each row's western half loses what its eastern half gains, in another order.
TEST_F(Landscape, IsostasyEqualFractionsConserveColumn)
{
    const landscape::SphereMetric metric(64, 32, radius);
    std::vector<float> erosion(metric.size(), 0.0f), deposition(metric.size(), 0.0f);
    double moved = 0.0;
    for (int y = 0; y < metric.height; y++)
    {
        for (int x = 0; x < metric.width / 2; x++)
        {
            const float thickness = static_cast<float>(1 + (x * 7 + y * 3) % 23) * 1.25f;
            erosion[metric.index(x, y)] = thickness;
            deposition[metric.index(metric.width - 1 - x, y)] = thickness;
            moved += thickness * metric.area_m2(metric.index(x, y));
        }
    }
    landscape::IsostasyParameters parameters;
    parameters.reboundfraction = rebound;
    parameters.loadingfraction = rebound;
    std::vector<float> elevation(metric.size(), 0.0f);
    landscape::applyisostasy(parameters, erosion, deposition, elevation);
    double response = 0.0, lifted = 0.0;
    for (std::size_t cell = 0; cell < metric.size(); cell++)
    {
        response += static_cast<double>(elevation[cell]) * metric.area_m2(cell);
        lifted += std::max(0.0, static_cast<double>(elevation[cell])) * metric.area_m2(cell);
    }
    EXPECT_NEAR(lifted / moved, rebound, 1.0e-6);
    EXPECT_LE(std::abs(response) / moved, 1.0e-9);

    // Unequal fractions do not conserve it: the response is (rebound - loading) x moved.
    elevation.assign(metric.size(), 0.0f);
    landscape::applyisostasy(defaults(), erosion, deposition, elevation);
    response = 0.0;
    for (std::size_t cell = 0; cell < metric.size(); cell++)
        response += static_cast<double>(elevation[cell]) * metric.area_m2(cell);
    EXPECT_NEAR(response / moved, rebound - loading, 1.0e-6);
}

TEST_F(Landscape, IsostasyRejectsSizeMismatch)
{
    const std::vector<float> erosion(4, 1.0f), deposition(3, 0.0f);
    std::vector<float> elevation(4, 0.0f);
    EXPECT_THROW(landscape::applyisostasy(defaults(), erosion, deposition, elevation), std::invalid_argument);
}
