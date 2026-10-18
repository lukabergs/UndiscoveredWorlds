// Society on the cube (docs/architecture/sphere-mesh.md section 12, W2-S):
// sites, routes, polities and trade on a cube planet of N 16 with synthetic land.

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <numeric>
#include <set>
#include <vector>

#include <gtest/gtest.h>

#include "core/climate_seasons.hpp"
#include "core/mesh/cell_grid.hpp"
#include "simulations/society/social_generation.hpp"
#include "simulations/society/social_generation_internal.hpp"
#include "world/planet.hpp"

namespace
{
constexpr int facecells = 16;

int faceof(mesh::CellId c)
{
    return static_cast<int>(c / static_cast<mesh::CellId>(facecells * facecells));
}

// A cube planet with land where keep(c) holds, a mild climate and a moist forest everywhere.
template <typename Keep>
std::unique_ptr<planet> cubeworld(Keep&& keep)
{
    auto world = std::make_unique<planet>();
    world->setseed(1);
    world->setgrid(mesh::GridKind::cube, 4 * facecells);
    world->setsealevel(0);

    for (mesh::CellId c = 0; c < world->cellcount(); c++)
    {
        world->setnomfloat(c, keep(c) ? 300.0f : -2000.0f);
        world->setjantempfloat(c, 12.0f);
        world->setjultempfloat(c, 16.0f);
        world->setjanrainfloat(c, 80.0f);
        world->setjulrainfloat(c, 80.0f);
        world->setbiome(c, biometropicalmoistforest);
        world->setriverdir(c, 0);
    }

    return world;
}

// Uniform land on every face, with the same small seas (nine single cells per
// face, away from the face edges) so that coasts and harbours rank the sites.
std::unique_ptr<planet> socialcube()
{
    std::unique_ptr<planet> world = cubeworld([](mesh::CellId c)
    {
        const int i = static_cast<int>(c % facecells);
        const int j = static_cast<int>((c / facecells) % facecells);
        return !(i % 5 == 2 && j % 5 == 2);
    });
    SocialGenerationOptions options;
    options.enabled = true;
    generatesocialworld(*world, options);
    return world;
}
}

TEST(SocietyCube, RunsWithFiniteOutputs)
{
    const std::unique_ptr<planet> world = socialcube();
    const std::vector<Settlement>& settlements = world->settlements();
    ASSERT_GE(settlements.size(), 8u);

    for (const Settlement& s : settlements)
    {
        ASSERT_LT(s.cell, world->cellcount());
        EXPECT_FALSE(world->sea(s.cell));
        EXPECT_GE(s.x, 0);
        EXPECT_LT(s.x, 4 * facecells);
        EXPECT_GE(s.y, 0);
        EXPECT_LT(s.y, 2 * facecells);
        EXPECT_GT(s.urbanPopulation + s.ruralPopulation, 0);
        EXPECT_TRUE(std::isfinite(s.carryingCapacity));
        EXPECT_TRUE(std::isfinite(s.infrastructure));
        EXPECT_TRUE(std::isfinite(s.marketStrength));
        EXPECT_GE(s.polityId, 0);
    }

    ASSERT_FALSE(world->routeedges().empty());
    for (const RouteEdge& r : world->routeedges())
    {
        EXPECT_TRUE(std::isfinite(r.cost) && r.cost >= 1.0f);
        EXPECT_TRUE(std::isfinite(r.capacity) && r.capacity > 0.0f);
        EXPECT_TRUE(std::isfinite(r.traffic) && r.traffic > 0.0f);
    }

    ASSERT_FALSE(world->polities().empty());
    for (const Polity& p : world->polities())
    {
        EXPECT_TRUE(std::isfinite(p.cohesion));
        EXPECT_TRUE(std::isfinite(p.wealth));
        EXPECT_GT(p.population, 0);
    }

    // Every land cell is reached by a settlement's influence (the table flood crosses face edges).
    int owned = 0;
    for (mesh::CellId c = 0; c < world->cellcount(); c++)
        owned += world->ownersettlementid(c) >= 0 ? 1 : 0;
    EXPECT_GT(owned, static_cast<int>(world->cellcount()) / 2);
}

TEST(SocietyCube, RouteCrossesFaceEdge)
{
    const std::unique_ptr<planet> world = socialcube();
    const mesh::CellGrid& grid = world->grid();
    const std::vector<Settlement>& settlements = world->settlements();

    bool crossfaceroute = false;
    for (const RouteEdge& r : world->routeedges())
        crossfaceroute = crossfaceroute || faceof(settlements[r.fromSettlementId].cell) != faceof(settlements[r.toSettlementId].cell);
    EXPECT_TRUE(crossfaceroute);

    // The traffic of some route continues across a face edge.
    bool crossing = false;
    for (mesh::CellId c = 0; c < grid.cellcount() && !crossing; c++)
    {
        if (world->routetraffic(c) <= 0)
            continue;
        for (int slot = 1; slot <= 8; slot++)
        {
            const int n = grid.neighbour(c, slot);
            if (n >= 0 && faceof(static_cast<mesh::CellId>(n)) != faceof(c) && world->routetraffic(static_cast<mesh::CellId>(n)) > 0)
                crossing = true;
        }
    }
    EXPECT_TRUE(crossing);
}

TEST(SocietyCube, SitesDrawnByArea)
{
    // Equal-scored candidates are taken in an order drawn by area: the first
    // 600 of the 1536 cells of uniform land fall evenly on the six faces.
    const std::shared_ptr<const mesh::CellGrid> grid = mesh::makecubegrid(facecells, 6371000.0);
    const std::vector<double> keys = socialgen::sitedrawkeys(*grid, 1);
    ASSERT_EQ(keys.size(), grid->cellcount());

    std::vector<mesh::CellId> order(grid->cellcount());
    std::iota(order.begin(), order.end(), 0u);
    std::sort(order.begin(), order.end(), [&](mesh::CellId a, mesh::CellId b) { return keys[a] > keys[b]; });

    const int taken = 600;
    std::array<int, 6> perface{};
    for (int k = 0; k < taken; k++)
        perface[faceof(order[k])]++;
    for (int f = 0; f < 6; f++)
    {
        EXPECT_GE(perface[f], 70) << "face " << f;
        EXPECT_LE(perface[f], 130) << "face " << f;
    }

    // On uniform land the settlements spread over the sphere instead of filling the first face.
    const std::unique_ptr<planet> world = socialcube();
    std::set<int> faces;
    for (const Settlement& s : world->settlements())
        faces.insert(faceof(s.cell));
    EXPECT_GE(faces.size(), 4u);
}
