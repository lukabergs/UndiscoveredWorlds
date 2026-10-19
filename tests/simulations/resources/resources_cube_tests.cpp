// The resource layers on the cube (docs/architecture/sphere-mesh.md section 12, W2-R):
// mineral and marine reserves, gold deposits and placers on a cube planet of N 16
// with synthetic land, sea and deposits.

#include <array>
#include <cmath>
#include <cstdint>
#include <memory>
#include <vector>

#include <gtest/gtest.h>

#include "core/climate_seasons.hpp"
#include "core/deterministic_random.hpp"
#include "core/mesh/cell_grid.hpp"
#include "simulations/geology/mineral_systems/mineral_deposit.hpp"
#include "simulations/resources/gold_deposits.hpp"
#include "simulations/resources/marine_resources.hpp"
#include "simulations/resources/mineral_resources.hpp"
#include "simulations/resources/physical_layers_internal.hpp"
#include "world/planet.hpp"

namespace
{
constexpr int facecells = 16;
constexpr double kPi = 3.14159265358979323846;
constexpr int depositcount = 600;

// A deterministic uniform double in [0, 1) for draw k of a seed.
double uniform(std::uint64_t seed, int k)
{
    return static_cast<double>(deterministiccombine(seed, k) >> 11) * (1.0 / 9007199254740992.0);
}

int faceof(mesh::CellId c)
{
    return static_cast<int>(c / static_cast<mesh::CellId>(facecells * facecells));
}

// Land rising toward +X (lon 0), sea where the centre's x is below -0.4 and in nine
// single-cell seas per face; a mild climate; erosion potential on land; a volcano and
// a convergent boundary on one land cell of each face (the volcanic and ore reserves read
// discs of one or two legacy pixels, which on these 625 km cells hold only the cell itself);
// primary deposits at points uniform
// on the sphere (u by longitude, v by latitude as the world fractions), each with
// liberated gold for the placers. On lat-lon (64 x 32) the same rules, i and j read from
// the id as if it were a cube cell.
std::unique_ptr<planet> resourcecube(mesh::GridKind kind = mesh::GridKind::cube)
{
    auto world = std::make_unique<planet>();
    world->setseed(1);
    world->setgrid(kind, 4 * facecells);
    world->setsealevel(0);
    const mesh::CellGrid& grid = world->grid();

    for (mesh::CellId c = 0; c < world->cellcount(); c++)
    {
        const mesh::Vec3 p = grid.centre(c);
        const int i = static_cast<int>(c % facecells);
        const int j = static_cast<int>((c / facecells) % facecells);
        const bool sea = p.x < -0.4 || (i % 5 == 2 && j % 5 == 2);
        world->setnomfloat(c, sea ? -2000.0f : static_cast<float>(50.0 + 2000.0 * (p.x + 0.4)));
        world->setjantempfloat(c, 12.0f);
        world->setjultempfloat(c, 16.0f);
        world->setjanrainfloat(c, 80.0f);
        world->setjulrainfloat(c, 80.0f);
        world->setriverdir(c, 0);
        world->seterosionpotential(c, sea ? 0 : 50);
        world->setdepositionpotential(c, sea ? 0 : 30);
        for (int season = 0; season < CLIMATESEASONCOUNT; season++)
            world->setseasonalsst(season, c, static_cast<int>(std::lround(20.0 * std::cos(grid.latitude_rad(c)))));
        if (i == 8 && j == 8)
        {
            world->setvolcano(c, 1500.0f);
            world->settectonicboundarytype(c, BoundaryType::convergent);
            world->setgeologicregime(c, GeologicRegime::convergent_arc);
        }
    }

    const std::uint64_t seed = deterministiccontextseed(1, 0x5243);
    std::vector<mineralsystems::Deposit>& deposits = world->mineraldeposits();
    for (int k = 0; k < depositcount; k++)
    {
        mineralsystems::Deposit deposit;
        deposit.id = static_cast<std::uint32_t>(k + 1);
        deposit.system = mineralsystems::DepositSystem::orogenic_gold;
        deposit.u = uniform(seed, 2 * k);
        deposit.v = std::acos(1.0 - 2.0 * uniform(seed, 2 * k + 1)) / kPi;  // colatitude uniform by area
        deposit.oreMt = 2.0;
        deposit.goldGradeGt = 5.0;
        deposit.recentLiberatedGoldT = 4.0;
        deposits.push_back(deposit);
    }
    return world;
}

void generatelayers(planet& world)
{
    physical_layers::generate_mineral_resources(world);
    physical_layers::generate_gold_deposits(world);
    physical_layers::generate_marine_resources(world, world.continentalshelves());
}
}

TEST(ResourcesCube, LayersRunWithFiniteOutputs)
{
    const std::unique_ptr<planet> world = resourcecube();
    generatelayers(*world);

    int volcanic = 0, ore = 0, placer = 0, fishery = 0, coastalfishery = 0;
    for (mesh::CellId c = 0; c < world->cellcount(); c++)
    {
        for (const int score : { world->volcanicreserve(c), world->metalorereserve(c), world->placerreserve(c),
                 world->evaporitereserve(c), world->fisheryreserve(c) })
        {
            EXPECT_GE(score, 0);
            EXPECT_LE(score, 100);
        }
        volcanic += world->volcanicreserve(c) > 0 ? 1 : 0;
        ore += world->metalorereserve(c) > 0 ? 1 : 0;
        placer += world->placerreserve(c) > 0 ? 1 : 0;
        fishery += world->sea(c) && world->fisheryreserve(c) > 0 ? 1 : 0;
        coastalfishery += !world->sea(c) && world->fisheryreserve(c) > 0 ? 1 : 0;
    }
    EXPECT_GT(volcanic, 0);
    EXPECT_GT(ore, 0);
    EXPECT_GT(placer, 0);
    EXPECT_GT(fishery, 0);
    EXPECT_GT(coastalfishery, 0);

    int placers = 0;
    for (const mineralsystems::Deposit& d : world->mineraldeposits())
    {
        ASSERT_TRUE(d.surface.evaluated);
        ASSERT_LT(d.surface.cell, world->cellcount());
        EXPECT_GE(d.surface.worldX, 0);
        EXPECT_LT(d.surface.worldX, 4 * facecells);
        EXPECT_GE(d.surface.worldY, 0);
        EXPECT_LT(d.surface.worldY, 2 * facecells);
        EXPECT_TRUE(std::isfinite(d.surface.elevationM));
        EXPECT_TRUE(std::isfinite(d.surface.localReliefM));
        EXPECT_TRUE(std::isfinite(d.containedGoldT()));
        if (d.system == mineralsystems::DepositSystem::placer_gold)
        {
            ++placers;
            EXPECT_FALSE(world->sea(d.surface.cell));
        }
    }
    EXPECT_GT(placers, 0);

    // The liberated gold is conserved over rivers, lakes, the sea and the placers.
    const mineralsystems::MineralBudget& budget = world->mineralbudget();
    EXPECT_NEAR(budget.placerInjectedGoldT, 4.0 * depositcount, 1e-6);
    EXPECT_NEAR(budget.placerDepositGoldT + budget.placerDispersedGoldT + budget.placerLakeGoldT + budget.placerExportedGoldT,
        budget.placerInjectedGoldT, 1e-6);
}

TEST(ResourcesCube, DepositsSpreadByArea)
{
    // Deposits at points uniform on the sphere land evenly on the six faces.
    const std::unique_ptr<planet> world = resourcecube();
    physical_layers::generate_gold_deposits(*world);

    std::array<int, 6> perface{};
    for (const mineralsystems::Deposit& d : world->mineraldeposits())
    {
        if (d.system == mineralsystems::DepositSystem::placer_gold)
            continue;
        ASSERT_LT(d.surface.cell, world->cellcount());
        perface[faceof(d.surface.cell)]++;
        // The cell is the one containing the deposit's point.
        const double latitude = kPi * (0.5 - d.v);
        const double longitude = kPi * (2.0 * d.u - 1.0);
        EXPECT_EQ(d.surface.cell, world->grid().locate(mesh::unitvector(latitude, longitude)));
    }
    for (int f = 0; f < 6; f++)
    {
        EXPECT_GE(perface[f], 60) << "face " << f;
        EXPECT_LE(perface[f], 140) << "face " << f;
    }
}

TEST(ResourcesCube, DrainageCrossesFaceEdges)
{
    // The receivers follow the neighbour table across face edges.
    const std::unique_ptr<planet> world = resourcecube();
    const std::vector<std::int32_t> receiver = physical_layers::cellqueries::drainagereceivers(*world);
    ASSERT_EQ(receiver.size(), world->cellcount());

    int crossing = 0;
    for (mesh::CellId c = 0; c < world->cellcount(); c++)
    {
        if (receiver[c] < 0)
            continue;
        EXPECT_LT(world->nomfloat(static_cast<mesh::CellId>(receiver[c])), world->nomfloat(c));
        crossing += faceof(static_cast<mesh::CellId>(receiver[c])) != faceof(c) ? 1 : 0;
    }
    EXPECT_GT(crossing, 0);
}

TEST(ResourcesCube, PlacerCellIsTheTrappingCell)
{
    // A placer lies in the cell where the routing trapped its gold: the cell of its point is
    // its surface cell, and that cell is downstream (drainagereceivers) of every primary
    // deposit it names as a source. On the cube and on lat-lon (a uniform point of the cell).
    for (const mesh::GridKind kind : { mesh::GridKind::cube, mesh::GridKind::latlon })
    {
        const std::unique_ptr<planet> world = resourcecube(kind);
        physical_layers::generate_gold_deposits(*world);
        const mesh::CellGrid& grid = world->grid();
        const std::vector<std::int32_t> receiver = physical_layers::cellqueries::drainagereceivers(*world);

        std::vector<mesh::CellId> primarycell(depositcount + 1, mesh::nocell);
        for (const mineralsystems::Deposit& d : world->mineraldeposits())
            if (d.system != mineralsystems::DepositSystem::placer_gold && d.id <= static_cast<std::uint32_t>(depositcount))
                primarycell[d.id] = d.surface.cell;

        int placers = 0;
        int checkedsources = 0;
        for (const mineralsystems::Deposit& d : world->mineraldeposits())
        {
            if (d.system != mineralsystems::DepositSystem::placer_gold)
                continue;
            ++placers;
            ASSERT_LT(d.surface.cell, world->cellcount());
            const double latitude = kPi * (0.5 - d.v);
            const double longitude = kPi * (2.0 * d.u - 1.0);
            EXPECT_EQ(d.surface.cell, grid.locate(mesh::unitvector(latitude, longitude))) << "placer " << d.id;
            if (kind == mesh::GridKind::latlon)
            {
                EXPECT_EQ(d.surface.worldX, grid.column(d.surface.cell));
                EXPECT_EQ(d.surface.worldY, grid.row(d.surface.cell));
            }

            EXPECT_FALSE(d.sources.empty()) << "placer " << d.id;
            for (const mineralsystems::SourceShare& share : d.sources)
            {
                ASSERT_LE(share.depositId, static_cast<std::uint32_t>(depositcount));
                std::int32_t c = static_cast<std::int32_t>(primarycell[share.depositId]);
                ASSERT_GE(c, 0);
                bool reached = false;
                for (std::size_t steps = 0; c >= 0 && steps <= world->cellcount(); ++steps)
                {
                    if (static_cast<mesh::CellId>(c) == d.surface.cell)
                    {
                        reached = true;
                        break;
                    }
                    c = receiver[static_cast<std::size_t>(c)];
                }
                EXPECT_TRUE(reached) << "placer " << d.id << " source " << share.depositId;
                ++checkedsources;
            }
        }
        EXPECT_GT(placers, 0);
        EXPECT_GT(checkedsources, 0);
    }
}
