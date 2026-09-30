#include "gtest/gtest.h"

#include "placer_routing.hpp"

#include <cmath>

using namespace mineralsystems;

namespace
{
DrainageNetwork emptynetwork(int width, int height, float stepkm, float transportkm)
{
    DrainageNetwork network;
    network.width = width;
    network.height = height;
    const std::size_t count = network.cellCount();
    network.receiver.assign(count, -1);
    network.sea.assign(count, 0);
    network.sink.assign(count, 0);
    network.stepLengthKm.assign(count, stepkm);
    network.transportLengthKm.assign(count, transportkm);
    return network;
}

// Each row drains eastward into a sea cell in the last column.
DrainageNetwork eastwardrivers(int width, int height, float stepkm, float transportkm)
{
    DrainageNetwork network = emptynetwork(width, height, stepkm, transportkm);
    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width - 1; ++x)
            network.receiver[static_cast<std::size_t>(y * width + x)] = y * width + x + 1;
        network.sea[static_cast<std::size_t>(y * width + width - 1)] = 1;
    }
    return network;
}
}

TEST(PlacerRouting, PlacersRequireAConnectedSource)
{
    const DrainageNetwork network = eastwardrivers(10, 2, 20.0f, 30.0f);
    const PlacerRouting routing = routeplacergold(network, { { 7, 0, 10.0 } });

    for (int x = 0; x < 10; ++x)
        EXPECT_EQ(routing.riverGoldT[static_cast<std::size_t>(10 + x)], 0.0) << "second catchment, x " << x;
    EXPECT_GT(routing.riverGoldT[0], routing.riverGoldT[1]);
    EXPECT_GT(routing.riverGoldT[1], 0.0);

    MineralBudget budget;
    const auto placers = makeplacerdeposits(routing, network, PlacerParameters{}, 600.0, 100, 1, budget);
    ASSERT_FALSE(placers.empty());
    for (const Deposit& placer : placers)
    {
        EXPECT_LT(placer.v, 0.5); // First row only.
        ASSERT_EQ(placer.sources.size(), 1u);
        EXPECT_EQ(placer.sources[0].depositId, 7u);
    }
    EXPECT_TRUE(routeplacergold(network, {}).contributions.empty());
}

TEST(PlacerRouting, RoutingConservesGoldAcrossSinksCyclesAndOutlets)
{
    DrainageNetwork network = emptynetwork(6, 1, 50.0f, 40.0f);
    network.receiver = { 1, 2, -1, 4, 3, 4 }; // 0->1->2 dead end; 3<->4 cycle; 5 feeds the cycle.
    network.sink[2] = 0;
    DrainageNetwork lake = network;
    lake.sink[1] = 1;
    DrainageNetwork sea = network;
    sea.sea[1] = 1;

    const std::vector<PlacerSource> sources{ { 1, 0, 3.0 }, { 2, 5, 4.5 }, { 3, 3, 1.25 } };
    for (const DrainageNetwork* candidate : { &network, &lake, &sea })
    {
        const PlacerRouting routing = routeplacergold(*candidate, sources);
        EXPECT_NEAR(routing.injectedT, 8.75, 1e-12);
        EXPECT_NEAR(routing.riverT + routing.lakeT + routing.exportedT, routing.injectedT, 1e-12);
        double cells = 0.0;
        for (double gold : routing.riverGoldT)
            cells += gold;
        EXPECT_NEAR(cells, routing.riverT, 1e-12);

        MineralBudget budget;
        const auto placers = makeplacerdeposits(routing, *candidate, PlacerParameters{}, 0.0, 1, 2, budget);
        double deposited = 0.0;
        for (const Deposit& placer : placers)
            deposited += placer.containedGoldT();
        EXPECT_NEAR(deposited, budget.placerDepositGoldT, 1e-9);
        EXPECT_NEAR(budget.placerDepositGoldT + budget.placerDispersedGoldT + budget.placerLakeGoldT
            + budget.placerExportedGoldT, budget.placerInjectedGoldT, 1e-12);
    }
    EXPECT_NEAR(routeplacergold(lake, sources).lakeT, 3.0 * std::exp(-50.0 / 40.0), 1e-12);
    EXPECT_GT(routeplacergold(sea, sources).exportedT, 0.0);
}

TEST(PlacerRouting, TrappingDependsOnDistanceNotCellCount)
{
    // A 600 km channel resolved with 30 km and 15 km steps.
    const double transportkm = 50.0;
    double within150[2] = {};
    double exported[2] = {};
    const int cells[2] = { 20, 40 };
    for (int grid = 0; grid < 2; ++grid)
    {
        const float step = static_cast<float>(600.0 / cells[grid]);
        DrainageNetwork network = eastwardrivers(cells[grid] + 1, 1, step, static_cast<float>(transportkm));
        const PlacerRouting routing = routeplacergold(network, { { 1, 0, 1.0 } });
        for (int x = 0; x < cells[grid] / 4; ++x)
            within150[grid] += routing.riverGoldT[static_cast<std::size_t>(x)];
        exported[grid] = routing.exportedT;
    }
    EXPECT_NEAR(within150[0], 1.0 - std::exp(-150.0 / transportkm), 1e-9);
    EXPECT_NEAR(within150[1], within150[0], 1e-9);
    EXPECT_NEAR(exported[1], exported[0], 1e-12);
    EXPECT_NEAR(exported[0], std::exp(-600.0 / transportkm), 1e-12);
}

TEST(PlacerRouting, PlacerDepositsNameTheirLargestSources)
{
    const DrainageNetwork network = eastwardrivers(4, 1, 100.0f, 10.0f);
    const PlacerRouting routing = routeplacergold(network, { { 21, 1, 1.0 }, { 22, 1, 3.0 } });
    MineralBudget budget;
    const auto placers = makeplacerdeposits(routing, network, PlacerParameters{}, 500.0, 40, 9, budget);

    ASSERT_EQ(placers.size(), 1u);
    const Deposit& placer = placers[0];
    EXPECT_EQ(placer.id, 40u);
    EXPECT_EQ(placer.system, DepositSystem::placer_gold);
    EXPECT_EQ(placer.origin, DepositOrigin::redistributed);
    ASSERT_EQ(placer.sources.size(), 2u);
    EXPECT_EQ(placer.sources[0].depositId, 22u);
    EXPECT_NEAR(placer.sources[0].goldTonnes + placer.sources[1].goldTonnes, placer.containedGoldT(), 1e-9);
    EXPECT_NEAR(placer.oreMt * placer.goldGradeGt, placer.containedGoldT(), 1e-12);

    Deposit source;
    source.id = 22;
    source.system = DepositSystem::epithermal_gold;
    const std::vector<Deposit> all{ source, placer };
    const std::string text = describedeposit(placer, &all);
    EXPECT_NE(text.find("from #22 (Epithermal gold)"), std::string::npos);
    EXPECT_NE(text.find("from #21"), std::string::npos);
}
