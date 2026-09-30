#include "gtest/gtest.h"

#include "frame_fixture.hpp"
#include "mineral_history.hpp"

#include <cmath>
#include <map>

using namespace mineralsystems;

namespace
{
constexpr double kCircumferenceKm = 40030.0;

// Two continental plates meeting at a stationary collision belt.
FrameFixture beltworld(std::uint32_t width)
{
    FrameFixture fixture(width, width / 2);
    fixture.frame.worldWidthKm = kCircumferenceKm;
    BoundarySummary belt;
    belt.boundaryType = boundarycode::convergent;
    belt.geologicRegime = regimecode::continent_collision;
    belt.convergenceScore = 80;
    belt.ageMyr = 50.0f;
    belt.leftPlate = 0;
    belt.rightPlate = 1;
    fixture.frame.boundaries = { BoundarySummary{}, belt };
    fixture.frame.plates = { { 1 }, { 2 } };
    for (std::uint32_t y = 0; y < fixture.frame.height; ++y)
    {
        for (std::uint32_t x = 0; x < width; ++x)
        {
            const std::size_t cell = fixture.cell(x, y);
            fixture.plate[cell] = x < width / 2 ? 0 : 1;
            fixture.nearest[cell] = 1;
            fixture.distance[cell] = static_cast<std::uint16_t>(std::abs(static_cast<int>(x) - static_cast<int>(width / 2)));
        }
    }
    return fixture;
}

double sumgold(const std::vector<Deposit>& deposits, bool remaining)
{
    double total = 0.0;
    for (const Deposit& deposit : deposits)
        total += remaining ? (deposit.fate == DepositFate::preserved ? deposit.remainingGoldT() : 0.0) : deposit.containedGoldT();
    return total;
}
}

TEST(MineralHistory, DepositsPersistAndTravelAfterTheirBoundaryDies)
{
    FrameFixture fixture = beltworld(128);
    MineralParameters parameters;
    parameters.inheritedProvinces = false;
    MineralHistory history(parameters, 5);
    history.begin(fixture.view());

    for (int update = 1; update <= 20; ++update)
    {
        fixture.frame.timeMyr = update;
        history.step(fixture.view());
    }
    const std::vector<Deposit> formed = history.deposits();
    const std::vector<Marker> formedmarkers = history.markers();
    ASSERT_GT(formed.size(), 10u);

    // The boundary stops and the plate carries its crust eastward.
    fixture.frame.boundaries[1] = BoundarySummary{};
    fixture.motion(0.5f);
    for (int update = 21; update <= 70; ++update)
    {
        fixture.frame.timeMyr = update;
        history.step(fixture.view());
    }

    ASSERT_EQ(history.deposits().size(), formed.size()); // Nothing forms without a boundary.
    for (std::size_t i = 0; i < formed.size(); ++i)
    {
        const Deposit& deposit = history.deposits()[i];
        EXPECT_EQ(deposit.fate, DepositFate::preserved);
        EXPECT_EQ(deposit.formation.geologicRegime, regimecode::continent_collision);
        EXPECT_LE(deposit.formation.timeMyr, 20.0);
        const double moved = std::fmod(history.markers()[i].x - formedmarkers[i].x + 128.0, 128.0);
        EXPECT_NEAR(moved, 25.0, 1e-9);
    }
}

TEST(MineralHistory, ExhumationErodesTheOreBodyUniformly)
{
    Deposit deposit;
    deposit.oreMt = 2.0;
    deposit.goldGradeGt = 5.0; // 10 t.
    deposit.formationDepthM = 3000.0;
    deposit.verticalExtentM = 1000.0;

    EXPECT_DOUBLE_EQ(exhumedeposit(deposit, 2500.0), 0.0);
    EXPECT_EQ(exposureclass(deposit.topDepthM()), Exposure::buried);
    EXPECT_NEAR(exhumedeposit(deposit, 750.0), 2.5, 1e-12);
    EXPECT_NEAR(deposit.remainingGoldT(), 7.5, 1e-12);
    EXPECT_EQ(exposureclass(deposit.topDepthM()), Exposure::outcropping);
    EXPECT_NEAR(exhumedeposit(deposit, 5000.0), 7.5, 1e-12);
    EXPECT_DOUBLE_EQ(deposit.erodedFraction, 1.0);
    EXPECT_DOUBLE_EQ(exhumedeposit(deposit, 100.0), 0.0);
}

TEST(MineralHistory, DenudationRisesSteeplyWithElevation)
{
    const ErosionParameters erosion;
    EXPECT_DOUBLE_EQ(denudationrate(-200.0, erosion), 0.0);
    EXPECT_NEAR(denudationrate(300.0, erosion), 0.135, 1e-12); // Craton.
    EXPECT_NEAR(denudationrate(2000.0, erosion), 40.0, 1e-9);
    EXPECT_NEAR(denudationrate(4000.0, erosion), 320.0, 1e-9); // High orogen.
    EXPECT_DOUBLE_EQ(denudationrate(20000.0, erosion), erosion.maximumRateMMyr);
}

TEST(MineralHistory, DenudationUsesCurrentEngineElevation)
{
    FrameFixture fixture = beltworld(64);
    MineralParameters parameters;
    parameters.formation.provincesPerMkm2 = 5.0;
    MineralHistory history(parameters, 1);
    history.begin(fixture.view());
    ASSERT_FALSE(history.deposits().empty());
    std::vector<double> initial;
    for (const auto& deposit : history.deposits()) initial.push_back(deposit.exhumationM);
    for (float& elevation : fixture.elevation) elevation = 300.0f;
    fixture.frame.timeMyr = 1;
    history.step(fixture.view());
    for (std::size_t i = 0; i < initial.size(); ++i)
        if (history.deposits()[i].fate == DepositFate::preserved)
            EXPECT_NEAR(history.deposits()[i].exhumationM - initial[i], 0.135, 1e-8);
    std::vector<double> after;
    for (const auto& deposit : history.deposits()) after.push_back(deposit.exhumationM);
    for (float& elevation : fixture.elevation) elevation = -100.0f;
    fixture.frame.timeMyr = 2;
    history.step(fixture.view());
    for (std::size_t i = 0; i < after.size(); ++i)
        EXPECT_DOUBLE_EQ(history.deposits()[i].exhumationM, after[i]);
}

TEST(MineralHistory, GoldBudgetBalancesThroughFormationAndErosion)
{
    FrameFixture fixture = beltworld(128);
    // High ground over the belt exhumes deposits within the run.
    for (std::uint32_t y = 0; y < fixture.frame.height; ++y)
        for (std::uint32_t x = 56; x < 72; ++x)
            fixture.elevation[fixture.cell(x, y)] = 3000.0f;
    MineralParameters parameters;
    parameters.formation.provincesPerMkm2 = 1.0;
    MineralHistory history(parameters, 9);
    history.begin(fixture.view());
    fixture.motion(0.2f);
    for (int update = 1; update <= 150; ++update)
    {
        fixture.frame.timeMyr = update;
        history.step(fixture.view());
    }

    const MineralBudget& budget = history.budget();
    double formed = budget.inheritedGoldT;
    for (double gold : budget.formedGoldT)
        formed += gold;
    const double contained = sumgold(history.deposits(), false);
    ASSERT_GT(contained, 0.0);
    ASSERT_GT(budget.erodedGoldT, 0.0);
    ASSERT_GT(budget.erodedCount, 0);
    EXPECT_NEAR(formed, contained, 1e-9 * contained);
    EXPECT_NEAR(sumgold(history.deposits(), true) + budget.erodedGoldT + budget.subductedGoldT, contained, 1e-9 * contained);

    std::map<DepositFate, int> fates;
    for (const Deposit& deposit : history.deposits())
    {
        ++fates[deposit.fate];
        if (deposit.fate == DepositFate::eroded)
            EXPECT_DOUBLE_EQ(deposit.remainingGoldT(), 0.0);
        EXPECT_GE(deposit.recentLiberatedGoldT, 0.0);
    }
    EXPECT_EQ(fates[DepositFate::eroded], budget.erodedCount);
}

TEST(MineralHistory, RepeatedFramesNeitherErodeNorFormDeposits)
{
    FrameFixture fixture = beltworld(64);
    for (float& elevation : fixture.elevation)
        elevation = 2000.0f;
    MineralParameters parameters;
    parameters.formation.provincesPerMkm2 = 5.0;
    MineralHistory history(parameters, 3);
    fixture.frame.timeMyr = 10.0;
    history.begin(fixture.view());
    const std::size_t count = history.deposits().size();
    const double exhumation = count > 0 ? history.deposits()[0].exhumationM : 0.0;

    history.step(fixture.view());

    EXPECT_EQ(history.laststep().dtMyr, 0.0);
    EXPECT_EQ(history.laststep().markers.moved, 0u);
    EXPECT_EQ(history.deposits().size(), count);
    if (count > 0)
        EXPECT_EQ(history.deposits()[0].exhumationM, exhumation);
}

TEST(MineralHistory, DescriptionExplainsProvenanceAndQuantity)
{
    Deposit deposit;
    deposit.id = 12;
    deposit.system = DepositSystem::porphyry_copper_gold;
    deposit.oreMt = 200.0;
    deposit.goldGradeGt = 0.3;
    deposit.copperGradePct = 0.5;
    deposit.formationDepthM = 2500.0;
    deposit.verticalExtentM = 1000.0;
    deposit.exhumationM = 2400.0;
    deposit.formation.timeMyr = 412.0;
    deposit.formation.geologicRegime = regimecode::convergent_arc;
    deposit.formation.hostCrustClass = crustcode::continental;
    deposit.formation.boundaryDistanceKm = 150.0f;
    deposit.formation.boundarySegmentId = 17;

    const std::string text = describedeposit(deposit);
    EXPECT_NE(text.find("Porphyry copper-gold #12"), std::string::npos);
    EXPECT_NE(text.find("60.0 t gold"), std::string::npos);
    EXPECT_NE(text.find("1000 kt copper"), std::string::npos);
    EXPECT_NE(text.find("subduction-arc boundary (segment 17)"), std::string::npos);
    EXPECT_NE(text.find("continental crust"), std::string::npos);
    EXPECT_NE(text.find("shallow: top at 100 m depth"), std::string::npos);
}
