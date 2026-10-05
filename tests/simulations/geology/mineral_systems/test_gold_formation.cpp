#include "gtest/gtest.h"

#include "frame_fixture.hpp"
#include "simulations/geology/mineral_systems/gold_formation.hpp"

#include <cmath>
#include <cstdlib>

using namespace mineralsystems;

namespace
{
constexpr double kCircumferenceKm = 40030.0;

BoundarySummary maturecollision()
{
    BoundarySummary boundary;
    boundary.boundaryType = boundarycode::convergent;
    boundary.geologicRegime = regimecode::continent_collision;
    boundary.convergenceScore = 60;
    boundary.ageMyr = 60.0f;
    return boundary;
}

// Plate 0 lies west of x = width / 4 and plate 1 east of it; boundary segment 1
// separates them within a band of rows, and distances are measured in native
// cells. Continental crust stands in a box east of the boundary, and across it
// onto plate 0 when the plates collide. Elsewhere lies deep ocean floor.
FrameFixture boxworld(std::uint32_t width, BoundarySummary boundary, bool collision = true)
{
    FrameFixture fixture(width, width / 2);
    fixture.frame.worldWidthKm = kCircumferenceKm;
    boundary.leftPlate = 0;
    boundary.rightPlate = 1;
    fixture.frame.boundaries = { BoundarySummary{}, boundary };
    const std::uint32_t west = width / 4;
    const std::uint32_t east = width / 2;
    const std::uint32_t start = collision ? west - width / 8 : west;
    for (std::uint32_t y = 0; y < fixture.frame.height; ++y)
    {
        const bool band = y >= fixture.frame.height * 3 / 8 && y < fixture.frame.height * 5 / 8;
        for (std::uint32_t x = 0; x < width; ++x)
        {
            const std::size_t cell = fixture.cell(x, y);
            const bool inside = band && x >= start && x < east;
            fixture.plate[cell] = x < west ? 0 : 1;
            fixture.crust[cell] = inside ? crustcode::continental : crustcode::oceanic;
            fixture.elevation[cell] = inside ? 500.0f : -5000.0f;
            fixture.continental[cell] = inside ? 35000.0f : 0.0f;
            fixture.nearest[cell] = band ? 1 : 0;
            fixture.distance[cell] = static_cast<std::uint16_t>(x >= west ? x - west : west - 1 - x);
        }
    }
    fixture.frame.plates = { { 1 }, { 2 } };
    return fixture;
}

double totalgold(const FormationSample& sample)
{
    double total = 0.0;
    for (const NewDeposit& item : sample.deposits)
        total += item.deposit.containedGoldT();
    return total;
}
}

TEST(GoldFormation, SettingKernelsUsePhysicalDistance)
{
    FrameFixture fixture = boxworld(256, maturecollision());
    FormationParameters parameters;
    const SettingField field(fixture.view(), parameters);
    const double cellKm = fixture.frame.cellKm();
    const std::uint32_t y = fixture.frame.height / 2;
    const std::uint32_t west = 64;

    // Mature with saturated activity: the kernel is evaluated at cell centres,
    // half a cell beyond the boundary cell's distance.
    const auto kernel = [&](std::uint32_t cells)
    {
        return std::exp(-std::pow((cells + 0.5) * cellKm / parameters.orogenicWidthKm, 2.0));
    };
    EXPECT_NEAR(field.at(fixture.cell(west, y)).orogenic, kernel(0), 1e-6);
    const auto widthcells = static_cast<std::uint32_t>(std::lround(parameters.orogenicWidthKm / cellKm));
    EXPECT_NEAR(field.at(fixture.cell(west + widthcells, y)).orogenic, kernel(widthcells), 1e-6);
    EXPECT_FLOAT_EQ(field.at(fixture.cell(west, y)).arc, 0.0f); // Collision zones host no arc.
    EXPECT_NEAR(field.at(fixture.cell(west - 1, y)).orogenic, kernel(0), 1e-6); // Both plates deform.
    EXPECT_FLOAT_EQ(field.at(fixture.cell(0, y)).orogenic, 0.0f); // Deep ocean floor.
}

TEST(GoldFormation, SubductionMarginsHostArcsAwayFromTheTrench)
{
    BoundarySummary margin;
    margin.boundaryType = boundarycode::convergent;
    margin.geologicRegime = regimecode::convergent_arc;
    margin.convergenceScore = 100;
    margin.ageMyr = 5.0f; // Immature: no orogenic gold yet.
    FrameFixture fixture = boxworld(512, margin, false);
    FormationParameters parameters;
    const SettingField field(fixture.view(), parameters);
    const std::uint32_t y = fixture.frame.height / 2;
    const std::uint32_t west = 128;
    const auto peak = static_cast<std::uint32_t>(std::lround(parameters.arcPeakKm / fixture.frame.cellKm() - 0.5));
    ASSERT_GT(peak, 0u);

    EXPECT_GT(field.at(fixture.cell(west + peak, y)).arc, 0.85f);
    EXPECT_GT(field.at(fixture.cell(west + peak, y)).arc, field.at(fixture.cell(west, y)).arc);
    EXPECT_FLOAT_EQ(field.at(fixture.cell(west + peak, y)).orogenic, 0.0f);
}

TEST(GoldFormation, ConvergentMarginsAreClassifiedFromCrustColumns)
{
    const std::uint32_t y = 64;
    const std::uint32_t west = 64;
    const FormationParameters parameters;

    // Labelled a collision, but only plate 1 brings continental crust: plate 0
    // subducts beneath an Andean-type margin. An emergent cell of thin crust
    // on plate 0 is host crust, yet hosts neither setting.
    FrameFixture andean = boxworld(256, maturecollision(), false);
    andean.elevation[andean.cell(west - 1, y)] = 0.0f;
    const SettingField margin(andean.view(), parameters);
    EXPECT_EQ(margin.regime(1), regimecode::convergent_arc);
    EXPECT_GT(margin.at(andean.cell(west + 1, y)).arc, 0.0f);
    EXPECT_GT(margin.at(andean.cell(west + 1, y)).orogenic, 0.0f); // Accretionary orogen.
    EXPECT_FLOAT_EQ(margin.at(andean.cell(west - 1, y)).arc, 0.0f);
    EXPECT_FLOAT_EQ(margin.at(andean.cell(west - 1, y)).orogenic, 0.0f);

    // Buoyant crust on both edges collides, whatever the engine regime.
    BoundarySummary labelledarc = maturecollision();
    labelledarc.geologicRegime = regimecode::convergent_arc;
    FrameFixture collision = boxworld(256, labelledarc, true);
    const SettingField belt(collision.view(), parameters);
    EXPECT_EQ(belt.regime(1), regimecode::continent_collision);
    EXPECT_GT(belt.at(collision.cell(west - 1, y)).orogenic, margin.at(andean.cell(west, y)).orogenic);
    EXPECT_FLOAT_EQ(belt.at(collision.cell(west, y)).arc, 0.0f);
}

TEST(GoldFormation, OlderOceanFloorSubductsBetweenOceanicPlates)
{
    FrameFixture fixture = boxworld(256, maturecollision(), false);
    const std::uint32_t y = 64;
    const std::uint32_t west = 64;
    // No continental crust: emergent arc crust stands on both edges.
    for (float& crust : fixture.continental)
        crust = 0.0f;
    fixture.elevation[fixture.cell(west - 1, y)] = 200.0f;
    const auto ages = [&](float westMyr, float eastMyr)
    {
        for (std::size_t cell = 0; cell < fixture.age.size(); ++cell)
            fixture.age[cell] = fixture.plate[cell] == 0 ? westMyr : eastMyr;
        return SettingField(fixture.view(), FormationParameters{});
    };

    const SettingField eastover = ages(120.0f, 30.0f);
    EXPECT_EQ(eastover.regime(1), regimecode::convergent_arc);
    EXPECT_GT(eastover.at(fixture.cell(west + 1, y)).arc, 0.0f);
    EXPECT_FLOAT_EQ(eastover.at(fixture.cell(west - 1, y)).arc, 0.0f);

    const SettingField westover = ages(30.0f, 120.0f);
    EXPECT_GT(westover.at(fixture.cell(west - 1, y)).arc, 0.0f);
    EXPECT_FLOAT_EQ(westover.at(fixture.cell(west + 1, y)).arc, 0.0f);
}

TEST(GoldFormation, RefiningTheGridPreservesExpectedEndowment)
{
    FormationParameters parameters;
    double activearea[2] = {};
    double gold[2] = {};
    int count[2] = {};
    const std::uint32_t widths[2] = { 128, 256 };
    for (int grid = 0; grid < 2; ++grid)
    {
        FrameFixture fixture = boxworld(widths[grid], maturecollision());
        // Long intervals make the realized counts large enough to compare.
        for (std::uint64_t seed = 1; seed <= 20; ++seed)
        {
            const FormationSample sample = sampleformation(fixture.view(), parameters, 200.0, seed);
            activearea[grid] = sample.orogenicAreaKm2;
            gold[grid] += totalgold(sample);
            for (const NewDeposit& item : sample.deposits)
                count[grid] += item.deposit.system == DepositSystem::orogenic_gold;
        }
    }

    // The kernel integral changes only by discretisation.
    EXPECT_NEAR(activearea[1] / activearea[0], 1.0, 0.03);
    const double expected = parameters.orogenic.ratePerMkm2Myr * 1e-6 * activearea[0] * 200.0 * 20.0;
    ASSERT_GT(expected, 500.0);
    for (int grid = 0; grid < 2; ++grid)
        EXPECT_NEAR(count[grid], expected, 5.0 * std::sqrt(expected)) << "grid " << widths[grid];
    // Lognormal sizes are heavy tailed; totals agree to sampling error.
    EXPECT_NEAR(gold[1] / gold[0], 1.0, 0.35);
}

TEST(GoldFormation, FixedSeedsReproduceDeposits)
{
    FrameFixture fixture = boxworld(128, maturecollision());
    FormationParameters parameters;
    const FormationSample first = sampleformation(fixture.view(), parameters, 100.0, 42);
    const FormationSample second = sampleformation(fixture.view(), parameters, 100.0, 42);
    const FormationSample other = sampleformation(fixture.view(), parameters, 100.0, 43);

    ASSERT_FALSE(first.deposits.empty());
    ASSERT_EQ(first.deposits.size(), second.deposits.size());
    for (std::size_t i = 0; i < first.deposits.size(); ++i)
    {
        EXPECT_EQ(first.deposits[i].x, second.deposits[i].x);
        EXPECT_EQ(first.deposits[i].y, second.deposits[i].y);
        EXPECT_EQ(first.deposits[i].deposit.oreMt, second.deposits[i].deposit.oreMt);
        EXPECT_EQ(first.deposits[i].deposit.goldGradeGt, second.deposits[i].deposit.goldGradeGt);
        EXPECT_EQ(first.deposits[i].deposit.formation.timeMyr, second.deposits[i].deposit.formation.timeMyr);
    }
    EXPECT_NE(totalgold(first), totalgold(other));
}

TEST(GoldFormation, DeepOceanFloorHostsNoPrimaryGold)
{
    FrameFixture fixture = boxworld(64, maturecollision());
    for (float& elevation : fixture.elevation)
        elevation = -4000.0f;
    const FormationSample sample = sampleformation(fixture.view(), FormationParameters{}, 1000.0, 7);
    EXPECT_EQ(sample.orogenicAreaKm2, 0.0);
    EXPECT_TRUE(sample.deposits.empty());
}

TEST(GoldFormation, FormationRecordsTheSetting)
{
    FrameFixture fixture = boxworld(128, maturecollision());
    fixture.frame.timeMyr = 300.0;
    const FormationSample sample = sampleformation(fixture.view(), FormationParameters{}, 50.0, 3);
    ASSERT_FALSE(sample.deposits.empty());
    for (const NewDeposit& item : sample.deposits)
    {
        const FormationRecord& record = item.deposit.formation;
        EXPECT_EQ(item.deposit.origin, DepositOrigin::simulated);
        EXPECT_EQ(record.geologicRegime, regimecode::continent_collision);
        EXPECT_EQ(record.boundarySegmentId, 1u);
        EXPECT_EQ(record.hostCrustClass, crustcode::continental);
        EXPECT_GE(record.timeMyr, 250.0);
        EXPECT_LE(record.timeMyr, 300.0);
        EXPECT_GT(record.settingIntensity, 0.0f);
        EXPECT_GT(item.deposit.containedGoldT(), 0.0);
    }
}

TEST(GoldFormation, InheritedProvincesLieOnContinentalCrust)
{
    FrameFixture fixture = boxworld(256, BoundarySummary{});
    FormationParameters parameters;
    parameters.provincesPerMkm2 = 20.0; // Dense enough to sample several provinces.
    std::uint32_t provinces = 0;
    const auto deposits = sampleinheritedprovinces(fixture.view(), parameters, 11, provinces);

    ASSERT_GT(provinces, 0u);
    ASSERT_FALSE(deposits.empty());
    for (const NewDeposit& item : deposits)
    {
        const auto cell = fixture.cell(static_cast<std::uint32_t>(item.x), static_cast<std::uint32_t>(item.y));
        EXPECT_EQ(fixture.crust[cell], crustcode::continental);
        EXPECT_EQ(item.deposit.origin, DepositOrigin::inherited);
        EXPECT_GE(item.deposit.formation.provinceId, 1u);
        EXPECT_LE(item.deposit.formation.provinceId, provinces);
        EXPECT_LT(item.deposit.formation.timeMyr, 0.0);
        EXPECT_GE(item.deposit.topDepthM(), parameters.inheritedTopDepthM.minimum - 1e-9);
        EXPECT_LE(item.deposit.topDepthM(), parameters.inheritedTopDepthM.maximum + 1e-9);
    }
}
