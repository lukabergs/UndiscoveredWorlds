#include "gtest/gtest.h"
#include "simulations/geology/plate_tectonics/lithosphere.hpp"
#include "simulations/geology/mineral_systems/mineral_history.hpp"
#include "simulations/geology/plate_tectonics/tectonic_api.hpp"
#include "simulations/geology/mineral_systems/tectonic_frame_reader.hpp"
#include "simulations/geology/plate_tectonics/tectonic_recipe.hpp"
#include <cmath>

using namespace mineralsystems;

TEST(TectonicMineralIntegration, ReadsColumnDisplacementAndDeterministicTopologyIdentity)
{
    auto scenario = tectonics::scenario::baseline();
    scenario.seed = 1;
    scenario.sphere_face_cells = 16;  // the frame is the id-order view of the cube, N x 6N
    scenario.evolution_duration_myr = 12;
    lithosphere a(scenario), b(scenario);
    // The frame's width in km: N mean cell sides (world_tectonics.cpp's choice until W2-R).
    const double pi = 3.14159265358979323846;
    TectonicFrameReader reader(scenario.radius_km * std::sqrt(2 * pi / 3));
    const double meanarea = 4 * pi * scenario.radius_km * scenario.radius_km /
                            (6.0 * scenario.sphere_face_cells * scenario.sphere_face_cells);
    MineralHistory minerals({}, 1);
    minerals.begin(reader.read(&a));
    EXPECT_EQ(reader.read(&a).materialShift, nullptr);
    while (!a.isFinished())
    {
        a.update(); b.update();
        const auto& frame = reader.read(&a);
        ASSERT_NE(frame.materialShift, nullptr);
        EXPECT_DOUBLE_EQ(frame.timeMyr, a.getTimeMyr());
        EXPECT_NEAR(frame.cellAreaKm2(), meanarea, meanarea * 1e-12);
        ASSERT_EQ(frame.plates.size(), a.getTopologyDiagnostics().stable_ids.size());
        for (std::uint32_t p = 0; p < frame.plates.size(); ++p)
        {
            EXPECT_EQ(frame.plates[p].uid, a.getTopologyDiagnostics().stable_ids[p]);
            EXPECT_EQ(frame.plates[p].uid, tectonics_get_plate_uid(&b, p));
        }
        ASSERT_NE(frame.continentalCrustM, nullptr);
        const auto& columns = a.getColumnModel().cells();
        for (std::size_t cell = 0; cell < frame.cellCount(); ++cell)
        {
            ASSERT_FLOAT_EQ(frame.continentalCrustM[cell], static_cast<float>(columns[cell].continent));
            ASSERT_DOUBLE_EQ(frame.materialShift[cell].x, a.getMaterialDisplacement()[cell].x);
            ASSERT_DOUBLE_EQ(frame.materialShift[cell].y, a.getMaterialDisplacement()[cell].y);
            ASSERT_TRUE(std::isfinite(frame.elevationM[cell]));
        }
        minerals.step(frame);
        for (const auto& marker : minerals.markers())
        {
            ASSERT_TRUE(std::isfinite(marker.x) && std::isfinite(marker.y));
            ASSERT_GE(marker.x, 0); ASSERT_LT(marker.x, frame.width);
            ASSERT_GE(marker.y, 0); ASSERT_LT(marker.y, frame.height);
        }
    }
    minerals.finish(reader.read(&a));
    ASSERT_FALSE(minerals.deposits().empty());
    double supplied = minerals.budget().inheritedGoldT;
    for (double gold : minerals.budget().formedGoldT) supplied += gold;
    double remaining = 0;
    for (const auto& deposit : minerals.deposits())
        if (deposit.fate == DepositFate::preserved) remaining += deposit.remainingGoldT();
    EXPECT_NEAR(supplied, remaining + minerals.budget().erodedGoldT + minerals.budget().subductedGoldT,
                supplied * 1e-10);
}
