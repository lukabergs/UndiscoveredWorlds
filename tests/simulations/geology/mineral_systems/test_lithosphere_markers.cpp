#include "gtest/gtest.h"
#include "simulations/geology/plate_tectonics/lithosphere.hpp"
#include "simulations/geology/mineral_systems/mineral_history.hpp"
#include "simulations/geology/plate_tectonics/tectonic_api.hpp"
#include "simulations/geology/mineral_systems/tectonic_frame_reader.hpp"
#include "simulations/geology/plate_tectonics/tectonic_recipe.hpp"
#include <algorithm>
#include <cmath>
#include <string>

using namespace mineralsystems;

TEST(TectonicMineralIntegration, ReadsColumnDisplacementAndDeterministicTopologyIdentity)
{
    auto scenario = tectonics::scenario::baseline();
    scenario.seed = 1;
    scenario.sphere_face_cells = 16;  // the frame is the engine's cube
    scenario.evolution_duration_myr = 12;
    lithosphere a(scenario), b(scenario);
    TectonicFrameReader reader;
    const double pi = 3.14159265358979323846;
    const double surface = 4 * pi * scenario.radius_km * scenario.radius_km;
    const double spacing_m = a.getGrid().cellspacing_m();
    MineralHistory minerals({}, 1);
    minerals.begin(reader.read(&a));
    EXPECT_EQ(reader.read(&a).materialShift, nullptr);
    while (!a.isFinished())
    {
        a.update(); b.update();
        const auto& frame = reader.read(&a);
        ASSERT_NE(frame.materialShift, nullptr);
        EXPECT_DOUBLE_EQ(frame.timeMyr, a.getTimeMyr());
        // True areas: the cells sum to the sphere; the spacing is the mean cell side.
        ASSERT_EQ(frame.grid->kind(), mesh::GridKind::cube);
        ASSERT_EQ(frame.cellCount(), a.getGrid().cellcount());
        double area = 0;
        for (std::size_t cell = 0; cell < frame.cellCount(); ++cell) area += frame.cellAreaKm2(cell);
        EXPECT_NEAR(area, surface, surface * 1e-9);
        EXPECT_NEAR(frame.spacingKm * 1000.0, spacing_m, spacing_m * 1e-12);
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
            // The published (east, north) cell spacings as a tangent vector in metres.
            const auto c = static_cast<mesh::CellId>(cell);
            const mesh::Vec3 shift = frame.materialShift[cell];
            ASSERT_NEAR(mesh::dot(shift, frame.grid->centre(c)), 0.0, 1e-6 * (1.0 + mesh::length(shift)));
            ASSERT_NEAR(mesh::dot(shift, frame.grid->east(c)), a.getMaterialDisplacement()[cell].x * spacing_m, 1e-6);
            ASSERT_NEAR(mesh::dot(shift, frame.grid->north(c)), a.getMaterialDisplacement()[cell].y * spacing_m, 1e-6);
            ASSERT_TRUE(std::isfinite(frame.elevationM[cell]));
        }
        minerals.step(frame);
        for (const auto& marker : minerals.markers())
        {
            ASSERT_NEAR(mesh::length(marker.position), 1.0, 1e-9);
        }
    }
    minerals.finish(reader.read(&a));
    ASSERT_FALSE(minerals.deposits().empty());
    // Positions as fractions of the sphere's longitude and colatitude.
    for (const auto& deposit : minerals.deposits())
    {
        EXPECT_GE(deposit.u, 0.0); EXPECT_LT(deposit.u, 1.0);
        EXPECT_GE(deposit.v, 0.0); EXPECT_LE(deposit.v, 1.0);
    }
    double supplied = minerals.budget().inheritedGoldT;
    for (double gold : minerals.budget().formedGoldT) supplied += gold;
    double remaining = 0;
    for (const auto& deposit : minerals.deposits())
        if (deposit.fate == DepositFate::preserved) remaining += deposit.remainingGoldT();
    EXPECT_NEAR(supplied, remaining + minerals.budget().erodedGoldT + minerals.budget().subductedGoldT,
                supplied * 1e-10);
}

// One engine update (T8 E5 review): away from plate boundaries, beyond the reach of the
// distributed deformation, a cell's published displacement is its plate's rigid motion, so
// shift - (omega dt) x R c is no larger than the (vanishing) deformation offset there. A unit
// error in omega, dt or the spacing of the published displacement shows as a residual of the
// order of the rigid motion itself.
TEST(TectonicMineralIntegration, DisplacementAwayFromBoundariesIsThePlateRotation)
{
    auto scenario = tectonics::scenario::baseline();
    scenario.seed = 1;
    scenario.sphere_face_cells = 16;
    scenario.evolution_duration_myr = 4;
    lithosphere engine(scenario);
    engine.update();
    TectonicFrameReader reader;
    const TectonicFrame& frame = reader.read(&engine);
    ASSERT_NE(frame.materialShift, nullptr);
    ASSERT_EQ(frame.plateRotation.size(), frame.plates.size());
    const double radius_m = frame.grid->radius_m();
    const double spacing_km = frame.spacingKm;
    // The deformation's reach: twice its peak distance, plus two cells.
    const double reach_km = 2.0 * scenario.deformation_peak_km + 2.0 * spacing_km;
    int interior = 0;
    double largestresidual = 0.0, largestrigid = 0.0;
    for (std::size_t cell = 0; cell < frame.cellCount(); ++cell)
    {
        const std::uint32_t owner = frame.plateIndex[cell];
        if (owner >= frame.plateRotation.size() || frame.boundaryDistance[cell] * spacing_km <= reach_km)
            continue;
        const mesh::Vec3 rigid = mesh::cross(frame.plateRotation[owner], radius_m * frame.grid->centre(static_cast<mesh::CellId>(cell)));
        const double residual = mesh::length(frame.materialShift[cell] - rigid);
        ASSERT_LE(residual, 1e-3 * mesh::length(rigid) + 1e-9 * spacing_km * 1000.0) << "cell " << cell;
        largestresidual = std::max(largestresidual, residual);
        largestrigid = std::max(largestrigid, mesh::length(rigid));
        ++interior;
    }
    RecordProperty("interior_cells", interior);
    RecordProperty("largest_residual_m", std::to_string(largestresidual));
    RecordProperty("largest_rigid_m", std::to_string(largestrigid));
    EXPECT_GT(interior, 50);
    EXPECT_GT(largestrigid, 1.0);  // the plates move (metres per update)
}
