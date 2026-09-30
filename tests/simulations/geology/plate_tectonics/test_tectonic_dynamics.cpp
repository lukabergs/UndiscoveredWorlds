#include "gtest/gtest.h"
#include "tectonic_dynamics.hpp"

#include "tectonic_provenance.hpp"
#include "lithosphere.hpp"
#include "plate.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>

namespace {
using tectonics::dynamics::Contact;
using tectonics::dynamics::Plate;
using tectonics::dynamics::solve;

void check_budget(const tectonics::dynamics::Solution& result) {
    EXPECT_LT(result.diagnostics.relative_force_residual, 1e-12);
    EXPECT_GE(result.diagnostics.dissipation_w, 0.0);
    EXPECT_NEAR(result.diagnostics.driving_power_w, result.diagnostics.dissipation_w,
                std::max(1.0, result.diagnostics.dissipation_w) * 1e-12);
}
}

TEST(TectonicDynamics, BasalDragHasAnalyticalEquilibriumAndNoMinimumSpeed) {
    const auto moving = solve({Plate{1e13, {1e-9, -2e-9}}}, {});
    EXPECT_NEAR(moving.velocity_m_s[0].x, 1e-9, 1e-23);
    EXPECT_NEAR(moving.velocity_m_s[0].y, -2e-9, 1e-23);
    check_budget(moving);
    const auto stopped = solve({Plate{1e13, {0, 0}}}, {});
    EXPECT_DOUBLE_EQ(stopped.velocity_m_s[0].x, 0);
    EXPECT_DOUBLE_EQ(stopped.velocity_m_s[0].y, 0);
    check_budget(stopped);
    EXPECT_THROW(solve({Plate{0, {0, 0}}}, {}), std::runtime_error);
}

TEST(TectonicDynamics, OceanRidgeSpreadsWithBalancedOpposingForces) {
    const std::vector<Plate> plates = {{1e13, {}}, {1e13, {}}};
    const Contact ridge{0, 1, {1, 0}, 1e6, false, false, 0, 0, -1e-9};
    const auto result = solve(plates, {ridge});
    const double expected_speed = 2e12 * 1e6 /
        (tectonics::dynamics::kBasalDrag * 1e13 + 2 * tectonics::dynamics::kOceanBoundaryDrag * 1e6);
    EXPECT_NEAR(result.velocity_m_s[0].x, -expected_speed, 1e-23);
    EXPECT_NEAR(result.velocity_m_s[1].x, expected_speed, 1e-23);
    EXPECT_DOUBLE_EQ(result.diagnostics.ridge_length_km, 1000);
    check_budget(result);
    auto split = ridge;
    split.length_m /= 2;
    const auto refined = solve(plates, {split, split});
    EXPECT_NEAR(refined.velocity_m_s[0].x, result.velocity_m_s[0].x, 1e-23);
}

TEST(TectonicDynamics, OldOceanSlabPullsTowardTrenchWithoutPullingContinentalMaterialDown) {
    const std::vector<Plate> plates = {{1e13, {}}, {1e13, {}}};
    Contact trench{0, 1, {1, 0}, 1e6, false, true, 80, 200, 1e-9};
    const auto old = solve(plates, {trench});
    EXPECT_GT(old.velocity_m_s[0].x, old.velocity_m_s[1].x);
    EXPECT_GT(old.velocity_m_s[1].x, 0); // Coupled by boundary drag, not slab buoyancy.
    trench.age_a_myr = 5;
    const auto young = solve(plates, {trench});
    EXPECT_GT(old.velocity_m_s[0].x, young.velocity_m_s[0].x);
    trench.continental_a = true;
    const auto continent = solve(plates, {trench});
    EXPECT_DOUBLE_EQ(continent.velocity_m_s[0].x, 0);
    EXPECT_DOUBLE_EQ(continent.diagnostics.subduction_length_km, 0);
    check_budget(old);
}

TEST(TectonicDynamics, ContinentalResistanceSlowsConvergenceAndDissipatesWork) {
    const std::vector<Plate> plates = {{1e13, {1e-9, 0}}, {1e13, {-1e-9, 0}}};
    const Contact collision{0, 1, {1, 0}, 1e6, true, true, 200, 200, 1e-9};
    const auto free = solve(plates, {}), resisted = solve(plates, {collision});
    EXPECT_GT(resisted.velocity_m_s[0].x, 0);
    EXPECT_LT(resisted.velocity_m_s[0].x, free.velocity_m_s[0].x);
    EXPECT_NEAR(resisted.velocity_m_s[0].x, -resisted.velocity_m_s[1].x, 1e-23);
    check_budget(resisted);
}

TEST(TectonicDynamics, ThermalBuoyancyScalingPreservesForceAndPowerBalance) {
    const std::vector<Plate> plates={{1e13,{}},{1e13,{}}};
    Contact trench{0,1,{1,0},1e6,false,true,80,200,1e-9};
    const auto full=solve(plates,{trench});
    trench.slab_pull_scale=.5;
    const auto half=solve(plates,{trench});
    EXPECT_NEAR(half.velocity_m_s[0].x,.5*full.velocity_m_s[0].x,1e-23);
    check_budget(half);
    trench.slab_pull_scale=0;
    const auto hot=solve(plates,{trench});
    EXPECT_DOUBLE_EQ(hot.velocity_m_s[0].x,0);
    EXPECT_DOUBLE_EQ(hot.velocity_m_s[1].x,0);
    trench.slab_pull_scale=-1;
    EXPECT_THROW(solve(plates,{trench}),std::runtime_error);
}

TEST(TectonicDynamics, MountainLoadOpposesConvergenceAndCanDriveSpreadingWithoutAHeightCap) {
    const std::vector<Plate> plates = {{1e13, {1e-9, 0}}, {1e13, {-1e-9, 0}}};
    Contact collision{0,1,{1,0},1e6,true,true,200,200,1e-9};
    const auto flat = solve(plates,{collision});
    collision.mountain_force_n_m = 3e12;
    const auto mountain = solve(plates,{collision});
    EXPECT_GT(mountain.velocity_m_s[0].x,0);
    EXPECT_LT(mountain.velocity_m_s[0].x,flat.velocity_m_s[0].x);
    EXPECT_NEAR(mountain.velocity_m_s[0].x,-mountain.velocity_m_s[1].x,1e-23);
    collision.mountain_force_n_m = 8e12;
    const auto spreading = solve(plates,{collision});
    EXPECT_LT(spreading.velocity_m_s[0].x,0);
    EXPECT_DOUBLE_EQ(spreading.diagnostics.mountain_boundary_force_n,8e18);
    check_budget(mountain); check_budget(spreading);
}

TEST(TectonicDynamics, CoastalMountainRootResistsShorteningWithoutChangingNetExternalForce) {
    const std::vector<Plate> plates={{1e13,{1e-9,0}},{1e13,{-1e-9,0}}};
    Contact margin{0,1,{1,0},1e6,false,true,80,200,1e-9};
    const auto flat=solve(plates,{margin});
    margin.mountain_force_n_m=3e12;
    const auto loaded=solve(plates,{margin});
    EXPECT_LT(loaded.velocity_m_s[0].x-loaded.velocity_m_s[1].x,
              flat.velocity_m_s[0].x-flat.velocity_m_s[1].x);
    const double delta=3e18/(tectonics::dynamics::kBasalDrag*1e13+2*tectonics::dynamics::kOceanBoundaryDrag*1e6);
    EXPECT_NEAR(loaded.velocity_m_s[0].x,flat.velocity_m_s[0].x-delta,1e-23);
    EXPECT_NEAR(loaded.velocity_m_s[1].x,flat.velocity_m_s[1].x+delta,1e-23);
    EXPECT_DOUBLE_EQ(loaded.diagnostics.mountain_boundary_force_n,3e18);
    margin.a=1;margin.b=0;margin.normal.x=-1;
    margin.continental_a=true;margin.continental_b=false;
    std::swap(margin.age_a_myr,margin.age_b_myr);
    const auto reversed=solve(plates,{margin});
    EXPECT_NEAR(reversed.velocity_m_s[0].x,loaded.velocity_m_s[0].x,1e-23);
    check_budget(loaded);check_budget(reversed);
}

TEST(TectonicDynamics, PrescribedMantleFlowIsSeededPeriodicAndSmoothInPhysicalTime) {
    using tectonics::dynamics::mantle_flow;
    const auto a = mantle_flow(.2, .7, 50, 1, 4);
    const auto same = mantle_flow(1.2, -.3, 650, 1, 4);
    const auto next = mantle_flow(.2, .7, 50.001, 1, 4);
    const auto other = mantle_flow(.2, .7, 50, 2, 4);
    EXPECT_NEAR(a.x, same.x, 1e-23);
    EXPECT_NEAR(a.y, same.y, 1e-23);
    EXPECT_LT(std::hypot(a.x - next.x, a.y - next.y), 1e-13);
    EXPECT_GT(std::hypot(a.x - other.x, a.y - other.y), 1e-13);
}

TEST(TectonicDynamics, RidgeAndSlabForcesVanishContinuouslyAsNormalMotionStops) {
    using namespace tectonics::dynamics;
    std::vector<Plate> plates(2, {1e12, {0, 0}});
    Contact edge{0, 1, {1, 0}, 1e6, false, false, 80, 10};
    edge.activation_speed_m_s = .1 / 3.15576e9;
    for (double sign : {-1., 1.}) {
        edge.closing_m_s = sign * edge.activation_speed_m_s * 1e-3;
        const auto first = solve(plates, {edge});
        edge.closing_m_s *= .1;
        const auto smaller = solve(plates, {edge});
        EXPECT_NEAR(smaller.diagnostics.max_speed_cm_yr / first.diagnostics.max_speed_cm_yr, .1, 1e-6);
        EXPECT_LT(smaller.diagnostics.relative_force_residual, 1e-12);
        EXPECT_NEAR(smaller.diagnostics.driving_power_w, smaller.diagnostics.dissipation_w, 1e-6);
    }
    edge.closing_m_s = 0;
    EXPECT_DOUBLE_EQ(solve(plates, {edge}).diagnostics.max_speed_cm_yr, 0);
    edge.activation_speed_m_s = -1;
    EXPECT_THROW(solve(plates, {edge}), std::runtime_error);
}

TEST(TectonicDynamics, BoundaryActivityUsesPhysicalSpeedIndependentOfTimestep) {
    std::vector<float> height(8, 1.1f);
    std::vector<uint32_t> owners = {0, 0, 1, 1, 0, 0, 1, 1};
    std::vector<tectonics::contract::PlateKinematics> plates(2);
    plates[0].linear_velocity_x = .02f;
    plates[1].linear_velocity_x = -.02f;
    tectonics::provenance::Inputs inputs{4, 2, height.data(), owners.data(), plates.data(), 2, nullptr, 25.0};
    std::vector<uint8_t> convergence(8), divergence(8), shear(8), regime(8);
    tectonics::provenance::compute_maps(inputs, convergence.data(), divergence.data(), shear.data(), regime.data());
    EXPECT_EQ(convergence[1], 20U); // Relative convergence = 1 cm/year.
    EXPECT_EQ(regime[1], static_cast<uint8_t>(tectonics::contract::GeologicRegime::ContinentCollision));
    const auto before = convergence;
    plates[0].linear_velocity_x *= .25f;
    plates[1].linear_velocity_x *= .25f;
    inputs.velocity_to_cm_yr *= 4;
    tectonics::provenance::compute_maps(inputs, convergence.data(), divergence.data(), shear.data(), regime.data());
    EXPECT_EQ(before, convergence);
}

TEST(TectonicDynamics, ContinuousStationaryWorldPreservesCrustUntilDurationEnds) {
    tectonics::scenario::Scenario s;
    s.seed = 1; s.width = 32; s.height = 16; s.plate_count = 1;
    s.evolution_duration_myr = 6; 
    s.mantle_flow_speed_cm_yr = 0; s.thermal_subsidence_strength = 0;
    s.delta_time_myr = .5; s.erosion_strength = 0; lithosphere world(s);
    const auto datum = world.getSeaLevelMeters();
    std::vector<uint16_t> terrain(32 * 16, datum + 1500);
    world.importMetricHeightMap(terrain.data(), datum);
    const std::vector<float> initial(world.getTopography(), world.getTopography() + terrain.size());
    const std::vector<float> age(world.getCrustAgeMyrMap(), world.getCrustAgeMyrMap() + terrain.size());
    int calls = 0;
    while (!world.isFinished() && calls < 20) { world.update(); ++calls; }
    ASSERT_EQ(calls, 12);
    EXPECT_DOUBLE_EQ(world.getTimeMyr(), 6);
    EXPECT_EQ(world.getCycleCount(), 0U);
    EXPECT_EQ(std::memcmp(initial.data(), world.getTopography(), initial.size() * sizeof(float)), 0);
    for (size_t i = 0; i < age.size(); ++i) EXPECT_FLOAT_EQ(world.getCrustAgeMyrMap()[i], age[i] + 6);
    world.update(); // Finished simulation is a no-op, including force solve.
    EXPECT_DOUBLE_EQ(world.getTimeMyr(), 6);
}

TEST(TectonicDynamics, ContinuousDrivenWorldReplaysAndRetainsMotionAt200Myr) {
    tectonics::scenario::Scenario s;
    s.seed = 1; s.width = 64; s.height = 32; s.plate_count = 8;
    s.evolution_duration_myr = 200; 
    lithosphere a(s), b(s);
    for (int i = 0; i < 200; ++i) {
        ASSERT_FALSE(a.isFinished()) << i;
        a.update(); b.update();
        ASSERT_LT(a.getMotionDiagnostics().relative_force_residual, 1e-12) << i;
        ASSERT_EQ(a.getCycleCount(), 0U) << i;
    }
    EXPECT_TRUE(a.isFinished());
    EXPECT_DOUBLE_EQ(a.getTimeMyr(), 200);
    EXPECT_GT(a.getMotionDiagnostics().mean_speed_cm_yr, .01);
    EXPECT_EQ(std::memcmp(a.getTopography(), b.getTopography(), 64 * 32 * sizeof(float)), 0);
    EXPECT_EQ(std::memcmp(a.getAgeMap(), b.getAgeMap(), 64 * 32 * sizeof(uint32_t)), 0);
}

TEST(TectonicDynamics, DrivenContinentalCollisionRaisesInitiallyFlatTerrain) {
    tectonics::scenario::Scenario s;
    s.seed = 1; s.width = 64; s.height = 32; s.plate_count = 8;
    s.evolution_duration_myr = 30; 
    s.world_width_km = 4000; s.delta_time_myr = .25;
    s.erosion_strength = 0; auto disabled = s; disabled.mantle_flow_speed_cm_yr = 0; lithosphere active(s), control(disabled);
    const auto datum = active.getSeaLevelMeters();
    std::vector<uint16_t> plateau(64 * 32, datum + 1500);
    active.importMetricHeightMap(plateau.data(), datum);
    control.importMetricHeightMap(plateau.data(), datum);
    while (!active.isFinished()) { active.update(); control.update(); }
    const auto peak = [datum](const lithosphere& world) {
        return static_cast<int>(TopographyCodec::internal_to_meters(
            *std::max_element(world.getTopography(), world.getTopography() + 64 * 32), datum)) - datum;
    };
    EXPECT_GT(peak(active), 1600);
    EXPECT_LE(peak(control), 1501);
}

TEST(TectonicDynamics, InlandPressureTransmitsRootsButStopsAtOtherPlates) {
    using tectonics::dynamics::inland_pressure;
    constexpr uint32_t w=8,h=4;
    std::vector<uint32_t> owners(w*h,0);
    std::vector<double> pressure(w*h,2e12);
    EXPECT_NEAR(inland_pressure(w,h,100000,owners,pressure,3,0,-1,500000),2e12,1e-3);
    pressure[3]=0;pressure[1]=1e13;
    EXPECT_DOUBLE_EQ(inland_pressure(w,h,100000,owners,pressure,3,0,-1,0),0);
    EXPECT_GT(inland_pressure(w,h,100000,owners,pressure,3,0,-1,500000),2e12);
    owners[2]=1;
    EXPECT_DOUBLE_EQ(inland_pressure(w,h,100000,owners,pressure,3,0,-1,500000),0);
    owners.assign(w*h,0);pressure.assign(w*h,0);pressure[7]=1e13;
    EXPECT_GT(inland_pressure(w,h,100000,owners,pressure,0,0,-1,300000),0);
    EXPECT_THROW(inland_pressure(w,h,100000,owners,pressure,0,0,0,300000),std::runtime_error);
}

TEST(TectonicDynamics, BoundaryHistoryBalancesItsReturnedVelocity) {
    using namespace tectonics::dynamics;
    std::vector<Plate> plates{{2e13,{1e-9,0}},{3e13,{-1e-9,0}}};
    Contact contact{0,1,{1,0},2e6,false,false,80,20};
    contact.activation_speed_m_s=.1/3.15576e9;
    contact.descending_side=0;
    contact.ridge_activity=.2;contact.slab_activity=.8;contact.activity_response=1./31;
    const auto solution=solve_consistent(plates,{contact},{{0,0},{0,0}});
    contact.closing_m_s=solution.velocity_m_s[0].x-solution.velocity_m_s[1].x;
    const auto check=solve(plates,{contact});
    EXPECT_LT(solution.diagnostics.velocity_consistency_cm_yr,1e-4);
    for (size_t p=0;p<plates.size();++p) {
        EXPECT_NEAR(check.velocity_m_s[p].x,solution.velocity_m_s[p].x,1e-4/3.15576e9);
        EXPECT_NEAR(check.velocity_m_s[p].y,solution.velocity_m_s[p].y,1e-4/3.15576e9);
    }
    // A reversed surface velocity does not instantaneously detach an old slab.
    contact.closing_m_s=-1e-9;
    EXPECT_GT(solve(plates,{contact}).diagnostics.subduction_length_km,0);
}
