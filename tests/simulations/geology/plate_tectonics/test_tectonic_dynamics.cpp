#include "gtest/gtest.h"
#include "simulations/geology/plate_tectonics/tectonic_dynamics.hpp"

#include "simulations/geology/plate_tectonics/tectonic_provenance.hpp"
#include "simulations/geology/plate_tectonics/lithosphere.hpp"

#include "core/mesh/cube_geometry.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <random>

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

// ---- The sphere (W2-E, step E2C) ----
namespace {
using tectonics::dynamics::MantleSample;
using tectonics::dynamics::SphereContact;
using tectonics::dynamics::SpherePlate;
using tectonics::grid::TectonicGrid;
using Vec3 = mesh::Vec3;
constexpr double kEarthRadiusKm = 6371.0;
constexpr double kCmYr = 3.15576e9;  // cm/yr per m/s

// A random partition: count distinct source cells from mt19937(seed), each cell to its
// nearest source over the edge graph.
std::vector<uint32_t> random_partition(const TectonicGrid& grid, uint32_t count, uint32_t seed) {
    std::mt19937 engine(seed);
    std::vector<tectonics::grid::CellId> sources;
    while (sources.size() < count) {
        const auto c = static_cast<tectonics::grid::CellId>(engine() % grid.cellcount());
        if (std::find(sources.begin(), sources.end(), c) == sources.end()) sources.push_back(c);
    }
    return grid.dijkstra(sources).source;
}

std::vector<SphereContact> plate_contacts(const TectonicGrid& grid, const std::vector<uint32_t>& owners) {
    std::vector<SphereContact> contacts;
    for (const auto& edge : grid.edges())
        if (owners[edge.a] != owners[edge.b])
            contacts.push_back(tectonics::dynamics::sphere_contact(grid, edge, owners[edge.a], owners[edge.b]));
    return contacts;
}

// Cyclic Jacobi eigenvalues of a symmetric n x n matrix (scaled to its largest diagonal).
std::vector<double> eigenvalues(std::vector<double> a, size_t n) {
    double scale = 0;
    for (size_t i = 0; i < n; ++i) scale = std::max(scale, std::abs(a[i * n + i]));
    for (auto& value : a) value /= scale;
    for (int sweep = 0; sweep < 100; ++sweep) {
        double off = 0;
        for (size_t p = 0; p < n; ++p) for (size_t q = p + 1; q < n; ++q) off += a[p * n + q] * a[p * n + q];
        if (off < 1e-40) break;
        for (size_t p = 0; p < n; ++p) for (size_t q = p + 1; q < n; ++q) {
            const double apq = a[p * n + q];
            if (std::abs(apq) < 1e-300) continue;
            const double theta = (a[q * n + q] - a[p * n + p]) / (2 * apq);
            const double t = (theta >= 0 ? 1.0 : -1.0) / (std::abs(theta) + std::sqrt(theta * theta + 1));
            const double c = 1 / std::sqrt(t * t + 1), s = t * c;
            for (size_t k = 0; k < n; ++k) {
                const double akp = a[k * n + p], akq = a[k * n + q];
                a[k * n + p] = c * akp - s * akq;
                a[k * n + q] = s * akp + c * akq;
            }
            for (size_t k = 0; k < n; ++k) {
                const double apk = a[p * n + k], aqk = a[q * n + k];
                a[p * n + k] = c * apk - s * aqk;
                a[q * n + k] = s * apk + c * aqk;
            }
        }
    }
    std::vector<double> result(n);
    for (size_t i = 0; i < n; ++i) result[i] = a[i * n + i] * scale;
    std::sort(result.begin(), result.end());
    return result;
}

double relative_difference(const Vec3& a, const Vec3& b) {
    return mesh::length(a - b) / std::max(mesh::length(a), mesh::length(b));
}
} // namespace

TEST(TectonicDynamics, SphereSolidBodyMantleRotationGivesThatRotationForEveryPlate) {
    using namespace tectonics::dynamics;
    const auto grid = TectonicGrid::cube(8, kEarthRadiusKm);
    const double radius = grid->mesh()->radius_m();
    const auto owners = random_partition(*grid, 5, 1);
    const Vec3 rotation{1e-15, -2e-15, 3e-15};  // about 7.5 cm/yr at the equator of the axis
    std::vector<MantleSample> mantle(grid->cellcount());
    for (tectonics::grid::CellId c = 0; c < grid->cellcount(); ++c) {
        const Vec3 p = grid->position(c);
        mantle[c] = {mesh::cross(rotation, radius * p), mesh::dot(rotation, p)};
    }
    const auto plates = sphere_plates(*grid, owners, 5, mantle);
    auto contacts = plate_contacts(*grid, owners);
    ASSERT_FALSE(contacts.empty());
    for (auto& contact : contacts) {
        contact.closing_m_s = sphere_closing_m_s(contact, rotation, rotation, radius);
        EXPECT_EQ(contact.closing_m_s, 0.0);
    }
    const auto solution = solve_sphere(plates, contacts, radius);
    for (const auto& w : solution.rotation_rad_s) EXPECT_LT(relative_difference(w, rotation), 1e-10);
    EXPECT_LT(solution.diagnostics.relative_force_residual, 1e-12);
    EXPECT_NEAR(solution.diagnostics.driving_power_w, solution.diagnostics.dissipation_w,
                solution.diagnostics.dissipation_w * 1e-10);
}

TEST(TectonicDynamics, SphereSymmetricPlatePairGetsEqualAndOppositeTorques) {
    using namespace tectonics::dynamics;
    // Plate b is plate a turned by pi about the contact's midpoint x = +X, which maps the
    // contact normal +Y to -Y: the same contact seen from b.
    const double radius = kEarthRadiusKm * 1000;
    const auto turned = [](const Vec3& v) { return Vec3{v.x, -v.y, -v.z}; };
    std::vector<SpherePlate> plates(2);
    for (int i = 0; i < 4; ++i) for (int j = 0; j < 3; ++j) {
        const Vec3 c = mesh::normalise({1.0, -0.05 - 0.1 * i, 0.08 * (j - 1) + 0.03 * i});
        add_sphere_cell(plates[0], c, 4e11 + 1e10 * j, radius, {});
        add_sphere_cell(plates[1], turned(c), 4e11 + 1e10 * j, radius, {});
    }
    SphereContact contact;
    contact.a = 0; contact.b = 1;
    contact.position = {1, 0, 0};
    contact.normal = {0, 1, 0};
    contact.length_m = 1e6;
    contact.closing_m_s = -1e-9;           // spreading: ridge push on both flanks
    contact.mountain_force_n_m = 3e12;
    const auto system = sphere_system(plates, {contact}, radius);
    for (int k = 0; k < 3; ++k) EXPECT_EQ(system.rhs[k], -system.rhs[3 + k]);
    EXPECT_GT(std::abs(system.rhs[2]), 0);
    for (size_t i = 0; i < 6; ++i) for (size_t j = 0; j < 6; ++j)
        EXPECT_NEAR(system.matrix[i * 6 + j], system.matrix[j * 6 + i], std::abs(system.matrix[i * 6 + i]) * 1e-15);
    const auto solution = solve_sphere(plates, {contact}, radius);
    EXPECT_LT(relative_difference(solution.rotation_rad_s[1], turned(solution.rotation_rad_s[0])), 1e-10);
    // Spreading: the flanks separate along the normal.
    EXPECT_LT(sphere_closing_m_s(contact, solution.rotation_rad_s[0], solution.rotation_rad_s[1], radius), 0);
    EXPECT_DOUBLE_EQ(solution.diagnostics.ridge_length_km, 1000);
    EXPECT_DOUBLE_EQ(solution.diagnostics.mountain_boundary_force_n, 3e18);
    EXPECT_NEAR(solution.diagnostics.driving_power_w, solution.diagnostics.dissipation_w,
                solution.diagnostics.dissipation_w * 1e-12);
    // Slab pull acts on the sinking plate alone.
    contact.closing_m_s = 1e-9;
    contact.mountain_force_n_m = 0;
    contact.age_a_myr = 80; contact.age_b_myr = 20;
    const auto trench = sphere_system(plates, {contact}, radius);
    EXPECT_GT(std::abs(trench.rhs[2]), 0);
    for (int k = 3; k < 6; ++k) EXPECT_EQ(trench.rhs[k], 0);
    EXPECT_DOUBLE_EQ(trench.diagnostics.subduction_length_km, 1000);
}

TEST(TectonicDynamics, SphereOneCellPlateKeepsTheTorqueBalancePositiveDefinite) {
    using namespace tectonics::dynamics;
    const auto grid = TectonicGrid::cube(16, kEarthRadiusKm);
    const double radius = grid->mesh()->radius_m();
    auto owners = random_partition(*grid, 12, 1);
    owners[700] = 12;  // a plate of one cell inside another
    const SphereMantleFlow flow(1, 4);
    std::vector<MantleSample> mantle(grid->cellcount());
    for (tectonics::grid::CellId c = 0; c < grid->cellcount(); ++c) mantle[c] = flow.at(grid->position(c), 0, radius);
    const auto plates = sphere_plates(*grid, owners, 13, mantle);
    // Its drag alone: lever arm in two directions, the disc's spin moment in the third.
    const auto lone = eigenvalues({plates[12].drag_m4.begin(), plates[12].drag_m4.end()}, 3);
    const double area = grid->area_m2(700);
    EXPECT_NEAR(lone[0], area * area / (2 * 3.14159265358979323846), lone[0] * 1e-9);
    EXPECT_NEAR(lone[2], area * radius * radius, lone[2] * 1e-9);
    const auto contacts = plate_contacts(*grid, owners);
    const auto system = sphere_system(plates, contacts, radius);
    const auto values = eigenvalues(system.matrix, 39);
    EXPECT_GT(values.front(), 0);
    std::cout << "one-cell plate: condition " << values.back() / values.front() << "\n";
    const auto solution = solve_sphere(plates, contacts, radius);
    EXPECT_LT(solution.diagnostics.relative_force_residual, 1e-12);
    // Without contacts it moves with the mantle sample of its cell.
    const auto alone = solve_sphere({plates[12]}, {}, radius);
    const Vec3 c = grid->position(700), w = alone.rotation_rad_s[0];
    EXPECT_LT(relative_difference(mesh::cross(w, radius * c), mantle[700].velocity_m_s), 1e-9);
    EXPECT_NEAR(mesh::dot(w, c), mantle[700].spin_rad_s, std::abs(mantle[700].spin_rad_s) * 1e-9);
}

TEST(TectonicDynamics, SphereMantleFlowIsTangentContinuousSeededPeriodicAndSmoothInPhysicalTime) {
    using tectonics::dynamics::SphereMantleFlow;
    const double radius = kEarthRadiusKm * 1000, speed = 4 / kCmYr;
    const SphereMantleFlow flow(1, 4), again(1, 4), other(2, 4);
    const Vec3 p = mesh::normalise({.3, -.5, .8});
    const auto a = flow.at(p, 50, radius);
    EXPECT_EQ(a.velocity_m_s.x, again.at(p, 50, radius).velocity_m_s.x);
    EXPECT_EQ(a.spin_rad_s, again.at(p, 50, radius).spin_rad_s);
    EXPECT_LT(mesh::length(a.velocity_m_s - flow.at(p, 650, radius).velocity_m_s), speed * 1e-9);
    EXPECT_LT(mesh::length(a.velocity_m_s - flow.at(p, 50.001, radius).velocity_m_s), speed * 1e-3);
    EXPECT_GT(mesh::length(a.velocity_m_s - other.at(p, 50, radius).velocity_m_s), speed * 1e-2);
    EXPECT_GT(mesh::length(a.velocity_m_s - flow.at(p, 350, radius).velocity_m_s), speed * 1e-2);
    // Tangent, and the torus field's RMS speed, over the cube's cells and the 600 Myr cycle.
    const auto grid = TectonicGrid::cube(16, kEarthRadiusKm);
    double square = 0, area = 0;
    for (int t = 0; t < 6; ++t) for (tectonics::grid::CellId c = 0; c < grid->cellcount(); ++c) {
        const Vec3 q = grid->position(c);
        const auto v = flow.at(q, 100.0 * t + 25, radius).velocity_m_s;
        EXPECT_LT(std::abs(mesh::dot(v, q)), speed * 1e-12);
        square += grid->area_m2(c) * mesh::dot(v, v);
        area += grid->area_m2(c);
    }
    const double rms = std::sqrt(square / area) / (std::sqrt(.6225) * speed);
    std::cout << "mantle flow RMS / torus RMS on the N 16 cube: " << rms << "\n";
    EXPECT_NEAR(rms, 1, .05);
    // Continuous across face edges and corners: points a hair apart on different faces.
    const double hair = 1e-7;
    for (const Vec3 corner : {Vec3{1, 1, 1}, Vec3{-1, 1, -1}, Vec3{1, 0, 1}, Vec3{0, -1, -1}}) {
        const Vec3 x = mesh::normalise(corner);
        const auto centre = flow.at(x, 50, radius).velocity_m_s;
        std::vector<int> faces;
        for (const Vec3 axis : {Vec3{corner.x, 0, 0}, Vec3{0, corner.y, 0}, Vec3{0, 0, corner.z}}) {
            if (mesh::length(axis) == 0) continue;
            const Vec3 q = mesh::normalise(x + hair * axis);
            faces.push_back(mesh::cubeposition(q).face);
            EXPECT_LT(mesh::length(flow.at(q, 50, radius).velocity_m_s - centre), 100 * hair * speed);
        }
        std::sort(faces.begin(), faces.end());
        EXPECT_TRUE(std::adjacent_find(faces.begin(), faces.end()) == faces.end());  // one point on each face
    }
    // Spin is half the radial vorticity: the circulation around a small circle.
    for (const Vec3 centre : {p, mesh::normalise({1, 1, 1}), mesh::normalise({-.2, .9, -.1})}) {
        const Vec3 e1 = mesh::normalise(mesh::cross(centre, {0.3, 0.1, 1})), e2 = mesh::cross(centre, e1);
        const double angle = 1e-3;
        const int steps = 128;
        double circulation = 0;
        for (int k = 0; k < steps; ++k) {
            const double theta = 2 * 3.14159265358979323846 * (k + .5) / steps;
            const Vec3 q = std::cos(angle) * centre + std::sin(angle) * (std::cos(theta) * e1 + std::sin(theta) * e2);
            const Vec3 tangent = -std::sin(theta) * e1 + std::cos(theta) * e2;
            circulation += mesh::dot(flow.at(q, 50, radius).velocity_m_s, tangent);
        }
        circulation *= 2 * 3.14159265358979323846 * radius * std::sin(angle) / steps;
        const double cap = 2 * 3.14159265358979323846 * radius * radius * (1 - std::cos(angle));
        EXPECT_NEAR(flow.at(centre, 50, radius).spin_rad_s, circulation / (2 * cap), 1e-3 * speed / radius);
    }
}

TEST(TectonicDynamics, SphereTwelvePlatePartitionRecordsConditionAndSpinShare) {
    using namespace tectonics::dynamics;
    const auto grid = TectonicGrid::cube(16, kEarthRadiusKm);
    const double radius = grid->mesh()->radius_m();
    const auto owners = random_partition(*grid, 12, 1);
    const SphereMantleFlow flow(1, 4);
    std::vector<MantleSample> mantle(grid->cellcount());
    for (tectonics::grid::CellId c = 0; c < grid->cellcount(); ++c) mantle[c] = flow.at(grid->position(c), 0, radius);
    const auto plates = sphere_plates(*grid, owners, 12, mantle);
    // Each plate's motion under its basal drag alone, for the boundaries' initial activity
    // (the engine's rule: activity follows the mantle forcing, memory 30 Myr at 1 Myr steps).
    std::vector<Vec3> free(12);
    for (uint32_t p = 0; p < 12; ++p) free[p] = solve_sphere({plates[p]}, {}, radius).rotation_rad_s[0];
    auto contacts = plate_contacts(*grid, owners);
    for (auto& contact : contacts) {
        contact.continental_a = contact.a < 3;
        contact.continental_b = contact.b < 3;
        contact.age_a_myr = 10 + 15.0 * contact.a;
        contact.age_b_myr = 10 + 15.0 * contact.b;
        contact.activation_speed_m_s = .1 / kCmYr;
        const double closing = sphere_closing_m_s(contact, free[contact.a], free[contact.b], radius);
        const double strength = std::abs(closing) / std::hypot(closing, contact.activation_speed_m_s);
        contact.ridge_activity = closing < 0 ? strength : 0;
        contact.slab_activity = closing > 0 ? strength : 0;
        contact.activity_response = 1.0 / 31;
    }
    const auto solution = solve_sphere_consistent(plates, contacts, std::vector<Vec3>(12), radius);
    EXPECT_LT(solution.diagnostics.velocity_consistency_cm_yr, 1e-4);
    EXPECT_LT(solution.diagnostics.relative_force_residual, 1e-12);
    EXPECT_GT(solution.diagnostics.mean_speed_cm_yr, .1);
    EXPECT_GT(solution.diagnostics.ridge_length_km, 0);
    EXPECT_GT(solution.diagnostics.subduction_length_km, 0);
    const auto values = eigenvalues(sphere_system(plates, contacts, radius).matrix, 36);
    ASSERT_GT(values.front(), 0);
    const double condition = values.back() / values.front();
    double mean_share = 0, max_share = 0;
    for (uint32_t p = 0; p < 12; ++p) {
        const Vec3 w = solution.rotation_rad_s[p];
        const double share = std::abs(mesh::dot(w, sphere_centroid(plates[p]))) / mesh::length(w);
        mean_share += share / 12;
        max_share = std::max(max_share, share);
    }
    std::cout << "12 plates, N 16, seed 1: condition " << condition << ", spin share mean " << mean_share
              << " max " << max_share << ", mean speed " << solution.diagnostics.mean_speed_cm_yr
              << " cm/yr, iterations " << solution.diagnostics.nonlinear_iterations << "\n";
    RecordProperty("condition", std::to_string(condition));
    RecordProperty("spin_share_mean", std::to_string(mean_share));
    RecordProperty("spin_share_max", std::to_string(max_share));
    EXPECT_LT(condition, 1e12);
}

TEST(TectonicDynamics, SphereInlandPressureWalksAcrossFaceEdgesAndStopsAtOtherPlates) {
    using tectonics::dynamics::inland_pressure;
    using tectonics::grid::east;
    const auto grid = TectonicGrid::cube(8, kEarthRadiusKm);
    std::vector<uint32_t> owners(grid->cellcount(), 0);
    std::vector<double> pressure(grid->cellcount(), 2e12);
    EXPECT_NEAR(inland_pressure(*grid, owners, pressure, 0, east, 3e6), 2e12, 1e-3);
    EXPECT_NEAR(inland_pressure(*grid, owners, pressure, 0, east, 1e9), 2e12, 1e-3);  // once around
    // A cell at a face's east edge: the root lies beyond it, on the next face.
    tectonics::grid::CellId cell = 0;
    while (mesh::cubeposition(grid->position(grid->ahead(cell, east).cell)).face ==
           mesh::cubeposition(grid->position(cell)).face) ++cell;
    const auto beyond = grid->ahead(cell, east);
    pressure.assign(grid->cellcount(), 0);
    pressure[beyond.cell] = 1e13;
    const double one = inland_pressure(*grid, owners, pressure, cell, east, 3e6);
    EXPECT_GT(one, 0);
    // The walk goes on straight in the next face's frame.
    pressure[grid->ahead(beyond.cell, beyond.direction).cell] = 1e13;
    EXPECT_GT(inland_pressure(*grid, owners, pressure, cell, east, 3e6), one);
    EXPECT_DOUBLE_EQ(inland_pressure(*grid, owners, pressure, cell, east, 0), 0);
    owners[beyond.cell] = 1;
    EXPECT_DOUBLE_EQ(inland_pressure(*grid, owners, pressure, cell, east, 3e6), 0);
    EXPECT_THROW(inland_pressure(*grid, owners, pressure, cell, 4, 3e6), std::runtime_error);
    const auto torus = TectonicGrid::planar(8, 4, 800);
    EXPECT_THROW(inland_pressure(*torus, std::vector<uint32_t>(32), std::vector<double>(32), 0, east, 3e6),
                 std::runtime_error);
}
