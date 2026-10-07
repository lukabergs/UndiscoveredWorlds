#include "gtest/gtest.h"
#include "simulations/geology/plate_tectonics/tectonic_columns.hpp"
#include "simulations/geology/plate_tectonics/tectonic_thermal.hpp"
#include "simulations/geology/plate_tectonics/lithosphere.hpp"
#include "validation/tectonics/tectonic_snapshots.hpp"
#include <algorithm>
#include <cmath>
#include <filesystem>

using namespace tectonics::columns;

TEST(TectonicColumns, AiryRootsAndWaterLoadingHaveAnalyticalElevations) {
    for (double z : {-5000., -500., 0., 1500., 8000., 12000.}) {
        const auto c = seed(z, true, 100);
        EXPECT_NEAR(elevation_m(c), z, 1e-10);
        EXPECT_GT(c.continent, z);
    }
    const auto ridge = seed(-100, false, 0);
    EXPECT_DOUBLE_EQ(ridge.ocean, 7000);
    EXPECT_DOUBLE_EQ(elevation_m(ridge), -2500);
    EXPECT_NEAR(elevation_m(seed(-100, false, 80)), -2500 - tectonics::thermal::subsidence_m(80), 1e-10);
}

TEST(TectonicColumns, StationaryCoolingIsIndependentOfTimestepAndDoesNotRemoveCrust) {
    Model a, b;
    const std::vector<Cell> initial(32, seed(-5000, false, 20));
    a.initialize(8, 4, 100000, initial); b.initialize(8, 4, 100000, initial);
    const std::vector<Displacement> velocity(32);
    for (int i = 0; i < 20; ++i) a.advance(velocity, 1, 0, 1);
    for (int i = 0; i < 80; ++i) b.advance(velocity, .25, 0, 1);
    EXPECT_DOUBLE_EQ(a.cells()[0].ocean, 7000);
    EXPECT_DOUBLE_EQ(age_myr(a.cells()[0]), 40);
    EXPECT_NEAR(elevation_m(a.cells()[0]), elevation_m(b.cells()[0]), 1e-10);
    EXPECT_NEAR(a.budget().relative_mass_residual, 0, 1e-14);
}

TEST(TectonicColumns, ConvergenceThickensContinentsWithoutInventingOrDeletingMass) {
    Model world;
    world.initialize(32, 8, 100000, std::vector<Cell>(256, seed(0, true, 100)));
    std::vector<Displacement> v(256);
    for (size_t i = 0; i < v.size(); ++i) v[i].x = i % 32 < 16 ? .05 : -.05;
    for (int i = 0; i < 80; ++i) world.advance(v, .25, 0, 1);
    double peak = 0;
    for (const auto& c : world.cells()) peak = std::max(peak, elevation_m(c));
    EXPECT_GT(peak, 5000); // A root can support high mountains; no 4 km ceiling.
    EXPECT_LT(std::abs(world.budget().relative_mass_residual), 1e-12);
    EXPECT_DOUBLE_EQ(world.budget().subducted_ocean_kg, 0);
    EXPECT_DOUBLE_EQ(world.budget().generated_ocean_kg, 0);
}

TEST(TectonicColumns, OceanOpeningAndSubductionAreExplicitMassSourcesAndSinks) {
    Model world;
    world.initialize(32, 8, 100000, std::vector<Cell>(256, seed(-4000, false, 80)));
    std::vector<Displacement> v(256);
    for (size_t i = 0; i < v.size(); ++i) v[i].x = i % 32 < 16 ? .05 : -.05;
    for (int i = 0; i < 20; ++i) world.advance(v, 1, 0, 1);
    EXPECT_GT(world.budget().generated_ocean_kg, 0);
    EXPECT_GT(world.budget().subducted_ocean_kg, 0);
    EXPECT_LT(std::abs(world.budget().relative_mass_residual), 1e-12);
    for (const auto& c : world.cells()) EXPECT_NEAR(c.ocean, 7000, 1e-9);
}

TEST(TectonicColumns, ErodedRockReappearsAsSedimentWithEqualMass) {
    auto cells = std::vector<Cell>(128, seed(0, true, 100));
    for (size_t i = 0; i < cells.size(); ++i) if (i % 16 == 8) cells[i] = seed(6000, true, 100);
    Model world;
    world.initialize(16, 8, 100000, cells);
    for (int i = 0; i < 20; ++i) world.advance(std::vector<Displacement>(128), 1, 1, 1);
    EXPECT_GT(world.budget().eroded_kg, 0);
    EXPECT_NEAR(world.budget().deposited_kg / world.budget().eroded_kg, 1, 1e-12);
    EXPECT_LT(std::abs(world.budget().relative_mass_residual), 1e-12);
    EXPECT_LT(elevation_m(world.cells()[8]), 6000);
    EXPECT_GT(world.cells()[7].sediment, 0);
}

TEST(TectonicColumns, FractionalTransportWrapsAndPreservesVolumeAndAge) {
    Model uninitialized;
    EXPECT_THROW(uninitialized.advance({}, 1, 0, 1), std::runtime_error);
    auto cells = std::vector<Cell>(32, seed(0, true, 100));
    cells[7] = seed(0, true, 250); // Uniform elevation isolates advection from root flow.
    Model world;
    world.initialize(8, 4, 100000, cells);
    world.advance(std::vector<Displacement>(32, {.25, 0}), 1, 0, 0);
    EXPECT_NEAR(world.cells()[0].continent, .75 * cells[0].continent + .25 * cells[7].continent, 1e-9);
    EXPECT_LT(std::abs(world.budget().relative_mass_residual), 1e-12);
    EXPECT_DOUBLE_EQ(age_myr(world.cells()[0]), .75 * 101 + .25 * 251);
    EXPECT_THROW(world.advance(std::vector<Displacement>(32, {1.1, 0}), 1, 0, 1), std::runtime_error);
}

TEST(TectonicColumns, RootsSpreadWithoutErosionAndPreserveMass) {
    auto cells = std::vector<Cell>(128, seed(0, true, 100));
    for (size_t i = 0; i < cells.size(); ++i) if (i % 16 == 8) cells[i] = seed(12000, true, 100);
    Model world;
    world.initialize(16, 8, 100000, cells);
    for (int i = 0; i < 20; ++i) world.advance(std::vector<Displacement>(128), 1, 0, 1);
    EXPECT_LT(elevation_m(world.cells()[8]), 6000);
    EXPECT_GT(elevation_m(world.cells()[7]), 0);
    EXPECT_GT(world.budget().redistributed_continent_kg, 0);
    EXPECT_DOUBLE_EQ(world.budget().eroded_kg, 0);
    EXPECT_LT(std::abs(world.budget().relative_mass_residual), 1e-12);
}

TEST(TectonicColumns, GentleInteriorReliefSurvivesRootFlowOver200Myr) {
    constexpr uint32_t w = 32, h = 8;
    constexpr double pi = 3.141592653589793;
    std::vector<Cell> initial(w * h);
    for (size_t i = 0; i < initial.size(); ++i)
        initial[i] = seed(700 + 300 * std::sin(2 * pi * (i % w) / 8), true, 100);
    Model world; world.initialize(w, h, 156250, initial);
    for (int i = 0; i < 200; ++i) world.advance(std::vector<Displacement>(w * h), 1, 0, 1);
    double amplitude = 0;
    for (size_t i = 0; i < initial.size(); ++i)
        amplitude += (elevation_m(world.cells()[i]) - 700) * std::sin(2 * pi * (i % w) / 8);
    const double retained = amplitude / (w * h * 150);
    RecordProperty("root_flow_relief_retained", retained);
    EXPECT_GT(retained, .95);
    EXPECT_LT(std::abs(world.budget().relative_mass_residual), 1e-12);
}

TEST(TectonicColumns, VanishingContinentalTracesHaveZeroRootConductance) {
    auto tiny = seed(0, true, 100);
    tiny.continent = 1e-80;
    tiny.continental_age_volume = tiny.continent * 100;
    Model world; world.initialize(8, 4, 156250, std::vector<Cell>(32, tiny));
    for (int i = 0; i < 5; ++i) world.advance(std::vector<Displacement>(32), 1, 0, 0);
    for (const auto& c : world.cells()) {
        EXPECT_DOUBLE_EQ(c.continent, tiny.continent);
        EXPECT_NEAR(age_myr(c), 105, 1e-12);
    }
    EXPECT_DOUBLE_EQ(world.budget().redistributed_continent_kg, 0);
}

TEST(TectonicColumns, RootFlowAndErosionConvergeWhenTimestepIsReduced) {
    auto cells = std::vector<Cell>(128, seed(0, true, 100));
    for (size_t i = 0; i < cells.size(); ++i) if (i % 16 == 8) cells[i] = seed(6000, true, 100);
    const auto run = [&](double dt) {
        Model world; world.initialize(16, 8, 100000, cells);
        for (int i = 0; i < static_cast<int>(10 / dt); ++i)
            world.advance(std::vector<Displacement>(128), dt, 1, 1);
        EXPECT_LT(std::abs(world.budget().relative_mass_residual), 1e-12);
        return world.cells();
    };
    // Stay below the adaptive root-flow stability step; halving a large outer
    // step need not halve its internally selected effective timestep.
    const auto coarse = run(.125), medium = run(.0625), fine = run(.015625);
    double coarse_error = 0, medium_error = 0;
    for (size_t i = 0; i < cells.size(); ++i) {
        coarse_error += std::abs(elevation_m(coarse[i]) - elevation_m(fine[i]));
        medium_error += std::abs(elevation_m(medium[i]) - elevation_m(fine[i]));
    }
    EXPECT_LT(medium_error, coarse_error * .8);
}

TEST(TectonicColumns, InitializationKeepsLowContinentalLandAboveWater) {
    tectonics::scenario::Scenario s;
    s.seed = 1; s.width = 32; s.height = 16; s.plate_count = 4;
    s.evolution_duration_myr = 1;  
    lithosphere world(s);
    const auto datum = world.getSeaLevelMeters();
    std::vector<uint16_t> lowland(512, datum + 200);
    world.importMetricHeightMap(lowland.data(), datum);
    for (const auto& c : world.getColumnModel().cells()) {
        EXPECT_DOUBLE_EQ(c.continental_area, 1);
        EXPECT_NEAR(elevation_m(c), 200, .1);
    }
}

TEST(TectonicColumns, Integrated200MyrReplayClosesMassBudgetAndPreservesContinents) {
    tectonics::scenario::Scenario s;
    s.seed = 1; s.width = 64; s.height = 32; s.plate_count = 8;
    s.evolution_duration_myr = 200;  
    lithosphere a(s), b(s);
    double initial_continent = 0;
    double initial_max_age = 0;
    for (const auto& c : a.getColumnModel().cells()) {
        initial_continent += c.continent;
        initial_max_age = std::max(initial_max_age, age_myr(c));
    }
    for (int i = 0; i < 200; ++i) {
        a.update(); b.update();
        ASSERT_LT(std::abs(a.getColumnModel().budget().relative_mass_residual), 1e-10) << i;
        for (const auto& c : a.getColumnModel().cells())
            ASSERT_LE(age_myr(c), initial_max_age + i + 1 + 1e-6) << "update " << i;
    }
    EXPECT_TRUE(a.isFinished());
    EXPECT_DOUBLE_EQ(a.getTimeMyr(), 200);
    EXPECT_GT(initial_continent, 0);
    for (size_t i = 0; i < a.getColumnModel().cells().size(); ++i) {
        EXPECT_DOUBLE_EQ(a.getColumnModel().cells()[i].continent, b.getColumnModel().cells()[i].continent);
        EXPECT_FLOAT_EQ(a.getTopography()[i], b.getTopography()[i]);
    }
}

TEST(TectonicColumns, ResolvedRiftsCreateIndependentPlatesAndReplayWithoutLosingMaterial) {
    // Stronger forcing is a controlled failure fixture, not the benchmark recipe.
    tectonics::scenario::Scenario s;
    s.seed=1; s.width=64; s.height=32; s.plate_count=4;
    s.evolution_duration_myr=100; s.delta_time_myr=.5;
      s.mantle_flow_speed_cm_yr=20;
    lithosphere a(s), b(s);
    for (int step=0;step<200;++step) {
        a.update(); b.update();
        ASSERT_LT(std::abs(a.getColumnModel().budget().relative_mass_residual),1e-10);
        ASSERT_EQ(a.getTopologyDiagnostics().stable_ids,b.getTopologyDiagnostics().stable_ids);
    }
    const auto& events=a.getTopologyDiagnostics().events;
    ASSERT_GE(events.size(),1U);
    EXPECT_EQ(events[0].kind,tectonics::topology::EventKind::Split);
    EXPECT_GT(a.getTopologyDiagnostics().stable_ids.size(),s.plate_count);
    EXPECT_EQ(events.size(),b.getTopologyDiagnostics().events.size());
    for (size_t i=0;i<events.size();++i) {
        EXPECT_DOUBLE_EQ(events[i].time_myr,b.getTopologyDiagnostics().events[i].time_myr);
        EXPECT_EQ(events[i].result,b.getTopologyDiagnostics().events[i].result);
    }
    for (uint32_t i=0;i<s.width*s.height;++i) {
        EXPECT_FLOAT_EQ(a.getTopography()[i],b.getTopography()[i]);
        EXPECT_EQ(a.getPlatesMap()[i],b.getPlatesMap()[i]);
        EXPECT_DOUBLE_EQ(a.getColumnModel().cells()[i].continental_age_volume,
                         b.getColumnModel().cells()[i].continental_age_volume);
    }
    RecordProperty("split_events",static_cast<int>(events.size()));
}

TEST(TectonicColumns, ExcessMountainForceMatchesAiryColumnUnitsAndVanishesOnUnthickenedCrust) {
    EXPECT_DOUBLE_EQ(mountain_force_n_m(seed(0,true,100)),0);
    EXPECT_DOUBLE_EQ(mountain_force_n_m(seed(-4000,false,100)),0);
    const auto mountain=seed(5000,true,100);
    const double expected=.5*9.80665*2800*(1-2800./3300)*(68000.*68000-35000.*35000);
    EXPECT_NEAR(mountain_force_n_m(mountain),expected,expected*1e-14);
    EXPECT_GT(mountain_force_n_m(seed(10000,true,100)),mountain_force_n_m(mountain));
}

// Every field of two snapshots is equal.
static bool same_fields(const tectonics::contract::Snapshot& a, const tectonics::contract::Snapshot& b) {
    return a.continental_crust_m == b.continental_crust_m && a.oceanic_crust_m == b.oceanic_crust_m &&
        a.sediment_m == b.sediment_m && a.rift_damage == b.rift_damage && a.plume_influence == b.plume_influence &&
        a.heightmap == b.heightmap && a.plate_id == b.plate_id && a.crust_age_steps == b.crust_age_steps &&
        a.crust_age_myr == b.crust_age_myr && a.crust_thickness == b.crust_thickness && a.crust_class == b.crust_class &&
        a.uplift_tendency == b.uplift_tendency && a.subsidence_tendency == b.subsidence_tendency &&
        a.accumulated_strain == b.accumulated_strain && a.boundary_type == b.boundary_type &&
        a.boundary_distance == b.boundary_distance && a.boundary_segment_id == b.boundary_segment_id &&
        a.nearest_boundary_id == b.nearest_boundary_id && a.deforming_region_id == b.deforming_region_id &&
        a.deforming_region_type == b.deforming_region_type && a.deformation_rate == b.deformation_rate &&
        a.deformation_velocity_x == b.deformation_velocity_x && a.deformation_velocity_y == b.deformation_velocity_y &&
        a.convergence_score == b.convergence_score && a.divergence_score == b.divergence_score &&
        a.shear_score == b.shear_score && a.geologic_regime == b.geologic_regime;
}

TEST(TectonicColumns, PhysicalLayersSurviveBundleReadAndAffectComparison) {
    tectonics::scenario::Scenario s;
    s.seed = 1; s.width = 32; s.height = 16; s.plate_count = 4;
    s.evolution_duration_myr = 1;  
    lithosphere world(s); world.update();
    auto snapshot = tectonics::contract::capture_snapshot(world);
    snapshot.column_budget.noise_redistributed_kg=123456;
    snapshot.column_budget.redistributed_sediment_kg=654321;
    snapshot.topology_diagnostics.events = {{tectonics::topology::EventKind::Split,1,0,0,4,128},
        {tectonics::topology::EventKind::Fusion,2,1,2,5,200}};
    const auto folder = std::filesystem::temp_directory_path() / "uw-isostatic-column-test";
    tectonics::tooling::write_snapshot_bundle(folder, snapshot, {});
    auto loaded = tectonics::tooling::load_snapshot_bundle(folder).snapshot;
    EXPECT_EQ(loaded.continental_crust_m, snapshot.continental_crust_m);
    EXPECT_EQ(loaded.oceanic_crust_m, snapshot.oceanic_crust_m);
    EXPECT_EQ(loaded.sediment_m, snapshot.sediment_m);
    EXPECT_DOUBLE_EQ(loaded.column_budget.mass_kg, snapshot.column_budget.mass_kg);
    EXPECT_DOUBLE_EQ(loaded.column_budget.noise_redistributed_kg,123456);
    EXPECT_DOUBLE_EQ(loaded.column_budget.redistributed_sediment_kg,654321);
    EXPECT_TRUE(loaded.topology_diagnostics.enabled);
    EXPECT_EQ(loaded.topology_diagnostics.reassigned_cells, snapshot.topology_diagnostics.reassigned_cells);
    EXPECT_EQ(loaded.topology_diagnostics.max_detached_components, snapshot.topology_diagnostics.max_detached_components);
    EXPECT_EQ(loaded.topology_diagnostics.stable_ids, snapshot.topology_diagnostics.stable_ids);
    ASSERT_EQ(loaded.topology_diagnostics.events.size(),2U);
    EXPECT_EQ(loaded.topology_diagnostics.events[0].kind,tectonics::topology::EventKind::Split);
    EXPECT_EQ(loaded.topology_diagnostics.events[1].kind,tectonics::topology::EventKind::Fusion);
    EXPECT_EQ(loaded.topology_diagnostics.events[1].result,5U);
    EXPECT_DOUBLE_EQ(loaded.topology_diagnostics.events[1].time_myr,2);
    EXPECT_TRUE(same_fields(snapshot, loaded));
    loaded.continental_crust_m[0] += 10;
    EXPECT_FALSE(same_fields(snapshot, loaded));
}

TEST(TectonicColumns, SedimentLoadContributesToLayeredGpeAndConservativeRootFlow) {
    for (double z:{-1000.,0.,1000.,5000.}) {
        const auto c=seed(z,true,40);
        EXPECT_NEAR(mountain_force_n_m(c,true),mountain_force_n_m(c),1);
    }
    auto column=seed(0,true,40);column.sediment=27500;
    EXPECT_NEAR(elevation_m(column,0),5000,1e-9);
    EXPECT_DOUBLE_EQ(mountain_force_n_m(column),0);
    EXPECT_NEAR(mountain_force_n_m(column,true),.5*9.80665*1.15875e12,1);
    constexpr uint32_t w=32,h=8;
    std::vector<Cell> cells(w*h,seed(0,true,40));
    for (size_t i=0;i<cells.size();++i)
        cells[i].sediment=50000*std::pow(.5+.5*std::cos(6.283185307179586*(i%w)/w),4);
    Model loaded,old;loaded.initialize(w,h,50000,cells);old.initialize(w,h,50000,cells);
    const std::vector<Displacement> stopped(cells.size());
    for (int i=0;i<10;++i) {
        loaded.advance(stopped,1,0,0,0,nullptr,{}, {},true);
        old.advance(stopped,1,0,0,0);
    }
    double initial_s=0,final_s=0,initial_c=0,final_c=0,initial_age=0,final_age=0,peak=-1e9,old_peak=-1e9;
    for (size_t i=0;i<cells.size();++i) {
        const auto& c=loaded.cells()[i];
        initial_s+=cells[i].sediment;final_s+=c.sediment;
        initial_c+=cells[i].continent;final_c+=c.continent;
        initial_age+=cells[i].continental_age_volume;final_age+=c.continental_age_volume;
        EXPECT_GE(c.sediment,0);EXPECT_GE(c.continent,0);
        peak=std::max(peak,elevation_m(c,0));old_peak=std::max(old_peak,elevation_m(old.cells()[i],0));
    }
    EXPECT_LT(peak,old_peak);
    EXPECT_NEAR(initial_s,final_s,initial_s*1e-12);
    EXPECT_NEAR(initial_c,final_c,initial_c*1e-12);
    EXPECT_NEAR(initial_age+10*initial_c,final_age,final_age*1e-12);
    EXPECT_GT(loaded.budget().redistributed_sediment_kg,0);
    EXPECT_DOUBLE_EQ(loaded.budget().eroded_kg,0);
    EXPECT_LT(std::abs(loaded.budget().relative_mass_residual),1e-12);
}

TEST(TectonicColumns, SedimentLoadedFlowConvergesWithSmallerTimesteps) {
    std::vector<Cell> cells(128,seed(0,true,40));
    for (size_t i=0;i<cells.size();++i) cells[i].sediment=i%16==8?25000:1000;
    const auto run=[&](double dt) {
        Model model;model.initialize(16,8,100000,cells);
        for (int i=0;i<static_cast<int>(10/dt);++i)
            model.advance(std::vector<Displacement>(128),dt,0,0,0,nullptr,{}, {},true);
        EXPECT_LT(std::abs(model.budget().relative_mass_residual),1e-12);
        return model.cells();
    };
    const auto coarse=run(.125),medium=run(.0625),fine=run(.015625);
    double coarse_error=0,medium_error=0;
    for (size_t i=0;i<cells.size();++i) {
        coarse_error+=std::abs(elevation_m(coarse[i],0)-elevation_m(fine[i],0));
        medium_error+=std::abs(elevation_m(medium[i],0)-elevation_m(fine[i],0));
    }
    EXPECT_LT(medium_error,.8*coarse_error);
}

// The rock uplift rate (erosion plan D1) counts the tectonic sub-steps only: the
// diffusion and the rivers change the columns, never the rate; a plume's basalt on a
// continental column raises it by 400/3300 of its thickness (Airy, land), so 330 m
// over 1 Myr is a step rate of 40 m/Myr, kept as an exponential mean with tau = 5 Myr.
TEST(TectonicColumns, RockUpliftRateCountsTectonicChangeOnly) {
    tectonics::margins::Parameters still; // No root flow: the roots would carry the change sideways.
    still.root_mobility = 0;
    {
        auto cells = std::vector<Cell>(128, seed(500, true, 100));
        for (size_t i = 0; i < cells.size(); ++i) if (i % 16 == 8) cells[i] = seed(3000, true, 100);
        auto rivers = still;
        rivers.enabled = true; // routeSediment with the default river erosion; no forcing, no thermal change
        Model world;
        world.initialize(16, 8, 100000, cells);
        for (int i = 0; i < 3; ++i) world.advance(std::vector<Displacement>(128), 1, 1, 1, 1, nullptr, rivers);
        EXPECT_GT(world.budget().eroded_kg, 0);
        EXPECT_GT(world.budget().river_eroded_kg, 0);
        EXPECT_LT(elevation_m(world.cells()[8]), 3000);
        ASSERT_EQ(world.upliftRate().size(), 128u);
        for (double rate : world.upliftRate()) EXPECT_NEAR(rate, 0, 1e-12);
    }
    tectonics::plumes::Source plume;
    plume.basalt_m.assign(32, 0);
    plume.influence.assign(32, 0);
    plume.basalt_m[9] = 330;
    Model world;
    world.initialize(8, 4, 100000, std::vector<Cell>(32, seed(500, true, 100)));
    world.advance(std::vector<Displacement>(32), 1, 0, 1, 1, nullptr, still, {}, false, &plume, 5);
    const double weight = -std::expm1(-1.0 / 5); // 1 - e^(-dt/tau)
    EXPECT_NEAR(world.upliftRate()[9], 40 * weight, 1e-12);
    EXPECT_NEAR(world.upliftRate()[9], 7.2507698768807245, 1e-12);
    for (size_t i = 0; i < 32; ++i) if (i != 9) EXPECT_NEAR(world.upliftRate()[i], 0, 1e-12);
    // A step without tectonic change decays the mean by e^(-dt/tau).
    world.advance(std::vector<Displacement>(32), 1, 0, 1, 1, nullptr, still, {}, false, nullptr, 5);
    EXPECT_NEAR(world.upliftRate()[9], 40 * weight * std::exp(-1.0 / 5), 1e-12);
    EXPECT_THROW(world.advance(std::vector<Displacement>(32), 1, 0, 1, 1, nullptr, still, {}, false, nullptr, 0),
        std::runtime_error);
}

// The measurement is Eulerian: plate motion over a slope changes the elevation of a
// fixed cell by -v . grad z. Land columns 500 + 100 x m (x = 0..15, 100 km cells, so
// grad z = 1 m/km) move east by v = 0.25 cell (25 km) in 1 Myr with nothing else
// acting. The remap's parabolic reconstruction is exact on linear data, so a cell
// whose stencil (x - 3 .. x + 2) stays off the wrap's jump gets the mean of the
// ramp shifted by 0.25 cell: the step rate is -v . grad z = -25 km/Myr x 1 m/km =
// -25 m/Myr and r = -25 (1 - e^(-1/5)) = -4.5317311730504545. A start point taken
// after transportColumns would give 0.
TEST(TectonicColumns, RockUpliftRateIsEulerianAdvectionOverASlope) {
    tectonics::margins::Parameters still;
    still.root_mobility = 0;
    std::vector<Cell> cells(128);
    for (size_t i = 0; i < cells.size(); ++i) cells[i] = seed(500 + 100. * (i % 16), true, 100);
    Model world;
    world.initialize(16, 8, 100000, cells);
    world.advance(std::vector<Displacement>(128, {.25, 0}), 1, 0, 1, 1, nullptr, still, {}, false, nullptr, 5);
    for (size_t i = 0; i < 128; ++i) {
        if (i % 16 < 3 || i % 16 > 13) continue; // stencils that see the 1500 m jump across the wrap
        EXPECT_NEAR(world.upliftRate()[i], -25 * -std::expm1(-1.0 / 5), 1e-9) << "cell " << i;
        EXPECT_NEAR(world.upliftRate()[i], -4.5317311730504545, 1e-9) << "cell " << i;
    }
}

// Thermal subsidence of an ageing ocean floor is tectonic: a 20 Myr old ocean column
// (thermal strength 1, nothing else acting) ages to 21 Myr in one step and sinks by
// subsidence_m(21) - subsidence_m(20) = 39.50061759809455 m (Parsons & Sclater law
// of tectonic_thermal.hpp), so r = -39.50061759809455 (1 - e^(-1/5)) = -7.160247204961215.
TEST(TectonicColumns, RockUpliftRateCountsThermalSubsidence) {
    Model world;
    world.initialize(8, 4, 100000, std::vector<Cell>(32, seed(-5000, false, 20)));
    world.advance(std::vector<Displacement>(32), 1, 0, 1, 1, nullptr, {}, {}, false, nullptr, 5);
    const double sinking = tectonics::thermal::subsidence_m(21) - tectonics::thermal::subsidence_m(20);
    EXPECT_NEAR(sinking, 39.50061759809455, 1e-9);
    for (double rate : world.upliftRate()) {
        EXPECT_NEAR(rate, -sinking * -std::expm1(-1.0 / 5), 1e-9);
        EXPECT_NEAR(rate, -7.160247204961215, 1e-9);
    }
}
