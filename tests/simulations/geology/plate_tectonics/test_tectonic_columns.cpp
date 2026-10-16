#include "gtest/gtest.h"
#include "simulations/geology/plate_tectonics/tectonic_columns.hpp"
#include "simulations/geology/plate_tectonics/tectonic_thermal.hpp"
#include "simulations/geology/plate_tectonics/lithosphere.hpp"
#include "validation/tectonics/tectonic_snapshots.hpp"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <limits>
#include <utility>

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

TEST(TectonicColumns, ContinentalFractionSeedsShelvesToTheRequestedArea) {
    tectonics::scenario::Scenario s;
    s.seed = 1; s.width = 128; s.height = 64; s.plate_count = 8; s.sea_level = .65f;
    s.evolution_duration_myr = 1; s.boundary_conditioned_terrain = true;
    const auto share = [](lithosphere& world) {
        double continental = 0;
        const auto& cells = world.getColumnModel().cells();
        for (const auto& c : cells) continental += c.continental_area >= .5 ? 1 : 0;
        return continental / cells.size();
    };
    lithosphere legacy(s);
    s.continental_fraction = .40;
    lithosphere shelved(s);
    EXPECT_GT(share(legacy), .5);
    EXPECT_NEAR(share(shelved), .40, .01);
    // The land share is unchanged: the shelf lies below sea level.
    double legacy_land = 0, shelved_land = 0;
    for (const auto& c : legacy.getColumnModel().cells()) legacy_land += elevation_m(c) > 0;
    for (const auto& c : shelved.getColumnModel().cells()) shelved_land += elevation_m(c) > 0;
    EXPECT_NEAR(shelved_land, legacy_land, .01 * legacy.getColumnModel().cells().size());
    // With the baseline's 2000 m margin taper the slope continues to -3500 m, so the
    // continental area summed over the columns, not only the cells over one half, is f.
    s.initial_margin_width_m = 2000;
    lithosphere tapered(s);
    double area = 0;
    for (const auto& c : tapered.getColumnModel().cells()) area += c.continental_area;
    EXPECT_NEAR(area / tapered.getColumnModel().cells().size(), .40, .01);
    s.initial_margin_width_m = 0;
    s.continental_fraction = .35; // 1 - sea_level: no shelf, a 2.5 km cliff at the coast
    EXPECT_THROW(lithosphere{s}, std::runtime_error);
    s.continental_fraction = .30;
    EXPECT_THROW(lithosphere{s}, std::runtime_error);
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

// The measurement is in the plate's frame (erosion plan H2): plate motion over a slope
// is not uplift. Land columns 500 + 100 x m (x = 0..15, 100 km cells, so grad z =
// 1 m/km) move east by v = 0.25 cell (25 km) in 1 Myr with nothing else acting. The
// elevation before the step is carried by the same remap as the columns, and land
// elevation is linear in them, so the rate is 0 up to rounding in every cell, also
// where the remap sees the 1500 m jump across the wrap (tolerance 1e-9 m/Myr). A
// fixed-cell (Eulerian) rate would be -v . grad z = -25 m/Myr, r = -4.53.
TEST(TectonicColumns, RockUpliftRateIsZeroForPlateMotionOverASlope) {
    tectonics::margins::Parameters still;
    still.root_mobility = 0;
    std::vector<Cell> cells(128);
    for (size_t i = 0; i < cells.size(); ++i) cells[i] = seed(500 + 100. * (i % 16), true, 100);
    Model world;
    world.initialize(16, 8, 100000, cells);
    world.advance(std::vector<Displacement>(128, {.25, 0}), 1, 0, 1, 1, nullptr, still, {}, false, nullptr, 5);
    for (size_t i = 0; i < 128; ++i) EXPECT_NEAR(world.upliftRate()[i], 0, 1e-9) << "cell " << i;
}

// A tectonic change on a moving plate is measured as the change alone. On the same
// ramp moving east by 0.25 cell, a plume adds 330 m of basalt under cell 9 of every
// row: 40 m on a land column (as in RockUpliftRateCountsTectonicChangeOnly), so the
// step rate there is 40 m/Myr and r = 40 (1 - e^(-1/5)) = 7.2507698768807245, not
// the Eulerian 40 - 25 = 15 m/Myr; every other cell stays at 0.
TEST(TectonicColumns, RockUpliftRateOnAMovingPlateIsTheTectonicChange) {
    tectonics::margins::Parameters still;
    still.root_mobility = 0;
    std::vector<Cell> cells(128);
    for (size_t i = 0; i < cells.size(); ++i) cells[i] = seed(500 + 100. * (i % 16), true, 100);
    tectonics::plumes::Source plume;
    plume.basalt_m.assign(128, 0);
    plume.influence.assign(128, 0);
    for (size_t row = 0; row < 8; ++row) plume.basalt_m[row * 16 + 9] = 330;
    Model world;
    world.initialize(16, 8, 100000, cells);
    world.advance(std::vector<Displacement>(128, {.25, 0}), 1, 0, 1, 1, nullptr, still, {}, false, &plume, 5);
    for (size_t i = 0; i < 128; ++i) {
        const double expected = i % 16 == 9 ? 40 * -std::expm1(-1.0 / 5) : 0;
        EXPECT_NEAR(world.upliftRate()[i], expected, 1e-9) << "cell " << i;
    }
    EXPECT_NEAR(world.upliftRate()[9], 7.2507698768807245, 1e-9);
}

// Pure translation of a curved surface (no other process): land columns
// z = 1000 + 500 sin(2 pi x / 16) + 300 cos(2 pi y / 8) on a 16 x 8 torus move by
// (0.16, 0.1) cell in 1 Myr, at subdivision 1 and at the baseline's 4. The carried
// elevation moves with the engine's own remap, so the mean is 0 up to rounding
// (largest 9.1e-13 m/Myr at subdivision 1, 1.0e-12 at 4; tolerance 1e-9). A bilinear
// sample at x - u (the first H2 method) leaves about a (1 - a) z'' / 2 per axis: at
// subdivision 1 a largest step rate of 12.9 m/Myr, a mean of 2.34 m/Myr after one
// step (w = 0.181). At subdivision 4 the subcells start as uniform 4 x 4 blocks,
// which a shift within one sub-step moves exactly as bilinear interpolation does.
TEST(TectonicColumns, RockUpliftRateIsZeroForPlateMotionOverACurvedSurface) {
    tectonics::margins::Parameters still;
    still.root_mobility = 0;
    const double pi = std::acos(-1.0);
    const auto z = [&](double x, double y) {
        return 1000 + 500 * std::sin(2 * pi * x / 16) + 300 * std::cos(2 * pi * y / 8);
    };
    std::vector<Cell> cells(128);
    for (size_t i = 0; i < cells.size(); ++i)
        cells[i] = seed(z(static_cast<double>(i % 16), static_cast<double>(i / 16)), true, 100);
    for (uint32_t subdivision : {1u, 4u}) {
        Model world;
        world.initialize(16, 8, 100000, cells, subdivision);
        const Displacement u{.16, .1};
        world.advance(std::vector<Displacement>(128, u), 1, 0, 1, 1, nullptr, still, {}, false, nullptr, 5);
        const double weight = -std::expm1(-1.0 / 5);
        double residual = 0, bilinear = 0;
        for (size_t i = 0; i < 128; ++i) {
            residual = std::max(residual, std::abs(world.upliftRate()[i]));
            const double sx = i % 16 - u.x, sy = i / 16 - u.y;
            const double fx = std::floor(sx), fy = std::floor(sy), ax = sx - fx, ay = sy - fy;
            const double sampled = (1 - ay) * ((1 - ax) * z(fx, fy) + ax * z(fx + 1, fy))
                + ay * ((1 - ax) * z(fx, fy + 1) + ax * z(fx + 1, fy + 1));
            bilinear = std::max(bilinear, std::abs((elevation_m(world.cells()[i]) - sampled) * weight));
        }
        EXPECT_LT(residual, 1e-9) << "subdivision " << subdivision;
        if (subdivision == 1) EXPECT_GT(bilinear, 2); // measured 2.34 m/Myr (a step rate of 12.9 m/Myr)
    }
}

// The running mean follows the rock, across the x and y wraps. Uniform land 500 m
// moves by (0.25, 0.25) cell per 1 Myr step; a plume adds 330 m of basalt (40 m) under
// cell (15, 7) in the first step only: r1 = 40 w there (w = 1 - e^(-1/5)), 0
// elsewhere. The second step has no tectonic change in the rock's frame, so its mean
// is (1 - w) r1 = e^(-1/5) r1 spread where the rock went. The limiter flattens the
// single-cell peak (subdivision 1) and the 4 x 4 block of subcells (subdivision 4,
// two sub-steps of half a subcell), so both move as constants: (15, 7) keeps
// 0.75 x 0.75 = 0.5625 of it, (0, 7) east across the x wrap and (15, 0) south across
// the y wrap get 0.1875, (0, 0) gets 0.0625. A mean read at the destination instead
// of the source, a wrong y sign or a missing wrap puts these shares elsewhere.
TEST(TectonicColumns, RockUpliftRateMeanFollowsTheRockAcrossTheWrap) {
    tectonics::margins::Parameters still;
    still.root_mobility = 0;
    tectonics::plumes::Source plume;
    plume.basalt_m.assign(128, 0);
    plume.influence.assign(128, 0);
    plume.basalt_m[7 * 16 + 15] = 330;
    for (uint32_t subdivision : {1u, 4u}) {
        Model world;
        world.initialize(16, 8, 100000, std::vector<Cell>(128, seed(500, true, 100)), subdivision);
        const std::vector<Displacement> u(128, {.25, .25});
        world.advance(u, 1, 0, 1, 1, nullptr, still, {}, false, &plume, 5);
        const double r1 = 40 * -std::expm1(-1.0 / 5);
        for (size_t i = 0; i < 128; ++i)
            EXPECT_NEAR(world.upliftRate()[i], i == 127 ? r1 : 0, 1e-9) << "subdivision " << subdivision << " cell " << i;
        world.advance(u, 1, 0, 1, 1, nullptr, still, {}, false, nullptr, 5);
        const double carried = r1 * std::exp(-1.0 / 5);
        EXPECT_NEAR(carried, 5.9364282816937015, 1e-12);
        for (size_t i = 0; i < 128; ++i) {
            double share = 0;
            if (i == 7 * 16 + 15) share = .5625;
            if (i == 7 * 16 + 0) share = .1875;
            if (i == 0 * 16 + 15) share = .1875;
            if (i == 0 * 16 + 0) share = .0625;
            EXPECT_NEAR(world.upliftRate()[i], share * carried, 1e-9) << "subdivision " << subdivision << " cell " << i;
        }
    }
}

// Diverging plates at subdivision 4: land z = 1000 + 500 sin(2 pi x / 16) (range
// R = 1000 m), the west half (x < 8) moves west and the east half east by 0.2 cell per
// 1 Myr step, so plates part between x = 7 and 8 and across the wrap; the opening
// stretches the boundary columns (real subsidence). With a window of 1e-3 Myr the
// mean is the step rate, so the carried source elevation is z_after - rate dt. Over
// six steps it stays within [low - 2R, high + 2R] of the elevations before each step
// (measured: within [low - 1353, high - 7.6] m). At the boundary columns (0, 7, 8, 15)
// it differs from the cell's own elevation before the step (the rock that stayed)
// by 4.9 m in the first step and by up to 1353 m later, once the boundary columns
// (thinned by a fifth each step) lie km below the sea level and their subcells
// differ strongly (the subcell elevations are not linear in the columns there): the
// measured rate at x = 7 then falls short of the cell's own elevation change by 23,
// 46 and 70 percent in the second to fourth steps (largest boundary rate 2327 m/Myr). Off the boundary (x = 2..5, 10..13) the rate stays at 0 (measured
// 4.0e-11 m/Myr), the translation of one varying column field.
TEST(TectonicColumns, RockUpliftRateAtDivergentBoundariesIsBounded) {
    tectonics::margins::Parameters still;
    still.root_mobility = 0;
    const double pi = std::acos(-1.0);
    std::vector<Cell> cells(128);
    for (size_t i = 0; i < cells.size(); ++i)
        cells[i] = seed(1000 + 500 * std::sin(2 * pi * static_cast<double>(i % 16) / 16), true, 100);
    std::vector<Displacement> u(128);
    for (size_t i = 0; i < u.size(); ++i) u[i].x = i % 16 < 8 ? -.2 : .2;
    Model world;
    world.initialize(16, 8, 100000, cells, 4);
    double boundary = 0, interior = 0, lowest = 1e300, highest = -1e300, boundary_rate = 0;
    for (int step = 0; step < 6; ++step) {
        std::vector<double> before(128);
        for (size_t i = 0; i < 128; ++i) before[i] = elevation_m(world.cells()[i]);
        const auto [low, high] = std::minmax_element(before.begin(), before.end());
        const double range = *high - *low;
        world.advance(u, 1, 0, 1, 1, nullptr, still, {}, false, nullptr, 1e-3);
        for (size_t i = 0; i < 128; ++i) {
            const double rate = world.upliftRate()[i];
            const double source = elevation_m(world.cells()[i]) - rate;
            EXPECT_GE(source, *low - 2 * range) << "step " << step << " cell " << i;
            EXPECT_LE(source, *high + 2 * range) << "step " << step << " cell " << i;
            lowest = std::min(lowest, source - *low);
            highest = std::max(highest, source - *high);
            const size_t x = i % 16;
            if (x == 0 || x == 7 || x == 8 || x == 15) {
                boundary = std::max(boundary, std::abs(source - before[i]));
                boundary_rate = std::max(boundary_rate, std::abs(rate));
            } else if ((x >= 2 && x <= 5) || (x >= 10 && x <= 13)) interior = std::max(interior, std::abs(rate));
        }
    }
    EXPECT_LT(interior, 1e-9);
    EXPECT_LT(boundary, 1400); // measured 1353 m
    EXPECT_GT(boundary_rate, 0);
    EXPECT_GT(lowest, -2000); EXPECT_LT(highest, 0);
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

// Cumulative shortening (multiscale plan S5.1): the export of a strain tensor. The
// shortening is minus the least principal component times the cell width; the
// direction is that component's axis in [0, pi) from x toward y. Simple shear xy = g/2
// has components +-g/2 and shortens along (1, -1) (3 pi / 4) for g > 0, along (1, 1)
// for g < 0; extension alone shortens nothing; equal components have direction 0.
TEST(TectonicColumns, ShorteningAndDirectionOfExplicitStrains) {
    const double pi = std::acos(-1.0);
    EXPECT_NEAR(shortening_m({-.1, 0, 0}, 1000), 100, 1e-12);
    EXPECT_EQ(shortening_direction({-.1, 0, 0}), 0);
    EXPECT_NEAR(shortening_m({0, -.1, 0}, 1000), 100, 1e-12);
    EXPECT_NEAR(shortening_direction({0, -.1, 0}), pi / 2, 1e-15);
    EXPECT_NEAR(shortening_m({0, 0, .05}, 1000), 50, 1e-12);
    EXPECT_NEAR(shortening_direction({0, 0, .05}), 3 * pi / 4, 1e-15);
    EXPECT_NEAR(shortening_direction({0, 0, -.05}), pi / 4, 1e-15);
    EXPECT_EQ(shortening_m({.1, .2, 0}, 1000), 0);
    EXPECT_NEAR(shortening_direction({.1, .2, 0}), 0, 1e-15); // the least extension's axis
    EXPECT_NEAR(shortening_m({-.1, -.1, 0}, 1000), 100, 1e-12);
    EXPECT_EQ(shortening_direction({-.1, -.1, 0}), 0);
    EXPECT_EQ(shortening_m({}, 1000), 0);
    EXPECT_EQ(shortening_direction({}), 0);
}

namespace {
// Land columns at rest except for the displacement: no erosion, no roots, no forcing.
Model strained_world(uint32_t subdivision) {
    Model world;
    world.initialize(16, 8, 100000, std::vector<Cell>(128, seed(500, true, 100)), subdivision);
    return world;
}
void advance_still(Model& world, const std::vector<Displacement>& u, const tectonics::margins::Forcing* forcing = nullptr) {
    tectonics::margins::Parameters still;
    still.root_mobility = 0;
    world.advance(u, 1, 0, 1, 1, forcing, still, {}, false, nullptr, 5);
}
// Periodic triangle profile of period 2 n: v at 0, -v at n, slope -v / (n / 2) between.
double triangle(int x, int n, double v) {
    const int p = ((x % (2 * n)) + 2 * n) % (2 * n);
    return p <= n ? v - 2 * v * p / n : -v + 2 * v * (p - n) / n;
}
}

// A uniform convergence along x: u_x a periodic triangle (0.2 at x = 4, -0.2 at x = 12,
// slope -0.05 cell per cell between), so each update shortens the rock of x = 5..11 by
// 0.05 along x: 5 km per 100 km cell, direction 0, nothing along y. Rock arriving from
// the kinks (x = 4, 12: derivative 0) mixes into x = 5 and 11, so x = 6..10 after one
// update, x = 7..9 after two (10 km); the diverging half x = 13..3 shortens nothing.
// Along y the same with period 8 (0.1 at y = 2, -0.1 at y = 6): row 4 shortens by 5 km
// per update at direction pi / 2. Subdivisions 1 and 4.
TEST(TectonicColumns, ShorteningOfAUniformConvergence) {
    const double pi = std::acos(-1.0);
    for (uint32_t subdivision : {1u, 4u}) {
        std::vector<Displacement> u(128), v(128);
        for (int i = 0; i < 128; ++i) {
            u[i].x = triangle(i % 16 - 4, 8, .2);
            v[i].y = triangle(i / 16 - 2, 4, .1);
        }
        EXPECT_NEAR(u[5].x - u[4].x, -.05, 1e-15);
        EXPECT_NEAR(v[3 * 16].y - v[2 * 16].y, -.05, 1e-15);
        auto world = strained_world(subdivision);
        advance_still(world, u);
        for (size_t i = 0; i < 128; ++i) {
            const auto x = i % 16;
            const auto& s = world.strain()[i];
            if (x >= 6 && x <= 10) {
                EXPECT_NEAR(s.xx, -.05, 1e-12) << "subdivision " << subdivision << " cell " << i;
                EXPECT_NEAR(shortening_m(s, 100000), 5000, 1e-7);
                EXPECT_NEAR(shortening_direction(s), 0, 1e-12);
            }
            if (x >= 14 || x <= 2) EXPECT_EQ(shortening_m(s, 100000), 0) << "cell " << i;
            EXPECT_NEAR(s.yy, 0, 1e-12);
            EXPECT_NEAR(s.xy, 0, 1e-12);
        }
        advance_still(world, u);
        for (size_t x = 7; x <= 9; ++x) for (size_t y = 0; y < 8; ++y)
            EXPECT_NEAR(shortening_m(world.strain()[y * 16 + x], 100000), 10000, 1e-7)
                << "subdivision " << subdivision << " x " << x << " y " << y;
        auto column = strained_world(subdivision);
        advance_still(column, v);
        for (size_t x = 0; x < 16; ++x) {
            const auto& s = column.strain()[4 * 16 + x];
            EXPECT_NEAR(shortening_m(s, 100000), 5000, 1e-7) << "subdivision " << subdivision << " x " << x;
            EXPECT_NEAR(shortening_direction(s), pi / 2, 1e-12);
        }
    }
}

// The tensor moves with the rock. u_x = 0.1 for x <= 4, 0 at x = 5, -0.1 for x >= 6
// strains column 5 alone (its two one-sided differences are -0.1; elsewhere one of
// them is 0); converging rock dilutes it. A second update translates everything by
// 0.25 cell east (no strain, no spin): the sum over a row is conserved and the
// shortening moves east by the translation (its centroid by 0.25 cell), keeping its
// direction 0. Subdivision 1.
TEST(TectonicColumns, ShorteningIsCarriedWithTheRock) {
    std::vector<Displacement> u(128);
    for (size_t i = 0; i < 128; ++i) u[i].x = i % 16 <= 4 ? .1 : i % 16 == 5 ? 0 : -.1;
    auto world = strained_world(1);
    advance_still(world, u);
    const auto row = [&](size_t y) {
        double total = 0, moment = 0;
        for (size_t x = 0; x < 16; ++x) {
            const double amount = shortening_m(world.strain()[y * 16 + x], 100000);
            total += amount; moment += amount * x;
        }
        return std::pair{total, moment / total};
    };
    const auto [first, before] = row(3);
    EXPECT_GT(shortening_m(world.strain()[3 * 16 + 5], 100000), 8000);
    for (size_t x : {0u, 1u, 2u, 3u, 9u, 10u, 11u, 12u, 13u, 14u, 15u})
        EXPECT_EQ(shortening_m(world.strain()[3 * 16 + x], 100000), 0) << "x " << x;
    advance_still(world, std::vector<Displacement>(128, {.25, 0}));
    for (size_t y = 0; y < 8; ++y) {
        const auto [total, centroid] = row(y);
        EXPECT_NEAR(total, first, 1e-6) << "row " << y;
        EXPECT_NEAR(centroid, before + .25, 1e-9) << "row " << y;
    }
    for (const auto& s : world.strain()) if (shortening_m(s, 100000) > 1) EXPECT_NEAR(shortening_direction(s), 0, 1e-12);
}

// The tensor turns with the plate. After one update of the uniform convergence along x
// (x = 6..10 shortened by 0.05 at direction 0), a rigid rotation by w = 0.05 rad about
// (8, 4) (u = w (-(y - 4), x - 8), with its jumps at the wrap seams inside the plate)
// adds no strain and turns the axis by w: x = 7..9 keep 5 km at direction 0.05 in
// every row, also next to the seam between rows 7 and 0, where the smaller one-sided
// difference ignores the jump of 8 w.
TEST(TectonicColumns, ShorteningTurnsWithAPlateRotation) {
    const double w = .05;
    for (uint32_t subdivision : {1u, 4u}) {
        std::vector<Displacement> u(128), rotation(128);
        for (int i = 0; i < 128; ++i) {
            u[i].x = triangle(i % 16 - 4, 8, .2);
            rotation[i] = {-w * (i / 16 - 4), w * (i % 16 - 8)};
        }
        auto world = strained_world(subdivision);
        advance_still(world, u);
        advance_still(world, rotation);
        for (size_t x = 7; x <= 9; ++x) for (size_t y = 0; y < 8; ++y) {
            const auto& s = world.strain()[y * 16 + x];
            EXPECT_NEAR(shortening_m(s, 100000), 5000, 1e-7) << "subdivision " << subdivision << " x " << x << " y " << y;
            EXPECT_NEAR(shortening_direction(s), w, 1e-12) << "subdivision " << subdivision << " x " << x << " y " << y;
        }
    }
}

// Plate edges: a jump of velocity to another plate is not strain. Plate 0 (x = 0..7)
// moves east by 0.2 cell and slows to 0.1 at x = 6 and 0 at x = 7; plate 1 (x = 8..15)
// moves east by 0.1. With the owners, x = 6 and 7 take the in-plate difference -0.1
// (10 km; x = 6 keeps 9 / 11 of it, mixed with unstrained rock arriving from x = 5),
// plate 1 and the rest of plate 0 nothing, though the boundary opens between x = 7
// and 8 and closes across the wrap. Without the owners (one plate) the smaller
// difference at x = 7 is the mean of -0.1 and +0.1, 0: only x = 6 is strained, and
// part of its rock reaches x = 7. Rigid plates converging shorten nothing.
TEST(TectonicColumns, ShorteningIgnoresJumpsBetweenPlates) {
    tectonics::margins::Forcing plates;
    plates.thermal_age_myr.assign(128, -1);
    plates.trench_m.assign(128, 0);
    plates.arc_receiver.assign(128, std::numeric_limits<uint32_t>::max());
    plates.owners.resize(128);
    std::vector<Displacement> u(128), rigid(128);
    for (size_t i = 0; i < 128; ++i) {
        const size_t x = i % 16;
        plates.owners[i] = x < 8 ? 0 : 1;
        u[i].x = x < 6 ? .2 : x == 6 ? .1 : x == 7 ? 0 : .1;
        rigid[i].x = x < 8 ? .2 : -.2;
    }
    auto owned = strained_world(1), single = strained_world(1), converging = strained_world(1);
    advance_still(owned, u, &plates);
    advance_still(single, u);
    advance_still(converging, rigid, &plates);
    for (size_t y = 0; y < 8; ++y) {
        const auto at = [&](const Model& world, size_t x) { return shortening_m(world.strain()[y * 16 + x], 100000); };
        EXPECT_NEAR(at(owned, 7), 10000, 1e-6) << "row " << y;
        EXPECT_NEAR(at(owned, 6), 10000. * 9 / 11, 1e-6) << "row " << y; // 0.9 of its rock stays, 0.2 of x = 5's arrives
        for (size_t x : {0u, 1u, 2u, 3u, 4u, 8u, 9u, 10u, 11u, 12u, 13u, 14u, 15u})
            EXPECT_NEAR(at(owned, x), 0, 1e-6) << "x " << x; // rounding of the carried shift
        EXPECT_NEAR(at(single, 7), 10000. / 11, 1e-6) << "row " << y; // 0.1 of x = 6's rock, coverage 1.1
        for (size_t x = 0; x < 16; ++x) EXPECT_EQ(at(converging, x), 0) << "x " << x;
    }
}

// New crust has no strain. Ocean floor is first shortened by 0.05 along y in row 4 and
// partly in rows 3 and 5 (the triangle profile of ShorteningOfAUniformConvergence, one
// plate), across every column, the future ridge included. Then two rigid plates part
// for 50 updates at subdivision 4: x = 0..7 moves west and x = 8..15 east by 0.2 cell,
// so a ridge opens between x = 7 and 8 and the floor subducts across the wrap. The
// crust the ridge creates carries no strain and the rock that stayed is diluted by it,
// so the ridge columns and the floor they made (x = 4..11, all of it created in the
// last 20 updates) keep almost none of the 5 km: below 1 m at the ridge, at most
// 28.3 m (0.6 percent) in the new floor, where the remap smears in the older rock that
// left the ridge first. Without the zero strain of the created floor both keep 5 km.
TEST(TectonicColumns, ShorteningOfNewCrustAtADivergentBoundaryIsZero) {
    tectonics::margins::Forcing plates;
    plates.thermal_age_myr.assign(128, -1);
    plates.trench_m.assign(128, 0);
    plates.arc_receiver.assign(128, std::numeric_limits<uint32_t>::max());
    plates.owners.resize(128);
    std::vector<Displacement> squeeze(128), apart(128);
    for (size_t i = 0; i < 128; ++i) {
        plates.owners[i] = i % 16 < 8 ? 0 : 1;
        squeeze[i].y = triangle(static_cast<int>(i / 16) - 2, 4, .1);
        apart[i].x = i % 16 < 8 ? -.2 : .2;
    }
    Model world;
    world.initialize(16, 8, 100000, std::vector<Cell>(128, seed(-4000, false, 20)), 4);
    advance_still(world, squeeze);
    for (size_t x = 0; x < 16; ++x) ASSERT_NEAR(shortening_m(world.strain()[4 * 16 + x], 100000), 5000, 1e-6);
    for (int step = 0; step < 50; ++step) advance_still(world, apart, &plates);
    double ridge = 0, floor = 0;
    for (size_t i = 0; i < 128; ++i) {
        const auto x = i % 16;
        const double amount = shortening_m(world.strain()[i], 100000);
        if (x == 7 || x == 8) ridge = std::max(ridge, amount);
        else if (x >= 4 && x <= 11) floor = std::max(floor, amount);
    }
    EXPECT_LT(ridge, 1) << "largest shortening at the ridge, m";
    EXPECT_LT(floor, 40) << "largest shortening of the new floor, m"; // measured 28.3 m
}

// Over several updates of a varied displacement (two plates, subdivision 4) the
// shortening stays finite and non-negative and the direction in [0, pi).
TEST(TectonicColumns, ShorteningIsNonNegativeAndItsDirectionInRange) {
    const double pi = std::acos(-1.0);
    tectonics::margins::Forcing plates;
    plates.thermal_age_myr.assign(128, -1);
    plates.trench_m.assign(128, 0);
    plates.arc_receiver.assign(128, std::numeric_limits<uint32_t>::max());
    plates.owners.resize(128);
    std::vector<Displacement> u(128);
    for (size_t i = 0; i < 128; ++i) {
        const double x = static_cast<double>(i % 16), y = static_cast<double>(i / 16);
        plates.owners[i] = (i % 16 + i / 16) % 16 < 8 ? 0 : 1;
        u[i] = {.3 * std::sin(2 * pi * x / 16 + y), .2 * std::cos(2 * pi * y / 8 - x / 3)};
    }
    auto world = strained_world(4);
    double largest = 0;
    for (int step = 0; step < 6; ++step) {
        advance_still(world, u, &plates);
        for (const auto& s : world.strain()) {
            const double amount = shortening_m(s, 100000), direction = shortening_direction(s);
            ASSERT_TRUE(std::isfinite(amount));
            EXPECT_GE(amount, 0);
            EXPECT_GE(direction, 0);
            EXPECT_LT(direction, pi);
            largest = std::max(largest, amount);
        }
    }
    EXPECT_GT(largest, 0);
}
