#ifndef NOMINMAX
#define NOMINMAX
#endif

#include "gtest/gtest.h"
#include "lithosphere.hpp"
#include "platecapi.hpp"
#include "tectonic_pipeline_support.hpp"

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

platec::scenario::Scenario make_phase7_scenario(long seed) {
    platec::scenario::Scenario scenario;
    scenario.seed = seed;
    scenario.width = 96;
    scenario.height = 48;
    scenario.sea_level = 0.65f;
    scenario.erosion_period = 0;
    scenario.folding_ratio = 0.02f;
    scenario.aggregation_overlap_abs = 1000000;
    scenario.aggregation_overlap_rel = 0.33f;
    scenario.cycle_count = 4;
    scenario.plate_count = 4;
    scenario.cycle_step_limit = 600;
    scenario.delta_time_myr = 1.0;
    return scenario;
}

platec::scenario::Scenario make_phase7_restart_scenario(long seed) {
    platec::scenario::Scenario scenario = make_phase7_scenario(seed);
    scenario.cycle_step_limit = 1;
    scenario.cycle_count = 3;
    return scenario;
}

platec::contract::Snapshot run_snapshot(long seed, uint32_t updates) {
    lithosphere simulation(make_phase7_scenario(seed));
    for (uint32_t i = 0; i < updates; ++i) {
        if (simulation.isFinished()) {
            ADD_FAILURE() << "simulation finished before requested update count";
            break;
        }
        simulation.update();
    }
    return platec::contract::capture_snapshot(simulation);
}

platec::contract::Snapshot run_snapshot(const platec::scenario::Scenario& scenario,
                                        uint32_t updates) {
    lithosphere simulation(scenario);
    for (uint32_t i = 0; i < updates; ++i) {
        if (simulation.isFinished()) {
            ADD_FAILURE() << "simulation finished before requested update count";
            break;
        }
        simulation.update();
    }
    return platec::contract::capture_snapshot(simulation);
}

fs::path temp_root() {
    return fs::temp_directory_path() / "plate-tectonics-phase7-tooling";
}

std::string read_text_file(const fs::path& path) {
    std::ifstream input(path, std::ios::binary);
    std::ostringstream contents;
    contents << input.rdbuf();
    return contents.str();
}

platec::contract::Snapshot make_snapshot_with_junction() {
    platec::contract::Snapshot snapshot;
    snapshot.width = 2;
    snapshot.height = 2;
    snapshot.sea_level_m = 1000;
    snapshot.run_scenario = make_phase7_scenario(4242);
    snapshot.run_scenario.width = snapshot.width;
    snapshot.run_scenario.height = snapshot.height;
    snapshot.run_scenario.plate_count = 3;

    const size_t cell_count = snapshot.cell_count();
    snapshot.heightmap.assign(cell_count, 0.5f);
    snapshot.plate_id = {0U, 1U, 2U, 0U};
    snapshot.crust_age_steps.assign(cell_count, 0U);
    snapshot.crust_age_myr.assign(cell_count, 0.0f);
    snapshot.crust_thickness.assign(cell_count, 1.0f);
    snapshot.crust_class.assign(
        cell_count, static_cast<uint8_t>(platec::contract::CrustClass::Oceanic));
    snapshot.uplift_tendency.assign(cell_count, 0.0f);
    snapshot.subsidence_tendency.assign(cell_count, 0.0f);
    snapshot.accumulated_strain.assign(cell_count, 0.0f);
    snapshot.boundary_type.assign(
        cell_count, static_cast<uint8_t>(platec::contract::BoundaryType::Divergent));
    snapshot.boundary_distance.assign(cell_count, 0U);
    snapshot.boundary_segment_id.assign(cell_count, platec::contract::kNoBoundaryId);
    snapshot.nearest_boundary_id.assign(cell_count, platec::contract::kNoBoundaryId);
    snapshot.deforming_region_id.assign(cell_count, platec::contract::kNoDeformingRegionId);
    snapshot.deforming_region_type.assign(
        cell_count, static_cast<uint8_t>(platec::contract::DeformingRegionType::None));
    snapshot.deformation_rate.assign(cell_count, 0.0f);
    snapshot.deformation_velocity_x.assign(cell_count, 0.0f);
    snapshot.deformation_velocity_y.assign(cell_count, 0.0f);
    snapshot.convergence_score.assign(cell_count, static_cast<uint8_t>(0));
    snapshot.divergence_score.assign(cell_count, static_cast<uint8_t>(255));
    snapshot.shear_score.assign(cell_count, static_cast<uint8_t>(0));
    snapshot.geologic_regime.assign(
        cell_count, static_cast<uint8_t>(platec::contract::GeologicRegime::MidOceanRidge));
    snapshot.plates.resize(3);

    platec::contract::TectonicJunction junction;
    junction.id = 7U;
    junction.plate_ids[0] = 0U;
    junction.plate_ids[1] = 1U;
    junction.plate_ids[2] = 2U;
    junction.boundary_segment_ids[0] = 10U;
    junction.boundary_segment_ids[1] = 20U;
    junction.boundary_segment_ids[2] = 30U;
    junction.arm_types[0] = platec::contract::JunctionArmType::Ridge;
    junction.arm_types[1] = platec::contract::JunctionArmType::Fault;
    junction.arm_types[2] = platec::contract::JunctionArmType::Trench;
    junction.arm_strike_x[0] = 1.0f;
    junction.arm_strike_y[1] = 1.0f;
    junction.arm_strike_x[2] = 0.70710678f;
    junction.arm_strike_y[2] = 0.70710678f;
    junction.arm_orientation_quality[0] = 0.9f;
    junction.arm_orientation_quality[1] = 0.8f;
    junction.arm_orientation_quality[2] = 0.7f;
    junction.candidate_cell_count = 3U;
    junction.centroid_x = 0.5f;
    junction.centroid_y = 1.5f;
    junction.junction_velocity_x = 0.25f;
    junction.junction_velocity_y = -0.5f;
    junction.velocity_closure_error = 0.125f;
    junction.velocity_quality_score = 0.8f;
    junction.velocity_constraint_error = 0.25f;
    junction.geometry_quality_score = 0.7f;
    junction.stability = platec::contract::JunctionStability::ConditionallyStable;
    junction.stability_score = 0.6f;
    junction.quality_score = 0.5f;
    snapshot.tectonic_junctions.push_back(junction);

    return snapshot;
}

} // namespace

TEST(Phase7Tooling, SnapshotBundleRoundTripsCoreRasterFields) {
    const platec::contract::Snapshot snapshot = run_snapshot(12345, 6);
    const fs::path output_dir = temp_root() / "roundtrip";
    std::error_code ec;
    fs::remove_all(output_dir, ec);

    platec::tooling::SnapshotBundleOptions options;
    options.label = "update000006";
    options.source = "test";
    options.include_update_count = true;
    options.update_count = 6;
    platec::tooling::write_snapshot_bundle(output_dir, snapshot, options);

    ASSERT_TRUE(fs::exists(output_dir / "manifest.json"));
    ASSERT_TRUE(fs::exists(output_dir / "heightmap.f32"));
    ASSERT_TRUE(fs::exists(output_dir / "boundary_segments.json"));
    ASSERT_TRUE(fs::exists(output_dir / "deforming_regions.json"));
    ASSERT_TRUE(fs::exists(output_dir / "tectonic_junctions.json"));

    const platec::tooling::LoadedSnapshotBundle loaded =
        platec::tooling::load_snapshot_bundle(output_dir);
    EXPECT_EQ(loaded.schema, "tectonic-snapshot-bundle/v1");
    EXPECT_EQ(loaded.label, "update000006");
    EXPECT_EQ(loaded.source, "test");
    EXPECT_EQ(loaded.snapshot.schema_version, snapshot.schema_version);
    EXPECT_EQ(loaded.snapshot.width, snapshot.width);
    EXPECT_EQ(loaded.snapshot.height, snapshot.height);
    EXPECT_EQ(loaded.snapshot.iteration_count, snapshot.iteration_count);
    EXPECT_EQ(loaded.snapshot.cycle_count, snapshot.cycle_count);
    EXPECT_EQ(loaded.snapshot.time_origin_step, snapshot.time_origin_step);
    EXPECT_EQ(loaded.snapshot.sea_level_m, snapshot.sea_level_m);
    EXPECT_EQ(loaded.snapshot.time_myr, snapshot.time_myr);
    EXPECT_EQ(loaded.snapshot.delta_time_myr, snapshot.delta_time_myr);
    EXPECT_EQ(loaded.snapshot.heightmap, snapshot.heightmap);
    EXPECT_EQ(loaded.snapshot.plate_id, snapshot.plate_id);
    EXPECT_EQ(loaded.snapshot.crust_age_steps, snapshot.crust_age_steps);
    EXPECT_EQ(loaded.snapshot.crust_age_myr, snapshot.crust_age_myr);
    EXPECT_EQ(loaded.snapshot.boundary_type, snapshot.boundary_type);
    EXPECT_EQ(loaded.snapshot.boundary_segment_id, snapshot.boundary_segment_id);
    EXPECT_EQ(loaded.snapshot.deforming_region_id, snapshot.deforming_region_id);
    EXPECT_EQ(loaded.declared_boundary_segment_count, snapshot.boundary_segments.size());
    EXPECT_EQ(loaded.declared_deforming_region_count, snapshot.deforming_regions.size());
    EXPECT_EQ(loaded.declared_tectonic_junction_count, snapshot.tectonic_junctions.size());

    fs::remove_all(output_dir, ec);
}

TEST(Phase7Tooling, TectonicJunctionJsonIncludesDiagnostics) {
    const platec::contract::Snapshot snapshot = make_snapshot_with_junction();
    const fs::path output_dir = temp_root() / "junction-json";
    std::error_code ec;
    fs::remove_all(output_dir, ec);

    platec::tooling::SnapshotBundleOptions options;
    options.label = "junction-diagnostics";
    options.source = "test";
    platec::tooling::write_snapshot_bundle(output_dir, snapshot, options);

    const std::string junction_json = read_text_file(output_dir / "tectonic_junctions.json");
    const std::vector<std::string> required_fields = {
        "\"id\"",
        "\"plate_ids\"",
        "\"boundary_segment_ids\"",
        "\"arm_types\"",
        "\"arm_strikes\"",
        "\"velocity_closure_error\"",
        "\"velocity_quality_score\"",
        "\"velocity_constraint_error\"",
        "\"geometry_quality_score\"",
        "\"stability\"",
        "\"stability_score\"",
        "\"quality_score\"",
    };
    for (const std::string& field : required_fields) {
        EXPECT_NE(junction_json.find(field), std::string::npos) << field;
    }

    fs::remove_all(output_dir, ec);
}

TEST(Phase7Tooling, SummaryExposesBoundaryAndDeformingStats) {
    const platec::contract::Snapshot snapshot = run_snapshot(12345, 6);
    const platec::tooling::SnapshotSummary summary =
        platec::tooling::summarize_snapshot(snapshot);

    EXPECT_EQ(summary.contract_schema_version, platec::contract::kPhase6SchemaVersion);
    EXPECT_EQ(summary.width, snapshot.width);
    EXPECT_EQ(summary.height, snapshot.height);
    EXPECT_GT(summary.active_convergence_cells, 0U);
    EXPECT_GT(summary.boundary.segment_count, 0U);
    EXPECT_GT(summary.boundary.active_cell_count, 0U);
    EXPECT_GT(summary.boundary.nearest_mapped_cell_count, 0U);
    EXPECT_GT(summary.deforming.region_count, 0U);
    EXPECT_GT(summary.deforming.active_cell_count, 0U);

    const std::string summary_json = platec::tooling::snapshot_summary_json(summary);
    const std::string boundary_json = platec::tooling::boundary_stats_json(summary);
    EXPECT_NE(summary_json.find("\"tectonic-snapshot-summary/v1\""), std::string::npos);
    EXPECT_NE(summary_json.find("\"segment_count\""), std::string::npos);
    EXPECT_NE(boundary_json.find("\"tectonic-boundary-stats/v1\""), std::string::npos);
    EXPECT_NE(boundary_json.find("\"deforming\""), std::string::npos);
}

TEST(Phase7Tooling, CApiTectonicJunctionsMatchCapturedSnapshot) {
    const platec::scenario::Scenario scenario = make_phase7_scenario(12345);
    void* simulation = platec_api_create_from_scenario(&scenario);
    ASSERT_NE(simulation, nullptr);

    while (platec_api_is_finished(simulation) == 0U) {
        platec_api_step(simulation);
    }

    const auto& litho = *static_cast<const lithosphere*>(simulation);
    const platec::contract::Snapshot snapshot = platec::contract::capture_snapshot(litho);
    const uint32_t api_count = platec_api_get_tectonic_junction_count(simulation);
    const auto* api_junctions = platec_api_get_tectonic_junctions(simulation);

    ASSERT_EQ(api_count, snapshot.tectonic_junctions.size());
    ASSERT_GT(api_count, 0U);
    ASSERT_NE(api_junctions, nullptr);

    const auto& expected = snapshot.tectonic_junctions.front();
    const auto& actual = api_junctions[0];
    EXPECT_EQ(actual.id, expected.id);
    EXPECT_EQ(actual.plate_ids[0], expected.plate_ids[0]);
    EXPECT_EQ(actual.plate_ids[1], expected.plate_ids[1]);
    EXPECT_EQ(actual.plate_ids[2], expected.plate_ids[2]);
    EXPECT_EQ(actual.boundary_segment_ids[0], expected.boundary_segment_ids[0]);
    EXPECT_EQ(actual.boundary_segment_ids[1], expected.boundary_segment_ids[1]);
    EXPECT_EQ(actual.boundary_segment_ids[2], expected.boundary_segment_ids[2]);
    EXPECT_EQ(actual.arm_types[0], expected.arm_types[0]);
    EXPECT_EQ(actual.arm_types[1], expected.arm_types[1]);
    EXPECT_EQ(actual.arm_types[2], expected.arm_types[2]);
    EXPECT_FLOAT_EQ(actual.centroid_x, expected.centroid_x);
    EXPECT_FLOAT_EQ(actual.centroid_y, expected.centroid_y);
    EXPECT_FLOAT_EQ(actual.velocity_quality_score, expected.velocity_quality_score);
    EXPECT_EQ(actual.stability, expected.stability);
    EXPECT_FLOAT_EQ(actual.quality_score, expected.quality_score);

    platec_api_destroy(simulation);
}

TEST(Phase7Tooling, ComparisonDetectsSameAndDifferentSnapshots) {
    const platec::contract::Snapshot snapshot_a = run_snapshot(12345, 6);
    const platec::contract::Snapshot snapshot_b = run_snapshot(67890, 6);

    const platec::tooling::SnapshotComparison identical =
        platec::tooling::compare_snapshots(snapshot_a, snapshot_a);
    EXPECT_TRUE(identical.dimensions_match);
    EXPECT_TRUE(identical.exact_match);
    EXPECT_EQ(identical.heightmap.changed_cells, 0U);
    EXPECT_EQ(identical.plate_id.changed_cells, 0U);
    EXPECT_EQ(identical.geologic_regime.changed_cells, 0U);

    const platec::tooling::SnapshotComparison different =
        platec::tooling::compare_snapshots(snapshot_a, snapshot_b);
    EXPECT_TRUE(different.dimensions_match);
    EXPECT_FALSE(different.exact_match);
    EXPECT_GT(different.heightmap.changed_cells, 0U);
    EXPECT_GT(different.plate_id.changed_cells, 0U);
    EXPECT_GT(different.geologic_regime.changed_cells, 0U);

    const std::string comparison_json =
        platec::tooling::snapshot_comparison_json(different);
    EXPECT_NE(comparison_json.find("\"tectonic-snapshot-comparison/v1\""), std::string::npos);
    EXPECT_NE(comparison_json.find("\"heightmap\""), std::string::npos);
}

TEST(Phase7Tooling, RestartedSnapshotsKeepWholeRunTimeInBundlesAndSummaries) {
    const platec::contract::Snapshot snapshot =
        run_snapshot(make_phase7_restart_scenario(12345), 3);

    EXPECT_EQ(snapshot.cycle_count, 1U);
    EXPECT_DOUBLE_EQ(snapshot.time_myr, 2.0);

    const platec::tooling::SnapshotSummary summary =
        platec::tooling::summarize_snapshot(snapshot);
    EXPECT_EQ(summary.cycle_count, 1U);
    EXPECT_DOUBLE_EQ(summary.time_myr, 2.0);

    const fs::path output_dir = temp_root() / "restart-time";
    std::error_code ec;
    fs::remove_all(output_dir, ec);

    platec::tooling::SnapshotBundleOptions options;
    options.label = "time2_000Myr";
    options.source = "test";
    options.include_update_count = true;
    options.update_count = 3;
    platec::tooling::write_snapshot_bundle(output_dir, snapshot, options);

    const platec::tooling::LoadedSnapshotBundle loaded =
        platec::tooling::load_snapshot_bundle(output_dir);
    EXPECT_EQ(loaded.snapshot.cycle_count, 1U);
    EXPECT_DOUBLE_EQ(loaded.snapshot.time_myr, 2.0);

    fs::remove_all(output_dir, ec);
}
