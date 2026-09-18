#ifndef TECTONIC_PIPELINE_SUPPORT_HPP
#define TECTONIC_PIPELINE_SUPPORT_HPP

#include "tectonic_contract.hpp"

#include <array>
#include <cstdint>
#include <filesystem>
#include <string>

namespace platec::tooling {

struct SnapshotBundleOptions {
    std::string schema = "tectonic-snapshot-bundle/v1";
    std::string label = "snapshot";
    std::string source = "runtime";
    bool include_update_count = false;
    uint32_t update_count = 0;
};

struct BoundaryStatsSummary {
    uint32_t segment_count = 0;
    uint32_t convergent_segment_count = 0;
    uint32_t divergent_segment_count = 0;
    uint32_t transform_segment_count = 0;
    uint32_t passive_margin_segment_count = 0;
    uint32_t active_cell_count = 0;
    uint32_t convergent_cell_count = 0;
    uint32_t divergent_cell_count = 0;
    uint32_t transform_cell_count = 0;
    uint32_t passive_margin_cell_count = 0;
    uint32_t nearest_mapped_cell_count = 0;
    double max_segment_age_myr = 0.0;
    double mean_segment_age_myr = 0.0;
    double mean_segment_persistence_steps = 0.0;
};

struct DeformingRegionStatsSummary {
    uint32_t region_count = 0;
    uint32_t continental_rift_region_count = 0;
    uint32_t diffuse_collision_region_count = 0;
    uint32_t active_cell_count = 0;
    uint32_t continental_rift_cell_count = 0;
    uint32_t diffuse_collision_cell_count = 0;
    double max_region_age_myr = 0.0;
    double mean_region_age_myr = 0.0;
    double mean_region_persistence_steps = 0.0;
    double mean_deformation_rate = 0.0;
};

struct SnapshotSummary {
    uint32_t contract_schema_version = 0;
    uint32_t width = 0;
    uint32_t height = 0;
    uint32_t iteration_count = 0;
    uint32_t cycle_count = 0;
    uint32_t time_origin_step = 0;
    uint16_t sea_level_m = 0;
    double time_myr = 0.0;
    double delta_time_myr = 0.0;
    uint32_t plate_count = 0;
    uint32_t tectonic_junction_count = 0;
    std::array<uint32_t, 5> tectonic_junction_stability_counts{};
    std::array<uint32_t, 8> regime_counts{};
    std::array<uint32_t, 4> crust_class_counts{};
    uint32_t active_convergence_cells = 0;
    uint32_t active_divergence_cells = 0;
    uint32_t active_shear_cells = 0;
    BoundaryStatsSummary boundary;
    DeformingRegionStatsSummary deforming;
    std::string heightmap_hash;
    std::string plate_id_hash;
    std::string geologic_regime_hash;
    std::string boundary_segment_id_hash;
    std::string deforming_region_id_hash;
    std::string tectonic_junctions_hash;
};

struct LoadedSnapshotBundle {
    platec::contract::Snapshot snapshot;
    std::string schema;
    std::string label;
    std::string source;
    uint32_t declared_boundary_segment_count = 0;
    uint32_t declared_deforming_region_count = 0;
    uint32_t declared_tectonic_junction_count = 0;
};

struct NumericFieldDiff {
    uint64_t changed_cells = 0;
    double max_abs_diff = 0.0;
    double mean_abs_diff = 0.0;
    double rmse = 0.0;
};

struct ExactFieldDiff {
    uint64_t changed_cells = 0;
};

struct SnapshotComparison {
    bool dimensions_match = false;
    bool exact_match = false;
    uint32_t width = 0;
    uint32_t height = 0;
    uint32_t contract_schema_version_a = 0;
    uint32_t contract_schema_version_b = 0;
    double time_myr_a = 0.0;
    double time_myr_b = 0.0;
    int64_t iteration_count_delta = 0;
    int64_t cycle_count_delta = 0;
    int64_t plate_count_delta = 0;
    int64_t boundary_segment_count_delta = 0;
    int64_t deforming_region_count_delta = 0;
    int64_t tectonic_junction_count_delta = 0;
    NumericFieldDiff heightmap;
    ExactFieldDiff plate_id;
    ExactFieldDiff crust_age_steps;
    NumericFieldDiff crust_age_myr;
    NumericFieldDiff crust_thickness;
    ExactFieldDiff crust_class;
    NumericFieldDiff uplift_tendency;
    NumericFieldDiff subsidence_tendency;
    NumericFieldDiff accumulated_strain;
    ExactFieldDiff boundary_type;
    NumericFieldDiff boundary_distance;
    ExactFieldDiff boundary_segment_id;
    ExactFieldDiff nearest_boundary_id;
    ExactFieldDiff deforming_region_id;
    ExactFieldDiff deforming_region_type;
    NumericFieldDiff deformation_rate;
    NumericFieldDiff deformation_velocity_x;
    NumericFieldDiff deformation_velocity_y;
    ExactFieldDiff convergence_score;
    ExactFieldDiff divergence_score;
    ExactFieldDiff shear_score;
    ExactFieldDiff geologic_regime;
};

void write_snapshot_bundle(const std::filesystem::path& output_dir,
                           const platec::contract::Snapshot& snapshot,
                           const SnapshotBundleOptions& options);
LoadedSnapshotBundle load_snapshot_bundle(const std::filesystem::path& input_dir);
SnapshotSummary summarize_snapshot(const platec::contract::Snapshot& snapshot);
SnapshotComparison compare_snapshots(const platec::contract::Snapshot& a,
                                     const platec::contract::Snapshot& b,
                                     uint32_t boundary_segment_count_a = 0,
                                     uint32_t boundary_segment_count_b = 0,
                                     uint32_t deforming_region_count_a = 0,
                                     uint32_t deforming_region_count_b = 0,
                                     uint32_t tectonic_junction_count_a = 0,
                                     uint32_t tectonic_junction_count_b = 0);
std::string snapshot_summary_json(const SnapshotSummary& summary);
std::string boundary_stats_json(const SnapshotSummary& summary);
std::string snapshot_comparison_json(const SnapshotComparison& comparison);

} // namespace platec::tooling

#endif
