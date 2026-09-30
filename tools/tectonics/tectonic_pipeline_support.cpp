#include "tectonic_pipeline_support.hpp"

#include "topography_codec.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace tectonics::tooling {
namespace {

[[noreturn]] void fail(const std::string& message) {
    throw std::runtime_error(message);
}

std::string json_escape(const std::string& value) {
    std::string escaped;
    escaped.reserve(value.size());
    for (char ch : value) {
        switch (ch) {
        case '\\':
            escaped += "\\\\";
            break;
        case '\"':
            escaped += "\\\"";
            break;
        case '\n':
            escaped += "\\n";
            break;
        case '\r':
            escaped += "\\r";
            break;
        case '\t':
            escaped += "\\t";
            break;
        default:
            escaped += ch;
            break;
        }
    }
    return escaped;
}

template <typename T>
void write_binary_file(const fs::path& path, const std::vector<T>& values) {
    std::ofstream output(path, std::ios::binary);
    if (!output) {
        fail("failed to open " + path.string());
    }
    output.write(reinterpret_cast<const char*>(values.data()),
                 static_cast<std::streamsize>(values.size() * sizeof(T)));
    if (!output) {
        fail("failed to write " + path.string());
    }
}

void write_text_file(const fs::path& path, const std::string& contents) {
    std::ofstream output(path, std::ios::binary);
    if (!output) {
        fail("failed to open " + path.string());
    }
    output << contents;
    if (!output) {
        fail("failed to write " + path.string());
    }
}

std::string read_text_file(const fs::path& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        fail("failed to open " + path.string());
    }
    std::ostringstream contents;
    contents << input.rdbuf();
    if (!input.good() && !input.eof()) {
        fail("failed to read " + path.string());
    }
    return contents.str();
}

template <typename T>
std::vector<T> read_binary_file(const fs::path& path, size_t expected_count) {
    std::ifstream input(path, std::ios::binary | std::ios::ate);
    if (!input) {
        fail("failed to open " + path.string());
    }

    const std::streamsize size = input.tellg();
    const std::streamsize expected_size =
        static_cast<std::streamsize>(expected_count * sizeof(T));
    if (size != expected_size) {
        fail("unexpected size for " + path.string());
    }

    input.seekg(0, std::ios::beg);
    std::vector<T> values(expected_count);
    input.read(reinterpret_cast<char*>(values.data()), expected_size);
    if (!input) {
        fail("failed to read " + path.string());
    }
    return values;
}

template <typename T>
std::string fnv1a_hash(const std::vector<T>& values) {
    constexpr uint64_t kOffset = 14695981039346656037ull;
    constexpr uint64_t kPrime = 1099511628211ull;

    uint64_t hash = kOffset;
    const auto* bytes = reinterpret_cast<const uint8_t*>(values.data());
    const size_t byte_count = values.size() * sizeof(T);
    for (size_t i = 0; i < byte_count; ++i) {
        hash ^= static_cast<uint64_t>(bytes[i]);
        hash *= kPrime;
    }

    std::ostringstream stream;
    stream << std::hex << std::setfill('0') << std::setw(16) << hash;
    return stream.str();
}

std::string fnv1a_hash(const std::string& value) {
    constexpr uint64_t kOffset = 14695981039346656037ull;
    constexpr uint64_t kPrime = 1099511628211ull;

    uint64_t hash = kOffset;
    for (unsigned char byte : value) {
        hash ^= static_cast<uint64_t>(byte);
        hash *= kPrime;
    }

    std::ostringstream stream;
    stream << std::hex << std::setfill('0') << std::setw(16) << hash;
    return stream.str();
}

template <typename T>
struct RangeStats {
    T min;
    T max;
    double mean = 0.0;
};

template <typename T>
RangeStats<T> compute_stats(const std::vector<T>& values) {
    if (values.empty()) {
        fail("cannot compute stats for empty vector");
    }

    RangeStats<T> stats{};
    stats.min = values.front();
    stats.max = values.front();
    long double total = 0.0;
    for (const T value : values) {
        stats.min = std::min(stats.min, value);
        stats.max = std::max(stats.max, value);
        total += static_cast<long double>(value);
    }
    stats.mean = static_cast<double>(total / static_cast<long double>(values.size()));
    return stats;
}

std::vector<uint16_t> to_metric_heightmap(const tectonics::contract::Snapshot& snapshot) {
    std::vector<uint16_t> values(snapshot.cell_count());
    for (size_t i = 0; i < snapshot.cell_count(); ++i) {
        values[i] = TopographyCodec::internal_to_meters(snapshot.heightmap[i], snapshot.sea_level_m);
    }
    return values;
}

std::array<uint32_t, 8> regime_counts(const tectonics::contract::Snapshot& snapshot) {
    std::array<uint32_t, 8> counts{};
    for (uint8_t regime : snapshot.geologic_regime) {
        if (regime < counts.size()) {
            ++counts[regime];
        }
    }
    return counts;
}

std::array<uint32_t, 4> crust_class_counts(const tectonics::contract::Snapshot& snapshot) {
    std::array<uint32_t, 4> counts{};
    for (uint8_t crust_class : snapshot.crust_class) {
        if (crust_class < counts.size()) {
            ++counts[crust_class];
        }
    }
    return counts;
}

std::vector<uint32_t> plate_cell_counts(const tectonics::contract::Snapshot& snapshot) {
    std::vector<uint32_t> counts(snapshot.plates.size(), 0);
    for (uint32_t plate_id : snapshot.plate_id) {
        if (plate_id < counts.size()) {
            ++counts[plate_id];
        }
    }
    return counts;
}

std::string build_boundary_segments_json(const tectonics::contract::Snapshot& snapshot) {
    std::ostringstream output;
    output << std::setprecision(17);
    output << "[\n";
    for (size_t i = 0; i < snapshot.boundary_segments.size(); ++i) {
        const auto& segment = snapshot.boundary_segments[i];
        output << "  {\n";
        output << "    \"id\": " << segment.id << ",\n";
        output << "    \"left_plate_id\": " << segment.left_plate_id << ",\n";
        output << "    \"right_plate_id\": " << segment.right_plate_id << ",\n";
        output << "    \"cell_count\": " << segment.cell_count << ",\n";
        output << "    \"persistence_steps\": " << segment.persistence_steps << ",\n";
        output << "    \"centroid\": [" << segment.centroid_x << ", " << segment.centroid_y
               << "],\n";
        output << "    \"length_cells\": " << segment.length_cells << ",\n";
        output << "    \"average_normal_motion\": " << segment.average_normal_motion << ",\n";
        output << "    \"average_shear_motion\": " << segment.average_shear_motion << ",\n";
        output << "    \"average_convergence_score\": "
               << static_cast<uint32_t>(segment.average_convergence_score) << ",\n";
        output << "    \"average_divergence_score\": "
               << static_cast<uint32_t>(segment.average_divergence_score) << ",\n";
        output << "    \"average_shear_score\": "
               << static_cast<uint32_t>(segment.average_shear_score) << ",\n";
        output << "    \"boundary_type\": \""
               << tectonics::contract::boundary_type_name(segment.boundary_type) << "\",\n";
        output << "    \"geologic_regime\": \""
               << tectonics::contract::geologic_regime_name(segment.geologic_regime) << "\",\n";
        output << "    \"age_myr\": " << segment.age_myr << "\n";
        output << "  }" << (i + 1 < snapshot.boundary_segments.size() ? "," : "") << "\n";
    }
    output << "]\n";
    return output.str();
}

std::string build_deforming_regions_json(const tectonics::contract::Snapshot& snapshot) {
    std::ostringstream output;
    output << std::setprecision(17);
    output << "[\n";
    for (size_t i = 0; i < snapshot.deforming_regions.size(); ++i) {
        const auto& region = snapshot.deforming_regions[i];
        output << "  {\n";
        output << "    \"id\": " << region.id << ",\n";
        output << "    \"boundary_segment_id\": " << region.boundary_segment_id << ",\n";
        output << "    \"primary_plate_id\": " << region.primary_plate_id << ",\n";
        output << "    \"secondary_plate_id\": " << region.secondary_plate_id << ",\n";
        output << "    \"cell_count\": " << region.cell_count << ",\n";
        output << "    \"persistence_steps\": " << region.persistence_steps << ",\n";
        output << "    \"centroid\": [" << region.centroid_x << ", " << region.centroid_y
               << "],\n";
        output << "    \"average_deformation_rate\": " << region.average_deformation_rate
               << ",\n";
        output << "    \"average_interpolated_velocity\": ["
               << region.average_interpolated_velocity_x << ", "
               << region.average_interpolated_velocity_y << "],\n";
        output << "    \"average_normal_motion\": " << region.average_normal_motion << ",\n";
        output << "    \"average_shear_motion\": " << region.average_shear_motion << ",\n";
        output << "    \"type\": \""
               << tectonics::contract::deforming_region_type_name(region.type) << "\",\n";
        output << "    \"age_myr\": " << region.age_myr << "\n";
        output << "  }" << (i + 1 < snapshot.deforming_regions.size() ? "," : "") << "\n";
    }
    output << "]\n";
    return output.str();
}

std::string build_tectonic_junctions_json(const tectonics::contract::Snapshot& snapshot) {
    std::ostringstream output;
    output << std::setprecision(17);
    output << "[\n";
    for (size_t i = 0; i < snapshot.tectonic_junctions.size(); ++i) {
        const auto& junction = snapshot.tectonic_junctions[i];
        output << "  {\n";
        output << "    \"id\": " << junction.id << ",\n";
        output << "    \"plate_ids\": [" << junction.plate_ids[0] << ", "
               << junction.plate_ids[1] << ", " << junction.plate_ids[2] << "],\n";
        output << "    \"plate_pairs\": [[" << junction.plate_ids[0] << ", "
               << junction.plate_ids[1] << "], [" << junction.plate_ids[0] << ", "
               << junction.plate_ids[2] << "], [" << junction.plate_ids[1] << ", "
               << junction.plate_ids[2] << "]],\n";
        output << "    \"boundary_segment_ids\": [" << junction.boundary_segment_ids[0]
               << ", " << junction.boundary_segment_ids[1] << ", "
               << junction.boundary_segment_ids[2] << "],\n";
        output << "    \"arm_types\": [\""
               << tectonics::contract::junction_arm_type_name(junction.arm_types[0]) << "\", \""
               << tectonics::contract::junction_arm_type_name(junction.arm_types[1]) << "\", \""
               << tectonics::contract::junction_arm_type_name(junction.arm_types[2])
               << "\"],\n";
        output << "    \"arm_strikes\": [[" << junction.arm_strike_x[0] << ", "
               << junction.arm_strike_y[0] << "], [" << junction.arm_strike_x[1] << ", "
               << junction.arm_strike_y[1] << "], [" << junction.arm_strike_x[2] << ", "
               << junction.arm_strike_y[2] << "]],\n";
        output << "    \"arm_orientation_quality\": ["
               << junction.arm_orientation_quality[0] << ", "
               << junction.arm_orientation_quality[1] << ", "
               << junction.arm_orientation_quality[2] << "],\n";
        output << "    \"configuration\": \""
               << tectonics::contract::junction_arm_type_code(junction.arm_types[0]) << "-"
               << tectonics::contract::junction_arm_type_code(junction.arm_types[1]) << "-"
               << tectonics::contract::junction_arm_type_code(junction.arm_types[2])
               << "\",\n";
        output << "    \"candidate_cell_count\": " << junction.candidate_cell_count << ",\n";
        output << "    \"centroid\": [" << junction.centroid_x << ", " << junction.centroid_y
               << "],\n";
        output << "    \"junction_velocity\": [" << junction.junction_velocity_x << ", "
               << junction.junction_velocity_y << "],\n";
        output << "    \"velocity_closure_error\": " << junction.velocity_closure_error << ",\n";
        output << "    \"velocity_quality_score\": " << junction.velocity_quality_score << ",\n";
        output << "    \"velocity_constraint_error\": "
               << junction.velocity_constraint_error << ",\n";
        output << "    \"geometry_quality_score\": "
               << junction.geometry_quality_score << ",\n";
        output << "    \"stability\": \""
               << tectonics::contract::junction_stability_name(junction.stability) << "\",\n";
        output << "    \"stability_score\": " << junction.stability_score << ",\n";
        output << "    \"quality_score\": " << junction.quality_score << "\n";
        output << "  }" << (i + 1 < snapshot.tectonic_junctions.size() ? "," : "") << "\n";
    }
    output << "]\n";
    return output.str();
}

void append_scenario_json(std::ostream& output, const tectonics::scenario::Scenario& scenario,
                          const std::string& indent) {
    output << indent << "{\n";
    output << indent << "  \"version\": " << scenario.version << ",\n";
    output << indent << "  \"seed\": " << scenario.seed << ",\n";
    output << indent << "  \"width\": " << scenario.width << ",\n";
    output << indent << "  \"height\": " << scenario.height << ",\n";
    output << indent << "  \"sea_level\": " << scenario.sea_level << ",\n";
    output << indent << "  \"plate_count\": " << scenario.plate_count << ",\n";
    output << indent << "  \"erosion_strength\": " << scenario.erosion_strength << ",\n";
    output << indent << "  \"sea_level_m\": " << scenario.sea_level_m << ",\n";
    output << indent << "  \"initial_ocean_depth_m\": " << scenario.initial_ocean_depth_m
           << ",\n";
    output << indent << "  \"initial_land_relief_m\": " << scenario.initial_land_relief_m
           << ",\n";
    output << indent << "  \"evolution_duration_myr\": " << scenario.evolution_duration_myr << ",\n";
    output << indent << "  \"motion_model\": \"" << "force-balance" << "\",\n";
    output << indent << "  \"surface_model\": \"" << "isostatic-columns" << "\",\n";
    output << indent << "  \"transport_subdivision\": " << scenario.transport_subdivision << ",\n";
    output << indent << "  \"deformation_peak_km\": " << scenario.deformation_peak_km << ",\n";
    output << indent << "  \"deformation_variation\": " << scenario.deformation_variation << ",\n";
    output << indent << "  \"subduction_coupling\": " << scenario.subduction_coupling << ",\n";
    output << indent << "  \"coastal_gpe_strength\": " << scenario.coastal_gpe_strength << ",\n";
    output << indent << "  \"boundary_force_ramp_cm_yr\": " << scenario.boundary_force_ramp_cm_yr << ",\n";
    output << indent << "  \"mountain_force_range_km\": " << scenario.mountain_force_range_km << ",\n";
    output << indent << "  \"initial_margin_width_m\": " << scenario.initial_margin_width_m << ",\n";
    output << indent << "  \"initial_boundary_influence\": " << scenario.initial_boundary_influence << ",\n";
    output << indent << "  \"consistent_motion\": " << (scenario.consistent_motion?"true":"false") << ",\n";
    output << indent << "  \"boundary_memory_myr\": " << scenario.boundary_memory_myr << ",\n";
    output << indent << "  \"noise_pass_amplitude_m\": " << scenario.experimental_noise.amplitude_m << ",\n";
    output << indent << "  \"noise_pass_wavelength_km\": " << scenario.experimental_noise.wavelength_km << ",\n";
    output << indent << "  \"noise_pass_period_myr\": " << scenario.experimental_noise.period_myr << ",\n";
    output << indent << "  \"noise_pass_coast_taper_m\": " << scenario.experimental_noise.coast_taper_m << ",\n";
    output << indent << "  \"ocean_crust_variation_m\": " << scenario.ocean_crust_variation_m << ",\n";
    output << indent << "  \"ocean_noise_wavelength_km\": " << scenario.ocean_noise_wavelength_km << ",\n";
    output << indent << "  \"seamount_crust_m\": " << scenario.seamount_crust_m << ",\n";
    output << indent << "  \"seamount_radius_km\": " << scenario.seamount_radius_km << ",\n";
    output << indent << "  \"seamount_spacing_km\": " << scenario.seamount_spacing_km << ",\n";
    output << indent << "  \"seamount_roughness\": " << scenario.seamount_roughness << ",\n";
    output << indent << "  \"submarine_erosion_ratio\": " << scenario.submarine_erosion_ratio << ",\n";
    output << indent << "  \"margin_processes\": " << (scenario.margin_processes.enabled ? "true" : "false") << ",\n";
    output << indent << "  \"boundary_conditioned_terrain\": " << (scenario.boundary_conditioned_terrain ? "true" : "false") << ",\n";
    output << indent << "  \"sediment_load_feedback\": " << (scenario.sediment_load_feedback ? "true" : "false") << ",\n";
    output << indent << "  \"trench_load_n_m\": " << scenario.margin_processes.trench_load_n_m << ",\n";
    output << indent << "  \"flexural_length_km\": " << scenario.margin_processes.flexural_length_km << ",\n";
    output << indent << "  \"trench_response_myr\": " << scenario.margin_processes.trench_response_myr << ",\n";
    output << indent << "  \"process_variation\": " << scenario.margin_processes.process_variation << ",\n";
    output << indent << "  \"process_wavelength_km\": " << scenario.margin_processes.process_wavelength_km << ",\n";
    output << indent << "  \"process_period_myr\": " << scenario.margin_processes.process_period_myr << ",\n";
    output << indent << "  \"arc_distance_km\": " << scenario.margin_processes.arc_distance_km << ",\n";
    output << indent << "  \"arc_melt_ratio\": " << scenario.margin_processes.arc_melt_ratio << ",\n";
    output << indent << "  \"river_erosion\": " << scenario.margin_processes.river_erosion << ",\n";
    output << indent << "  \"sediment_subduction_fraction\": " << scenario.margin_processes.sediment_subduction_fraction << ",\n";
    output << indent << "  \"progressive_rifting\": " << scenario.topology_processes.progressive_rifting << ",\n";
    output << indent << "  \"rift_yield_n_m\": " << scenario.topology_processes.yield_n_m << ",\n";
    output << indent << "  \"rift_damage_time_myr\": " << scenario.topology_processes.damage_time_myr << ",\n";
    output << indent << "  \"rift_strain_rate_myr\": " << scenario.topology_processes.strain_rate_myr << ",\n";
    output << indent << "  \"thermal_age_rate\": " << scenario.margin_processes.thermal_age_rate << ",\n";
    output << indent << "  \"thermal_force_feedback\": " << (scenario.thermal_force_feedback?"true":"false") << ",\n";
    output << indent << "  \"wave_erosion_m_myr\": " << scenario.margin_processes.wave_erosion_m_myr << ",\n";
    output << indent << "  \"wave_height_scale_m\": " << scenario.margin_processes.wave_height_scale_m << ",\n";
    output << indent << "  \"alluvial_transport\": " << scenario.margin_processes.alluvial_transport << ",\n";
    output << indent << "  \"river_outlet_sea_level\": " << (scenario.margin_processes.river_outlet_sea_level?"true":"false") << ",\n";
    output << indent << "  \"wave_retreat_km_myr\": " << scenario.margin_processes.wave_retreat_km_myr << ",\n";
    output << indent << "  \"root_mobility\": " << scenario.margin_processes.root_mobility << ",\n";
    output << indent << "  \"sediment_travel_km\": " << scenario.margin_processes.sediment_travel_km << ",\n";
    output << indent << "  \"rift_breakup_thickness_m\": " << scenario.margin_processes.rift_breakup_thickness_m << ",\n";
    output << indent << "  \"world_width_km\": " << scenario.world_width_km << ",\n";
    output << indent << "  \"mantle_flow_speed_cm_yr\": " << scenario.mantle_flow_speed_cm_yr << ",\n";
    output << indent << "  \"thermal_subsidence_strength\": " << scenario.thermal_subsidence_strength << ",\n";
    output << indent << "  \"delta_time_myr\": " << scenario.delta_time_myr << "\n";
    output << indent << "}";
}

void append_summary_json(std::ostream& output, const SnapshotSummary& summary,
                         const std::string& indent) {
    output << indent << "{\n";
    output << indent << "  \"contract_schema_version\": " << summary.contract_schema_version
           << ",\n";
    output << indent << "  \"width\": " << summary.width << ",\n";
    output << indent << "  \"height\": " << summary.height << ",\n";
    output << indent << "  \"iteration_count\": " << summary.iteration_count << ",\n";
    output << indent << "  \"cycle_count\": " << summary.cycle_count << ",\n";
    output << indent << "  \"time_origin_step\": " << summary.time_origin_step << ",\n";
    output << indent << "  \"sea_level_m\": " << summary.sea_level_m << ",\n";
    output << indent << "  \"time_myr\": " << summary.time_myr << ",\n";
    output << indent << "  \"delta_time_myr\": " << summary.delta_time_myr << ",\n";
    output << indent << "  \"plate_count\": " << summary.plate_count << ",\n";
    output << indent << "  \"tectonic_junction_count\": "
           << summary.tectonic_junction_count << ",\n";
    output << indent << "  \"tectonic_junction_stability_counts\": {\n";
    for (size_t i = 0; i < summary.tectonic_junction_stability_counts.size(); ++i) {
        output << indent << "    \""
               << tectonics::contract::junction_stability_name(
                      static_cast<tectonics::contract::JunctionStability>(i))
               << "\": " << summary.tectonic_junction_stability_counts[i]
               << (i + 1 < summary.tectonic_junction_stability_counts.size() ? ",\n" : "\n");
    }
    output << indent << "  },\n";
    output << indent << "  \"active_convergence_cells\": " << summary.active_convergence_cells
           << ",\n";
    output << indent << "  \"active_divergence_cells\": " << summary.active_divergence_cells
           << ",\n";
    output << indent << "  \"active_shear_cells\": " << summary.active_shear_cells << ",\n";
    output << indent << "  \"regime_counts\": {\n";
    for (size_t i = 0; i < summary.regime_counts.size(); ++i) {
        output << indent << "    \""
               << tectonics::contract::geologic_regime_name(
                      static_cast<tectonics::contract::GeologicRegime>(i))
               << "\": " << summary.regime_counts[i]
               << (i + 1 < summary.regime_counts.size() ? ",\n" : "\n");
    }
    output << indent << "  },\n";
    output << indent << "  \"crust_class_counts\": {\n";
    for (size_t i = 0; i < summary.crust_class_counts.size(); ++i) {
        output << indent << "    \""
               << tectonics::contract::crust_class_name(
                      static_cast<tectonics::contract::CrustClass>(i))
               << "\": " << summary.crust_class_counts[i]
               << (i + 1 < summary.crust_class_counts.size() ? ",\n" : "\n");
    }
    output << indent << "  },\n";
    output << indent << "  \"boundary\": {\n";
    output << indent << "    \"segment_count\": " << summary.boundary.segment_count << ",\n";
    output << indent << "    \"convergent_segment_count\": "
           << summary.boundary.convergent_segment_count << ",\n";
    output << indent << "    \"divergent_segment_count\": "
           << summary.boundary.divergent_segment_count << ",\n";
    output << indent << "    \"transform_segment_count\": "
           << summary.boundary.transform_segment_count << ",\n";
    output << indent << "    \"passive_margin_segment_count\": "
           << summary.boundary.passive_margin_segment_count << ",\n";
    output << indent << "    \"active_cell_count\": " << summary.boundary.active_cell_count
           << ",\n";
    output << indent << "    \"convergent_cell_count\": "
           << summary.boundary.convergent_cell_count << ",\n";
    output << indent << "    \"divergent_cell_count\": "
           << summary.boundary.divergent_cell_count << ",\n";
    output << indent << "    \"transform_cell_count\": "
           << summary.boundary.transform_cell_count << ",\n";
    output << indent << "    \"passive_margin_cell_count\": "
           << summary.boundary.passive_margin_cell_count << ",\n";
    output << indent << "    \"nearest_mapped_cell_count\": "
           << summary.boundary.nearest_mapped_cell_count << ",\n";
    output << indent << "    \"max_segment_age_myr\": "
           << summary.boundary.max_segment_age_myr << ",\n";
    output << indent << "    \"mean_segment_age_myr\": "
           << summary.boundary.mean_segment_age_myr << ",\n";
    output << indent << "    \"mean_segment_persistence_steps\": "
           << summary.boundary.mean_segment_persistence_steps << "\n";
    output << indent << "  },\n";
    output << indent << "  \"deforming\": {\n";
    output << indent << "    \"region_count\": " << summary.deforming.region_count << ",\n";
    output << indent << "    \"continental_rift_region_count\": "
           << summary.deforming.continental_rift_region_count << ",\n";
    output << indent << "    \"diffuse_collision_region_count\": "
           << summary.deforming.diffuse_collision_region_count << ",\n";
    output << indent << "    \"active_cell_count\": " << summary.deforming.active_cell_count
           << ",\n";
    output << indent << "    \"continental_rift_cell_count\": "
           << summary.deforming.continental_rift_cell_count << ",\n";
    output << indent << "    \"diffuse_collision_cell_count\": "
           << summary.deforming.diffuse_collision_cell_count << ",\n";
    output << indent << "    \"max_region_age_myr\": "
           << summary.deforming.max_region_age_myr << ",\n";
    output << indent << "    \"mean_region_age_myr\": "
           << summary.deforming.mean_region_age_myr << ",\n";
    output << indent << "    \"mean_region_persistence_steps\": "
           << summary.deforming.mean_region_persistence_steps << ",\n";
    output << indent << "    \"mean_deformation_rate\": "
           << summary.deforming.mean_deformation_rate << "\n";
    output << indent << "  },\n";
    output << indent << "  \"hashes\": {\n";
    output << indent << "    \"heightmap\": \"" << summary.heightmap_hash << "\",\n";
    output << indent << "    \"plate_id\": \"" << summary.plate_id_hash << "\",\n";
    output << indent << "    \"geologic_regime\": \"" << summary.geologic_regime_hash
           << "\",\n";
    output << indent << "    \"boundary_segment_id\": \""
           << summary.boundary_segment_id_hash << "\",\n";
    output << indent << "    \"deforming_region_id\": \""
           << summary.deforming_region_id_hash << "\",\n";
    output << indent << "    \"tectonic_junctions\": \"" << summary.tectonic_junctions_hash
           << "\"\n";
    output << indent << "  }\n";
    output << indent << "}";
}

void write_bundle_manifest(const fs::path& path, const tectonics::contract::Snapshot& snapshot,
                           const SnapshotBundleOptions& options) {
    const SnapshotSummary summary = summarize_snapshot(snapshot);
    const std::vector<uint16_t> metric_heightmap = to_metric_heightmap(snapshot);
    const RangeStats<float> height_stats = compute_stats(snapshot.heightmap);
    const RangeStats<uint16_t> metric_stats = compute_stats(metric_heightmap);
    const RangeStats<uint32_t> age_steps_stats = compute_stats(snapshot.crust_age_steps);
    const RangeStats<float> age_myr_stats = compute_stats(snapshot.crust_age_myr);
    const std::vector<uint32_t> cell_counts = plate_cell_counts(snapshot);
    const std::string boundary_segments_json = build_boundary_segments_json(snapshot);
    const std::string deforming_regions_json = build_deforming_regions_json(snapshot);
    const std::string tectonic_junctions_json = build_tectonic_junctions_json(snapshot);

    std::ofstream output(path);
    if (!output) {
        fail("failed to open " + path.string());
    }

    output << std::setprecision(17);
    output << "{\n";
    output << "  \"schema\": \"" << json_escape(options.schema) << "\",\n";
    output << "  \"label\": \"" << json_escape(options.label) << "\",\n";
    output << "  \"source\": \"" << json_escape(options.source) << "\",\n";
    output << "  \"contract_schema_version\": " << snapshot.schema_version << ",\n";
    if (options.include_update_count) {
        output << "  \"update_count\": " << options.update_count << ",\n";
    }
    output << "  \"seed\": " << snapshot.run_scenario.seed << ",\n";
    output << "  \"width\": " << snapshot.width << ",\n";
    output << "  \"height\": " << snapshot.height << ",\n";
    output << "  \"iteration_count\": " << snapshot.iteration_count << ",\n";
    output << "  \"cycle_count\": " << snapshot.cycle_count << ",\n";
    output << "  \"time_origin_step\": " << snapshot.time_origin_step << ",\n";
    output << "  \"time_myr\": " << snapshot.time_myr << ",\n";
    output << "  \"delta_time_myr\": " << snapshot.delta_time_myr << ",\n";
    output << "  \"sea_level_m\": " << snapshot.sea_level_m << ",\n";
    output << "  \"run_scenario\": ";
    append_scenario_json(output, snapshot.run_scenario, "  ");
    const auto& motion = snapshot.motion_diagnostics;
    output << ",\n  \"motion_diagnostics\": {\n";
    output << "    \"enabled\": " << (motion.enabled ? "true" : "false") << ",\n";
    output << "    \"relative_force_residual\": " << motion.relative_force_residual << ",\n";
    output << "    \"nonlinear_iterations\": " << motion.nonlinear_iterations << ",\n";
    output << "    \"velocity_consistency_cm_yr\": " << motion.velocity_consistency_cm_yr << ",\n";
    output << "    \"max_velocity_consistency_cm_yr\": " << motion.max_velocity_consistency_cm_yr << ",\n";
    output << "    \"persistent_face_steps\": " << motion.persistent_face_steps << ",\n";
    output << "    \"reversing_face_steps\": " << motion.reversing_face_steps << ",\n";
    output << "    \"strong_face_steps\": " << motion.strong_face_steps << ",\n";
    output << "    \"strong_reversing_face_steps\": " << motion.strong_reversing_face_steps << ",\n";
    output << "    \"driving_power_w\": " << motion.driving_power_w << ",\n";
    output << "    \"dissipation_w\": " << motion.dissipation_w << ",\n";
    output << "    \"max_speed_cm_yr\": " << motion.max_speed_cm_yr << ",\n";
    output << "    \"mean_speed_cm_yr\": " << motion.mean_speed_cm_yr << ",\n";
    output << "    \"ridge_length_km\": " << motion.ridge_length_km << ",\n";
    output << "    \"subduction_length_km\": " << motion.subduction_length_km << ",\n";
    output << "    \"continental_contact_length_km\": " << motion.continental_contact_length_km << ",\n";
    output << "    \"mountain_boundary_force_n\": " << motion.mountain_boundary_force_n << "\n  }";
    const auto& budget = snapshot.column_budget;
    output << ",\n  \"column_budget\": {\n";
    output << "    \"enabled\": " << (budget.enabled ? "true" : "false") << ",\n";
    output << "    \"initial_mass_kg\": " << budget.initial_mass_kg << ",\n";
    output << "    \"mass_kg\": " << budget.mass_kg << ",\n";
    output << "    \"generated_ocean_kg\": " << budget.generated_ocean_kg << ",\n";
    output << "    \"generated_arc_kg\": " << budget.generated_arc_kg << ",\n";
    output << "    \"river_eroded_kg\": " << budget.river_eroded_kg << ",\n";
    output << "    \"subducted_ocean_kg\": " << budget.subducted_ocean_kg << ",\n";
    output << "    \"subducted_sediment_kg\": " << budget.subducted_sediment_kg << ",\n";
    output << "    \"eroded_kg\": " << budget.eroded_kg << ",\n";
    output << "    \"deposited_kg\": " << budget.deposited_kg << ",\n";
    output << "    \"redistributed_continent_kg\": " << budget.redistributed_continent_kg << ",\n";
    output << "    \"redistributed_sediment_kg\": " << budget.redistributed_sediment_kg << ",\n";
    output << "    \"noise_redistributed_kg\": " << budget.noise_redistributed_kg << ",\n";
    output << "    \"relative_mass_residual\": " << budget.relative_mass_residual << "\n  }";
    const auto& topology = snapshot.topology_diagnostics;
    if (topology.enabled) {
        output << ",\n  \"topology_diagnostics\": {\n    \"enabled\": true,\n";
        output << "    \"reassigned_cells\": " << topology.reassigned_cells << ",\n";
        output << "    \"max_detached_components\": " << topology.max_detached_components << ",\n";
        output << "    \"event_model\": \"traction-rift-welding-v1\",\n";
        output << "    \"max_rift_damage\": " << topology.max_rift_damage << ",\n";
        output << "    \"max_tension_n_m\": " << topology.max_tension_n_m << ",\n";
        output << "    \"damaged_rift_cells\": " << topology.damaged_rift_cells << ",\n";
        output << "    \"max_weld_time_myr\": " << topology.max_weld_time_myr << ",\n";
        output << "    \"stable_ids\": [";
        for (size_t i = 0; i < topology.stable_ids.size(); ++i)
            output << (i ? ", " : "") << topology.stable_ids[i];
        output << "],\n    \"events\": [";
        for (size_t i = 0; i < topology.events.size(); ++i) {
            const auto& e = topology.events[i];
            output << (i ? ",\n" : "\n") << "      {\"kind\": \""
                << (e.kind == tectonics::topology::EventKind::Split ? "split" : "fusion")
                << "\", \"time_myr\": " << e.time_myr << ", \"parent_a\": " << e.parent_a
                << ", \"parent_b\": " << e.parent_b << ", \"result\": " << e.result
                << ", \"cells\": " << e.cells << "}";
        }
        output << "]\n  }";
    }
    output << ",\n";
    output << "  \"files\": {\n";
    if (budget.enabled) {
        output << "    \"continental_crust_m\": \"continental_crust_m.f32\",\n";
        output << "    \"oceanic_crust_m\": \"oceanic_crust_m.f32\",\n";
        output << "    \"sediment_m\": \"sediment_m.f32\",\n";
        output << "    \"rift_damage\": \"rift_damage.f32\",\n";
    }
    output << "    \"heightmap\": \"heightmap.f32\",\n";
    output << "    \"plate_id\": \"plate_id.u32\",\n";
    output << "    \"crust_age_steps\": \"crust_age_steps.u32\",\n";
    output << "    \"crust_age_myr\": \"crust_age_myr.f32\",\n";
    output << "    \"crust_thickness\": \"crust_thickness.f32\",\n";
    output << "    \"crust_class\": \"crust_class.u8\",\n";
    output << "    \"uplift_tendency\": \"uplift_tendency.f32\",\n";
    output << "    \"subsidence_tendency\": \"subsidence_tendency.f32\",\n";
    output << "    \"accumulated_strain\": \"accumulated_strain.f32\",\n";
    output << "    \"boundary_type\": \"boundary_type.u8\",\n";
    output << "    \"boundary_distance\": \"boundary_distance.u16\",\n";
    output << "    \"boundary_segment_id\": \"boundary_segment_id.u32\",\n";
    output << "    \"nearest_boundary_id\": \"nearest_boundary_id.u32\",\n";
    output << "    \"deforming_region_id\": \"deforming_region_id.u32\",\n";
    output << "    \"deforming_region_type\": \"deforming_region_type.u8\",\n";
    output << "    \"deformation_rate\": \"deformation_rate.f32\",\n";
    output << "    \"deformation_velocity_x\": \"deformation_velocity_x.f32\",\n";
    output << "    \"deformation_velocity_y\": \"deformation_velocity_y.f32\",\n";
    output << "    \"convergence_score\": \"convergence_score.u8\",\n";
    output << "    \"divergence_score\": \"divergence_score.u8\",\n";
    output << "    \"shear_score\": \"shear_score.u8\",\n";
    output << "    \"geologic_regime\": \"geologic_regime.u8\",\n";
    output << "    \"boundary_segments\": \"boundary_segments.json\",\n";
    output << "    \"deforming_regions\": \"deforming_regions.json\",\n";
    output << "    \"tectonic_junctions\": \"tectonic_junctions.json\"\n";
    output << "  },\n";
    output << "  \"hashes\": {\n";
    if (budget.enabled) {
        output << "    \"continental_crust_m\": \"" << fnv1a_hash(snapshot.continental_crust_m) << "\",\n";
        output << "    \"oceanic_crust_m\": \"" << fnv1a_hash(snapshot.oceanic_crust_m) << "\",\n";
        output << "    \"sediment_m\": \"" << fnv1a_hash(snapshot.sediment_m) << "\",\n";
        output << "    \"rift_damage\": \"" << fnv1a_hash(snapshot.rift_damage) << "\",\n";
    }
    output << "    \"heightmap\": \"" << fnv1a_hash(snapshot.heightmap) << "\",\n";
    output << "    \"plate_id\": \"" << fnv1a_hash(snapshot.plate_id) << "\",\n";
    output << "    \"crust_age_steps\": \"" << fnv1a_hash(snapshot.crust_age_steps)
           << "\",\n";
    output << "    \"crust_age_myr\": \"" << fnv1a_hash(snapshot.crust_age_myr) << "\",\n";
    output << "    \"crust_thickness\": \"" << fnv1a_hash(snapshot.crust_thickness)
           << "\",\n";
    output << "    \"crust_class\": \"" << fnv1a_hash(snapshot.crust_class) << "\",\n";
    output << "    \"uplift_tendency\": \"" << fnv1a_hash(snapshot.uplift_tendency)
           << "\",\n";
    output << "    \"subsidence_tendency\": \""
           << fnv1a_hash(snapshot.subsidence_tendency) << "\",\n";
    output << "    \"accumulated_strain\": \"" << fnv1a_hash(snapshot.accumulated_strain)
           << "\",\n";
    output << "    \"boundary_type\": \"" << fnv1a_hash(snapshot.boundary_type)
           << "\",\n";
    output << "    \"boundary_distance\": \"" << fnv1a_hash(snapshot.boundary_distance)
           << "\",\n";
    output << "    \"boundary_segment_id\": \"" << fnv1a_hash(snapshot.boundary_segment_id)
           << "\",\n";
    output << "    \"nearest_boundary_id\": \"" << fnv1a_hash(snapshot.nearest_boundary_id)
           << "\",\n";
    output << "    \"deforming_region_id\": \"" << fnv1a_hash(snapshot.deforming_region_id)
           << "\",\n";
    output << "    \"deforming_region_type\": \""
           << fnv1a_hash(snapshot.deforming_region_type) << "\",\n";
    output << "    \"deformation_rate\": \"" << fnv1a_hash(snapshot.deformation_rate)
           << "\",\n";
    output << "    \"deformation_velocity_x\": \""
           << fnv1a_hash(snapshot.deformation_velocity_x) << "\",\n";
    output << "    \"deformation_velocity_y\": \""
           << fnv1a_hash(snapshot.deformation_velocity_y) << "\",\n";
    output << "    \"convergence_score\": \"" << fnv1a_hash(snapshot.convergence_score)
           << "\",\n";
    output << "    \"divergence_score\": \"" << fnv1a_hash(snapshot.divergence_score)
           << "\",\n";
    output << "    \"shear_score\": \"" << fnv1a_hash(snapshot.shear_score) << "\",\n";
    output << "    \"geologic_regime\": \"" << fnv1a_hash(snapshot.geologic_regime)
           << "\",\n";
    output << "    \"boundary_segments\": \"" << fnv1a_hash(boundary_segments_json)
           << "\",\n";
    output << "    \"deforming_regions\": \"" << fnv1a_hash(deforming_regions_json)
           << "\",\n";
    output << "    \"tectonic_junctions\": \"" << fnv1a_hash(tectonic_junctions_json)
           << "\"\n";
    output << "  },\n";
    output << "  \"heightmap_stats\": {\n";
    output << "    \"internal_min\": " << height_stats.min << ",\n";
    output << "    \"internal_max\": " << height_stats.max << ",\n";
    output << "    \"internal_mean\": " << height_stats.mean << ",\n";
    output << "    \"metric_min\": " << metric_stats.min << ",\n";
    output << "    \"metric_max\": " << metric_stats.max << ",\n";
    output << "    \"metric_mean\": " << metric_stats.mean << "\n";
    output << "  },\n";
    output << "  \"crust_age_steps_stats\": {\n";
    output << "    \"min\": " << age_steps_stats.min << ",\n";
    output << "    \"max\": " << age_steps_stats.max << ",\n";
    output << "    \"mean\": " << age_steps_stats.mean << "\n";
    output << "  },\n";
    output << "  \"crust_age_myr_stats\": {\n";
    output << "    \"min\": " << age_myr_stats.min << ",\n";
    output << "    \"max\": " << age_myr_stats.max << ",\n";
    output << "    \"mean\": " << age_myr_stats.mean << "\n";
    output << "  },\n";
    output << "  \"summary\": ";
    append_summary_json(output, summary, "  ");
    output << ",\n";
    output << "  \"boundary_segment_count\": " << snapshot.boundary_segments.size() << ",\n";
    output << "  \"deforming_region_count\": " << snapshot.deforming_regions.size() << ",\n";
    output << "  \"tectonic_junction_count\": " << snapshot.tectonic_junctions.size()
           << ",\n";
    output << "  \"plate_cell_counts\": [";
    for (size_t i = 0; i < cell_counts.size(); ++i) {
        output << cell_counts[i] << (i + 1 < cell_counts.size() ? ", " : "");
    }
    output << "],\n";
    output << "  \"plates\": [\n";
    for (size_t i = 0; i < snapshot.plates.size(); ++i) {
        const auto& plate = snapshot.plates[i];
        output << "    {\n";
        output << "      \"plate_id\": " << i << ",\n";
        output << "      \"unit_velocity\": [" << plate.unit_velocity_x << ", "
               << plate.unit_velocity_y << "],\n";
        output << "      \"linear_velocity\": [" << plate.linear_velocity_x << ", "
               << plate.linear_velocity_y << "],\n";
        output << "      \"angular_velocity\": " << plate.angular_velocity << ",\n";
        output << "      \"mass_center\": [" << plate.mass_center_x << ", "
               << plate.mass_center_y << "]\n";
        output << "    }" << (i + 1 < snapshot.plates.size() ? "," : "") << "\n";
    }
    output << "  ]\n";
    output << "}\n";
    if (!output) {
        fail("failed to write " + path.string());
    }
}

size_t find_key(const std::string& text, const std::string& key) {
    const std::string needle = "\"" + key + "\"";
    const size_t key_pos = text.find(needle);
    if (key_pos == std::string::npos) {
        fail("missing key in manifest: " + key);
    }
    const size_t colon_pos = text.find(':', key_pos + needle.size());
    if (colon_pos == std::string::npos) {
        fail("invalid manifest for key: " + key);
    }
    return colon_pos + 1;
}

size_t find_optional_key_after(const std::string& text, const std::string& key,
                               size_t start_pos) {
    const std::string needle = "\"" + key + "\"";
    const size_t key_pos = text.find(needle, start_pos);
    if (key_pos == std::string::npos) {
        return std::string::npos;
    }
    const size_t colon_pos = text.find(':', key_pos + needle.size());
    if (colon_pos == std::string::npos) {
        fail("invalid manifest for key: " + key);
    }
    return colon_pos + 1;
}

void skip_json_whitespace(const std::string& text, size_t& pos) {
    while (pos < text.size() &&
           (text[pos] == ' ' || text[pos] == '\n' || text[pos] == '\r' || text[pos] == '\t')) {
        ++pos;
    }
}

std::string parse_json_string(const std::string& text, size_t pos) {
    skip_json_whitespace(text, pos);
    if (pos >= text.size() || text[pos] != '"') {
        fail("expected string in manifest");
    }

    ++pos;
    std::string value;
    while (pos < text.size()) {
        const char ch = text[pos++];
        if (ch == '"') {
            return value;
        }
        if (ch == '\\') {
            if (pos >= text.size()) {
                fail("unterminated escape in manifest");
            }
            const char escaped = text[pos++];
            switch (escaped) {
            case '\\':
            case '"':
            case '/':
                value.push_back(escaped);
                break;
            case 'n':
                value.push_back('\n');
                break;
            case 'r':
                value.push_back('\r');
                break;
            case 't':
                value.push_back('\t');
                break;
            default:
                fail("unsupported escape in manifest");
            }
            continue;
        }
        value.push_back(ch);
    }
    fail("unterminated string in manifest");
}

double parse_json_double(const std::string& text, size_t pos) {
    skip_json_whitespace(text, pos);
    char* end = nullptr;
    const double value = std::strtod(text.c_str() + pos, &end);
    if (end == text.c_str() + pos || !std::isfinite(value)) {
        fail("invalid number in manifest");
    }
    return value;
}

uint32_t parse_json_u32(const std::string& text, size_t pos) {
    skip_json_whitespace(text, pos);
    char* end = nullptr;
    const unsigned long value = std::strtoul(text.c_str() + pos, &end, 10);
    if (end == text.c_str() + pos ||
        value > static_cast<unsigned long>(std::numeric_limits<uint32_t>::max())) {
        fail("invalid integer in manifest");
    }
    return static_cast<uint32_t>(value);
}

NumericFieldDiff compare_numeric_vectors(const std::vector<float>& a, const std::vector<float>& b,
                                        double epsilon = 1e-6) {
    NumericFieldDiff diff{};
    long double total_abs = 0.0;
    long double total_sq = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        const double delta = static_cast<double>(a[i]) - static_cast<double>(b[i]);
        const double abs_delta = std::fabs(delta);
        if (abs_delta > epsilon) {
            ++diff.changed_cells;
        }
        diff.max_abs_diff = std::max(diff.max_abs_diff, abs_delta);
        total_abs += abs_delta;
        total_sq += delta * delta;
    }
    const long double denom = static_cast<long double>(a.size());
    diff.mean_abs_diff = static_cast<double>(total_abs / denom);
    diff.rmse = static_cast<double>(std::sqrt(total_sq / denom));
    return diff;
}

template <typename T>
NumericFieldDiff compare_numeric_vectors(const std::vector<T>& a, const std::vector<T>& b) {
    NumericFieldDiff diff{};
    long double total_abs = 0.0;
    long double total_sq = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        const double delta = static_cast<double>(a[i]) - static_cast<double>(b[i]);
        const double abs_delta = std::fabs(delta);
        if (abs_delta > 0.0) {
            ++diff.changed_cells;
        }
        diff.max_abs_diff = std::max(diff.max_abs_diff, abs_delta);
        total_abs += abs_delta;
        total_sq += delta * delta;
    }
    const long double denom = static_cast<long double>(a.size());
    diff.mean_abs_diff = static_cast<double>(total_abs / denom);
    diff.rmse = static_cast<double>(std::sqrt(total_sq / denom));
    return diff;
}

template <typename T>
ExactFieldDiff compare_exact_vectors(const std::vector<T>& a, const std::vector<T>& b) {
    ExactFieldDiff diff{};
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] != b[i]) {
            ++diff.changed_cells;
        }
    }
    return diff;
}

bool all_exact_zero(const ExactFieldDiff& diff) {
    return diff.changed_cells == 0;
}

bool all_numeric_zero(const NumericFieldDiff& diff) {
    return diff.changed_cells == 0 && diff.max_abs_diff == 0.0;
}

} // namespace

void write_snapshot_bundle(const fs::path& output_dir, const tectonics::contract::Snapshot& snapshot,
                           const SnapshotBundleOptions& options) {
    fs::create_directories(output_dir);

    write_binary_file(output_dir / "heightmap.f32", snapshot.heightmap);
    write_binary_file(output_dir / "plate_id.u32", snapshot.plate_id);
    write_binary_file(output_dir / "crust_age_steps.u32", snapshot.crust_age_steps);
    write_binary_file(output_dir / "crust_age_myr.f32", snapshot.crust_age_myr);
    write_binary_file(output_dir / "crust_thickness.f32", snapshot.crust_thickness);
    if (snapshot.column_budget.enabled) {
        write_binary_file(output_dir / "continental_crust_m.f32", snapshot.continental_crust_m);
        write_binary_file(output_dir / "oceanic_crust_m.f32", snapshot.oceanic_crust_m);
        write_binary_file(output_dir / "sediment_m.f32", snapshot.sediment_m);
        write_binary_file(output_dir / "rift_damage.f32", snapshot.rift_damage);
    }
    write_binary_file(output_dir / "crust_class.u8", snapshot.crust_class);
    write_binary_file(output_dir / "uplift_tendency.f32", snapshot.uplift_tendency);
    write_binary_file(output_dir / "subsidence_tendency.f32", snapshot.subsidence_tendency);
    write_binary_file(output_dir / "accumulated_strain.f32", snapshot.accumulated_strain);
    write_binary_file(output_dir / "boundary_type.u8", snapshot.boundary_type);
    write_binary_file(output_dir / "boundary_distance.u16", snapshot.boundary_distance);
    write_binary_file(output_dir / "boundary_segment_id.u32", snapshot.boundary_segment_id);
    write_binary_file(output_dir / "nearest_boundary_id.u32", snapshot.nearest_boundary_id);
    write_binary_file(output_dir / "deforming_region_id.u32", snapshot.deforming_region_id);
    write_binary_file(output_dir / "deforming_region_type.u8", snapshot.deforming_region_type);
    write_binary_file(output_dir / "deformation_rate.f32", snapshot.deformation_rate);
    write_binary_file(output_dir / "deformation_velocity_x.f32", snapshot.deformation_velocity_x);
    write_binary_file(output_dir / "deformation_velocity_y.f32", snapshot.deformation_velocity_y);
    write_binary_file(output_dir / "convergence_score.u8", snapshot.convergence_score);
    write_binary_file(output_dir / "divergence_score.u8", snapshot.divergence_score);
    write_binary_file(output_dir / "shear_score.u8", snapshot.shear_score);
    write_binary_file(output_dir / "geologic_regime.u8", snapshot.geologic_regime);
    write_text_file(output_dir / "boundary_segments.json", build_boundary_segments_json(snapshot));
    write_text_file(output_dir / "deforming_regions.json", build_deforming_regions_json(snapshot));
    write_text_file(output_dir / "tectonic_junctions.json",
                    build_tectonic_junctions_json(snapshot));
    write_bundle_manifest(output_dir / "manifest.json", snapshot, options);
}

LoadedSnapshotBundle load_snapshot_bundle(const fs::path& input_dir) {
    LoadedSnapshotBundle bundle;
    const std::string manifest = read_text_file(input_dir / "manifest.json");

    bundle.schema = parse_json_string(manifest, find_key(manifest, "schema"));
    bundle.label = parse_json_string(manifest, find_key(manifest, "label"));
    bundle.source = parse_json_string(manifest, find_key(manifest, "source"));
    bundle.snapshot.schema_version =
        parse_json_u32(manifest, find_key(manifest, "contract_schema_version"));
    bundle.snapshot.run_scenario.seed =
        static_cast<long>(parse_json_u32(manifest, find_key(manifest, "seed")));
    bundle.snapshot.width = parse_json_u32(manifest, find_key(manifest, "width"));
    bundle.snapshot.height = parse_json_u32(manifest, find_key(manifest, "height"));
    bundle.snapshot.iteration_count = parse_json_u32(manifest, find_key(manifest, "iteration_count"));
    bundle.snapshot.cycle_count = parse_json_u32(manifest, find_key(manifest, "cycle_count"));
    bundle.snapshot.time_origin_step =
        parse_json_u32(manifest, find_key(manifest, "time_origin_step"));
    bundle.snapshot.time_myr = parse_json_double(manifest, find_key(manifest, "time_myr"));
    bundle.snapshot.delta_time_myr =
        parse_json_double(manifest, find_key(manifest, "delta_time_myr"));
    bundle.snapshot.sea_level_m =
        static_cast<uint16_t>(parse_json_u32(manifest, find_key(manifest, "sea_level_m")));
    bundle.declared_boundary_segment_count =
        parse_json_u32(manifest, find_key(manifest, "boundary_segment_count"));
    const size_t deforming_count_pos = find_key(manifest, "deforming_region_count");
    bundle.declared_deforming_region_count =
        parse_json_u32(manifest, deforming_count_pos);
    const size_t junction_count_pos =
        find_optional_key_after(manifest, "tectonic_junction_count", deforming_count_pos);
    if (junction_count_pos != std::string::npos) {
        bundle.declared_tectonic_junction_count =
            parse_json_u32(manifest, junction_count_pos);
    }

    const size_t cell_count = bundle.snapshot.cell_count();
    const auto topology_start = find_optional_key_after(manifest, "topology_diagnostics", 0);
    if (topology_start != std::string::npos) {
        auto& topology = bundle.snapshot.topology_diagnostics;
        topology.enabled = true;
        topology.reassigned_cells = std::stoull(manifest.substr(
            find_optional_key_after(manifest, "reassigned_cells", topology_start)));
        topology.max_detached_components = parse_json_u32(manifest,
            find_optional_key_after(manifest, "max_detached_components", topology_start));
        const auto ids = find_optional_key_after(manifest, "stable_ids", topology_start);
        if (ids != std::string::npos) {
            const auto start = manifest.find('[', ids), end = manifest.find(']', start);
            if (start == std::string::npos || end == std::string::npos) fail("Invalid plate identities");
            auto values = manifest.substr(start + 1, end - start - 1);
            std::replace(values.begin(), values.end(), ',', ' ');
            std::istringstream stream(values);
            uint32_t id;
            while (stream >> id) topology.stable_ids.push_back(id);
            if (!stream.eof()) fail("Invalid plate identity value");
            topology.max_rift_damage = parse_json_double(manifest, find_key(manifest, "max_rift_damage"));
            if (find_optional_key_after(manifest, "max_tension_n_m", 0) != std::string::npos)
                topology.max_tension_n_m = parse_json_double(manifest, find_key(manifest, "max_tension_n_m"));
            if (find_optional_key_after(manifest, "damaged_rift_cells", 0) != std::string::npos)
                topology.damaged_rift_cells = static_cast<uint32_t>(parse_json_double(manifest, find_key(manifest, "damaged_rift_cells")));
            topology.max_weld_time_myr = parse_json_double(manifest, find_key(manifest, "max_weld_time_myr"));
            const auto events = find_optional_key_after(manifest, "events", end);
            const auto array_start = manifest.find('[', events), array_end = manifest.find(']', array_start);
            if (array_start == std::string::npos || array_end == std::string::npos) fail("Invalid plate events");
            auto pos = manifest.find('{', array_start);
            while (pos != std::string::npos && pos < array_end) {
                const auto finish = manifest.find('}', pos);
                if (finish == std::string::npos || finish > array_end) fail("Invalid plate event");
                const auto object = manifest.substr(pos, finish - pos + 1);
                const auto kind = parse_json_string(object, find_key(object, "kind"));
                if (kind != "split" && kind != "fusion") fail("Invalid plate event kind");
                topology.events.push_back({kind == "split" ? tectonics::topology::EventKind::Split : tectonics::topology::EventKind::Fusion,
                    parse_json_double(object, find_key(object, "time_myr")),
                    parse_json_u32(object, find_key(object, "parent_a")), parse_json_u32(object, find_key(object, "parent_b")),
                    parse_json_u32(object, find_key(object, "result")), parse_json_u32(object, find_key(object, "cells"))});
                pos = manifest.find('{', finish + 1);
            }
        }
    }
    if (find_optional_key_after(manifest, "continental_crust_m", 0) != std::string::npos) {
        bundle.snapshot.continental_crust_m = read_binary_file<float>(input_dir / "continental_crust_m.f32", cell_count);
        bundle.snapshot.oceanic_crust_m = read_binary_file<float>(input_dir / "oceanic_crust_m.f32", cell_count);
        bundle.snapshot.sediment_m = read_binary_file<float>(input_dir / "sediment_m.f32", cell_count);
        if (find_optional_key_after(manifest, "rift_damage", 0) != std::string::npos)
            bundle.snapshot.rift_damage = read_binary_file<float>(input_dir / "rift_damage.f32", cell_count);
        auto& budget = bundle.snapshot.column_budget;
        budget.enabled = true;
        const auto start = find_key(manifest, "column_budget");
        const auto value = [&](const char* name) {
            return parse_json_double(manifest, find_optional_key_after(manifest, name, start));
        };
        budget.initial_mass_kg = value("initial_mass_kg"); budget.mass_kg = value("mass_kg");
        budget.generated_ocean_kg = value("generated_ocean_kg");
        if (find_optional_key_after(manifest, "generated_arc_kg", start) != std::string::npos)
            budget.generated_arc_kg = value("generated_arc_kg");
        if (find_optional_key_after(manifest, "river_eroded_kg", start) != std::string::npos)
            budget.river_eroded_kg = value("river_eroded_kg");
        budget.subducted_ocean_kg = value("subducted_ocean_kg");
        if (find_optional_key_after(manifest, "subducted_sediment_kg", 0) != std::string::npos)
            budget.subducted_sediment_kg = value("subducted_sediment_kg");
        budget.eroded_kg = value("eroded_kg"); budget.deposited_kg = value("deposited_kg");
        budget.relative_mass_residual = value("relative_mass_residual");
        if (find_optional_key_after(manifest, "redistributed_continent_kg", start) != std::string::npos)
            budget.redistributed_continent_kg = value("redistributed_continent_kg");
        if (find_optional_key_after(manifest, "redistributed_sediment_kg", start) != std::string::npos)
            budget.redistributed_sediment_kg = value("redistributed_sediment_kg");
        if (find_optional_key_after(manifest, "noise_redistributed_kg", start) != std::string::npos)
            budget.noise_redistributed_kg = value("noise_redistributed_kg");
    }
    bundle.snapshot.heightmap = read_binary_file<float>(input_dir / "heightmap.f32", cell_count);
    bundle.snapshot.plate_id = read_binary_file<uint32_t>(input_dir / "plate_id.u32", cell_count);
    bundle.snapshot.crust_age_steps =
        read_binary_file<uint32_t>(input_dir / "crust_age_steps.u32", cell_count);
    bundle.snapshot.crust_age_myr =
        read_binary_file<float>(input_dir / "crust_age_myr.f32", cell_count);
    bundle.snapshot.crust_thickness =
        read_binary_file<float>(input_dir / "crust_thickness.f32", cell_count);
    bundle.snapshot.crust_class =
        read_binary_file<uint8_t>(input_dir / "crust_class.u8", cell_count);
    bundle.snapshot.uplift_tendency =
        read_binary_file<float>(input_dir / "uplift_tendency.f32", cell_count);
    bundle.snapshot.subsidence_tendency =
        read_binary_file<float>(input_dir / "subsidence_tendency.f32", cell_count);
    bundle.snapshot.accumulated_strain =
        read_binary_file<float>(input_dir / "accumulated_strain.f32", cell_count);
    bundle.snapshot.boundary_type =
        read_binary_file<uint8_t>(input_dir / "boundary_type.u8", cell_count);
    bundle.snapshot.boundary_distance =
        read_binary_file<uint16_t>(input_dir / "boundary_distance.u16", cell_count);
    bundle.snapshot.boundary_segment_id =
        read_binary_file<uint32_t>(input_dir / "boundary_segment_id.u32", cell_count);
    bundle.snapshot.nearest_boundary_id =
        read_binary_file<uint32_t>(input_dir / "nearest_boundary_id.u32", cell_count);
    bundle.snapshot.deforming_region_id =
        read_binary_file<uint32_t>(input_dir / "deforming_region_id.u32", cell_count);
    bundle.snapshot.deforming_region_type =
        read_binary_file<uint8_t>(input_dir / "deforming_region_type.u8", cell_count);
    bundle.snapshot.deformation_rate =
        read_binary_file<float>(input_dir / "deformation_rate.f32", cell_count);
    bundle.snapshot.deformation_velocity_x =
        read_binary_file<float>(input_dir / "deformation_velocity_x.f32", cell_count);
    bundle.snapshot.deformation_velocity_y =
        read_binary_file<float>(input_dir / "deformation_velocity_y.f32", cell_count);
    bundle.snapshot.convergence_score =
        read_binary_file<uint8_t>(input_dir / "convergence_score.u8", cell_count);
    bundle.snapshot.divergence_score =
        read_binary_file<uint8_t>(input_dir / "divergence_score.u8", cell_count);
    bundle.snapshot.shear_score =
        read_binary_file<uint8_t>(input_dir / "shear_score.u8", cell_count);
    bundle.snapshot.geologic_regime =
        read_binary_file<uint8_t>(input_dir / "geologic_regime.u8", cell_count);
    return bundle;
}

SnapshotSummary summarize_snapshot(const tectonics::contract::Snapshot& snapshot) {
    SnapshotSummary summary{};
    summary.contract_schema_version = snapshot.schema_version;
    summary.width = snapshot.width;
    summary.height = snapshot.height;
    summary.iteration_count = snapshot.iteration_count;
    summary.cycle_count = snapshot.cycle_count;
    summary.time_origin_step = snapshot.time_origin_step;
    summary.sea_level_m = snapshot.sea_level_m;
    summary.time_myr = snapshot.time_myr;
    summary.delta_time_myr = snapshot.delta_time_myr;
    summary.plate_count = static_cast<uint32_t>(snapshot.plates.size());
    summary.tectonic_junction_count =
        static_cast<uint32_t>(snapshot.tectonic_junctions.size());
    summary.regime_counts = regime_counts(snapshot);
    summary.crust_class_counts = crust_class_counts(snapshot);
    summary.heightmap_hash = fnv1a_hash(snapshot.heightmap);
    summary.plate_id_hash = fnv1a_hash(snapshot.plate_id);
    summary.geologic_regime_hash = fnv1a_hash(snapshot.geologic_regime);
    summary.boundary_segment_id_hash = fnv1a_hash(snapshot.boundary_segment_id);
    summary.deforming_region_id_hash = fnv1a_hash(snapshot.deforming_region_id);
    summary.tectonic_junctions_hash = fnv1a_hash(build_tectonic_junctions_json(snapshot));

    for (const tectonics::contract::TectonicJunction& junction : snapshot.tectonic_junctions) {
        const size_t index = static_cast<size_t>(junction.stability);
        if (index < summary.tectonic_junction_stability_counts.size()) {
            ++summary.tectonic_junction_stability_counts[index];
        }
    }

    for (size_t i = 0; i < snapshot.cell_count(); ++i) {
        if (snapshot.convergence_score[i] > 0U) {
            ++summary.active_convergence_cells;
        }
        if (snapshot.divergence_score[i] > 0U) {
            ++summary.active_divergence_cells;
        }
        if (snapshot.shear_score[i] > 0U) {
            ++summary.active_shear_cells;
        }

        switch (static_cast<tectonics::contract::BoundaryType>(snapshot.boundary_type[i])) {
        case tectonics::contract::BoundaryType::Convergent:
            ++summary.boundary.active_cell_count;
            ++summary.boundary.convergent_cell_count;
            break;
        case tectonics::contract::BoundaryType::Divergent:
            ++summary.boundary.active_cell_count;
            ++summary.boundary.divergent_cell_count;
            break;
        case tectonics::contract::BoundaryType::Transform:
            ++summary.boundary.active_cell_count;
            ++summary.boundary.transform_cell_count;
            break;
        case tectonics::contract::BoundaryType::PassiveMargin:
            ++summary.boundary.active_cell_count;
            ++summary.boundary.passive_margin_cell_count;
            break;
        case tectonics::contract::BoundaryType::None:
        default:
            break;
        }

        if (snapshot.nearest_boundary_id[i] != tectonics::contract::kNoBoundaryId) {
            ++summary.boundary.nearest_mapped_cell_count;
        }

        switch (static_cast<tectonics::contract::DeformingRegionType>(snapshot.deforming_region_type[i])) {
        case tectonics::contract::DeformingRegionType::ContinentalRift:
            ++summary.deforming.active_cell_count;
            ++summary.deforming.continental_rift_cell_count;
            break;
        case tectonics::contract::DeformingRegionType::DiffuseCollision:
            ++summary.deforming.active_cell_count;
            ++summary.deforming.diffuse_collision_cell_count;
            break;
        case tectonics::contract::DeformingRegionType::None:
        default:
            break;
        }
    }

    summary.boundary.segment_count = static_cast<uint32_t>(snapshot.boundary_segments.size());
    if (!snapshot.boundary_segments.empty()) {
        long double age_total = 0.0;
        long double persistence_total = 0.0;
        for (const tectonics::contract::BoundarySegment& segment : snapshot.boundary_segments) {
            age_total += segment.age_myr;
            persistence_total += segment.persistence_steps;
            summary.boundary.max_segment_age_myr =
                std::max(summary.boundary.max_segment_age_myr, segment.age_myr);
            switch (segment.boundary_type) {
            case tectonics::contract::BoundaryType::Convergent:
                ++summary.boundary.convergent_segment_count;
                break;
            case tectonics::contract::BoundaryType::Divergent:
                ++summary.boundary.divergent_segment_count;
                break;
            case tectonics::contract::BoundaryType::Transform:
                ++summary.boundary.transform_segment_count;
                break;
            case tectonics::contract::BoundaryType::PassiveMargin:
                ++summary.boundary.passive_margin_segment_count;
                break;
            case tectonics::contract::BoundaryType::None:
            default:
                break;
            }
        }
        const long double denom = static_cast<long double>(snapshot.boundary_segments.size());
        summary.boundary.mean_segment_age_myr = static_cast<double>(age_total / denom);
        summary.boundary.mean_segment_persistence_steps =
            static_cast<double>(persistence_total / denom);
    }

    summary.deforming.region_count = static_cast<uint32_t>(snapshot.deforming_regions.size());
    if (!snapshot.deforming_regions.empty()) {
        long double age_total = 0.0;
        long double persistence_total = 0.0;
        long double deformation_total = 0.0;
        for (const tectonics::contract::DeformingRegion& region : snapshot.deforming_regions) {
            age_total += region.age_myr;
            persistence_total += region.persistence_steps;
            deformation_total += region.average_deformation_rate;
            summary.deforming.max_region_age_myr =
                std::max(summary.deforming.max_region_age_myr, region.age_myr);
            switch (region.type) {
            case tectonics::contract::DeformingRegionType::ContinentalRift:
                ++summary.deforming.continental_rift_region_count;
                break;
            case tectonics::contract::DeformingRegionType::DiffuseCollision:
                ++summary.deforming.diffuse_collision_region_count;
                break;
            case tectonics::contract::DeformingRegionType::None:
            default:
                break;
            }
        }
        const long double denom = static_cast<long double>(snapshot.deforming_regions.size());
        summary.deforming.mean_region_age_myr = static_cast<double>(age_total / denom);
        summary.deforming.mean_region_persistence_steps =
            static_cast<double>(persistence_total / denom);
        summary.deforming.mean_deformation_rate =
            static_cast<double>(deformation_total / denom);
    }

    return summary;
}

SnapshotComparison compare_snapshots(const tectonics::contract::Snapshot& a,
                                     const tectonics::contract::Snapshot& b,
                                     uint32_t boundary_segment_count_a,
                                     uint32_t boundary_segment_count_b,
                                     uint32_t deforming_region_count_a,
                                     uint32_t deforming_region_count_b,
                                     uint32_t tectonic_junction_count_a,
                                     uint32_t tectonic_junction_count_b) {
    SnapshotComparison comparison{};
    comparison.dimensions_match = a.width == b.width && a.height == b.height &&
                                  a.cell_count() == b.cell_count();
    comparison.width = comparison.dimensions_match ? a.width : 0;
    comparison.height = comparison.dimensions_match ? a.height : 0;
    comparison.contract_schema_version_a = a.schema_version;
    comparison.contract_schema_version_b = b.schema_version;
    comparison.time_myr_a = a.time_myr;
    comparison.time_myr_b = b.time_myr;
    comparison.iteration_count_delta =
        static_cast<int64_t>(a.iteration_count) - static_cast<int64_t>(b.iteration_count);
    comparison.cycle_count_delta =
        static_cast<int64_t>(a.cycle_count) - static_cast<int64_t>(b.cycle_count);
    comparison.plate_count_delta =
        static_cast<int64_t>(a.plates.size()) - static_cast<int64_t>(b.plates.size());
    comparison.boundary_segment_count_delta =
        static_cast<int64_t>(boundary_segment_count_a != 0 ? boundary_segment_count_a
                                                           : a.boundary_segments.size()) -
        static_cast<int64_t>(boundary_segment_count_b != 0 ? boundary_segment_count_b
                                                           : b.boundary_segments.size());
    comparison.deforming_region_count_delta =
        static_cast<int64_t>(deforming_region_count_a != 0 ? deforming_region_count_a
                                                           : a.deforming_regions.size()) -
        static_cast<int64_t>(deforming_region_count_b != 0 ? deforming_region_count_b
                                                           : b.deforming_regions.size());
    comparison.tectonic_junction_count_delta =
        static_cast<int64_t>(tectonic_junction_count_a != 0 ? tectonic_junction_count_a
                                                            : a.tectonic_junctions.size()) -
        static_cast<int64_t>(tectonic_junction_count_b != 0 ? tectonic_junction_count_b
                                                            : b.tectonic_junctions.size());

    if (!comparison.dimensions_match) {
        return comparison;
    }

    comparison.heightmap = compare_numeric_vectors(a.heightmap, b.heightmap);
    comparison.plate_id = compare_exact_vectors(a.plate_id, b.plate_id);
    comparison.crust_age_steps = compare_exact_vectors(a.crust_age_steps, b.crust_age_steps);
    comparison.crust_age_myr = compare_numeric_vectors(a.crust_age_myr, b.crust_age_myr);
    comparison.crust_thickness = compare_numeric_vectors(a.crust_thickness, b.crust_thickness);
    comparison.physical_columns_comparable = !a.continental_crust_m.empty() && !b.continental_crust_m.empty();
    if (comparison.physical_columns_comparable) {
        comparison.continental_crust_m = compare_numeric_vectors(a.continental_crust_m, b.continental_crust_m);
        comparison.oceanic_crust_m = compare_numeric_vectors(a.oceanic_crust_m, b.oceanic_crust_m);
        comparison.sediment_m = compare_numeric_vectors(a.sediment_m, b.sediment_m);
    }
    comparison.crust_class = compare_exact_vectors(a.crust_class, b.crust_class);
    comparison.uplift_tendency = compare_numeric_vectors(a.uplift_tendency, b.uplift_tendency);
    comparison.subsidence_tendency =
        compare_numeric_vectors(a.subsidence_tendency, b.subsidence_tendency);
    comparison.accumulated_strain =
        compare_numeric_vectors(a.accumulated_strain, b.accumulated_strain);
    comparison.boundary_type = compare_exact_vectors(a.boundary_type, b.boundary_type);
    comparison.boundary_distance = compare_numeric_vectors(a.boundary_distance, b.boundary_distance);
    comparison.boundary_segment_id =
        compare_exact_vectors(a.boundary_segment_id, b.boundary_segment_id);
    comparison.nearest_boundary_id =
        compare_exact_vectors(a.nearest_boundary_id, b.nearest_boundary_id);
    comparison.deforming_region_id =
        compare_exact_vectors(a.deforming_region_id, b.deforming_region_id);
    comparison.deforming_region_type =
        compare_exact_vectors(a.deforming_region_type, b.deforming_region_type);
    comparison.deformation_rate = compare_numeric_vectors(a.deformation_rate, b.deformation_rate);
    comparison.deformation_velocity_x =
        compare_numeric_vectors(a.deformation_velocity_x, b.deformation_velocity_x);
    comparison.deformation_velocity_y =
        compare_numeric_vectors(a.deformation_velocity_y, b.deformation_velocity_y);
    comparison.convergence_score =
        compare_exact_vectors(a.convergence_score, b.convergence_score);
    comparison.divergence_score =
        compare_exact_vectors(a.divergence_score, b.divergence_score);
    comparison.shear_score = compare_exact_vectors(a.shear_score, b.shear_score);
    comparison.geologic_regime = compare_exact_vectors(a.geologic_regime, b.geologic_regime);

    comparison.exact_match =
        a.continental_crust_m == b.continental_crust_m && a.oceanic_crust_m == b.oceanic_crust_m &&
        a.sediment_m == b.sediment_m && a.rift_damage == b.rift_damage &&
        comparison.contract_schema_version_a == comparison.contract_schema_version_b &&
        comparison.time_myr_a == comparison.time_myr_b &&
        comparison.iteration_count_delta == 0 &&
        comparison.cycle_count_delta == 0 &&
        comparison.plate_count_delta == 0 &&
        comparison.boundary_segment_count_delta == 0 &&
        comparison.deforming_region_count_delta == 0 &&
        comparison.tectonic_junction_count_delta == 0 &&
        all_numeric_zero(comparison.heightmap) &&
        all_exact_zero(comparison.plate_id) &&
        all_exact_zero(comparison.crust_age_steps) &&
        all_numeric_zero(comparison.crust_age_myr) &&
        all_numeric_zero(comparison.crust_thickness) &&
        all_exact_zero(comparison.crust_class) &&
        all_numeric_zero(comparison.uplift_tendency) &&
        all_numeric_zero(comparison.subsidence_tendency) &&
        all_numeric_zero(comparison.accumulated_strain) &&
        all_exact_zero(comparison.boundary_type) &&
        all_numeric_zero(comparison.boundary_distance) &&
        all_exact_zero(comparison.boundary_segment_id) &&
        all_exact_zero(comparison.nearest_boundary_id) &&
        all_exact_zero(comparison.deforming_region_id) &&
        all_exact_zero(comparison.deforming_region_type) &&
        all_numeric_zero(comparison.deformation_rate) &&
        all_numeric_zero(comparison.deformation_velocity_x) &&
        all_numeric_zero(comparison.deformation_velocity_y) &&
        all_exact_zero(comparison.convergence_score) &&
        all_exact_zero(comparison.divergence_score) &&
        all_exact_zero(comparison.shear_score) &&
        all_exact_zero(comparison.geologic_regime);

    return comparison;
}

std::string snapshot_summary_json(const SnapshotSummary& summary) {
    std::ostringstream output;
    output << std::setprecision(17);
    output << "{\n";
    output << "  \"schema\": \"tectonic-snapshot-summary/v1\",\n";
    output << "  \"summary\": ";
    append_summary_json(output, summary, "  ");
    output << "\n}\n";
    return output.str();
}

std::string boundary_stats_json(const SnapshotSummary& summary) {
    std::ostringstream output;
    output << std::setprecision(17);
    output << "{\n";
    output << "  \"schema\": \"tectonic-boundary-stats/v1\",\n";
    output << "  \"time_myr\": " << summary.time_myr << ",\n";
    output << "  \"iteration_count\": " << summary.iteration_count << ",\n";
    output << "  \"cycle_count\": " << summary.cycle_count << ",\n";
    output << "  \"tectonic_junction_count\": " << summary.tectonic_junction_count << ",\n";
    output << "  \"tectonic_junction_stability_counts\": {\n";
    for (size_t i = 0; i < summary.tectonic_junction_stability_counts.size(); ++i) {
        output << "    \""
               << tectonics::contract::junction_stability_name(
                      static_cast<tectonics::contract::JunctionStability>(i))
               << "\": " << summary.tectonic_junction_stability_counts[i]
               << (i + 1 < summary.tectonic_junction_stability_counts.size() ? ",\n" : "\n");
    }
    output << "  },\n";
    output << "  \"boundary\": {\n";
    output << "    \"segment_count\": " << summary.boundary.segment_count << ",\n";
    output << "    \"convergent_segment_count\": "
           << summary.boundary.convergent_segment_count << ",\n";
    output << "    \"divergent_segment_count\": "
           << summary.boundary.divergent_segment_count << ",\n";
    output << "    \"transform_segment_count\": "
           << summary.boundary.transform_segment_count << ",\n";
    output << "    \"passive_margin_segment_count\": "
           << summary.boundary.passive_margin_segment_count << ",\n";
    output << "    \"active_cell_count\": " << summary.boundary.active_cell_count << ",\n";
    output << "    \"convergent_cell_count\": "
           << summary.boundary.convergent_cell_count << ",\n";
    output << "    \"divergent_cell_count\": "
           << summary.boundary.divergent_cell_count << ",\n";
    output << "    \"transform_cell_count\": "
           << summary.boundary.transform_cell_count << ",\n";
    output << "    \"passive_margin_cell_count\": "
           << summary.boundary.passive_margin_cell_count << ",\n";
    output << "    \"nearest_mapped_cell_count\": "
           << summary.boundary.nearest_mapped_cell_count << ",\n";
    output << "    \"max_segment_age_myr\": "
           << summary.boundary.max_segment_age_myr << ",\n";
    output << "    \"mean_segment_age_myr\": "
           << summary.boundary.mean_segment_age_myr << ",\n";
    output << "    \"mean_segment_persistence_steps\": "
           << summary.boundary.mean_segment_persistence_steps << "\n";
    output << "  },\n";
    output << "  \"deforming\": {\n";
    output << "    \"region_count\": " << summary.deforming.region_count << ",\n";
    output << "    \"continental_rift_region_count\": "
           << summary.deforming.continental_rift_region_count << ",\n";
    output << "    \"diffuse_collision_region_count\": "
           << summary.deforming.diffuse_collision_region_count << ",\n";
    output << "    \"active_cell_count\": " << summary.deforming.active_cell_count << ",\n";
    output << "    \"continental_rift_cell_count\": "
           << summary.deforming.continental_rift_cell_count << ",\n";
    output << "    \"diffuse_collision_cell_count\": "
           << summary.deforming.diffuse_collision_cell_count << ",\n";
    output << "    \"max_region_age_myr\": "
           << summary.deforming.max_region_age_myr << ",\n";
    output << "    \"mean_region_age_myr\": "
           << summary.deforming.mean_region_age_myr << ",\n";
    output << "    \"mean_region_persistence_steps\": "
           << summary.deforming.mean_region_persistence_steps << ",\n";
    output << "    \"mean_deformation_rate\": "
           << summary.deforming.mean_deformation_rate << "\n";
    output << "  }\n";
    output << "}\n";
    return output.str();
}

std::string snapshot_comparison_json(const SnapshotComparison& comparison) {
    const auto append_numeric = [](std::ostream& output, const NumericFieldDiff& diff,
                                   const std::string& indent) {
        output << indent << "{\n";
        output << indent << "  \"changed_cells\": " << diff.changed_cells << ",\n";
        output << indent << "  \"max_abs_diff\": " << diff.max_abs_diff << ",\n";
        output << indent << "  \"mean_abs_diff\": " << diff.mean_abs_diff << ",\n";
        output << indent << "  \"rmse\": " << diff.rmse << "\n";
        output << indent << "}";
    };
    const auto append_exact = [](std::ostream& output, const ExactFieldDiff& diff,
                                 const std::string& indent) {
        output << indent << "{\n";
        output << indent << "  \"changed_cells\": " << diff.changed_cells << "\n";
        output << indent << "}";
    };

    std::ostringstream output;
    output << std::setprecision(17);
    output << "{\n";
    output << "  \"schema\": \"tectonic-snapshot-comparison/v1\",\n";
    output << "  \"dimensions_match\": " << (comparison.dimensions_match ? "true" : "false")
           << ",\n";
    output << "  \"exact_match\": " << (comparison.exact_match ? "true" : "false") << ",\n";
    output << "  \"width\": " << comparison.width << ",\n";
    output << "  \"height\": " << comparison.height << ",\n";
    output << "  \"contract_schema_version_a\": " << comparison.contract_schema_version_a
           << ",\n";
    output << "  \"contract_schema_version_b\": " << comparison.contract_schema_version_b
           << ",\n";
    output << "  \"time_myr_a\": " << comparison.time_myr_a << ",\n";
    output << "  \"time_myr_b\": " << comparison.time_myr_b << ",\n";
    output << "  \"iteration_count_delta\": " << comparison.iteration_count_delta << ",\n";
    output << "  \"cycle_count_delta\": " << comparison.cycle_count_delta << ",\n";
    output << "  \"plate_count_delta\": " << comparison.plate_count_delta << ",\n";
    output << "  \"boundary_segment_count_delta\": "
           << comparison.boundary_segment_count_delta << ",\n";
    output << "  \"deforming_region_count_delta\": "
           << comparison.deforming_region_count_delta << ",\n";
    output << "  \"tectonic_junction_count_delta\": "
           << comparison.tectonic_junction_count_delta << ",\n";
    output << "  \"fields\": {\n";
    output << "    \"heightmap\": ";
    append_numeric(output, comparison.heightmap, "    ");
    output << ",\n    \"plate_id\": ";
    append_exact(output, comparison.plate_id, "    ");
    output << ",\n    \"crust_age_steps\": ";
    append_exact(output, comparison.crust_age_steps, "    ");
    output << ",\n    \"crust_age_myr\": ";
    append_numeric(output, comparison.crust_age_myr, "    ");
    output << ",\n    \"crust_thickness\": ";
    append_numeric(output, comparison.crust_thickness, "    ");
    output << ",\n    \"physical_columns_comparable\": " << (comparison.physical_columns_comparable ? "true" : "false");
    if (comparison.physical_columns_comparable) {
        output << ",\n    \"continental_crust_m\": "; append_numeric(output, comparison.continental_crust_m, "    ");
        output << ",\n    \"oceanic_crust_m\": "; append_numeric(output, comparison.oceanic_crust_m, "    ");
        output << ",\n    \"sediment_m\": "; append_numeric(output, comparison.sediment_m, "    ");
    }
    output << ",\n    \"crust_class\": ";
    append_exact(output, comparison.crust_class, "    ");
    output << ",\n    \"uplift_tendency\": ";
    append_numeric(output, comparison.uplift_tendency, "    ");
    output << ",\n    \"subsidence_tendency\": ";
    append_numeric(output, comparison.subsidence_tendency, "    ");
    output << ",\n    \"accumulated_strain\": ";
    append_numeric(output, comparison.accumulated_strain, "    ");
    output << ",\n    \"boundary_type\": ";
    append_exact(output, comparison.boundary_type, "    ");
    output << ",\n    \"boundary_distance\": ";
    append_numeric(output, comparison.boundary_distance, "    ");
    output << ",\n    \"boundary_segment_id\": ";
    append_exact(output, comparison.boundary_segment_id, "    ");
    output << ",\n    \"nearest_boundary_id\": ";
    append_exact(output, comparison.nearest_boundary_id, "    ");
    output << ",\n    \"deforming_region_id\": ";
    append_exact(output, comparison.deforming_region_id, "    ");
    output << ",\n    \"deforming_region_type\": ";
    append_exact(output, comparison.deforming_region_type, "    ");
    output << ",\n    \"deformation_rate\": ";
    append_numeric(output, comparison.deformation_rate, "    ");
    output << ",\n    \"deformation_velocity_x\": ";
    append_numeric(output, comparison.deformation_velocity_x, "    ");
    output << ",\n    \"deformation_velocity_y\": ";
    append_numeric(output, comparison.deformation_velocity_y, "    ");
    output << ",\n    \"convergence_score\": ";
    append_exact(output, comparison.convergence_score, "    ");
    output << ",\n    \"divergence_score\": ";
    append_exact(output, comparison.divergence_score, "    ");
    output << ",\n    \"shear_score\": ";
    append_exact(output, comparison.shear_score, "    ");
    output << ",\n    \"geologic_regime\": ";
    append_exact(output, comparison.geologic_regime, "    ");
    output << "\n  }\n";
    output << "}\n";
    return output.str();
}

} // namespace tectonics::tooling
