#include "lithosphere.hpp"
#include "tectonic_contract.hpp"
#include "topography_codec.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
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

namespace {

struct Params {
    long seed = 12345;
    uint32_t width = 96;
    uint32_t height = 48;
    fs::path output_dir;
    float sea_level = 0.65f;
    uint32_t erosion_period = 0;
    float folding_ratio = 0.02f;
    uint32_t aggregation_overlap_abs = 1000000;
    float aggregation_overlap_rel = 0.33f;
    uint32_t cycle_count = 4;
    uint32_t plate_count = 4;
    float erosion_strength = 1.0f;
    float crust_rotation_strength = 0.20f;
    float rotation_strength = 1.0f;
    float subduction_strength = 1.0f;
    int32_t sea_level_m_override = TopographyCodec::kNoSeaLevelOverride;
    uint16_t initial_min_height_m = TopographyCodec::kDefaultInitialMinHeightMeters;
    uint16_t initial_max_height_m = TopographyCodec::kDefaultInitialMaxHeightMeters;
    uint32_t cycle_step_limit = 600;
    float divergent_carve_strength = 0.015f;
    double delta_time_myr = platec::scenario::kDefaultDeltaTimeMyr;
};

[[noreturn]] void fail(const std::string& message) {
    throw std::runtime_error(message);
}

uint32_t parse_u32(const char* flag, const char* value) {
    char* end = nullptr;
    const unsigned long parsed = std::strtoul(value, &end, 10);
    if (end == value || *end != '\0' ||
        parsed > static_cast<unsigned long>(std::numeric_limits<uint32_t>::max())) {
        fail(std::string("invalid value for ") + flag + ": " + value);
    }
    return static_cast<uint32_t>(parsed);
}

long parse_long(const char* flag, const char* value) {
    char* end = nullptr;
    const long parsed = std::strtol(value, &end, 10);
    if (end == value || *end != '\0') {
        fail(std::string("invalid value for ") + flag + ": " + value);
    }
    return parsed;
}

float parse_float(const char* flag, const char* value) {
    char* end = nullptr;
    const float parsed = std::strtof(value, &end);
    if (end == value || *end != '\0' || !std::isfinite(parsed)) {
        fail(std::string("invalid value for ") + flag + ": " + value);
    }
    return parsed;
}

double parse_double(const char* flag, const char* value) {
    char* end = nullptr;
    const double parsed = std::strtod(value, &end);
    if (end == value || *end != '\0' || !std::isfinite(parsed)) {
        fail(std::string("invalid value for ") + flag + ": " + value);
    }
    return parsed;
}

void print_help() {
    std::printf("contract_fixture --output DIR [options]\n");
    std::printf("  --seed N\n");
    std::printf("  --width N --height N\n");
    std::printf("  --sea-level X\n");
    std::printf("  --erosion-period N\n");
    std::printf("  --folding-ratio X\n");
    std::printf("  --aggregation-overlap-abs N\n");
    std::printf("  --aggregation-overlap-rel X\n");
    std::printf("  --cycles N\n");
    std::printf("  --plates N\n");
    std::printf("  --cycle-step-limit N\n");
    std::printf("  --delta-time-myr X\n");
}

Params parse_args(int argc, char** argv) {
    Params params;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            print_help();
            std::exit(0);
        }
        if (i + 1 >= argc) {
            fail("missing value for " + arg);
        }

        const char* value = argv[++i];
        if (arg == "--output") {
            params.output_dir = value;
        } else if (arg == "--seed") {
            params.seed = parse_long("--seed", value);
        } else if (arg == "--width") {
            params.width = parse_u32("--width", value);
        } else if (arg == "--height") {
            params.height = parse_u32("--height", value);
        } else if (arg == "--sea-level") {
            params.sea_level = parse_float("--sea-level", value);
        } else if (arg == "--erosion-period") {
            params.erosion_period = parse_u32("--erosion-period", value);
        } else if (arg == "--folding-ratio") {
            params.folding_ratio = parse_float("--folding-ratio", value);
        } else if (arg == "--aggregation-overlap-abs") {
            params.aggregation_overlap_abs = parse_u32("--aggregation-overlap-abs", value);
        } else if (arg == "--aggregation-overlap-rel") {
            params.aggregation_overlap_rel = parse_float("--aggregation-overlap-rel", value);
        } else if (arg == "--cycles") {
            params.cycle_count = parse_u32("--cycles", value);
        } else if (arg == "--plates") {
            params.plate_count = parse_u32("--plates", value);
        } else if (arg == "--cycle-step-limit") {
            params.cycle_step_limit = parse_u32("--cycle-step-limit", value);
        } else if (arg == "--delta-time-myr") {
            params.delta_time_myr = parse_double("--delta-time-myr", value);
        } else {
            fail("unknown argument: " + arg);
        }
    }

    if (params.output_dir.empty()) {
        fail("--output is required");
    }

    return params;
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

void write_text_file(const fs::path& path, const std::string& contents) {
    std::ofstream output(path);
    if (!output) {
        fail("failed to open " + path.string());
    }
    output << contents;
    if (!output) {
        fail("failed to write " + path.string());
    }
}

template <typename T>
struct RangeStats {
    T min;
    T max;
    double mean;
};

template <typename T>
RangeStats<T> compute_stats(const std::vector<T>& values) {
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

std::vector<uint16_t> to_metric_heightmap(const platec::contract::Snapshot& snapshot) {
    std::vector<uint16_t> values(snapshot.cell_count());
    for (size_t i = 0; i < snapshot.cell_count(); ++i) {
        values[i] = TopographyCodec::internal_to_meters(snapshot.heightmap[i], snapshot.sea_level_m);
    }
    return values;
}

std::array<uint32_t, 8> regime_counts(const platec::contract::Snapshot& snapshot) {
    std::array<uint32_t, 8> counts{};
    for (uint8_t regime : snapshot.geologic_regime) {
        if (regime < counts.size()) {
            ++counts[regime];
        }
    }
    return counts;
}

std::vector<uint32_t> plate_cell_counts(const platec::contract::Snapshot& snapshot) {
    std::vector<uint32_t> counts(snapshot.plates.size(), 0);
    for (uint32_t plate_id : snapshot.plate_id) {
        if (plate_id < counts.size()) {
            ++counts[plate_id];
        }
    }
    return counts;
}

std::string build_boundary_segments_json(const platec::contract::Snapshot& snapshot) {
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
               << platec::contract::boundary_type_name(segment.boundary_type) << "\",\n";
        output << "    \"geologic_regime\": \""
               << platec::contract::geologic_regime_name(segment.geologic_regime) << "\",\n";
        output << "    \"age_myr\": " << segment.age_myr << "\n";
        output << "  }" << (i + 1 < snapshot.boundary_segments.size() ? "," : "") << "\n";
    }
    output << "]\n";
    return output.str();
}

std::string build_deforming_regions_json(const platec::contract::Snapshot& snapshot) {
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
        output << "    \"average_deformation_rate\": " << region.average_deformation_rate << ",\n";
        output << "    \"average_interpolated_velocity\": ["
               << region.average_interpolated_velocity_x << ", "
               << region.average_interpolated_velocity_y << "],\n";
        output << "    \"average_normal_motion\": " << region.average_normal_motion << ",\n";
        output << "    \"average_shear_motion\": " << region.average_shear_motion << ",\n";
        output << "    \"type\": \""
               << platec::contract::deforming_region_type_name(region.type) << "\",\n";
        output << "    \"age_myr\": " << region.age_myr << "\n";
        output << "  }" << (i + 1 < snapshot.deforming_regions.size() ? "," : "") << "\n";
    }
    output << "]\n";
    return output.str();
}

std::string build_tectonic_junctions_json(const platec::contract::Snapshot& snapshot) {
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
               << platec::contract::junction_arm_type_name(junction.arm_types[0]) << "\", \""
               << platec::contract::junction_arm_type_name(junction.arm_types[1]) << "\", \""
               << platec::contract::junction_arm_type_name(junction.arm_types[2])
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
               << platec::contract::junction_arm_type_code(junction.arm_types[0]) << "-"
               << platec::contract::junction_arm_type_code(junction.arm_types[1]) << "-"
               << platec::contract::junction_arm_type_code(junction.arm_types[2])
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
               << platec::contract::junction_stability_name(junction.stability) << "\",\n";
        output << "    \"stability_score\": " << junction.stability_score << ",\n";
        output << "    \"quality_score\": " << junction.quality_score << "\n";
        output << "  }" << (i + 1 < snapshot.tectonic_junctions.size() ? "," : "") << "\n";
    }
    output << "]\n";
    return output.str();
}

void write_manifest(const fs::path& path, const Params& params,
                    const platec::contract::Snapshot& snapshot) {
    const std::vector<uint16_t> metric_heightmap = to_metric_heightmap(snapshot);
    const RangeStats<float> height_stats = compute_stats(snapshot.heightmap);
    const RangeStats<uint32_t> age_stats = compute_stats(snapshot.crust_age_steps);
    const RangeStats<uint16_t> metric_stats = compute_stats(metric_heightmap);
    const auto regimes = regime_counts(snapshot);
    const std::vector<uint32_t> plate_counts = plate_cell_counts(snapshot);

    std::ofstream output(path);
    if (!output) {
        fail("failed to open " + path.string());
    }

    output << "{\n";
    output << std::setprecision(17);
    output << "  \"schema\": \"tectonic-contract-fixture/v" << snapshot.schema_version << "\",\n";
    output << "  \"schema_version\": " << snapshot.schema_version << ",\n";
    output << "  \"seed\": " << params.seed << ",\n";
    output << "  \"width\": " << snapshot.width << ",\n";
    output << "  \"height\": " << snapshot.height << ",\n";
    output << "  \"iteration_count\": " << snapshot.iteration_count << ",\n";
    output << "  \"cycle_count\": " << snapshot.cycle_count << ",\n";
    output << "  \"time_origin_step\": " << snapshot.time_origin_step << ",\n";
    output << "  \"time_myr\": " << snapshot.time_myr << ",\n";
    output << "  \"delta_time_myr\": " << snapshot.delta_time_myr << ",\n";
    output << "  \"sea_level_m\": " << snapshot.sea_level_m << ",\n";
    output << "  \"parameters\": {\n";
    output << "    \"sea_level\": " << params.sea_level << ",\n";
    output << "    \"erosion_period\": " << params.erosion_period << ",\n";
    output << "    \"folding_ratio\": " << params.folding_ratio << ",\n";
    output << "    \"aggregation_overlap_abs\": " << params.aggregation_overlap_abs << ",\n";
    output << "    \"aggregation_overlap_rel\": " << params.aggregation_overlap_rel << ",\n";
    output << "    \"cycle_count\": " << params.cycle_count << ",\n";
    output << "    \"plate_count\": " << params.plate_count << ",\n";
    output << "    \"cycle_step_limit\": " << params.cycle_step_limit << ",\n";
    output << "    \"delta_time_myr\": " << params.delta_time_myr << "\n";
    output << "  },\n";
    output << "  \"scenario\": {\n";
    output << "    \"version\": " << snapshot.run_scenario.version << ",\n";
    output << "    \"seed\": " << snapshot.run_scenario.seed << ",\n";
    output << "    \"width\": " << snapshot.run_scenario.width << ",\n";
    output << "    \"height\": " << snapshot.run_scenario.height << ",\n";
    output << "    \"sea_level\": " << snapshot.run_scenario.sea_level << ",\n";
    output << "    \"cycle_count\": " << snapshot.run_scenario.cycle_count << ",\n";
    output << "    \"plate_count\": " << snapshot.run_scenario.plate_count << ",\n";
    output << "    \"delta_time_myr\": " << snapshot.run_scenario.delta_time_myr << "\n";
    output << "  },\n";
    output << "  \"files\": {\n";
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
    output << "    \"heightmap\": \"" << fnv1a_hash(snapshot.heightmap) << "\",\n";
    output << "    \"plate_id\": \"" << fnv1a_hash(snapshot.plate_id) << "\",\n";
    output << "    \"crust_age_steps\": \"" << fnv1a_hash(snapshot.crust_age_steps) << "\",\n";
    output << "    \"crust_age_myr\": \"" << fnv1a_hash(snapshot.crust_age_myr) << "\",\n";
    output << "    \"crust_thickness\": \"" << fnv1a_hash(snapshot.crust_thickness) << "\",\n";
    output << "    \"crust_class\": \"" << fnv1a_hash(snapshot.crust_class) << "\",\n";
    output << "    \"uplift_tendency\": \"" << fnv1a_hash(snapshot.uplift_tendency) << "\",\n";
    output << "    \"subsidence_tendency\": \"" << fnv1a_hash(snapshot.subsidence_tendency) << "\",\n";
    output << "    \"accumulated_strain\": \"" << fnv1a_hash(snapshot.accumulated_strain) << "\",\n";
    output << "    \"boundary_type\": \"" << fnv1a_hash(snapshot.boundary_type) << "\",\n";
    output << "    \"boundary_distance\": \"" << fnv1a_hash(snapshot.boundary_distance) << "\",\n";
    output << "    \"boundary_segment_id\": \"" << fnv1a_hash(snapshot.boundary_segment_id) << "\",\n";
    output << "    \"nearest_boundary_id\": \"" << fnv1a_hash(snapshot.nearest_boundary_id) << "\",\n";
    output << "    \"deforming_region_id\": \"" << fnv1a_hash(snapshot.deforming_region_id) << "\",\n";
    output << "    \"deforming_region_type\": \"" << fnv1a_hash(snapshot.deforming_region_type) << "\",\n";
    output << "    \"deformation_rate\": \"" << fnv1a_hash(snapshot.deformation_rate) << "\",\n";
    output << "    \"deformation_velocity_x\": \"" << fnv1a_hash(snapshot.deformation_velocity_x) << "\",\n";
    output << "    \"deformation_velocity_y\": \"" << fnv1a_hash(snapshot.deformation_velocity_y) << "\",\n";
    output << "    \"convergence_score\": \"" << fnv1a_hash(snapshot.convergence_score) << "\",\n";
    output << "    \"divergence_score\": \"" << fnv1a_hash(snapshot.divergence_score) << "\",\n";
    output << "    \"shear_score\": \"" << fnv1a_hash(snapshot.shear_score) << "\",\n";
    output << "    \"geologic_regime\": \"" << fnv1a_hash(snapshot.geologic_regime) << "\",\n";
    output << "    \"boundary_segments\": \""
           << fnv1a_hash(build_boundary_segments_json(snapshot)) << "\",\n";
    output << "    \"deforming_regions\": \""
           << fnv1a_hash(build_deforming_regions_json(snapshot)) << "\",\n";
    output << "    \"tectonic_junctions\": \""
           << fnv1a_hash(build_tectonic_junctions_json(snapshot)) << "\"\n";
    output << "  },\n";
    output << "  \"heightmap_stats\": {\n";
    output << "    \"internal_min\": " << height_stats.min << ",\n";
    output << "    \"internal_max\": " << height_stats.max << ",\n";
    output << "    \"internal_mean\": " << height_stats.mean << ",\n";
    output << "    \"metric_min\": " << metric_stats.min << ",\n";
    output << "    \"metric_max\": " << metric_stats.max << ",\n";
    output << "    \"metric_mean\": " << metric_stats.mean << "\n";
    output << "  },\n";
    output << "  \"crust_age_stats\": {\n";
    output << "    \"min\": " << age_stats.min << ",\n";
    output << "    \"max\": " << age_stats.max << ",\n";
    output << "    \"mean\": " << age_stats.mean << "\n";
    output << "  },\n";
    output << "  \"regime_counts\": {\n";
    for (size_t i = 0; i < regimes.size(); ++i) {
        output << "    \""
               << platec::contract::geologic_regime_name(
                      static_cast<platec::contract::GeologicRegime>(i))
               << "\": " << regimes[i];
        output << (i + 1 < regimes.size() ? ",\n" : "\n");
    }
    output << "  },\n";
    output << "  \"boundary_segment_count\": " << snapshot.boundary_segments.size() << ",\n";
    output << "  \"deforming_region_count\": " << snapshot.deforming_regions.size() << ",\n";
    output << "  \"tectonic_junction_count\": " << snapshot.tectonic_junctions.size() << ",\n";
    output << "  \"plate_cell_counts\": [";
    for (size_t i = 0; i < plate_counts.size(); ++i) {
        output << plate_counts[i];
        output << (i + 1 < plate_counts.size() ? ", " : "");
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
}

} // namespace

int main(int argc, char** argv) {
    try {
        const Params params = parse_args(argc, argv);
        fs::create_directories(params.output_dir);

        platec::scenario::Scenario scenario;
        scenario.seed = params.seed;
        scenario.width = params.width;
        scenario.height = params.height;
        scenario.sea_level = params.sea_level;
        scenario.erosion_period = params.erosion_period;
        scenario.folding_ratio = params.folding_ratio;
        scenario.aggregation_overlap_abs = params.aggregation_overlap_abs;
        scenario.aggregation_overlap_rel = params.aggregation_overlap_rel;
        scenario.cycle_count = params.cycle_count;
        scenario.plate_count = params.plate_count;
        scenario.erosion_strength = params.erosion_strength;
        scenario.crust_rotation_strength = params.crust_rotation_strength;
        scenario.rotation_strength = params.rotation_strength;
        scenario.subduction_strength = params.subduction_strength;
        scenario.sea_level_m = params.sea_level_m_override;
        scenario.initial_min_height_m = params.initial_min_height_m;
        scenario.initial_max_height_m = params.initial_max_height_m;
        scenario.cycle_step_limit = params.cycle_step_limit;
        scenario.divergent_carve_strength = params.divergent_carve_strength;
        scenario.delta_time_myr = params.delta_time_myr;

        lithosphere simulation(scenario);

        while (!simulation.isFinished()) {
            simulation.update();
        }
        const platec::contract::Snapshot snapshot =
            platec::contract::capture_snapshot(simulation);

        write_binary_file(params.output_dir / "heightmap.f32", snapshot.heightmap);
        write_binary_file(params.output_dir / "plate_id.u32", snapshot.plate_id);
        write_binary_file(params.output_dir / "crust_age_steps.u32", snapshot.crust_age_steps);
        write_binary_file(params.output_dir / "crust_age_myr.f32", snapshot.crust_age_myr);
        write_binary_file(params.output_dir / "crust_thickness.f32", snapshot.crust_thickness);
        write_binary_file(params.output_dir / "crust_class.u8", snapshot.crust_class);
        write_binary_file(params.output_dir / "uplift_tendency.f32", snapshot.uplift_tendency);
        write_binary_file(params.output_dir / "subsidence_tendency.f32",
                          snapshot.subsidence_tendency);
        write_binary_file(params.output_dir / "accumulated_strain.f32",
                          snapshot.accumulated_strain);
        write_binary_file(params.output_dir / "boundary_type.u8", snapshot.boundary_type);
        write_binary_file(params.output_dir / "boundary_distance.u16",
                          snapshot.boundary_distance);
        write_binary_file(params.output_dir / "boundary_segment_id.u32",
                          snapshot.boundary_segment_id);
        write_binary_file(params.output_dir / "nearest_boundary_id.u32",
                          snapshot.nearest_boundary_id);
        write_binary_file(params.output_dir / "deforming_region_id.u32",
                          snapshot.deforming_region_id);
        write_binary_file(params.output_dir / "deforming_region_type.u8",
                          snapshot.deforming_region_type);
        write_binary_file(params.output_dir / "deformation_rate.f32",
                          snapshot.deformation_rate);
        write_binary_file(params.output_dir / "deformation_velocity_x.f32",
                          snapshot.deformation_velocity_x);
        write_binary_file(params.output_dir / "deformation_velocity_y.f32",
                          snapshot.deformation_velocity_y);
        write_binary_file(params.output_dir / "convergence_score.u8",
                          snapshot.convergence_score);
        write_binary_file(params.output_dir / "divergence_score.u8",
                          snapshot.divergence_score);
        write_binary_file(params.output_dir / "shear_score.u8", snapshot.shear_score);
        write_binary_file(params.output_dir / "geologic_regime.u8",
                          snapshot.geologic_regime);
        write_text_file(params.output_dir / "boundary_segments.json",
                        build_boundary_segments_json(snapshot));
        write_text_file(params.output_dir / "deforming_regions.json",
                        build_deforming_regions_json(snapshot));
        write_text_file(params.output_dir / "tectonic_junctions.json",
                        build_tectonic_junctions_json(snapshot));
        write_manifest(params.output_dir / "manifest.json", params, snapshot);

        std::printf("fixture written to %s\n", params.output_dir.string().c_str());
        return 0;
    } catch (const std::exception& exception) {
        std::fprintf(stderr, "error: %s\n", exception.what());
        return 1;
    }
}
