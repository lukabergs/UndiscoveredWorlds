#ifndef NOMINMAX
#define NOMINMAX
#endif

#include "lithosphere.hpp"
#include "tectonic_pipeline_support.hpp"
#include "tectonic_scenario.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

enum class TargetMode {
    Final,
    UpdateCount,
};

struct TargetSelector {
    TargetMode mode = TargetMode::Final;
    uint32_t update_count = 0;
    double requested_time_myr = 0.0;
    bool came_from_time = false;
};

struct BatchExportParams {
    platec::scenario::Scenario scenario;
    fs::path output_dir;
    bool include_initial = true;
    bool include_final = true;
    uint32_t checkpoint_every_updates = 0;
    std::vector<uint32_t> explicit_updates;
    std::vector<double> explicit_times_myr;
};

struct InspectParams {
    platec::scenario::Scenario scenario;
    TargetSelector target;
    bool has_output = false;
    fs::path output_path;
    bool has_bundle_output = false;
    fs::path bundle_output_dir;
};

struct CompareParams {
    bool bundle_mode = false;
    fs::path bundle_a;
    fs::path bundle_b;
    platec::scenario::Scenario shared_scenario;
    long seed_a = 0;
    long seed_b = 0;
    bool has_seed_a = false;
    bool has_seed_b = false;
    TargetSelector target;
    bool has_output = false;
    fs::path output_path;
};

struct ExportedSnapshotRecord {
    std::string label;
    uint32_t update_count = 0;
    uint32_t iteration_count = 0;
    uint32_t cycle_count = 0;
    double time_myr = 0.0;
    std::string directory;
};

struct RequestedCheckpoint {
    uint32_t update_count = 0;
    std::string label;
    bool from_time = false;
    double requested_time_myr = 0.0;
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

int32_t parse_i32(const char* flag, const char* value) {
    char* end = nullptr;
    const long parsed = std::strtol(value, &end, 10);
    if (end == value || *end != '\0' ||
        parsed < static_cast<long>(std::numeric_limits<int32_t>::min()) ||
        parsed > static_cast<long>(std::numeric_limits<int32_t>::max())) {
        fail(std::string("invalid value for ") + flag + ": " + value);
    }
    return static_cast<int32_t>(parsed);
}

uint16_t parse_u16(const char* flag, const char* value) {
    const uint32_t parsed = parse_u32(flag, value);
    if (parsed > static_cast<uint32_t>(std::numeric_limits<uint16_t>::max())) {
        fail(std::string("invalid value for ") + flag + ": " + value);
    }
    return static_cast<uint16_t>(parsed);
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

platec::scenario::Scenario default_scenario() {
    platec::scenario::Scenario scenario;
    scenario.seed = 12345;
    return scenario;
}

void print_general_help() {
    std::printf("tectonic_pipeline <command> [options]\n");
    std::printf("Commands:\n");
    std::printf("  batch-export   Export structured snapshot bundles for a run\n");
    std::printf("  inspect        Inspect a seed/scenario at a specific checkpoint\n");
    std::printf("  boundary-stats Dump boundary and deforming-region stats at a checkpoint\n");
    std::printf("  compare        Compare two seeds or two exported snapshot bundles\n");
}

void print_common_scenario_help() {
    std::printf("Common scenario options:\n");
    std::printf("  --seed N\n");
    std::printf("  --width N --height N\n");
    std::printf("  --sea-level X\n");
    std::printf("  --erosion-period N\n");
    std::printf("  --folding-ratio X\n");
    std::printf("  --aggregation-overlap-abs N\n");
    std::printf("  --aggregation-overlap-rel X\n");
    std::printf("  --cycles N\n");
    std::printf("  --plates N\n");
    std::printf("  --erosion-strength X\n");
    std::printf("  --crust-rotation-strength X\n");
    std::printf("  --rotation-strength X\n");
    std::printf("  --subduction-strength X\n");
    std::printf("  --sea-level-m N\n");
    std::printf("  --initial-min-height-m N\n");
    std::printf("  --initial-max-height-m N\n");
    std::printf("  --cycle-step-limit N\n");
    std::printf("  --divergent-carve X\n");
    std::printf("  --delta-time-myr X\n");
}

void print_batch_help() {
    std::printf("tectonic_pipeline batch-export --output DIR [options]\n");
    std::printf("  --checkpoint-every-updates N\n");
    std::printf("  --checkpoint-at-update N\n");
    std::printf("  --checkpoint-at-time-myr X\n");
    std::printf("  --no-initial\n");
    std::printf("  --no-final\n");
    print_common_scenario_help();
}

void print_inspect_help(const char* command_name) {
    std::printf("tectonic_pipeline %s [target options] [scenario options]\n", command_name);
    std::printf("Target options:\n");
    std::printf("  --final\n");
    std::printf("  --update N\n");
    std::printf("  --time-myr X\n");
    if (std::string(command_name) == "inspect") {
        std::printf("  --output FILE\n");
        std::printf("  --bundle-output DIR\n");
    } else {
        std::printf("  --output FILE\n");
    }
    print_common_scenario_help();
}

void print_compare_help() {
    std::printf("tectonic_pipeline compare [compare mode] [target options] [scenario options]\n");
    std::printf("Compare modes:\n");
    std::printf("  --seed-a N --seed-b N\n");
    std::printf("  --bundle-a DIR --bundle-b DIR\n");
    std::printf("Target options for live seed compare:\n");
    std::printf("  --final\n");
    std::printf("  --update N\n");
    std::printf("  --time-myr X\n");
    std::printf("Other options:\n");
    std::printf("  --output FILE\n");
    print_common_scenario_help();
}

uint32_t updates_for_time_myr(double time_myr, double delta_time_myr) {
    if (time_myr < 0.0) {
        fail("--time-myr must be >= 0");
    }
    const double raw = time_myr / delta_time_myr;
    const double rounded = std::round(raw);
    if (std::fabs(raw - rounded) > 1e-9) {
        std::ostringstream message;
        message << "--time-myr " << time_myr
                << " is not aligned to delta_time_myr=" << delta_time_myr;
        fail(message.str());
    }
    return static_cast<uint32_t>(rounded);
}

std::string checkpoint_label_from_update(uint32_t update_count) {
    std::ostringstream output;
    output << "update" << std::setfill('0') << std::setw(6) << update_count;
    return output.str();
}

std::string checkpoint_label_from_time(double time_myr) {
    std::ostringstream output;
    output << std::fixed << std::setprecision(3) << "time" << time_myr << "Myr";
    return output.str();
}

std::string sanitize_label(const std::string& value) {
    std::string sanitized;
    sanitized.reserve(value.size());
    for (char ch : value) {
        if ((ch >= 'a' && ch <= 'z') || (ch >= 'A' && ch <= 'Z') ||
            (ch >= '0' && ch <= '9')) {
            sanitized.push_back(ch);
        } else {
            sanitized.push_back('_');
        }
    }
    return sanitized;
}

fs::path snapshot_directory_name(size_t index, const std::string& label) {
    std::ostringstream output;
    output << std::setfill('0') << std::setw(3) << index << "_" << sanitize_label(label);
    return fs::path(output.str());
}

void write_output_or_stdout(const std::string& contents, bool has_output, const fs::path& output_path) {
    if (!has_output) {
        std::printf("%s", contents.c_str());
        return;
    }

    if (!output_path.parent_path().empty()) {
        fs::create_directories(output_path.parent_path());
    }

    std::ofstream output(output_path, std::ios::binary);
    if (!output) {
        fail("failed to open " + output_path.string());
    }
    output << contents;
    if (!output) {
        fail("failed to write " + output_path.string());
    }
}

void parse_common_scenario_arg(platec::scenario::Scenario& scenario, const std::string& arg,
                               const char* value) {
    if (arg == "--seed") {
        scenario.seed = parse_long("--seed", value);
    } else if (arg == "--width") {
        scenario.width = parse_u32("--width", value);
    } else if (arg == "--height") {
        scenario.height = parse_u32("--height", value);
    } else if (arg == "--sea-level") {
        scenario.sea_level = parse_float("--sea-level", value);
    } else if (arg == "--erosion-period") {
        scenario.erosion_period = parse_u32("--erosion-period", value);
    } else if (arg == "--folding-ratio") {
        scenario.folding_ratio = parse_float("--folding-ratio", value);
    } else if (arg == "--aggregation-overlap-abs") {
        scenario.aggregation_overlap_abs = parse_u32("--aggregation-overlap-abs", value);
    } else if (arg == "--aggregation-overlap-rel") {
        scenario.aggregation_overlap_rel = parse_float("--aggregation-overlap-rel", value);
    } else if (arg == "--cycles") {
        scenario.cycle_count = parse_u32("--cycles", value);
    } else if (arg == "--plates") {
        scenario.plate_count = parse_u32("--plates", value);
    } else if (arg == "--erosion-strength") {
        scenario.erosion_strength = parse_float("--erosion-strength", value);
    } else if (arg == "--crust-rotation-strength") {
        scenario.crust_rotation_strength = parse_float("--crust-rotation-strength", value);
    } else if (arg == "--rotation-strength") {
        scenario.rotation_strength = parse_float("--rotation-strength", value);
    } else if (arg == "--subduction-strength") {
        scenario.subduction_strength = parse_float("--subduction-strength", value);
    } else if (arg == "--sea-level-m") {
        scenario.sea_level_m = parse_i32("--sea-level-m", value);
    } else if (arg == "--initial-min-height-m") {
        scenario.initial_min_height_m = parse_u16("--initial-min-height-m", value);
    } else if (arg == "--initial-max-height-m") {
        scenario.initial_max_height_m = parse_u16("--initial-max-height-m", value);
    } else if (arg == "--cycle-step-limit") {
        scenario.cycle_step_limit = parse_u32("--cycle-step-limit", value);
    } else if (arg == "--divergent-carve") {
        scenario.divergent_carve_strength = parse_float("--divergent-carve", value);
    } else if (arg == "--delta-time-myr") {
        scenario.delta_time_myr = parse_double("--delta-time-myr", value);
    } else {
        fail("unknown argument: " + arg);
    }
}

TargetSelector parse_target_selector(const platec::scenario::Scenario& scenario,
                                     const std::vector<std::string>& args, size_t& index) {
    TargetSelector target;
    bool seen_target_flag = false;

    while (index < args.size()) {
        const std::string& arg = args[index];
        if (arg == "--final") {
            if (seen_target_flag) {
                fail("multiple target flags provided");
            }
            target.mode = TargetMode::Final;
            seen_target_flag = true;
            ++index;
        } else if (arg == "--update") {
            if (seen_target_flag || index + 1 >= args.size()) {
                fail("invalid --update usage");
            }
            target.mode = TargetMode::UpdateCount;
            target.update_count = parse_u32("--update", args[index + 1].c_str());
            seen_target_flag = true;
            index += 2;
        } else if (arg == "--time-myr") {
            if (seen_target_flag || index + 1 >= args.size()) {
                fail("invalid --time-myr usage");
            }
            target.requested_time_myr = parse_double("--time-myr", args[index + 1].c_str());
            target.update_count =
                updates_for_time_myr(target.requested_time_myr, scenario.delta_time_myr);
            target.mode = TargetMode::UpdateCount;
            target.came_from_time = true;
            seen_target_flag = true;
            index += 2;
        } else {
            break;
        }
    }

    return target;
}

std::vector<std::string> collect_args(int argc, char** argv, int start) {
    std::vector<std::string> args;
    args.reserve(static_cast<size_t>(argc - start));
    for (int i = start; i < argc; ++i) {
        args.emplace_back(argv[i]);
    }
    return args;
}

platec::contract::Snapshot run_to_target(const platec::scenario::Scenario& scenario,
                                         const TargetSelector& target,
                                         uint32_t& resolved_update_count) {
    lithosphere simulation(scenario);
    resolved_update_count = 0;

    if (target.mode == TargetMode::Final) {
        while (!simulation.isFinished()) {
            simulation.update();
            ++resolved_update_count;
        }
        return platec::contract::capture_snapshot(simulation);
    }

    while (resolved_update_count < target.update_count && !simulation.isFinished()) {
        simulation.update();
        ++resolved_update_count;
    }
    if (resolved_update_count < target.update_count) {
        std::ostringstream message;
        message << "simulation finished at update " << resolved_update_count << " (time "
                << simulation.getTimeMyr() << " Myr) before requested target "
                << target.update_count;
        fail(message.str());
    }
    return platec::contract::capture_snapshot(simulation);
}

void write_batch_manifest(const fs::path& path, const BatchExportParams& params,
                          const std::vector<ExportedSnapshotRecord>& records,
                          const std::vector<RequestedCheckpoint>& missing_targets) {
    std::ofstream output(path);
    if (!output) {
        fail("failed to open " + path.string());
    }

    output << std::setprecision(17);
    output << "{\n";
    output << "  \"schema\": \"tectonic-batch-export/v1\",\n";
    output << "  \"output_dir\": \"" << json_escape(path.parent_path().string()) << "\",\n";
    output << "  \"scenario\": {\n";
    output << "    \"seed\": " << params.scenario.seed << ",\n";
    output << "    \"width\": " << params.scenario.width << ",\n";
    output << "    \"height\": " << params.scenario.height << ",\n";
    output << "    \"sea_level\": " << params.scenario.sea_level << ",\n";
    output << "    \"erosion_period\": " << params.scenario.erosion_period << ",\n";
    output << "    \"folding_ratio\": " << params.scenario.folding_ratio << ",\n";
    output << "    \"aggregation_overlap_abs\": " << params.scenario.aggregation_overlap_abs << ",\n";
    output << "    \"aggregation_overlap_rel\": " << params.scenario.aggregation_overlap_rel << ",\n";
    output << "    \"cycle_count\": " << params.scenario.cycle_count << ",\n";
    output << "    \"plate_count\": " << params.scenario.plate_count << ",\n";
    output << "    \"erosion_strength\": " << params.scenario.erosion_strength << ",\n";
    output << "    \"crust_rotation_strength\": " << params.scenario.crust_rotation_strength << ",\n";
    output << "    \"rotation_strength\": " << params.scenario.rotation_strength << ",\n";
    output << "    \"subduction_strength\": " << params.scenario.subduction_strength << ",\n";
    output << "    \"sea_level_m\": " << params.scenario.sea_level_m << ",\n";
    output << "    \"initial_min_height_m\": " << params.scenario.initial_min_height_m << ",\n";
    output << "    \"initial_max_height_m\": " << params.scenario.initial_max_height_m << ",\n";
    output << "    \"cycle_step_limit\": " << params.scenario.cycle_step_limit << ",\n";
    output << "    \"divergent_carve_strength\": " << params.scenario.divergent_carve_strength << ",\n";
    output << "    \"delta_time_myr\": " << params.scenario.delta_time_myr << "\n";
    output << "  },\n";
    output << "  \"targets\": {\n";
    output << "    \"include_initial\": " << (params.include_initial ? "true" : "false") << ",\n";
    output << "    \"include_final\": " << (params.include_final ? "true" : "false") << ",\n";
    output << "    \"checkpoint_every_updates\": " << params.checkpoint_every_updates << ",\n";
    output << "    \"explicit_updates\": [";
    for (size_t i = 0; i < params.explicit_updates.size(); ++i) {
        output << params.explicit_updates[i]
               << (i + 1 < params.explicit_updates.size() ? ", " : "");
    }
    output << "],\n";
    output << "    \"explicit_times_myr\": [";
    for (size_t i = 0; i < params.explicit_times_myr.size(); ++i) {
        output << params.explicit_times_myr[i]
               << (i + 1 < params.explicit_times_myr.size() ? ", " : "");
    }
    output << "]\n";
    output << "  },\n";
    output << "  \"snapshots\": [\n";
    for (size_t i = 0; i < records.size(); ++i) {
        const ExportedSnapshotRecord& record = records[i];
        output << "    {\n";
        output << "      \"label\": \"" << json_escape(record.label) << "\",\n";
        output << "      \"update_count\": " << record.update_count << ",\n";
        output << "      \"iteration_count\": " << record.iteration_count << ",\n";
        output << "      \"cycle_count\": " << record.cycle_count << ",\n";
        output << "      \"time_myr\": " << record.time_myr << ",\n";
        output << "      \"directory\": \"" << json_escape(record.directory) << "\"\n";
        output << "    }" << (i + 1 < records.size() ? "," : "") << "\n";
    }
    output << "  ],\n";
    output << "  \"missing_targets\": [\n";
    for (size_t i = 0; i < missing_targets.size(); ++i) {
        const RequestedCheckpoint& missing = missing_targets[i];
        output << "    {\n";
        output << "      \"label\": \"" << json_escape(missing.label) << "\",\n";
        output << "      \"update_count\": " << missing.update_count;
        if (missing.from_time) {
            output << ",\n      \"requested_time_myr\": " << missing.requested_time_myr << "\n";
        } else {
            output << "\n";
        }
        output << "    }" << (i + 1 < missing_targets.size() ? "," : "") << "\n";
    }
    output << "  ]\n";
    output << "}\n";
    if (!output) {
        fail("failed to write " + path.string());
    }
}

BatchExportParams parse_batch_export_args(const std::vector<std::string>& args) {
    BatchExportParams params;
    params.scenario = default_scenario();

    for (size_t i = 0; i < args.size(); ++i) {
        const std::string& arg = args[i];
        if (arg == "--help" || arg == "-h") {
            print_batch_help();
            std::exit(0);
        }
        if (arg == "--no-initial") {
            params.include_initial = false;
            continue;
        }
        if (arg == "--no-final") {
            params.include_final = false;
            continue;
        }
        if (i + 1 >= args.size()) {
            fail("missing value for " + arg);
        }

        const char* value = args[i + 1].c_str();
        if (arg == "--output") {
            params.output_dir = value;
        } else if (arg == "--checkpoint-every-updates") {
            params.checkpoint_every_updates =
                parse_u32("--checkpoint-every-updates", value);
        } else if (arg == "--checkpoint-at-update") {
            params.explicit_updates.push_back(parse_u32("--checkpoint-at-update", value));
        } else if (arg == "--checkpoint-at-time-myr") {
            params.explicit_times_myr.push_back(
                parse_double("--checkpoint-at-time-myr", value));
        } else {
            parse_common_scenario_arg(params.scenario, arg, value);
        }
        ++i;
    }

    if (params.output_dir.empty()) {
        fail("--output is required");
    }
    if (params.scenario.delta_time_myr <= 0.0) {
        fail("--delta-time-myr must be > 0");
    }
    return params;
}

InspectParams parse_inspect_args(const std::vector<std::string>& args, bool boundary_stats_mode) {
    InspectParams params;
    params.scenario = default_scenario();

    for (size_t index = 0; index < args.size();) {
        const std::string& arg = args[index];
        if (arg == "--help" || arg == "-h") {
            print_inspect_help(boundary_stats_mode ? "boundary-stats" : "inspect");
            std::exit(0);
        }
        if (arg == "--final" || arg == "--update" || arg == "--time-myr") {
            params.target = parse_target_selector(params.scenario, args, index);
            continue;
        }
        if (index + 1 >= args.size()) {
            fail("missing value for " + arg);
        }

        const char* value = args[index + 1].c_str();
        if (arg == "--output") {
            params.has_output = true;
            params.output_path = value;
        } else if (!boundary_stats_mode && arg == "--bundle-output") {
            params.has_bundle_output = true;
            params.bundle_output_dir = value;
        } else {
            parse_common_scenario_arg(params.scenario, arg, value);
        }
        index += 2;
    }

    if (params.scenario.delta_time_myr <= 0.0) {
        fail("--delta-time-myr must be > 0");
    }
    if (params.target.came_from_time) {
        params.target.update_count =
            updates_for_time_myr(params.target.requested_time_myr, params.scenario.delta_time_myr);
    }
    return params;
}

CompareParams parse_compare_args(const std::vector<std::string>& args) {
    CompareParams params;
    params.shared_scenario = default_scenario();

    for (size_t index = 0; index < args.size();) {
        const std::string& arg = args[index];
        if (arg == "--help" || arg == "-h") {
            print_compare_help();
            std::exit(0);
        }
        if (arg == "--final" || arg == "--update" || arg == "--time-myr") {
            params.target = parse_target_selector(params.shared_scenario, args, index);
            continue;
        }
        if (index + 1 >= args.size()) {
            fail("missing value for " + arg);
        }

        const char* value = args[index + 1].c_str();
        if (arg == "--bundle-a") {
            params.bundle_mode = true;
            params.bundle_a = value;
        } else if (arg == "--bundle-b") {
            params.bundle_mode = true;
            params.bundle_b = value;
        } else if (arg == "--seed-a") {
            params.has_seed_a = true;
            params.seed_a = parse_long("--seed-a", value);
        } else if (arg == "--seed-b") {
            params.has_seed_b = true;
            params.seed_b = parse_long("--seed-b", value);
        } else if (arg == "--output") {
            params.has_output = true;
            params.output_path = value;
        } else {
            parse_common_scenario_arg(params.shared_scenario, arg, value);
        }
        index += 2;
    }

    if (params.bundle_mode) {
        if (params.bundle_a.empty() || params.bundle_b.empty()) {
            fail("--bundle-a and --bundle-b are both required in bundle compare mode");
        }
        return params;
    }

    if (!params.has_seed_a || !params.has_seed_b) {
        fail("compare requires either --bundle-a/--bundle-b or --seed-a/--seed-b");
    }
    if (params.shared_scenario.delta_time_myr <= 0.0) {
        fail("--delta-time-myr must be > 0");
    }
    if (params.target.came_from_time) {
        params.target.update_count = updates_for_time_myr(params.target.requested_time_myr,
                                                          params.shared_scenario.delta_time_myr);
    }
    return params;
}

std::string target_label(const TargetSelector& target, uint32_t resolved_update_count,
                         double resolved_time_myr) {
    if (target.mode == TargetMode::Final) {
        return "final";
    }
    if (target.came_from_time) {
        return checkpoint_label_from_time(target.requested_time_myr);
    }
    if (target.update_count == 0 && resolved_time_myr == 0.0) {
        return "initial";
    }
    return checkpoint_label_from_update(resolved_update_count);
}

int command_batch_export(const std::vector<std::string>& args) {
    const BatchExportParams params = parse_batch_export_args(args);
    fs::create_directories(params.output_dir);
    const fs::path snapshots_dir = params.output_dir / "snapshots";
    std::error_code ec;
    fs::remove_all(snapshots_dir, ec);
    ec.clear();
    fs::remove(params.output_dir / "batch_manifest.json", ec);
    fs::create_directories(snapshots_dir);

    std::vector<RequestedCheckpoint> explicit_targets;
    explicit_targets.reserve(params.explicit_updates.size() + params.explicit_times_myr.size());
    for (uint32_t update_count : params.explicit_updates) {
        explicit_targets.push_back({update_count, checkpoint_label_from_update(update_count), false, 0.0});
    }
    for (double time_myr : params.explicit_times_myr) {
        explicit_targets.push_back({updates_for_time_myr(time_myr, params.scenario.delta_time_myr),
                                    checkpoint_label_from_time(time_myr), true, time_myr});
    }
    std::stable_sort(explicit_targets.begin(), explicit_targets.end(),
                     [](const RequestedCheckpoint& lhs, const RequestedCheckpoint& rhs) {
                         return lhs.update_count < rhs.update_count;
                     });

    lithosphere simulation(params.scenario);
    uint32_t update_count = 0;
    size_t export_index = 0;
    size_t explicit_index = 0;
    std::vector<ExportedSnapshotRecord> records;

    const auto export_snapshot = [&](const std::string& label, uint32_t current_update) {
        const platec::contract::Snapshot snapshot = platec::contract::capture_snapshot(simulation);
        const fs::path relative_dir = fs::path("snapshots") / snapshot_directory_name(export_index, label);
        platec::tooling::SnapshotBundleOptions options;
        options.label = label;
        options.source = "tectonic_pipeline.batch-export";
        options.include_update_count = true;
        options.update_count = current_update;
        platec::tooling::write_snapshot_bundle(params.output_dir / relative_dir, snapshot, options);
        records.push_back({label, current_update, snapshot.iteration_count, snapshot.cycle_count,
                           snapshot.time_myr, relative_dir.generic_string()});
        ++export_index;
    };

    if (params.include_initial) {
        export_snapshot("initial", 0);
    }
    while (explicit_index < explicit_targets.size() &&
           explicit_targets[explicit_index].update_count == 0) {
        export_snapshot(explicit_targets[explicit_index].label, 0);
        ++explicit_index;
    }

    while (!simulation.isFinished()) {
        simulation.update();
        ++update_count;

        if (params.checkpoint_every_updates != 0 &&
            update_count % params.checkpoint_every_updates == 0) {
            export_snapshot(checkpoint_label_from_update(update_count), update_count);
        }

        while (explicit_index < explicit_targets.size() &&
               explicit_targets[explicit_index].update_count == update_count) {
            export_snapshot(explicit_targets[explicit_index].label, update_count);
            ++explicit_index;
        }
    }

    if (params.include_final) {
        export_snapshot("final", update_count);
    }

    std::vector<RequestedCheckpoint> missing_targets;
    while (explicit_index < explicit_targets.size()) {
        missing_targets.push_back(explicit_targets[explicit_index]);
        ++explicit_index;
    }

    write_batch_manifest(params.output_dir / "batch_manifest.json", params, records, missing_targets);
    std::printf("batch export written to %s\n", params.output_dir.string().c_str());
    return 0;
}

int command_inspect(const std::vector<std::string>& args, bool boundary_stats_mode) {
    const InspectParams params = parse_inspect_args(args, boundary_stats_mode);
    uint32_t resolved_update_count = 0;
    const platec::contract::Snapshot snapshot =
        run_to_target(params.scenario, params.target, resolved_update_count);
    const platec::tooling::SnapshotSummary summary = platec::tooling::summarize_snapshot(snapshot);

    if (!boundary_stats_mode && params.has_bundle_output) {
        platec::tooling::SnapshotBundleOptions options;
        options.label = target_label(params.target, resolved_update_count, snapshot.time_myr);
        options.source = "tectonic_pipeline.inspect";
        options.include_update_count = true;
        options.update_count = resolved_update_count;
        platec::tooling::write_snapshot_bundle(params.bundle_output_dir, snapshot, options);
    }

    const std::string json = boundary_stats_mode
        ? platec::tooling::boundary_stats_json(summary)
        : platec::tooling::snapshot_summary_json(summary);
    write_output_or_stdout(json, params.has_output, params.output_path);
    return 0;
}

int command_compare(const std::vector<std::string>& args) {
    const CompareParams params = parse_compare_args(args);

    platec::tooling::SnapshotComparison comparison;
    if (params.bundle_mode) {
        const platec::tooling::LoadedSnapshotBundle bundle_a =
            platec::tooling::load_snapshot_bundle(params.bundle_a);
        const platec::tooling::LoadedSnapshotBundle bundle_b =
            platec::tooling::load_snapshot_bundle(params.bundle_b);
        comparison = platec::tooling::compare_snapshots(
            bundle_a.snapshot, bundle_b.snapshot, bundle_a.declared_boundary_segment_count,
            bundle_b.declared_boundary_segment_count, bundle_a.declared_deforming_region_count,
            bundle_b.declared_deforming_region_count,
            bundle_a.declared_tectonic_junction_count,
            bundle_b.declared_tectonic_junction_count);
    } else {
        platec::scenario::Scenario scenario_a = params.shared_scenario;
        platec::scenario::Scenario scenario_b = params.shared_scenario;
        scenario_a.seed = params.seed_a;
        scenario_b.seed = params.seed_b;

        uint32_t resolved_a = 0;
        uint32_t resolved_b = 0;
        const platec::contract::Snapshot snapshot_a =
            run_to_target(scenario_a, params.target, resolved_a);
        const platec::contract::Snapshot snapshot_b =
            run_to_target(scenario_b, params.target, resolved_b);
        if (resolved_a != resolved_b && params.target.mode != TargetMode::Final) {
            fail("live compare resolved different checkpoint updates");
        }
        comparison = platec::tooling::compare_snapshots(snapshot_a, snapshot_b);
    }

    write_output_or_stdout(platec::tooling::snapshot_comparison_json(comparison),
                           params.has_output, params.output_path);
    return 0;
}

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 2) {
            print_general_help();
            return 1;
        }

        const std::string command = argv[1];
        const std::vector<std::string> args = collect_args(argc, argv, 2);
        if (command == "--help" || command == "-h" || command == "help") {
            print_general_help();
            return 0;
        }
        if (command == "batch-export") {
            return command_batch_export(args);
        }
        if (command == "inspect") {
            return command_inspect(args, false);
        }
        if (command == "boundary-stats") {
            return command_inspect(args, true);
        }
        if (command == "compare") {
            return command_compare(args);
        }

        fail("unknown command: " + command);
    } catch (const std::exception& exception) {
        std::fprintf(stderr, "error: %s\n", exception.what());
        return 1;
    }
}
