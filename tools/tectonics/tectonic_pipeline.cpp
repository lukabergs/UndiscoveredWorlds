#ifndef NOMINMAX
#define NOMINMAX
#endif

#include "tectonic_api.hpp"
#include "mineral_history.hpp"
#include "tectonic_frame_reader.hpp"
#include "tectonic_pipeline_support.hpp"
#include "tectonic_recipe.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <optional>
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
    tectonics::scenario::Scenario scenario = tectonics::scenario::baseline();
    fs::path output_dir;
    bool include_initial = true;
    bool include_final = true;
    uint32_t checkpoint_every_updates = 0;
    std::vector<uint32_t> explicit_updates;
    std::vector<double> explicit_times_myr;
    bool mineral_history = false;
};

struct InspectParams {
    tectonics::scenario::Scenario scenario = tectonics::scenario::baseline();
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
    tectonics::scenario::Scenario shared_scenario = tectonics::scenario::baseline();
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

tectonics::scenario::Scenario default_scenario() {
    tectonics::scenario::Scenario scenario = tectonics::scenario::baseline();
    scenario.seed = 1;
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
    std::printf("  --plates N\n");
    std::printf("  --erosion-strength X\n");
    std::printf("  --thermal-subsidence-strength X (0 disables conductive plate cooling)\n");
    std::printf("  --sea-level-m N\n");
    std::printf("  --initial-ocean-depth-m N (positive depth below sea level)\n");
    std::printf("  --initial-land-relief-m N (relief above sea level)\n");
    std::printf("  --evolution-duration-myr X (continuous duration)\n");
    std::printf("  --world-width-km X --mantle-flow-speed-cm-yr X\n");
    std::printf("  --transport-subdivision 1|2|4 --deformation-peak-km X (0 disables)\n");
    std::printf("  --deformation-variation X (0..0.5) --subduction-coupling X (0..1)\n");
    std::printf("  --coastal-gpe-strength X (0..1; mountain-root resistance at active margins)\n");
    std::printf("  --boundary-force-ramp-cm-yr X (nonnegative; 0 retains abrupt force activation)\n");
    std::printf("  --consistent-motion on|off (nonlinear activity balance and shared slab polarity)\n");
    std::printf("  --boundary-memory-myr X (positive; attached boundary activity response)\n");
    std::printf("  --mountain-force-range-km X (nonnegative; 0 samples only contact-cell roots)\n");
    std::printf("  --ocean-crust-variation-m X (0..7000 exclusive) --ocean-noise-wavelength-km X\n");
    std::printf("  --seamount-crust-m X --seamount-radius-km X --seamount-spacing-km X\n");
    std::printf("  --submarine-erosion-ratio X (0..1; relative to land transport)\n");
    std::printf("  --seamount-roughness X (0..1; coherent shape modulation)\n");
    std::printf("  --margin-processes on|off --trench-load-n-m X --flexural-length-km X\n");
    std::printf("  --arc-distance-km X --arc-melt-ratio X --river-erosion X --sediment-travel-km X\n");
    std::printf("  --thermal-age-rate X\n");
    std::printf("  --thermal-force-feedback on|off; --initial-boundary-influence X\n");
    std::printf("  --wave-erosion-m-myr X\n");
    std::printf("  --wave-retreat-km-myr X (nonnegative; subgrid coastal cliff retreat)\n");
    std::printf("  --alluvial-transport X (nonnegative; 0 retains land sediment bypass)\n");
    std::printf("  --river-outlet-sea-level on|off (independent of alluvial deposition; default off)\n");
    std::printf("  --wave-height-scale-m X\n");
    std::printf("  --root-mobility X\n");
    std::printf("  --progressive-rifting on|off; --rift-yield-n-m X; --rift-damage-time-myr X; --rift-strain-rate-myr X\n");
    std::printf("  --rift-breakup-thickness-m X (onset thickness for oceanization during extension)\n");
    std::printf("  --boundary-conditioned-terrain on|off (preserve seeded relief and condition opening contacts)\n");
    std::printf("  --initial-margin-width-m X (0..4000; smooth initial crust mixture around -2500 m)\n");
    std::printf("  --sediment-load-feedback on|off (continental sediment contributes to root forces and flow)\n");
    std::printf("  --noise-pass-amplitude-m X (nonnegative; synthetic conservative detail experiment)\n");
    std::printf("  --noise-pass-wavelength-km X --noise-pass-period-myr X (positive)\n");
    std::printf("  --noise-pass-coast-taper-m X (nonnegative; smooth suppression near sea level, default 500 m)\n");
    std::printf("  --plume-count N (0..64; mantle plumes fixed in the mantle frame, default 0 = off)\n");
    std::printf("  --plume-tail-flux-km3-yr X --plume-tail-radius-km X (basalt per plume after birth; Gaussian sigma)\n");
    std::printf("  --plume-head-volume-km3 X --plume-head-radius-km X --plume-head-duration-myr X (head pulse at birth; 0 = none)\n");
    std::printf("  --plume-birth-window-myr X (births spread over [0, X]) --plume-min-spacing-km X\n");
    std::printf("  --plume-swell-age-myr X --plume-swell-radius-km X (reheat ocean floor toward X Myr; 0 = none)\n");
    std::printf("  --plume-weakening X (0..1 rock-strength reduction; needs --progressive-rifting on)\n");
    std::printf("  --trench-response-myr X (time constant of local bending response)\n");
    std::printf("  --process-variation X --process-wavelength-km X --process-period-myr X\n");
    std::printf("  --delta-time-myr X\n");
}

void print_batch_help() {
    std::printf("tectonic_pipeline batch-export --output DIR [options]\n");
    std::printf("  --checkpoint-every-updates N\n");
    std::printf("  --checkpoint-at-update N\n");
    std::printf("  --checkpoint-at-time-myr X\n");
    std::printf("  --no-initial\n");
    std::printf("  --no-final\n");
    std::printf("  --mineral-history on|off (gold deposits: mineral_deposits.csv, mineral_budget.json)\n");
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
    if (!std::isfinite(time_myr) || time_myr < 0.0 ||
        !std::isfinite(delta_time_myr) || delta_time_myr <= 0.0) {
        fail("time must be finite and >= 0; --delta-time-myr must be finite and > 0");
    }
    const double raw = time_myr / delta_time_myr;
    if (!std::isfinite(raw) || raw > std::numeric_limits<uint32_t>::max()) {
        fail("requested time exceeds the supported update count");
    }
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

void parse_common_scenario_arg(tectonics::scenario::Scenario& scenario, const std::string& arg,
                               const char* value) {
    if (arg == "--seed") {
        scenario.seed = parse_long("--seed", value);
    } else if (arg == "--width") {
        scenario.width = parse_u32("--width", value);
    } else if (arg == "--height") {
        scenario.height = parse_u32("--height", value);
    } else if (arg == "--sea-level") {
        scenario.sea_level = parse_float("--sea-level", value);
    } else if (arg == "--plates") {
        scenario.plate_count = parse_u32("--plates", value);
    } else if (arg == "--erosion-strength") {
        scenario.erosion_strength = parse_float("--erosion-strength", value);
    } else if (arg == "--thermal-subsidence-strength") {
        scenario.thermal_subsidence_strength = parse_float("--thermal-subsidence-strength", value);
    } else if (arg == "--sea-level-m") {
        scenario.sea_level_m = parse_i32("--sea-level-m", value);
    } else if (arg == "--initial-ocean-depth-m") {
        scenario.initial_ocean_depth_m = parse_u16("--initial-ocean-depth-m", value);
    } else if (arg == "--initial-land-relief-m") {
        scenario.initial_land_relief_m = parse_u16("--initial-land-relief-m", value);
    } else if (arg == "--evolution-duration-myr") {
        scenario.evolution_duration_myr = parse_double("--evolution-duration-myr", value);
    } else if (arg == "--world-width-km") {
        scenario.world_width_km = parse_double("--world-width-km", value);
    } else if (arg == "--mantle-flow-speed-cm-yr") {
        scenario.mantle_flow_speed_cm_yr = parse_double("--mantle-flow-speed-cm-yr", value);
    } else if (arg == "--transport-subdivision") {
        scenario.transport_subdivision = parse_u32("--transport-subdivision", value);
    } else if (arg == "--deformation-peak-km") {
        scenario.deformation_peak_km = parse_double("--deformation-peak-km", value);
    } else if (arg == "--deformation-variation") {
        scenario.deformation_variation = parse_double("--deformation-variation", value);
    } else if (arg == "--subduction-coupling") {
        scenario.subduction_coupling = parse_double("--subduction-coupling", value);
    } else if (arg == "--coastal-gpe-strength") {
        scenario.coastal_gpe_strength = parse_double("--coastal-gpe-strength", value);
    } else if (arg == "--boundary-force-ramp-cm-yr") {
        scenario.boundary_force_ramp_cm_yr = parse_double("--boundary-force-ramp-cm-yr", value);
    } else if (arg == "--consistent-motion") {
        if (std::string(value) != "on" && std::string(value) != "off") throw std::runtime_error("consistent-motion expects on or off");
        scenario.consistent_motion = std::string(value) == "on";
    } else if (arg == "--boundary-memory-myr") {
        scenario.boundary_memory_myr=parse_double(arg.c_str(),value);
    } else if (arg == "--mountain-force-range-km") {
        scenario.mountain_force_range_km = parse_double(arg.c_str(), value);
    } else if (arg == "--ocean-crust-variation-m") {
        scenario.ocean_crust_variation_m = parse_double("--ocean-crust-variation-m", value);
    } else if (arg == "--ocean-noise-wavelength-km") {
        scenario.ocean_noise_wavelength_km = parse_double("--ocean-noise-wavelength-km", value);
    } else if (arg == "--seamount-crust-m") {
        scenario.seamount_crust_m = parse_double("--seamount-crust-m", value);
    } else if (arg == "--seamount-radius-km") {
        scenario.seamount_radius_km = parse_double("--seamount-radius-km", value);
    } else if (arg == "--seamount-spacing-km") {
        scenario.seamount_spacing_km = parse_double("--seamount-spacing-km", value);
    } else if (arg == "--submarine-erosion-ratio") {
        scenario.submarine_erosion_ratio = parse_double("--submarine-erosion-ratio", value);
    } else if (arg == "--seamount-roughness") {
        scenario.seamount_roughness = parse_double("--seamount-roughness", value);
    } else if (arg == "--boundary-conditioned-terrain") {
        const std::string mode(value);
        if (mode!="on" && mode!="off") fail("--boundary-conditioned-terrain requires on or off");
        scenario.boundary_conditioned_terrain=mode=="on";
    } else if (arg == "--noise-pass-amplitude-m") {
        scenario.experimental_noise.amplitude_m=parse_double(arg.c_str(),value);
    } else if (arg == "--noise-pass-wavelength-km") {
        scenario.experimental_noise.wavelength_km=parse_double(arg.c_str(),value);
    } else if (arg == "--noise-pass-period-myr") {
        scenario.experimental_noise.period_myr=parse_double(arg.c_str(),value);
    } else if (arg == "--noise-pass-coast-taper-m") {
        scenario.experimental_noise.coast_taper_m=parse_double(arg.c_str(),value);
    } else if (arg == "--plume-count") {
        scenario.plume_processes.count=parse_u32("--plume-count",value);
    } else if (arg == "--plume-tail-flux-km3-yr") {
        scenario.plume_processes.tail_flux_km3_yr=parse_double(arg.c_str(),value);
    } else if (arg == "--plume-tail-radius-km") {
        scenario.plume_processes.tail_radius_km=parse_double(arg.c_str(),value);
    } else if (arg == "--plume-head-volume-km3") {
        scenario.plume_processes.head_volume_km3=parse_double(arg.c_str(),value);
    } else if (arg == "--plume-head-radius-km") {
        scenario.plume_processes.head_radius_km=parse_double(arg.c_str(),value);
    } else if (arg == "--plume-head-duration-myr") {
        scenario.plume_processes.head_duration_myr=parse_double(arg.c_str(),value);
    } else if (arg == "--plume-birth-window-myr") {
        scenario.plume_processes.birth_window_myr=parse_double(arg.c_str(),value);
    } else if (arg == "--plume-swell-age-myr") {
        scenario.plume_processes.swell_age_myr=parse_double(arg.c_str(),value);
    } else if (arg == "--plume-swell-radius-km") {
        scenario.plume_processes.swell_radius_km=parse_double(arg.c_str(),value);
    } else if (arg == "--plume-weakening") {
        scenario.plume_processes.weakening=parse_double(arg.c_str(),value);
    } else if (arg == "--plume-min-spacing-km") {
        scenario.plume_processes.min_spacing_km=parse_double(arg.c_str(),value);
    } else if (arg == "--sediment-load-feedback") {
        const std::string mode(value);
        if (mode!="on" && mode!="off") fail("--sediment-load-feedback requires on or off");
        scenario.sediment_load_feedback=mode=="on";
    } else if (arg == "--initial-margin-width-m") {
        scenario.initial_margin_width_m=parse_double(arg.c_str(),value);
    } else if (arg == "--margin-processes") {
        const std::string mode(value);
        if (mode!="on" && mode!="off") fail("--margin-processes requires on or off");
        scenario.margin_processes.enabled=mode=="on";
    } else if (arg == "--trench-load-n-m") {
        scenario.margin_processes.trench_load_n_m=parse_double(arg.c_str(),value);
    } else if (arg == "--flexural-length-km") {
        scenario.margin_processes.flexural_length_km=parse_double(arg.c_str(),value);
    } else if (arg == "--trench-response-myr") {
        scenario.margin_processes.trench_response_myr=parse_double(arg.c_str(),value);
    } else if (arg == "--process-variation") {
        scenario.margin_processes.process_variation=parse_double(arg.c_str(),value);
    } else if (arg == "--process-wavelength-km") {
        scenario.margin_processes.process_wavelength_km=parse_double(arg.c_str(),value);
    } else if (arg == "--process-period-myr") {
        scenario.margin_processes.process_period_myr=parse_double(arg.c_str(),value);
    } else if (arg == "--arc-distance-km") {
        scenario.margin_processes.arc_distance_km=parse_double(arg.c_str(),value);
    } else if (arg == "--arc-melt-ratio") {
        scenario.margin_processes.arc_melt_ratio=parse_double(arg.c_str(),value);
    } else if (arg == "--thermal-age-rate") {
        scenario.margin_processes.thermal_age_rate=parse_double(arg.c_str(),value);
    } else if (arg == "--thermal-force-feedback") {
        if (std::string(value)!="on" && std::string(value)!="off") fail("thermal-force-feedback requires on/off");
        scenario.thermal_force_feedback=std::string(value)=="on";
    } else if (arg == "--wave-erosion-m-myr") {
        scenario.margin_processes.wave_erosion_m_myr=parse_double(arg.c_str(),value);
    } else if (arg == "--wave-height-scale-m") {
        scenario.margin_processes.wave_height_scale_m=parse_double(arg.c_str(),value);
    } else if (arg == "--wave-retreat-km-myr") {
        scenario.margin_processes.wave_retreat_km_myr=parse_double(arg.c_str(),value);
    } else if (arg == "--root-mobility") {
        scenario.margin_processes.root_mobility=parse_double(arg.c_str(),value);
    } else if (arg == "--progressive-rifting") {
        if (std::string(value)!="on" && std::string(value)!="off") fail("progressive-rifting requires on/off");
        scenario.topology_processes.progressive_rifting=std::string(value)=="on";
    } else if (arg == "--rift-yield-n-m") {
        scenario.topology_processes.yield_n_m=parse_double(arg.c_str(),value);
    } else if (arg == "--rift-damage-time-myr") {
        scenario.topology_processes.damage_time_myr=parse_double(arg.c_str(),value);
    } else if (arg == "--rift-strain-rate-myr") {
        scenario.topology_processes.strain_rate_myr=parse_double(arg.c_str(),value);
    } else if (arg == "--sediment-subduction-fraction") {
        scenario.margin_processes.sediment_subduction_fraction=parse_double(arg.c_str(),value);
    } else if (arg == "--initial-boundary-influence") {
        scenario.initial_boundary_influence=parse_double(arg.c_str(),value);
    } else if (arg == "--river-erosion") {
        scenario.margin_processes.river_erosion=parse_double(arg.c_str(),value);
    } else if (arg == "--alluvial-transport") {
        scenario.margin_processes.alluvial_transport=parse_double(arg.c_str(),value);
    } else if (arg == "--river-outlet-sea-level") {
        if (std::string(value)!="on" && std::string(value)!="off") fail("river-outlet-sea-level requires on/off");
        scenario.margin_processes.river_outlet_sea_level=std::string(value)=="on";
    } else if (arg == "--sediment-travel-km") {
        scenario.margin_processes.sediment_travel_km=parse_double(arg.c_str(),value);
    } else if (arg == "--rift-breakup-thickness-m") {
        scenario.margin_processes.rift_breakup_thickness_m=parse_double(arg.c_str(),value);
    } else if (arg == "--delta-time-myr") {
        scenario.delta_time_myr = parse_double("--delta-time-myr", value);
    } else {
        fail("unknown argument: " + arg);
    }
}

TargetSelector parse_target_selector(const std::vector<std::string>& args, size_t& index) {
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
            // Resolve after all scenario flags, including a later timestep flag.
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

// Owns one engine instance through the public tectonics API.
class Simulation {
public:
    explicit Simulation(const tectonics::scenario::Scenario& scenario)
        : handle_(tectonics_create_from_scenario(&scenario)) {
        if (handle_ == nullptr) throw std::runtime_error("Failed to create tectonics simulation");
    }
    ~Simulation() { tectonics_destroy(handle_); }
    Simulation(const Simulation&) = delete;
    Simulation& operator=(const Simulation&) = delete;

    void update() { tectonics_step(handle_); }
    bool isFinished() const { return tectonics_is_finished(handle_) != 0; }
    double getTimeMyr() const { return tectonics_get_time_myr(handle_); }
    tectonics::contract::Snapshot snapshot() const { return tectonics_capture_snapshot(handle_); }
    void* handle() const { return handle_; }

private:
    void* handle_;
};

tectonics::contract::Snapshot run_to_target(const tectonics::scenario::Scenario& scenario,
                                         const TargetSelector& target,
                                         uint32_t& resolved_update_count) {
    Simulation simulation(scenario);
    resolved_update_count = 0;

    if (target.mode == TargetMode::Final) {
        while (!simulation.isFinished()) {
            simulation.update();
            ++resolved_update_count;
        }
        return simulation.snapshot();
    }

    const auto target_pending = [&]() {
        return target.came_from_time
            ? simulation.getTimeMyr() + 1e-9 < target.requested_time_myr
            : resolved_update_count < target.update_count;
    };
    while (target_pending() && !simulation.isFinished()) {
        simulation.update();
        ++resolved_update_count;
    }
    if (target_pending()) {
        std::ostringstream message;
        message << "simulation finished at update " << resolved_update_count << " (time "
                << simulation.getTimeMyr() << " Myr) before requested target "
                << (target.came_from_time ? target.requested_time_myr : target.update_count)
                << (target.came_from_time ? " Myr" : " updates");
        fail(message.str());
    }
    return simulation.snapshot();
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
    output << "    \"plate_count\": " << params.scenario.plate_count << ",\n";
    output << "    \"erosion_strength\": " << params.scenario.erosion_strength << ",\n";
    output << "    \"sea_level_m\": " << params.scenario.sea_level_m << ",\n";
    output << "    \"initial_ocean_depth_m\": " << params.scenario.initial_ocean_depth_m << ",\n";
    output << "    \"initial_land_relief_m\": " << params.scenario.initial_land_relief_m << ",\n";
    output << "    \"evolution_duration_myr\": " << params.scenario.evolution_duration_myr << ",\n";
    output << "    \"motion_model\": \"" << "force-balance" << "\",\n";
    output << "    \"surface_model\": \"" << "isostatic-columns" << "\",\n";
    output << "    \"transport_subdivision\": " << params.scenario.transport_subdivision << ",\n";
    output << "    \"deformation_peak_km\": " << params.scenario.deformation_peak_km << ",\n";
    output << "    \"deformation_variation\": " << params.scenario.deformation_variation << ",\n";
    output << "    \"subduction_coupling\": " << params.scenario.subduction_coupling << ",\n";
    output << "    \"coastal_gpe_strength\": " << params.scenario.coastal_gpe_strength << ",\n";
    output << "    \"boundary_force_ramp_cm_yr\": " << params.scenario.boundary_force_ramp_cm_yr << ",\n";
    output << "    \"mountain_force_range_km\": " << params.scenario.mountain_force_range_km << ",\n";
    output << "    \"initial_margin_width_m\": " << params.scenario.initial_margin_width_m << ",\n";
    output << "    \"initial_boundary_influence\": " << params.scenario.initial_boundary_influence << ",\n";
    output << "    \"consistent_motion\": " << (params.scenario.consistent_motion?"true":"false") << ",\n";
    output << "    \"boundary_memory_myr\": " << params.scenario.boundary_memory_myr << ",\n";
    output << "    \"noise_pass_amplitude_m\": " << params.scenario.experimental_noise.amplitude_m << ",\n";
    output << "    \"noise_pass_wavelength_km\": " << params.scenario.experimental_noise.wavelength_km << ",\n";
    output << "    \"noise_pass_period_myr\": " << params.scenario.experimental_noise.period_myr << ",\n";
    output << "    \"noise_pass_coast_taper_m\": " << params.scenario.experimental_noise.coast_taper_m << ",\n";
    const auto& plume = params.scenario.plume_processes;
    output << "    \"plume_count\": " << plume.count << ",\n";
    output << "    \"plume_tail_flux_km3_yr\": " << plume.tail_flux_km3_yr << ",\n";
    output << "    \"plume_tail_radius_km\": " << plume.tail_radius_km << ",\n";
    output << "    \"plume_head_volume_km3\": " << plume.head_volume_km3 << ",\n";
    output << "    \"plume_head_radius_km\": " << plume.head_radius_km << ",\n";
    output << "    \"plume_head_duration_myr\": " << plume.head_duration_myr << ",\n";
    output << "    \"plume_birth_window_myr\": " << plume.birth_window_myr << ",\n";
    output << "    \"plume_swell_age_myr\": " << plume.swell_age_myr << ",\n";
    output << "    \"plume_swell_radius_km\": " << plume.swell_radius_km << ",\n";
    output << "    \"plume_weakening\": " << plume.weakening << ",\n";
    output << "    \"plume_min_spacing_km\": " << plume.min_spacing_km << ",\n";
    output << "    \"ocean_crust_variation_m\": " << params.scenario.ocean_crust_variation_m << ",\n";
    output << "    \"ocean_noise_wavelength_km\": " << params.scenario.ocean_noise_wavelength_km << ",\n";
    output << "    \"seamount_crust_m\": " << params.scenario.seamount_crust_m << ",\n";
    output << "    \"seamount_radius_km\": " << params.scenario.seamount_radius_km << ",\n";
    output << "    \"seamount_spacing_km\": " << params.scenario.seamount_spacing_km << ",\n";
    output << "    \"seamount_roughness\": " << params.scenario.seamount_roughness << ",\n";
    output << "    \"submarine_erosion_ratio\": " << params.scenario.submarine_erosion_ratio << ",\n";
    output << "    \"margin_processes\": " << (params.scenario.margin_processes.enabled ? "true" : "false") << ",\n";
    output << "    \"boundary_conditioned_terrain\": " << (params.scenario.boundary_conditioned_terrain ? "true" : "false") << ",\n";
    output << "    \"sediment_load_feedback\": " << (params.scenario.sediment_load_feedback ? "true" : "false") << ",\n";
    output << "    \"trench_load_n_m\": " << params.scenario.margin_processes.trench_load_n_m << ",\n";
    output << "    \"flexural_length_km\": " << params.scenario.margin_processes.flexural_length_km << ",\n";
    output << "    \"trench_response_myr\": " << params.scenario.margin_processes.trench_response_myr << ",\n";
    output << "    \"process_variation\": " << params.scenario.margin_processes.process_variation << ",\n";
    output << "    \"process_wavelength_km\": " << params.scenario.margin_processes.process_wavelength_km << ",\n";
    output << "    \"process_period_myr\": " << params.scenario.margin_processes.process_period_myr << ",\n";
    output << "    \"arc_distance_km\": " << params.scenario.margin_processes.arc_distance_km << ",\n";
    output << "    \"arc_melt_ratio\": " << params.scenario.margin_processes.arc_melt_ratio << ",\n";
    output << "    \"river_erosion\": " << params.scenario.margin_processes.river_erosion << ",\n";
    output << "    \"sediment_subduction_fraction\": " << params.scenario.margin_processes.sediment_subduction_fraction << ",\n";
    output << "    \"progressive_rifting\": " << params.scenario.topology_processes.progressive_rifting << ",\n";
    output << "    \"rift_yield_n_m\": " << params.scenario.topology_processes.yield_n_m << ",\n";
    output << "    \"rift_damage_time_myr\": " << params.scenario.topology_processes.damage_time_myr << ",\n";
    output << "    \"rift_strain_rate_myr\": " << params.scenario.topology_processes.strain_rate_myr << ",\n";
    output << "    \"thermal_age_rate\": " << params.scenario.margin_processes.thermal_age_rate << ",\n";
    output << "    \"thermal_force_feedback\": " << (params.scenario.thermal_force_feedback?"true":"false") << ",\n";
    output << "    \"wave_erosion_m_myr\": " << params.scenario.margin_processes.wave_erosion_m_myr << ",\n";
    output << "    \"wave_height_scale_m\": " << params.scenario.margin_processes.wave_height_scale_m << ",\n";
    output << "    \"alluvial_transport\": " << params.scenario.margin_processes.alluvial_transport << ",\n";
    output << "    \"river_outlet_sea_level\": " << (params.scenario.margin_processes.river_outlet_sea_level?"true":"false") << ",\n";
    output << "    \"wave_retreat_km_myr\": " << params.scenario.margin_processes.wave_retreat_km_myr << ",\n";
    output << "    \"root_mobility\": " << params.scenario.margin_processes.root_mobility << ",\n";
    output << "    \"sediment_travel_km\": " << params.scenario.margin_processes.sediment_travel_km << ",\n";
    output << "    \"rift_breakup_thickness_m\": " << params.scenario.margin_processes.rift_breakup_thickness_m << ",\n";
    output << "    \"world_width_km\": " << params.scenario.world_width_km << ",\n";
    output << "    \"mantle_flow_speed_cm_yr\": " << params.scenario.mantle_flow_speed_cm_yr << ",\n";
    output << "    \"thermal_subsidence_strength\": " << params.scenario.thermal_subsidence_strength << ",\n";
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
        } else if (arg == "--mineral-history") {
            if (std::string(value) != "on" && std::string(value) != "off")
                fail("--mineral-history requires on or off");
            params.mineral_history = std::string(value) == "on";
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
            params.target = parse_target_selector(args, index);
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
            params.target = parse_target_selector(args, index);
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
    fs::remove(params.output_dir / "mineral_deposits.csv", ec);
    fs::remove(params.output_dir / "mineral_budget.json", ec);
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

    Simulation simulation(params.scenario);
    uint32_t update_count = 0;
    size_t export_index = 0;
    std::vector<ExportedSnapshotRecord> records;
    mineralsystems::TectonicFrameReader mineral_frames(params.scenario.world_width_km);
    std::optional<mineralsystems::MineralHistory> minerals;
    if (params.mineral_history) {
        minerals.emplace(mineralsystems::MineralParameters{}, mineralsystems::mineralseed(params.scenario.seed));
        minerals->begin(mineral_frames.read(simulation.handle()));
    }

    const auto export_snapshot = [&](const std::string& label, uint32_t current_update) {
        const tectonics::contract::Snapshot snapshot = simulation.snapshot();
        const fs::path relative_dir = fs::path("snapshots") / snapshot_directory_name(export_index, label);
        tectonics::tooling::SnapshotBundleOptions options;
        options.label = label;
        options.source = "tectonic_pipeline.batch-export";
        options.include_update_count = true;
        options.update_count = current_update;
        tectonics::tooling::write_snapshot_bundle(params.output_dir / relative_dir, snapshot, options);
        records.push_back({label, current_update, snapshot.iteration_count,
                           snapshot.time_myr, relative_dir.generic_string()});
        ++export_index;
    };

    if (params.include_initial) {
        export_snapshot("initial", 0);
    }
    const auto export_reached_targets = [&]() {
        for (auto it = explicit_targets.begin(); it != explicit_targets.end();) {
            const bool reached = it->from_time
                ? simulation.getTimeMyr() + 1e-9 >= it->requested_time_myr
                : update_count >= it->update_count;
            if (reached) {
                export_snapshot(it->label, update_count);
                it = explicit_targets.erase(it);
            } else {
                ++it;
            }
        }
    };
    export_reached_targets();

    while (!simulation.isFinished()) {
        simulation.update();
        ++update_count;
        if (minerals) {
            minerals->step(mineral_frames.read(simulation.handle()));
        }

        if (params.checkpoint_every_updates != 0 &&
            update_count % params.checkpoint_every_updates == 0) {
            export_snapshot(checkpoint_label_from_update(update_count), update_count);
        }

        export_reached_targets();
    }

    if (params.include_final) {
        export_snapshot("final", update_count);
    }

    write_batch_manifest(params.output_dir / "batch_manifest.json", params, records, explicit_targets);
    if (minerals) {
        minerals->finish(mineral_frames.read(simulation.handle()));
        std::ofstream deposits(params.output_dir / "mineral_deposits.csv");
        mineralsystems::writedepositscsv(deposits, minerals->deposits());
        std::ofstream budget(params.output_dir / "mineral_budget.json");
        budget << mineralsystems::mineralbudgetjson(minerals->budget());
        if (!deposits || !budget) {
            fail("failed to write mineral history outputs");
        }
    }
    std::printf("batch export written to %s\n", params.output_dir.string().c_str());
    return 0;
}

int command_inspect(const std::vector<std::string>& args, bool boundary_stats_mode) {
    const InspectParams params = parse_inspect_args(args, boundary_stats_mode);
    uint32_t resolved_update_count = 0;
    const tectonics::contract::Snapshot snapshot =
        run_to_target(params.scenario, params.target, resolved_update_count);
    const tectonics::tooling::SnapshotSummary summary = tectonics::tooling::summarize_snapshot(snapshot);

    if (!boundary_stats_mode && params.has_bundle_output) {
        tectonics::tooling::SnapshotBundleOptions options;
        options.label = target_label(params.target, resolved_update_count, snapshot.time_myr);
        options.source = "tectonic_pipeline.inspect";
        options.include_update_count = true;
        options.update_count = resolved_update_count;
        tectonics::tooling::write_snapshot_bundle(params.bundle_output_dir, snapshot, options);
    }

    const std::string json = boundary_stats_mode
        ? tectonics::tooling::boundary_stats_json(summary)
        : tectonics::tooling::snapshot_summary_json(summary);
    write_output_or_stdout(json, params.has_output, params.output_path);
    return 0;
}

int command_compare(const std::vector<std::string>& args) {
    const CompareParams params = parse_compare_args(args);

    tectonics::tooling::SnapshotComparison comparison;
    if (params.bundle_mode) {
        const tectonics::tooling::LoadedSnapshotBundle bundle_a =
            tectonics::tooling::load_snapshot_bundle(params.bundle_a);
        const tectonics::tooling::LoadedSnapshotBundle bundle_b =
            tectonics::tooling::load_snapshot_bundle(params.bundle_b);
        comparison = tectonics::tooling::compare_snapshots(
            bundle_a.snapshot, bundle_b.snapshot, bundle_a.declared_boundary_segment_count,
            bundle_b.declared_boundary_segment_count, bundle_a.declared_deforming_region_count,
            bundle_b.declared_deforming_region_count,
            bundle_a.declared_tectonic_junction_count,
            bundle_b.declared_tectonic_junction_count);
    } else {
        tectonics::scenario::Scenario scenario_a = params.shared_scenario;
        tectonics::scenario::Scenario scenario_b = params.shared_scenario;
        scenario_a.seed = params.seed_a;
        scenario_b.seed = params.seed_b;

        uint32_t resolved_a = 0;
        uint32_t resolved_b = 0;
        const tectonics::contract::Snapshot snapshot_a =
            run_to_target(scenario_a, params.target, resolved_a);
        const tectonics::contract::Snapshot snapshot_b =
            run_to_target(scenario_b, params.target, resolved_b);
        if (resolved_a != resolved_b && params.target.mode != TargetMode::Final) {
            fail("live compare resolved different checkpoint updates");
        }
        comparison = tectonics::tooling::compare_snapshots(snapshot_a, snapshot_b);
    }

    write_output_or_stdout(tectonics::tooling::snapshot_comparison_json(comparison),
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
