// Scenario options of the tectonics engine: every key of configs/tectonics-baseline.json
// (and the engine's other scenario members) as --key value for world_pipeline;
// range checks stay in the engine.

#include "cli/tectonic_options.hpp"

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <stdexcept>

namespace
{
[[noreturn]] void fail(const std::string& message)
{
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

bool parse_scenario_arg(tectonics::scenario::Scenario& scenario, const std::string& arg,
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
    } else if (arg == "--river-implicit") {
        if (std::string(value)!="on" && std::string(value)!="off") fail("river-implicit requires on/off");
        scenario.margin_processes.river_implicit=std::string(value)=="on";
    } else if (arg == "--alluvial-transport") {
        scenario.margin_processes.alluvial_transport=parse_double(arg.c_str(),value);
    } else if (arg == "--river-outlet-sea-level") {
        if (std::string(value)!="on" && std::string(value)!="off") fail("river-outlet-sea-level requires on/off");
        scenario.margin_processes.river_outlet_sea_level=std::string(value)=="on";
    } else if (arg == "--sediment-travel-km") {
        scenario.margin_processes.sediment_travel_km=parse_double(arg.c_str(),value);
    } else if (arg == "--rift-breakup-thickness-m") {
        scenario.margin_processes.rift_breakup_thickness_m=parse_double(arg.c_str(),value);
    } else if (arg == "--continental-fraction") {
        scenario.continental_fraction=parse_double(arg.c_str(),value);
    } else if (arg == "--continental-root-kg-m3") {
        scenario.continental_root_kg_m3=parse_double(arg.c_str(),value);
    } else if (arg == "--continental-min-thickness-m") {
        scenario.margin_processes.continental_min_thickness_m=parse_double(arg.c_str(),value);
    } else if (arg == "--continental-max-thickness-m") {
        scenario.margin_processes.continental_max_thickness_m=parse_double(arg.c_str(),value);
    } else if (arg == "--orogenic-collapse-km2-myr") {
        scenario.margin_processes.orogenic_collapse_km2_myr=parse_double(arg.c_str(),value);
    } else if (arg == "--arc-continent-fraction") {
        scenario.margin_processes.arc_continent_fraction=parse_double(arg.c_str(),value);
    } else if (arg == "--subduction-erosion-ratio") {
        scenario.margin_processes.subduction_erosion_ratio=parse_double(arg.c_str(),value);
    } else if (arg == "--continental-growth-control") {
        if (std::string(value)!="on" && std::string(value)!="off") fail("continental-growth-control requires on/off");
        scenario.margin_processes.continental_growth_control=std::string(value)=="on";
    } else if (arg == "--continental-growth-per-gyr") {
        scenario.margin_processes.continental_growth_per_gyr=parse_double(arg.c_str(),value);
    } else if (arg == "--delta-time-myr") {
        scenario.delta_time_myr = parse_double("--delta-time-myr", value);
    } else if (arg == "--uplift-rate-window-myr") {
        scenario.uplift_rate_window_myr=parse_double(arg.c_str(),value);
    } else {
        return false;
    }
    return true;
}
}

bool parsetectonicoption(tectonics::scenario::Scenario& scenario, const std::string& argument, const std::string& value)
{
    // The grid and the seed belong to the run (--resolution, --seed), not to the scenario.
    if (argument == "--seed" || argument == "--width" || argument == "--height")
        return false;
    return parse_scenario_arg(scenario, argument, value.c_str());
}

void printtectonicoptionsusage(std::ostream& out)
{
    out << "  Tectonic scenario options (--key value, defaults from configs/tectonics-baseline.json):\n";
    const char* const flags[] = {
        "--sea-level",
        "--plates",
        "--erosion-strength",
        "--thermal-subsidence-strength",
        "--sea-level-m",
        "--initial-ocean-depth-m",
        "--initial-land-relief-m",
        "--evolution-duration-myr",
        "--world-width-km",
        "--mantle-flow-speed-cm-yr",
        "--transport-subdivision",
        "--deformation-peak-km",
        "--deformation-variation",
        "--subduction-coupling",
        "--coastal-gpe-strength",
        "--boundary-force-ramp-cm-yr",
        "--consistent-motion",
        "--boundary-memory-myr",
        "--mountain-force-range-km",
        "--ocean-crust-variation-m",
        "--ocean-noise-wavelength-km",
        "--seamount-crust-m",
        "--seamount-radius-km",
        "--seamount-spacing-km",
        "--submarine-erosion-ratio",
        "--seamount-roughness",
        "--boundary-conditioned-terrain",
        "--noise-pass-amplitude-m",
        "--noise-pass-wavelength-km",
        "--noise-pass-period-myr",
        "--noise-pass-coast-taper-m",
        "--plume-count",
        "--plume-tail-flux-km3-yr",
        "--plume-tail-radius-km",
        "--plume-head-volume-km3",
        "--plume-head-radius-km",
        "--plume-head-duration-myr",
        "--plume-birth-window-myr",
        "--plume-swell-age-myr",
        "--plume-swell-radius-km",
        "--plume-weakening",
        "--plume-min-spacing-km",
        "--sediment-load-feedback",
        "--initial-margin-width-m",
        "--margin-processes",
        "--trench-load-n-m",
        "--flexural-length-km",
        "--trench-response-myr",
        "--process-variation",
        "--process-wavelength-km",
        "--process-period-myr",
        "--arc-distance-km",
        "--arc-melt-ratio",
        "--thermal-age-rate",
        "--thermal-force-feedback",
        "--wave-erosion-m-myr",
        "--wave-height-scale-m",
        "--wave-retreat-km-myr",
        "--root-mobility",
        "--progressive-rifting",
        "--rift-yield-n-m",
        "--rift-damage-time-myr",
        "--rift-strain-rate-myr",
        "--sediment-subduction-fraction",
        "--initial-boundary-influence",
        "--river-erosion",
        "--river-implicit",
        "--alluvial-transport",
        "--river-outlet-sea-level",
        "--sediment-travel-km",
        "--rift-breakup-thickness-m",
        "--continental-fraction",
        "--continental-root-kg-m3",
        "--continental-min-thickness-m",
        "--continental-max-thickness-m",
        "--orogenic-collapse-km2-myr",
        "--arc-continent-fraction",
        "--subduction-erosion-ratio",
        "--continental-growth-control",
        "--continental-growth-per-gyr",
        "--delta-time-myr",
        "--uplift-rate-window-myr",
    };
    std::string line = "   ";
    for (const char* flag : flags)
    {
        if (line.size() + std::string(flag).size() + 1 > 100)
        {
            out << line << '\n';
            line = "   ";
        }
        line += ' ';
        line += flag;
    }
    out << line << '\n';
}
