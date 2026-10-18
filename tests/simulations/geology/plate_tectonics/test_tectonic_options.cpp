#include "gtest/gtest.h"
#include "simulations/geology/plate_tectonics/lithosphere.hpp"
#include "cli/tectonic_options.hpp"

#include <stdexcept>
#include <string>
#include <vector>

namespace {

// The baseline recipe without its Earth-like set: those members at the engine's
// struct defaults, so that each case below sets only the inputs of the rule it probes.
tectonics::scenario::Scenario without_earth_set() {
    tectonics::scenario::Scenario scenario = tectonics::scenario::baseline();
    const tectonics::scenario::Scenario defaults{};
    scenario.continental_fraction = defaults.continental_fraction;
    scenario.continental_root_kg_m3 = defaults.continental_root_kg_m3;
    scenario.continental_root_max_thickness_m = defaults.continental_root_max_thickness_m;
    scenario.continental_root_delamination_width_m = defaults.continental_root_delamination_width_m;
    auto& margins = scenario.margin_processes;
    const auto& legacy = defaults.margin_processes;
    margins.rift_breakup_thickness_m = legacy.rift_breakup_thickness_m;
    margins.continental_min_thickness_m = legacy.continental_min_thickness_m;
    margins.sediment_subduction_fraction = legacy.sediment_subduction_fraction;
    margins.subduction_erosion_ratio = legacy.subduction_erosion_ratio;
    margins.continental_growth_control = legacy.continental_growth_control;
    margins.continental_growth_per_gyr = legacy.continental_growth_per_gyr;
    margins.continental_growth_continent_only = legacy.continental_growth_continent_only;
    margins.continental_max_thickness_m = legacy.continental_max_thickness_m;
    margins.orogenic_collapse_km2_myr = legacy.orogenic_collapse_km2_myr;
    margins.orogenic_collapse_onset_thickness_m = legacy.orogenic_collapse_onset_thickness_m;
    margins.orogenic_collapse_below_onset_km2_myr = legacy.orogenic_collapse_below_onset_km2_myr;
    return scenario;
}

// The engine as world_pipeline starts it: the baseline recipe (or the given start),
// then the scenario options of the command line (--key value), on a small grid.
void start(const std::vector<std::string>& options,
           tectonics::scenario::Scenario scenario = tectonics::scenario::baseline()) {
    scenario.width = 32;
    scenario.height = 16;
    for (size_t i = 0; i + 1 < options.size(); i += 2) {
        if (!parsetectonicoption(scenario, options[i], options[i + 1])) {
            throw std::runtime_error("not a scenario option: " + options[i]);
        }
    }
    lithosphere engine(scenario);
}

// A case of a single rule: the options on the baseline without its Earth-like set.
void probe(const std::vector<std::string>& options) {
    start(options, without_earth_set());
}

} // namespace

TEST(TectonicOptions, BaselineAndExplicitValuesStart) {
    EXPECT_NO_THROW(start({}));
    EXPECT_NO_THROW(probe({"--delta-time-myr", ".5", "--evolution-duration-myr", "2", "--transport-subdivision", "2"}));
    EXPECT_NO_THROW(probe({"--continental-max-thickness-m", "40000", "--orogenic-collapse-km2-myr", "3000"}));
    EXPECT_NO_THROW(probe({"--continental-root-kg-m3", "246"}));
    EXPECT_NO_THROW(probe({"--continental-root-kg-m3", "246", "--initial-margin-width-m", "2800"}));
    EXPECT_NO_THROW(probe({"--continental-root-kg-m3", "246", "--continental-root-max-thickness-m", "50000"}));
    EXPECT_NO_THROW(probe({"--continental-root-max-thickness-m", "0"}));
    EXPECT_NO_THROW(probe({"--continental-root-kg-m3", "246", "--continental-root-max-thickness-m", "50000",
                           "--continental-root-delamination-width-m", "10000"}));
    EXPECT_NO_THROW(probe({"--continental-root-delamination-width-m", "0"}));
    EXPECT_NO_THROW(probe({"--continental-max-thickness-m", "34000", "--orogenic-collapse-km2-myr", "30000",
                           "--orogenic-collapse-onset-thickness-m", "60000"}));
    EXPECT_NO_THROW(probe({"--continental-max-thickness-m", "34000", "--orogenic-collapse-km2-myr", "30000",
                           "--orogenic-collapse-onset-thickness-m", "34000"}));
    EXPECT_NO_THROW(probe({"--orogenic-collapse-onset-thickness-m", "0"}));
    EXPECT_NO_THROW(probe({"--continental-max-thickness-m", "34000", "--orogenic-collapse-km2-myr", "30000",
                           "--orogenic-collapse-onset-thickness-m", "65000", "--orogenic-collapse-margin-cap", "on",
                           "--orogenic-collapse-below-onset-km2-myr", "10000"}));
    EXPECT_NO_THROW(probe({"--orogenic-collapse-margin-cap", "off", "--orogenic-collapse-below-onset-km2-myr", "0"}));
}

// A value the option parser or the engine refuses stops the run before it starts.
TEST(TectonicOptions, InvalidScenarioValuesAreRejected) {
    const std::vector<std::vector<std::string>> cases = {
        {"--delta-time-myr", "0"},
        {"--evolution-duration-myr", "0"},
        {"--evolution-duration-myr", "nan"},
        {"--evolution-duration-myr", "-1"},
        {"--evolution-duration-myr", "1e20"},
        {"--world-width-km", "0"},
        {"--mantle-flow-speed-cm-yr", "-1"},
        {"--coastal-gpe-strength", "1.01"},
        {"--boundary-force-ramp-cm-yr", "-1"},
        {"--initial-margin-width-m", "4001"},
        {"--sediment-load-feedback", "maybe"},
        {"--noise-pass-wavelength-km", "0"},
        {"--transport-subdivision", "3"},
        {"--deformation-variation", ".6"},
        {"--subduction-coupling", "1.1"},
        {"--submarine-erosion-ratio", "1.1"},
        {"--thermal-subsidence-strength", "-1"},
        {"--initial-ocean-depth-m", "65535"},
        {"--plume-count", "65"},
        {"--plume-tail-radius-km", "0"},
        {"--plume-weakening", "1.5"},
        {"--plume-count", "2", "--plume-weakening", ".5"}, // weakening needs progressive rifting
        {"--continental-max-thickness-m", "12000"}, // not above the rift onset
        {"--continental-max-thickness-m", "-1"},
        {"--orogenic-collapse-km2-myr", "3000"}, // needs a maximum thickness
        {"--continental-fraction", ".40", "--initial-ocean-depth-m", "3500"}, // slope foot at 3500 m
        {"--continental-root-kg-m3", "500"}, // no crust buoyancy left
        {"--continental-root-kg-m3", "-1"},
        // At root 246 a continental column has no crust below 3916 m: the slope foot
        // 2500 m + width / 2 needs a width under 2832 m.
        {"--continental-root-kg-m3", "246", "--initial-margin-width-m", "3000"},
        {"--continental-growth-continent-only", "maybe"},
        {"--continental-root-max-thickness-m", "30000"}, // the cap must lie above 35 km
        {"--continental-root-max-thickness-m", "35000"},
        {"--continental-root-max-thickness-m", "-1"},
        {"--continental-root-delamination-width-m", "10000"}, // a ramp needs a cap
        {"--continental-root-max-thickness-m", "50000", "--continental-root-delamination-width-m", "-1"},
        {"--continental-root-max-thickness-m", "50000", "--continental-root-delamination-width-m", "inf"},
        {"--orogenic-collapse-onset-thickness-m", "60000"}, // an onset needs a maximum thickness
        {"--continental-max-thickness-m", "40000", "--orogenic-collapse-onset-thickness-m", "30000"}, // below it
        {"--continental-max-thickness-m", "40000", "--orogenic-collapse-onset-thickness-m", "-1"},
        {"--continental-max-thickness-m", "40000", "--orogenic-collapse-onset-thickness-m", "nan"},
        {"--orogenic-collapse-margin-cap", "maybe"},
        {"--continental-max-thickness-m", "40000", "--orogenic-collapse-km2-myr", "3000",
         "--orogenic-collapse-onset-thickness-m", "65000", "--orogenic-collapse-below-onset-km2-myr", "3001"}, // creep > D
        {"--continental-max-thickness-m", "40000", "--orogenic-collapse-onset-thickness-m", "65000"}, // no collapse
        {"--continental-max-thickness-m", "40000", "--orogenic-collapse-onset-thickness-m", "65000",
         "--orogenic-collapse-margin-cap", "on"},
        {"--continental-max-thickness-m", "40000", "--orogenic-collapse-margin-cap", "on"}, // needs an onset
        {"--continental-max-thickness-m", "40000", "--orogenic-collapse-below-onset-km2-myr", "10000"},
        {"--continental-max-thickness-m", "40000", "--orogenic-collapse-onset-thickness-m", "65000",
         "--orogenic-collapse-below-onset-km2-myr", "-1"},
    };
    for (const auto& options : cases) {
        EXPECT_ANY_THROW(probe(options)) << options[0] << ' ' << options[1];
    }
}
