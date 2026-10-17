#include "gtest/gtest.h"
#include "simulations/geology/plate_tectonics/lithosphere.hpp"
#include "cli/tectonic_options.hpp"

#include <stdexcept>
#include <string>
#include <vector>

namespace {

// The engine as world_pipeline starts it: the baseline recipe, then the scenario
// options of the command line (--key value), on a small grid.
void start(const std::vector<std::string>& options) {
    tectonics::scenario::Scenario scenario = tectonics::scenario::baseline();
    scenario.width = 32;
    scenario.height = 16;
    for (size_t i = 0; i + 1 < options.size(); i += 2) {
        if (!parsetectonicoption(scenario, options[i], options[i + 1])) {
            throw std::runtime_error("not a scenario option: " + options[i]);
        }
    }
    lithosphere engine(scenario);
}

} // namespace

TEST(TectonicOptions, BaselineAndExplicitValuesStart) {
    EXPECT_NO_THROW(start({}));
    EXPECT_NO_THROW(start({"--delta-time-myr", ".5", "--evolution-duration-myr", "2", "--transport-subdivision", "2"}));
    EXPECT_NO_THROW(start({"--continental-max-thickness-m", "40000", "--orogenic-collapse-km2-myr", "3000"}));
    EXPECT_NO_THROW(start({"--continental-root-kg-m3", "246"}));
    EXPECT_NO_THROW(start({"--continental-root-kg-m3", "246", "--initial-margin-width-m", "2800"}));
    EXPECT_NO_THROW(start({"--continental-root-kg-m3", "246", "--continental-root-max-thickness-m", "50000"}));
    EXPECT_NO_THROW(start({"--continental-root-max-thickness-m", "0"}));
    EXPECT_NO_THROW(start({"--continental-root-kg-m3", "246", "--continental-root-max-thickness-m", "50000",
                           "--continental-root-delamination-width-m", "10000"}));
    EXPECT_NO_THROW(start({"--continental-root-delamination-width-m", "0"}));
    EXPECT_NO_THROW(start({"--continental-max-thickness-m", "34000", "--orogenic-collapse-km2-myr", "30000",
                           "--orogenic-collapse-onset-thickness-m", "60000"}));
    EXPECT_NO_THROW(start({"--continental-max-thickness-m", "34000", "--orogenic-collapse-km2-myr", "30000",
                           "--orogenic-collapse-onset-thickness-m", "34000"}));
    EXPECT_NO_THROW(start({"--orogenic-collapse-onset-thickness-m", "0"}));
    EXPECT_NO_THROW(start({"--continental-max-thickness-m", "34000", "--orogenic-collapse-km2-myr", "30000",
                           "--orogenic-collapse-onset-thickness-m", "65000", "--orogenic-collapse-margin-cap", "on",
                           "--orogenic-collapse-below-onset-km2-myr", "10000"}));
    EXPECT_NO_THROW(start({"--orogenic-collapse-margin-cap", "off", "--orogenic-collapse-below-onset-km2-myr", "0"}));
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
        EXPECT_ANY_THROW(start(options)) << options[0] << ' ' << options[1];
    }
}
