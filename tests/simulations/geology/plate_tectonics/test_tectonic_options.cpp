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
    };
    for (const auto& options : cases) {
        EXPECT_ANY_THROW(start(options)) << options[0] << ' ' << options[1];
    }
}
