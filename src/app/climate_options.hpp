#pragma once

// Climate solver options shared by world_pipeline and climate_pipeline.

#include <ostream>
#include <string>
#include <vector>

#include "climate_context.hpp"

struct ClimateRunOptions
{
    int climateResolution = 128;
    int hydrologyResolution = 0; // Zero keeps the 2x climate default.
    int couplingIterations = 0;  // Zero uses the tuning default.
    int oceanHeatStepHours = 24;
    bool reuseOceanCirculation = true;
    bool krylovOceanCirculation = true;
    bool gatherOceanHeatFluxes = true;
    int oceanHeatWorkers = 0;
    climateocean::CoastalScheme oceanCoastalScheme = climateocean::CoastalScheme::FaceForcing;
    std::string oceanCaptureDirectory;
    climateatmosphere::StationarySolver atmosphereSolver = climateatmosphere::StationarySolver::Zonal;
    std::string atmosphereCaptureDirectory;
    bool calibrateTemperature = true;
    bool transientWeather = false;
    bool mixedLayerTropics = false;
    climatecontext::TransportMethod transport = climatecontext::TransportMethod::Ffsl;
    std::vector<float> genericLevelPressuresHpa; // experimental diagnostic levels
    float genericLevelReversalHpa = 760.0f;
};

enum class ClimateOptionResult { notclimate, parsed, invalid };

// True for the --climate-* options, which all take a value.
bool isclimateoption(const std::string& argument);
// Parses one option and its value; invalid values report a reason.
ClimateOptionResult parseclimateoption(const std::string& argument, const std::string& value, ClimateRunOptions& options,
    std::string* reason = nullptr);
void applyclimateoptions(const ClimateRunOptions& options, climatecontext::Context& context);
void printclimateoptionsusage(std::ostream& out);
