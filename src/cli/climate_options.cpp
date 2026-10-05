#include "cli/climate_options.hpp"

#include <stdexcept>

#ifndef ARRAYWIDTH
#define ARRAYWIDTH 2048
#endif

using namespace std;

bool isclimateoption(const string& argument)
{
    return argument == "--climate-resolution" || argument == "--climate-hydrology-resolution" ||
        argument == "--climate-coupling-iterations" || argument == "--climate-ocean-step-hours" ||
        argument == "--climate-ocean-heat-kernel" || argument == "--climate-ocean-heat-workers" ||
        argument == "--climate-ocean-coastal-scheme" || argument == "--climate-ocean-capture" ||
        argument == "--climate-atmosphere-solver" || argument == "--climate-atmosphere-capture" ||
        argument == "--climate-ocean-circulation-solver" || argument == "--climate-ocean-circulation-cache" ||
        argument == "--climate-mode" || argument == "--climate-tropical-closure" || argument == "--climate-transport" ||
        argument == "--climate-temperature-mode" || argument == "--climate-generic-levels" ||
        argument == "--climate-generic-reversal-hpa";
}

ClimateOptionResult parseclimateoption(const string& argument, const string& value, ClimateRunOptions& options, string* reason)
{
    if (!isclimateoption(argument))
        return ClimateOptionResult::notclimate;

    try
    {
        if (argument == "--climate-generic-levels")
        {
            // Comma-separated pressures strictly between 500 and 1000 hPa.
            options.genericLevelPressuresHpa.clear();
            size_t start = 0;
            while (start <= value.size())
            {
                const size_t end = value.find(',', start);
                const float pressure = stof(value.substr(start, end == string::npos ? string::npos : end - start));
                if (!(pressure > 500.0f && pressure < 1000.0f))
                    throw invalid_argument("climate-generic-levels takes pressures between 500 and 1000 hPa");
                options.genericLevelPressuresHpa.push_back(pressure);
                if (end == string::npos)
                    break;
                start = end + 1;
            }
            if (options.genericLevelPressuresHpa.size() > 8)
                throw invalid_argument("climate-generic-levels takes at most eight pressures");
        }
        else if (argument == "--climate-generic-reversal-hpa")
        {
            options.genericLevelReversalHpa = stof(value);
            if (!(options.genericLevelReversalHpa > 500.0f && options.genericLevelReversalHpa < 1000.0f))
                throw invalid_argument("climate-generic-reversal-hpa must be between 500 and 1000");
        }
        else if (argument == "--climate-resolution")
        {
            options.climateResolution = stoi(value);
            if (options.climateResolution < 8 || options.climateResolution % 4 != 0 || options.climateResolution > ARRAYWIDTH / 2)
                throw invalid_argument("climate-resolution must be a multiple of 4 between 8 and half ARRAYWIDTH");
        }
        else if (argument == "--climate-hydrology-resolution")
        {
            options.hydrologyResolution = stoi(value);
            if (options.hydrologyResolution < 8 || options.hydrologyResolution % 4 != 0 || options.hydrologyResolution > ARRAYWIDTH)
                throw invalid_argument("climate-hydrology-resolution must be a multiple of 4 between 8 and ARRAYWIDTH");
        }
        else if (argument == "--climate-coupling-iterations")
        {
            options.couplingIterations = stoi(value);
            if (options.couplingIterations < 1 || options.couplingIterations > 100)
                throw invalid_argument("climate-coupling-iterations must be between 1 and 100");
        }
        else if (argument == "--climate-ocean-step-hours")
        {
            options.oceanHeatStepHours = stoi(value);
            if (options.oceanHeatStepHours < 1 || options.oceanHeatStepHours > 24 || 24 % options.oceanHeatStepHours != 0)
                throw invalid_argument("climate-ocean-step-hours must be a positive divisor of 24");
        }
        else if (argument == "--climate-ocean-heat-kernel")
        {
            if (value != "gather" && value != "legacy") throw invalid_argument("climate-ocean-heat-kernel");
            options.gatherOceanHeatFluxes = value == "gather";
        }
        else if (argument == "--climate-ocean-heat-workers")
        {
            options.oceanHeatWorkers = stoi(value);
            if (options.oceanHeatWorkers < 0 || options.oceanHeatWorkers > 32)
                throw invalid_argument("climate-ocean-heat-workers must be between 0 and 32");
        }
        else if (argument == "--climate-ocean-coastal-scheme")
        {
            using Scheme = climateocean::CoastalScheme;
            if (value != "legacy" && value != "face-drag" && value != "face-forcing")
                throw invalid_argument("climate-ocean-coastal-scheme");
            options.oceanCoastalScheme = value == "legacy" ? Scheme::Legacy :
                value == "face-drag" ? Scheme::FaceDrag : Scheme::FaceForcing;
        }
        else if (argument == "--climate-ocean-capture")
            options.oceanCaptureDirectory = value;
        else if (argument == "--climate-atmosphere-solver")
        {
            using Solver = climateatmosphere::StationarySolver;
            if (value != "zonal" && value != "jacobi" && value != "legacy") throw invalid_argument("climate-atmosphere-solver");
            options.atmosphereSolver = value == "zonal" ? Solver::Zonal : value == "jacobi" ? Solver::ProjectedJacobi : Solver::LegacyJacobi;
        }
        else if (argument == "--climate-atmosphere-capture")
            options.atmosphereCaptureDirectory = value;
        else if (argument == "--climate-ocean-circulation-solver")
        {
            if (value != "krylov" && value != "relaxation") throw invalid_argument("climate-ocean-circulation-solver");
            options.krylovOceanCirculation = value == "krylov";
        }
        else if (argument == "--climate-ocean-circulation-cache")
        {
            if (value != "reuse" && value != "recompute") throw invalid_argument("climate-ocean-circulation-cache");
            options.reuseOceanCirculation = value == "reuse";
        }
        else if (argument == "--climate-mode")
        {
            if (value != "baseline" && value != "transient") throw invalid_argument("climate-mode");
            options.transientWeather = value == "transient";
        }
        else if (argument == "--climate-tropical-closure")
        {
            if (value != "mixed-layer" && value != "legacy") throw invalid_argument("climate-tropical-closure");
            options.mixedLayerTropics = value == "mixed-layer";
        }
        else if (argument == "--climate-transport")
        {
            if (value != "ffsl" && value != "mpdata") throw invalid_argument("climate-transport");
            options.transport = value == "ffsl" ? climatecontext::TransportMethod::Ffsl : climatecontext::TransportMethod::Mpdata;
        }
        else if (argument == "--climate-temperature-mode")
        {
            if (value != "calibrated" && value != "radiative") throw invalid_argument("climate-temperature-mode");
            options.calibrateTemperature = value == "calibrated";
        }
    }
    catch (const exception& error)
    {
        if (reason != nullptr)
            *reason = error.what();
        return ClimateOptionResult::invalid;
    }
    return ClimateOptionResult::parsed;
}

void applyclimateoptions(const ClimateRunOptions& options, climatecontext::Context& context)
{
    context.resolution = { options.climateResolution, options.climateResolution,
        options.hydrologyResolution ? options.hydrologyResolution : 2 * options.climateResolution, options.climateResolution / 2 };
    context.calibrateMeanTemperature = options.calibrateTemperature;
    context.transientWeather = options.transientWeather;
    context.mixedLayerTropics = options.mixedLayerTropics;
    context.transportMethod = options.transport;
    context.oceanHeatStepHours = options.oceanHeatStepHours;
    context.reuseOceanCirculation = options.reuseOceanCirculation;
    context.krylovOceanCirculation = options.krylovOceanCirculation;
    context.gatherOceanHeatFluxes = options.gatherOceanHeatFluxes;
    context.oceanHeatWorkers = options.oceanHeatWorkers;
    context.oceanCoastalScheme = options.oceanCoastalScheme;
    context.oceanCaptureDirectory = options.oceanCaptureDirectory;
    context.atmosphereSolver = options.atmosphereSolver;
    context.atmosphereCaptureDirectory = options.atmosphereCaptureDirectory;
    context.couplingIterationLimit = options.couplingIterations;
    context.genericLevelPressuresHpa = options.genericLevelPressuresHpa;
    context.genericLevelReversalHpa = options.genericLevelReversalHpa;
}

void printclimateoptionsusage(ostream& out)
{
    out << "  --climate-resolution <width>  Atmosphere/ocean width; moisture uses twice, weather half (default 128).\n";
    out << "  --climate-coupling-iterations <1..100>  Coupled annual update limit (default 24).\n";
    out << "  --climate-mode <baseline|transient>  Seasonal baseline or sampled weather (default baseline).\n";
    out << "  --climate-tropical-closure <mixed-layer|legacy>  Tropical circulation (default legacy; mixed-layer experimental).\n";
    out << "  --climate-ocean-step-hours <hours>  Maximum ocean heat step; positive divisor of 24 (default 24).\n";
    out << "  --climate-hydrology-resolution <width>  Native moisture grid width (default 2x climate resolution).\n";
    out << "  --climate-ocean-circulation-cache <reuse|recompute>  Seasonal spin-up circulation (default reuse).\n";
    out << "  --climate-ocean-circulation-solver <krylov|relaxation>  Ocean momentum solver (default krylov).\n";
    out << "  --climate-ocean-heat-kernel <gather|legacy>  Heat transport execution (default gather).\n";
    out << "  --climate-ocean-heat-workers <0..32>  Maximum heat workers; 0 automatic, 1 serial.\n";
    out << "  --climate-ocean-coastal-scheme <legacy|face-drag|face-forcing>  Depth/forcing discretization (default face-forcing).\n";
    out << "  --climate-ocean-capture <directory>  Optional exact seasonal ocean forcing fixtures.\n";
    out << "  --climate-atmosphere-solver <zonal|jacobi|legacy>  Pressure preconditioner (default zonal).\n";
    out << "  --climate-atmosphere-capture <directory>  Optional exact stationary-solver input fixtures.\n";
    out << "  --climate-transport <ffsl|mpdata>  Conservative moisture transport (default ffsl).\n";
    out << "  --climate-temperature-mode <calibrated|radiative>  Mean-temperature calibration (default calibrated).\n";
    out << "  --climate-generic-levels <hPa[,hPa...]>  Experimental diagnostic levels from one rule set (500-1000 hPa).\n";
    out << "  --climate-generic-reversal-hpa <hPa>  Heat-low depth of those levels' projection (default 760).\n";
}
