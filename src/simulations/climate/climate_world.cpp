#include "simulations/climate/climate_world.hpp"

#include <cmath>
#include <utility>

#include "core/climate_grid.hpp"
#include "simulations/climate/climate_hydrology.hpp"
#include "simulations/climate/climate_weather.hpp"
#include "simulations/climate/climate_parameters.hpp"

void climateworld::resize(int width, int height)
{
    inputs_.width = width;
    inputs_.height = height;
    levels_.clear();
    diagnostics_.reset();
    // A new grid discards the solver state and keeps the run's options (cli/climate_options.cpp):
    // a procedural world is resized after its options were applied.
    climatecontext::Context options;
    options.resolution = climatecontext_.resolution;
    options.calibrateMeanTemperature = climatecontext_.calibrateMeanTemperature;
    options.transientWeather = climatecontext_.transientWeather;
    options.mixedLayerTropics = climatecontext_.mixedLayerTropics;
    options.genericLevelPressuresHpa = climatecontext_.genericLevelPressuresHpa;
    options.genericLevelReversalHpa = climatecontext_.genericLevelReversalHpa;
    options.transportMethod = climatecontext_.transportMethod;
    options.oceanHeatStepHours = climatecontext_.oceanHeatStepHours;
    options.reuseOceanCirculation = climatecontext_.reuseOceanCirculation;
    options.krylovOceanCirculation = climatecontext_.krylovOceanCirculation;
    options.gatherOceanHeatFluxes = climatecontext_.gatherOceanHeatFluxes;
    options.oceanHeatWorkers = climatecontext_.oceanHeatWorkers;
    options.oceanCoastalScheme = climatecontext_.oceanCoastalScheme;
    options.oceanCaptureDirectory = climatecontext_.oceanCaptureDirectory;
    options.atmosphereSolver = climatecontext_.atmosphereSolver;
    options.atmosphereCaptureDirectory = climatecontext_.atmosphereCaptureDirectory;
    options.couplingIterationLimit = climatecontext_.couplingIterationLimit;
    climatecontext_ = std::move(options);
    const int cellcount = (inputs_.width + 1) * (inputs_.height + 1);
    jantempmap.resize(inputs_.width + 1, inputs_.height + 1);
    jultempmap.resize(inputs_.width + 1, inputs_.height + 1);
    climatemap.resize(inputs_.width + 1, inputs_.height + 1);
    biomemap.resize(inputs_.width + 1, inputs_.height + 1);
    janrainmap.resize(inputs_.width + 1, inputs_.height + 1);
    julrainmap.resize(inputs_.width + 1, inputs_.height + 1);
    janmountainrainmap.resize(inputs_.width + 1, inputs_.height + 1);
    julmountainrainmap.resize(inputs_.width + 1, inputs_.height + 1);
    janmountainraindirmap.resize(inputs_.width + 1, inputs_.height + 1);
    julmountainraindirmap.resize(inputs_.width + 1, inputs_.height + 1);
    seaicemap.resize(inputs_.width + 1, inputs_.height + 1);
    windmap.resize(inputs_.width + 1, inputs_.height + 1);
    // The save format stores subtropical-high rows for the full 2048-column capacity.
    horselats.resize(ARRAYWIDTH, 6);

    auto resizefield = [cellcount](std::array<std::vector<short>, CLIMATESEASONCOUNT>& field)
    {
        for (std::vector<short>& season : field)
            season.assign(cellcount, 0);
    };
    auto resizefloatfield = [cellcount](
        std::array<std::vector<float>, CLIMATESEASONCOUNT>& field)
    {
        for (std::vector<float>& season : field)
            season.assign(cellcount, 0.0f);
    };

    resizefloatfield(seasonaltempmaps);
    resizefloatfield(seasonalrainmaps);
    resizefloatfield(seasonalpressuremaps);
    resizefloatfield(seasonaluwindmaps);
    resizefloatfield(seasonalvwindmaps);
    resizefloatfield(seasonalupperheightmaps);
    resizefloatfield(seasonalupperuwindmaps);
    resizefloatfield(seasonaluppervwindmaps);
    resizefloatfield(seasonalverticalvelocitymaps);
    resizefloatfield(seasonalcurrentumaps);
    resizefloatfield(seasonalcurrentvmaps);
    resizefloatfield(seasonalsstmaps);
    resizefloatfield(seasonalevaporationmaps);
    resizefloatfield(seasonalmaritimeinfluencemaps);
    resizefloatfield(seasonalmaritimethermalanomalymaps);
    resizefield(seasonalmaritimefetchmaps);
    resizefloatfield(seasonalmoisturemaps);
    resizefloatfield(seasonalconvergencemaps);
    resizefloatfield(seasonalupliftmaps);
    resizefloatfield(seasonalsubsidencemaps);
}

std::vector<short> climateworld::legacyfield(const std::vector<float>& values, float scale)
{
    std::vector<short> encoded(values.size());
    for (std::size_t index = 0; index < values.size(); index++)
        encoded[index] = static_cast<short>(legacyvalue(values[index], scale));
    return encoded;
}

std::vector<short> climateworld::legacyfield(const worldgrid<float>& values)
{
    std::vector<short> encoded(values.size());
    for (std::size_t index = 0; index < values.size(); index++)
        encoded[index] = static_cast<short>(legacyvalue(values.data()[index]));
    return encoded;
}

const std::vector<std::uint8_t>& climateworld::weatheranomalystate() const
{
    UW_TABLE_READ(weather_state);
    return weatheranomalystatebytes;
}

void climateworld::setweatheranomalystate(std::vector<std::uint8_t> state)
{
    weatheranomalystatebytes = std::move(state);
}

void climateworld::syncseasonalclimatefromlegacy()
{
    bool hasexplicittransitiontemp = false;
    bool hasexplicittransitionrain = false;

    for (const float value : seasonaltempmaps[seasonapril])
    {
        if (value != 0.0f)
        {
            hasexplicittransitiontemp = true;
            break;
        }
    }

    if (hasexplicittransitiontemp == false)
    {
        for (const float value : seasonaltempmaps[seasonoctober])
        {
            if (value != 0.0f)
            {
                hasexplicittransitiontemp = true;
                break;
            }
        }
    }

    for (const float value : seasonalrainmaps[seasonapril])
    {
        if (value != 0.0f)
        {
            hasexplicittransitionrain = true;
            break;
        }
    }

    if (hasexplicittransitionrain == false)
    {
        for (const float value : seasonalrainmaps[seasonoctober])
        {
            if (value != 0.0f)
            {
                hasexplicittransitionrain = true;
                break;
            }
        }
    }

    for (int y = 0; y <= inputs_.height; y++)
    {
        for (int x = 0; x <= inputs_.width; x++)
        {
            const int index = seasonalclimateindex(x, y);

            const float januarytemp = jantempmap[x][y];
            const float julytemp = jultempmap[x][y];
            const float januaryrain = janrainmap[x][y];
            const float julyrain = julrainmap[x][y];

            seasonaltempmaps[seasonjanuary][index] = jantempmap[x][y];
            seasonaltempmaps[seasonjuly][index] = jultempmap[x][y];
            seasonalrainmaps[seasonjanuary][index] = janrainmap[x][y];
            seasonalrainmaps[seasonjuly][index] = julrainmap[x][y];

            float summertemp = julytemp;
            float wintertemp = januarytemp;

            if (inputs_.perihelion == 1)
            {
                summertemp = januarytemp;
                wintertemp = julytemp;
            }

            const float winterstrength = 0.5f + inputs_.eccentricity * 0.5f;
            const float summerstrength = 1.0f - winterstrength;
            float transitiontemp = summertemp * summerstrength + wintertemp * winterstrength;

            float fourseason = inputs_.tilt * 0.294592f - 2.45428f;
            float lat = static_cast<float>(y);

            if (y > inputs_.height / 2.0f)
                lat = static_cast<float>(inputs_.height - y);

            const float fourseasonstrength = lat / (static_cast<float>(inputs_.height) / 2.0f);
            const float transitiontempdiff = (fourseason * fourseasonstrength) / 2.0f;

            if (hasexplicittransitiontemp == false)
            {
                seasonaltempmaps[seasonapril][index] = transitiontemp + transitiontempdiff;
                seasonaltempmaps[seasonoctober][index] = transitiontemp + transitiontempdiff;
            }

            if (hasexplicittransitionrain == false)
            {
                float apriljanfactor = 0.5f;
                float apriljulfactor = 0.5f;

                if (julytemp > januarytemp && julyrain > januaryrain && julyrain > 0.0f)
                {
                    const float monsoonfactor = 1.0f - januaryrain / julyrain;
                    apriljanfactor = monsoonfactor * 0.9f;
                    apriljulfactor = 1.0f - apriljanfactor;
                }

                if (julytemp < januarytemp && julyrain < januaryrain && januaryrain > 0.0f)
                {
                    const float monsoonfactor = 1.0f - julyrain / januaryrain;
                    apriljanfactor = monsoonfactor * 0.7f;
                    apriljulfactor = 1.0f - apriljanfactor;
                }

                seasonalrainmaps[seasonapril][index] = januaryrain * apriljanfactor + julyrain * apriljulfactor;

                float octoberjanfactor = 0.5f;
                float octoberjulfactor = 0.5f;

                if (julytemp > januarytemp && julyrain > januaryrain && julyrain > 0.0f)
                {
                    const float monsoonfactor = 1.0f - januaryrain / julyrain;
                    octoberjulfactor = monsoonfactor * 0.7f;
                    octoberjanfactor = 1.0f - octoberjulfactor;
                }

                if (julytemp < januarytemp && julyrain < januaryrain && januaryrain > 0.0f)
                {
                    const float monsoonfactor = 1.0f - julyrain / januaryrain;
                    octoberjulfactor = monsoonfactor * 0.9f;
                    octoberjanfactor = 1.0f - octoberjulfactor;
                }

                seasonalrainmaps[seasonoctober][index] = januaryrain * octoberjanfactor + julyrain * octoberjulfactor;
            }
        }
    }
}

bool climateworld::advanceweather(float elapsedSeconds, int horizontalCells)
{
    climateweather::ShallowWaterState state;
    if (!weatheranomalystatebytes.empty() && !climateweather::deserializeState(weatheranomalystatebytes, state))
        return false;
    const int columns = horizontalCells > 0 ? std::clamp(horizontalCells, 8, 128) :
        state.columns > 0 ? state.columns : tuning::climate::moistureadvection::weatherAnomalyHorizontalCells;
    if (state.layers.empty()) state = climateweather::makeState(columns, columns / 2, 1, static_cast<std::uint64_t>(seed()));
    else if (columns != state.columns) state = climateweather::resampleState(state, columns);
    climateweather::ShallowWaterConfig config;
    config.layerCount = static_cast<int>(state.layers.size());
    config.rotationDirection = rotation() ? 1.0f : -1.0f;
    config.gravityMetresPerSecondSquared *= std::max(0.05f, gravity());
    config.lowerMeanDepthMetres = tuning::climate::moistureadvection::weatherAnomalyEquivalentDepthMetres;
    config.lowerDragTimeSeconds = tuning::climate::moistureadvection::weatherAnomalyDragTimeDays * 86400.0f;
    config.heightRelaxationTimeSeconds = tuning::climate::moistureadvection::weatherAnomalyRelaxationTimeDays * 86400.0f;
    config.stochasticHeightForcingMetresPerSecond = tuning::climate::moistureadvection::weatherAnomalyStochasticHeightMps;
    if (!std::isfinite(elapsedSeconds) || elapsedSeconds < 0.0f) return false;
    float remaining = elapsedSeconds;
    // Refresh the seasonally interpolated jets at least once per simulated day.
    while (remaining > 0.0f)
    {
        const float step = std::min(remaining, 86400.0f);
        const float season = static_cast<float>(std::fmod((state.elapsedSeconds + step * 0.5) / (365.25 * 86400.0), 1.0) * 4.0);
        const int first = static_cast<int>(season), second = (first + 1) % CLIMATESEASONCOUNT;
        const float fraction = season - first;
        const auto background = [&](bool south, bool upper)
        {
            std::vector<float> raster((width() + 1) * (height() + 1));
            for (int y = 0; y <= height(); ++y)
                for (int x = 0; x <= width(); ++x)
                {
                    const auto wind = [&](int s) { return upper ? (south ? seasonaluppervwindfloat(s, x, y) : seasonalupperuwindfloat(s, x, y)) :
                        (south ? seasonalvwindfloat(s, x, y) : seasonaluwindfloat(s, x, y)); };
                    raster[y * (width() + 1) + x] = wind(first) * (1.0f - fraction) + wind(second) * fraction;
                }
            return climategrid::remapField(
                width() + 1, height() + 1, raster, state.columns, state.rows);
        };
        climateweather::ShallowWaterForcing forcing;
        for (int layer = 0; layer < config.layerCount; ++layer)
        {
            forcing.backgroundEastWindMps.push_back(background(false, layer > 0));
            forcing.backgroundSouthWindMps.push_back(background(true, layer > 0));
        }
        const auto diagnostics = climateweather::advance(state, config, forcing, step);
        if (!diagnostics.finite || !diagnostics.bounded) return false;
        remaining -= step;
    }
    weatheranomalystatebytes = climateweather::serializeState(state);
    return true;
}

void climateworld::recordprocessdiagnostics(int season, const climatehydrology::SeasonalProcessFields& fields)
{
    if (!validseasonindex(season))
        return;
    climatediagnostics& diagnostics = diagnostics_.ensure();
    if (!diagnostics.hasprocess || season == 0)
        diagnostics.process = {};
    diagnostics.hasprocess = true;
    diagnostics.process[season] = fields;
}

void climateworld::recordoceandiagnostics(int season, int columns, int rows, const climateocean::OceanState& state)
{
    if (!validseasonindex(season))
        return;
    climatediagnostics& diagnostics = diagnostics_.ensure();
    if (!diagnostics.hasocean || season == 0)
        diagnostics.ocean = {};
    diagnostics.hasocean = true;
    diagnostics.oceancolumns = columns;
    diagnostics.oceanrows = rows;
    diagnostics.ocean[season] = state;
}

void climateworld::recordweatherdiagnostics(int season, int columns, int rows,
    const climateweather::WeatherStatistics& statistics, int layer)
{
    if (!validseasonindex(season) || layer < 0 || layer > 1 || columns < 3 || rows < 2)
        return;
    climatediagnostics& diagnostics = diagnostics_.ensure();
    if (!diagnostics.hasweather || diagnostics.weathercolumns != columns || diagnostics.weatherrows != rows ||
        (season == 0 && layer == 0))
        diagnostics.weather = {};
    diagnostics.hasweather = true;
    diagnostics.weathercolumns = columns;
    diagnostics.weatherrows = rows;
    diagnostics.weather[layer][season] = statistics;
}

void climateworld::updatelevels()
{
    const climatecontext::WindFields& winds = climatecontext_.winds;
    const int columns = inputs_.width + 1;
    const int rows = inputs_.height + 1;
    const std::size_t cells = static_cast<std::size_t>(columns) * static_cast<std::size_t>(rows);
    const auto sized = [&](const climatecontext::FloatGrid& grid)
    {
        return static_cast<int>(grid.size()) == columns && !grid.empty() && static_cast<int>(grid[0].size()) == rows;
    };

    levels_.assign(3 + winds.genericlevels.size(), atmospherelevel{});
    for (std::size_t index = 0; index < winds.genericlevels.size(); ++index)
    {
        atmospherelevel& level = levels_[3 + index];
        const float pressure = winds.genericlevels[index].pressureHpa;
        level.name = "generic-" + std::to_string(static_cast<int>(std::lround(pressure))) + "hPa";
        level.pressurehpa = pressure;
    }
    levels_[0].name = "10m";
    levels_[1].name = "850hPa";
    levels_[1].pressurehpa = 850.0f;
    levels_[2].name = "500hPa";
    levels_[2].pressurehpa = 500.0f;

    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
    {
        for (atmospherelevel& level : levels_)
        {
            level.eastwind[season].assign(cells, 0.0f);
            level.northwind[season].assign(cells, 0.0f);
            level.available[season].assign(cells, 0);
        }
        levels_[2].geopotentialheight[season].assign(cells, 0.0f);
        levels_[2].verticalvelocity[season].assign(cells, 0.0f);

        if (!winds.populated[season] || !sized(winds.surfacewindu[season]) || !sized(winds.upperwindu[season]))
            continue;
        const bool lowlevel = winds.hasLowLevelProfile[season] && sized(winds.lowlevelwindu[season]) && sized(winds.lowlevelavailable[season]);
        const bool seasonalfields = seasonalupperheightmaps[season].size() == cells && seasonalverticalvelocitymaps[season].size() == cells;

        for (int y = 0; y < rows; y++)
        {
            for (int x = 0; x < columns; x++)
            {
                const std::size_t cell = static_cast<std::size_t>(y) * columns + x;
                // The solver works in an east/south basis; levels store northward wind.
                levels_[0].eastwind[season][cell] = winds.surfacewindu[season][x][y];
                levels_[0].northwind[season][cell] = -winds.surfacewindv[season][x][y];
                levels_[0].available[season][cell] = 1;
                if (lowlevel)
                {
                    levels_[1].eastwind[season][cell] = winds.lowlevelwindu[season][x][y];
                    levels_[1].northwind[season][cell] = -winds.lowlevelwindv[season][x][y];
                    levels_[1].available[season][cell] = winds.lowlevelavailable[season][x][y] > 0.0f ? 1 : 0;
                }
                levels_[2].eastwind[season][cell] = winds.upperwindu[season][x][y];
                levels_[2].northwind[season][cell] = -winds.upperwindv[season][x][y];
                levels_[2].available[season][cell] = 1;
                if (seasonalfields)
                {
                    levels_[2].geopotentialheight[season][cell] = seasonalupperheightmaps[season][cell];
                    levels_[2].verticalvelocity[season][cell] = seasonalverticalvelocitymaps[season][cell];
                }
                for (std::size_t index = 0; index < winds.genericlevels.size(); ++index)
                {
                    const auto& generic = winds.genericlevels[index];
                    if (!generic.solved[season] || !sized(generic.windu[season]))
                        continue;
                    levels_[3 + index].eastwind[season][cell] = generic.windu[season][x][y];
                    levels_[3 + index].northwind[season][cell] = -generic.windv[season][x][y];
                    levels_[3 + index].available[season][cell] = generic.available[season][x][y] > 0.5f ? 1 : 0;
                }
            }
        }
    }
}
