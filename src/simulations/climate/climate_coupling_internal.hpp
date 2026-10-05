#pragma once

// Shared by the sources of the coupled climate processes (pressure, winds, ocean
// currents, sea surface temperature, rainfall and coupling): grid aliases and
// small helpers on the world grid. Not for use outside simulations/climate.

#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <vector>

#include "core/classes.hpp"
#include "core/climate_grid.hpp"
#include "core/grid_coordinates.hpp"
#include "io/climate_experiment_input.hpp"
#include "simulations/climate/climate_atmosphere.hpp"
#include "simulations/climate/climate_context.hpp"
#include "simulations/climate/climate_energy.hpp"
#include "simulations/climate/climate_world.hpp"
#include "world/generation_tuning.hpp"
#include "simulations/climate/climate_parameters.hpp"

using floatgrid = std::vector<std::vector<float>>;
using hydrologyfloatgrid = climatecontext::RowMajorGrid<float>;
using hydrologybytegrid = climatecontext::RowMajorGrid<signed char>;

inline const climateio::experiment::Inputs& referenceexperiment()
{
    static const climateio::experiment::Inputs inputs([] {
        const char* directory = std::getenv("UW_CLIMATE_EXPERIMENT_INPUT");
        if (directory && *directory)
            std::cout << "REFERENCE-FORCED DIAGNOSTIC: " << directory << '\n';
        return directory && *directory ? std::filesystem::path(directory) : std::filesystem::path{};
    }());
    return inputs;
}

inline float worldpermanentlandicefraction(const climateworld& world, int x, int y)
{
    if (!tuning::climate::energybalance::couplePermanentLandIce || world.sea(x,y)) return 0.0f;
    float warmest = world.seasonaltempfloat(0,x,y);
    for (int season = 1; season < CLIMATESEASONCOUNT; ++season)
        warmest = std::max(warmest, world.seasonaltempfloat(season,x,y));
    return climateenergy::permanentLandIceFraction(warmest);
}

inline int wrapx(int x, int width)
{
    if (x < 0 || x > width)
        return wrap(x, width);

    return x;
}

inline float latitudeforrow(int y, int height)
{
    constexpr double radiansToDegrees = 57.2957795130823208768;
    return static_cast<float>(
        climategrid::latitudeCentreRadians(y, height + 1) * radiansToDegrees);
}

inline double rowareaweight(int y, int height)
{
    if (height <= 0)
        return 0.0;

    constexpr double pi = 3.14159265358979323846;
    const int rows = height + 1;
    const double halfstep = 0.5 * pi / static_cast<double>(rows);
    return climategrid::latitudeBandMeasure(y, rows) /
        (2.0 * std::sin(halfstep));
}

inline float pressuresurfacetemperature(climateworld& world, int season, int x, int y)
{
    auto& climate = world.climatesimulation();
    if (tuning::climate::oceancurrents::oneWayDiagnosticsOnly)
        return world.seasonaltempfloat(season, x, y);
    if (world.sea(x, y) == 1 && climate.hasOcean &&
        climate.oceanSst[season].size() == static_cast<std::size_t>(world.width() + 1))
        return climate.oceanSkin[season][x][y];
    if (world.sea(x, y) == 1 && world.seasonalsstfloat(season, x, y) != 0.0f)
        return world.seasonalsstfloat(season, x, y);

    return world.seasonaltempfloat(season, x, y);
}

// climate_pressure.cpp
void smoothallfield(climateworld& world, floatgrid& field, int iterations);
