// climate.sea_surface_temperature: the sea surface temperature of the first guess.

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <vector>

#include "simulations/climate/climate_atmosphere.hpp"
#include "simulations/climate/climate_coupling.hpp"
#include "simulations/climate/climate_energy.hpp"
#include "core/climate_grid.hpp"
#include "simulations/climate/climate_hydrology.hpp"
#include "simulations/climate/climate_transport.hpp"
#include "simulations/climate/climate_ocean_dynamics.hpp"
#include "simulations/climate/detail/ocean_thermal_remap.hpp"
#include "simulations/climate/detail/ocean_replay.hpp"
#include "simulations/climate/detail/surface_wind_remap.hpp"
#include "simulations/climate/detail/surface_drag_mobility.hpp"
#include "simulations/climate/detail/katabatic_drainage.hpp"
#include "simulations/climate/detail/process_time_step.hpp"
#include "simulations/climate/detail/upper_thermal_projection.hpp"
#include "simulations/climate/detail/circulation_capture.hpp"
#include "simulations/climate/detail/mountain_wave_height.hpp"
#include "simulations/climate/detail/seasonal_overturning.hpp"
#include "simulations/climate/detail/thermal_trough.hpp"
#include "simulations/climate/detail/surface_momentum_advection.hpp"
#include "simulations/climate/detail/surface_850_exchange.hpp"
#include "simulations/climate/detail/moisture_carrier.hpp"
#include "simulations/climate/detail/lifted_parcel.hpp"
#include "simulations/climate/detail/surface_evaporation_humidity.hpp"
#include "simulations/climate/detail/frozen_surface_drag.hpp"
#include "simulations/climate/detail/vegetation_drag.hpp"
#include "io/climate_surface_input.hpp"
#include "io/climate_experiment_input.hpp"
#include "simulations/climate/climate_physics.hpp"
#include "simulations/climate/climate_weather.hpp"
#include "simulations/climate/climate_parameters.hpp"
#include "simulations/climate/climate_world.hpp"
#include "simulations/climate/climate_classification.hpp"
#include "simulations/climate/climate_fields.hpp"
#include "core/classes.hpp"
#include "simulations/climate/climate_elevation.hpp"
#include "core/deterministic_random.hpp"
#include "pipeline/generation_progress.hpp"
#include "core/grid_coordinates.hpp"
#include "simulations/terrain/land_distance.hpp"
#include "core/parallel_rows.hpp"
#include "core/procedural_noise.hpp"
#include "pipeline/world_generation_debug.hpp"
#include "simulations/climate/climate_coupling_internal.hpp"

using namespace std;

void createsurfacetemperaturemap(climateworld& world)
{
    auto& climate = world.climatesimulation();
    const int width = world.width();
    const int height = world.height();
    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
    {
        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int x = 0; x <= width; x++)
                {
                    if (!world.basebelowsea(x, y))
                    {
                        world.setseasonalsstfloat(season, x, y, 0.0f);
                        world.setseasonalevaporation(season, x, y, 0);
                        world.setseasonalmoisture(season, x, y, 0);
                        continue;
                    }
                    const float sst = tuning::climate::oceancurrents::oneWayDiagnosticsOnly
                        ? world.seasonaltempfloat(season, x, y)
                        : climate.oceanSkin[season][x][y];
                    const float eastwind = world.seasonaluwindfloat(season, x, y);
                    const float southwind = world.seasonalvwindfloat(season, x, y);
                    const float evaporation = climatephysics::bulkAerodynamicEvaporationMm(
                        sst,
                        0.0f,
                        std::max(1.0f, std::hypot(eastwind, southwind)),
                        0.75f,
                        tuning::climate::circulation::secondsPerDay,
                        tuning::climate::atmosphere::oceanMomentumDragCoefficient);
                    world.setseasonalevaporationfloat(season, x, y, evaporation);
                    world.setseasonalmoisturefloat(season, x, y, evaporation);
                }
            }
        });
    }
}
