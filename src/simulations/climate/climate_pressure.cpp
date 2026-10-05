// climate.pressure: surface pressure and the subtropical belts diagnosed from it.

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
#include "world/generation_tuning.hpp"
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

void smoothallfield(climateworld& world, floatgrid& field, int iterations)
{
    const int width = world.width();
    const int height = world.height();
    floatgrid scratch = field;

    for (int iteration = 0; iteration < iterations; iteration++)
    {
        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int x = 0; x <= width; x++)
                {
                    float total = 0.0f;
                    float weighttotal = 0.0f;

                    for (int dy = -1; dy <= 1; dy++)
                    {
                        const int yy = y + dy;

                        if (yy < 0 || yy > height)
                            continue;

                        for (int dx = -1; dx <= 1; dx++)
                        {
                            const int xx = wrapx(x + dx, width);
                            const float weight = (dx == 0 && dy == 0) ? 2.0f : 1.0f;
                            total += field[xx][yy] * weight;
                            weighttotal += weight;
                        }
                    }

                    scratch[x][y] = (weighttotal > 0.0f) ? total / weighttotal : field[x][y];
                }
            }
        });

        field.swap(scratch);
    }
}

void createpressuremap(climateworld& world)
{
    auto& climate = world.climatesimulation();
    climate.resetDynamics();
    climatephysics::clearClimateCouplingDiagnostics();
    const int width = world.width();
    const int height = world.height();
    const float thicknessresponse = climateatmosphere::hypsometricHeightResponseMetresPerKelvin(
        tuning::climate::circulation::surfaceReferencePressurePa,
        tuning::climate::circulation::upperReferencePressurePa);

    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
    {
        floatgrid surface(width + 1, vector<float>(height + 1, 0.0f));

        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int x = 0; x <= width; x++)
                    surface[x][y] = pressuresurfacetemperature(world, season, x, y);
            }
        });

        double temperaturesum = 0.0;
        double weighttotal = 0.0;

        for (int y = 0; y <= height; y++)
        {
            constexpr double pi = 3.14159265358979323846;
            const double weight = std::max(
                0.0,
                std::cos(static_cast<double>(latitudeforrow(y, height)) * pi / 180.0));

            for (int x = 0; x <= width; x++)
            {
                temperaturesum += weight * static_cast<double>(surface[x][y]);
                weighttotal += weight;
            }
        }

        const float globalmean = weighttotal > 0.0 ?
            static_cast<float>(temperaturesum / weighttotal) : 0.0f;
        floatgrid stationarytemperature = surface;
        smoothallfield(
            world,
            stationarytemperature,
            tuning::climate::circulation::thermalHeightSmoothingIterations);
        floatgrid pressure(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid upperheight(width + 1, vector<float>(height + 1, 0.0f));
        vector<float> zonaltemperature(height + 1, 0.0f);
        vector<float> stationaryzonaltemperature(height + 1, 0.0f);

        for (int y = 0; y <= height; y++)
        {
            double rowsum = 0.0;
            double stationaryrowsum = 0.0;

            for (int x = 0; x <= width; x++)
            {
                rowsum += static_cast<double>(surface[x][y]);
                stationaryrowsum += static_cast<double>(stationarytemperature[x][y]);
            }

            zonaltemperature[y] = static_cast<float>(rowsum / static_cast<double>(width + 1));
            stationaryzonaltemperature[y] = static_cast<float>(
                stationaryrowsum / static_cast<double>(width + 1));
        }

        for (int iteration = 0;
             iteration < tuning::climate::circulation::zonalTemperatureSmoothingIterations;
             iteration++)
        {
            vector<float> smoothed = zonaltemperature;

            for (int y = 0; y <= height; y++)
            {
                float total = 2.0f * zonaltemperature[y];
                float weights = 2.0f;

                if (y > 0)
                {
                    total += zonaltemperature[y - 1];
                    weights += 1.0f;
                }

                if (y < height)
                {
                    total += zonaltemperature[y + 1];
                    weights += 1.0f;
                }

                smoothed[y] = total / weights;
            }

            zonaltemperature.swap(smoothed);
        }

        const int thermalrow = static_cast<int>(std::distance(
            zonaltemperature.begin(),
            std::max_element(zonaltemperature.begin(), zonaltemperature.end())));
        const int polarcaprows = std::max(1, height / 12);
        double polartemperaturetotal = 0.0;
        int polartemperaturecount = 0;

        for (int y = 0; y <= polarcaprows; y++)
        {
            polartemperaturetotal += zonaltemperature[y];
            polartemperaturetotal += zonaltemperature[height - y];
            polartemperaturecount += 2;
        }

        const float polartemperature = polartemperaturecount > 0 ?
            static_cast<float>(polartemperaturetotal / polartemperaturecount) : globalmean;
        const float temperaturecontrast = std::max(
            1.0f,
            zonaltemperature[thermalrow] - polartemperature);
        const float gravity = tuning::climate::atmosphere::gravityMetresPerSecondSquared *
            std::max(0.05f, world.gravity());
        const float hadleyhalfwidth = std::clamp(
            climateatmosphere::heldHouHadleyEdgeLatitudeDegrees(
                temperaturecontrast,
                tuning::climate::circulation::troposphereHeightMetres,
                tuning::climate::circulation::referenceTemperatureK,
                gravity,
                tuning::climate::atmosphere::rotationRatePerSecond,
                tuning::climate::moistureadvection::referencePlanetRadiusMetres),
            tuning::climate::circulation::minimumHadleyHalfWidthDegrees,
            tuning::climate::circulation::maximumHadleyHalfWidthDegrees);
        const float maximumthermalequatorshift = hadleyhalfwidth *
            tuning::climate::circulation::maximumThermalEquatorShiftHadleyFraction;
        const float thermalequator = std::clamp(
            latitudeforrow(thermalrow, height),
            -maximumthermalequatorshift,
            maximumthermalequatorshift);
        const float pressureamplitude = std::clamp(
            tuning::climate::circulation::surfaceReferencePressurePa / 100.0f *
                temperaturecontrast / tuning::climate::circulation::referenceTemperatureK *
                tuning::climate::circulation::overturningMassRedistributionEfficiency,
            tuning::climate::circulation::minimumOverturningPressureAmplitudeHpa,
            tuning::climate::circulation::maximumOverturningPressureAmplitudeHpa);
        vector<float> zonalpressure(height + 1, 0.0f);
        double pressuretotal = 0.0;
        double pressureweighttotal = 0.0;

        for (int y = 0; y <= height; y++)
        {
            constexpr double pi = 3.14159265358979323846;
            const float latitude = latitudeforrow(y, height);
            const double areaweight = std::max(
                0.0,
                std::cos(static_cast<double>(latitude) * pi / 180.0));
            zonalpressure[y] = climateatmosphere::axisymmetricOverturningPressureAnomalyHpa(
                latitude,
                thermalequator,
                hadleyhalfwidth,
                pressureamplitude);
            pressuretotal += areaweight * zonalpressure[y];
            pressureweighttotal += areaweight;
        }

        const float meanpressure = pressureweighttotal > 0.0 ?
            static_cast<float>(pressuretotal / pressureweighttotal) : 0.0f;

        for (float& rowpressure : zonalpressure)
            rowpressure -= meanpressure;

        vector<float> localstationarytemperature(
            static_cast<size_t>(width + 1) * static_cast<size_t>(height + 1),
            0.0f);
        for (int y = 0; y <= height; y++)
        {
            for (int x = 0; x <= width; x++)
            {
                localstationarytemperature[
                    static_cast<size_t>(y) * (width + 1) + x] =
                    stationarytemperature[x][y] - stationaryzonaltemperature[y];
            }
        }
        cout
            << "Overturning pressure season " << season
            << ": thermal_equator=" << thermalequator
            << " C_contrast=" << temperaturecontrast
            << " hadley_half_width=" << hadleyhalfwidth
            << " pressure_amplitude_hpa=" << pressureamplitude << '\n';

        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int x = 0; x <= width; x++)
                {
                    pressure[x][y] = zonalpressure[y];
                    upperheight[x][y] = 0.0f;
                }
            }
        });

        smoothallfield(world, upperheight, tuning::climate::circulation::thermalHeightSmoothingIterations);

        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int x = 0; x <= width; x++)
                {
                    world.setseasonalpressurefloat(season, x, y, pressure[x][y]);
                    world.setseasonalupperheightfloat(season, x, y, upperheight[x][y]);
                }
            }
        });

    }
}

void updatehorsebeltsfrompressure(climateworld& world)
{
    const int width = world.width();
    const int height = world.height();

    auto rowfromlatitude = [height](float latitude)
    {
        const int rows = height + 1;
        return std::clamp(static_cast<int>(std::round(
            (90.0f - latitude) * static_cast<float>(rows) / 180.0f - 0.5f)),
            0, height);
    };

    const int defaultnorthouter = rowfromlatitude(32.0f);
    const int defaultnorthinner = rowfromlatitude(26.0f);
    const int defaultsouthinner = rowfromlatitude(-26.0f);
    const int defaultsouthouter = rowfromlatitude(-32.0f);
    const int minband = std::max(2, height / 60);

    vector<int> northouter(width + 1, defaultnorthouter);
    vector<int> northinner(width + 1, defaultnorthinner);
    vector<int> southinner(width + 1, defaultsouthinner);
    vector<int> southouter(width + 1, defaultsouthouter);

    for (int x = 0; x <= width; x++)
    {
        int northpeakpoleward = height;
        int northpeakequatorward = 0;
        int southpeakequatorward = height;
        int southpeakpoleward = 0;
        bool foundnorth = false;
        bool foundsouth = false;

        for (int season = 0; season < CLIMATESEASONCOUNT; season++)
        {
            float bestnorth = -1000000.0f;
            float bestsouth = -1000000.0f;
            int bestnorthy = defaultnorthinner;
            int bestsouthy = defaultsouthinner;

            for (int y = 0; y <= height; y++)
            {
                const float latitude = latitudeforrow(y, height);
                float pressure = 0.0f;

                for (int dx = -1; dx <= 1; dx++)
                {
                    const int xx = wrapx(x + dx, width);
                    pressure += world.seasonalpressurefloat(season, xx, y);
                }

                pressure = pressure / 3.0f;

                if (latitude >= 15.0f && latitude <= 45.0f)
                {
                    if (pressure > bestnorth)
                    {
                        bestnorth = pressure;
                        bestnorthy = y;
                    }
                }

                if (latitude <= -15.0f && latitude >= -45.0f)
                {
                    if (pressure > bestsouth)
                    {
                        bestsouth = pressure;
                        bestsouthy = y;
                    }
                }
            }

            northpeakpoleward = std::min(northpeakpoleward, bestnorthy);
            northpeakequatorward = std::max(northpeakequatorward, bestnorthy);
            southpeakequatorward = std::min(southpeakequatorward, bestsouthy);
            southpeakpoleward = std::max(southpeakpoleward, bestsouthy);
            foundnorth = true;
            foundsouth = true;
        }

        if (foundnorth)
        {
            if (northpeakequatorward - northpeakpoleward < minband)
            {
                const int centre = (northpeakequatorward + northpeakpoleward) / 2;
                northpeakpoleward = std::max(rowfromlatitude(45.0f), centre - minband / 2);
                northpeakequatorward = std::min(rowfromlatitude(15.0f), northpeakpoleward + minband);
            }

            northouter[x] = northpeakpoleward;
            northinner[x] = northpeakequatorward;
        }

        if (foundsouth)
        {
            if (southpeakpoleward - southpeakequatorward < minband)
            {
                const int centre = (southpeakpoleward + southpeakequatorward) / 2;
                southpeakequatorward = std::max(rowfromlatitude(-15.0f), centre - minband / 2);
                southpeakpoleward = std::min(rowfromlatitude(-45.0f), southpeakequatorward + minband);
            }

            southinner[x] = southpeakequatorward;
            southouter[x] = southpeakpoleward;
        }
    }

    auto smoothband = [width](vector<int>& band)
    {
        vector<int> smoothed = band;

        for (int x = 0; x <= width; x++)
        {
            int total = 0;
            int count = 0;

            for (int dx = -2; dx <= 2; dx++)
            {
                total += band[wrapx(x + dx, width)];
                count++;
            }

            smoothed[x] = total / count;
        }

        band.swap(smoothed);
    };

    smoothband(northouter);
    smoothband(northinner);
    smoothband(southinner);
    smoothband(southouter);

    for (int x = 0; x <= width; x++)
    {
        if (northinner[x] <= northouter[x])
            northinner[x] = std::min(height, northouter[x] + minband);

        if (southouter[x] <= southinner[x])
            southouter[x] = std::min(height, southinner[x] + minband);

        world.sethorse(x, 1, northouter[x]);
        world.sethorse(x, 2, northinner[x]);
        world.sethorse(x, 3, southinner[x]);
        world.sethorse(x, 4, southouter[x]);
    }
}
