// climate.rainfall_and_coupling: moisture transport and rainfall, and the coupled
// wind, pressure, ocean and rainfall iteration.

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

namespace
{
void smoothconvergencefootprint(
    const hydrologyfloatgrid& source,
    hydrologyfloatgrid& destination,
    hydrologyfloatgrid& scratch,
    const vector<float>& rowareas,
    int width,
    int height,
    float mixingfraction,
    int passes)
{
    const float mixing = std::clamp(mixingfraction, 0.0f, 1.0f);

    for (int pass = 0; pass < std::max(1, passes); pass++)
    {
        const hydrologyfloatgrid& current = pass == 0 ? source : destination;

        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int x = 0; x <= width; x++)
                {
                    scratch[x][y] =
                        (1.0f - mixing) * current[x][y] +
                        0.5f * mixing *
                            (current[wrapx(x - 1, width)][y] +
                                current[wrapx(x + 1, width)][y]);
                }
            }
        });

        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                const float centrearea = rowareas[y];

                for (int x = 0; x <= width; x++)
                {
                    const float centre = scratch[x][y];
                    float smoothed = centre;

                    if (y > 0)
                    {
                        smoothed += mixing * rowareas[y - 1] /
                            std::max(1.0e-6f, centrearea + rowareas[y - 1]) *
                            (scratch[x][y - 1] - centre);
                    }

                    if (y < height)
                    {
                        smoothed += mixing * rowareas[y + 1] /
                            std::max(1.0e-6f, centrearea + rowareas[y + 1]) *
                            (scratch[x][y + 1] - centre);
                    }

                    destination[x][y] = smoothed;
                }
            }
        });
    }
}

void exchangeeddymoisture(
    float& first,
    float& second,
    double firstarea,
    double secondarea,
    float exchangefraction);

void damphighwavenumbers(
    hydrologyfloatgrid& field,
    hydrologyfloatgrid& laplacian,
    hydrologyfloatgrid& next,
    const vector<float>& rowareas,
    const vector<float>& meridionalfactors,
    int width,
    int height,
    float strength,
    int passes,
    bool nonnegative)
{
    (void)laplacian;
    (void)nonnegative;
    const float damping = std::clamp(strength, 0.0f, 0.20f);

    for (int pass = 0; pass < std::max(0, passes); pass++)
    {
        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int x = 0; x <= width; x++)
                {
                    const float zonalmean = 0.5f * (
                        field[wrapx(x - 1, width)][y] +
                        field[wrapx(x + 1, width)][y]);
                    next[x][y] = field[x][y] +
                        damping * (zonalmean - field[x][y]);
                }
            }
        });

        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int x = 0; x <= width; x++)
                    field[x][y] = next[x][y];
            }
        });

        parallelforrows(0, width, [&](int startcolumn, int endcolumn)
        {
            for (int x = startcolumn; x <= endcolumn; x++)
            {
                for (int parity = 0; parity < 2; parity++)
                {
                    for (int y = parity; y < height; y += 2)
                    {
                        const float edgefactor = std::min(
                            meridionalfactors[y],
                            meridionalfactors[y + 1]);
                        exchangeeddymoisture(
                            field[x][y],
                            field[x][y + 1],
                            rowareas[y],
                            rowareas[y + 1],
                            0.5f * damping * edgefactor);
                    }
                }
            }
        });
    }
}

float samplehydrologygrid(
    const hydrologyfloatgrid& field,
    const climategrid::ConservativeRemap& remap,
    int targetX,
    int targetY)
{
    return climategrid::sampleConservative(remap, targetX, targetY,
        [&](int x, int y) { return field[x][y]; });
}

float eddyexchangefraction(
    float diffusivitym2s,
    double edgelengthmetres,
    double celldistancemetres,
    double firstcellareametres2,
    double secondcellareametres2,
    float timestepseconds)
{
    if (diffusivitym2s <= 0.0f || edgelengthmetres <= 0.0 || celldistancemetres <= 0.0 ||
        firstcellareametres2 <= 0.0 || secondcellareametres2 <= 0.0 || timestepseconds <= 0.0f)
    {
        return 0.0f;
    }

    const double coupling = static_cast<double>(diffusivitym2s) *
        edgelengthmetres / celldistancemetres;
    const double decayrate = coupling *
        (1.0 / firstcellareametres2 + 1.0 / secondcellareametres2);
    return static_cast<float>(std::clamp(
        1.0 - std::exp(-decayrate * static_cast<double>(timestepseconds)),
        0.0,
        1.0));
}

void exchangeeddymoisture(
    float& first,
    float& second,
    double firstarea,
    double secondarea,
    float exchangefraction)
{
    if (exchangefraction <= 0.0f || firstarea <= 0.0 || secondarea <= 0.0)
        return;

    const double equilibrium =
        (static_cast<double>(first) * firstarea + static_cast<double>(second) * secondarea) /
        (firstarea + secondarea);
    first = static_cast<float>(first + (equilibrium - first) * exchangefraction);
    second = static_cast<float>(second + (equilibrium - second) * exchangefraction);
}

void createeddyexchangefractions(
    int width,
    int height,
    const vector<float>& latitudes,
    const vector<float>& areaweights,
    const vector<float>& meridionaltransportfactors,
    const hydrologyfloatgrid& surfacewindu,
    const hydrologyfloatgrid& surfacewindv,
    const hydrologyfloatgrid& upperwindu,
    const hydrologyfloatgrid& upperwindv,
    float timestepseconds,
    hydrologyfloatgrid& zonalfractions,
    hydrologyfloatgrid& meridionalfractions)
{
    constexpr double pi = 3.14159265358979323846;
    const double radius = tuning::climate::moistureadvection::referencePlanetRadiusMetres;
    const double equatorialcellwidth = 2.0 * pi * radius / static_cast<double>(width + 1);
    const double meridionalcellheight = pi * radius / static_cast<double>(height + 1);
    const double equatorialcellarea = equatorialcellwidth * meridionalcellheight;
    hydrologyfloatgrid diffusivity(width + 1, height + 1, 0.0f);
    vector<double> rowareas(height + 1, 0.0);
    vector<double> zonalcelldistances(height + 1, 0.0);

    for (int y = 0; y <= height; y++)
    {
        constexpr double degreestoradians = pi / 180.0;
        const float latitude = latitudes[y];
        rowareas[y] = equatorialcellarea * areaweights[y];
        zonalcelldistances[y] = equatorialcellwidth *
            std::max(0.05, std::abs(std::cos(static_cast<double>(latitude) * degreestoradians)));

        for (int x = 0; x <= width; x++)
        {
            diffusivity[x][y] = std::min(
                tuning::climate::moistureadvection::maximumTransientEddyDiffusivityM2S,
                tuning::climate::moistureadvection::backgroundMoistureDiffusivityM2S +
                    climatephysics::transientEddyDiffusivityM2S(
                        surfacewindu[x][y],
                        surfacewindv[x][y],
                        upperwindu[x][y],
                        upperwindv[x][y],
                        latitude,
                        tuning::climate::moistureadvection::transientEddyMixingLengthMetres,
                        tuning::climate::moistureadvection::maximumTransientEddyDiffusivityM2S,
                        tuning::climate::moistureadvection::transientEddyMinimumLatitudeDegrees,
                        tuning::climate::moistureadvection::transientEddyFullStrengthLatitudeDegrees));
        }
    }

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int y = startrow; y <= endrow; y++)
        {
            for (int x = 0; x <= width; x++)
            {
                const int east = x < width ? x + 1 : 0;
                const float edgediffusivity = 0.5f * (diffusivity[x][y] + diffusivity[east][y]);
                zonalfractions[x][y] = eddyexchangefraction(
                    edgediffusivity,
                    meridionalcellheight,
                    zonalcelldistances[y],
                    rowareas[y],
                    rowareas[y],
                    timestepseconds);
            }
        }
    });

    parallelforrows(0, width, [&](int startcolumn, int endcolumn)
    {
        for (int x = startcolumn; x <= endcolumn; x++)
        {
            for (int y = 0; y < height; y++)
            {
                const float edgediffusivity = 0.5f * (diffusivity[x][y] + diffusivity[x][y + 1]);
                const double edgelength = equatorialcellwidth *
                    0.5 * (areaweights[y] + areaweights[y + 1]);
                meridionalfractions[x][y] = eddyexchangefraction(
                    edgediffusivity,
                    edgelength,
                    meridionalcellheight,
                    rowareas[y],
                    rowareas[y + 1],
                    timestepseconds) *
                    std::min(
                        meridionaltransportfactors[y],
                        meridionaltransportfactors[y + 1]);
            }
        }
    });
}

void applyeddyexchange(
    int width,
    int height,
    const vector<float>& areaweights,
    hydrologyfloatgrid& moisture,
    const hydrologyfloatgrid& zonalfractions,
    const hydrologyfloatgrid& meridionalfractions)
{
    vector<double> rowareas(height + 1, 0.0);

    for (int y = 0; y <= height; y++)
        rowareas[y] = areaweights[y];

    auto applyzonal = [&](bool halfstep)
    {
        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int parity = 0; parity < 2; parity++)
                {
                    for (int x = parity; x < width; x += 2)
                    {
                        float fraction = zonalfractions[x][y];

                        if (halfstep)
                            fraction = 1.0f - std::sqrt(std::max(0.0f, 1.0f - fraction));

                        exchangeeddymoisture(
                            moisture[x][y],
                            moisture[x + 1][y],
                            rowareas[y],
                            rowareas[y],
                            fraction);
                    }
                }

                float wrapfraction = zonalfractions[width][y];

                if (halfstep)
                    wrapfraction = 1.0f - std::sqrt(std::max(0.0f, 1.0f - wrapfraction));

                exchangeeddymoisture(
                    moisture[width][y],
                    moisture[0][y],
                    rowareas[y],
                    rowareas[y],
                    wrapfraction);
            }
        });
    };

    applyzonal(true);

    parallelforrows(0, width, [&](int startcolumn, int endcolumn)
    {
        for (int parity = 0; parity < 2; parity++)
        {
            for (int y = 1 + parity; y < height - 1; y += 2)
            {
                for (int x = startcolumn; x <= endcolumn; x++)
                {
                    exchangeeddymoisture(
                        moisture[x][y],
                        moisture[x][y + 1],
                        rowareas[y],
                        rowareas[y + 1],
                        meridionalfractions[x][y]);
                }
            }
        }
    });

    applyzonal(true);
}
}

void convergeclimatecoupling(climateworld& world, vector<vector<int>>& fractal)
{
    auto& climate = world.climatesimulation();
    using namespace tuning::climate::circulation;
    const int iterationLimit = climate.couplingIterationLimit > 0
        ? climate.couplingIterationLimit : maximumCouplingIterations;
    if (!climate.hasProcesses) return;
    const auto snapshot = [&]()
    {
        std::array<vector<double>, 4> values;
        for (int season = 0; season < CLIMATESEASONCOUNT; ++season)
        {
            const auto& fields = climate.processes[season];
            for (int y = 0; y < fields.rows; ++y)
                for (int x = 0; x < fields.columns; ++x)
                {
                    const int cell = y * fields.columns + x;
                    const int xx = std::min(world.width(), static_cast<int>((x + 0.5) * (world.width() + 1) / fields.columns));
                    const int yy = std::min(world.height(), static_cast<int>((y + 0.5) * (world.height() + 1) / fields.rows));
                    const double weight = std::sqrt(std::max(0.0, std::cos(3.141592653589793 * ((y + 0.5) / fields.rows - 0.5))));
                    values[0].push_back(weight * climate.winds.surfacewindu[season][xx][yy]);
                    values[0].push_back(weight * climate.winds.surfacewindv[season][xx][yy]);
                    values[0].push_back(weight * climate.winds.upperwindu[season][xx][yy]);
                    values[0].push_back(weight * climate.winds.upperwindv[season][xx][yy]);
                    if (world.sea(xx, yy)) values[1].push_back(weight * pressuresurfacetemperature(world, season, xx, yy));
                    for (int layer = 0; layer < 2; ++layer)
                        values[2].push_back(weight * (fields.radiativeHeatingWm2[layer][cell] + fields.latentHeatingWm2[layer][cell] +
                            (layer == 0 ? fields.sensibleHeatingWm2[cell] : 0.0f)));
                    values[2].push_back(weight * fields.surfaceNetHeatingWm2[cell]);
                    values[3].push_back(weight * world.seasonalrainfloat(season, xx, yy));
                }
        }
        return values;
    };
    const auto relative = [](const vector<double>& old, const vector<double>& now)
    {
        if (old.size() != now.size()) return std::numeric_limits<double>::infinity();
        double change = 0.0, scale = 0.0;
        for (std::size_t i = 0; i < old.size(); ++i)
        {
            if (!std::isfinite(old[i]) || !std::isfinite(now[i])) return std::numeric_limits<double>::infinity();
            change += std::pow(now[i] - old[i], 2);
            scale += std::max(old[i] * old[i], now[i] * now[i]);
        }
        return std::sqrt(change / std::max(static_cast<double>(std::max<std::size_t>(1, old.size())), scale));
    };
    auto old = snapshot();
    auto previousforcing = climate.processes;
    vector<double> previousinput, previousresidual;
    double relaxation = couplingHeatingRelaxation;
    for (int iteration = 1; iteration <= iterationLimit; ++iteration)
    {
        if (iteration > 1)
        {
            auto residual = snapshot()[2];
            for (std::size_t i = 0; i < residual.size(); ++i) residual[i] -= previousinput[i];
            relaxation = climatephysics::adaptiveCouplingRelaxation(relaxation, previousresidual, residual,
                couplingMinimumHeatingRelaxation, couplingHeatingRelaxation);
            previousresidual = std::move(residual);
            for (int season = 0; season < CLIMATESEASONCOUNT; ++season)
            {
                auto& fields = climate.processes[season];
                const auto blend = [&](vector<float>& now, const vector<float>& previous) {
                    for (std::size_t i = 0; i < now.size(); ++i)
                        now[i] = previous[i] + static_cast<float>(relaxation) * (now[i] - previous[i]); };
                for (int layer = 0; layer < 2; ++layer)
                {
                    blend(fields.radiativeHeatingWm2[layer], previousforcing[season].radiativeHeatingWm2[layer]);
                    blend(fields.latentHeatingWm2[layer], previousforcing[season].latentHeatingWm2[layer]);
                }
                blend(fields.sensibleHeatingWm2, previousforcing[season].sensibleHeatingWm2);
                blend(fields.surfaceRadiativeHeatingWm2, previousforcing[season].surfaceRadiativeHeatingWm2);
                blend(fields.surfaceLatentHeatLossWm2, previousforcing[season].surfaceLatentHeatLossWm2);
                blend(fields.surfaceNetHeatingWm2, previousforcing[season].surfaceNetHeatingWm2);
                blend(fields.surfaceTemperatureC, previousforcing[season].surfaceTemperatureC);
            }
        }
        previousforcing = climate.processes;
        // Test the generated heating against the forcing actually supplied to
        // this iteration, not against the previous unrelaxed diagnostic.
        old[2] = snapshot()[2];
        previousinput = old[2];
        createvectorwindmap(world);
        updatehorsebeltsfrompressure(world);
        createoceancurrentmap(world);
        createsurfacetemperaturemap(world);
        refreshadvectedrainfall(world, fractal);
        const auto now = snapshot();
        climatephysics::ClimateCouplingDiagnostics d;
        d.iteration = iteration;
        d.heatingRelaxation = iteration == 1 ? 1.0 : relaxation;
        d.windChange = relative(old[0], now[0]); d.sstChange = relative(old[1], now[1]);
        d.heatingChange = relative(old[2], now[2]); d.rainfallChange = relative(old[3], now[3]);
        std::array<double, 3> heatingChange{}, heatingScale{};
        std::size_t peakHeatingIndex = 0;
        double peakHeatingChange = 0.0;
        for (std::size_t i = 0; i < old[2].size(); ++i)
        {
            const double change = now[2][i] - old[2][i];
            heatingChange[i % 3] += change * change;
            heatingScale[i % 3] += std::max(old[2][i] * old[2][i], now[2][i] * now[2][i]);
            if (std::abs(change) > std::abs(peakHeatingChange))
            {
                peakHeatingChange = change;
                peakHeatingIndex = i;
            }
        }
        const int forcingColumns = std::max(1, climate.processes[0].columns);
        const int forcingRows = std::max(1, climate.processes[0].rows);
        const std::size_t peakCell = peakHeatingIndex / 3;
        d.innerSolvesAccepted = std::all_of(climate.oceanAccepted.begin(), climate.oceanAccepted.end(), [](bool accepted) { return accepted; }) &&
            std::all_of(climate.atmosphereAccepted.begin(), climate.atmosphereAccepted.end(), [](bool accepted) { return accepted; }) &&
            climatephysics::lastHydrologySpinupDiagnostics().converged && climate.oceanStorage.periodic;
        d.converged = climatephysics::climateCouplingConverged(d, minimumCouplingIterations, couplingRelativeTolerance);
        climatephysics::appendClimateCouplingDiagnostics(d);
        std::cout << "Climate coupling iteration=" << iteration << " wind_change=" << d.windChange << " sst_change=" << d.sstChange
            << " heating_change=" << d.heatingChange << " rainfall_change=" << d.rainfallChange
            << " heating_relaxation=" << d.heatingRelaxation
            << " heating_lower=" << std::sqrt(heatingChange[0] / std::max(1.0, heatingScale[0]))
            << " heating_upper=" << std::sqrt(heatingChange[1] / std::max(1.0, heatingScale[1]))
            << " heating_surface=" << std::sqrt(heatingChange[2] / std::max(1.0, heatingScale[2]))
            << " heating_peak_component=" << peakHeatingIndex % 3
            << " heating_peak_season=" << peakCell / (forcingColumns * forcingRows)
            << " heating_peak_lat=" << 90.0 - 180.0 * (peakCell / forcingColumns % forcingRows + 0.5) / forcingRows
            << " heating_peak_lon=" << -180.0 + 360.0 * (peakCell % forcingColumns + 0.5) / forcingColumns
            << " heating_peak_area_scaled_wm2=" << peakHeatingChange
            << " inner_accepted=" << d.innerSolvesAccepted << " converged=" << d.converged << '\n';
        if (d.converged) break;
        old = now;
    }
    const auto& diagnostics = climatephysics::lastClimateCouplingDiagnostics();
    if (!diagnostics.empty() && !diagnostics.back().converged)
        std::cout << "Climate coupling did not converge within " << iterationLimit
            << " iterations; inspect the climate_coupling diagnostics.\n";
}

void createadvectedrainfall(climateworld& world, vector<vector<int>>& inland, vector<vector<int>>& fractal)
{
    auto& climate = world.climatesimulation();
    (void)inland;
    (void)fractal;
    const bool diurnalCycle = climate.transientWeather &&
        tuning::climate::moistureadvection::enablePrescribedCoastalDiurnalCycle;
    const bool prescribedSynoptic = climate.transientWeather &&
        tuning::climate::moistureadvection::enablePrescribedSynopticPerturbations;
    const int outputwidth = world.width();
    const int outputheight = world.height();
    const climatehydrology::ClimateGridDimensions climatedimensions =
        climatehydrology::climateGridDimensions(
            outputwidth + 1,
            outputheight + 1,
            climate.resolution.hydrologyColumns);
    const int width = climatedimensions.columns - 1;
    const int height = climatedimensions.rows - 1;
    if (tuning::climate::moistureadvection::oceanEvaporationUseBoundaryHumidity)
        std::cout << "Ocean evaporation humidity=lower-reservoir proxy\n";
    const auto* referenceqa = referenceexperiment().field(width + 1, "qa");
    const auto* referenceevapmonthly = referenceexperiment().field(width + 1, "evap_m");
    const auto* referenceevap = referenceevapmonthly ? referenceevapmonthly : referenceexperiment().field(width + 1, "evap");
    const auto* referenceradiation = referenceexperiment().field(width + 1, "rad");
    const auto* referenceuppermodeeast = referenceexperiment().field(width + 1, "ut500");
    const auto* referenceuppermodenorth = referenceexperiment().field(width + 1, "vt500");
    if (bool(referenceuppermodeeast) != bool(referenceuppermodenorth))
        throw std::runtime_error("Reference upper transport mode requires both wind components");
    // Diagnostic ledgers are separate from prognostic stores and saved-world fields.
    std::array<std::vector<std::array<double, 5>>, CLIMATESEASONCOUNT> referenceledger;
    if (referenceqa || referenceevap || referenceradiation)
        for (auto& quarter : referenceledger) quarter.resize((width + 1) * (height + 1));
    std::array<std::vector<std::array<double, 3>>, CLIMATESEASONCOUNT> uppertransportledger;
    if (referenceuppermodeeast)
        for (auto& quarter : uppertransportledger) quarter.resize((width + 1) * (height + 1));
    const int sealevel = world.sealevel();
    const float gravitymultiplier = std::max(0.05f, world.gravity());
    std::array<std::vector<float>, 12> monthlyrain;
    const bool hasfloatingwinds = climate.winds.completeFor(world.width(), world.height());
    const int baselineiterations = tuning::climate::moistureadvection::iterations;
    const int referenceiterations = tuning::climateresolution::scaleDistance(
        baselineiterations, width, height);
    const float condensationconversiontimeseconds =
        tuning::climate::moistureadvection::condensationConversionTimeDays *
        tuning::climate::circulation::secondsPerDay;
    vector<float> rowareaweights(height + 1, 0.0f);
    vector<float> climatelatitudes(height + 1, 0.0f);
    vector<float> polartransportfactors(height + 1, 1.0f);
    vector<float> polarconvectionfactors(height + 1, 1.0f);

    for (int y = 0; y <= height; y++)
    {
        climatelatitudes[y] = climatehydrology::climateCellLatitudeDegrees(
            y,
            height + 1);
        rowareaweights[y] = climatehydrology::climateCellAreaWeight(
            y,
            height + 1);
        // Cell-face polar boundaries close transport; do not suppress valid
        // meridional transport over whole high-latitude bands.
        polarconvectionfactors[y] = climatehydrology::polarTaperFactor(
            climatelatitudes[y],
            tuning::climate::moistureadvection::polarConvectionTaperStartDegrees,
            tuning::climate::moistureadvection::polarConvectionTaperEndDegrees);
    }

    const auto inputremap = climategrid::makeConservativeRemap(
        outputwidth + 1, outputheight + 1, width + 1, height + 1);
    auto aggregateworldfield = [&](int x, int y, auto accessor)
    {
        double total = 0.0;
        for (const auto& yw : inputremap.latitude[y])
        {
            for (const auto& xw : inputremap.longitude[x])
                total += xw.fraction * yw.fraction * accessor(xw.source, yw.source);
        }
        return static_cast<float>(total);
    };
    auto aggregateworldmaximum = [&](int x, int y, auto accessor)
    {
        float maximum = 0.0f;
        for (const auto& yw : inputremap.latitude[y])
        {
            for (const auto& xw : inputremap.longitude[x])
                maximum = std::max(maximum, static_cast<float>(accessor(xw.source, yw.source)));
        }

        return maximum;
    };

    hydrologyfloatgrid landfraction(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid meanelevation(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid peakelevation(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid icefraction(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid landpermanentice(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid coastnormalu(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid coastnormalv(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid coastalstrength(width + 1, height + 1, 0.0f);
    hydrologybytegrid seacell(width + 1, height + 1, 0);

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int y = startrow; y <= endrow; y++)
        {
            for (int x = 0; x <= width; x++)
            {
                landfraction[x][y] = aggregateworldfield(
                    x,
                    y,
                    [&](int xx, int yy) { return world.sea(xx, yy) == 0 ? 1.0f : 0.0f; });
                meanelevation[x][y] = aggregateworldfield(
                    x,
                    y,
                    [&](int xx, int yy)
                    {
                        return std::max(0.0f, world.mapfloat(xx, yy) - static_cast<float>(sealevel));
                    });
                peakelevation[x][y] = aggregateworldmaximum(
                    x,
                    y,
                    [&](int xx, int yy)
                    {
                        return std::max(0.0f, world.mapfloat(xx, yy) - static_cast<float>(sealevel));
                    });
                icefraction[x][y] = aggregateworldfield(
                    x,
                    y,
                    [&](int xx, int yy) { return world.seaice(xx, yy) != 0 ? 1.0f : 0.0f; });
                seacell[x][y] = landfraction[x][y] < 0.5f ? 1 : 0;
                landpermanentice[x][y] = aggregateworldfield(x,y,
                    [&](int xx, int yy) { return worldpermanentlandicefraction(world,xx,yy); });
            }
        }
    });

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int y = startrow; y <= endrow; y++)
        {
            const int north = std::max(0, y - 1);
            const int south = std::min(height, y + 1);

            for (int x = 0; x <= width; x++)
            {
                const float gradienteast = 0.5f *
                    (landfraction[wrapx(x + 1, width)][y] -
                        landfraction[wrapx(x - 1, width)][y]);
                const float gradientsouth = 0.5f *
                    (landfraction[x][south] - landfraction[x][north]);
                const float magnitude = std::sqrt(
                    gradienteast * gradienteast + gradientsouth * gradientsouth);

                if (magnitude > 1.0e-5f)
                {
                    coastnormalu[x][y] = gradienteast / magnitude;
                    coastnormalv[x][y] = gradientsouth / magnitude;
                }

                coastalstrength[x][y] = std::clamp(
                    magnitude * 2.0f + 4.0f * landfraction[x][y] *
                        (1.0f - landfraction[x][y]),
                    0.0f,
                    1.0f);
            }
        }
    });

    std::cout
        << "Climate hydrology grid output=" << outputwidth + 1 << 'x' << outputheight + 1
        << " internal=" << width + 1 << 'x' << height + 1
        << " weather_phases="
        << (diurnalCycle ? tuning::climate::moistureadvection::weatherPhaseCount : 1)
        << " transport=" << (climate.transportMethod == climatecontext::TransportMethod::Ffsl ? "ffsl" : "mpdata")
        << " forcing=" << (climate.transientWeather ? "transient" : "daily-mean")
        << " circulation_background=" << (climate.transientWeather ? "updated" : "radiative-seasonal")
        << " remap=continuous-conservative"
        << " years_per_update=" << tuning::climate::moistureadvection::yearsPerCouplingUpdate
        << std::endl;

    hydrologyfloatgrid soilmoisture(
        width + 1,
        height + 1,
        tuning::climate::moistureadvection::landSoilMoistureCapacity *
            tuning::climate::moistureadvection::initialSoilMoistureFraction);
    climatecontext::RowMajorGrid<double> landSurface(width + 1, height + 1, 0.0);
    climatecontext::RowMajorGrid<double> landDeep(width + 1, height + 1, 0.0);
    climatecontext::RowMajorGrid<double> landInputEnergy(width + 1, height + 1, 0.0);
    for (int y = 0; y <= height; ++y)
        for (int x = 0; x <= width; ++x)
        {
            const double initial = aggregateworldfield(x, y, [&](int xx, int yy) {
                double total = 0.0;
                for (int season = 0; season < CLIMATESEASONCOUNT; ++season)
                    total += world.seasonaltempfloat(season, xx, yy) / 4.0;
                return total;
            });
            landSurface[x][y] = landDeep[x][y] = initial;
        }
    hydrologyfloatgrid snowwater(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid boundarymoisture(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid freemoisture(width + 1, height + 1, 0.0f);
    // Continue spin-up across coupled iterations, but never across worlds.
    // Reinitializing these stores discarded the very equilibrium being solved.
    if (climate.storage.columns == width + 1 && climate.storage.rows == height + 1)
    {
        boundarymoisture = climate.storage.boundary;
        freemoisture = climate.storage.free;
        soilmoisture = climate.storage.soil;
        snowwater = climate.storage.snow;
        landSurface = climate.storage.landSurfaceC;
        landDeep = climate.storage.landDeepC;
    }
    else
        climate.hydrologyYears = 0;
    hydrologyfloatgrid rawannualrain(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid quarterrain(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid quarterconvergence(width + 1, height + 1, 0.0f);
    std::array<hydrologyfloatgrid, CLIMATESEASONCOUNT> storedseasonalrain;
    std::array<hydrologyfloatgrid, CLIMATESEASONCOUNT> storedseasonalmoisture;
    std::array<hydrologyfloatgrid, CLIMATESEASONCOUNT> storedseasonalconvergence;

    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
    {
        storedseasonalrain[season] = hydrologyfloatgrid(width + 1, height + 1, 0.0f);
        storedseasonalmoisture[season] = hydrologyfloatgrid(width + 1, height + 1, 0.0f);
        storedseasonalconvergence[season] = hydrologyfloatgrid(width + 1, height + 1, 0.0f);
    }
    std::array<climatephysics::WaterBudget, CLIMATESEASONCOUNT> finalareaweightedbudgets{};
    std::array<climatephysics::WaterBudget, CLIMATESEASONCOUNT> finalbudgets{};
    std::array<climatephysics::CondensationActivityDiagnostics, CLIMATESEASONCOUNT>
        finalcondensationactivity{};
    std::array<climatephysics::PrecipitationProcessDiagnostics, CLIMATESEASONCOUNT>
        finalprocessdiagnostics{};

    struct monthlysolverresult
    {
        climatephysics::WaterBudget budget;
        climatephysics::WaterBudget areaweightedbudget;
        climatephysics::CondensationActivityDiagnostics condensation;
        climatephysics::PrecipitationProcessDiagnostics processes;
        std::vector<std::uint8_t> weatherstate;
    };
    std::vector<std::uint8_t> finalweatherstate;
    climateweather::ShallowWaterState continuedweather;
    std::array<vector<climateweather::ShallowWaterState>, CLIMATESEASONCOUNT> seasonalweather;
    std::array<vector<double>, CLIMATESEASONCOUNT> seasonalweights;
    std::array<climatehydrology::SeasonalProcessFields, CLIMATESEASONCOUNT> seasonalprocess;
    hydrologyfloatgrid sensibleheat(width + 1, height + 1, 0.0f);

    hydrologyfloatgrid totalrain(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid totalseaevaporation(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid totallandevaporation(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid runoff(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid uplift(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid descent(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid dynamicvertical(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid saturationcapacity(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid boundarycapacity(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid surfacepressure(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid temperature(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid surfacewindu(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid surfacewindv(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid upperwindu(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid upperwindv(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid surfacewindspeed(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid surfaceexchangeincrementu(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid surfaceexchangeincrementv(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid phaseexchangewindu(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid phaseexchangewindv(width + 1, height + 1, 0.0f);
    const float drainageTransportMean = tuning::climate::atmosphere::katabaticShallowTenMetreDiagnostic ?
        climateatmosphere::detail::coldJetTransportFraction(
            tuning::climate::atmosphere::katabaticShallowDepthMetres,
            tuning::climate::atmosphere::surfaceBoundaryLayerMomentumDepthMetres) :
        climateatmosphere::detail::katabaticLayerMeanFraction(
        tuning::climate::atmosphere::katabaticVelocityDecayHeightMetres,
        tuning::climate::atmosphere::surfaceBoundaryLayerMomentumDepthMetres);
    hydrologyfloatgrid phasewindu(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid phasewindv(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid phaseupperwindu(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid phaseupperwindv(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid previousphasewindu, previousphasewindv, previousphaseupperwindu, previousphaseupperwindv;
    hydrologyfloatgrid phasetemperature(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid phasesurfacetemperature(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid phasedynamicvertical(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid phasesaturationcapacity(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid phaseboundarycapacity(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid persistentheating(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid cloudmemory(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid surfaceevaporation(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid nextboundarymoisture(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid nextfreemoisture(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid boundaryfluxtendency(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid freefluxtendency(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid convectiveconvergence(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid convergencescratch(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid highpasslaplacian(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid highpassscratch(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid totalconvergence(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid zonaleddyexchange(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid meridionaleddyexchange(width + 1, height + 1, 0.0f);

    double transportSeconds = 0.0, weatherSeconds = 0.0;
    std::size_t transportCalls = 0, transportCellSubsteps = 0;
    int maximumTransportSubsteps = 0;
    auto runsolver = [&](int month)
    {
        const climatehydrology::CalendarMonth calendar =
            climatehydrology::calendarMonth(month);
        monthlyrain[month].assign((width + 1) * (height + 1), 0.0f);
        const int quarter = month / 3;
        auto& exchanges = seasonalprocess[quarter];
        int startday = 0;
        for (int m = 0; m < month; ++m) startday += climatehydrology::calendarMonth(m).days;
        const int stepsperday = std::max({
            static_cast<int>(std::ceil(1.0f / tuning::climate::moistureadvection::maximumClimateStepDays)),
            climate.transientWeather
                ? static_cast<int>(std::ceil(static_cast<float>(referenceiterations) / baselineiterations)) : 1,
            diurnalCycle
                ? tuning::climate::moistureadvection::weatherPhaseCount : 1});
        const int iterations = calendar.days * stepsperday;
        const float timestepseconds =
            static_cast<float>(calendar.days) *
            tuning::climate::moistureadvection::advectionTimeStepSeconds /
            static_cast<float>(iterations);
        totalrain.fill(0.0f);
        totalseaevaporation.fill(0.0f);
        totallandevaporation.fill(0.0f);
        runoff.fill(0.0f);
        totalconvergence.fill(0.0f);
        climatephysics::WaterBudget budget;
        climatephysics::WaterBudget areaweightedbudget;
        vector<climatephysics::CondensationActivityDiagnostics> rowcondensation(height + 1);
        vector<climatephysics::PrecipitationProcessDiagnostics> rowprocesses(height + 1);

        if (month % 3 == 0)
        {
            quarterrain.fill(0.0f);
            quarterconvergence.fill(0.0f);
        }

        for (int y = 0; y <= height; y++)
        {
            const double areaweight = rowareaweights[y];

            for (int x = 0; x <= width; x++)
            {
                const float atmosphericwater =
                    boundarymoisture[x][y] + freemoisture[x][y];
                budget.initialAtmosphericStorage += atmosphericwater;
                areaweightedbudget.initialAtmosphericStorage +=
                    areaweight * atmosphericwater;

                if (seacell[x][y] == 0)
                {
                    budget.initialSoilStorage += soilmoisture[x][y];
                    areaweightedbudget.initialSoilStorage += areaweight * soilmoisture[x][y];
                    budget.initialSnowStorage += snowwater[x][y];
                    areaweightedbudget.initialSnowStorage += areaweight * snowwater[x][y];
                }
            }
        }

        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int x = 0; x <= width; x++)
                {
                    auto aggregatefield = [&](int season, auto accessor)
                    {
                        return aggregateworldfield(
                            x,
                            y,
                            [&](int xx, int yy) { return accessor(season, xx, yy); });
                    };
                    auto interpolatefield = [&](auto accessor)
                    {
                        return climatehydrology::interpolateSeasonal(
                            aggregatefield(calendar.firstSeason, accessor),
                            aggregatefield(calendar.secondSeason, accessor),
                            calendar.interpolation);
                    };
                    auto interpolatecachedwind = [&](
                        const std::array<floatgrid, CLIMATESEASONCOUNT>& field)
                    {
                        return climatehydrology::interpolateSeasonal(
                            aggregateworldfield(
                                x,
                                y,
                                [&](int xx, int yy)
                                {
                                    return field[calendar.firstSeason][xx][yy];
                                }),
                            aggregateworldfield(
                                x,
                                y,
                                [&](int xx, int yy)
                                {
                                    return field[calendar.secondSeason][xx][yy];
                                }),
                            calendar.interpolation);
                    };
                    const float u = hasfloatingwinds
                        ? interpolatecachedwind(climate.winds.surfacewindu)
                        : interpolatefield(
                            [&](int season, int xx, int yy)
                            {
                                return world.seasonaluwindfloat(season, xx, yy);
                            });
                    const float v = hasfloatingwinds
                        ? interpolatecachedwind(climate.winds.surfacewindv)
                        : interpolatefield(
                            [&](int season, int xx, int yy)
                            {
                                return world.seasonalvwindfloat(season, xx, yy);
                            });
                    float upperu = hasfloatingwinds
                        ? interpolatecachedwind(climate.winds.upperwindu)
                        : interpolatefield(
                            [&](int season, int xx, int yy)
                            {
                                return world.seasonalupperuwindfloat(season, xx, yy);
                            });
                    float upperv = hasfloatingwinds
                        ? interpolatecachedwind(climate.winds.upperwindv)
                        : interpolatefield(
                            [&](int season, int xx, int yy)
                            {
                                return world.seasonaluppervwindfloat(season, xx, yy);
                            });
                    const float upperwinddeparturemultiplier = std::clamp(
                        (tuning::climate::moistureadvection::
                            polarUpperWindBlendTaperEndDegrees -
                            std::fabs(climatelatitudes[y])) /
                            std::max(
                                0.1f,
                                tuning::climate::moistureadvection::
                                    polarUpperWindBlendTaperEndDegrees -
                                tuning::climate::moistureadvection::
                                    polarUpperWindBlendTaperStartDegrees),
                        0.0f,
                        1.0f) *
                        tuning::climate::moistureadvection::
                            freeTroposphereUpperWindDepartureFraction;
                    const bool sea = seacell[x][y] != 0;
                    if (sea && climate.hasOcean && !tuning::climate::oceancurrents::oneWayDiagnosticsOnly)
                        icefraction[x][y] = interpolatefield([&](int season, int xx, int yy)
                            { return climate.oceanIce[season][xx][yy]; });
                    const float localtemperature = interpolatefield(
                        [&](int season, int xx, int yy)
                        {
                            return pressuresurfacetemperature(world, season, xx, yy);
                        });
                    const float elevation = meanelevation[x][y];
                    const bool hasLowLevel = hasfloatingwinds &&
                        climate.winds.hasLowLevelProfile[calendar.firstSeason] &&
                        climate.winds.hasLowLevelProfile[calendar.secondSeason];
                    const float lowAvailability = hasLowLevel
                        ? interpolatecachedwind(climate.winds.lowlevelavailable) : 0.0f;
                    const float lowWeight = climateatmosphere::detail::lowLevelMoistureWeight(
                        tuning::climate::moistureadvection::freeTroposphereLowLevelWindFraction,
                        lowAvailability, elevation, localtemperature + 273.15f, 9.80665f * gravitymultiplier);
                    if (referenceuppermodeeast)
                    {
                        // Replace only the 500 hPa component of the existing moisture
                        // carrier blend. Momentum diagnosis and surface exchange retain
                        // their own winds; the lower/upper moisture partition is unchanged.
                        const float height500 = 287.05f * (localtemperature + 273.15f) /
                            (9.80665f * gravitymultiplier) * std::log(100000.0f / 50000.0f);
                        const float missing = std::numeric_limits<float>::quiet_NaN();
                        const float east = climateio::experiment::sample(referenceuppermodeeast,
                            calendar.firstSeason, calendar.secondSeason, calendar.interpolation, x, y, missing);
                        const float north = climateio::experiment::sample(referenceuppermodenorth,
                            calendar.firstSeason, calendar.secondSeason, calendar.interpolation, x, y, missing);
                        if (elevation < height500 && std::isfinite(east) && std::isfinite(north))
                        {
                            auto& ledger = uppertransportledger[quarter][y * (width + 1) + x];
                            const double seconds = calendar.days * 86400.0;
                            ledger[0] += seconds;
                            ledger[1] += seconds * upperwinddeparturemultiplier * (1.0f - lowWeight) * (east - upperu);
                            ledger[2] += seconds * upperwinddeparturemultiplier * (1.0f - lowWeight) *
                                (north + upperv) * polartransportfactors[y];
                            upperu = east;
                            upperv = -north;
                        }
                    }
                    const float pressureanomaly = interpolatefield(
                        [&](int season, int xx, int yy)
                        {
                            return world.seasonalpressurefloat(season, xx, yy);
                        });

                    temperature[x][y] = localtemperature;
                    surfacepressure[x][y] = std::max(
                        100.0f,
                        climatephysics::surfacePressureHpa(
                            elevation,
                            gravitymultiplier) + pressureanomaly);
                    const float boundarytransportwindfraction = 1.0f -
                        landfraction[x][y] *
                            (1.0f - tuning::climate::moistureadvection::
                                landBoundaryLayerTransportWindFraction);
                    const auto drainage = [&](const auto& field)
                    {
                        if (!hasfloatingwinds || !tuning::climate::atmosphere::separateKatabaticTransport) return 0.0f;
                        return interpolatefield([&](int s, int xx, int yy)
                        { return climate.winds.hasDrainageProfile[s] ? field[s][xx][yy] : 0.0f; });
                    };
                    const float du = drainage(climate.winds.drainagewindu);
                    const float dv = drainage(climate.winds.drainagewindv);
                    if (lowWeight > 0.0f && lowAvailability > 0.0f)
                    {
                        // Normalize masked remapping near level/terrain intersections.
                        // The availability weight then makes the transition continuous.
                        upperu = climateatmosphere::detail::moistureCarrierMode(upperu,
                            interpolatecachedwind(climate.winds.lowlevelwindu) / lowAvailability, lowWeight);
                        upperv = climateatmosphere::detail::moistureCarrierMode(upperv,
                            interpolatecachedwind(climate.winds.lowlevelwindv) / lowAvailability, lowWeight);
                    }
                    surfacewindu[x][y] = (u - (1.0f - drainageTransportMean) * du) * boundarytransportwindfraction;
                    surfacewindv[x][y] = (v - (1.0f - drainageTransportMean) * dv) * boundarytransportwindfraction;
                    // The profile vanishes above the boundary layer. Preserve
                    // the carrier/upper blend, including the existing polar taper.
                    upperwindu[x][y] = (u - du) +
                        (upperu - (u - du)) * upperwinddeparturemultiplier;
                    upperwindv[x][y] = (v - dv) +
                        (upperv - (v - dv)) * upperwinddeparturemultiplier;
                    surfaceexchangeincrementu[x][y] = (1.0f - drainageTransportMean) * du * boundarytransportwindfraction;
                    surfaceexchangeincrementv[x][y] = (1.0f - drainageTransportMean) * dv * boundarytransportwindfraction;
                    saturationcapacity[x][y] =
                        climatephysics::saturationColumnWaterAtPressure(
                            localtemperature,
                            surfacepressure[x][y],
                            gravitymultiplier);
                    boundarycapacity[x][y] = saturationcapacity[x][y] *
                        tuning::climate::moistureadvection::boundaryLayerCapacityFraction;
                    surfacewindspeed[x][y] = std::max(
                        tuning::climate::moistureadvection::minimumSurfaceWind,
                        std::sqrt(u * u + v * v));
                    auto interpolateorographicextreme = [&](auto accessor)
                    {
                        auto retainextreme = [&](int season)
                        {
                            const float mean = aggregatefield(season, accessor);
                            const float maximum = aggregateworldmaximum(
                                x,
                                y,
                                [&](int xx, int yy) { return accessor(season, xx, yy); });
                            return mean + (maximum - mean) *
                                tuning::climate::moistureadvection::
                                    subgridOrographicExtremeRetention;
                        };
                        return climatehydrology::interpolateSeasonal(
                            retainextreme(calendar.firstSeason),
                            retainextreme(calendar.secondSeason),
                            calendar.interpolation);
                    };
                    uplift[x][y] = interpolateorographicextreme(
                        [&](int season, int xx, int yy)
                        {
                            return world.seasonalupliftfloat(season, xx, yy);
                        });
                    descent[x][y] = interpolateorographicextreme(
                        [&](int season, int xx, int yy)
                        {
                            return world.seasonalsubsidencefloat(season, xx, yy);
                        });
                    dynamicvertical[x][y] = interpolatefield(
                        [&](int season, int xx, int yy)
                        {
                            return world.seasonalverticalvelocityfloat(season, xx, yy);
                        });
                }
            }
        });

        createeddyexchangefractions(
            width,
            height,
            climatelatitudes,
            rowareaweights,
            polartransportfactors,
            surfacewindu,
            surfacewindv,
            upperwindu,
            upperwindv,
            timestepseconds * static_cast<float>(
                tuning::climate::moistureadvection::freeTroposphereTransportCadence),
            zonaleddyexchange,
            meridionaleddyexchange);

        auto advectlayer = [&] (
            const hydrologyfloatgrid& source,
            const hydrologyfloatgrid& windu,
            const hydrologyfloatgrid& windv,
            hydrologyfloatgrid& destination,
            float transporttimestepseconds)
        {
            const int columns = width + 1;
            const int rows = height + 1;
            const std::size_t cellcount = static_cast<std::size_t>(columns) * rows;
            const bool boundary = &source == &boundarymoisture;
            const auto& startu = boundary ? previousphasewindu : previousphaseupperwindu;
            const auto& startv = boundary ? previousphasewindv : previousphaseupperwindv;
            climatehydrology::MpdataOptions options;
            options.endZonalWindMps = windu.values();
            options.endMeridionalWindMps = windv.values();
            options.maximumCourantPerSubstep =
                tuning::climate::moistureadvection::transportMaximumCourant;
            options.correctivePasses =
                tuning::climate::moistureadvection::transportCorrectivePasses;
            options.monotone = true;
            std::vector<float> transported;
            const auto transportStart = std::chrono::steady_clock::now();
            climatehydrology::SphericalTracerTransportDiagnostics flux;
            if (climate.transportMethod == climatecontext::TransportMethod::Ffsl)
            {
                climatehydrology::FfslOptions ffsl;
                ffsl.endZonalWindMps = windu.values(); ffsl.endMeridionalWindMps = windv.values();
                flux = climatehydrology::advectSphericalTracerFfsl(columns, rows, source.values(),
                    startu.values(), startv.values(), transporttimestepseconds,
                    tuning::climate::moistureadvection::referencePlanetRadiusMetres, ffsl, transported);
            }
            else flux = climatehydrology::advectSphericalTracerMpdata(
                columns,
                rows,
                source.values(),
                startu.values(),
                startv.values(),
                transporttimestepseconds,
                tuning::climate::moistureadvection::referencePlanetRadiusMetres,
                options,
                transported);
            transportSeconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - transportStart).count();
            ++transportCalls;
            transportCellSubsteps += flux.cellSubsteps;
            maximumTransportSubsteps = std::max(maximumTransportSubsteps, flux.substeps);
            const int layer = &source == &boundarymoisture ? 0 : 1;
            for (std::size_t cell = 0; cell < cellcount; ++cell)
            {
                exchanges.eastIntegratedFlux[layer][cell] += flux.eastIntegratedFlux[cell];
                exchanges.southIntegratedFlux[layer][cell] += flux.southIntegratedFlux[cell];
            }
            destination.values() = std::move(transported);
        };

        const int weathercolumns = std::min(
            width + 1,
            climate.resolution.weatherColumns);
        const int weatherrows = std::max(2, weathercolumns / 2);
        climateweather::ShallowWaterConfig weatherconfig;
        weatherconfig.layerCount = 1;
        weatherconfig.rotationDirection = world.rotation() ? 1.0f : -1.0f;
        weatherconfig.gravityMetresPerSecondSquared *= std::max(0.05f, world.gravity());
        weatherconfig.lowerMeanDepthMetres =
            tuning::climate::moistureadvection::weatherAnomalyEquivalentDepthMetres;
        weatherconfig.lowerDragTimeSeconds =
            tuning::climate::moistureadvection::weatherAnomalyDragTimeDays *
            tuning::climate::circulation::secondsPerDay;
        weatherconfig.heightRelaxationTimeSeconds =
            tuning::climate::moistureadvection::weatherAnomalyRelaxationTimeDays *
            tuning::climate::circulation::secondsPerDay;
        weatherconfig.stochasticHeightForcingMetresPerSecond =
            tuning::climate::moistureadvection::weatherAnomalyStochasticHeightMps;
        auto weatherstate = climateweather::makeState(
            weathercolumns,
            weatherrows,
            1,
            static_cast<std::uint64_t>(world.seed()) ^
                (static_cast<std::uint64_t>(month + 1) * 0x9e3779b97f4a7c15ULL));
        if (continuedweather.columns == weathercolumns && continuedweather.rows == weatherrows)
            weatherstate = continuedweather;
        climateweather::ShallowWaterForcing weatherforcing;
        const auto weatherbackground = [&](const hydrologyfloatgrid& field)
        {
            vector<float> flat((width + 1) * (height + 1));
            for (int y = 0; y <= height; ++y)
                for (int x = 0; x <= width; ++x) flat[y * (width + 1) + x] = field[x][y];
            return climategrid::remapField(
                width + 1, height + 1, flat, weathercolumns, weatherrows);
        };
        weatherforcing.backgroundEastWindMps = {weatherbackground(surfacewindu)};
        weatherforcing.backgroundSouthWindMps = {weatherbackground(surfacewindv)};
        bool evolvingweather =
            climate.transientWeather && tuning::climate::moistureadvection::enableEvolvingWeatherAnomalies;
        const auto sampleweather = [&](const std::vector<float>& field, int x, int y)
        {
            const float sourcex = (static_cast<float>(x) + 0.5f) * weathercolumns /
                static_cast<float>(width + 1) - 0.5f;
            const float sourcey = std::clamp((static_cast<float>(y) + 0.5f) * weatherrows /
                static_cast<float>(height + 1) - 0.5f, 0.0f, static_cast<float>(weatherrows - 1));
            const int basex = static_cast<int>(std::floor(sourcex));
            const int x0 = ((basex % weathercolumns) + weathercolumns) % weathercolumns;
            const int x1 = (x0 + 1) % weathercolumns;
            const int y0 = std::clamp(
                static_cast<int>(std::floor(sourcey)), 0, weatherrows - 1);
            const int y1 = std::min(weatherrows - 1, y0 + 1);
            const float fractionx = sourcex - std::floor(sourcex);
            const float fractiony = std::clamp(
                sourcey - std::floor(sourcey), 0.0f, 1.0f);
            const auto at = [&](int xx, int yy)
            {
                return field[static_cast<std::size_t>(yy) * weathercolumns + xx];
            };
            const float north = at(x0, y0) * (1.0f - fractionx) +
                at(x1, y0) * fractionx;
            const float south = at(x0, y1) * (1.0f - fractionx) +
                at(x1, y1) * fractionx;
            return north * (1.0f - fractiony) + south * fractiony;
        };

        auto climatesample = climateweather::makeState(width + 1, height + 1, 2);
        const auto firstStepStart = std::chrono::steady_clock::now();
        for (int iteration = 0; iteration < iterations; iteration++)
        {
            const auto startweather = weatherstate;
            exchanges.durationSeconds += timestepseconds;
            const float dayofyear = startday + (iteration + 0.5f) * calendar.days / iterations;
            const auto solar = climateenergy::solarForcing(dayofyear, world.tilt(), world.eccentricity(), world.perihelion());
            const float solardeclination = solar.declinationRadians, solardistance = solar.distanceFactor;
            if (evolvingweather)
            {
                const auto weatherStart = std::chrono::steady_clock::now();
                const auto weatherdiagnostics = climateweather::advance(
                    weatherstate,
                    weatherconfig,
                    weatherforcing,
                    timestepseconds);
                evolvingweather = weatherdiagnostics.finite && weatherdiagnostics.bounded;
                weatherSeconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - weatherStart).count();
            }
            if (month == 0 && iteration == 0 && isworldgendebugrunactive())
                std::cout << "Climate first-step weather seconds=" << std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - firstStepStart).count() << std::endl;
            const auto buildphase = [&](const climateweather::ShallowWaterState& phaseweather, double phaseseconds)
            {
            parallelforrows(0, height, [&](int startrow, int endrow)
            {
                for (int y = startrow; y <= endrow; y++)
                {
                    constexpr float pi = 3.14159265358979323846f;
                    const float latitude = climatelatitudes[y] * pi / 180.0f;
                    const float latitudetaper = std::max(0.0f, std::cos(latitude));

                    for (int x = 0; x <= width; x++)
                    {
                        const climatehydrology::WeatherPhase weatherphase =
                            climatehydrology::weatherPhaseAtTime(
                                phaseseconds, -180.0f + 360.0f * (x + 0.5f) / (width + 1),
                                prescribedSynoptic
                                    ? tuning::climate::moistureadvection::synopticWindRotationDegrees : 0.0f,
                                tuning::climate::moistureadvection::daytimeLandTemperatureAnomalyC,
                                tuning::climate::moistureadvection::nighttimeLandTemperatureAnomalyC,
                                tuning::climate::moistureadvection::daytimeSeaTemperatureAnomalyC,
                                tuning::climate::moistureadvection::nighttimeSeaTemperatureAnomalyC);
                        const float rotationcos = std::cos(weatherphase.windRotationRadians);
                        const float rotationsin = std::sin(weatherphase.windRotationRadians);
                        const float longitude = 2.0f * pi * static_cast<float>(x) /
                            static_cast<float>(width + 1);
                        const float synopticu = evolvingweather
                            ? std::clamp(
                                sampleweather(phaseweather.layers[0].eastWindMps, x, y),
                                -tuning::climate::moistureadvection::
                                    maximumWeatherAnomalyWindMps,
                                tuning::climate::moistureadvection::
                                    maximumWeatherAnomalyWindMps)
                            : (prescribedSynoptic ? 1.0f : 0.0f) * tuning::climate::moistureadvection::synopticWindPerturbationMps *
                                std::sin(
                                    2.0f * longitude + weatherphase.synopticPhaseRadians) *
                                latitudetaper;
                        const float synopticv = evolvingweather
                            ? std::clamp(
                                sampleweather(phaseweather.layers[0].southWindMps, x, y),
                                -tuning::climate::moistureadvection::
                                    maximumWeatherAnomalyWindMps,
                                tuning::climate::moistureadvection::
                                    maximumWeatherAnomalyWindMps)
                            : (prescribedSynoptic ? 1.0f : 0.0f) * tuning::climate::moistureadvection::synopticWindPerturbationMps *
                                std::cos(
                                    3.0f * longitude - weatherphase.synopticPhaseRadians) *
                                latitudetaper;
                        const float coastalmagnitude = weatherphase.coastalDirection >= 0.0f
                            ? tuning::climate::moistureadvection::daytimeOnshoreCoastalWindMps
                            : tuning::climate::moistureadvection::nighttimeOffshoreCoastalWindMps;
                        const float coastalwind = tuning::climate::moistureadvection::enablePrescribedCoastalDiurnalCycle
                            ? coastalstrength[x][y] * (diurnalCycle ? weatherphase.coastalDirection * coastalmagnitude
                                : climatehydrology::dailyMeanCycle(
                                    tuning::climate::moistureadvection::daytimeOnshoreCoastalWindMps,
                                    -tuning::climate::moistureadvection::nighttimeOffshoreCoastalWindMps)) : 0.0f;
                        const float rotatedu =
                            surfacewindu[x][y] * rotationcos - surfacewindv[x][y] * rotationsin;
                        const float rotatedv =
                            surfacewindu[x][y] * rotationsin + surfacewindv[x][y] * rotationcos;
                        phasewindu[x][y] = rotatedu + synopticu +
                            coastalwind * coastnormalu[x][y];
                        phasewindv[x][y] = (
                            rotatedv + synopticv + coastalwind * coastnormalv[x][y]) *
                            polartransportfactors[y];
                        phaseupperwindu[x][y] =
                            upperwindu[x][y] * rotationcos - upperwindv[x][y] * rotationsin +
                            0.5f * synopticu;
                        phaseupperwindv[x][y] = (
                            upperwindu[x][y] * rotationsin + upperwindv[x][y] * rotationcos +
                            0.5f * synopticv) * polartransportfactors[y];
                        const float diurnalanomaly = tuning::climate::moistureadvection::enablePrescribedCoastalDiurnalCycle
                            ? (diurnalCycle ? weatherphase.seaTemperatureAnomalyC * (1.0f - landfraction[x][y]) +
                                weatherphase.landTemperatureAnomalyC * landfraction[x][y]
                                : climatehydrology::dailyMeanCycle(tuning::climate::moistureadvection::daytimeSeaTemperatureAnomalyC,
                                    tuning::climate::moistureadvection::nighttimeSeaTemperatureAnomalyC) * (1.0f - landfraction[x][y]) +
                                    climatehydrology::dailyMeanCycle(tuning::climate::moistureadvection::daytimeLandTemperatureAnomalyC,
                                        tuning::climate::moistureadvection::nighttimeLandTemperatureAnomalyC) * landfraction[x][y]) : 0.0f;
                        // Atmospheric latent/cloud adjustments do not directly
                        // change the independently integrated surface or ocean ice skin.
                        phasesurfacetemperature[x][y] = seacell[x][y]
                            ? temperature[x][y] + diurnalanomaly
                            : static_cast<float>(landSurface[x][y]);
                        phasetemperature[x][y] = temperature[x][y] + diurnalanomaly +
                            persistentheating[x][y] -
                            cloudmemory[x][y] *
                                tuning::climate::moistureadvection::maximumCloudRadiativeCoolingC;
                        phasedynamicvertical[x][y] =
                            dynamicvertical[x][y] * polartransportfactors[y] +
                            (tuning::climate::circulation::enableLaggedDiabaticCoupling ? 0.0f : persistentheating[x][y]) *
                                tuning::climate::moistureadvection::heatingVerticalVelocityResponse *
                                polarconvectionfactors[y];
                        phasesaturationcapacity[x][y] =
                            climatephysics::saturationColumnWaterAtPressure(
                                phasetemperature[x][y],
                                surfacepressure[x][y],
                                gravitymultiplier);
                        phaseboundarycapacity[x][y] = phasesaturationcapacity[x][y] *
                            tuning::climate::moistureadvection::boundaryLayerCapacityFraction;
                        // Scalar surface exchange retains the full local jet;
                        // its layer-mean transport velocity may be much smaller.
                        phaseexchangewindu[x][y] = phasewindu[x][y] +
                            surfaceexchangeincrementu[x][y] * rotationcos - surfaceexchangeincrementv[x][y] * rotationsin;
                        phaseexchangewindv[x][y] = phasewindv[x][y] +
                            (surfaceexchangeincrementu[x][y] * rotationsin + surfaceexchangeincrementv[x][y] * rotationcos) *
                            polartransportfactors[y];
                        surfacewindspeed[x][y] = std::max(
                            tuning::climate::moistureadvection::minimumSurfaceWind,
                            std::sqrt(
                                phaseexchangewindu[x][y] * phaseexchangewindu[x][y] +
                                phaseexchangewindv[x][y] * phaseexchangewindv[x][y]));
                    }
                }
            });

            };
            buildphase(startweather, startday * 86400.0 + iteration * timestepseconds);
            previousphasewindu = phasewindu; previousphasewindv = phasewindv;
            previousphaseupperwindu = phaseupperwindu; previousphaseupperwindv = phaseupperwindv;
            buildphase(weatherstate, startday * 86400.0 + (iteration + 1) * timestepseconds);
            climatesample.elapsedSeconds = startday * 86400.0 + (iteration / stepsperday + 0.5) * 86400.0;
            for (int y = 0; y <= height; ++y)
                for (int x = 0; x <= width; ++x)
                {
                    const int cell = y * (width + 1) + x;
                    climatesample.layers[0].eastWindMps[cell] += 0.5f * (previousphasewindu[x][y] + phasewindu[x][y]) / stepsperday;
                    climatesample.layers[0].southWindMps[cell] += 0.5f * (previousphasewindv[x][y] + phasewindv[x][y]) / stepsperday;
                    climatesample.layers[1].eastWindMps[cell] += 0.5f * (previousphaseupperwindu[x][y] + phaseupperwindu[x][y]) / stepsperday;
                    climatesample.layers[1].southWindMps[cell] += 0.5f * (previousphaseupperwindv[x][y] + phaseupperwindv[x][y]) / stepsperday;
                }
            if ((iteration + 1) % stepsperday == 0)
            {
                seasonalweather[quarter].push_back(climatesample);
                seasonalweights[quarter].push_back(86400.0);
                for (auto& layer : climatesample.layers)
                {
                    std::fill(layer.eastWindMps.begin(), layer.eastWindMps.end(), 0.0f);
                    std::fill(layer.southWindMps.begin(), layer.southWindMps.end(), 0.0f);
                }
            }

            parallelforrows(0, height, [&](int startrow, int endrow)
            {
                for (int y = startrow; y <= endrow; y++)
                {
                    for (int x = 0; x <= width; x++)
                    {
                        const bool sea = seacell[x][y] != 0;
                        const float localtemperature = phasetemperature[x][y];
                        const float surfacetemperature = phasesurfacetemperature[x][y];
                        const float retained = std::max(
                            0.0f,
                            boundarymoisture[x][y] + freemoisture[x][y]);
                        // The humidity deficit uses surface saturation and the
                        // actual vapour store, independently of parcel cooling.
                        const float surfacecapacity = climatephysics::saturationColumnWaterAtPressure(
                            surfacetemperature, surfacepressure[x][y], gravitymultiplier);
                        float originalhumidity = std::clamp(retained / surfacecapacity, 0.0f, 1.0f);
                        if (sea && tuning::climate::moistureadvection::oceanEvaporationUseBoundaryHumidity)
                            originalhumidity = climatephysics::detail::boundarySurfaceRelativeHumidity(
                                boundarymoisture[x][y], surfacecapacity,
                                tuning::climate::moistureadvection::boundaryLayerCapacityFraction);
                        float relativehumidity = originalhumidity;
                        if (sea && referenceqa)
                        {
                            const float qa = climateio::experiment::sample(referenceqa,
                                calendar.firstSeason, calendar.secondSeason, calendar.interpolation,
                                x, y, std::numeric_limits<float>::quiet_NaN());
                            if (std::isfinite(qa))
                            {
                                const float saturation = climatephysics::saturationSpecificHumidity(
                                    surfacetemperature, surfacepressure[x][y]);
                                relativehumidity = std::clamp(qa / std::max(1.0e-8f, saturation), 0.0f, 1.0f);
                                auto& ledger = referenceledger[quarter][y * (width + 1) + x];
                                ledger[0] += timestepseconds;
                                ledger[1] += qa * timestepseconds;
                                ledger[2] += qa >= saturation ? timestepseconds : 0.0;
                            }
                        }
                        const float airtemperature = localtemperature;
                        // Scalar turbulent fluxes need wind variance missing
                        // from the daily baseline mean vector. Leave transport
                        // and stress vectors unchanged; transient runs opt out.
                        const float exchangespeed = sea && !climate.transientWeather
                            ? climatephysics::surfaceExchangeWindSpeedMps(phaseexchangewindu[x][y], phaseexchangewindv[x][y],
                                tuning::climate::moistureadvection::minimumSurfaceWind,
                                tuning::climate::moistureadvection::baselineOceanUnresolvedWindRmsMps)
                            : surfacewindspeed[x][y];
                        const float frozencover = sea ? std::clamp(icefraction[x][y], 0.0f, 1.0f)
                            : std::max(landpermanentice[x][y], std::clamp(snowwater[x][y] /
                                tuning::climate::moistureadvection::fullSnowCoverWaterEquivalentMm, 0.0f, 1.0f));
                        const float roughnesslength = sea
                            ? tuning::climate::moistureadvection::oceanRoughnessLengthMetres
                            : tuning::climate::moistureadvection::landRoughnessLengthMetres;
                        const float icecoefficient = climatephysics::neutralDragCoefficient(
                            tuning::climate::moistureadvection::iceRoughnessLengthMetres,
                            tuning::climate::moistureadvection::exchangeReferenceHeightMetres);
                        // Exchange coefficients are area-weighted over open and
                        // frozen surface; trace snow must not switch a whole cell.
                        const float neutralcoefficient =
                            (1.0f - frozencover) * climatephysics::neutralDragCoefficient(
                                roughnesslength,
                                tuning::climate::moistureadvection::exchangeReferenceHeightMetres) +
                            frozencover * icecoefficient;
                        const float referencecoefficient =
                            climatephysics::neutralDragCoefficient(
                                tuning::climate::moistureadvection::oceanRoughnessLengthMetres,
                                tuning::climate::moistureadvection::exchangeReferenceHeightMetres);
                        const float stabilitymultiplier =
                            climatephysics::bulkRichardsonExchangeMultiplier(
                                surfacetemperature,
                                airtemperature,
                                exchangespeed,
                                gravitymultiplier,
                                tuning::climate::moistureadvection::exchangeReferenceHeightMetres,
                                tuning::climate::moistureadvection::minimumStabilityExchangeMultiplier,
                                tuning::climate::moistureadvection::maximumStabilityExchangeMultiplier);
                        const float transfercoefficient =
                            tuning::climate::moistureadvection::surfaceExchangeCoefficient *
                            neutralcoefficient / std::max(0.0001f, referencecoefficient) *
                            stabilitymultiplier * (sea
                                ? tuning::climate::moistureadvection::oceanScalarExchangeFraction : 1.0f);
                        sensibleheat[x][y] = tuning::climate::atmosphere::surfaceAirDensityKgM3 * 1004.0f *
                            transfercoefficient * exchangespeed * (surfacetemperature - airtemperature);
                        const float potentialevaporation =
                            climatephysics::bulkAerodynamicEvaporationMmAtPressure(
                            surfacetemperature,
                            surfacepressure[x][y],
                            exchangespeed,
                            relativehumidity,
                            timestepseconds,
                            transfercoefficient);
                        if (sea && referenceqa)
                            referenceledger[quarter][y * (width + 1) + x][3] += potentialevaporation -
                                climatephysics::bulkAerodynamicEvaporationMmAtPressure(surfacetemperature,
                                    surfacepressure[x][y], exchangespeed, originalhumidity, timestepseconds, transfercoefficient);
                        float seaevaporation = 0.0f;
                        float landevaporation = 0.0f;

                        if (sea)
                        {
                            seaevaporation = potentialevaporation;

                            seaevaporation *= 1.0f - icefraction[x][y] *
                                (1.0f - tuning::climate::moistureadvection::seaIceFactor);

                            if (referenceevap)
                            {
                                const float rate = referenceevapmonthly ? referenceevapmonthly->value(month, x, y) :
                                    climateio::experiment::sample(referenceevap,
                                    calendar.firstSeason, calendar.secondSeason, calendar.interpolation,
                                    x, y, std::numeric_limits<float>::quiet_NaN());
                                if (std::isfinite(rate) && rate != -9999.9f)
                                {
                                    // The reference is already area/time-mean evaporation, including ice.
                                    // Native condensation is unsupported; negative imposed rates are recorded.
                                    auto& ledger = referenceledger[quarter][y * (width + 1) + x];
                                    const float imposed = std::max(0.0f, rate) * timestepseconds / 86400.0f;
                                    ledger[0] += referenceqa ? 0.0 : timestepseconds;
                                    ledger[1] += referenceqa ? 0.0 : rate * timestepseconds;
                                    ledger[2] += rate < 0.0f ? timestepseconds : 0.0;
                                    ledger[3] += imposed - seaevaporation;
                                    seaevaporation = imposed;
                                }
                            }

                            totalseaevaporation[x][y] += seaevaporation;
                        }
                        else
                        {
                            const float melt = climatehydrology::snowMeltAmount(
                                snowwater[x][y],
                                localtemperature,
                                timestepseconds,
                                tuning::climate::moistureadvection::degreeDaySnowMeltMmPerDegreeC);
                            snowwater[x][y] -= melt;
                            const float meltstorageavailable = std::max(
                                0.0f,
                                tuning::climate::moistureadvection::landSoilMoistureCapacity -
                                    soilmoisture[x][y]);
                            const float meltinfiltration = std::min(
                                meltstorageavailable,
                                melt * tuning::climate::moistureadvection::landInfiltrationFraction);
                            soilmoisture[x][y] += meltinfiltration;
                            runoff[x][y] += melt - meltinfiltration;

                            if (snowwater[x][y] > 0.0f)
                            {
                                landevaporation = std::min(
                                    snowwater[x][y],
                                    potentialevaporation *
                                        tuning::climate::moistureadvection::snowSublimationResistance);
                                snowwater[x][y] -= landevaporation;
                            }
                            else if (soilmoisture[x][y] > 0.0f)
                            {
                                const float wateravailability =
                                    climatehydrology::soilMoistureStress(
                                    soilmoisture[x][y],
                                    tuning::climate::moistureadvection::landSoilMoistureCapacity,
                                    tuning::climate::moistureadvection::soilMoistureCriticalFraction,
                                    tuning::climate::moistureadvection::soilMoistureStressExponent);
                                landevaporation = std::min(
                                    soilmoisture[x][y],
                                    potentialevaporation * wateravailability *
                                        tuning::climate::moistureadvection::landSurfaceResistance);
                                soilmoisture[x][y] -= landevaporation;
                            }

                            totallandevaporation[x][y] += landevaporation;
                        }

                        surfaceevaporation[x][y] = seaevaporation + landevaporation;
                        boundarymoisture[x][y] = std::max(0.0f, boundarymoisture[x][y]) +
                            surfaceevaporation[x][y];
                        freemoisture[x][y] = std::max(0.0f, freemoisture[x][y]);
                    }
                }
            });

            advectlayer(
                boundarymoisture,
                phasewindu,
                phasewindv,
                nextboundarymoisture,
                timestepseconds);
            const int uppertimestepcount = iteration %
                tuning::climate::moistureadvection::freeTroposphereTransportCadence + 1;
            const bool updateuppertransport = uppertimestepcount ==
                    tuning::climate::moistureadvection::freeTroposphereTransportCadence ||
                iteration == iterations - 1;

            if (updateuppertransport)
            {
                advectlayer(
                    freemoisture,
                    phaseupperwindu,
                    phaseupperwindv,
                    nextfreemoisture,
                    timestepseconds * static_cast<float>(uppertimestepcount));
                applyeddyexchange(
                    width,
                    height,
                    rowareaweights,
                    nextfreemoisture,
                    zonaleddyexchange,
                    meridionaleddyexchange);
            }
            else
            {
                nextfreemoisture = freemoisture;
            }

            if (month == 0 && iteration == 0 && isworldgendebugrunactive())
                std::cout << "Climate first-step through transport seconds=" << std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - firstStepStart).count() << std::endl;
            damphighwavenumbers(
                nextboundarymoisture,
                highpasslaplacian,
                highpassscratch,
                rowareaweights,
                polartransportfactors,
                width,
                height,
                tuning::climate::moistureadvection::vapourHighWavenumberDamping * timestepseconds / 86400.0f,
                tuning::climate::moistureadvection::highWavenumberDampingPasses,
                true);
            damphighwavenumbers(
                nextfreemoisture,
                highpasslaplacian,
                highpassscratch,
                rowareaweights,
                polartransportfactors,
                width,
                height,
                tuning::climate::moistureadvection::vapourHighWavenumberDamping * timestepseconds / 86400.0f,
                tuning::climate::moistureadvection::highWavenumberDampingPasses,
                true);

            parallelforrows(0, height, [&](int startrow, int endrow)
            {
                for (int y = startrow; y <= endrow; y++)
                {
                    for (int x = 0; x <= width; x++)
                    {

                        boundaryfluxtendency[x][y] =
                            nextboundarymoisture[x][y] - boundarymoisture[x][y];
                        freefluxtendency[x][y] =
                            nextfreemoisture[x][y] - freemoisture[x][y];
                        totalconvergence[x][y] +=
                            (boundaryfluxtendency[x][y] + freefluxtendency[x][y]) /
                            std::max(0.05f, phasesaturationcapacity[x][y]);
                    }
                }
            });

            smoothconvergencefootprint(
                boundaryfluxtendency,
                convectiveconvergence,
                convergencescratch,
                rowareaweights,
                width,
                height,
                tuning::climate::moistureadvection::convectiveConvergenceMixingFraction,
                tuning::climate::moistureadvection::convectiveConvergenceSmoothingPasses);
            damphighwavenumbers(
                convectiveconvergence,
                highpasslaplacian,
                highpassscratch,
                rowareaweights,
                polartransportfactors,
                width,
                height,
                tuning::climate::moistureadvection::convergenceHighWavenumberDamping,
                tuning::climate::moistureadvection::highWavenumberDampingPasses,
                false);

            parallelforrows(0, height, [&](int startrow, int endrow)
            {
                for (int y = startrow; y <= endrow; y++)
                {
                    for (int x = 0; x <= width; x++)
                    {
                        const bool sea = seacell[x][y] != 0;
                        const float localtemperature = phasetemperature[x][y];
                        const float dynamiccooling = std::clamp(
                            std::max(0.0f, phasedynamicvertical[x][y]) *
                                tuning::climate::moistureadvection::dynamicVerticalCooling,
                            0.0f,
                            tuning::climate::moistureadvection::maximumParcelTemperatureAdjustment);
                        const float dynamicwarming = std::clamp(
                            std::max(0.0f, -phasedynamicvertical[x][y]) *
                                tuning::climate::moistureadvection::dynamicSubsidenceWarming,
                            0.0f,
                            tuning::climate::moistureadvection::maximumParcelTemperatureAdjustment);
                        const float nonorographicparceltemperature =
                            localtemperature - dynamiccooling + dynamicwarming -
                            tuning::climate::moistureadvection::freeTroposphereSaturationLapseFraction *
                                tuning::climate::moistureadvection::freeTroposphereEnvironmentalLapseC *
                                polarconvectionfactors[y];
                        const float terraincooling = std::clamp(
                            uplift[x][y] * polarconvectionfactors[y] *
                                tuning::climate::moistureadvection::topographicUpliftCoolingCPerMetre,
                            0.0f,
                            tuning::climate::moistureadvection::maximumParcelTemperatureAdjustment);
                        const float terrainwarming = std::clamp(
                            descent[x][y] *
                                tuning::climate::moistureadvection::topographicDescentWarmingCPerMetre,
                            0.0f,
                            tuning::climate::moistureadvection::maximumParcelTemperatureAdjustment);
                        const float terrainadjustedparceltemperature =
                            nonorographicparceltemperature - terraincooling + terrainwarming;
                        const float nonorographiccapacity =
                            climatephysics::saturationColumnWaterAtPressure(
                                nonorographicparceltemperature,
                                surfacepressure[x][y],
                                gravitymultiplier) *
                            tuning::climate::moistureadvection::freeTroposphereCapacityFraction;
                        const float terrainadjustedcapacity =
                            climatephysics::saturationColumnWaterAtPressure(
                                terrainadjustedparceltemperature,
                                surfacepressure[x][y],
                                gravitymultiplier) *
                            tuning::climate::moistureadvection::freeTroposphereCapacityFraction;
                        const float stratiformrelativehumidity = std::clamp(
                            tuning::climate::moistureadvection::stratiformCriticalRelativeHumidity,
                            0.0f,
                            1.0f);
                        const float stratiformnonorographiccapacity =
                            nonorographiccapacity * stratiformrelativehumidity;
                        const float stratiformterrainadjustedcapacity =
                            terrainadjustedcapacity * stratiformrelativehumidity;
                        const float boundaryrelativehumidity = std::clamp(
                            nextboundarymoisture[x][y] /
                                std::max(0.05f, phaseboundarycapacity[x][y]),
                            0.0f,
                            1.0f);
                        const float freerelativehumidity = std::clamp(
                            nextfreemoisture[x][y] /
                                std::max(0.05f, nonorographiccapacity),
                            0.0f,
                            1.0f);
                        const float cloudfraction = std::max(
                            climatehydrology::diagnosticCloudFraction(
                                boundaryrelativehumidity,
                                tuning::climate::moistureadvection::cloudOnsetRelativeHumidity),
                            climatehydrology::diagnosticCloudFraction(
                                freerelativehumidity,
                                tuning::climate::moistureadvection::cloudOnsetRelativeHumidity));
                        const float dayfraction = timestepseconds /
                            tuning::climate::circulation::secondsPerDay;
                        const float backgroundexchange = -std::expm1(
                            -dayfraction /
                            tuning::climate::moistureadvection::backgroundVerticalExchangeTimeDays);
                        const float normalizedvertical = std::clamp(
                            phasedynamicvertical[x][y] /
                                tuning::climate::circulation::maximumVerticalVelocity,
                            -1.0f,
                            1.0f);
                        const float convectivesupplyfraction = std::clamp(
                            std::max(
                                0.0f,
                                convectiveconvergence[x][y] + surfaceevaporation[x][y] +
                                    std::max(0.0f, freefluxtendency[x][y]) *
                                        tuning::climate::moistureadvection::elevatedMoistureAccessionFraction) /
                                std::max(0.05f, phaseboundarycapacity[x][y]),
                            0.0f,
                            1.0f);
                        const float environmentalfreetemperature =
                            temperature[x][y] -
                            tuning::climate::moistureadvection::freeTroposphereEnvironmentalLapseC *
                                polarconvectionfactors[y] +
                            dynamicwarming;
                        const float liftmetres = tuning::climate::moistureadvection::convectiveParcelLiftMetres *
                            polarconvectionfactors[y];
                        const float liftedparceltemperature = liftmetres > 0.0f ?
                            climatehydrology::detail::liftedParcelTemperatureC(
                                localtemperature, boundaryrelativehumidity, surfacepressure[x][y],
                                surfacepressure[x][y] * std::exp(-9.80665f * std::max(0.05f, gravitymultiplier) *
                                    liftmetres / (287.05f * std::max(150.0f, localtemperature + 273.15f)))) :
                            localtemperature;
                        const float parcelbuoyancy = liftedparceltemperature - environmentalfreetemperature;
                        const float verticalwindshear = std::sqrt(
                            std::pow(phaseupperwindu[x][y] - phasewindu[x][y], 2.0f) +
                            std::pow(phaseupperwindv[x][y] - phasewindv[x][y], 2.0f));
                        const float shallowexchange =
                            climatehydrology::shallowConvectionExchangeFraction(
                                boundaryrelativehumidity,
                                freerelativehumidity,
                                parcelbuoyancy,
                                verticalwindshear,
                                timestepseconds,
                                tuning::climate::moistureadvection::shallowConvectionMixingTimeDays,
                                tuning::climate::moistureadvection::shallowConvectionHumidityOnset,
                                tuning::climate::moistureadvection::shallowConvectionFullHumidity,
                                tuning::climate::moistureadvection::shallowConvectionFullShearMps,
                                tuning::climate::moistureadvection::maximumShallowExchangeFraction);
                        const float dryexchange =
                            climatehydrology::dryConvectionExchangeFraction(
                                parcelbuoyancy,
                                timestepseconds,
                                tuning::climate::moistureadvection::dryConvectionMixingTimeDays,
                                tuning::climate::moistureadvection::dryConvectionActivationBuoyancyC,
                                tuning::climate::moistureadvection::dryConvectionFullStrengthBuoyancyC,
                                tuning::climate::moistureadvection::maximumDryExchangeFraction);
                        const float upwardfraction = std::clamp(
                            backgroundexchange +
                                climatehydrology::detail::dailyFractionForTimeStep(
                                    std::max(0.0f, normalizedvertical) *
                                        tuning::climate::moistureadvection::ascentVerticalExchangeFraction,
                                    timestepseconds) +
                                cloudfraction * convectivesupplyfraction *
                                    tuning::climate::moistureadvection::convectiveVerticalExchangeFraction +
                                shallowexchange + dryexchange,
                            0.0f,
                            tuning::climate::moistureadvection::maximumVerticalExchangeFraction);
                        const float downwardfraction = std::clamp(
                            backgroundexchange +
                                climatehydrology::detail::dailyFractionForTimeStep(
                                    std::max(0.0f, -normalizedvertical) *
                                        tuning::climate::moistureadvection::subsidenceVerticalExchangeFraction,
                                    timestepseconds) +
                                0.5f * shallowexchange + dryexchange,
                            0.0f,
                            tuning::climate::moistureadvection::maximumVerticalExchangeFraction);
                        const climatehydrology::MoistureLayerExchange layerexchange =
                            climatehydrology::exchangeMoistureLayers(
                                nextboundarymoisture[x][y],
                                nextfreemoisture[x][y],
                                upwardfraction,
                                downwardfraction);
                        const float availablemoisture =
                            layerexchange.boundaryLayerMm + layerexchange.freeTroposphereMm;
                        const float condensablewater = std::max(
                            0.0f,
                            layerexchange.freeTroposphereMm -
                                stratiformterrainadjustedcapacity);
                        const climatehydrology::PrecipitationPartition partition =
                            climatehydrology::partitionTwoLayerPrecipitation(
                            layerexchange.boundaryLayerMm,
                            layerexchange.freeTroposphereMm,
                            phaseboundarycapacity[x][y],
                            stratiformnonorographiccapacity,
                            stratiformterrainadjustedcapacity,
                            convectiveconvergence[x][y],
                            surfaceevaporation[x][y],
                            liftedparceltemperature,
                            environmentalfreetemperature,
                            updateuppertransport
                                ? timestepseconds * static_cast<float>(uppertimestepcount)
                                : 0.0f,
                            condensationconversiontimeseconds,
                            tuning::climate::moistureadvection::kuoCriticalRelativeHumidity,
                            tuning::climate::moistureadvection::convectiveConversionEfficiency,
                            tuning::climate::moistureadvection::convectiveActivationBuoyancyC,
                            tuning::climate::moistureadvection::convectiveFullStrengthBuoyancyC,
                            tuning::climate::moistureadvection::moistAdjustmentIterations,
                            tuning::climate::moistureadvection::latentHeatingCPerMillimetre,
                            tuning::climate::moistureadvection::saturationCapacityTemperatureSensitivityPerC,
                            freefluxtendency[x][y],
                            tuning::climate::moistureadvection::elevatedMoistureAccessionFraction,
                            tuning::climate::moistureadvection::kuoHumidityExponent);
                        const float condensate = partition.totalMm();
                        const float boundaryaftercondensation = std::max(
                            0.0f,
                            layerexchange.boundaryLayerMm - partition.convectiveMm);
                        const float freeaftercondensation = std::max(
                            0.0f,
                            layerexchange.freeTroposphereMm -
                                partition.stratiformMm - partition.orographicMm);
                        const float falloutrelativehumidity = std::clamp(
                            boundaryaftercondensation /
                                std::max(0.05f, phaseboundarycapacity[x][y]),
                            0.0f,
                            1.0f);
                        const climatehydrology::FallingPrecipitation falling =
                            climatehydrology::processFallingPrecipitation(
                                condensate,
                                localtemperature,
                                falloutrelativehumidity,
                                tuning::climate::moistureadvection::maximumFallingPrecipitationReevaporationFraction,
                                std::max(
                                    0.0f,
                                    phaseboundarycapacity[x][y] - boundaryaftercondensation),
                                tuning::climate::moistureadvection::allSnowTemperatureC,
                                tuning::climate::moistureadvection::allRainTemperatureC);
                        const float precipitation = falling.surfaceTotalMm();
                        const int cell = y * (width + 1) + x;
                        climateatmosphere::ColumnHeatingInput heatinput;
                        const double falloutfraction = condensate > 0.0f ? precipitation / condensate : 0.0;
                        exchanges.surfacePrecipitationMm[0][cell] += partition.stratiformMm * falloutfraction;
                        exchanges.surfacePrecipitationMm[1][cell] += partition.orographicMm * falloutfraction;
                        exchanges.surfacePrecipitationMm[2][cell] += partition.convectiveMm * falloutfraction;
                        exchanges.relativeHumiditySeconds[0][cell] += boundaryrelativehumidity * timestepseconds;
                        exchanges.relativeHumiditySeconds[1][cell] += freerelativehumidity * timestepseconds;
                        exchanges.upliftMetreSeconds[cell] += uplift[x][y] * timestepseconds;
                        exchanges.descentMetreSeconds[cell] += descent[x][y] * timestepseconds;
                        exchanges.terrainCoolingCSeconds[cell] += terraincooling * timestepseconds;
                        exchanges.terrainWarmingCSeconds[cell] += terrainwarming * timestepseconds;
                        heatinput.incomingSolarWm2 = climateenergy::dailyMeanInsolationWm2(climatelatitudes[y], solardeclination, solardistance);
                        const float frozencover = seacell[x][y] ? std::clamp(icefraction[x][y], 0.0f, 1.0f)
                            : std::clamp(snowwater[x][y] /
                                tuning::climate::moistureadvection::fullSnowCoverWaterEquivalentMm, 0.0f, 1.0f);
                        const double openalbedo = seacell[x][y] ? tuning::climate::energybalance::oceanAlbedo
                            : tuning::climate::energybalance::landAlbedo;
                        const double frozenalbedo = seacell[x][y] ? tuning::climate::energybalance::seaIceAlbedo
                            : tuning::climate::energybalance::snowAlbedo;
                        heatinput.surfaceAlbedo = openalbedo + frozencover * (frozenalbedo - openalbedo);
                        if (!seacell[x][y] && landpermanentice[x][y] > 0.0f)
                            heatinput.surfaceAlbedo = climateatmosphere::detail::frozenLandAlbedo(
                                frozencover, landpermanentice[x][y], static_cast<float>(openalbedo),
                                static_cast<float>(frozenalbedo), static_cast<float>(tuning::climate::energybalance::permanentIceAlbedo));
                        heatinput.surfaceTemperatureK = phasesurfacetemperature[x][y] + 273.15;
                        heatinput.airTemperatureK = {localtemperature + 273.15,
                            environmentalfreetemperature + 273.15};
                        heatinput.longwaveOpticalDepth = {tuning::climate::circulation::lowerLongwaveOpticalDepth,
                            tuning::climate::circulation::upperLongwaveOpticalDepth};
                        heatinput.shortwaveOpticalDepth = {tuning::climate::circulation::lowerShortwaveOpticalDepth,
                            tuning::climate::circulation::upperShortwaveOpticalDepth};
                        heatinput.sensibleHeatingWm2 = sensibleheat[x][y];
                        heatinput.condensationMm = {partition.convectiveMm, partition.stratiformMm + partition.orographicMm};
                        heatinput.reevaporationMm = falling.reevaporatedMm;
                        heatinput.surfaceEvaporationMm = surfaceevaporation[x][y];
                        heatinput.accumulationSeconds = timestepseconds;
                        auto heat = climateatmosphere::diagnoseColumnHeating(heatinput);
                        if (seacell[x][y] && referenceradiation)
                        {
                            const float imposed = referenceradiation->value(month, x, y);
                            if (std::isfinite(imposed) && imposed != -9999.9f)
                            {
                                const double external = imposed - heat.surfaceRadiativeWm2;
                                heat.surfaceRadiativeWm2 = imposed;
                                heat.surfaceNetHeatingWm2 += external;
                                referenceledger[quarter][cell][4] += external * timestepseconds;
                                // Native column residual remains the residual after subtracting this
                                // explicitly recorded external surface-energy supply. Atmospheric
                                // radiation is deliberately left free in this surface-only experiment.
                            }
                        }
                        for (int layer = 0; layer < 2; ++layer)
                        {
                            exchanges.radiativeHeatingWm2[layer][cell] += static_cast<float>(heat.radiativeWm2[layer] * timestepseconds);
                            exchanges.latentHeatingWm2[layer][cell] += static_cast<float>(heat.latentWm2[layer] * timestepseconds);
                        }
                        exchanges.sensibleHeatingWm2[cell] += sensibleheat[x][y] * timestepseconds;
                        exchanges.surfaceTemperatureC[cell] += phasesurfacetemperature[x][y] * timestepseconds;
                        exchanges.landSnowCoverFraction[cell] += (seacell[x][y] ? 0.0f : frozencover) * timestepseconds;
                        exchanges.permanentLandIceFraction[cell] += (seacell[x][y] ? 0.0f : landpermanentice[x][y]) * timestepseconds;
                        exchanges.surfaceNetHeatingWm2[cell] += static_cast<float>(heat.surfaceNetHeatingWm2 * timestepseconds);
                        exchanges.surfaceRadiativeHeatingWm2[cell] += static_cast<float>(heat.surfaceRadiativeWm2 * timestepseconds);
                        exchanges.surfaceLatentHeatLossWm2[cell] += 2.5e6f * surfaceevaporation[x][y];
                        exchanges.columnEnergyResidualWm2[cell] += static_cast<float>(heat.closureResidualWm2 * timestepseconds);
                        if (!seacell[x][y])
                        {
                            const auto thermal = climateenergy::stepLandThermalState(
                                {landSurface[x][y], landDeep[x][y]}, heat.surfaceNetHeatingWm2,
                                tuning::climate::energybalance::landHeatCapacityJm2K,
                                tuning::climate::energybalance::landDeepLayerHeatCapacityJm2K,
                                tuning::climate::energybalance::landDeepLayerCouplingWm2K, timestepseconds);
                            landSurface[x][y] = thermal.surfaceC;
                            landDeep[x][y] = thermal.deepC;
                            landInputEnergy[x][y] += heat.surfaceNetHeatingWm2 * timestepseconds;
                        }
                        exchanges.columnWaterMm[cell] += (boundaryaftercondensation + falling.reevaporatedMm + freeaftercondensation) * timestepseconds;
                        exchanges.ascentHpaPerDay[cell] += phasedynamicvertical[x][y] * timestepseconds;
                        boundarymoisture[x][y] =
                            boundaryaftercondensation + falling.reevaporatedMm;
                        freemoisture[x][y] = freeaftercondensation;
                        const float heatingdecay = std::exp(
                            -timestepseconds /
                            (tuning::climate::moistureadvection::persistentHeatingDecayTimeDays *
                                tuning::climate::circulation::secondsPerDay));
                        persistentheating[x][y] = std::clamp(
                            persistentheating[x][y] * heatingdecay +
                                precipitation *
                                    tuning::climate::moistureadvection::persistentHeatingRetentionPerMillimetreC,
                            0.0f,
                            tuning::climate::moistureadvection::maximumPersistentHeatingC);
                        const float clouddecay = std::exp(
                            -timestepseconds /
                            (tuning::climate::moistureadvection::cloudMemoryTimeDays *
                                tuning::climate::circulation::secondsPerDay));
                        cloudmemory[x][y] = std::clamp(
                            cloudmemory[x][y] * clouddecay +
                                cloudfraction * (1.0f - clouddecay),
                            0.0f,
                            1.0f);
                        const double areaweight = rowareaweights[y];
                        climatephysics::CondensationActivityDiagnostics& activity =
                            rowcondensation[y];
                        climatephysics::PrecipitationProcessDiagnostics& processes =
                            rowprocesses[y];
                        activity.cellStepAreaWeight += areaweight;
                        activity.atmosphericWaterAreaWeighted +=
                            areaweight * static_cast<double>(availablemoisture);

                        if (precipitation > 0.0f)
                        {
                            activity.activeCellStepAreaWeight += areaweight;
                            activity.activeAtmosphericWaterAreaWeighted +=
                                areaweight * static_cast<double>(availablemoisture);
                            activity.excessWaterAreaWeighted +=
                                areaweight * static_cast<double>(
                                    std::max(condensablewater, precipitation));
                        }

                        activity.precipitationAreaWeighted +=
                            areaweight * static_cast<double>(precipitation);
                        processes.stratiformPrecipitation +=
                            areaweight * static_cast<double>(
                                condensate > 0.0f
                                    ? partition.stratiformMm * precipitation / condensate
                                    : 0.0f);
                        processes.orographicPrecipitation +=
                            areaweight * static_cast<double>(
                                condensate > 0.0f
                                    ? partition.orographicMm * precipitation / condensate
                                    : 0.0f);
                        processes.convectivePrecipitation +=
                            areaweight * static_cast<double>(
                                condensate > 0.0f
                                    ? partition.convectiveMm * precipitation / condensate
                                    : 0.0f);
                        processes.reevaporatedPrecipitation +=
                            areaweight * static_cast<double>(falling.reevaporatedMm);
                        processes.snowfall +=
                            areaweight * static_cast<double>(falling.snowMm);
                        processes.upwardMoistureTransfer +=
                            areaweight * static_cast<double>(layerexchange.upwardTransferMm);
                        processes.downwardMoistureTransfer +=
                            areaweight * static_cast<double>(layerexchange.downwardTransferMm);
                        processes.cloudFractionAreaWeighted += areaweight * cloudfraction;

                        const float combinedfluxtendency =
                            boundaryfluxtendency[x][y] + freefluxtendency[x][y];

                        if (combinedfluxtendency >= 0.0f)
                        {
                            processes.positiveMoistureFluxConvergence +=
                                areaweight * static_cast<double>(combinedfluxtendency);
                        }
                        else
                        {
                            processes.negativeMoistureFluxConvergence +=
                                areaweight * static_cast<double>(combinedfluxtendency);
                        }
                        totalrain[x][y] += precipitation;

                        if (sea == false)
                        {
                            const climatehydrology::SnowAccumulation snowaccumulation =
                                climatehydrology::accumulateSnowfall(
                                    snowwater[x][y],
                                    falling.snowMm,
                                    tuning::climate::moistureadvection::maximumSnowStorageMm);
                            snowwater[x][y] = snowaccumulation.storageMm;
                            const float storageavailable = std::max(
                                0.0f,
                                tuning::climate::moistureadvection::landSoilMoistureCapacity - soilmoisture[x][y]);
                            const float infiltration = std::min(
                                storageavailable,
                                falling.rainMm *
                                    tuning::climate::moistureadvection::landInfiltrationFraction);

                            soilmoisture[x][y] += infiltration;
                            runoff[x][y] += falling.rainMm - infiltration +
                                snowaccumulation.overflowMm;
                        }
                    }
                }
            });
        }

        climatephysics::CondensationActivityDiagnostics condensationactivity;
        climatephysics::PrecipitationProcessDiagnostics processdiagnostics;

        for (const climatephysics::CondensationActivityDiagnostics& row : rowcondensation)
        {
            condensationactivity.cellStepAreaWeight += row.cellStepAreaWeight;
            condensationactivity.activeCellStepAreaWeight += row.activeCellStepAreaWeight;
            condensationactivity.atmosphericWaterAreaWeighted +=
                row.atmosphericWaterAreaWeighted;
            condensationactivity.activeAtmosphericWaterAreaWeighted +=
                row.activeAtmosphericWaterAreaWeighted;
            condensationactivity.excessWaterAreaWeighted += row.excessWaterAreaWeighted;
            condensationactivity.precipitationAreaWeighted +=
                row.precipitationAreaWeighted;
        }

        for (const climatephysics::PrecipitationProcessDiagnostics& row : rowprocesses)
        {
            processdiagnostics.stratiformPrecipitation +=
                row.stratiformPrecipitation;
            processdiagnostics.orographicPrecipitation +=
                row.orographicPrecipitation;
            processdiagnostics.convectivePrecipitation +=
                row.convectivePrecipitation;
            processdiagnostics.reevaporatedPrecipitation +=
                row.reevaporatedPrecipitation;
            processdiagnostics.snowfall += row.snowfall;
            processdiagnostics.upwardMoistureTransfer +=
                row.upwardMoistureTransfer;
            processdiagnostics.downwardMoistureTransfer +=
                row.downwardMoistureTransfer;
            processdiagnostics.cloudFractionAreaWeighted +=
                row.cloudFractionAreaWeighted;
            processdiagnostics.positiveMoistureFluxConvergence +=
                row.positiveMoistureFluxConvergence;
            processdiagnostics.negativeMoistureFluxConvergence +=
                row.negativeMoistureFluxConvergence;
        }

        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int x = 0; x <= width; x++)
                {
                    const float rawrain = totalrain[x][y] *
                        tuning::climate::moistureadvection::rainfallScale;
                    monthlyrain[month][y * (width + 1) + x] = rawrain;
                    rawannualrain[x][y] += rawrain;
                    quarterrain[x][y] += rawrain;
                    quarterconvergence[x][y] +=
                        totalconvergence[x][y] / static_cast<float>(iterations);


                    if (month % 3 == 2)
                    {
                        const float physicalrain = std::max(
                            0.0f,
                            quarterrain[x][y] / 3.0f);
                        storedseasonalrain[quarter][x][y] = physicalrain;
                        storedseasonalmoisture[quarter][x][y] =
                            seasonalprocess[quarter].columnWaterMm[y * (width + 1) + x] /
                            static_cast<float>(seasonalprocess[quarter].durationSeconds);
                        storedseasonalconvergence[quarter][x][y] = quarterconvergence[x][y] / 3.0f;
                    }
                }
            }
        });

        for (int y = 0; y <= height; y++)
        {
            const double areaweight = rowareaweights[y];

            for (int x = 0; x <= width; x++)
            {
                if (seacell[x][y] != 0)
                {
                    budget.oceanEvaporation += totalseaevaporation[x][y];
                    budget.oceanPrecipitation += totalrain[x][y];
                    areaweightedbudget.oceanEvaporation += areaweight * totalseaevaporation[x][y];
                    areaweightedbudget.oceanPrecipitation += areaweight * totalrain[x][y];
                }
                else
                {
                    budget.landEvaporation += totallandevaporation[x][y];
                    budget.landPrecipitation += totalrain[x][y];
                    budget.runoff += runoff[x][y];
                    budget.soilStorage += soilmoisture[x][y];
                    budget.snowStorage += snowwater[x][y];
                    areaweightedbudget.landEvaporation += areaweight * totallandevaporation[x][y];
                    areaweightedbudget.landPrecipitation += areaweight * totalrain[x][y];
                    areaweightedbudget.runoff += areaweight * runoff[x][y];
                    areaweightedbudget.soilStorage += areaweight * soilmoisture[x][y];
                    areaweightedbudget.snowStorage += areaweight * snowwater[x][y];
                }

                const float atmosphericwater =
                    boundarymoisture[x][y] + freemoisture[x][y];
                budget.atmosphericStorage += atmosphericwater;
                areaweightedbudget.atmosphericStorage += areaweight * atmosphericwater;
            }
        }

        continuedweather = weatherstate;
        return monthlysolverresult{
            budget,
            areaweightedbudget,
            condensationactivity,
            processdiagnostics,
            climateweather::serializeState(weatherstate)
        };
    };

    climatephysics::HydrologySpinupDiagnostics diagnostics;

    for (int cycle = 0; cycle < tuning::climate::moistureadvection::yearsPerCouplingUpdate; cycle++)
    {
        transportSeconds = weatherSeconds = 0.0;
        transportCalls = transportCellSubsteps = 0;
        maximumTransportSubsteps = 0;
        // Replay the same evolving year while hydrological storage spins up.
        // Months within a year share one dynamical state, not independent noise.
        continuedweather = {};
        for (auto& samples : seasonalweather) samples.clear();
        for (auto& weights : seasonalweights) weights.clear();
        for (auto& fields : seasonalprocess)
        {
            fields = {};
            fields.columns = width + 1; fields.rows = height + 1;
            const int cells = fields.columns * fields.rows;
            for (int layer = 0; layer < 2; ++layer)
            {
                fields.eastIntegratedFlux[layer].assign(cells, 0.0);
                fields.southIntegratedFlux[layer].assign(cells, 0.0);
                fields.radiativeHeatingWm2[layer].assign(cells, 0.0f);
                fields.latentHeatingWm2[layer].assign(cells, 0.0f);
            }
            fields.sensibleHeatingWm2.assign(cells, 0.0f);
            fields.surfaceNetHeatingWm2.assign(cells, 0.0f);
            fields.surfaceRadiativeHeatingWm2.assign(cells, 0.0f);
            fields.surfaceLatentHeatLossWm2.assign(cells, 0.0f);
            fields.surfaceTemperatureC.assign(cells, 0.0f);
            fields.landSnowCoverFraction.assign(cells, 0.0f);
            fields.permanentLandIceFraction.assign(cells, 0.0f);
            fields.columnWaterMm.assign(cells, 0.0f);
            fields.ascentHpaPerDay.assign(cells, 0.0f);
            fields.columnEnergyResidualWm2.assign(cells, 0.0f);
            for (auto& values : fields.surfacePrecipitationMm) values.assign(cells, 0.0);
            for (auto& values : fields.relativeHumiditySeconds) values.assign(cells, 0.0);
            fields.upliftMetreSeconds.assign(cells, 0.0);
            fields.descentMetreSeconds.assign(cells, 0.0);
            fields.terrainCoolingCSeconds.assign(cells, 0.0);
            fields.terrainWarmingCSeconds.assign(cells, 0.0);
        }
        const auto cyclestart = std::chrono::steady_clock::now();
        const hydrologyfloatgrid initialboundarymoisture = boundarymoisture;
        const hydrologyfloatgrid initialfreemoisture = freemoisture;
        const hydrologyfloatgrid initialsoilmoisture = soilmoisture;
        const hydrologyfloatgrid initialsnowwater = snowwater;
        const auto initialLandSurface = landSurface;
        const auto initialLandDeep = landDeep;
        landInputEnergy.fill(0.0);

        rawannualrain.fill(0.0f);

        finalbudgets = {};
        finalareaweightedbudgets = {};
        finalcondensationactivity = {};
        finalprocessdiagnostics = {};

        auto appendbudget = [](
            climatephysics::WaterBudget& quarterbudget,
            const climatephysics::WaterBudget& monthbudget,
            bool firstmonth)
        {
            if (firstmonth)
            {
                quarterbudget.initialAtmosphericStorage =
                    monthbudget.initialAtmosphericStorage;
                quarterbudget.initialSoilStorage = monthbudget.initialSoilStorage;
                quarterbudget.initialSnowStorage = monthbudget.initialSnowStorage;
            }

            quarterbudget.oceanEvaporation += monthbudget.oceanEvaporation;
            quarterbudget.landEvaporation += monthbudget.landEvaporation;
            quarterbudget.oceanPrecipitation += monthbudget.oceanPrecipitation;
            quarterbudget.landPrecipitation += monthbudget.landPrecipitation;
            quarterbudget.runoff += monthbudget.runoff;
            quarterbudget.atmosphericStorage = monthbudget.atmosphericStorage;
            quarterbudget.soilStorage = monthbudget.soilStorage;
            quarterbudget.snowStorage = monthbudget.snowStorage;
        };

        for (int month = 0; month < climatehydrology::monthCount; month++)
        {
            const monthlysolverresult result = runsolver(month);
            if (month == climatehydrology::monthCount - 1)
                finalweatherstate = result.weatherstate;
            const int quarter = month / 3;
            appendbudget(finalbudgets[quarter], result.budget, month % 3 == 0);
            appendbudget(
                finalareaweightedbudgets[quarter],
                result.areaweightedbudget,
                month % 3 == 0);

            auto& activity = finalcondensationactivity[quarter];
            activity.cellStepAreaWeight += result.condensation.cellStepAreaWeight;
            activity.activeCellStepAreaWeight +=
                result.condensation.activeCellStepAreaWeight;
            activity.atmosphericWaterAreaWeighted +=
                result.condensation.atmosphericWaterAreaWeighted;
            activity.activeAtmosphericWaterAreaWeighted +=
                result.condensation.activeAtmosphericWaterAreaWeighted;
            activity.excessWaterAreaWeighted +=
                result.condensation.excessWaterAreaWeighted;
            activity.precipitationAreaWeighted +=
                result.condensation.precipitationAreaWeighted;

            auto& processes = finalprocessdiagnostics[quarter];
            processes.stratiformPrecipitation +=
                result.processes.stratiformPrecipitation;
            processes.orographicPrecipitation +=
                result.processes.orographicPrecipitation;
            processes.convectivePrecipitation +=
                result.processes.convectivePrecipitation;
            processes.reevaporatedPrecipitation +=
                result.processes.reevaporatedPrecipitation;
            processes.snowfall += result.processes.snowfall;
            processes.upwardMoistureTransfer +=
                result.processes.upwardMoistureTransfer;
            processes.downwardMoistureTransfer +=
                result.processes.downwardMoistureTransfer;
            processes.cloudFractionAreaWeighted +=
                result.processes.cloudFractionAreaWeighted;
            processes.positiveMoistureFluxConvergence +=
                result.processes.positiveMoistureFluxConvergence;
            processes.negativeMoistureFluxConvergence +=
                result.processes.negativeMoistureFluxConvergence;
        }

        double atmosphericabsolutechange = 0.0;
        double atmosphericreferencestorage = 0.0;
        double soilabsolutechange = 0.0;
        double soilreferencestorage = 0.0;
        double snowabsolutechange = 0.0;
        double snowreferencestorage = 0.0;
        double snowcoverchangedarea = 0.0;
        double landarea = 0.0, landChange = 0.0, landResidual = 0.0;
        double atmosphericstorage = 0.0;
        double soilstorage = 0.0;
        double snowstorage = 0.0;

        for (int y = 0; y <= height; y++)
        {
            const double areaweight = rowareaweights[y];

            for (int x = 0; x <= width; x++)
            {
                const double atmosphericwater =
                    boundarymoisture[x][y] + freemoisture[x][y];
                const double initialatmosphericwater =
                    initialboundarymoisture[x][y] + initialfreemoisture[x][y];
                atmosphericabsolutechange += areaweight * std::abs(
                    atmosphericwater - initialatmosphericwater);
                atmosphericreferencestorage += areaweight * 0.5 *
                    (atmosphericwater + initialatmosphericwater);
                atmosphericstorage += areaweight * atmosphericwater;

                if (seacell[x][y] == 0)
                {
                    landarea += areaweight;
                    const double ds = landSurface[x][y] - initialLandSurface[x][y];
                    const double dd = landDeep[x][y] - initialLandDeep[x][y];
                    landChange += areaweight * std::max(ds * ds, dd * dd);
                    landResidual += areaweight * std::abs(
                        tuning::climate::energybalance::landHeatCapacityJm2K * ds +
                        tuning::climate::energybalance::landDeepLayerHeatCapacityJm2K * dd - landInputEnergy[x][y]);
                    soilabsolutechange += areaweight * std::abs(static_cast<double>(
                        soilmoisture[x][y] - initialsoilmoisture[x][y]));
                    soilreferencestorage += areaweight * 0.5 * static_cast<double>(
                        soilmoisture[x][y] + initialsoilmoisture[x][y]);
                    soilstorage += areaweight * soilmoisture[x][y];
                    snowabsolutechange += areaweight * std::abs(static_cast<double>(
                        snowwater[x][y] - initialsnowwater[x][y]));
                    snowreferencestorage += areaweight * 0.5 * static_cast<double>(
                        snowwater[x][y] + initialsnowwater[x][y]);
                    snowstorage += areaweight * snowwater[x][y];
                    const bool hassnow = snowwater[x][y] >=
                        tuning::climate::moistureadvection::snowCoverConvergenceThresholdMm;
                    const bool initiallyhadsnow = initialsnowwater[x][y] >=
                        tuning::climate::moistureadvection::snowCoverConvergenceThresholdMm;

                    if (hassnow != initiallyhadsnow)
                        snowcoverchangedarea += areaweight;
                }
            }
        }

        diagnostics.cyclesCompleted = ++climate.hydrologyYears;
        diagnostics.relativeAtmosphericStorageChange = atmosphericabsolutechange /
            std::max(1.0, atmosphericreferencestorage);
        diagnostics.relativeSoilStorageChange = soilabsolutechange /
            std::max(1.0, soilreferencestorage);
        diagnostics.relativeSnowStorageChange = snowabsolutechange /
            std::max(1.0, snowreferencestorage);
        diagnostics.relativeSnowCoverChange = snowcoverchangedarea /
            std::max(1.0, landarea);
        diagnostics.relativeStorageChange = std::max(
            diagnostics.relativeAtmosphericStorageChange,
            std::max(
                diagnostics.relativeSoilStorageChange,
                std::max(diagnostics.relativeSnowCoverChange, diagnostics.relativeSnowStorageChange)));
        diagnostics.atmosphericStorage = atmosphericstorage;
        diagnostics.soilStorage = soilstorage;
        diagnostics.snowStorage = snowstorage;
        climate.landAnnualDriftK = std::sqrt(landChange / std::max(1.0e-30, landarea));
        climate.landEnergyResidualJm2 = landResidual / std::max(1.0e-30, landarea);
        diagnostics.converged = climate.landAnnualDriftK <= 0.05 && diagnostics.cyclesCompleted >=
                tuning::climate::moistureadvection::minimumSpinupCycles &&
            diagnostics.relativeStorageChange <=
                tuning::climate::moistureadvection::spinupRelativeStorageTolerance;
        const double cycleelapsedseconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - cyclestart).count();

        std::cout
            << "Climate hydrology spinup cycle=" << diagnostics.cyclesCompleted
            << " relative_storage_change=" << diagnostics.relativeStorageChange
            << " relative_atmospheric_change="
            << diagnostics.relativeAtmosphericStorageChange
            << " relative_soil_change=" << diagnostics.relativeSoilStorageChange
            << " relative_snow_change=" << diagnostics.relativeSnowStorageChange
            << " relative_snow_cover_change=" << diagnostics.relativeSnowCoverChange
            << " atmospheric_storage=" << diagnostics.atmosphericStorage
            << " soil_storage=" << diagnostics.soilStorage
            << " snow_storage=" << diagnostics.snowStorage
            << " land_drift_k=" << climate.landAnnualDriftK
            << " land_energy_residual_jm2=" << climate.landEnergyResidualJm2
            << " elapsed_seconds=" << cycleelapsedseconds
            << " transport_seconds=" << transportSeconds
            << " weather_seconds=" << weatherSeconds
            << " transport_calls=" << transportCalls
            << " transport_cell_substeps=" << transportCellSubsteps
            << " transport_max_substeps=" << maximumTransportSubsteps
            << " converged=" << (diagnostics.converged ? 1 : 0)
            << std::endl;

        if (diagnostics.converged)
            break;
    }

    const auto outputremap = climategrid::makeSmoothRemap(
        width + 1, height + 1, outputwidth + 1, outputheight + 1);
    auto& monthly = climate.monthly;
    monthly.columns = outputwidth + 1; monthly.rows = outputheight + 1;
    for (int month = 0; month < 12; ++month)
        monthly.rainfallMm[month] = climategrid::remapField(outputremap, monthlyrain[month], width + 1);
    for (auto& field : monthly.referenceRainfallMm) field.resize((outputwidth + 1) * (outputheight + 1));

    hydrologyfloatgrid outputrawannualrain(
        outputwidth + 1,
        outputheight + 1,
        0.0f);
    hydrologyfloatgrid outputrawmaximummonthlyrain(
        outputwidth + 1,
        outputheight + 1,
        0.0f);
    hydrologybytegrid outputrawmaximumrainmonth(
        outputwidth + 1,
        outputheight + 1,
        static_cast<signed char>(-1));
    vector<float> outputrowareaweights(outputheight + 1, 0.0f);

    for (int y = 0; y <= outputheight; y++)
        outputrowareaweights[y] = static_cast<float>(rowareaweight(y, outputheight));

    parallelforrows(0, outputheight, [&](int startrow, int endrow)
    {
        for (int y = startrow; y <= endrow; y++)
        {
            for (int x = 0; x <= outputwidth; x++)
            {
                outputrawannualrain[x][y] = samplehydrologygrid(
                    rawannualrain, outputremap, x, y);
                for (int month = 0; month < 12; ++month)
                {
                    const float rain = monthly.rainfallMm[month][y * (outputwidth + 1) + x];
                    if (rain > outputrawmaximummonthlyrain[x][y])
                    {
                        outputrawmaximummonthlyrain[x][y] = rain;
                        outputrawmaximumrainmonth[x][y] = static_cast<signed char>(month);
                    }
                }

                for (int season = 0; season < CLIMATESEASONCOUNT; season++)
                {
                    const float seasonalrain = std::max(
                        0.0f,
                        samplehydrologygrid(
                            storedseasonalrain[season], outputremap, x, y));
                    const float seasonalmoisture = samplehydrologygrid(
                        storedseasonalmoisture[season], outputremap, x, y);
                    const float seasonalconvergence = samplehydrologygrid(
                        storedseasonalconvergence[season], outputremap, x, y);
                    world.setseasonalrainfloat(season, x, y, seasonalrain);
                    monthly.referenceRainfallMm[season][y * (outputwidth + 1) + x] = seasonalrain;
                    world.setseasonalmoisturefloat(season, x, y, seasonalmoisture);
                    world.setseasonalconvergencefloat(season, x, y, seasonalconvergence);

                    if (season == seasonjanuary)
                        world.setjanrainfloat(x, y, seasonalrain);
                    else if (season == seasonjuly)
                        world.setjulrainfloat(x, y, seasonalrain);
                }
            }
        }
    });

    struct distributionaccumulator
    {
        climatephysics::PrecipitationDistributionScope metrics;
        long long rawzeros = 0;
        long long storedzeros = 0;
        long long belowone = 0;
        double areaweightedrawzeros = 0.0;
        double areaweightedstoredzeros = 0.0;
        double areaweightedbelowone = 0.0;
        double precipitationtotal = 0.0;
        double areaweightedprecipitationtotal = 0.0;
        vector<float> precipitation;
    };

    distributionaccumulator landdistribution;
    distributionaccumulator oceandistribution;

    for (int y = 0; y <= outputheight; y++)
    {
        const double areaweight = outputrowareaweights[y];

        for (int x = 0; x <= outputwidth; x++)
        {
            distributionaccumulator& distribution = world.sea(x, y) == 1
                ? oceandistribution
                : landdistribution;
            const float meanmonthlyrain = outputrawannualrain[x][y] /
                static_cast<float>(climatehydrology::monthCount);
            float storedseasonalrain = 0.0f;

            for (int season = 0; season < CLIMATESEASONCOUNT; season++)
                storedseasonalrain += world.seasonalrainfloat(season, x, y);

            distribution.metrics.cells++;
            distribution.metrics.areaWeight += areaweight;
            distribution.precipitationtotal += meanmonthlyrain;
            distribution.areaweightedprecipitationtotal += areaweight * meanmonthlyrain;
            if (outputrawmaximummonthlyrain[x][y] >
                distribution.metrics.maximumMonthlyPrecipitationMm)
            {
                const int maximummonth = static_cast<int>(outputrawmaximumrainmonth[x][y]);
                distribution.metrics.maximumMonthlyPrecipitationMm =
                    static_cast<double>(outputrawmaximummonthlyrain[x][y]);
                distribution.metrics.maximumMonthlyPrecipitationX = x;
                distribution.metrics.maximumMonthlyPrecipitationY = y;
                distribution.metrics.maximumMonthlyPrecipitationMonth = maximummonth;
                distribution.metrics.roundedPrecipitationAtMaximumMm =
                    maximummonth >= 0
                        ? static_cast<int>(std::round(outputrawmaximummonthlyrain[x][y]))
                        : 0;
            }
            distribution.precipitation.push_back(meanmonthlyrain);

            if (outputrawannualrain[x][y] <= 1.0e-6f)
            {
                distribution.rawzeros++;
                distribution.areaweightedrawzeros += areaweight;
            }

            if (storedseasonalrain == 0.0f)
            {
                distribution.storedzeros++;
                distribution.areaweightedstoredzeros += areaweight;
            }

            if (meanmonthlyrain < 1.0f)
            {
                distribution.belowone++;
                distribution.areaweightedbelowone += areaweight;
            }
        }
    }

    auto finalisedistribution = [](distributionaccumulator& distribution)
    {
        auto& metrics = distribution.metrics;
        const double inversecellcount = metrics.cells > 0
            ? 1.0 / static_cast<double>(metrics.cells)
            : 0.0;
        const double inverseareaweight = metrics.areaWeight > 0.0
            ? 1.0 / metrics.areaWeight
            : 0.0;

        metrics.rawZeroFraction = distribution.rawzeros * inversecellcount;
        metrics.storedZeroFraction = distribution.storedzeros * inversecellcount;
        metrics.belowOneMillimetreFraction = distribution.belowone * inversecellcount;
        metrics.areaWeightedRawZeroFraction = distribution.areaweightedrawzeros * inverseareaweight;
        metrics.areaWeightedStoredZeroFraction = distribution.areaweightedstoredzeros * inverseareaweight;
        metrics.areaWeightedBelowOneMillimetreFraction = distribution.areaweightedbelowone * inverseareaweight;
        metrics.meanMonthlyPrecipitationMm = distribution.precipitationtotal * inversecellcount;
        metrics.areaWeightedMeanMonthlyPrecipitationMm =
            distribution.areaweightedprecipitationtotal * inverseareaweight;

        std::sort(distribution.precipitation.begin(), distribution.precipitation.end());
        const size_t wettestcount = static_cast<size_t>(std::ceil(
            static_cast<double>(distribution.precipitation.size()) * 0.10));
        double wettesttotal = 0.0;

        for (size_t index = distribution.precipitation.size() - wettestcount;
            index < distribution.precipitation.size(); index++)
        {
            wettesttotal += distribution.precipitation[index];
        }

        if (distribution.precipitationtotal > 0.0)
            metrics.wettestTenPercentShare = wettesttotal / distribution.precipitationtotal;
    };

    finalisedistribution(landdistribution);
    finalisedistribution(oceandistribution);

    climatephysics::PrecipitationDistributionDiagnostics precipitationdiagnostics;
    precipitationdiagnostics.land = landdistribution.metrics;
    precipitationdiagnostics.ocean = oceandistribution.metrics;
    climatephysics::setLastPrecipitationDistributionDiagnostics(precipitationdiagnostics);
    climatephysics::setLastHydrologySpinupDiagnostics(diagnostics);

    std::cout
        << "Climate precipitation distribution land_raw_zero="
        << precipitationdiagnostics.land.areaWeightedRawZeroFraction
        << " land_stored_zero="
        << precipitationdiagnostics.land.areaWeightedStoredZeroFraction
        << " land_below_1mm="
        << precipitationdiagnostics.land.areaWeightedBelowOneMillimetreFraction
        << " land_wettest_10_percent_share="
        << precipitationdiagnostics.land.wettestTenPercentShare
        << '\n';

    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
    {
        const climatephysics::WaterBudget& budget = finalbudgets[season];
        const climatephysics::WaterBudget& areaweightedbudget = finalareaweightedbudgets[season];
        climatephysics::setLastWaterBudget(season, budget);
        climatephysics::setLastAreaWeightedWaterBudget(season, areaweightedbudget);
        climatephysics::setLastCondensationActivityDiagnostics(
            season,
            finalcondensationactivity[season]);
        climatephysics::setLastPrecipitationProcessDiagnostics(
            season,
            finalprocessdiagnostics[season]);
        std::cout
            << "Climate area-weighted water budget season=" << season
            << " initial_atmospheric_storage=" << areaweightedbudget.initialAtmosphericStorage
            << " initial_soil_storage=" << areaweightedbudget.initialSoilStorage
            << " initial_snow_storage=" << areaweightedbudget.initialSnowStorage
            << " ocean_evaporation=" << areaweightedbudget.oceanEvaporation
            << " land_evaporation=" << areaweightedbudget.landEvaporation
            << " ocean_precipitation=" << areaweightedbudget.oceanPrecipitation
            << " land_precipitation=" << areaweightedbudget.landPrecipitation
            << " runoff=" << areaweightedbudget.runoff
            << " atmospheric_storage=" << areaweightedbudget.atmosphericStorage
            << " soil_storage=" << areaweightedbudget.soilStorage
            << " snow_storage=" << areaweightedbudget.snowStorage
            << " residual=" << areaweightedbudget.residual()
            << " relative_residual=" << areaweightedbudget.relativeResidual()
            << '\n';
    }

    // Climate sampling does not initialize/advance the deferred gameplay state.
    climate.storage = { width + 1, height + 1, boundarymoisture, freemoisture, soilmoisture, snowwater, landSurface, landDeep };
    climate.hasProcesses = true;
    if (referenceqa || referenceevap || referenceradiation)
    {
        const auto path = referenceexperiment().directory() /
            ("humidity-use-year-" + std::to_string(climate.hydrologyYears) + ".csv");
        std::ofstream output(path);
        output << "season,x,y,seconds,qa_seconds,suppressed_evaporation_seconds,potential_evaporation_delta_mm,external_surface_radiation_jm2\n";
        // For evaporation-only cases qa_seconds stores rate*time and suppression records negative rates.
        output << std::setprecision(10);
        for (int season = 0; season < CLIMATESEASONCOUNT; ++season)
            for (int y = 0; y <= height; ++y)
                for (int x = 0; x <= width; ++x)
                {
                    const auto& ledger = referenceledger[season][y * (width + 1) + x];
                    output << season << ',' << x << ',' << y;
                    for (double value : ledger) output << ',' << value;
                    output << '\n';
                }
        if (!output) throw std::runtime_error("Failed to write reference humidity ledger");
    }
    if (referenceuppermodeeast)
    {
        const auto path = referenceexperiment().directory() /
            ("upper-transport-use-year-" + std::to_string(climate.hydrologyYears) + ".csv");
        std::ofstream output(path);
        output << "season,x,y,seconds,east_velocity_delta_m,north_velocity_delta_m\n" << std::setprecision(10);
        for (int season = 0; season < CLIMATESEASONCOUNT; ++season)
            for (int y = 0; y <= height; ++y)
                for (int x = 0; x <= width; ++x)
                {
                    output << season << ',' << x << ',' << y;
                    for (double value : uppertransportledger[season][y * (width + 1) + x]) output << ',' << value;
                    output << '\n';
                }
        if (!output) throw std::runtime_error("Failed to write reference upper transport ledger");
    }
    for (int season = 0; season < CLIMATESEASONCOUNT; ++season)
    {
        auto& fields = seasonalprocess[season];
        const auto normalize = [&](vector<float>& values) {
            for (float& value : values) value /= static_cast<float>(std::max(1.0, fields.durationSeconds)); };
        for (int layer = 0; layer < 2; ++layer)
        {
            normalize(fields.radiativeHeatingWm2[layer]); normalize(fields.latentHeatingWm2[layer]);
        }
        normalize(fields.sensibleHeatingWm2); normalize(fields.surfaceNetHeatingWm2); normalize(fields.surfaceTemperatureC);
        normalize(fields.landSnowCoverFraction);
        normalize(fields.permanentLandIceFraction);
        normalize(fields.surfaceRadiativeHeatingWm2); normalize(fields.surfaceLatentHeatLossWm2);
        normalize(fields.columnWaterMm); normalize(fields.ascentHpaPerDay); normalize(fields.columnEnergyResidualWm2);
        world.recordprocessdiagnostics(season, fields);
        climate.processes[season] = std::move(fields);
        for (int layer = 0; layer < 2; ++layer)
            world.recordweatherdiagnostics(season, width + 1, height + 1,
                climateweather::calculateStatistics(seasonalweather[season], layer, seasonalweights[season]), layer);
    }
}
