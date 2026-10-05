// climate.ocean_currents: the wind-driven ocean circulation on the world grid.

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

void createoceancurrentmap(climateworld& world)
{
    auto& climate = world.climatesimulation();
    const int outputcolumns = world.width() + 1;
    const int outputrows = world.height() + 1;
    const int sealevel = world.sealevel();
    const auto dimensions = climatehydrology::climateGridDimensions(
        outputcolumns,
        outputrows,
        climate.resolution.oceanColumns);
    const int columns = dimensions.columns;
    const int rows = dimensions.rows;
    const std::size_t outputcellcount = static_cast<std::size_t>(outputcolumns) * outputrows;
    const std::size_t cellcount = static_cast<std::size_t>(columns) * rows;

    const auto inputRemap = climategrid::makeConservativeRemap(outputcolumns, outputrows, columns, rows);
    const auto outputRemap = climategrid::makeSmoothRemap(columns, rows, outputcolumns, outputrows);
    auto reservoir = climate.oceanStorage;
    if (reservoir.columns != columns || reservoir.rows != rows)
        reservoir = {columns, rows};
    std::array<climateocean::OceanForcing, CLIMATESEASONCOUNT> forcingBySeason;
    std::array<climateocean::OceanConfig, CLIMATESEASONCOUNT> configs;
    std::array<std::vector<float>, CLIMATESEASONCOUNT> surfaceRadiation, surfaceLatentLoss, surfaceSensibleLoss;
    std::vector<float> outputland(outputcellcount, 0.0f);
    std::vector<std::uint8_t> outputlandmask(outputcellcount, 0);
    std::vector<float> hydrologyland(outputcellcount, 0.0f);
    for (int season = 0; season < CLIMATESEASONCOUNT; ++season)
    {
        auto& config = configs[season];
        std::vector<float> outputdepth(outputcellcount, 0.0f);
        std::vector<float> outputeastwind(outputcellcount, 0.0f);
        std::vector<float> outputsouthwind(outputcellcount, 0.0f);
        std::vector<float> outputtemperature(outputcellcount, 0.0f);
        std::vector<float> outputsst(outputcellcount, 0.0f);
        std::vector<float> outputice(outputcellcount, 0.0f), outputfluxreference(outputcellcount, 0.0f);
        for (int y = 0; y < outputrows; y++)
        {
            for (int x = 0; x < outputcolumns; x++)
            {
                const std::size_t cell = static_cast<std::size_t>(y) * outputcolumns + x;
                const bool land = world.nom(x, y) > sealevel;
                outputland[cell] = land ? 1.0f : 0.0f;
                outputlandmask[cell] = land ? 1 : 0;
                hydrologyland[cell] = world.sea(x, y) ? 0.0f : 1.0f;
                outputdepth[cell] = land
                    ? 0.0f
                    : std::max(50.0f, static_cast<float>(sealevel) - world.nomfloat(x, y));
                const bool cached = climate.winds.populated[season];
                outputeastwind[cell] = cached ? climate.winds.surfacewindu[season][x][y] : world.seasonaluwindfloat(season, x, y);
                outputsouthwind[cell] = cached ? climate.winds.surfacewindv[season][x][y] : world.seasonalvwindfloat(season, x, y);
                outputtemperature[cell] = world.seasonaltempfloat(season, x, y);
                // The base climate supplies ice-skin temperature, not supercooled
                // liquid. Prescribe a seasonal 1 m ice reservoir where frozen.
                const float freezing = config.freezingTemperatureC;
                outputsst[cell] = land ? outputtemperature[cell] : std::max(freezing, outputtemperature[cell]);
                outputice[cell] = !land && outputtemperature[cell] < freezing ? 1.0f : 0.0f;
                outputfluxreference[cell] = climate.hasOcean && !climate.oceanSkin[season].empty()
                    ? climate.oceanSkin[season][x][y] : outputtemperature[cell];
            }
        }

        const auto remapinput = [&](const std::vector<float>& field)
        {
            return climategrid::remapField(
                inputRemap, field, outputcolumns);
        };
        const auto remapwetinput = [&](const std::vector<float>& field)
        {
            return climateocean::detail::remapOceanThermalField(
                inputRemap, outputlandmask, field, std::vector<float>(cellcount, 0.0f));
        };
        auto& forcing = forcingBySeason[season];
        const std::vector<float> landfraction = remapinput(outputland);
        forcing.landMask.resize(cellcount, 0);
        for (std::size_t cell = 0; cell < cellcount; cell++)
            forcing.landMask[cell] = landfraction[cell] >= 0.5f ? 1 : 0;
        // These intensive fields describe the water portion of a mixed cell.
        // Land depth placeholders and land winds/temperatures are not donors.
        forcing.bathymetryMetres = remapwetinput(outputdepth);
        // Solve the wind-driven gyre above a quiescent abyss. This is a
        // fixed active-layer approximation, not a prognostic thermocline.
        // The diagnostic depth now reports this dynamical column; world
        // bathymetry is retained. A 300 m cap preserves the thermal column.
        if (tuning::climate::oceancurrents::windDrivenActiveDepthMetres > 0.0f)
            for (float& depth : forcing.bathymetryMetres)
                depth = std::min(depth, tuning::climate::oceancurrents::windDrivenActiveDepthMetres);
        forcing.eastWindMps = remapwetinput(outputeastwind);
        forcing.southWindMps = remapwetinput(outputsouthwind);
        forcing.atmosphericTemperatureC = remapwetinput(outputtemperature);
        forcing.initialSstC = remapwetinput(outputsst);
        for (std::size_t cell = 0; cell < cellcount; ++cell)
            if (!forcing.landMask[cell]) forcing.initialSstC[cell] = std::max(config.freezingTemperatureC, forcing.initialSstC[cell]);
        forcing.initialIceThicknessMetres = remapwetinput(outputice);
        if (reservoir.sstC.empty())
        {
            reservoir.sstC = forcing.initialSstC;
            reservoir.iceMetres = forcing.initialIceThicknessMetres;
            reservoir.deepTemperatureC = forcing.initialSstC;
            for (auto& value : reservoir.deepTemperatureC)
                value = std::max(config.freezingTemperatureC, value - config.deepWaterTemperatureContrastK);
        }
        forcing.initialSstC = reservoir.sstC;
        forcing.initialIceThicknessMetres = reservoir.iceMetres;
        forcing.deepWaterTemperatureC = reservoir.deepTemperatureC;
        forcing.initialStorageDepthMetres = reservoir.storageDepthMetres;
        forcing.initialReservoirTemperatureC = reservoir.reservoirTemperatureC;
        forcing.surfaceHeatFluxReferenceTemperatureC = remapwetinput(outputfluxreference);
        if (climate.hasProcesses && climate.processes[season].durationSeconds > 0.0)
        {
            const auto& fields = climate.processes[season];
            // Match the binary sea-cell diagnosis used to compute these fluxes.
            // Remap flux and its linearization temperature with the SAME wet
            // support; mixing land fluxes changes the ocean boundary condition.
            const auto processlandfraction = climategrid::remapField(
                outputcolumns, outputrows, hydrologyland, fields.columns, fields.rows);
            std::vector<std::uint8_t> processland(processlandfraction.size());
            for (std::size_t cell = 0; cell < processland.size(); ++cell)
                processland[cell] = processlandfraction[cell] >= 0.5f ? 1 : 0;
            const auto processRemap = climategrid::makeConservativeRemap(
                fields.columns, fields.rows, columns, rows);
            std::vector<float> fallbackflux(cellcount);
            for (std::size_t cell = 0; cell < cellcount; ++cell)
                fallbackflux[cell] = tuning::climate::oceancurrents::surfaceHeatExchangeWm2K *
                    (forcing.atmosphericTemperatureC[cell] - forcing.surfaceHeatFluxReferenceTemperatureC[cell]);
            forcing.surfaceHeatFluxWm2 = climateocean::detail::remapOceanThermalField(
                processRemap, processland, fields.surfaceNetHeatingWm2, fallbackflux);
            forcing.surfaceHeatFluxReferenceTemperatureC = climateocean::detail::remapOceanThermalField(
                processRemap, processland, fields.surfaceTemperatureC, forcing.surfaceHeatFluxReferenceTemperatureC);
            const std::vector<float> missing(cellcount, 0.0f);
            surfaceRadiation[season] = climateocean::detail::remapOceanThermalField(
                processRemap, processland, fields.surfaceRadiativeHeatingWm2, missing);
            surfaceLatentLoss[season] = climateocean::detail::remapOceanThermalField(
                processRemap, processland, fields.surfaceLatentHeatLossWm2, missing);
            surfaceSensibleLoss[season] = climateocean::detail::remapOceanThermalField(
                processRemap, processland, fields.sensibleHeatingWm2, missing);
        }

        config.rotationDirection = world.rotation() ? 1.0f : -1.0f;
        config.rotationRatePerSecond = tuning::climate::atmosphere::rotationRatePerSecond;
        config.planetRadiusMetres = tuning::climate::atmosphere::referencePlanetRadiusMetres;
        config.airDensityKgM3 = tuning::climate::atmosphere::surfaceAirDensityKgM3;
        config.dragCoefficient = tuning::climate::atmosphere::oceanMomentumDragCoefficient;
        config.waterDensityKgM3 = tuning::climate::oceancurrents::waterDensityKgM3;
        config.barotropicDragPerSecond =
            tuning::climate::oceancurrents::barotropicDragPerSecond;
        config.minimumCoriolisPerSecond = tuning::climate::oceancurrents::ekmanDampingPerSecond;
        config.linearBottomDragMps = tuning::climate::oceancurrents::linearBottomDragMps;
        config.mixedLayerDepthMetres =
            tuning::climate::oceancurrents::mixedLayerDepthMetres;
        config.variableHeatStorage = tuning::climate::oceancurrents::variableHeatStorage;
        config.storageStratificationPerSecond2 = tuning::climate::oceancurrents::storageStratificationPerSecond2;
        config.storageWindMixingEfficiency = tuning::climate::oceancurrents::storageWindMixingEfficiency;
        config.minimumStorageDepthMetres = tuning::climate::oceancurrents::minimumStorageDepthMetres;
        config.heatDiffusivityM2S = tuning::climate::oceancurrents::heatDiffusivityM2S;
        config.surfaceHeatExchangeWm2K =
            tuning::climate::oceancurrents::surfaceHeatExchangeWm2K;
        config.streamfunctionIterations =
            tuning::climate::oceancurrents::streamfunctionIterations;
        config.heatStepsPerIteration =
            tuning::climate::oceancurrents::heatStepsPerCouplingIteration;
        config.couplingIterations = tuning::climate::oceancurrents::couplingIterations;
        config.underRelaxation =
            tuning::climate::oceancurrents::couplingUnderRelaxation;
        config.convergenceTolerance = tuning::climate::oceancurrents::couplingTolerance;
        config.sstWindFeedbackMpsPerK =
            tuning::climate::oceancurrents::sstWindFeedbackMpsPerK;
        config.maximumCurrentMps = tuning::climate::oceancurrents::maximumCurrentMps;
        // Atmosphere feedback is handled by the shared outer climate solve,
        // never by the ocean's legacy SST-to-pressure surrogate.
        config.oneWay = true;
        const int hours = std::clamp(climate.oceanHeatStepHours, 1, 24);
        const int stepsPerDay = (24 + hours - 1) / hours;
        config.oceanTimeStepSeconds = 86400.0f / stepsPerDay;
        config.reuseSeasonalCirculation = climate.reuseOceanCirculation;
        config.krylovCirculation = climate.krylovOceanCirculation;
        config.gatherHeatFluxes = climate.gatherOceanHeatFluxes;
        config.heatWorkers = climate.oceanHeatWorkers;
        config.coastalScheme = climate.oceanCoastalScheme;
        config.heatStepsPerIteration = 0;
        for (int month = season * 3; month < season * 3 + 3; ++month)
            config.heatStepsPerIteration += climatehydrology::calendarMonth(month).days * stepsPerDay;
    }
    if (!climate.oceanCaptureDirectory.empty())
        for (int season = 0; season < CLIMATESEASONCOUNT; ++season)
            std::cout << "Ocean capture path=" << climateocean::detail::OceanReplay{
                columns, rows, season, configs[season], forcingBySeason[season]}.write(climate.oceanCaptureDirectory).string() << '\n';
    auto periodic = climateocean::solvePeriodicOcean(columns, rows, forcingBySeason, configs);
    if (periodic.accepted)
    {
        reservoir.sstC = periodic.seasons.back().sstC;
        reservoir.iceMetres = periodic.seasons.back().iceThicknessMetres;
        if (configs[0].variableHeatStorage)
        {
            reservoir.storageDepthMetres = periodic.seasons.back().storageDepthMetres;
            reservoir.reservoirTemperatureC = periodic.seasons.back().reservoirTemperatureC;
        }
        reservoir.annualEnthalpyDriftK = periodic.annualEnthalpyDriftK;
        reservoir.periodic = periodic.converged;
        climate.oceanStorage = std::move(reservoir);
    }
    else climate.oceanStorage.periodic = false;
    std::cout << "Ocean spinup years=" << periodic.years << " annual_enthalpy_drift_k="
        << periodic.annualEnthalpyDriftK << " periodic=" << periodic.converged << '\n';
    std::cout << "Ocean profile columns=" << columns << " rows=" << rows
        << " circulation_seconds=" << periodic.circulationSeconds << " heat_seconds=" << periodic.heatSeconds
        << " zonal_seconds=" << periodic.zonalHeatSeconds
        << " flux_seconds=" << periodic.heatFluxSeconds << " update_seconds=" << periodic.heatUpdateSeconds
        << " heat_kernel=" << (climate.gatherOceanHeatFluxes ? "gather-v1" : "scatter-v1")
        << " heat_workers=" << periodic.seasons[0].heatWorkersUsed
        << " coastal_scheme=" << static_cast<int>(climate.oceanCoastalScheme)
        << " circulation_solves=" << periodic.circulationSolves << " heat_steps=" << periodic.heatSubsteps
        << " step_hours=" << climate.oceanHeatStepHours << " reuse=" << climate.reuseOceanCirculation
        << " solver=" << (climate.krylovOceanCirculation ? "bicgstab-sgs-v1" : "relaxation-v1") << '\n';
    for (int season = 0; season < CLIMATESEASONCOUNT; ++season)
    {
        auto& ocean = periodic.seasons[season];
        const auto& config = configs[season];
        const auto& forcing = forcingBySeason[season];
        climate.oceanAccepted[season] = climateocean::usableOceanState(ocean, cellcount);
        const auto& limit = ocean.stabilityLimit;
        if (limit.cell >= 0)
            std::cout << "Ocean timestep season=" << season << " columns=" << columns << " rows=" << rows
                << " cell=" << limit.cell << " longitude=" << (-180.0 + 360.0 * (limit.cell % columns + 0.5) / columns)
                << " latitude=" << (90.0 - 180.0 * (limit.cell / columns + 0.5) / rows)
                << " rate_s=" << limit.totalRate << " stable_seconds=" << (0.7 / limit.totalRate)
                << " actual_seconds=" << ocean.actualHeatStepSeconds << " zonal_implicit=" << config.implicitZonalHeatDiffusion
                << " zonal_advection_s=" << limit.zonalAdvectionRate << " meridional_advection_s=" << limit.meridionalAdvectionRate
                << " zonal_diffusion_s=" << limit.zonalDiffusionRate << " meridional_diffusion_s=" << limit.meridionalDiffusionRate
                << " vertical_s=" << limit.verticalRate << " surface_s=" << limit.surfaceRate << '\n';
        const auto& distribution = ocean.stabilityDistribution;
        std::cout << "Ocean rate distribution season=" << season << " wet_cells=" << distribution.wetCells
            << " wet_area_m2=" << distribution.wetAreaM2
            << " p50_rate_s=" << distribution.rateQuantiles[0] << " p95_rate_s=" << distribution.rateQuantiles[1]
            << " p99_rate_s=" << distribution.rateQuantiles[2]
            << " without_advection_rate_s=" << distribution.maximumRateWithoutAdvection << '\n';
        constexpr std::array<int, 5> diagnosticHours{1, 3, 6, 12, 24};
        for (std::size_t bin = 0; bin < diagnosticHours.size(); ++bin)
            std::cout << "Ocean rate threshold season=" << season << " hours=" << diagnosticHours[bin]
                << " cells=" << distribution.restrictedCells[bin]
                << " area_fraction=" << (distribution.restrictedAreaM2[bin] / std::max(1.0, distribution.wetAreaM2)) << '\n';
        constexpr std::array<const char*, 4> faceNames{"east", "west", "south", "north"};
        for (std::size_t rank = 0; rank < distribution.highestRates.size(); ++rank)
        {
            const auto& hot = distribution.highestRates[rank];
            if (hot.cell < 0) continue;
            std::cout << "Ocean hotspot season=" << season << " rank=" << rank + 1 << " cell=" << hot.cell
                << " longitude=" << (-180.0 + 360.0 * (hot.cell % columns + 0.5) / columns)
                << " latitude=" << (90.0 - 180.0 * (hot.cell / columns + 0.5) / rows)
                << " rate_s=" << hot.totalRate << " depth_m=" << hot.depthMetres
                << " zonal_advection_s=" << hot.zonalAdvectionRate << " meridional_advection_s=" << hot.meridionalAdvectionRate
                << " vertical_s=" << hot.verticalRate;
            for (std::size_t face = 0; face < faceNames.size(); ++face)
                std::cout << ' ' << faceNames[face] << "_mps=" << hot.faceVelocityMps[face]
                    << ' ' << faceNames[face] << "_barotropic_mps=" << hot.barotropicFaceVelocityMps[face]
                    << ' ' << faceNames[face] << "_ekman_mps=" << hot.faceVelocityMps[face] - hot.barotropicFaceVelocityMps[face]
                    << ' ' << faceNames[face] << "_depth_m=" << hot.faceDepthMetres[face];
            std::cout << '\n';
        }
        if (!climate.oceanAccepted[season])
        {
            ocean.converged = false;
            std::cout << "Ocean fallback: season=" << season << " iterations=" << ocean.couplingIterations
                << " residual=" << ocean.relativeResidual << " basin_residual=" << ocean.streamfunctionRelativeResidual << '\n';
            ocean.sstC = forcing.initialSstC;
            ocean.surfaceSkinTemperatureC = forcing.atmosphericTemperatureC;
            ocean.iceThicknessMetres = forcing.initialIceThicknessMetres;
            ocean.eastCurrentMps.assign(cellcount, 0.0f);
            ocean.southCurrentMps.assign(cellcount, 0.0f);
            ocean.ekmanUpwellingMps.assign(cellcount, 0.0f);
            ocean.coupledEastWindMps = forcing.eastWindMps;
            ocean.coupledSouthWindMps = forcing.southWindMps;
            ocean.coupledPressureAnomalyHpa.assign(cellcount, 0.0f);
        }

        if (!climate.oceanAccepted[season])
        {
            ocean.meanSstC = ocean.sstC;
            ocean.meanSurfaceSkinTemperatureC = ocean.surfaceSkinTemperatureC;
            ocean.meanIceCover = ocean.iceThicknessMetres;
            for (auto& h : ocean.meanIceCover) h = std::clamp(h / 0.10f, 0.0f, 1.0f);
        }
        world.recordoceandiagnostics(season, columns, rows, ocean);
        // Log the forcing actually used, before the next hydrology update.
        // Sequential groups of four seasons identify outer ocean calls.
        if (climate.oceanAccepted[season] && !surfaceRadiation[season].empty())
            for (int band = 0; band < 4; ++band)
            {
                std::array<double, 11> sums{};
                double area = 0.0, maximumClosure = 0.0;
                int cells = 0;
                for (int y = 0; y < rows; ++y)
                    for (int x = 0; x < columns; ++x)
                    {
                        const int cell = y * columns + x;
                        const double latitude = std::abs(90.0 - 180.0 * (y + 0.5) / rows);
                        if (forcing.landMask[cell] || forcing.bathymetryMetres[cell] <= 0.0f ||
                            (band == 1 && latitude >= 20.0) ||
                            (band == 2 && (latitude < 20.0 || latitude >= 60.0)) ||
                            (band == 3 && latitude < 60.0)) continue;
                        const double weight = climategrid::latitudeBandMeasure(y, rows);
                        const auto& budget = ocean.heatBudget[cell];
                        const std::array<double, 11> values = {
                            forcing.surfaceHeatFluxReferenceTemperatureC[cell], forcing.surfaceHeatFluxWm2[cell],
                            surfaceRadiation[season][cell], surfaceLatentLoss[season][cell], surfaceSensibleLoss[season][cell],
                            ocean.meanSstC[cell], ocean.meanSurfaceSkinTemperatureC[cell], budget.surfaceWm2,
                            budget.horizontalWm2 + budget.verticalWm2, budget.storageWm2, forcing.atmosphericTemperatureC[cell]};
                        for (std::size_t k = 0; k < sums.size(); ++k) sums[k] += weight * values[k];
                        maximumClosure = std::max(maximumClosure, std::abs(values[1] - values[2] + values[3] + values[4]));
                        area += weight; ++cells;
                    }
                if (area == 0.0) continue;
                constexpr const char* labels[] = {"reference_c", "reference_flux_wm2", "radiation_wm2", "latent_loss_wm2",
                    "sensible_loss_wm2", "sst_c", "skin_c", "applied_surface_wm2", "net_transport_wm2", "storage_wm2", "air_c"};
                std::ostringstream trace;
                trace << std::setprecision(10) << "Ocean applied season=" << season << " band=" << band << " cells=" << cells;
                for (std::size_t k = 0; k < sums.size(); ++k) trace << ' ' << labels[k] << '=' << sums[k] / area;
                trace << " maximum_forcing_closure_wm2=" << maximumClosure;
                trace << " zonal_diffusion=" << (config.implicitZonalHeatDiffusion ? "implicit-v1" : "explicit-v1")
                    << " heat_substeps=" << ocean.heatSubsteps << " explicit_heat_substeps=" << ocean.explicitHeatSubsteps
                    << " zonal_iterations=" << ocean.maximumZonalDiffusionIterations
                    << " zonal_residual_k=" << ocean.maximumZonalDiffusionResidualK;
                std::cout << trace.str() << '\n';
            }
        const auto remapoutput = [&](const std::vector<float>& field)
        {
            return climategrid::remapField(
                outputRemap, field, columns);
        };
        const std::vector<float> eastcurrent = remapoutput(ocean.eastCurrentMps);
        const std::vector<float> southcurrent = remapoutput(ocean.southCurrentMps);
        // These are temperatures/cover over water, not extensive whole-cell
        // quantities. Land placeholders must not cool or de-ice coastal water.
        // Small seas with no resolved wet donor retain the base thermal state.
        std::vector<float> fallbackSst(outputcellcount), fallbackSkin(outputcellcount), fallbackIce(outputcellcount);
        for (int y = 0; y < outputrows; ++y)
            for (int x = 0; x < outputcolumns; ++x)
            {
                const auto cell = y * outputcolumns + x;
                const float base = world.seasonaltempfloat(season, x, y);
                fallbackSst[cell] = std::max(config.freezingTemperatureC, base);
                fallbackSkin[cell] = base;
                fallbackIce[cell] = base < config.freezingTemperatureC ? 1.0f : 0.0f;
            }
        auto sst = climateocean::detail::remapOceanThermalField(
            outputRemap, forcing.landMask, ocean.meanSstC, fallbackSst);
        for (std::size_t cell = 0; cell < outputcellcount; ++cell)
            if (outputland[cell] == 0.0f) sst[cell] = std::max(config.freezingTemperatureC, sst[cell]);
        const auto skin = climateocean::detail::remapOceanThermalField(
            outputRemap, forcing.landMask, ocean.meanSurfaceSkinTemperatureC, fallbackSkin);
        const auto ice = climateocean::detail::remapOceanThermalField(
            outputRemap, forcing.landMask, ocean.meanIceCover, fallbackIce);
        climate.hasOcean = true;
        climate.oceanSst[season].assign(outputcolumns, vector<float>(outputrows));
        climate.oceanSkin[season].assign(outputcolumns, vector<float>(outputrows));
        climate.oceanIce[season].assign(outputcolumns, vector<float>(outputrows));
        for (int y = 0; y < outputrows; ++y)
            for (int x = 0; x < outputcolumns; ++x)
            {
                const auto cell = y * outputcolumns + x;
                climate.oceanSst[season][x][y] = sst[cell];
                climate.oceanSkin[season][x][y] = skin[cell];
                climate.oceanIce[season][x][y] = ice[cell];
            }
        for (int y = 0; y < outputrows; y++)
        {
            for (int x = 0; x < outputcolumns; x++)
            {
                const std::size_t cell = static_cast<std::size_t>(y) * outputcolumns + x;
                if (world.nom(x, y) > sealevel)
                {
                    world.setseasonalcurrentu(season, x, y, 0);
                    world.setseasonalcurrentv(season, x, y, 0);
                    world.setseasonalsst(season, x, y, 0);
                    continue;
                }
                world.setseasonalcurrentufloat(season, x, y, eastcurrent[cell]);
                world.setseasonalcurrentvfloat(season, x, y, southcurrent[cell]);
                world.setseasonalsstfloat(season, x, y, sst[cell]);
                // This is a one-way ocean solve. Its coarser forcing is never
                // written back over the accepted atmospheric wind/pressure.
            }
        }
        std::cout
            << "Coupled ocean season " << season
            << ": iterations=" << ocean.couplingIterations
            << " residual=" << ocean.relativeResidual
            << " streamfunction_residual=" << ocean.streamfunctionRelativeResidual
            << " converged=" << ocean.converged
            << " max_divergence=" << ocean.maximumTransportDivergenceMps
            << " heat_residual_j=" << ocean.heatBudgetResidualJ
            << " relative_heat_residual=" << ocean.relativeHeatBudgetResidual << '\n';
    }
    if (!tuning::climate::oceancurrents::oneWayDiagnosticsOnly)
        for (int y = 0; y < outputrows; ++y)
            for (int x = 0; x < outputcolumns; ++x)
            {
                int frozenSeasons = 0;
                for (int season = 0; season < CLIMATESEASONCOUNT; ++season)
                    frozenSeasons += climate.oceanIce[season][x][y] >= 0.5f;
                world.setseaice(x, y, world.sea(x, y)
                    ? (frozenSeasons == CLIMATESEASONCOUNT ? 2 : (frozenSeasons > 0 ? 1 : 0)) : 0);
            }
}
