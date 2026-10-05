// Diagnostics of the climate solvers that no world field carries: the coupling
// iterations, the energy and water budgets, precipitation statistics, the
// atmospheric budget and the precision of the circulation. One CSV table each,
// written at the end of a climate run into the diagnostics directory of
// configs/app.env. Comparisons with observed data are made from run bundles by
// scripts/evaluation.

#include "simulations/climate/climate_diagnostics.hpp"

#include "core/app_environment.hpp"
#include "core/climate_grid.hpp"
#include "core/climate_seasons.hpp"
#include "pipeline/world_generation_debug.hpp"
#include "simulations/climate/climate_atmosphere.hpp"
#include "simulations/climate/climate_coupling.hpp"
#include "simulations/climate/climate_energy.hpp"
#include "simulations/climate/climate_physics.hpp"
#include "simulations/climate/climate_world.hpp"
#include "simulations/climate/climate_parameters.hpp"

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <string>

using namespace std;

namespace
{
double gridcellareaweight(int row, int height)
{
    if (height <= 0)
        return 0.0;

    return climategrid::latitudeBandMeasure(row, height + 1);
}

float latitudeforrow(int row, int height)
{
    if (height <= 0)
        return 0.0f;

    return 90.0f - 180.0f * (static_cast<float>(row) + 0.5f) /
        static_cast<float>(height + 1);
}

filesystem::path climatediagnosticsdirectory()
{
    filesystem::path outputroot = getappenvironment().climateDiagnosticsDirectory;

    if (outputroot.empty())
        outputroot = filesystem::current_path();

    outputroot /= "seed_" + to_string(worldgenerationdebugseed());
    filesystem::create_directories(outputroot);

    return outputroot;
}

void writeatmosphericbudget(const filesystem::path& outputdir, climateworld& world)
{
    constexpr array<const char*, CLIMATESEASONCOUNT> seasonnames = { "january", "april", "july", "october" };
    const int width = world.width();
    const int height = world.height();
    ofstream outfile(outputdir / "climate_atmosphere_budget.csv");

    if (!outfile.is_open())
        return;

    outfile << "season,area_weighted_mean_pressure_anomaly_hpa,rms_pressure_anomaly_hpa,rms_zonal_mean_pressure_anomaly_hpa,rms_stationary_pressure_anomaly_hpa,pressure_limiter_fraction,area_weighted_mean_surface_divergence_per_day,area_weighted_mean_upper_divergence_per_day,area_weighted_mean_vertical_velocity_hpa_per_day,rms_surface_wind_m_s,surface_wind_limiter_fraction,rms_upper_wind_m_s,rms_vertical_velocity_hpa_per_day,vertical_limiter_fraction,rms_column_divergence_per_day,mean_orographic_uplift_m,maximum_orographic_uplift_m,uplift_above_500m_fraction,mean_orographic_descent_m,maximum_orographic_descent_m,descent_above_500m_fraction\n";
    outfile << fixed << setprecision(9);

    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
    {
        double weighttotal = 0.0;
        double pressuretotal = 0.0;
        double pressuresquared = 0.0;
        double zonalpressuresquared = 0.0;
        double stationarypressuresquared = 0.0;
        double pressurelimitedweight = 0.0;
        double surfacedivergencetotal = 0.0;
        double upperdivergencetotal = 0.0;
        double verticaltotal = 0.0;
        double surfacewindsquared = 0.0;
        double upperwindsquared = 0.0;
        double verticalsquared = 0.0;
        double columndivergencesquared = 0.0;
        double surfacewindlimitedweight = 0.0;
        double verticallimitedweight = 0.0;
        double upliftmetrestotal = 0.0;
        double descentmetrestotal = 0.0;
        double upliftabove500weight = 0.0;
        double descentabove500weight = 0.0;
        double maximumupliftmetres = 0.0;
        double maximumdescentmetres = 0.0;
        vector<double> zonalpressure(height + 1, 0.0);

        for (int y = 0; y <= height; y++)
        {
            for (int x = 0; x <= width; x++)
                zonalpressure[y] += world.seasonalpressurefloat(season, x, y);

            zonalpressure[y] /= static_cast<double>(width + 1);
        }

        for (int y = 0; y <= height; y++)
        {
            const int ynorth = y > 0 ? y - 1 : y;
            const int ysouth = y < height ? y + 1 : y;
            const double latitude = latitudeforrow(y, height);
            const double centrecosine = (std::max)(
                0.02,
                std::abs(std::cos(latitude * 3.14159265358979323846 / 180.0)));
            const double northcosine = (std::max)(
                0.0,
                std::cos(latitudeforrow(ynorth, height) * 3.14159265358979323846 / 180.0));
            const double southcosine = (std::max)(
                0.0,
                std::cos(latitudeforrow(ysouth, height) * 3.14159265358979323846 / 180.0));
            const double weight = gridcellareaweight(y, height);
            const auto spacing = climateatmosphere::cellSpacingMetres(
                latitude,
                width + 1,
                height + 1,
                tuning::climate::atmosphere::referencePlanetRadiusMetres);

            for (int x = 0; x <= width; x++)
            {
                const int xwest = x > 0 ? x - 1 : width - 1;
                const int xeast = x < width ? x + 1 : 1;
                const double surfaceu = world.seasonaluwindfloat(season, x, y);
                const double surfacev = world.seasonalvwindfloat(season, x, y);
                const double upperu = world.seasonalupperuwindfloat(season, x, y);
                const double upperv = world.seasonaluppervwindfloat(season, x, y);
                const double surfacedivergence =
                    ((world.seasonaluwindfloat(season, xeast, y) - world.seasonaluwindfloat(season, xwest, y)) /
                        (2.0 * spacing.zonalMetres) +
                        (world.seasonalvwindfloat(season, x, ysouth) * southcosine -
                            world.seasonalvwindfloat(season, x, ynorth) * northcosine) /
                        (2.0 * spacing.meridionalMetres * centrecosine)) *
                    tuning::climate::circulation::secondsPerDay;
                const double upperdivergence =
                    ((world.seasonalupperuwindfloat(season, xeast, y) - world.seasonalupperuwindfloat(season, xwest, y)) /
                        (2.0 * spacing.zonalMetres) +
                        (world.seasonaluppervwindfloat(season, x, ysouth) * southcosine -
                            world.seasonaluppervwindfloat(season, x, ynorth) * northcosine) /
                        (2.0 * spacing.meridionalMetres * centrecosine)) *
                    tuning::climate::circulation::secondsPerDay;
                const double verticalvelocity =
                    world.seasonalverticalvelocityfloat(season, x, y);
                const double columndivergence = surfacedivergence + upperdivergence;
                const double upliftmetres = world.seasonalupliftfloat(season, x, y);
                const double descentmetres = world.seasonalsubsidencefloat(season, x, y);

                weighttotal += weight;
                pressuretotal += weight * world.seasonalpressurefloat(season, x, y);
                pressuresquared += weight * world.seasonalpressurefloat(season, x, y) *
                    world.seasonalpressurefloat(season, x, y);
                zonalpressuresquared += weight * zonalpressure[y] * zonalpressure[y];
                const double stationarypressure =
                    world.seasonalpressurefloat(season, x, y) - zonalpressure[y];
                stationarypressuresquared +=
                    weight * stationarypressure * stationarypressure;
                surfacedivergencetotal += weight * surfacedivergence;
                upperdivergencetotal += weight * upperdivergence;
                verticaltotal += weight * verticalvelocity;
                surfacewindsquared += weight * (surfaceu * surfaceu + surfacev * surfacev);
                upperwindsquared += weight * (upperu * upperu + upperv * upperv);
                verticalsquared += weight * verticalvelocity * verticalvelocity;
                columndivergencesquared += weight * columndivergence * columndivergence;
                upliftmetrestotal += weight * upliftmetres;
                descentmetrestotal += weight * descentmetres;
                maximumupliftmetres = (std::max)(maximumupliftmetres, upliftmetres);
                maximumdescentmetres = (std::max)(maximumdescentmetres, descentmetres);
                upliftabove500weight += upliftmetres >= 500.0 ? weight : 0.0;
                descentabove500weight += descentmetres >= 500.0 ? weight : 0.0;

                if (std::abs(world.seasonalpressurefloat(season, x, y)) >=
                    tuning::climate::circulation::maximumSurfacePressureAnomalyHpa - 1.0f)
                {
                    pressurelimitedweight += weight;
                }

                if (std::sqrt(surfaceu * surfaceu + surfacev * surfacev) >=
                    tuning::climate::atmosphere::maxVectorWind - 1.0f)
                {
                    surfacewindlimitedweight += weight;
                }

                if (std::abs(verticalvelocity) >=
                    tuning::climate::circulation::maximumVerticalVelocity - 1.0f)
                {
                    verticallimitedweight += weight;
                }
            }
        }

        const double inverseweight = weighttotal > 0.0 ? 1.0 / weighttotal : 0.0;
        outfile
            << seasonnames[season] << ','
            << pressuretotal * inverseweight << ','
            << sqrt(pressuresquared * inverseweight) << ','
            << sqrt(zonalpressuresquared * inverseweight) << ','
            << sqrt(stationarypressuresquared * inverseweight) << ','
            << pressurelimitedweight * inverseweight << ','
            << surfacedivergencetotal * inverseweight << ','
            << upperdivergencetotal * inverseweight << ','
            << verticaltotal * inverseweight << ','
            << sqrt(surfacewindsquared * inverseweight) << ','
            << surfacewindlimitedweight * inverseweight << ','
            << sqrt(upperwindsquared * inverseweight) << ','
            << sqrt(verticalsquared * inverseweight) << ','
            << verticallimitedweight * inverseweight << ','
            << sqrt(columndivergencesquared * inverseweight) << ','
            << upliftmetrestotal * inverseweight << ','
            << maximumupliftmetres << ','
            << upliftabove500weight * inverseweight << ','
            << descentmetrestotal * inverseweight << ','
            << maximumdescentmetres << ','
            << descentabove500weight * inverseweight << '\n';
    }
}

void writecirculationprecision(const filesystem::path& outputdir)
{
    constexpr array<const char*, CLIMATESEASONCOUNT> seasonnames = {
        "january", "april", "july", "october"
    };
    ofstream output(outputdir / "climate_circulation_precision.csv");

    if (!output.is_open())
        return;

    output << "season,stage,float_pressure_gradient_rms_pa_m,rounded_pressure_gradient_rms_pa_m,pressure_gradient_difference_rms_pa_m,float_surface_wind_rms_m_s,rounded_surface_wind_rms_m_s,surface_wind_difference_rms_m_s,float_upper_height_gradient_rms,rounded_upper_height_gradient_rms,upper_height_gradient_difference_rms,float_upper_wind_rms_m_s,rounded_upper_wind_rms_m_s,upper_wind_difference_rms_m_s\n";
    output << fixed << setprecision(9);
    const auto& diagnostics = climateatmosphere::lastCirculationPrecisionDiagnostics();

    auto writerow = [&](int season, const char* stage,
        const climateatmosphere::CirculationPrecisionStageDiagnostics& values)
    {
        output
            << seasonnames[season] << ',' << stage << ','
            << values.areaWeightedFloatPressureGradientRmsPaPerMetre << ','
            << values.areaWeightedRoundedPressureGradientRmsPaPerMetre << ','
            << values.areaWeightedPressureGradientDifferenceRmsPaPerMetre << ','
            << values.areaWeightedFloatSurfaceWindRmsMetresPerSecond << ','
            << values.areaWeightedRoundedSurfaceWindRmsMetresPerSecond << ','
            << values.areaWeightedSurfaceWindDifferenceRmsMetresPerSecond << ','
            << values.areaWeightedFloatUpperHeightGradientRms << ','
            << values.areaWeightedRoundedUpperHeightGradientRms << ','
            << values.areaWeightedUpperHeightGradientDifferenceRms << ','
            << values.areaWeightedFloatUpperWindRmsMetresPerSecond << ','
            << values.areaWeightedRoundedUpperWindRmsMetresPerSecond << ','
            << values.areaWeightedUpperWindDifferenceRmsMetresPerSecond << '\n';
    };

    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
    {
        writerow(season, "base", diagnostics[season].base);
        writerow(season, "final", diagnostics[season].final);
    }
}
}

void writeclimatediagnostics(climateworld& world)
{
    if (isworldgendebugrunactive() == false)
        return;

    const filesystem::path outputdir = climatediagnosticsdirectory();

    ofstream couplingfile(outputdir / "climate_coupling.csv");
    couplingfile << "iteration,wind_change,sst_change,heating_change,rainfall_change,inner_accepted,converged,heating_relaxation\n";
    for (const auto& d : climatephysics::lastClimateCouplingDiagnostics())
        couplingfile << d.iteration << ',' << d.windChange << ',' << d.sstChange << ',' << d.heatingChange << ','
            << d.rainfallChange << ',' << d.innerSolvesAccepted << ',' << d.converged << ',' << d.heatingRelaxation << '\n';

    ofstream energybudgetfile(outputdir / "climate_energy_budget.csv");

    if (energybudgetfile.is_open())
    {
        const climateenergy::AnnualEnergyBudget& budget = climateenergy::lastAnnualEnergyBudget();
        energybudgetfile << "incoming_solar_wm2,absorbed_solar_wm2,outgoing_longwave_wm2,atmospheric_transport_wm2,storage_tendency_wm2,residual_wm2,calibrated_longwave_intercept_wm2,area_weighted_mean_temperature_c,area_weighted_land_temperature_c,area_weighted_ocean_temperature_c,area_weighted_land_elevation_cooling_c,area_weighted_permanent_land_ice_fraction,local_ice_coupling_iterations,local_ice_coupling_maximum_row_residual,local_ice_coupling_temperature_error_c,local_ice_coupling_converged,northern_50_70_warmest_temperature_c,southern_60_90_warmest_temperature_c,southern_60_90_icecap_thermal_fraction\n";
        energybudgetfile << fixed << setprecision(9)
            << budget.incomingSolarWm2 << ','
            << budget.absorbedSolarWm2 << ','
            << budget.outgoingLongwaveWm2 << ','
            << budget.atmosphericTransportWm2 << ','
            << budget.storageTendencyWm2 << ','
            << budget.residualWm2 << ','
            << budget.calibratedLongwaveInterceptWm2 << ','
            << budget.areaWeightedMeanTemperatureC << ','
            << budget.areaWeightedLandTemperatureC << ','
            << budget.areaWeightedOceanTemperatureC << ','
            << budget.areaWeightedLandElevationCoolingC << ','
            << budget.areaWeightedPermanentLandIceFraction << ','
            << budget.localIceCouplingIterations << ','
            << budget.localIceCouplingMaximumRowResidual << ','
            << budget.localIceCouplingTemperatureErrorC << ','
            << (budget.localIceCouplingConverged ? 1 : 0) << ','
            << budget.northern5070WarmestTemperatureC << ','
            << budget.southern6090WarmestTemperatureC << ','
            << budget.southern6090IcecapThermalFraction << '\n';
    }

    ofstream waterbudgetfile(outputdir / "climate_water_budget.csv");

    if (waterbudgetfile.is_open())
    {
        waterbudgetfile << "quarter,initial_atmospheric_storage,initial_soil_storage,initial_snow_storage,ocean_evaporation,land_evaporation,ocean_precipitation,land_precipitation,runoff,atmospheric_storage,soil_storage,snow_storage,residual,relative_residual\n";
        waterbudgetfile << fixed << setprecision(9);

        const auto& budgets = climatephysics::lastWaterBudgets();

        for (int season = 0; season < CLIMATESEASONCOUNT; season++)
        {
            const climatephysics::WaterBudget& budget = budgets[season];
            waterbudgetfile
                << season << ','
                << budget.initialAtmosphericStorage << ','
                << budget.initialSoilStorage << ','
                << budget.initialSnowStorage << ','
                << budget.oceanEvaporation << ','
                << budget.landEvaporation << ','
                << budget.oceanPrecipitation << ','
                << budget.landPrecipitation << ','
                << budget.runoff << ','
                << budget.atmosphericStorage << ','
                << budget.soilStorage << ','
                << budget.snowStorage << ','
                << budget.residual() << ','
                << budget.relativeResidual() << '\n';
        }
    }

    ofstream areaweightedwaterbudgetfile(outputdir / "climate_water_budget_area_weighted.csv");

    if (areaweightedwaterbudgetfile.is_open())
    {
        areaweightedwaterbudgetfile << "quarter,initial_atmospheric_storage,initial_soil_storage,initial_snow_storage,ocean_evaporation,land_evaporation,ocean_precipitation,land_precipitation,runoff,atmospheric_storage,soil_storage,snow_storage,residual,relative_residual\n";
        areaweightedwaterbudgetfile << fixed << setprecision(9);

        const auto& budgets = climatephysics::lastAreaWeightedWaterBudgets();

        for (int season = 0; season < CLIMATESEASONCOUNT; season++)
        {
            const climatephysics::WaterBudget& budget = budgets[season];
            areaweightedwaterbudgetfile
                << season << ','
                << budget.initialAtmosphericStorage << ','
                << budget.initialSoilStorage << ','
                << budget.initialSnowStorage << ','
                << budget.oceanEvaporation << ','
                << budget.landEvaporation << ','
                << budget.oceanPrecipitation << ','
                << budget.landPrecipitation << ','
                << budget.runoff << ','
                << budget.atmosphericStorage << ','
                << budget.soilStorage << ','
                << budget.snowStorage << ','
                << budget.residual() << ','
                << budget.relativeResidual() << '\n';
        }
    }

    ofstream precipitationdistributionfile(outputdir / "climate_precipitation_distribution.csv");

    if (precipitationdistributionfile.is_open())
    {
        precipitationdistributionfile << "scope,cells,area_weight,raw_zero_fraction,stored_zero_fraction,below_1mm_fraction,area_weighted_raw_zero_fraction,area_weighted_stored_zero_fraction,area_weighted_below_1mm_fraction,mean_monthly_precipitation_mm,area_weighted_mean_monthly_precipitation_mm,maximum_monthly_precipitation_mm,maximum_x,maximum_y,maximum_month,rounded_precipitation_at_maximum_mm,wettest_10_percent_share\n";
        precipitationdistributionfile << fixed << setprecision(9);
        const auto& diagnostics = climatephysics::lastPrecipitationDistributionDiagnostics();

        auto writerow = [&](const char* scope, const climatephysics::PrecipitationDistributionScope& values)
        {
            precipitationdistributionfile
                << scope << ','
                << values.cells << ','
                << values.areaWeight << ','
                << values.rawZeroFraction << ','
                << values.storedZeroFraction << ','
                << values.belowOneMillimetreFraction << ','
                << values.areaWeightedRawZeroFraction << ','
                << values.areaWeightedStoredZeroFraction << ','
                << values.areaWeightedBelowOneMillimetreFraction << ','
                << values.meanMonthlyPrecipitationMm << ','
                << values.areaWeightedMeanMonthlyPrecipitationMm << ','
                << values.maximumMonthlyPrecipitationMm << ','
                << values.maximumMonthlyPrecipitationX << ','
                << values.maximumMonthlyPrecipitationY << ','
                << values.maximumMonthlyPrecipitationMonth << ','
                << values.roundedPrecipitationAtMaximumMm << ','
                << values.wettestTenPercentShare << '\n';
        };

        writerow("land", diagnostics.land);
        writerow("ocean", diagnostics.ocean);
    }

    ofstream spinupfile(outputdir / "climate_hydrology_spinup.csv");

    if (spinupfile.is_open())
    {
        const climatephysics::HydrologySpinupDiagnostics& diagnostics =
            climatephysics::lastHydrologySpinupDiagnostics();
        spinupfile << "cycles_completed,converged,relative_storage_change,relative_atmospheric_storage_change,relative_soil_storage_change,relative_snow_storage_change,relative_snow_cover_change,atmospheric_storage,soil_storage,snow_storage\n";
        spinupfile << fixed << setprecision(9)
            << diagnostics.cyclesCompleted << ','
            << (diagnostics.converged ? 1 : 0) << ','
            << diagnostics.relativeStorageChange << ','
            << diagnostics.relativeAtmosphericStorageChange << ','
            << diagnostics.relativeSoilStorageChange << ','
            << diagnostics.relativeSnowStorageChange << ','
            << diagnostics.relativeSnowCoverChange << ','
            << diagnostics.atmosphericStorage << ','
            << diagnostics.soilStorage << ','
            << diagnostics.snowStorage << '\n';
    }

    ofstream condensationfile(outputdir / "climate_condensation_activity.csv");

    if (condensationfile.is_open())
    {
        constexpr array<const char*, CLIMATESEASONCOUNT> seasonnames = {
            "january-march", "april-june", "july-september", "october-december"
        };
        const auto& activity = climatephysics::lastCondensationActivityDiagnostics();
        const auto& budgets = climatephysics::lastAreaWeightedWaterBudgets();
        constexpr double representativeperioddays = 365.0 / CLIMATESEASONCOUNT;
        condensationfile << "season,cell_step_area_weight,active_cell_step_fraction,active_atmospheric_water_fraction,condensable_excess_water_fraction,precipitation_to_excess_fraction,atmospheric_water_residence_time_days\n";
        condensationfile << fixed << setprecision(9);

        for (int season = 0; season < CLIMATESEASONCOUNT; season++)
        {
            const double precipitation = budgets[season].oceanPrecipitation +
                budgets[season].landPrecipitation;
            condensationfile
                << seasonnames[season] << ','
                << activity[season].cellStepAreaWeight << ','
                << (activity[season].cellStepAreaWeight > 0.0 ?
                    activity[season].activeCellStepAreaWeight /
                        activity[season].cellStepAreaWeight : 0.0) << ','
                << (activity[season].atmosphericWaterAreaWeighted > 0.0 ?
                    activity[season].activeAtmosphericWaterAreaWeighted /
                        activity[season].atmosphericWaterAreaWeighted : 0.0) << ','
                << (activity[season].atmosphericWaterAreaWeighted > 0.0 ?
                    activity[season].excessWaterAreaWeighted /
                        activity[season].atmosphericWaterAreaWeighted : 0.0) << ','
                << (activity[season].excessWaterAreaWeighted > 0.0 ?
                    activity[season].precipitationAreaWeighted /
                        activity[season].excessWaterAreaWeighted : 0.0) << ','
                << (precipitation > 0.0 ?
                    representativeperioddays * budgets[season].atmosphericStorage /
                        precipitation : 0.0) << '\n';
        }
    }

    ofstream precipitationprocessfile(outputdir / "climate_precipitation_processes.csv");

    if (precipitationprocessfile.is_open())
    {
        constexpr array<const char*, CLIMATESEASONCOUNT> quarternames = {
            "january-march", "april-june", "july-september", "october-december"
        };
        const auto& processes = climatephysics::lastPrecipitationProcessDiagnostics();
        const auto& budgets = climatephysics::lastAreaWeightedWaterBudgets();
        precipitationprocessfile
            << "quarter,stratiform_precipitation_mm_area_weighted,orographic_precipitation_mm_area_weighted,convective_precipitation_mm_area_weighted,total_precipitation_mm_area_weighted,water_budget_precipitation_mm_area_weighted,component_budget_difference_mm_area_weighted,reevaporated_precipitation_mm_area_weighted,snowfall_mm_area_weighted,upward_moisture_transfer_mm_area_weighted,downward_moisture_transfer_mm_area_weighted,cloud_fraction_area_weighted,positive_moisture_flux_convergence_mm_area_weighted,negative_moisture_flux_convergence_mm_area_weighted,net_moisture_flux_convergence_mm_area_weighted\n";
        precipitationprocessfile << fixed << setprecision(9);

        for (int season = 0; season < CLIMATESEASONCOUNT; season++)
        {
            const auto& process = processes[season];
            const double processtotal = process.stratiformPrecipitation +
                process.orographicPrecipitation + process.convectivePrecipitation;
            const double budgettotal = budgets[season].oceanPrecipitation +
                budgets[season].landPrecipitation;
            precipitationprocessfile
                << quarternames[season] << ','
                << process.stratiformPrecipitation << ','
                << process.orographicPrecipitation << ','
                << process.convectivePrecipitation << ','
                << processtotal << ','
                << budgettotal << ','
                << processtotal - budgettotal << ','
                << process.reevaporatedPrecipitation << ','
                << process.snowfall << ','
                << process.upwardMoistureTransfer << ','
                << process.downwardMoistureTransfer << ','
                << process.cloudFractionAreaWeighted << ','
                << process.positiveMoistureFluxConvergence << ','
                << process.negativeMoistureFluxConvergence << ','
                << process.positiveMoistureFluxConvergence +
                    process.negativeMoistureFluxConvergence << '\n';
        }
    }

    writeatmosphericbudget(outputdir, world);
    writecirculationprecision(outputdir);
}
