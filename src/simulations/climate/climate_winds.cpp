// climate.winds: surface and upper winds, terrain vertical motion and the precision of the circulation.

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

namespace
{
float samplewrappedfield(const floatgrid& field, climateworld& world, float x, float y)
{
    const int width = world.width();
    const int height = world.height();
    const float span = static_cast<float>(width + 1);

    while (x < 0.0f)
        x += span;

    while (x >= span)
        x -= span;

    y = std::clamp(y, 0.0f, static_cast<float>(height));

    const int x0 = static_cast<int>(std::floor(x));
    const int y0 = static_cast<int>(std::floor(y));
    const int x1 = wrapx(x0 + 1, width);
    const int y1 = std::min(y0 + 1, height);
    const float fracx = x - static_cast<float>(x0);
    const float fracy = y - static_cast<float>(y0);
    const float v00 = field[x0][y0];
    const float v10 = field[x1][y0];
    const float v01 = field[x0][y1];
    const float v11 = field[x1][y1];

    return
        v00 * (1.0f - fracx) * (1.0f - fracy) +
        v10 * fracx * (1.0f - fracy) +
        v01 * (1.0f - fracx) * fracy +
        v11 * fracx * fracy;
}

void storeterrainverticalmotion(climateworld& world, int season, const floatgrid& macroterrain, const floatgrid& windu, const floatgrid& windv)
{
    const int width = world.width();
    const int height = world.height();

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int y = startrow; y <= endrow; y++)
        {
            const int ynorth = (y > 0) ? y - 1 : y;
            const int ysouth = (y < height) ? y + 1 : y;

            for (int x = 0; x <= width; x++)
            {
                const float u = windu[x][y];
                const float v = windv[x][y];
                const float speed = std::sqrt(u * u + v * v);

                if (speed < tuning::climate::atmosphere::topographyMinimumWindSpeed)
                {
                    world.setseasonaluplift(season, x, y, 0);
                    world.setseasonalsubsidence(season, x, y, 0);
                    continue;
                }

                const int xwest = wrapx(x - 1, width);
                const int xeast = wrapx(x + 1, width);
                const float latitude = latitudeforrow(y, height);
                const auto spacing = climateatmosphere::cellSpacingMetres(
                    latitude,
                    width + 1,
                    height + 1,
                    tuning::climate::atmosphere::referencePlanetRadiusMetres);
                const float elevationchangeeast =
                    (macroterrain[xeast][y] - macroterrain[xwest][y]) / 2.0f;
                const float elevationchangesouth =
                    (macroterrain[x][ysouth] - macroterrain[x][ynorth]) / 2.0f;
                const float gradmag = std::sqrt(
                    elevationchangeeast * elevationchangeeast +
                    elevationchangesouth * elevationchangesouth);
                const float terrainhere = macroterrain[x][y];

                if (gradmag < tuning::climate::atmosphere::topographyMinimumRelief && terrainhere < tuning::climate::atmosphere::topographyMinimumRelief)
                {
                    world.setseasonaluplift(season, x, y, 0);
                    world.setseasonalsubsidence(season, x, y, 0);
                    continue;
                }

                const float slopeeast = elevationchangeeast / spacing.zonalMetres;
                const float slopesouth = elevationchangesouth / spacing.meridionalMetres;
                const float verticalvelocity = u * slopeeast + v * slopesouth;
                const float directioneast = u / speed;
                const float directionsouth = v / speed;
                float windwardrise = 0.0f;
                float leewarddrop = 0.0f;

                for (int step = 1;
                    step <= tuning::climate::atmosphere::topographyLookaheadDistance;
                    step++)
                {
                    const float samplex = static_cast<float>(x) +
                        directioneast * static_cast<float>(step);
                    const float sampley = static_cast<float>(y) +
                        directionsouth * static_cast<float>(step);
                    const float elevationchange =
                        samplewrappedfield(macroterrain, world, samplex, sampley) - terrainhere;
                    windwardrise = std::max(windwardrise, elevationchange);
                    leewarddrop = std::max(leewarddrop, -elevationchange);
                }

                const float inversecellcrossingtime = std::sqrt(
                    (u / spacing.zonalMetres) * (u / spacing.zonalMetres) +
                    (v / spacing.meridionalMetres) * (v / spacing.meridionalMetres));
                const float cellcrossingtime = inversecellcrossingtime > 0.0f ?
                    1.0f / inversecellcrossingtime : 0.0f;
                float parceldisplacement = verticalvelocity * cellcrossingtime;

                if (parceldisplacement > 0.0f)
                    parceldisplacement = std::max(parceldisplacement, windwardrise);
                else if (parceldisplacement < 0.0f)
                    parceldisplacement = std::min(parceldisplacement, -leewarddrop);

                parceldisplacement = std::clamp(
                    parceldisplacement,
                    -tuning::climate::atmosphere::maximumOrographicParcelDisplacementMetres,
                    tuning::climate::atmosphere::maximumOrographicParcelDisplacementMetres);
                const float uplift = std::max(0.0f, parceldisplacement);
                const float subsidence = std::max(0.0f, -parceldisplacement);

                world.setseasonalupliftfloat(season, x, y, uplift);
                world.setseasonalsubsidencefloat(season, x, y, subsidence);
            }
        }
    });
}

climateatmosphere::CirculationPrecisionStageDiagnostics measurecirculationprecision(
    climateworld& world,
    const floatgrid& pressure,
    const floatgrid& upperheight,
    const floatgrid& continentality,
    const floatgrid& macroterrain)
{
    const int width = world.width();
    const int height = world.height();
    double weighttotal = 0.0;
    double floatpressuregradientsquared = 0.0;
    double roundedpressuregradientsquared = 0.0;
    double pressuregradientdifferencesquared = 0.0;
    double floatwindsquared = 0.0;
    double roundedwindsquared = 0.0;
    double winddifferencesquared = 0.0;
    double floatheightgradientsquared = 0.0;
    double roundedheightgradientsquared = 0.0;
    double heightgradientdifferencesquared = 0.0;
    double floatupperwindsquared = 0.0;
    double roundedupperwindsquared = 0.0;
    double upperwinddifferencesquared = 0.0;

    for (int y = 0; y <= height; y++)
    {
        const int ynorth = y > 0 ? y - 1 : y;
        const int ysouth = y < height ? y + 1 : y;
        const float latitude = latitudeforrow(y, height);
        const double weight = rowareaweight(y, height);
        const auto spacing = climateatmosphere::cellSpacingMetres(
            latitude,
            width + 1,
            height + 1,
            tuning::climate::atmosphere::referencePlanetRadiusMetres);

        for (int x = 0; x <= width; x++)
        {
            const int xwest = wrapx(x - 1, width);
            const int xeast = wrapx(x + 1, width);
            const float floatpressureeast =
                (pressure[xeast][y] - pressure[xwest][y]) *
                tuning::climate::atmosphere::pressurePascalsPerHectopascal /
                (2.0f * spacing.zonalMetres);
            const float floatpressurenorth =
                (pressure[x][ynorth] - pressure[x][ysouth]) *
                tuning::climate::atmosphere::pressurePascalsPerHectopascal /
                (2.0f * spacing.meridionalMetres);
            const float roundedpressureeast =
                (std::round(pressure[xeast][y]) - std::round(pressure[xwest][y])) *
                tuning::climate::atmosphere::pressurePascalsPerHectopascal /
                (2.0f * spacing.zonalMetres);
            const float roundedpressurenorth =
                (std::round(pressure[x][ynorth]) - std::round(pressure[x][ysouth])) *
                tuning::climate::atmosphere::pressurePascalsPerHectopascal /
                (2.0f * spacing.meridionalMetres);
            const float continental = std::clamp(continentality[x][y], 0.0f, 1.0f);
            const float relief = std::clamp(macroterrain[x][y] / 4500.0f, 0.0f, 1.0f);
            const float basedragcoefficient =
                tuning::climate::atmosphere::oceanMomentumDragCoefficient +
                continental * (tuning::climate::atmosphere::landMomentumDragCoefficient -
                    tuning::climate::atmosphere::oceanMomentumDragCoefficient);
            const float dragcoefficient = basedragcoefficient + relief *
                (tuning::climate::atmosphere::highReliefMomentumDragCoefficient -
                    basedragcoefficient);
            const auto floatwind = climateatmosphere::steadyQuadraticDragCoriolisWind(
                -floatpressureeast / tuning::climate::atmosphere::surfaceAirDensityKgM3,
                -floatpressurenorth / tuning::climate::atmosphere::surfaceAirDensityKgM3,
                latitude,
                dragcoefficient,
                tuning::climate::atmosphere::surfaceBoundaryLayerMomentumDepthMetres,
                tuning::climate::atmosphere::rotationRatePerSecond,
                world.rotation() ? 1.0f : -1.0f);
            const auto roundedwind = climateatmosphere::steadyQuadraticDragCoriolisWind(
                -roundedpressureeast / tuning::climate::atmosphere::surfaceAirDensityKgM3,
                -roundedpressurenorth / tuning::climate::atmosphere::surfaceAirDensityKgM3,
                latitude,
                dragcoefficient,
                tuning::climate::atmosphere::surfaceBoundaryLayerMomentumDepthMetres,
                tuning::climate::atmosphere::rotationRatePerSecond,
                world.rotation() ? 1.0f : -1.0f);
            const float floatheighteast =
                (upperheight[xeast][y] - upperheight[xwest][y]) /
                (2.0f * spacing.zonalMetres);
            const float floatheightnorth =
                (upperheight[x][ynorth] - upperheight[x][ysouth]) /
                (2.0f * spacing.meridionalMetres);
            const float roundedheighteast =
                (std::round(upperheight[xeast][y]) - std::round(upperheight[xwest][y])) /
                (2.0f * spacing.zonalMetres);
            const float roundedheightnorth =
                (std::round(upperheight[x][ynorth]) - std::round(upperheight[x][ysouth])) /
                (2.0f * spacing.meridionalMetres);
            const auto floatupperwind = climateatmosphere::steadyRayleighCoriolisWind(
                -tuning::climate::atmosphere::gravityMetresPerSecondSquared * floatheighteast,
                -tuning::climate::atmosphere::gravityMetresPerSecondSquared * floatheightnorth,
                latitude,
                tuning::climate::circulation::upperLayerDragTimeSeconds,
                tuning::climate::atmosphere::rotationRatePerSecond,
                world.rotation() ? 1.0f : -1.0f);
            const auto roundedupperwind = climateatmosphere::steadyRayleighCoriolisWind(
                -tuning::climate::atmosphere::gravityMetresPerSecondSquared * roundedheighteast,
                -tuning::climate::atmosphere::gravityMetresPerSecondSquared * roundedheightnorth,
                latitude,
                tuning::climate::circulation::upperLayerDragTimeSeconds,
                tuning::climate::atmosphere::rotationRatePerSecond,
                world.rotation() ? 1.0f : -1.0f);
            const double pressureeastdifference = roundedpressureeast - floatpressureeast;
            const double pressurenorthdifference = roundedpressurenorth - floatpressurenorth;
            const double windeastdifference = roundedwind.eastMetresPerSecond - floatwind.eastMetresPerSecond;
            const double windsouthdifference = roundedwind.southMetresPerSecond - floatwind.southMetresPerSecond;
            const double heighteastdifference = roundedheighteast - floatheighteast;
            const double heightnorthdifference = roundedheightnorth - floatheightnorth;
            const double upperwindeastdifference = roundedupperwind.eastMetresPerSecond - floatupperwind.eastMetresPerSecond;
            const double upperwindsouthdifference = roundedupperwind.southMetresPerSecond - floatupperwind.southMetresPerSecond;

            weighttotal += weight;
            floatpressuregradientsquared += weight *
                (floatpressureeast * floatpressureeast + floatpressurenorth * floatpressurenorth);
            roundedpressuregradientsquared += weight *
                (roundedpressureeast * roundedpressureeast + roundedpressurenorth * roundedpressurenorth);
            pressuregradientdifferencesquared += weight *
                (pressureeastdifference * pressureeastdifference + pressurenorthdifference * pressurenorthdifference);
            floatwindsquared += weight *
                (floatwind.eastMetresPerSecond * floatwind.eastMetresPerSecond +
                    floatwind.southMetresPerSecond * floatwind.southMetresPerSecond);
            roundedwindsquared += weight *
                (roundedwind.eastMetresPerSecond * roundedwind.eastMetresPerSecond +
                    roundedwind.southMetresPerSecond * roundedwind.southMetresPerSecond);
            winddifferencesquared += weight *
                (windeastdifference * windeastdifference + windsouthdifference * windsouthdifference);
            floatheightgradientsquared += weight *
                (floatheighteast * floatheighteast + floatheightnorth * floatheightnorth);
            roundedheightgradientsquared += weight *
                (roundedheighteast * roundedheighteast + roundedheightnorth * roundedheightnorth);
            heightgradientdifferencesquared += weight *
                (heighteastdifference * heighteastdifference + heightnorthdifference * heightnorthdifference);
            floatupperwindsquared += weight *
                (floatupperwind.eastMetresPerSecond * floatupperwind.eastMetresPerSecond +
                    floatupperwind.southMetresPerSecond * floatupperwind.southMetresPerSecond);
            roundedupperwindsquared += weight *
                (roundedupperwind.eastMetresPerSecond * roundedupperwind.eastMetresPerSecond +
                    roundedupperwind.southMetresPerSecond * roundedupperwind.southMetresPerSecond);
            upperwinddifferencesquared += weight *
                (upperwindeastdifference * upperwindeastdifference +
                    upperwindsouthdifference * upperwindsouthdifference);
        }
    }

    const double inverseweight = weighttotal > 0.0 ? 1.0 / weighttotal : 0.0;
    climateatmosphere::CirculationPrecisionStageDiagnostics diagnostics;
    diagnostics.areaWeightedFloatPressureGradientRmsPaPerMetre =
        std::sqrt(floatpressuregradientsquared * inverseweight);
    diagnostics.areaWeightedRoundedPressureGradientRmsPaPerMetre =
        std::sqrt(roundedpressuregradientsquared * inverseweight);
    diagnostics.areaWeightedPressureGradientDifferenceRmsPaPerMetre =
        std::sqrt(pressuregradientdifferencesquared * inverseweight);
    diagnostics.areaWeightedFloatSurfaceWindRmsMetresPerSecond =
        std::sqrt(floatwindsquared * inverseweight);
    diagnostics.areaWeightedRoundedSurfaceWindRmsMetresPerSecond =
        std::sqrt(roundedwindsquared * inverseweight);
    diagnostics.areaWeightedSurfaceWindDifferenceRmsMetresPerSecond =
        std::sqrt(winddifferencesquared * inverseweight);
    diagnostics.areaWeightedFloatUpperHeightGradientRms =
        std::sqrt(floatheightgradientsquared * inverseweight);
    diagnostics.areaWeightedRoundedUpperHeightGradientRms =
        std::sqrt(roundedheightgradientsquared * inverseweight);
    diagnostics.areaWeightedUpperHeightGradientDifferenceRms =
        std::sqrt(heightgradientdifferencesquared * inverseweight);
    diagnostics.areaWeightedFloatUpperWindRmsMetresPerSecond =
        std::sqrt(floatupperwindsquared * inverseweight);
    diagnostics.areaWeightedRoundedUpperWindRmsMetresPerSecond =
        std::sqrt(roundedupperwindsquared * inverseweight);
    diagnostics.areaWeightedUpperWindDifferenceRmsMetresPerSecond =
        std::sqrt(upperwinddifferencesquared * inverseweight);
    return diagnostics;
}
}

void createvectorwindmap(climateworld& world)
{
    auto& climate = world.climatesimulation();
    using namespace tuning::climate;
    const int width = world.width(), height = world.height();
    const int columns = climatehydrology::climateGridDimensions(width + 1, height + 1,
        climate.resolution.atmosphereColumns).columns;
    const int rows = columns / 2;
    if (columns < 6) return;
    const auto grid = climategrid::makeSphericalGrid(columns, rows, atmosphere::referencePlanetRadiusMetres);
    const auto input = climategrid::makeConservativeRemap(
        width + 1, height + 1, columns, rows);
    const auto sample = [&](auto accessor)
    {
        vector<float> result(columns * rows, 0.0f);
        for (int y = 0; y < rows; ++y)
            for (int x = 0; x < columns; ++x)
            {
                double sum = 0.0;
                for (const auto& wy : input.latitude[y])
                    for (const auto& wx : input.longitude[x])
                        sum += wy.fraction * wx.fraction * accessor(wx.source, wy.source);
                result[y * columns + x] = static_cast<float>(sum);
            }
        return result;
    };
    const auto output = [&](const vector<float>& field)
    {
        const auto flat = climategrid::remapSmoothField(
            columns, rows, field, width + 1, height + 1);
        floatgrid result(width + 1, vector<float>(height + 1));
        for (int y = 0; y <= height; ++y)
            for (int x = 0; x <= width; ++x) result[x][y] = flat[y * (width + 1) + x];
        return result;
    };
    const auto land = sample([&](int x, int y) { return world.nom(x, y) > world.sealevel() ? 1.0f : 0.0f; });
    vector<float> vegetationcontrast((width + 1) * (height + 1), 1.0f);
    if (atmosphere::vegetationDragContrastFraction > 0.0f)
    {
        const char* path = std::getenv("UW_CLIMATE_SURFACE_DRAG_INPUT");
        if (path && *path)
        {
            std::ifstream stream(path);
            vegetationcontrast = climateio::readSurfaceDragContrast(stream, width + 1, height + 1);
            std::cout << "Prescribed surface drag input=" << path << " strength="
                << atmosphere::vegetationDragContrastFraction << '\n';
        }
        else std::cout << "No prescribed vegetation input; retaining uniform land drag.\n";
    }
    auto nativevegetation = sample([&](int x, int y) {
        return world.nom(x,y) > world.sealevel() ? vegetationcontrast[y * (width + 1) + x] : 0.0f;
    });
    for (std::size_t cell = 0; cell < nativevegetation.size(); ++cell)
        nativevegetation[cell] = land[cell] > 0 ? nativevegetation[cell] / land[cell] : 1.0f;
    const auto terrain = sample([&](int x, int y) { return std::max(0.0f, world.nomfloat(x, y) - static_cast<float>(world.sealevel())); });
    const auto slopes = climateatmosphere::detail::terrainSlopes(grid, terrain,
        atmosphere::katabaticReflectedPolarGradient);
    const auto permanentice = sample([&](int x, int y) { return worldpermanentlandicefraction(world,x,y); });
    vector<float> persistentsnow(columns * rows, 0.0f);
    vector<float> persistentsnowoutput((width + 1) * (height + 1), 0.0f);
    if (atmosphere::persistentSnowDragFraction > 0.0f && climate.hasProcesses)
    {
        const auto& first = climate.processes[0];
        std::array<const vector<float>*, CLIMATESEASONCOUNT> quarters{};
        for (int season = 0; season < CLIMATESEASONCOUNT; ++season)
        {
            const auto& fields = climate.processes[season];
            if (fields.durationSeconds > 0.0 && fields.columns == first.columns && fields.rows == first.rows)
                quarters[season] = &fields.landSnowCoverFraction;
        }
        if (first.columns > 0 && first.rows > 0)
        {
            const auto cover = climateatmosphere::detail::persistentSnowCover(
                quarters, static_cast<std::size_t>(first.columns) * first.rows);
            persistentsnow = climategrid::remapField(first.columns, first.rows, cover, columns, rows);
            persistentsnowoutput = climategrid::remapSmoothField(first.columns, first.rows, cover, width + 1, height + 1);
        }
    }
    // Compare the reduced pressure mode's thermal columns at a common level.
    // This removes only the lapse offset imposed by the temperature module;
    // actual temperatures used for energy, snow and cold drainage are retained.
    const auto pressurelapseoffset = sample([&](int x, int y)
    {
        return world.sea(x,y) ? 0.0f : circulation::thermalPressureElevationCorrectionFraction *
            world.tempdecrease() * std::max(0.0f, world.mapfloat(x,y) - static_cast<float>(world.sealevel())) / 1000.0f;
    });
    const auto terrainoutput = output(terrain);
    const auto landoutput = output(land);
    // Preserve resolved world ridges for moisture uplift; the pressure solver
    // still uses its own conservative atmosphere-grid terrain above.
    floatgrid orographicterrain(width + 1, vector<float>(height + 1));
    for (int y = 0; y <= height; ++y)
        for (int x = 0; x <= width; ++x)
            orographicterrain[x][y] = std::max(0.0f, world.nomfloat(x, y) - static_cast<float>(world.sealevel()));
    const auto surfacedrag = [](float landFraction, float elevationMetres, float frozenLandFraction, float vegetationContrast)
    {
        return climateatmosphere::detail::vegetationSurfaceDragCoefficient(landFraction, elevationMetres,
            atmosphere::frozenSurfaceDragFraction * frozenLandFraction, atmosphere::oceanMomentumDragCoefficient,
            atmosphere::landMomentumDragCoefficient, atmosphere::highReliefMomentumDragCoefficient,
            atmosphere::frozenMomentumDragCoefficient, vegetationContrast, atmosphere::vegetationDragContrastFraction);
    };
    vector<float> parallelDrag;
    if (circulation::surfaceDragParallelMobility)
    {
        vector<float> localDrag((width + 1) * (height + 1));
        for (int y = 0; y <= height; ++y)
            for (int x = 0; x <= width; ++x)
            {
                const int cell = y * (width + 1) + x;
                localDrag[cell] = surfacedrag(world.nom(x, y) > world.sealevel() ? 1.0f : 0.0f,
                    std::max(0.0f, world.nomfloat(x, y) - static_cast<float>(world.sealevel())),
                    std::max(worldpermanentlandicefraction(world, x, y),
                        atmosphere::persistentSnowDragFraction * persistentsnowoutput[cell]), vegetationcontrast[cell]);
            }
        parallelDrag = climateatmosphere::detail::parallelSurfaceDragCoefficients(
            width + 1, height + 1, localDrag, columns, rows);
    }
    climate.winds.reset(world.width(), world.height());
    for (int season = 0; season < CLIMATESEASONCOUNT; ++season)
    {
        const auto temperature = sample([&](int x, int y) { return pressuresurfacetemperature(world, season, x, y); });
        const auto rain = sample([&](int x, int y) { return std::max(0.0f, world.seasonalrainfloat(season, x, y)); });
        vector<float> rowtemperature(rows, 0.0f), zonalpressure(rows, 0.0f);
        double globaltemperature = 0.0, totalarea = 0.0;
        for (int y = 0; y < rows; ++y)
        {
            for (int x = 0; x < columns; ++x) rowtemperature[y] += temperature[y * columns + x] / columns;
            globaltemperature += grid.cellAreasSquareMetres[y] * rowtemperature[y];
            totalarea += grid.cellAreasSquareMetres[y];
        }
        globaltemperature /= totalarea;
        // A stationary climate response needs a fixed seasonal background.
        // Re-selecting a warmest SST row on every coupling pass moves the entire
        // Hadley cell discontinuously when neighbouring warm rows exchange rank.
        // Local thermal pressure gradients below still respond to evolving SST.
        auto& radiative = climate.radiativeZonalTemperature[season];
        if (radiative.size() != rowtemperature.size()) radiative = rowtemperature;
        const auto& background = climate.transientWeather ? rowtemperature : radiative;
        double backgroundMean = 0.0;
        for (int y = 0; y < rows; ++y) backgroundMean += grid.cellAreasSquareMetres[y] * background[y] / totalarea;
        const int warmrow = static_cast<int>(std::max_element(background.begin(), background.end()) - background.begin());
        const float contrast = std::max(1.0f, background[warmrow] - 0.5f * (background.front() + background.back()));
        climateatmosphere::ModeSeparatedCirculationConfig config;
        config.adjustZonalSurfacePressure = circulation::adjustZonalSurfacePressure;
        config.adjustZonalUpperHeight = circulation::adjustZonalUpperHeight;
        config.boundaryLayer.enabled = climate.mixedLayerTropics;
        config.surfaceBoundaryLayerDepthMetres = atmosphere::surfaceBoundaryLayerMomentumDepthMetres;
        config.gravityMetresPerSecondSquared *= std::max(0.05f, world.gravity());
        config.rotationDirection = world.rotation() ? 1.0f : -1.0f;
        config.enabled = {circulation::enableZonalMode, circulation::enableStationaryMode,
            circulation::enableSurfaceMode, circulation::enableUpperMode};
        config.surfaceToUpperCoupling = circulation::surfaceToUpperModeCoupling;
        config.interlayerMomentumCoupling = circulation::stationaryInterlayerMomentumCoupling;
        config.upperDragTimeSeconds = circulation::upperLayerDragTimeSeconds;
        config.maximumIterations = circulation::stationaryWaveMaximumIterations;
        config.stationarySolver = climate.atmosphereSolver;
        config.stationaryCaptureDirectory = climate.atmosphereCaptureDirectory;
        config.relativeTolerance = circulation::stationaryWaveRelativeTolerance;
        const float hadleywidth = std::clamp(climateatmosphere::heldHouHadleyEdgeLatitudeDegrees(
            contrast, circulation::troposphereHeightMetres, static_cast<float>(backgroundMean + 273.15),
            config.gravityMetresPerSecondSquared, config.rotationRatePerSecond, config.planetRadiusMetres),
            circulation::minimumHadleyHalfWidthDegrees, circulation::maximumHadleyHalfWidthDegrees);
        const float thermalshift = hadleywidth * circulation::maximumThermalEquatorShiftHadleyFraction;
        const float thermalequator = std::clamp(static_cast<float>(grid.latitudeCentresRadians[warmrow] * 180.0 / 3.141592653589793),
            -thermalshift, thermalshift);
        const float confinedshift = hadleywidth * circulation::confinedThermalEquatorShiftHadleyFraction;
        const float confinedtrough = std::clamp(static_cast<float>(grid.latitudeCentresRadians[warmrow] * 180.0 / 3.141592653589793),
            -confinedshift, confinedshift);
        const float amplitude = std::clamp(circulation::surfaceReferencePressurePa / 100.0f *
            contrast / circulation::referenceTemperatureK * circulation::overturningMassRedistributionEfficiency,
            circulation::minimumOverturningPressureAmplitudeHpa, circulation::maximumOverturningPressureAmplitudeHpa);
        config.zonalUpperHeightMetres.resize(rows);
        const float hydrostatic = climateatmosphere::hypsometricHeightResponseMetresPerKelvin(
            circulation::surfaceReferencePressurePa, circulation::upperReferencePressurePa, 287.05f,
            config.gravityMetresPerSecondSquared);
        const auto upperthermal = climateatmosphere::detail::projectUpperThermalGradient(
            background, grid.cellAreasSquareMetres, circulation::upperThermalCoreGradientFraction,
            1.0f, 0.7f * hadleywidth, 1.7f * hadleywidth,
            circulation::upperThermalPolarGradientFraction,
            circulation::upperThermalPolarStartLatitude, circulation::upperThermalPolarEndLatitude,
            circulation::upperThermalJetGradientBoost);
        double pressuremean = 0.0;
        vector<float> continentalPressure(columns * rows, 0.0f);
        for (int y = 0; y < rows; ++y)
        {
            zonalpressure[y] = climateatmosphere::detail::confinedOverturningPressureHpa(
                static_cast<float>(grid.latitudeCentresRadians[y] * 180.0 / 3.141592653589793),
                thermalequator, hadleywidth, circulation::winterHadleyEdgeAnchorFraction, confinedtrough, amplitude);
            const float localcontrast = std::max(1.0f, background[warmrow] -
                (grid.latitudeCentresRadians[y] >= 0 ? background.front() : background.back()));
            zonalpressure[y] = climateatmosphere::detail::polarBranchPressureHpa(zonalpressure[y],
                static_cast<float>(grid.latitudeCentresRadians[y] * 180.0 / 3.141592653589793),
                thermalequator, hadleywidth, circulation::winterHadleyEdgeAnchorFraction,
                amplitude, localcontrast / contrast, circulation::polarThermalContrastFraction);
            zonalpressure[y] = climateatmosphere::detail::deepenSubpolarLowPressureHpa(zonalpressure[y],
                static_cast<float>(grid.latitudeCentresRadians[y] * 180.0 / 3.141592653589793),
                thermalequator, hadleywidth, circulation::winterHadleyEdgeAnchorFraction,
                amplitude, circulation::subpolarLowDeepeningFraction,
                circulation::subpolarDeepeningPreservePolarGradient);
            zonalpressure[y] = climateatmosphere::detail::hemisphericSubpolarPressureHpa(
                zonalpressure[y], static_cast<float>(grid.latitudeCentresRadians[y] * 180.0 / 3.141592653589793),
                thermalequator, hadleywidth, circulation::winterHadleyEdgeAnchorFraction,
                localcontrast / contrast, circulation::subpolarHemisphericContrastExponent);
            double continentalMean = 0.0;
            for (int x = 0; x < columns; ++x)
            {
                const int cell = y * columns + x;
                continentalPressure[cell] = climateatmosphere::detail::continentalSubpolarCorrectionHpa(
                    zonalpressure[y], static_cast<float>(grid.latitudeCentresRadians[y] * 180.0 / 3.141592653589793),
                    thermalequator, hadleywidth, circulation::winterHadleyEdgeAnchorFraction,
                    land[cell], circulation::continentalSubpolarLowSuppression);
                continentalMean += continentalPressure[cell] / columns;
            }
            // Carry the row mean in the zonal solve and the remainder in the
            // stationary solve. The global gauge below still conserves mass.
            zonalpressure[y] += static_cast<float>(continentalMean);
            for (int x = 0; x < columns; ++x)
                continentalPressure[y * columns + x] -= static_cast<float>(continentalMean);
            pressuremean += grid.cellAreasSquareMetres[y] * zonalpressure[y] / totalarea;
            config.zonalUpperHeightMetres[y] = hydrostatic * upperthermal[y];
        }
        for (float& value : zonalpressure) value -= static_cast<float>(pressuremean);

        // Stable dry/moist lapse-rate blend is diagnosed from simulated rain.
        double meanrain = 0.0;
        for (int y = 0; y < rows; ++y)
            for (int x = 0; x < columns; ++x)
                meanrain += rain[y * columns + x] * grid.cellAreasSquareMetres[y] / (totalarea * columns);
        const float humidity = std::clamp(static_cast<float>(meanrain) / lapse::rainfallHumidityScale, 0.0f, 1.0f);
        const float lapserate = (lapse::dryLapseRate + humidity * (lapse::moistLapseRate - lapse::dryLapseRate)) / 1000.0f;
        const float stability = std::max(0.001f, climateatmosphere::diagnoseBruntVaisalaFrequency(
            static_cast<float>(globaltemperature + 273.15), lapserate, config.gravityMetresPerSecondSquared));
        const auto surfaceparameters = climateatmosphere::diagnoseStationaryParameters(stability,
            circulation::surfaceModeDepthMetres, config.gravityMetresPerSecondSquared, config.planetRadiusMetres,
            config.rotationRatePerSecond, columns, rows, circulation::stationaryWaveNondimensionalDamping,
            circulation::diagnosedForcingScaleMetres);
        const auto upperparameters = climateatmosphere::diagnoseStationaryParameters(stability,
            circulation::upperModeDepthMetres, config.gravityMetresPerSecondSquared, config.planetRadiusMetres,
            config.rotationRatePerSecond, columns, rows, circulation::stationaryWaveNondimensionalDamping,
            circulation::diagnosedForcingScaleMetres);
        const float pressurepermetre = config.airDensityKgM3 * config.gravityMetresPerSecondSquared / 100.0f;
        config.surfaceEquivalentPressureDepthHpa = pressurepermetre * surfaceparameters.equivalentDepthMetres;
        config.upperEquivalentPressureDepthHpa = pressurepermetre * upperparameters.equivalentDepthMetres;
        config.surfaceDampingTimeSeconds = surfaceparameters.dampingTimeSeconds;
        config.upperDampingTimeSeconds = upperparameters.dampingTimeSeconds;
        config.maximumZonalWavenumber = surfaceparameters.maximumZonalWavenumber;
        config.maximumMeridionalWavenumber = surfaceparameters.maximumMeridionalWavenumber;
        config.upperMaximumZonalWavenumber = upperparameters.maximumZonalWavenumber;
        config.upperMaximumMeridionalWavenumber = upperparameters.maximumMeridionalWavenumber;
        config.surfaceDragTimesSeconds.resize(columns * rows);
        config.surfaceDragCoefficients.resize(columns * rows);
        for (int cell = 0; cell < columns * rows; ++cell)
        {
            const float drag = parallelDrag.empty() ? surfacedrag(land[cell], terrain[cell], std::max(permanentice[cell],
                atmosphere::persistentSnowDragFraction * persistentsnow[cell]), nativevegetation[cell]) : parallelDrag[cell];
            config.surfaceDragCoefficients[cell] = drag;
            config.surfaceDragTimesSeconds[cell] = atmosphere::surfaceBoundaryLayerMomentumDepthMetres /
                (drag * circulation::stationaryWaveLinearizationWindMps);
        }
        vector<float> zero(columns * rows, 0.0f);
        auto backgroundconfig = config;
        backgroundconfig.enabled = {true, false, true, true};
        const auto zonal = climateatmosphere::solveModeSeparatedCirculation(columns, rows, zonalpressure, zero, zero, backgroundconfig);
        auto surfaceforcing = climateatmosphere::mechanicalTopographicPressureForcingHpa(columns, rows,
            terrain, zonal.surfaceEastWindMps, zonal.surfaceSouthWindMps, 1.0f,
            circulation::mechanicalTopographicTerrainScaleMetres, circulation::mechanicalTopographicPressureAmplitudeHpa,
            circulation::mechanicalTopographicMinimumWindMps, circulation::mechanicalTopographicFullStrengthWindMps,
            circulation::mechanicalTopographicMinimumLatitudeDegrees, circulation::mechanicalTopographicFullStrengthLatitudeDegrees);
        const auto mechanicalforcing = surfaceforcing;
        auto thermaltemperature = temperature;
        if (circulation::surfaceThermalLandSkinFraction > 0.0f && climate.hasProcesses)
        {
            const auto& fields = climate.processes[season];
            if (fields.columns > 0 && fields.rows > 0 && fields.surfaceTemperatureC.size() ==
                static_cast<std::size_t>(fields.columns) * fields.rows)
            {
                const auto skin = climategrid::remapField(fields.columns, fields.rows,
                    fields.surfaceTemperatureC, columns, rows);
                const float fraction = std::clamp(circulation::surfaceThermalLandSkinFraction, 0.0f, 1.0f);
                for (int cell = 0; cell < columns * rows; ++cell)
                    if (std::isfinite(skin[cell]))
                        thermaltemperature[cell] += fraction * land[cell] * (skin[cell] - temperature[cell]);
            }
        }
        vector<float> thermalrowtemperature(rows, 0.0f);
        for (int cell = 0; cell < columns * rows; ++cell)
        {
            thermaltemperature[cell] += pressurelapseoffset[cell];
            thermalrowtemperature[cell / columns] += thermaltemperature[cell] / columns;
        }
        // Hydrostatic boundary-layer pressure response to the simulated surface
        // temperature, with anomalies tapering linearly to zero aloft. Net
        // column heating alone omits this temperature-gradient-driven flow
        // (Lindzen & Nigam, 1987). Keep the existing diabatic mode as well.
        for (int y = 0; y < rows; ++y)
            for (int x = 0; x < columns; ++x)
            {
                const int cell = y * columns + x;
                const float anomaly = thermaltemperature[cell] - thermalrowtemperature[y];
                // Blend integrated hydrostatic responses, not top pressures.
                // This retains the ocean response exactly and mixes coastline
                // tiles continuously as the continental column is deepened.
                surfaceforcing[cell] += continentalPressure[cell] +
                    (1.0f - land[cell]) * climateatmosphere::thermalModePressureAnomalyHpa(
                        circulation::surfaceThermalProfileFraction * anomaly, config.airDensityKgM3,
                        circulation::surfaceReferencePressurePa, circulation::surfaceThermalTopPressurePa) +
                    land[cell] * climateatmosphere::thermalModePressureAnomalyHpa(
                        circulation::surfaceThermalLandProfileFraction * anomaly, config.airDensityKgM3,
                        circulation::surfaceReferencePressurePa, circulation::surfaceThermalLandTopPressurePa);
            }
        vector<float> regionalPressure(columns * rows, 0.0f), regionalTrough(columns * rows, confinedtrough);
        if (circulation::regionalThermalTroughFraction > 0.0f)
        {
            const auto targets = climateatmosphere::detail::regionalThermalEquatorDegrees(
                columns, rows, thermaltemperature, grid.latitudeCentresRadians, confinedtrough, hadleywidth);
            for (int y = 0; y < rows; ++y)
            {
                double rowMean = 0.0;
                const float latitude = static_cast<float>(grid.latitudeCentresRadians[y] * 180.0 / 3.141592653589793);
                for (int x = 0; x < columns; ++x)
                {
                    const int cell = y * columns + x;
                    const float target = confinedtrough + circulation::regionalThermalTroughFraction * (targets[x] - confinedtrough);
                    regionalTrough[cell] = target;
                    const auto pressureAt = [&](float trough) {
                        return climateatmosphere::detail::confinedOverturningPressureHpa(latitude,
                            thermalequator, hadleywidth, circulation::winterHadleyEdgeAnchorFraction, trough, amplitude);
                    };
                    regionalPressure[cell] = pressureAt(target) - pressureAt(confinedtrough);
                    rowMean += regionalPressure[cell] / columns;
                }
                // This is a regional departure. The zonal pressure template and
                // its subtropical edges are not translated by local heating.
                for (int x = 0; x < columns; ++x)
                {
                    const int cell = y * columns + x;
                    regionalPressure[cell] -= static_cast<float>(rowMean);
                    surfaceforcing[cell] += regionalPressure[cell];
                }
            }
        }
        auto upperWaveTerrain = terrain;
        if (circulation::upperOrographicBlocking)
            for (int y = 0; y < rows; ++y)
            {
                float backgroundWind = 0.0f;
                for (int x = 0; x < columns; ++x)
                    backgroundWind += zonal.upperEastWindMps[y * columns + x] / columns;
                for (int x = 0; x < columns; ++x)
                {
                    const int cell = y * columns + x;
                    upperWaveTerrain[cell] = circulation::upperOrographicLaunchFraction *
                        climateatmosphere::detail::effectiveMountainWaveHeightMetres(
                        terrain[cell], backgroundWind, stability);
                }
            }
        config.upperOrographicHeightMetres = climateatmosphere::upperOrographicHeightForcing(columns, rows,
            upperWaveTerrain, zonal.upperEastWindMps, stability, hydrostatic * static_cast<float>(globaltemperature + 273.15),
            config.upperDampingTimeSeconds, config);
        constexpr std::array<float, CLIMATESEASONCOUNT> days = {15.0f, 105.0f, 196.0f, 288.0f};
        const auto solar = climateenergy::solarForcing(days[season], world.tilt(), world.eccentricity(), world.perihelion());
        const float declination = solar.declinationRadians, distance = solar.distanceFactor;
        vector<float> lowerheating(columns * rows), upperheating(columns * rows);
        for (int y = 0; y < rows; ++y)
            for (int x = 0; x < columns; ++x)
            {
                const int cell = y * columns + x;
                const float albedo = static_cast<float>(temperature[cell] < -2.0f ? energybalance::snowAlbedo :
                    energybalance::oceanAlbedo + land[cell] * (energybalance::landAlbedo - energybalance::oceanAlbedo));
                climateatmosphere::ColumnHeatingInput initial;
                initial.incomingSolarWm2 = climateenergy::dailyMeanInsolationWm2(
                    static_cast<float>(grid.latitudeCentresRadians[y] * 180.0 / 3.141592653589793), declination, distance);
                initial.surfaceAlbedo = albedo;
                initial.surfaceTemperatureK = temperature[cell] + 273.15;
                initial.airTemperatureK = {initial.surfaceTemperatureK,
                    initial.surfaceTemperatureK - moistureadvection::freeTroposphereEnvironmentalLapseC};
                initial.longwaveOpticalDepth = {circulation::lowerLongwaveOpticalDepth, circulation::upperLongwaveOpticalDepth};
                initial.shortwaveOpticalDepth = {circulation::lowerShortwaveOpticalDepth, circulation::upperShortwaveOpticalDepth};
                const auto heat = climateatmosphere::diagnoseColumnHeating(initial);
                lowerheating[cell] = static_cast<float>(heat.totalWm2[0]);
                upperheating[cell] = static_cast<float>(heat.totalWm2[1]);
            }
        if (climate.hasProcesses && climate.processes[season].durationSeconds > 0.0)
        {
            const auto& fields = climate.processes[season];
            for (int layer = 0; layer < 2; ++layer)
            {
                vector<float> net(fields.columns * fields.rows);
                for (std::size_t cell = 0; cell < net.size(); ++cell)
                    net[cell] = fields.radiativeHeatingWm2[layer][cell] +
                        circulation::stationaryLatentHeatProjectionFraction * fields.latentHeatingWm2[layer][cell] +
                        (layer == 0 ? fields.sensibleHeatingWm2[cell] : 0.0f);
                (layer == 0 ? lowerheating : upperheating) = climategrid::remapField(
                    fields.columns, fields.rows, net, columns, rows);
            }
        }
        // Q / (cp * dp/g) is a temperature tendency. Damping time and the
        // hypsometric response convert it to each mode's equilibrium forcing.
        const float columnmass = circulation::surfaceReferencePressurePa / config.gravityMetresPerSecondSquared;
        config.surfaceHeatingPressureResponseHpaPerWm2 = -pressurepermetre * hydrostatic * config.surfaceDampingTimeSeconds /
            (1004.0f * columnmass * moistureadvection::boundaryLayerCapacityFraction);
        config.upperHeatingHeightResponseMetresPerWm2 = hydrostatic * config.upperDampingTimeSeconds /
            (1004.0f * columnmass * moistureadvection::freeTroposphereCapacityFraction);
        config.upperStationaryHeatingWm2 = upperheating;
        auto result = climateatmosphere::solveModeSeparatedCirculation(columns, rows, zonalpressure,
            lowerheating, surfaceforcing, config);
        // The fixed 5 m/s pressure linearization overestimates stress in calm
        // forest cells. Update its positive drag rate from the diagnosed flow;
        // retain the full receiver roughness and the existing mass-flux closure.
        if (!climate.transientWeather && !config.boundaryLayer.enabled)
            for (int iteration = 0; iteration < circulation::surfacePressureDragIterations; ++iteration)
            {
                if (!result.surfaceStationarySolver.converged) break;
                double maximumChange = 0.0;
                for (int cell = 0; cell < columns * rows; ++cell)
                {
                    const double speed = std::max(circulation::surfacePressureDragMinimumWindMps,
                        std::hypot(result.surfaceEastWindMps[cell], result.surfaceSouthWindMps[cell]));
                    const double target = config.surfaceDragCoefficients[cell] * speed /
                        config.surfaceBoundaryLayerDepthMetres;
                    const double old = 1.0 / config.surfaceDragTimesSeconds[cell];
                    const double updated = old + circulation::surfacePressureDragRelaxation * (target - old);
                    maximumChange = std::max(maximumChange, std::abs(target - old) / std::max(target, old));
                    config.surfaceDragTimesSeconds[cell] = static_cast<float>(1.0 / updated);
                }
                result = climateatmosphere::solveModeSeparatedCirculation(columns, rows, zonalpressure,
                    lowerheating, surfaceforcing, config);
                std::cout << "Surface nonlinear drag season=" << season << " iteration=" << iteration
                    << " maximum_relative_rate_change=" << maximumChange << '\n';
            }
        climate.atmosphereAccepted[season] = result.surfaceStationarySolver.converged && result.upperStationarySolver.converged;
        if (climate.transientWeather && circulation::enablePrognosticClimateAtmosphere && config.enabled.stationary &&
            result.surfaceStationarySolver.converged && result.upperStationarySolver.converged)
        {
            climateweather::ShallowWaterConfig dynamics;
            dynamics.layerCount = 2;
            dynamics.rotationDirection = config.rotationDirection;
            dynamics.gravityMetresPerSecondSquared = config.gravityMetresPerSecondSquared;
            dynamics.lowerMeanDepthMetres = surfaceparameters.equivalentDepthMetres;
            dynamics.upperMeanDepthMetres = upperparameters.equivalentDepthMetres;
            dynamics.lowerDragTimeSeconds = config.surfaceDragTimeSeconds;
            dynamics.upperDragTimeSeconds = config.upperDragTimeSeconds;
            dynamics.heightRelaxationTimeSeconds = config.surfaceDampingTimeSeconds;
            auto state = climateweather::makeState(columns, rows, 2, static_cast<std::uint64_t>(world.seed()) + season + 1);
            climateweather::ShallowWaterForcing forcing;
            forcing.equilibriumHeightMetres = {zero, zero};
            forcing.backgroundEastWindMps = {zonal.surfaceEastWindMps, zonal.upperEastWindMps};
            forcing.backgroundSouthWindMps = {zonal.surfaceSouthWindMps, zonal.upperSouthWindMps};
            for (int cell = 0; cell < columns * rows; ++cell)
            {
                forcing.equilibriumHeightMetres[0][cell] = result.surfaceStationarySolver.equilibriumPressureAnomalyHpa[cell] / pressurepermetre;
                forcing.equilibriumHeightMetres[1][cell] = result.upperStationarySolver.equilibriumPressureAnomalyHpa[cell] / pressurepermetre;
            }
            bool stable = true;
            for (int step = 0; step < circulation::prognosticClimateSpinupSteps && stable; ++step)
            {
                const auto diagnostics = climateweather::advance(state, dynamics, forcing, circulation::prognosticClimateStepSeconds);
                stable = diagnostics.finite && diagnostics.bounded;
            }
            const auto samples = stable ? climateweather::generateWeatherSequence(state, dynamics, forcing,
                circulation::prognosticClimateSamples, circulation::prognosticClimateStepSeconds) : vector<climateweather::ShallowWaterState>{};
            if (samples.size() == static_cast<std::size_t>(circulation::prognosticClimateSamples))
            {
                const auto surface = climateweather::calculateStatistics(samples, 0);
                const auto upper = climateweather::calculateStatistics(samples, 1);
                for (int y = 0; y < rows; ++y)
                {
                    double surfacemean = 0.0, uppermean = 0.0;
                    for (int x = 0; x < columns; ++x)
                    {
                        surfacemean += surface.meanHeightAnomalyMetres[y * columns + x] / columns;
                        uppermean += upper.meanHeightAnomalyMetres[y * columns + x] / columns;
                    }
                    for (int x = 0; x < columns; ++x)
                    {
                        const int cell = y * columns + x;
                        if (config.enabled.surface)
                            result.surfacePressureAnomalyHpa[cell] += circulation::prognosticClimateBlend *
                                (pressurepermetre * (surface.meanHeightAnomalyMetres[cell] - static_cast<float>(surfacemean)) -
                                    result.surfaceStationarySolver.pressureAnomalyHpa[cell]);
                        if (config.enabled.upper)
                            result.upperHeightAnomalyMetres[cell] += circulation::prognosticClimateBlend *
                                (upper.meanHeightAnomalyMetres[cell] - static_cast<float>(uppermean) -
                                    result.upperStationarySolver.pressureAnomalyHpa[cell] / pressurepermetre);
                    }
                }
                climateatmosphere::diagnoseModeWinds(columns, rows, config, result);
            }
            else
            {
                const auto amplitude = [](const vector<float>& field)
                {
                    float maximum = 0.0f;
                    for (float value : field) maximum = std::max(maximum, std::abs(value));
                    return maximum;
                };
                std::cout << "Prognostic atmosphere fallback: season=" << season
                    << " samples=" << samples.size() << " stable_spinup=" << stable
                    << " surface_forcing_m=" << amplitude(forcing.equilibriumHeightMetres[0])
                    << " upper_forcing_m=" << amplitude(forcing.equilibriumHeightMetres[1])
                    << " background_upper_mps=" << amplitude(zonal.upperEastWindMps) << '\n';
            }
        }
        if (climate.mixedLayerTropics && config.enabled.surface)
        {
            // Calibrated projection into the reduced inversion mode. Applying
            // the former surface amplitude unchanged overdrives 850 hPa winds.
            constexpr float inversionOverturningProjection = 0.4f;
            config.inversionEquilibriumPressureHpa.resize(columns * rows);
            for (int y = 0; y < rows; ++y)
                for (int x = 0; x < columns; ++x)
                    config.inversionEquilibriumPressureHpa[y * columns + x] = config.enabled.zonal ?
                        inversionOverturningProjection * zonalpressure[y] : 0.0f;
            vector<float> columnheating(columns * rows);
            for (std::size_t i = 0; i < columnheating.size(); ++i) columnheating[i] = lowerheating[i] + upperheating[i];
            auto& layer = climate.boundaryLayer[season];
            layer = climateatmosphere::solveBoundaryLayer(columns, rows, temperature, columnheating, mechanicalforcing, config);
            const bool accepted = climateatmosphere::applyBoundaryLayer(layer, config, result);
            climate.atmosphereAccepted[season] = climate.atmosphereAccepted[season] && accepted;
            std::cout << "Tropical boundary layer season=" << season << " accepted=" << accepted
                << " inversion_hpa=" << config.boundaryLayer.inversionPressurePa / 100.0f
                << " overturning_projection=" << inversionOverturningProjection
                << " residual=" << layer.inversionSolver.relativeResidual
                << " iterations=" << layer.inversionSolver.iterations
                << " momentum_residual_mps2=" << result.maximumBoundaryLayerMomentumResidualMps2 << '\n';
        }
        const auto nativeSolvedPressure = result.surfacePressureAnomalyHpa;
        if (const auto* referencepressure = referenceexperiment().field(columns, "slp"))
        {
            for (int y = 0; y < rows; ++y)
                for (int x = 0; x < columns; ++x)
                    result.surfacePressureAnomalyHpa[y * columns + x] = climateio::experiment::sample(
                        referencepressure, season, season, 0.0f, x, y, result.surfacePressureAnomalyHpa[y * columns + x]);
            // Replace the diagnosed pressure, not its equilibrium forcing. This also
            // recomputes ascent from the resulting pressure fluxes before later corrections.
            climateatmosphere::diagnoseModeWinds(columns, rows, config, result);
            std::cout << "Reference surface pressure imposed at wind diagnosis season=" << season << '\n';
        }
        const auto preAdvectionEast = result.surfaceEastWindMps, preAdvectionSouth = result.surfaceSouthWindMps;
        if (circulation::surfaceMomentumAdvectionIterations > 0 && !climate.transientWeather &&
            !config.boundaryLayer.enabled && config.enabled.surface)
        {
            // Recover the existing effective forcing, including stationary
            // exchange, before adding bulk horizontal momentum transport.
            // Pressure remains the preceding inversion; this is a bounded
            // stationary correction, not a prognostic pressure closure.
            vector<float> forceEast(columns * rows), forceNorth(columns * rows);
            for (int y = 0; y < rows; ++y)
                for (int x = 0; x < columns; ++x)
                {
                    const int cell = y * columns + x;
                    const auto force = climateatmosphere::detail::surfaceMomentumForce(
                        {preAdvectionEast[cell], preAdvectionSouth[cell]}, config.surfaceDragCoefficients[cell],
                        latitudeforrow(y, rows - 1), config);
                    forceEast[cell] = force.eastMps2; forceNorth[cell] = force.northMps2;
                }
            const auto advected = climateatmosphere::detail::advectSurfaceMomentum(grid, forceEast, forceNorth,
                config.surfaceDragCoefficients, preAdvectionEast, preAdvectionSouth, config,
                circulation::surfaceMomentumAdvectionIterations);
            if (advected.east.size() == preAdvectionEast.size())
            {
                vector<float> du(columns * rows), dv(columns * rows);
                for (int cell = 0; cell < columns * rows; ++cell)
                {
                    du[cell] = advected.east[cell] - preAdvectionEast[cell];
                    dv[cell] = advected.south[cell] - preAdvectionSouth[cell];
                }
                const auto ascentchange = climateatmosphere::detail::surfaceAscentCorrectionHpaPerDay(
                    grid, du, dv, 0.5f * config.surfaceEquivalentPressureDepthHpa);
                for (std::size_t cell = 0; cell < ascentchange.size(); ++cell) result.ascentHpaPerDay[cell] += ascentchange[cell];
                result.surfaceEastWindMps = advected.east; result.surfaceSouthWindMps = advected.south;
                std::cout << "Surface momentum advection season=" << season << " iterations=" << advected.iterations
                    << " maximum_change_mps=" << advected.maximumChangeMps
                    << " relative_momentum_residual=" << advected.relativeMomentumResidual
                    << " forcing_work_wm2=" << advected.forcingWorkWm2 << " drag_work_wm2=" << advected.dragWorkWm2
                    << " advection_work_wm2=" << advected.advectionWorkWm2
                    << " upwind_mixing_work_wm2=" << advected.upwindMixingWorkWm2
                    << " ascent_integral=" << climategrid::areaWeightedIntegral(columns, rows, ascentchange) << '\n';
            }
        }
        vector<float> wind850East(columns * rows, 0.0f), wind850South(columns * rows, 0.0f),
            pressure850(columns * rows, 0.0f), available850(columns * rows, 0.0f), thermal850(columns * rows, 0.0f);
        const auto pre850SurfaceEast = result.surfaceEastWindMps, pre850SurfaceSouth = result.surfaceSouthWindMps;
        auto preExchange850East = wind850East, preExchange850South = wind850South;
        auto native850East = wind850East, native850South = wind850South;
        if (circulation::enableAuxiliary850Momentum && !climate.transientWeather &&
            !config.boundaryLayer.enabled && config.enabled.surface)
        {
            // Diagnose a distinct lower-free-tropospheric mode without applying
            // the alternative surface closure or replacing the 500 hPa mode.
            auto lowConfig = config;
            // Keep inversion mass adjustment independent of the surface mode.
            // Sharing its equivalent depth made a shallow surface-pressure
            // experiment also amplify the weakly damped 850 hPa circulation.
            const auto lowParameters = climateatmosphere::diagnoseStationaryParameters(stability,
                circulation::auxiliary850ModeDepthMetres, config.gravityMetresPerSecondSquared,
                config.planetRadiusMetres, config.rotationRatePerSecond, columns, rows,
                circulation::stationaryWaveNondimensionalDamping, circulation::diagnosedForcingScaleMetres);
            lowConfig.surfaceEquivalentPressureDepthHpa = pressurepermetre * lowParameters.equivalentDepthMetres;
            lowConfig.surfaceDampingTimeSeconds = lowParameters.dampingTimeSeconds;
            lowConfig.boundaryLayer.enabled = true;
            lowConfig.stationaryCaptureDirectory.clear();
            lowConfig.inversionEquilibriumPressureHpa.resize(columns * rows);
            vector<float> columnHeating(columns * rows);
            for (int cell = 0; cell < columns * rows; ++cell)
            {
                // Reduced-mode sensitivity only: actual column energy remains unchanged.
                columnHeating[cell] = circulation::auxiliary850HeatingResponseFraction *
                    (lowerheating[cell] + upperheating[cell]);
                // Stored regional warmth is distinct from instantaneous net
                // heating. Approximate an anomaly decaying toward 500 hPa;
                // use the same pressure-only elevation proxy as the surface.
                thermal850[cell] = climateatmosphere::thermalModePressureAnomalyHpa(
                    circulation::auxiliary850ThermalProfileFraction * (thermaltemperature[cell] - thermalrowtemperature[cell / columns]),
                    lowConfig.boundaryLayer.inversionDensityKgM3, lowConfig.boundaryLayer.inversionPressurePa, 50000.0f);
                lowConfig.inversionEquilibriumPressureHpa[cell] = circulation::auxiliary850OverturningProjection *
                    (circulation::auxiliary850ZonalPressureFraction * zonalpressure[cell / columns] +
                        regionalPressure[cell] + continentalPressure[cell]) + thermal850[cell];
            }
            const auto low = climateatmosphere::solveBoundaryLayer(columns, rows, temperature,
                columnHeating, mechanicalforcing, lowConfig);
            if (low.inversionSolver.converged)
            {
                pressure850 = low.inversionPressureHpa;
                for (int cell = 0; cell < columns * rows; ++cell)
                {
                    const float height850 = 287.05f * (temperature[cell] + 273.15f) /
                        config.gravityMetresPerSecondSquared * std::log(100000.0f / lowConfig.boundaryLayer.inversionPressurePa);
                    // Approximate pressure-surface availability, not cut-cell terrain dynamics.
                    available850[cell] = terrain[cell] < height850 ? 1.0f : 0.0f;
                    wind850East[cell] = available850[cell] * low.inversionEastWindMps[cell];
                    wind850South[cell] = available850[cell] * low.inversionSouthWindMps[cell];
                }
                native850East = wind850East; native850South = wind850South;
                const auto* reference850east = referenceexperiment().field(columns, "u850");
                const auto* reference850north = referenceexperiment().field(columns, "v850");
                if (bool(reference850east) != bool(reference850north))
                    throw std::runtime_error("Reference 850-hPa wind requires both components");
                if (reference850east)
                {
                    for (int y = 0; y < rows; ++y)
                        for (int x = 0; x < columns; ++x)
                        {
                            const int cell = y * columns + x;
                            if (available850[cell] == 0.0f) continue;
                            wind850East[cell] = climateio::experiment::sample(reference850east,
                                season, season, 0.0f, x, y, wind850East[cell]);
                            wind850South[cell] = -climateio::experiment::sample(reference850north,
                                season, season, 0.0f, x, y, -wind850South[cell]);
                        }
                    std::cout << "Reference 850-hPa driver imposed before two-way exchange season=" << season << '\n';
                }
                preExchange850East = wind850East; preExchange850South = wind850South;
                const float surfaceMass = config.airDensityKgM3 * config.surfaceBoundaryLayerDepthMetres;
                // The remainder of the nominal surface-to-850 pressure mass is
                // a bulk low-level reservoir represented by the 850 hPa wind.
                // This is not an additional prognostic moisture-transport layer.
                const float lowMass = (100000.0f - lowConfig.boundaryLayer.inversionPressurePa) /
                    config.gravityMetresPerSecondSquared - surfaceMass;
                if (circulation::surface850ExchangeRatePerSecond > 0.0f && lowMass > 0.0f)
                {
                    vector<float> du(columns * rows, 0.0f), dv(columns * rows, 0.0f);
                    double area = 0.0, loss = 0.0, surfaceWork = 0.0, lowWork = 0.0, maximumResidual = 0.0;
                    int failed = 0;
                    const float rayleigh = 1.0f / lowConfig.boundaryLayer.inversionDragTimeSeconds;
                    for (int y = 0; y < rows; ++y)
                        for (int x = 0; x < columns; ++x)
                        {
                            const int cell = y * columns + x;
                            const double cellArea = grid.cellAreasSquareMetres[y];
                            area += cellArea;
                            const float latitude = latitudeforrow(y, rows - 1);
                            const float exchange = circulation::surface850ExchangeRatePerSecond * available850[cell] *
                                climateatmosphere::boundaryLayerWeight(latitude, lowConfig.boundaryLayer);
                            if (exchange == 0.0f) continue;
                            const float f = climateatmosphere::coriolisParameterPerSecond(latitude,
                                config.rotationRatePerSecond, config.rotationDirection);
                            const auto surfaceForce = climateatmosphere::detail::surfaceMomentumForce(
                                {pre850SurfaceEast[cell], pre850SurfaceSouth[cell]}, config.surfaceDragCoefficients[cell], latitude, config);
                            const climateatmosphere::detail::SurfaceMomentumForce lowForce{
                                rayleigh * wind850East[cell] + f * wind850South[cell],
                                f * wind850East[cell] - rayleigh * wind850South[cell]};
                            const auto coupled = climateatmosphere::detail::coupleSurface850(surfaceForce, lowForce, f,
                                config.surfaceDragCoefficients[cell], config.surfaceBoundaryLayerDepthMetres,
                                config.airDensityKgM3, lowMass, rayleigh, exchange);
                            if (!coupled.valid) { ++failed; continue; }
                            result.surfaceEastWindMps[cell] = coupled.surface.eastMetresPerSecond;
                            result.surfaceSouthWindMps[cell] = coupled.surface.southMetresPerSecond;
                            wind850East[cell] = coupled.lowLevel.eastMetresPerSecond;
                            wind850South[cell] = coupled.lowLevel.southMetresPerSecond;
                            du[cell] = result.surfaceEastWindMps[cell] - pre850SurfaceEast[cell];
                            dv[cell] = result.surfaceSouthWindMps[cell] - pre850SurfaceSouth[cell];
                            loss += cellArea * coupled.mixingLossWm2;
                            surfaceWork += cellArea * (coupled.stressEastNm2 * result.surfaceEastWindMps[cell] -
                                coupled.stressNorthNm2 * result.surfaceSouthWindMps[cell]);
                            lowWork -= cellArea * (coupled.stressEastNm2 * wind850East[cell] -
                                coupled.stressNorthNm2 * wind850South[cell]);
                            maximumResidual = std::max(maximumResidual, coupled.maximumResidualMps2);
                        }
                    const auto ascentchange = climateatmosphere::detail::surfaceAscentCorrectionHpaPerDay(
                        grid, du, dv, 0.5f * config.surfaceEquivalentPressureDepthHpa);
                    for (std::size_t cell = 0; cell < ascentchange.size(); ++cell) result.ascentHpaPerDay[cell] += ascentchange[cell];
                    std::cout << "Surface 850 exchange season=" << season << " failed_cells=" << failed
                        << " low_mass_kgm2=" << lowMass << " maximum_momentum_residual_mps2=" << maximumResidual
                        << " mixing_loss_wm2=" << loss / area << " surface_work_wm2=" << surfaceWork / area
                        << " low_work_wm2=" << lowWork / area
                        << " work_balance_wm2=" << (loss + surfaceWork + lowWork) / area
                        << " ascent_integral=" << climategrid::areaWeightedIntegral(columns, rows, ascentchange) << '\n';
                }
            }
            std::cout << "Auxiliary 850 hPa season=" << season << " converged=" << low.inversionSolver.converged
                << " residual=" << low.inversionSolver.relativeResidual
                << " overturning_projection=" << circulation::auxiliary850OverturningProjection
                << " exchange_rate_per_s=" << circulation::surface850ExchangeRatePerSecond << '\n';
        }
        if (!climate.genericLevelPressuresHpa.empty() && !climate.transientWeather && config.enabled.surface)
        {
            // Experimental: every generic level uses the same rules, so that one level
            // class serves any pressure between the boundary layer and 500 hPa:
            //  - a troposphere-deep mode (the 850 hPa mode depth) and its diagnosed damping;
            //  - the free-tropospheric drag of the 850 hPa mode;
            //  - ideal-gas density at the level from the mean temperature and lapse rate;
            //  - the surface pressure pattern projected with ln(p / p_r) / ln(p_s / p_r),
            //    reversing at the heat-low depth p_r, plus the solver's own heating term.
            // Diagnostic only: nothing below reads these winds.
            const float reversalPa = 100.0f * climate.genericLevelReversalHpa;
            climate.winds.genericlevels.resize(climate.genericLevelPressuresHpa.size());
            const auto levelParameters = climateatmosphere::diagnoseStationaryParameters(stability,
                circulation::auxiliary850ModeDepthMetres, config.gravityMetresPerSecondSquared,
                config.planetRadiusMetres, config.rotationRatePerSecond, columns, rows,
                circulation::stationaryWaveNondimensionalDamping, circulation::diagnosedForcingScaleMetres);
            vector<float> levelHeating(columns * rows);
            for (int cell = 0; cell < columns * rows; ++cell)
                levelHeating[cell] = lowerheating[cell] + upperheating[cell];
            double meanTemperature = 0.0, meanWeight = 0.0;
            for (int y = 0; y < rows; ++y)
                for (int x = 0; x < columns; ++x)
                {
                    meanTemperature += grid.cellAreasSquareMetres[y] * temperature[y * columns + x];
                    meanWeight += grid.cellAreasSquareMetres[y];
                }
            meanTemperature /= std::max(1.0, meanWeight);
            for (std::size_t index = 0; index < climate.genericLevelPressuresHpa.size(); ++index)
            {
                auto& stored = climate.winds.genericlevels[index];
                const float levelPa = 100.0f * climate.genericLevelPressuresHpa[index];
                stored.pressureHpa = climate.genericLevelPressuresHpa[index];
                auto levelConfig = config;
                levelConfig.surfaceEquivalentPressureDepthHpa = pressurepermetre * levelParameters.equivalentDepthMetres;
                levelConfig.surfaceDampingTimeSeconds = levelParameters.dampingTimeSeconds;
                levelConfig.boundaryLayer.enabled = true;
                levelConfig.boundaryLayer.inversionPressurePa = levelPa;
                const double levelHeight = 287.05 * (meanTemperature + 273.15) / config.gravityMetresPerSecondSquared *
                    std::log(100000.0 / levelPa);
                levelConfig.boundaryLayer.inversionDensityKgM3 = static_cast<float>(levelPa /
                    (287.05 * (meanTemperature + 273.15 - 0.0065 * levelHeight)));
                levelConfig.stationaryCaptureDirectory.clear();
                const float projection = std::log(levelPa / reversalPa) / std::log(100000.0f / reversalPa);
                levelConfig.inversionEquilibriumPressureHpa.assign(columns * rows, 0.0f);
                for (int cell = 0; cell < columns * rows; ++cell)
                    levelConfig.inversionEquilibriumPressureHpa[cell] = projection *
                        (zonalpressure[cell / columns] + regionalPressure[cell] + continentalPressure[cell]);
                const auto level = climateatmosphere::solveBoundaryLayer(columns, rows, temperature,
                    levelHeating, mechanicalforcing, levelConfig);
                std::cout << "Generic level " << stored.pressureHpa << " hPa season=" << season
                    << " converged=" << level.inversionSolver.converged << " projection=" << projection
                    << " density=" << levelConfig.boundaryLayer.inversionDensityKgM3 << '\n';
                if (!level.inversionSolver.converged)
                    continue;
                vector<float> east(columns * rows), south(columns * rows), available(columns * rows);
                for (int cell = 0; cell < columns * rows; ++cell)
                {
                    const float height = 287.05f * (temperature[cell] + 273.15f) / config.gravityMetresPerSecondSquared *
                        std::log(100000.0f / levelPa);
                    available[cell] = terrain[cell] < height ? 1.0f : 0.0f;
                    east[cell] = available[cell] * level.inversionEastWindMps[cell];
                    south[cell] = available[cell] * level.inversionSouthWindMps[cell];
                }
                stored.windu[season] = output(east);
                stored.windv[season] = output(south);
                stored.available[season] = output(available);
                stored.solved[season] = true;
            }
        }
        vector<float> carrierEast, carrierSouth;
        const auto preDrainageEast = result.surfaceEastWindMps, preDrainageSouth = result.surfaceSouthWindMps;
        auto diagnosticSkin = temperature, diagnosticSnow = permanentice;
        vector<float> diagnosticDeficit(columns * rows, 0.0f), diagnosticForceEast(columns * rows, 0.0f),
            diagnosticForceNorth(columns * rows, 0.0f);
        vector<float> diagnosticRadiativeDeficit(columns * rows, 0.0f), diagnosticRadiation(columns * rows, 0.0f);
        const float drainageMean = atmosphere::katabaticShallowTenMetreDiagnostic ?
            climateatmosphere::detail::coldJetTransportFraction(
                atmosphere::katabaticShallowDepthMetres, config.surfaceBoundaryLayerDepthMetres) :
            atmosphere::separateKatabaticTransport ?
            climateatmosphere::detail::katabaticLayerMeanFraction(
                atmosphere::katabaticVelocityDecayHeightMetres, config.surfaceBoundaryLayerDepthMetres) : 1.0f;
        climate.winds.hasDrainageProfile[season] = false;
        if (atmosphere::katabaticDeficitFraction > 0.0f && !climate.transientWeather &&
            !config.boundaryLayer.enabled && config.enabled.surface && climate.hasProcesses)
        {
            const auto started = std::chrono::steady_clock::now();
            const auto& fields = climate.processes[season];
            const auto skin = climategrid::remapField(fields.columns, fields.rows,
                fields.surfaceTemperatureC, columns, rows);
            const auto snow = climategrid::remapField(fields.columns, fields.rows,
                fields.landSnowCoverFraction, columns, rows);
            const auto surfaceRadiation = climategrid::remapField(fields.columns, fields.rows,
                fields.surfaceRadiativeHeatingWm2, columns, rows);
            diagnosticRadiation = surfaceRadiation;
            diagnosticSkin = skin;
            for (int cell = 0; cell < columns * rows; ++cell) diagnosticSnow[cell] = std::max(snow[cell], permanentice[cell]);
            vector<float> deltau(columns * rows), deltav(columns * rows);
            if (atmosphere::separateKatabaticTransport)
            {
                carrierEast = result.surfaceEastWindMps;
                carrierSouth = result.surfaceSouthWindMps;
            }
            double area = 0.0, activeArea = 0.0, deficitArea = 0.0, speedChangeArea = 0.0, powerArea = 0.0;
            float maximumForce = 0.0f, maximumDeficit = 0.0f;
            for (int y = 0; y < rows; ++y)
                for (int x = 0; x < columns; ++x)
                {
                    const int cell = y * columns + x;
                    const auto skinDeficitFunction = atmosphere::katabaticSnowPatchTemperature ?
                        climateatmosphere::detail::snowPatchKatabaticDeficitK : climateatmosphere::detail::katabaticDeficitK;
                    const float skinDeficit = skinDeficitFunction(
                        temperature[cell], skin[cell], std::max(snow[cell], permanentice[cell]), atmosphere::katabaticDeficitFraction,
                        atmosphere::katabaticMaximumDeficitK);
                    const float radiativeDeficit = climateatmosphere::detail::radiativeKatabaticDeficitK(
                        surfaceRadiation[cell], skin[cell], std::max(snow[cell], permanentice[cell]), config.airDensityKgM3,
                        1004.0f, atmosphere::katabaticShallowTenMetreDiagnostic ? atmosphere::katabaticShallowDepthMetres :
                            atmosphere::katabaticVelocityDecayHeightMetres, atmosphere::katabaticColdLayerRenewalSeconds,
                        atmosphere::katabaticRadiativeCoolingFraction, atmosphere::katabaticMaximumDeficitK);
                    // These are alternative estimates of the same cold anomaly.
                    const float deficit = std::max(skinDeficit, radiativeDeficit);
                    diagnosticRadiativeDeficit[cell] = radiativeDeficit;
                    const auto katabatic = climateatmosphere::detail::katabaticForce(
                        slopes[cell].east, slopes[cell].north, temperature[cell] + 273.15f,
                        deficit, config.gravityMetresPerSecondSquared, atmosphere::katabaticMaximumSlope);
                    const float magnitude = std::hypot(katabatic.eastMps2, katabatic.northMps2);
                    diagnosticDeficit[cell] = deficit;
                    diagnosticForceEast[cell] = katabatic.eastMps2;
                    diagnosticForceNorth[cell] = katabatic.northMps2;
                    area += grid.cellAreasSquareMetres[y];
                    if (magnitude == 0.0f) continue;
                    const climateatmosphere::HorizontalWind original{result.surfaceEastWindMps[cell], result.surfaceSouthWindMps[cell]};
                    const auto force = climateatmosphere::detail::surfaceMomentumForce(original,
                        config.surfaceDragCoefficients[cell], latitudeforrow(y, rows - 1), config);
                    const auto wind = atmosphere::katabaticShallowTenMetreDiagnostic ?
                        climateatmosphere::detail::coldJetTenMetreWind(original, katabatic, latitudeforrow(y, rows - 1),
                            config.surfaceDragCoefficients[cell] * atmosphere::katabaticTenMetreStressFraction,
                            atmosphere::katabaticShallowDepthMetres, config.rotationRatePerSecond, config.rotationDirection) :
                        climateatmosphere::steadyQuadraticDragCoriolisWind(
                        force.eastMps2 + katabatic.eastMps2, force.northMps2 + katabatic.northMps2,
                        latitudeforrow(y, rows - 1), config.surfaceDragCoefficients[cell],
                        config.surfaceBoundaryLayerDepthMetres, config.rotationRatePerSecond, config.rotationDirection);
                    deltau[cell] = wind.eastMetresPerSecond - original.eastMetresPerSecond;
                    deltav[cell] = wind.southMetresPerSecond - original.southMetresPerSecond;
                    result.surfaceEastWindMps[cell] = wind.eastMetresPerSecond;
                    result.surfaceSouthWindMps[cell] = wind.southMetresPerSecond;
                    activeArea += grid.cellAreasSquareMetres[y];
                    deficitArea += grid.cellAreasSquareMetres[y] * deficit;
                    speedChangeArea += grid.cellAreasSquareMetres[y] *
                        (std::hypot(wind.eastMetresPerSecond, wind.southMetresPerSecond) -
                         std::hypot(original.eastMetresPerSecond, original.southMetresPerSecond));
                    const float powerDepth = atmosphere::katabaticShallowTenMetreDiagnostic ?
                        atmosphere::katabaticShallowDepthMetres : config.surfaceBoundaryLayerDepthMetres;
                    const float powerMean = atmosphere::katabaticShallowTenMetreDiagnostic ?
                        static_cast<float>(climateatmosphere::detail::coldJetMeanToTenMetres(powerDepth)) : 1.0f;
                    powerArea += grid.cellAreasSquareMetres[y] * config.airDensityKgM3 * powerDepth *
                        (katabatic.eastMps2 * (original.eastMetresPerSecond + powerMean * deltau[cell]) -
                         katabatic.northMps2 * (original.southMetresPerSecond + powerMean * deltav[cell]));
                    maximumForce = std::max(maximumForce, magnitude);
                    maximumDeficit = std::max(maximumDeficit, deficit);
                }
            for (float& value : deltau) value *= drainageMean;
            for (float& value : deltav) value *= drainageMean;
            const auto ascentchange = climateatmosphere::detail::surfaceAscentCorrectionHpaPerDay(
                grid, deltau, deltav, 0.5f * config.surfaceEquivalentPressureDepthHpa);
            for (std::size_t cell = 0; cell < ascentchange.size(); ++cell) result.ascentHpaPerDay[cell] += ascentchange[cell];
            std::cout << "Katabatic drainage season=" << season << " active_area_fraction=" << activeArea / area
                << " mean_deficit_k=" << deficitArea / std::max(1.0, activeArea) << " max_deficit_k=" << maximumDeficit
                << " mean_speed_change_mps=" << speedChangeArea / std::max(1.0, activeArea)
                << " max_force_mps2=" << maximumForce << " mean_power_wm2=" << powerArea / area
                << " transport_increment_fraction=" << drainageMean
                << " shallow_10m=" << atmosphere::katabaticShallowTenMetreDiagnostic
                << " cold_depth_m=" << (atmosphere::katabaticShallowTenMetreDiagnostic ?
                    atmosphere::katabaticShallowDepthMetres : atmosphere::katabaticVelocityDecayHeightMetres)
                << " ascent_integral=" << climategrid::areaWeightedIntegral(columns, rows, ascentchange)
                << " seconds=" << std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count() << '\n';
        }
        if (!config.stationaryCaptureDirectory.empty())
        {
            vector<float> thermal(surfaceforcing.size()), diabatic(surfaceforcing.size()), zonalfield(surfaceforcing.size());
            auto radiativePressure = lowerheating;
            vector<float> latentPressure(surfaceforcing.size(), 0.0f), sensiblePressure(surfaceforcing.size(), 0.0f);
            if (climate.hasProcesses && climate.processes[season].durationSeconds > 0.0)
            {
                const auto& fields = climate.processes[season];
                radiativePressure = climategrid::remapField(fields.columns, fields.rows, fields.radiativeHeatingWm2[0], columns, rows);
                latentPressure = climategrid::remapField(fields.columns, fields.rows, fields.latentHeatingWm2[0], columns, rows);
                sensiblePressure = climategrid::remapField(fields.columns, fields.rows, fields.sensibleHeatingWm2, columns, rows);
            }
            vector<float> linearizationSpeed(surfaceforcing.size());
            vector<float> upperZonalHeight(surfaceforcing.size()), upperHeatingHeight(surfaceforcing.size()),
                upperSurfaceHeight(surfaceforcing.size());
            for (int cell = 0; cell < columns * rows; ++cell)
            {
                linearizationSpeed[cell] = config.surfaceBoundaryLayerDepthMetres /
                    (config.surfaceDragCoefficients[cell] * config.surfaceDragTimesSeconds[cell]);
                thermal[cell] = surfaceforcing[cell] - mechanicalforcing[cell] - regionalPressure[cell] - continentalPressure[cell];
                diabatic[cell] = lowerheating[cell] * config.surfaceHeatingPressureResponseHpaPerWm2;
                radiativePressure[cell] *= config.surfaceHeatingPressureResponseHpaPerWm2;
                latentPressure[cell] *= circulation::stationaryLatentHeatProjectionFraction * config.surfaceHeatingPressureResponseHpaPerWm2;
                sensiblePressure[cell] *= config.surfaceHeatingPressureResponseHpaPerWm2;
                zonalfield[cell] = zonalpressure[cell / columns];
                upperZonalHeight[cell] = config.zonalUpperHeightMetres[cell / columns];
                upperHeatingHeight[cell] = upperheating[cell] * config.upperHeatingHeightResponseMetresPerWm2;
                upperSurfaceHeight[cell] = (diabatic[cell] + surfaceforcing[cell]) *
                    config.surfaceToUpperCoupling / pressurepermetre;
            }
            climateatmosphere::detail::captureCirculation(config.stationaryCaptureDirectory, columns, rows, season,
                {{"temperature_c", temperature}, {"land_fraction", land}, {"terrain_m", terrain},
                 {"pressure_linearization_speed_mps", linearizationSpeed},
                 {"permanent_ice_fraction", permanentice}, {"drag_coefficient", config.surfaceDragCoefficients},
                 {"persistent_snow_fraction", persistentsnow},
                 {"zonal_hpa", zonalfield}, {"mechanical_hpa", mechanicalforcing}, {"thermal_hpa", thermal},
                 {"continental_eddy_hpa", continentalPressure},
                 {"regional_trough_hpa", regionalPressure}, {"regional_trough_latitude", regionalTrough},
                 {"diabatic_hpa", diabatic}, {"filtered_equilibrium_hpa", result.surfaceStationarySolver.equilibriumPressureAnomalyHpa},
                 {"diabatic_radiative_hpa", radiativePressure}, {"diabatic_latent_hpa", latentPressure},
                 {"diabatic_sensible_hpa", sensiblePressure},
                 {"upper_zonal_height_m", upperZonalHeight}, {"upper_heating_height_m", upperHeatingHeight},
                 {"upper_orographic_height_m", config.upperOrographicHeightMetres},
                 {"upper_wave_terrain_m", upperWaveTerrain},
                 {"upper_surface_coupling_height_m", upperSurfaceHeight},
                 {"upper_filtered_equilibrium_hpa", result.upperStationarySolver.equilibriumPressureAnomalyHpa},
                 {"upper_solved_height_m", result.upperHeightAnomalyMetres},
                 {"upper_east_mps", result.upperEastWindMps}, {"upper_south_mps", result.upperSouthWindMps},
                 {"solved_hpa", result.surfacePressureAnomalyHpa}, {"east_mps", result.surfaceEastWindMps},
                 {"native_solved_hpa", nativeSolvedPressure},
                 {"south_mps", result.surfaceSouthWindMps}, {"thermal_input_c", thermaltemperature},
                 {"skin_c", diagnosticSkin}, {"snow_ice_fraction", diagnosticSnow}, {"drainage_deficit_k", diagnosticDeficit},
                 {"drainage_east_mps2", diagnosticForceEast}, {"drainage_north_mps2", diagnosticForceNorth},
                 {"pre_drainage_east_mps", preDrainageEast}, {"pre_drainage_south_mps", preDrainageSouth},
                 {"pre_advection_east_mps", preAdvectionEast}, {"pre_advection_south_mps", preAdvectionSouth},
                 {"low_level_850_east_mps", wind850East}, {"low_level_850_south_mps", wind850South},
                 {"low_level_850_pressure_hpa", pressure850}, {"low_level_850_available", available850},
                 {"low_level_850_thermal_hpa", thermal850},
                 {"pre_850_surface_east_mps", pre850SurfaceEast}, {"pre_850_surface_south_mps", pre850SurfaceSouth},
                 {"pre_exchange_850_east_mps", preExchange850East}, {"pre_exchange_850_south_mps", preExchange850South},
                 {"native_850_east_mps", native850East}, {"native_850_south_mps", native850South},
                 {"process_fields_available", vector<float>(columns * rows, climate.hasProcesses ? 1.0f : 0.0f)},
                 {"radiative_deficit_k", diagnosticRadiativeDeficit}, {"surface_radiative_wm2", diagnosticRadiation},
                 {"drainage_transport_fraction", vector<float>(columns * rows, drainageMean)},
                 {"cold_layer_depth_m", vector<float>(columns * rows, atmosphere::katabaticShallowTenMetreDiagnostic ?
                     atmosphere::katabaticShallowDepthMetres : atmosphere::katabaticVelocityDecayHeightMetres)},
                 {"shallow_10m_enabled", vector<float>(columns * rows, atmosphere::katabaticShallowTenMetreDiagnostic ? 1.0f : 0.0f)}});
        }
        auto pressure = output(result.surfacePressureAnomalyHpa);
        auto upperheight = output(result.upperHeightAnomalyMetres);
        auto windu = output(result.surfaceEastWindMps), windv = output(result.surfaceSouthWindMps);
        auto upperu = output(result.upperEastWindMps), upperv = output(result.upperSouthWindMps);
        auto vertical = output(result.ascentHpaPerDay);
        floatgrid carrieru, carrierv;
        if (!carrierEast.empty()) { carrieru = output(carrierEast); carrierv = output(carrierSouth); }
        if (!climate.transientWeather && !config.boundaryLayer.enabled && config.enabled.surface &&
            climategrid::validGlobalGridDimensions(width + 1, height + 1))
        {
            // Interpolate the native effective acceleration, then apply local
            // drag. Averaging wind and Cd separately changes quadratic stress.
            vector<float> nativeeast(columns * rows), nativenorth(columns * rows);
            vector<float> carrierForceEast(carrierEast.size()), carrierForceNorth(carrierEast.size());
            for (int y = 0; y < rows; ++y)
                for (int x = 0; x < columns; ++x)
                {
                    const int cell = y * columns + x;
                    const auto force = climateatmosphere::detail::surfaceMomentumForce(
                        {result.surfaceEastWindMps[cell], result.surfaceSouthWindMps[cell]},
                        config.surfaceDragCoefficients[cell], latitudeforrow(y, rows - 1), config);
                    nativeeast[cell] = force.eastMps2;
                    nativenorth[cell] = force.northMps2;
                    if (!carrierEast.empty())
                    {
                        const auto carrierForce = climateatmosphere::detail::surfaceMomentumForce(
                            {carrierEast[cell], carrierSouth[cell]}, config.surfaceDragCoefficients[cell],
                            latitudeforrow(y, rows - 1), config);
                        carrierForceEast[cell] = carrierForce.eastMps2;
                        carrierForceNorth[cell] = carrierForce.northMps2;
                    }
                }
            const auto eastforce = output(nativeeast), northforce = output(nativenorth);
            const auto carrierEastForce = carrierEast.empty() ? floatgrid{} : output(carrierForceEast);
            const auto carrierNorthForce = carrierEast.empty() ? floatgrid{} : output(carrierForceNorth);
            // Reconstruct background and buoyancy separately: reinterpreting a
            // shallow 10 m wind as a 300 m force would undo its momentum closure.
            const auto coldEastForce = atmosphere::katabaticShallowTenMetreDiagnostic ? output(diagnosticForceEast) : floatgrid{};
            const auto coldNorthForce = atmosphere::katabaticShallowTenMetreDiagnostic ? output(diagnosticForceNorth) : floatgrid{};
            const auto outputgrid = climategrid::makeSphericalGrid(
                width + 1, height + 1, config.planetRadiusMetres);
            vector<float> deltau((width + 1) * (height + 1)), deltav(deltau.size());
            for (int y = 0; y <= height; ++y)
                for (int x = 0; x <= width; ++x)
                {
                    const float drag = surfacedrag(world.nom(x, y) > world.sealevel() ? 1.0f : 0.0f,
                        std::max(0.0f, world.nomfloat(x, y) - static_cast<float>(world.sealevel())),
                        std::max(worldpermanentlandicefraction(world,x,y), atmosphere::persistentSnowDragFraction *
                            persistentsnowoutput[y * (width + 1) + x]), vegetationcontrast[y * (width + 1) + x]);
                    auto corrected = climateatmosphere::steadyQuadraticDragCoriolisWind(
                        eastforce[x][y], northforce[x][y], latitudeforrow(y, height), drag,
                        config.surfaceBoundaryLayerDepthMetres,
                        config.rotationRatePerSecond, config.rotationDirection);
                    const int cell = y * (width + 1) + x;
                    deltau[cell] = corrected.eastMetresPerSecond - windu[x][y];
                    deltav[cell] = corrected.southMetresPerSecond - windv[x][y];
                    if (!carrierEast.empty())
                    {
                        const auto carrier = climateatmosphere::steadyQuadraticDragCoriolisWind(
                            carrierEastForce[x][y], carrierNorthForce[x][y], latitudeforrow(y, height), drag,
                            config.surfaceBoundaryLayerDepthMetres, config.rotationRatePerSecond, config.rotationDirection);
                        if (atmosphere::katabaticShallowTenMetreDiagnostic)
                        {
                            corrected = climateatmosphere::detail::coldJetTenMetreWind(carrier,
                                {coldEastForce[x][y], coldNorthForce[x][y]}, latitudeforrow(y, height),
                                drag * atmosphere::katabaticTenMetreStressFraction, atmosphere::katabaticShallowDepthMetres,
                                config.rotationRatePerSecond, config.rotationDirection);
                            deltau[cell] = corrected.eastMetresPerSecond - windu[x][y];
                            deltav[cell] = corrected.southMetresPerSecond - windv[x][y];
                        }
                        // Correct convergence of the layer mean, while retaining
                        // the full local velocity increment in the surface map.
                        deltau[cell] -= (1.0f - drainageMean) * ((corrected.eastMetresPerSecond -
                            carrier.eastMetresPerSecond) - (windu[x][y] - carrieru[x][y]));
                        deltav[cell] -= (1.0f - drainageMean) * ((corrected.southMetresPerSecond -
                            carrier.southMetresPerSecond) - (windv[x][y] - carrierv[x][y]));
                        carrieru[x][y] = carrier.eastMetresPerSecond;
                        carrierv[x][y] = carrier.southMetresPerSecond;
                    }
                    windu[x][y] = corrected.eastMetresPerSecond;
                    windv[x][y] = corrected.southMetresPerSecond;
                }
            // Retain the legacy upper contribution and update lower convergence
            // with the same finite-volume mass-flux diagnosis as the atmosphere.
            const auto ascentchange = climateatmosphere::detail::surfaceAscentCorrectionHpaPerDay(
                outputgrid, deltau, deltav, 0.5f * config.surfaceEquivalentPressureDepthHpa);
            for (int y = 0; y <= height; ++y)
                for (int x = 0; x <= width; ++x)
                    vertical[x][y] += ascentchange[y * (width + 1) + x];
        }
        // No raster-scale pressure relaxation may erase the accepted modes.
        if (!carrierEast.empty())
        {
            auto bulku = windu, bulkv = windv;
            auto du = windu, dv = windv;
            for (int y = 0; y <= height; ++y)
                for (int x = 0; x <= width; ++x)
                {
                    du[x][y] -= carrieru[x][y]; dv[x][y] -= carrierv[x][y];
                    bulku[x][y] -= (1.0f - drainageMean) * du[x][y];
                    bulkv[x][y] -= (1.0f - drainageMean) * dv[x][y];
                }
            climate.winds.drainagewindu[season] = std::move(du);
            climate.winds.drainagewindv[season] = std::move(dv);
            climate.winds.hasDrainageProfile[season] = true;
            storeterrainverticalmotion(world, season, orographicterrain, bulku, bulkv);
        }
        else storeterrainverticalmotion(world, season, orographicterrain, windu, windv);
        climate.winds.surfacewindu[season] = windu;
        climate.winds.surfacewindv[season] = windv;
        climate.winds.upperwindu[season] = upperu;
        climate.winds.upperwindv[season] = upperv;
        climate.winds.lowlevelwindu[season] = output(wind850East);
        climate.winds.lowlevelwindv[season] = output(wind850South);
        climate.winds.lowlevelavailable[season] = output(available850);
        climate.winds.hasLowLevelProfile[season] = true;
        climate.winds.populated[season] = true;
        climateatmosphere::CirculationPrecisionDiagnostics precision;
        precision.base = measurecirculationprecision(world, pressure, upperheight, landoutput, terrainoutput);
        precision.final = precision.base;
        climateatmosphere::setLastCirculationPrecisionDiagnostics(season, precision);
        for (int y = 0; y <= height; ++y)
            for (int x = 0; x <= width; ++x)
            {
                const auto windstorage = [](float value) { return std::clamp(value,
                    -atmosphere::maxVectorWind, atmosphere::maxVectorWind); };
                world.setseasonalpressurefloat(season, x, y, pressure[x][y]);
                world.setseasonalupperheightfloat(season, x, y, upperheight[x][y]);
                world.setseasonaluwindfloat(season, x, y, windstorage(windu[x][y]));
                world.setseasonalvwindfloat(season, x, y, windstorage(windv[x][y]));
                world.setseasonalupperuwindfloat(season, x, y, windstorage(upperu[x][y]));
                world.setseasonaluppervwindfloat(season, x, y, windstorage(upperv[x][y]));
                world.setseasonalverticalvelocityfloat(season, x, y, std::clamp(vertical[x][y],
                    -circulation::maximumVerticalVelocity, circulation::maximumVerticalVelocity));
            }
        std::cout << "Mode circulation season " << season << ": grid=" << columns << 'x' << rows
            << " surface_residual=" << result.surfaceStationarySolver.relativeResidual
            << " upper_residual=" << result.upperStationarySolver.relativeResidual
            << " surface_converged=" << result.surfaceStationarySolver.converged
            << " upper_converged=" << result.upperStationarySolver.converged
            << " stagnated=" << (result.surfaceStationarySolver.stagnated || result.upperStationarySolver.stagnated)
            << " N=" << stability << " depth_m=" << surfaceparameters.equivalentDepthMetres
            << " hadley_deg=" << hadleywidth << " thermal_equator_deg=" << confinedtrough
            << " control_thermal_equator_deg=" << thermalequator
            << " adjusted_zonal_pressure=" << config.adjustZonalSurfacePressure
            << " adjusted_zonal_upper=" << config.adjustZonalUpperHeight
            << " bandwidth=" << config.maximumZonalWavenumber << '/' << config.maximumMeridionalWavenumber
            << " kinetic_jm2=" << result.areaWeightedKineticEnergyJm2
            << " drag_wm2=" << result.areaWeightedDragDissipationWm2
            << " mass_kgm2=" << result.areaWeightedMassAnomalyKgM2
            << " momentum_exchange_residual=" << result.maximumMomentumExchangeResidual
            << " row_transfer_hpa=" << result.maximumStationaryRowMeanHpa << '\n';
        for (int mode = 0; mode < 2; ++mode)
        {
            const auto& solver = mode ? result.upperStationarySolver : result.surfaceStationarySolver;
            std::cout << "Atmosphere solver season=" << season << " mode=" << (mode ? "upper" : "surface")
                << " columns=" << columns << " rows=" << rows
                << " method=" << static_cast<int>(config.stationarySolver)
                << " spatial=" << solver.usedSpatialPreconditioner << " iterations=" << solver.iterations
                << " iteration_limit=" << config.maximumIterations << " tolerance=" << config.relativeTolerance
                << " residual=" << solver.relativeResidual << " returned_residual=" << solver.returnedRelativeResidual
                << " operator_seconds=" << solver.operatorSeconds << " preconditioner_seconds=" << solver.preconditionerSeconds
                << " orthogonalization_seconds=" << solver.orthogonalizationSeconds << '\n';
        }
    }


    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int y = startrow; y <= endrow; y++)
        {
            for (int x = 0; x <= width; x++)
            {
                float averageu = 0.0f;
                float averageupperu = 0.0f;

                for (int season = 0; season < CLIMATESEASONCOUNT; season++)
                {
                    averageu += world.seasonaluwindfloat(season, x, y);
                    averageupperu += world.seasonalupperuwindfloat(season, x, y);
                }

                averageu = averageu / static_cast<float>(CLIMATESEASONCOUNT);
                averageupperu = averageupperu / static_cast<float>(CLIMATESEASONCOUNT);

                int scalarwind = static_cast<int>(std::round(averageu / tuning::climate::atmosphere::scalarWindDivisor));
                scalarwind = std::clamp(scalarwind, -10, 10);

                if (std::fabs(averageu) < tuning::climate::atmosphere::minimumScalarZonalWind || scalarwind == 0)
                {
                    const float transportu = averageu + (averageupperu - averageu) *
                        tuning::climate::moistureadvection::upperWindTransportFraction;

                    if (transportu > 0.25f)
                        world.setwind(x, y, 101);
                    else if (transportu < -0.25f)
                        world.setwind(x, y, 99);
                    else
                        world.setwind(x, y, 0);
                }
                else
                {
                    if (scalarwind > 0)
                        scalarwind = std::max(1, scalarwind);
                    else
                        scalarwind = std::min(-1, scalarwind);

                    world.setwind(x, y, scalarwind);
                }
            }
        }
    });

    for (int y = 0; y <= height; y++)
        world.setwind(0, y, world.wind(width, y));
}
