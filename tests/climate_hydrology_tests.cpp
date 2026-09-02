#include "climate_hydrology.hpp"
#include "climate_grid.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace
{
int failures = 0;

void expect(bool condition, const char* message)
{
    if (condition)
        return;

    std::cerr << "FAIL: " << message << '\n';
    failures++;
}

bool near(float first, float second, float tolerance = 1.0e-5f)
{
    return std::abs(first - second) <= tolerance;
}

double rotationError(int columns, int correctivePasses)
{
    constexpr double pi = 3.14159265358979323846;
    const int rows = columns / 2;
    const auto grid = climategrid::makeSphericalGrid(columns, rows, 1000.0);
    std::vector<float> tracer(columns * rows), u(tracer.size()), v(tracer.size(), 0.0f);
    for (int y = 0; y < rows; y++)
    {
        for (int x = 0; x < columns; x++)
        {
            tracer[grid.index(x, y)] = static_cast<float>(1.0 + 0.5 *
                std::cos(2.0 * pi * (x + 0.5) / columns) *
                std::cos(grid.latitudeCentresRadians[y]));
            u[grid.index(x, y)] = static_cast<float>(columns *
                grid.cellAreasSquareMetres[y] / grid.zonalFaceLengthsMetres[y]);
        }
    }
    std::vector<float> result;
    climatehydrology::MpdataOptions options;
    options.correctivePasses = correctivePasses;
    const auto diagnostics = climatehydrology::advectSphericalTracerMpdata(
        columns, rows, tracer, u, v, 1.0f, 1000.0f, options, result);
    double squaredError = 0.0, area = 0.0;
    for (std::size_t cell = 0; cell < tracer.size(); cell++)
    {
        const double weight = grid.cellAreasSquareMetres[cell / columns];
        const double error = result[cell] - tracer[cell];
        squaredError += weight * error * error;
        area += weight;
    }
    expect(diagnostics.minimumMixingRatio >= 0.0f &&
            std::abs(diagnostics.finalAreaWeightedMass / diagnostics.initialAreaWeightedMass - 1.0) < 2.0e-6,
        "full spherical rotation must preserve positivity and mass through the seam");
    return std::sqrt(squaredError / area);
}

void testVariableFlow()
{
    constexpr int columns = 32, rows = 16;
    constexpr double pi = 3.14159265358979323846;
    const auto grid = climategrid::makeSphericalGrid(columns, rows, 1000.0);
    std::vector<float> initial(columns * rows), u(initial.size()), v(initial.size());
    for (int y = 0; y < rows; y++)
        for (int x = 0; x < columns; x++)
            initial[grid.index(x, y)] = static_cast<float>(0.1 + std::exp(
                -std::pow((x - 15.0) / 5.0, 2) - std::pow((y - 7.5) / 3.0, 2)));
    const double initialMass = climategrid::areaWeightedIntegral(
        columns, rows, climategrid::LatitudeLayout::cellCentred, initial);
    auto tracer = initial;
    climatehydrology::MpdataOptions options;
    options.correctivePasses = 2;
    // Start at rest, accelerate through more than one cell per interval. The
    // integrated midpoint velocity must match a resolved sequence of steps.
    std::vector<float> accelerating(initial.size(), 0.0f), still(initial.size(), 0.0f), accelerated;
    for (int y = 0; y < rows; ++y)
        for (int x = 0; x < columns; ++x)
            accelerating[grid.index(x, y)] = static_cast<float>(4.0 * grid.cellAreasSquareMetres[y] / grid.zonalFaceLengthsMetres[y]);
    options.endZonalWindMps = accelerating;
    const auto acceleratingDiagnostics = climatehydrology::advectSphericalTracerMpdata(columns, rows,
        initial, still, still, 1.0f, 1000.0f, options, accelerated);
    expect(acceleratingDiagnostics.substeps >= 5 && accelerated != initial &&
        acceleratingDiagnostics.minimumMixingRatio >= 0.0f &&
        std::abs(acceleratingDiagnostics.finalAreaWeightedMass / acceleratingDiagnostics.initialAreaWeightedMass - 1.0) < 1.0e-6,
        "MPDATA must time-centre changing face winds with a complete-interval CFL bound");
    options.endZonalWindMps.clear();
    // Reversible, time-dependent sheared solid-body flow. Midpoint winds are
    // evaluated anew on every physical timestep; the second half retraces it.
    for (int step = 0; step < 80; step++)
    {
        for (int y = 0; y < rows; y++)
        {
            const double latitude = grid.latitudeCentresRadians[y];
            for (int x = 0; x < columns; x++)
            {
                u[grid.index(x, y)] = static_cast<float>(100.0 * std::cos(latitude) *
                    std::sin(2.0 * latitude) * std::cos(pi * (step + 0.5) / 80.0));
                v[grid.index(x, y)] = 0.0f;
            }
        }
        std::vector<float> next;
        climatehydrology::advectSphericalTracerMpdata(
            columns, rows, tracer, u, v, 0.5f, 1000.0f, options, next);
        tracer.swap(next);
    }
    double squaredError = 0.0;
    for (std::size_t cell = 0; cell < tracer.size(); cell++)
        squaredError += std::pow(tracer[cell] - initial[cell], 2);
    expect(std::sqrt(squaredError / tracer.size()) < 0.04,
        "reversible shear must recover tracer shape without substantial numerical blur");

    for (int sign : {-1, 1})
    {
        for (int y = 0; y < rows; y++)
        {
            for (int x = 0; x < columns; x++)
            {
                const double longitude = 2.0 * pi * (x + 0.5) / columns;
                u[grid.index(x, y)] = static_cast<float>(sign * 120.0 * std::sin(longitude));
                v[grid.index(x, y)] = static_cast<float>(sign * 60.0 *
                    std::sin(2.0 * grid.latitudeCentresRadians[y]));
            }
        }
        std::vector<float> next, repeat;
        climatehydrology::advectSphericalTracerMpdata(
            columns, rows, initial, u, v, 4.0f, 1000.0f, options, next);
        climatehydrology::advectSphericalTracerMpdata(
            columns, rows, initial, u, v, 4.0f, 1000.0f, options, repeat);
        const double finalMass = climategrid::areaWeightedIntegral(
            columns, rows, climategrid::LatitudeLayout::cellCentred, next);
        expect(next == repeat && std::all_of(next.begin(), next.end(), [](float q)
                { return std::isfinite(q) && q >= 0.0f; }) &&
                std::abs(finalMass / initialMass - 1.0) < 2.0e-6,
            "convergent and divergent face winds must remain deterministic, positive, and conservative");
    }
}
}

int main()
{
    const double coarseRotationError = rotationError(24, 1);
    const double fineRotationError = rotationError(48, 1);
    expect(fineRotationError < coarseRotationError * 0.75 &&
            fineRotationError < rotationError(48, 0),
        "MPDATA must converge under grid refinement and sharpen full rotations versus donor cell");
    testVariableFlow();
    const auto january = climatehydrology::calendarMonth(0);
    const auto march = climatehydrology::calendarMonth(2);
    const auto december = climatehydrology::calendarMonth(11);
    expect(january.firstSeason == 0 && january.secondSeason == 1,
        "January must begin at the January snapshot");
    expect(january.interpolation == 0.0f && january.days == 31,
        "January calendar metadata is wrong");
    expect(march.firstSeason == 0 && march.secondSeason == 1,
        "March must interpolate from January to April");
    expect(near(march.interpolation, 2.0f / 3.0f),
        "March interpolation phase is wrong");
    expect(december.firstSeason == 3 && december.secondSeason == 0,
        "December must interpolate periodically from October to January");
    expect(december.days == 31, "December day count is wrong");

    int annualDays = 0;
    for (int month = 0; month < climatehydrology::monthCount; month++)
        annualDays += climatehydrology::calendarMonth(month).days;
    expect(annualDays == 365, "hydrology calendar must span 365 days");
    expect(near(climatehydrology::interpolateSeasonal(10.0f, 16.0f, 1.0f / 3.0f), 12.0f),
        "seasonal interpolation is wrong");
    const auto climategrid = climatehydrology::climateGridDimensions(512, 257, 128);
    const auto oddRaster = climatehydrology::climateGridDimensions(63, 24, 128);
    expect(oddRaster.columns == 62 && oddRaster.rows == 31,
        "odd arbitrary output widths must still use an exactly 2:1 internal grid");
    expect(climategrid.columns == 128 && climategrid.rows == 64,
        "the reduced climate grid must use the W by W/2 internal contract");
    const auto nonstandardOutputGrid =
        climatehydrology::climateGridDimensions(600, 400, 64);
    expect(nonstandardOutputGrid.columns == 64 && nonstandardOutputGrid.rows == 32,
        "requested raster aspect ratios must not leak into the internal spherical grid");
    expect(climatehydrology::climateCellLatitudeDegrees(0, 64) < 90.0f &&
            climatehydrology::climateCellLatitudeDegrees(63, 64) > -90.0f,
        "the reduced climate grid must use finite-area cap cells instead of point poles");
    expect(near(
            climatehydrology::climateCellLatitudeDegrees(31, 64),
            -climatehydrology::climateCellLatitudeDegrees(32, 64)) &&
            climatehydrology::climateCellAreaWeight(0, 64) > 0.0f,
        "the even-row climate grid must remain equator-symmetric with positive cap area");
    expect(near(climatehydrology::polarTaperFactor(50.0f, 72.0f, 88.0f), 1.0f) &&
            near(climatehydrology::polarTaperFactor(90.0f, 72.0f, 88.0f), 0.0f) &&
            near(
                climatehydrology::polarTaperFactor(-80.0f, 72.0f, 88.0f),
                climatehydrology::polarTaperFactor(80.0f, 72.0f, 88.0f)),
        "polar transport tapering must be hemispherically symmetric and vanish at the poles");
    expect(climatehydrology::adjacentMeridionalTransportTargetRow(7, -8.0f, 64) == 6 &&
            climatehydrology::adjacentMeridionalTransportTargetRow(57, 8.0f, 64) == 58,
        "a meridional flux step must not leap across intervening latitude rows");
    expect(climatehydrology::adjacentMeridionalTransportTargetRow(0, -2.0f, 64) == 0 &&
            climatehydrology::adjacentMeridionalTransportTargetRow(64, 2.0f, 64) == 64,
        "meridional flux must remain inside the polar boundaries");

    constexpr int tracerColumns = 16;
    constexpr int tracerRows = tracerColumns / 2;
    constexpr float tracerRadius = 1000.0f;
    constexpr float tracerTimeStep = 1.0f;
    constexpr float pi = 3.14159265358979323846f;
    const size_t tracerCellCount = static_cast<size_t>(tracerColumns * tracerRows);
    const auto tracerIndex = [=](int x, int y)
    {
        return static_cast<size_t>(y * tracerColumns + x);
    };
    const auto tracerMass = [&](const std::vector<float>& field)
    {
        double mass = 0.0;
        for (int y = 0; y < tracerRows; y++)
        {
            const double area = climatehydrology::climateCellAreaWeight(y, tracerRows);
            for (int x = 0; x < tracerColumns; x++)
                mass += area * field[tracerIndex(x, y)];
        }
        return mass;
    };

    std::vector<float> tracer(tracerCellCount, 1.0f);
    std::vector<float> zonalWind(tracerCellCount, 8.0f);
    std::vector<float> meridionalWind(tracerCellCount, 0.0f);
    std::vector<float> transportedTracer;
    auto transportDiagnostics = climatehydrology::advectSphericalTracer(
        tracerColumns,
        tracerRows,
        tracer,
        zonalWind,
        meridionalWind,
        tracerTimeStep,
        tracerRadius,
        0.5f,
        48.0f,
        transportedTracer);
    bool uniformZonalTracer = transportedTracer.size() == tracerCellCount;
    for (float value : transportedTracer)
        uniformZonalTracer = uniformZonalTracer && near(value, 1.0f, 2.0e-5f);
    expect(uniformZonalTracer,
        "periodic zonal transport must preserve a uniform tracer");
    expect(std::abs(
            transportDiagnostics.finalAreaWeightedMass -
                transportDiagnostics.initialAreaWeightedMass) < 1.0e-4,
        "zonal tracer transport must conserve area-weighted mass");

    tracer.assign(tracerCellCount, 0.0f);
    zonalWind.assign(tracerCellCount, 0.0f);
    const int zonalSourceX = 4;
    const int equatorRow = tracerRows / 2;
    tracer[tracerIndex(zonalSourceX, equatorRow)] = 1.0f;
    const float sourceLatitude = climatehydrology::climateCellLatitudeDegrees(
        equatorRow, tracerRows) * pi / 180.0f;
    const float equatorialCellWidth = 2.0f * pi * tracerRadius *
        std::cos(sourceLatitude) / static_cast<float>(tracerColumns);
    for (int x = 0; x < tracerColumns; x++)
        zonalWind[tracerIndex(x, equatorRow)] = 2.25f * equatorialCellWidth;
    climatehydrology::MpdataOptions donorOptions;
    donorOptions.maximumCourantPerSubstep = 0.5f;
    donorOptions.correctivePasses = 0;
    std::vector<float> donorTracer;
    climatehydrology::advectSphericalTracerMpdata(
        tracerColumns,
        tracerRows,
        tracer,
        zonalWind,
        meridionalWind,
        tracerTimeStep,
        tracerRadius,
        donorOptions,
        donorTracer);
    climatehydrology::MpdataOptions mpdataOptions = donorOptions;
    mpdataOptions.correctivePasses = 1;
    transportDiagnostics = climatehydrology::advectSphericalTracerMpdata(
        tracerColumns,
        tracerRows,
        tracer,
        zonalWind,
        meridionalWind,
        tracerTimeStep,
        tracerRadius,
        mpdataOptions,
        transportedTracer);
    const auto zonalMoments = [&](const std::vector<float>& field)
    {
        double total = 0.0;
        double centroid = 0.0;
        for (int x = 0; x < tracerColumns; x++)
        {
            int offset = x - zonalSourceX;
            if (offset > tracerColumns / 2)
                offset -= tracerColumns;
            if (offset < -tracerColumns / 2)
                offset += tracerColumns;
            total += field[tracerIndex(x, equatorRow)];
            centroid += static_cast<double>(offset) * field[tracerIndex(x, equatorRow)];
        }
        centroid /= total;
        double variance = 0.0;
        for (int x = 0; x < tracerColumns; x++)
        {
            int offset = x - zonalSourceX;
            if (offset > tracerColumns / 2)
                offset -= tracerColumns;
            if (offset < -tracerColumns / 2)
                offset += tracerColumns;
            const double delta = static_cast<double>(offset) - centroid;
            variance += delta * delta * field[tracerIndex(x, equatorRow)] / total;
        }
        return std::pair<double, double>{ centroid, variance };
    };
    const auto donorMoments = zonalMoments(donorTracer);
    const auto mpdataMoments = zonalMoments(transportedTracer);
    expect(std::abs(mpdataMoments.first - 2.25) < 0.05,
        "super-CFL face transport must preserve the analytical zonal centroid");
    expect(mpdataMoments.second < donorMoments.second,
        "the MPDATA corrective pass must diffuse less than donor-cell transport");
    expect(std::abs(
            transportDiagnostics.finalAreaWeightedMass /
                transportDiagnostics.initialAreaWeightedMass - 1.0) < 1.0e-6,
        "multi-cell zonal flux transport must conserve tracer mass");

    tracer.assign(tracerCellCount, 0.0f);
    zonalWind.assign(tracerCellCount, 0.0f);
    meridionalWind.assign(tracerCellCount, 0.0f);
    for (int y = 1; y < tracerRows - 1; y++)
        meridionalWind[tracerIndex(5, y)] = 2.2f * pi * tracerRadius /
            static_cast<float>(tracerRows);
    const int meridionalSourceRow = 2;
    tracer[tracerIndex(5, meridionalSourceRow)] = 1.0f /
        climatehydrology::climateCellAreaWeight(meridionalSourceRow, tracerRows);
    transportDiagnostics = climatehydrology::advectSphericalTracer(
        tracerColumns,
        tracerRows,
        tracer,
        zonalWind,
        meridionalWind,
        tracerTimeStep,
        tracerRadius,
        0.5f,
        48.0f,
        transportedTracer);
    double meridionalCentroid = 0.0;
    double meridionalVariance = 0.0;
    const double transportedMass = tracerMass(transportedTracer);
    for (int y = 0; y < tracerRows; y++)
    {
        const double cellMass = climatehydrology::climateCellAreaWeight(y, tracerRows) *
            transportedTracer[tracerIndex(5, y)];
        meridionalCentroid += static_cast<double>(y) * cellMass;
    }
    meridionalCentroid /= transportedMass;
    for (int y = 0; y < tracerRows; y++)
    {
        const double cellMass = climatehydrology::climateCellAreaWeight(y, tracerRows) *
            transportedTracer[tracerIndex(5, y)];
        const double offset = static_cast<double>(y) - meridionalCentroid;
        meridionalVariance += offset * offset * cellMass / transportedMass;
    }
    expect(transportDiagnostics.substeps > 1 &&
            transportDiagnostics.maximumMultidimensionalCourant > 1.0f,
        "the spherical multidimensional CFL condition must select stable substeps");
    expect(meridionalCentroid > static_cast<double>(meridionalSourceRow) + 0.5 &&
            meridionalVariance > 0.0,
        "super-CFL meridional fluxes must traverse intervening rows without teleportation");
    expect(std::abs(
            transportDiagnostics.finalAreaWeightedMass /
                transportDiagnostics.initialAreaWeightedMass - 1.0) < 1.0e-6,
        "meridional transport must conserve area-weighted tracer mass");
    bool positiveTracer = true;
    for (float value : transportedTracer)
        positiveTracer = positiveTracer && value >= 0.0f;
    expect(positiveTracer, "tracer transport must remain positive");

    tracer.assign(tracerCellCount, 0.0f);
    meridionalWind.assign(tracerCellCount, -4.0f * pi * tracerRadius /
        static_cast<float>(tracerRows));
    tracer[tracerIndex(3, 0)] = 1.0f;
    climatehydrology::advectSphericalTracer(
        tracerColumns,
        tracerRows,
        tracer,
        zonalWind,
        meridionalWind,
        tracerTimeStep,
        tracerRadius,
        0.5f,
        48.0f,
        transportedTracer);
    expect(near(transportedTracer[tracerIndex(3, 0)], 1.0f) &&
            std::abs(tracerMass(transportedTracer) - tracerMass(tracer)) < 1.0e-6,
        "outward polar transport must neither leak nor jump across the cap");

    const auto daytime = climatehydrology::deterministicWeatherPhase(
        0, 3, 7.0f, 2.0f, -2.0f, 0.35f, -0.35f);
    const auto neutral = climatehydrology::deterministicWeatherPhase(
        1, 3, 7.0f, 2.0f, -2.0f, 0.35f, -0.35f);
    const auto nighttime = climatehydrology::deterministicWeatherPhase(
        2, 3, 7.0f, 2.0f, -2.0f, 0.35f, -0.35f);
    expect(daytime.coastalDirection > 0.0f && nighttime.coastalDirection < 0.0f &&
            near(neutral.coastalDirection, 0.0f),
        "weather phases must include deterministic onshore, neutral, and offshore flow");
    expect(near(
            daytime.landTemperatureAnomalyC + neutral.landTemperatureAnomalyC +
                nighttime.landTemperatureAnomalyC,
            0.0f),
        "weather-phase land temperature anomalies must not bias the monthly mean");

    expect(near(climatehydrology::kuoPrecipitationEfficiency(0.0f, 0.0f, 3.0f), 0.0f),
        "a dry column must have zero Kuo precipitation efficiency");
    expect(climatehydrology::kuoPrecipitationEfficiency(0.8f, 0.0f, 3.0f) >
            climatehydrology::kuoPrecipitationEfficiency(0.4f, 0.0f, 3.0f),
        "Kuo precipitation efficiency must increase continuously with humidity");
    expect(near(climatehydrology::convectiveBuoyancyEfficiency(0.0f, 0.0f, 6.0f), 0.0f) &&
            near(climatehydrology::convectiveBuoyancyEfficiency(6.0f, 0.0f, 6.0f), 1.0f),
        "convective buoyancy closure must span its configured activation range");
    const float shallowexchange = climatehydrology::shallowConvectionExchangeFraction(
        0.8f, 0.6f, 4.0f, 10.0f, 86400.0f, 2.0f, 0.55f, 0.85f, 20.0f, 0.30f);
    expect(shallowexchange > 0.0f && shallowexchange <= 0.30f,
        "humid unstable columns must undergo bounded shallow-convective mixing");
    expect(near(climatehydrology::dryConvectionExchangeFraction(
            8.0f, 86400.0f, 1.0f, 10.0f, 18.0f, 0.25f), 0.0f),
        "stable or weakly unstable columns must not undergo dry-convective mixing");
    expect(climatehydrology::dryConvectionExchangeFraction(
            18.0f, 86400.0f, 1.0f, 10.0f, 18.0f, 0.25f) > 0.0f,
        "strongly unstable columns must undergo dry-convective mixing");
    expect(near(climatehydrology::soilMoistureStress(0.0f, 80.0f, 0.5f, 0.5f), 0.0f),
        "empty soil must suppress evapotranspiration");
    expect(near(climatehydrology::soilMoistureStress(40.0f, 80.0f, 0.5f, 0.5f), 1.0f),
        "field-capacity soil must permit potential evapotranspiration");
    expect(
        climatehydrology::soilMoistureStress(20.0f, 80.0f, 0.5f, 0.5f) >
            climatehydrology::soilMoistureStress(10.0f, 80.0f, 0.5f, 0.5f),
        "evapotranspiration stress must increase monotonically with soil water");
    expect(near(climatehydrology::diagnosticCloudFraction(0.70f, 0.75f), 0.0f),
        "sub-cloud humidity must not create cloud cover");
    expect(climatehydrology::diagnosticCloudFraction(0.90f, 0.75f) > 0.0f,
        "humid sub-saturated air must be able to form cloud without raining");

    const auto adjustment = climatehydrology::moistSaturationAdjustment(
        30.0f, 20.0f, 10.0f, 86400.0f, 172800.0f, 2, 0.35f, 0.065f);
    expect(adjustment.condensedMm > 0.0f && adjustment.condensedMm < 10.0f,
        "moist adjustment must condense only part of a supersaturated column");
    expect(adjustment.adjustedTemperatureC > 10.0f,
        "condensation must warm the adjusted parcel");
    expect(near(adjustment.condensedMm + adjustment.remainingVapourMm, 30.0f),
        "moist adjustment must conserve column water");

    const auto exchange = climatehydrology::exchangeMoistureLayers(
        20.0f, 10.0f, 0.25f, 0.10f);
    expect(near(exchange.boundaryLayerMm, 16.0f) &&
            near(exchange.freeTroposphereMm, 14.0f),
        "vertical exchange must apply simultaneous up- and downward transfers");
    expect(near(exchange.boundaryLayerMm + exchange.freeTroposphereMm, 30.0f),
        "vertical exchange must conserve atmospheric water");

    const auto falling = climatehydrology::processFallingPrecipitation(
        10.0f, -5.0f, 0.25f, 0.4f, 8.0f, -1.0f, 2.0f);
    expect(falling.reevaporatedMm > 0.0f && near(falling.rainMm, 0.0f),
        "cold precipitation must sublimate and reach the surface as snow");
    expect(near(falling.surfaceTotalMm() + falling.reevaporatedMm, 10.0f),
        "falling precipitation processing must conserve water");
    expect(near(climatehydrology::snowMeltAmount(20.0f, 5.0f, 86400.0f, 3.0f), 15.0f),
        "degree-day snowmelt must respect temperature and duration");
    const auto snowaccumulation = climatehydrology::accumulateSnowfall(
        4995.0f, 12.0f, 5000.0f);
    expect(near(snowaccumulation.storageMm, 5000.0f) &&
            near(snowaccumulation.overflowMm, 7.0f),
        "snow above the finite glacier store must become runoff");
    expect(near(snowaccumulation.storageMm + snowaccumulation.overflowMm, 5007.0f),
        "snow accumulation and overflow must conserve water");

    const auto dry = climatehydrology::partitionPrecipitation(
        20.0f, 50.0f, 50.0f, -2.0f, 0.0f, 25.0f, 86400.0f,
        0.8f, 172800.0f, 0.6f, 0.75f, 8.0f, 24.0f);
    expect(near(dry.totalMm(), 0.0f),
        "a subcritical divergent column must not precipitate");

    const auto stratiform = climatehydrology::partitionPrecipitation(
        50.0f, 50.0f, 50.0f, 0.0f, 0.0f, 0.0f, 86400.0f,
        0.8f, 172800.0f, 0.6f, 0.75f, 8.0f, 24.0f);
    expect(stratiform.stratiformMm > 0.0f && near(stratiform.orographicMm, 0.0f),
        "uniform saturation excess must be stratiform");
    expect(near(stratiform.convectiveMm, 0.0f),
        "cold air must suppress the convective closure");

    const auto terrain = climatehydrology::partitionPrecipitation(
        38.0f, 50.0f, 40.0f, 0.0f, 0.0f, 20.0f, 86400.0f,
        0.8f, 172800.0f, 0.6f, 0.75f, 8.0f, 24.0f);
    expect(terrain.orographicMm > 0.0f,
        "terrain cooling must create a distinct orographic contribution");

    const auto convective = climatehydrology::partitionPrecipitation(
        38.0f, 50.0f, 50.0f, 4.0f, 0.0f, 30.0f, 86400.0f,
        0.8f, 172800.0f, 0.6f, 0.75f, 8.0f, 24.0f);
    expect(convective.convectiveMm > 0.0f,
        "warm convergent moisture supply must trigger convection");
    expect(convective.totalMm() <= 38.0f,
        "precipitation partition must conserve atmospheric water");

    const auto divergent = climatehydrology::partitionPrecipitation(
        38.0f, 50.0f, 50.0f, -4.0f, 0.0f, 30.0f, 86400.0f,
        0.8f, 172800.0f, 0.6f, 0.75f, 8.0f, 24.0f);
    expect(near(divergent.convectiveMm, 0.0f),
        "moisture divergence must not feed the convective closure");

    const auto evaporationFed = climatehydrology::partitionPrecipitation(
        38.0f, 50.0f, 50.0f, -1.0f, 3.0f, 30.0f, 86400.0f,
        0.8f, 172800.0f, 0.6f, 0.75f, 8.0f, 24.0f);
    expect(evaporationFed.convectiveMm > 0.0f,
        "surface evaporation plus net moisture supply must feed warm convection");

    const auto layered = climatehydrology::partitionTwoLayerPrecipitation(
        24.0f, 18.0f, 30.0f, 15.0f, 12.0f, 4.0f, 1.0f,
        28.0f, 12.0f, 86400.0f, 172800.0f, 0.6f, 0.75f,
        8.0f, 24.0f, 2, 0.35f, 0.065f);
    expect(layered.stratiformMm > 0.0f && layered.orographicMm > 0.0f,
        "a supersaturated lifted free troposphere must rain stratiformly and orographically");
    expect(layered.convectiveMm > 0.0f,
        "warm convergent boundary-layer moisture must rain convectively");
    expect(layered.stratiformMm + layered.orographicMm <= 18.0f + 1.0e-5f &&
            layered.convectiveMm <= 24.0f + 1.0e-5f,
        "two-layer precipitation must conserve each source reservoir");

    const auto elevatedFed = climatehydrology::partitionTwoLayerPrecipitation(
        24.0f, 10.0f, 30.0f, 15.0f, 15.0f, -1.0f, 0.0f,
        28.0f, 16.0f, 86400.0f, 172800.0f, 0.0f, 0.75f,
        0.0f, 6.0f, 2, 0.35f, 0.065f, 6.0f, 0.5f, 3.0f);
    expect(elevatedFed.convectiveMm > 0.0f,
        "elevated free-tropospheric moisture accession must be able to feed convection");

    for (int water = 0; water <= 100; water += 5)
    {
        const auto partition = climatehydrology::partitionPrecipitation(
            static_cast<float>(water), 50.0f, 42.0f, 6.0f, 2.0f, 28.0f,
            86400.0f, 0.8f, 172800.0f, 0.6f, 0.75f, 8.0f, 24.0f);
        expect(partition.stratiformMm >= 0.0f && partition.orographicMm >= 0.0f &&
                partition.convectiveMm >= 0.0f,
            "precipitation components must stay non-negative");
        expect(partition.totalMm() <= static_cast<float>(water) + 1.0e-5f,
            "precipitation must never remove more water than the column contains");
    }

    if (failures == 0)
        std::cout << "All climate hydrology tests passed\n";

    return failures == 0 ? 0 : 1;
}
