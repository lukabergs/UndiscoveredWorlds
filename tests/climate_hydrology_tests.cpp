#include "climate_hydrology.hpp"

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
}

int main()
{
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
    expect(climategrid.columns == 128 && climategrid.rows == 65,
        "the reduced climate grid must preserve the output aspect ratio and poles");
    expect(climatehydrology::climateCellLatitudeDegrees(0, 65) < 90.0f &&
            climatehydrology::climateCellLatitudeDegrees(64, 65) > -90.0f,
        "the reduced climate grid must use finite-area cap cells instead of point poles");
    expect(near(climatehydrology::climateCellLatitudeDegrees(32, 65), 0.0f) &&
            climatehydrology::climateCellAreaWeight(0, 65) > 0.0f,
        "the reduced climate grid must remain equator-symmetric with positive cap area");
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
