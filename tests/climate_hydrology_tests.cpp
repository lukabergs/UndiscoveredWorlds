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
    expect(near(climatehydrology::soilMoistureStress(0.0f, 80.0f, 0.5f, 0.5f), 0.0f),
        "empty soil must suppress evapotranspiration");
    expect(near(climatehydrology::soilMoistureStress(40.0f, 80.0f, 0.5f, 0.5f), 1.0f),
        "field-capacity soil must permit potential evapotranspiration");
    expect(
        climatehydrology::soilMoistureStress(20.0f, 80.0f, 0.5f, 0.5f) >
            climatehydrology::soilMoistureStress(10.0f, 80.0f, 0.5f, 0.5f),
        "evapotranspiration stress must increase monotonically with soil water");

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
