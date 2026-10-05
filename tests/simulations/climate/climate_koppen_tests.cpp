#include "simulations/climate/climate_koppen.hpp"

#include <cstdlib>
#include <iostream>
#include <gtest/gtest.h>

namespace
{
void expect(bool condition, const char* message)
{
    EXPECT_TRUE(condition) << message;
}
}

namespace
{

// A temperate northern-hemisphere year: monthly temperatures with 80 mm of rain each month.
class ClimateKoppen : public ::testing::Test
{
protected:
    std::array<float, 12> temperature = {5, 6, 9, 12, 16, 20, 24, 23, 19, 14, 9, 6};
    std::array<float, 12> rain;

    ClimateKoppen()
    {
        rain.fill(80.0f);
    }
};

}

TEST_F(ClimateKoppen, TemperateDrySummer)
{
    expect(climatekoppen::classifyMonthly(temperature, rain, true) == 15, "humid hot-summer climate must be Cfa");
    rain[6] = 5.0f;
    expect(climatekoppen::classifyMonthly(temperature, rain, true) == 9,
        "one dry summer month must produce Csa even when its quarterly mean is wet");
    std::array<float, 12> southernTemperature, southernRain;
    for (int month = 0; month < 12; ++month)
    {
        southernTemperature[month] = temperature[(month + 6) % 12];
        southernRain[month] = rain[(month + 6) % 12];
    }
    expect(climatekoppen::classifyMonthly(southernTemperature, southernRain, false) == 9,
        "hemisphere reversal must preserve the physical classification");
}

TEST_F(ClimateKoppen, TemperateDryWinter)
{
    rain[0] = 1.0f;
    expect(climatekoppen::classifyMonthly(temperature, rain, true) == 12, "dry winter must produce Cwa");
    temperature[0] = -1.0f;
    expect(climatekoppen::classifyMonthly(temperature, rain, true) == 22, "monthly C/D boundary must use 0 C");
}

TEST_F(ClimateKoppen, Tropical)
{
    temperature.fill(26.0f); rain.fill(200.0f);
    expect(climatekoppen::classifyMonthly(temperature, rain, true) == 1, "wet tropical climate must be Af");
    rain[0] = 50.0f;
    expect(climatekoppen::classifyMonthly(temperature, rain, true) == 2, "short tropical dry season must be Am");
}

TEST_F(ClimateKoppen, Arid)
{
    temperature.fill(26.0f);
    rain.fill(0.0f);
    expect(climatekoppen::classifyMonthly(temperature, rain, true) == 5, "hot arid climate must be BWh");
}

TEST_F(ClimateKoppen, Polar)
{
    temperature.fill(-5.0f); rain.fill(80.0f);
    expect(climatekoppen::classifyMonthly(temperature, rain, true) == 31, "year-round freezing climate must be EF");
    temperature[6] = 4.0f;
    expect(climatekoppen::classifyMonthly(temperature, rain, true) == 30, "positive polar summer must be ET");
}

TEST_F(ClimateKoppen, WinterDryRule)
{
    expect(
        climatekoppen::isWinterDry(9.9f, 100.0f),
        "winter precipitation below one tenth of the summer maximum must be dry-winter");
    expect(
        !climatekoppen::isWinterDry(10.0f, 100.0f),
        "the one-tenth boundary itself must not be dry-winter");
    expect(
        !climatekoppen::isWinterDry(25.0f, 100.0f),
        "the obsolete one-quarter temperate rule must not classify dry-winter");
}
