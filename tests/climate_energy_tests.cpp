#include "climate_energy.hpp"

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
}

int main()
{
    const float equinoxDeclination = climateenergy::solarDeclinationRadians(80.0f, 23.44f);
    expect(std::abs(equinoxDeclination) < 1.0e-5f, "day 80 must be an equinox in the orbital convention");

    const float equatorialEquinox = climateenergy::dailyMeanInsolationWm2(0.0f, equinoxDeclination);
    expect(std::abs(equatorialEquinox - 433.22f) < 0.2f, "equatorial equinox insolation must equal S0/pi");

    const float northernSummerDeclination = climateenergy::solarDeclinationRadians(172.0f, 23.44f);
    const float northernWinterDeclination = climateenergy::solarDeclinationRadians(355.0f, 23.44f);
    expect(
        climateenergy::dailyMeanInsolationWm2(65.0f, northernSummerDeclination) >
            climateenergy::dailyMeanInsolationWm2(65.0f, northernWinterDeclination),
        "northern high latitudes must receive more summer than winter energy");
    expect(
        climateenergy::dailyMeanInsolationWm2(90.0f, northernWinterDeclination) == 0.0f,
        "the winter pole must enter polar night");

    expect(
        climateenergy::orbitalDistanceFactor(3.0f, 0.1f, 0) >
            climateenergy::orbitalDistanceFactor(186.0f, 0.1f, 0),
        "solar flux must be greater at perihelion than aphelion");

    constexpr float previous = 14.0f;
    constexpr float absorbed = 300.0f;
    constexpr float intercept = 210.0f;
    constexpr float slope = 2.0f;
    constexpr float transport = 3.0f;
    constexpr float seconds = 86400.0f;
    const float land = climateenergy::implicitSlabTemperatureStep(
        previous, absorbed, intercept, slope, transport, 14.0f, 2.0e7f, seconds);
    const float ocean = climateenergy::implicitSlabTemperatureStep(
        previous, absorbed, intercept, slope, transport, 14.0f, 1.25e8f, seconds);
    expect(land > ocean && ocean > previous, "lower land heat capacity must produce a faster warming response");

    const float storage = 2.0e7f * (land - previous) / seconds;
    const float netFlux = absorbed - intercept - slope * land + transport * (14.0f - land);
    expect(std::abs(storage - netFlux) < 0.001f, "implicit slab step must close its energy equation");

    return failures == 0 ? 0 : 1;
}
