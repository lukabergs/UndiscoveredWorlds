#include "climate_energy.hpp"
#include "detail/land_ice.hpp"
#include "detail/tropical_cloud_albedo.hpp"

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
    for (double latitude : {0.0, 12.0, 30.0, 60.0, 90.0})
    {
        const auto cloud = climateenergy::detail::tropicalCloudAlbedo(latitude, 0.04, 12.0);
        expect(cloud >= 0.0 && cloud <= 0.04, "prescribed cloud albedo must remain bounded");
        expect(cloud == climateenergy::detail::tropicalCloudAlbedo(-latitude, 0.04, 12.0),
            "annual zonal cloud closure must preserve hemispheric symmetry");
        expect(climateenergy::detail::tropicalCloudAlbedo(latitude, 0.0, 12.0) == 0.0,
            "zero cloud amplitude must recover the previous shortwave forcing");
    }
    expect(climateenergy::detail::tropicalCloudAlbedo(30.0, 0.04, 12.0) < 0.002,
        "tropical cloud shading must become small in the subtropics");
    expect(climateenergy::detail::tropicalCloudAlbedo(60.0, 0.04, 12.0) < 1e-6,
        "tropical cloud shading must not directly cool polar regions");
    for (double dt : {3600.0, 28800.0, 86400.0, 365.0 * 86400.0})
    {
        const climateenergy::LandThermalState initial{25.0, -5.0};
        const auto next = climateenergy::stepLandThermalState(initial, 100.0, 7.5e6, 3.0e7, 1.5, dt);
        const double storage = 7.5e6 * (next.surfaceC - initial.surfaceC) + 3.0e7 * (next.deepC - initial.deepC);
        expect(std::abs(storage - 100.0 * dt) < 1.0e-6 * dt,
            "land surface and deep reservoirs must jointly conserve the supplied heat");
        const auto isolated = climateenergy::stepLandThermalState(initial, 0.0, 7.5e6, 3.0e7, 1.5, dt);
        expect(isolated.surfaceC < initial.surfaceC && isolated.deepC > initial.deepC && isolated.surfaceC >= isolated.deepC,
            "implicit soil exchange must relax a temperature contrast without overshooting");
    }
    for (int rows : {16, 64, 128})
    {
        std::vector<double> source(rows), response(rows, 0.5);
        for (int y = 0; y < rows; ++y) source[y] = std::cos((y + 0.5) * 3.141592653589793 / rows);
        const auto result = climateenergy::implicitMeridionalDiffusion(source, response, 0.6);
        double heat = 0.0, error = 0.0;
        for (int y = 0; y < rows; ++y)
        {
            const double area = std::cos(y * 3.141592653589793 / rows) - std::cos((y + 1) * 3.141592653589793 / rows);
            heat += area * result.heatConvergenceWm2[y];
            error = std::max(error, std::abs(result.temperatureC[y] - source[y] / 1.6));
            expect(std::abs(result.temperatureC[y] - source[y] - response[y] * result.heatConvergenceWm2[y]) < 1.0e-10,
                "implicit spherical transport must satisfy each local heat budget");
        }
        expect(std::abs(heat) < 1.0e-10, "polar-closed diffusion must exchange zero net global heat");
        expect(error < 1.0 / (rows * rows), "spherical dipole diffusion must approach its analytic decay rate");
    }
    for (float eccentricity : {0.0f, 0.3f, 0.8f, 0.95f})
    {
        double distanceMean = 0.0;
        constexpr int samples = 20000;
        for (int i = 0; i < samples; ++i)
            distanceMean += climateenergy::orbitalDistanceFactor(365.0f * (i + 0.5f) / samples, eccentricity, 0) / samples;
        expect(std::abs(distanceMean * std::sqrt(1.0 - eccentricity * eccentricity) - 1.0) < 2.0e-5,
            "time-averaged inverse-square orbital flux must satisfy Kepler's law");
        const auto equinox = climateenergy::solarForcing(80.0f, 23.44f, eccentricity, 0);
        expect(std::abs(equinox.declinationRadians) < 1.0e-6f, "eccentric orbits must retain the day-80 equinox");
        const auto first = climateenergy::solarForcing(123.0f, 23.44f, eccentricity, 1);
        const auto next = climateenergy::solarForcing(488.0f, 23.44f, eccentricity, 1);
        expect(std::abs(first.distanceFactor - next.distanceFactor) < 1.0e-6f &&
            std::abs(first.declinationRadians - next.declinationRadians) < 1.0e-6f,
            "orbital forcing must be periodic across years");
    }
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

    expect(
        climateenergy::permanentLandIceFraction(-14.0f) == 1.0f,
        "land whose warmest season reaches the cold transition must retain full ice cover");
    expect(
        std::abs(climateenergy::permanentLandIceFraction(-9.0f) - 0.5f) < 1.0e-6f,
        "the warm-season permanent-ice transition midpoint must have half coverage");
    expect(
        climateenergy::permanentLandIceFraction(-4.0f) == 0.0f,
        "land whose warmest season reaches the warm transition must lose permanent ice");

    // At 5 km this profile has a -17.25 C seasonal mean but a +5.5 C
    // summer. Cold annual/deep conditions must not bypass summer ice loss.
    constexpr double highlandOffset = -5.0 * 6.5;
    const std::array<double, CLIMATESEASONCOUNT> warmSummerHighland = {-6.0, 20.0, 38.0, 9.0};
    expect(climateenergy::detail::seasonalPermanentLandIceFraction(warmSummerHighland, highlandOffset) == 0.0f,
        "highlands with a thawing summer must lose permanent ice despite a cold annual temperature");
    expect(climateenergy::detail::seasonalPermanentLandIceFraction({38.0, 9.0, -6.0, 20.0}, highlandOffset) == 0.0f,
        "summer ice survival must follow the warmest season in either hemisphere");
    expect(climateenergy::detail::seasonalPermanentLandIceFraction({-12.0, 4.0, 18.5, 0.0}, highlandOffset) == 1.0f,
        "highlands whose summer stays below the cold transition must retain permanent ice");
    expect(std::abs(climateenergy::detail::seasonalPermanentLandIceFraction(
            {-12.0, 4.0, 23.5, 0.0}, highlandOffset) - 0.5f) < 1.0e-6f,
        "highland ice survival must retain the existing warm-season transition");
    expect(climateenergy::detail::seasonalPermanentLandIceFraction(
            {-44.5, -28.5, -9.0, -32.5}, 0.0) ==
        climateenergy::detail::seasonalPermanentLandIceFraction({-12.0, 4.0, 23.5, 0.0}, highlandOffset),
        "the local lapse correction must enter permanent ice diagnosis exactly once");

    return failures == 0 ? 0 : 1;
}
