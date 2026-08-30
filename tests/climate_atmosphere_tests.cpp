#include "climate_atmosphere.hpp"

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
    constexpr float earthRotation = 7.2921159e-5f;
    expect(
        std::abs(climateatmosphere::coriolisParameterPerSecond(0.0f, earthRotation)) < 1.0e-9f,
        "Coriolis acceleration must vanish at the equator");
    expect(
        climateatmosphere::coriolisParameterPerSecond(45.0f, earthRotation) > 0.0f &&
            climateatmosphere::coriolisParameterPerSecond(-45.0f, earthRotation) < 0.0f,
        "Coriolis sign must reverse across the equator");

    const float heightResponse = climateatmosphere::hypsometricHeightResponseMetresPerKelvin(
        100000.0f, 50000.0f);
    expect(
        std::abs(heightResponse - 20.27f) < 0.05f,
        "1000-to-500 hPa thickness must change by about 20.27 metres per kelvin");

    const auto spacing = climateatmosphere::cellSpacingMetres(0.0f, 2048, 1025, 6371000.0f);
    expect(
        spacing.zonalMetres > 19000.0f && spacing.zonalMetres < 20000.0f &&
            spacing.meridionalMetres > 19000.0f && spacing.meridionalMetres < 20000.0f,
        "Earth benchmark cells must be about 19.5 kilometres at the equator");

    const auto unrotated = climateatmosphere::steadyRayleighCoriolisWind(
        0.001f, 0.002f, 0.0f, 1000.0f, earthRotation);
    expect(
        std::abs(unrotated.eastMetresPerSecond - 1.0f) < 0.001f &&
            std::abs(unrotated.southMetresPerSecond + 2.0f) < 0.001f,
        "Rayleigh flow at the equator must follow the pressure-gradient acceleration");

    const auto northern = climateatmosphere::steadyRayleighCoriolisWind(
        0.0f, 0.001f, 45.0f, 86400.0f, earthRotation);
    const auto southern = climateatmosphere::steadyRayleighCoriolisWind(
        0.0f, 0.001f, -45.0f, 86400.0f, earthRotation);
    expect(
        northern.eastMetresPerSecond > 0.0f && southern.eastMetresPerSecond < 0.0f,
        "the same meridional height force must produce opposite zonal flow across the equator");

    return failures == 0 ? 0 : 1;
}
