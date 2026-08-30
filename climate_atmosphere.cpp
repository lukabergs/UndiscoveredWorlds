#include "climate_atmosphere.hpp"

#include <algorithm>
#include <cmath>

namespace climateatmosphere
{
namespace
{
constexpr float pi = 3.14159265358979323846f;
}

float coriolisParameterPerSecond(
    float latitudeDegrees,
    float rotationRatePerSecond,
    float rotationDirection)
{
    const float latitudeRadians = latitudeDegrees * pi / 180.0f;
    return 2.0f * rotationRatePerSecond * rotationDirection * std::sin(latitudeRadians);
}

float hypsometricHeightResponseMetresPerKelvin(
    float lowerPressurePa,
    float upperPressurePa,
    float dryAirGasConstant,
    float gravityMetresPerSecondSquared)
{
    if (lowerPressurePa <= upperPressurePa || upperPressurePa <= 0.0f ||
        dryAirGasConstant <= 0.0f || gravityMetresPerSecondSquared <= 0.0f)
    {
        return 0.0f;
    }

    return dryAirGasConstant / gravityMetresPerSecondSquared *
        std::log(lowerPressurePa / upperPressurePa);
}

CellSpacing cellSpacingMetres(
    float latitudeDegrees,
    int longitudeCells,
    int latitudeCells,
    float planetRadiusMetres)
{
    if (longitudeCells <= 0 || latitudeCells <= 1 || planetRadiusMetres <= 0.0f)
        return {};

    const float latitudeRadians = latitudeDegrees * pi / 180.0f;
    const float polarSafeCosine = std::max(0.02f, std::fabs(std::cos(latitudeRadians)));
    return {
        2.0f * pi * planetRadiusMetres * polarSafeCosine / static_cast<float>(longitudeCells),
        pi * planetRadiusMetres / static_cast<float>(latitudeCells - 1)
    };
}

HorizontalWind steadyRayleighCoriolisWind(
    float forceEastMetresPerSecondSquared,
    float forceNorthMetresPerSecondSquared,
    float latitudeDegrees,
    float dragTimeSeconds,
    float rotationRatePerSecond,
    float rotationDirection)
{
    if (dragTimeSeconds <= 0.0f)
        return {};

    const float dragRate = 1.0f / dragTimeSeconds;
    const float coriolis = coriolisParameterPerSecond(
        latitudeDegrees,
        rotationRatePerSecond,
        rotationDirection);
    const float denominator = dragRate * dragRate + coriolis * coriolis;

    if (denominator <= 0.0f)
        return {};

    const float east =
        (dragRate * forceEastMetresPerSecondSquared + coriolis * forceNorthMetresPerSecondSquared) /
        denominator;
    const float north =
        (dragRate * forceNorthMetresPerSecondSquared - coriolis * forceEastMetresPerSecondSquared) /
        denominator;
    return { east, -north };
}
}
