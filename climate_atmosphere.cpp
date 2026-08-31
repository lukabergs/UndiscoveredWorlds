#include "climate_atmosphere.hpp"

#include <algorithm>
#include <cmath>

namespace climateatmosphere
{
namespace
{
constexpr float pi = 3.14159265358979323846f;
std::array<CirculationPrecisionDiagnostics, CLIMATESEASONCOUNT> circulationPrecisionDiagnostics{};

float smoothstep(float value)
{
    const float clamped = std::clamp(value, 0.0f, 1.0f);
    return clamped * clamped * (3.0f - 2.0f * clamped);
}

float smoothinterpolate(float start, float end, float fraction)
{
    return start + (end - start) * smoothstep(fraction);
}
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

float heldHouHadleyEdgeLatitudeDegrees(
    float equatorToPoleTemperatureContrastK,
    float troposphereHeightMetres,
    float referenceTemperatureK,
    float gravityMetresPerSecondSquared,
    float rotationRatePerSecond,
    float planetRadiusMetres)
{
    if (equatorToPoleTemperatureContrastK <= 0.0f || troposphereHeightMetres <= 0.0f ||
        referenceTemperatureK <= 0.0f || gravityMetresPerSecondSquared <= 0.0f ||
        planetRadiusMetres <= 0.0f)
    {
        return 0.0f;
    }

    if (rotationRatePerSecond <= 0.0f)
        return 90.0f;

    const float radiusRotation = rotationRatePerSecond * planetRadiusMetres;
    const float edgeRadiansSquared =
        (5.0f / 3.0f) *
        gravityMetresPerSecondSquared * troposphereHeightMetres /
        (radiusRotation * radiusRotation) *
        equatorToPoleTemperatureContrastK / referenceTemperatureK;
    return std::min(90.0f, std::sqrt(std::max(0.0f, edgeRadiansSquared)) * 180.0f / pi);
}

float thermalSurfacePressureAnomalyHpa(
    float temperatureAnomalyK,
    float referencePressureHpa,
    float referenceTemperatureK,
    float massRedistributionEfficiency)
{
    if (referencePressureHpa <= 0.0f || referenceTemperatureK <= 0.0f ||
        massRedistributionEfficiency <= 0.0f)
    {
        return 0.0f;
    }

    return -referencePressureHpa * temperatureAnomalyK /
        referenceTemperatureK * massRedistributionEfficiency;
}

float axisymmetricOverturningPressureAnomalyHpa(
    float latitudeDegrees,
    float thermalEquatorLatitudeDegrees,
    float hadleyHalfWidthDegrees,
    float pressureAmplitudeHpa)
{
    if (hadleyHalfWidthDegrees <= 0.0f || pressureAmplitudeHpa <= 0.0f)
        return 0.0f;

    const bool north = latitudeDegrees >= thermalEquatorLatitudeDegrees;
    const float poleLatitude = north ? 90.0f : -90.0f;
    const float hemisphereSpan = std::fabs(poleLatitude - thermalEquatorLatitudeDegrees);
    const float distance = std::clamp(
        std::fabs(latitudeDegrees - thermalEquatorLatitudeDegrees),
        0.0f,
        hemisphereSpan);
    const float hadleyEdge = std::min(hadleyHalfWidthDegrees, hemisphereSpan * 0.80f);
    const float subpolarLow = hadleyEdge + 0.5f * (hemisphereSpan - hadleyEdge);

    if (distance <= hadleyEdge)
    {
        return pressureAmplitudeHpa * smoothinterpolate(
            -0.55f,
            1.0f,
            distance / std::max(0.001f, hadleyEdge));
    }

    if (distance <= subpolarLow)
    {
        return pressureAmplitudeHpa * smoothinterpolate(
            1.0f,
            -0.85f,
            (distance - hadleyEdge) / std::max(0.001f, subpolarLow - hadleyEdge));
    }

    return pressureAmplitudeHpa * smoothinterpolate(
        -0.85f,
        0.45f,
        (distance - subpolarLow) / std::max(0.001f, hemisphereSpan - subpolarLow));
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

void setLastCirculationPrecisionDiagnostics(
    int season,
    const CirculationPrecisionDiagnostics& diagnostics)
{
    if (season < 0 || season >= CLIMATESEASONCOUNT)
        return;

    circulationPrecisionDiagnostics[season] = diagnostics;
}

const std::array<CirculationPrecisionDiagnostics, CLIMATESEASONCOUNT>&
lastCirculationPrecisionDiagnostics()
{
    return circulationPrecisionDiagnostics;
}
}
