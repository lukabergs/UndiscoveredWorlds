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

std::vector<float> nonlocalThermalResponse(
    int longitudeCells,
    int latitudeCells,
    const std::vector<float>& localTemperatureAnomalyK,
    float tropicalLatitudeDegrees,
    float longitudinalReachDegrees,
    float meridionalReachDegrees,
    float rotationDirection)
{
    const size_t cellCount = static_cast<size_t>(std::max(0, longitudeCells)) *
        static_cast<size_t>(std::max(0, latitudeCells));
    if (longitudeCells <= 0 || latitudeCells <= 1 ||
        localTemperatureAnomalyK.size() != cellCount)
    {
        return localTemperatureAnomalyK;
    }

    const float longitudeDegreesPerCell = 360.0f / static_cast<float>(longitudeCells);
    const float latitudeDegreesPerCell = 180.0f / static_cast<float>(latitudeCells - 1);
    const int meridionalRadius = std::max(
        1,
        static_cast<int>(std::round(
            std::max(0.0f, meridionalReachDegrees) / latitudeDegreesPerCell)));
    const int longZonalRadius = std::max(
        1,
        static_cast<int>(std::round(
            std::max(0.0f, longitudinalReachDegrees) / longitudeDegreesPerCell)));
    const int shortZonalRadius = std::max(1, longZonalRadius / 4);
    const int symmetricZonalRadius = std::max(1, longZonalRadius / 3);
    const int rotationSign = rotationDirection >= 0.0f ? 1 : -1;

    const auto index = [longitudeCells](int x, int y)
    {
        return static_cast<size_t>(y) * longitudeCells + x;
    };
    const auto wrappedColumn = [longitudeCells](int x)
    {
        const int remainder = x % longitudeCells;
        return remainder < 0 ? remainder + longitudeCells : remainder;
    };

    std::vector<float> meridional(cellCount, 0.0f);
    for (int y = 0; y < latitudeCells; y++)
    {
        const int firstRow = std::max(0, y - meridionalRadius);
        const int lastRow = std::min(latitudeCells - 1, y + meridionalRadius);
        for (int x = 0; x < longitudeCells; x++)
        {
            double weighted = 0.0;
            double weightTotal = 0.0;
            for (int sourceY = firstRow; sourceY <= lastRow; sourceY++)
            {
                const float distance = static_cast<float>(std::abs(sourceY - y));
                const double weight = std::exp(
                    -2.0 * distance / static_cast<double>(meridionalRadius));
                weighted += weight * localTemperatureAnomalyK[index(x, sourceY)];
                weightTotal += weight;
            }
            meridional[index(x, y)] = weightTotal > 0.0
                ? static_cast<float>(weighted / weightTotal)
                : localTemperatureAnomalyK[index(x, y)];
        }
    }

    std::vector<float> response(cellCount, 0.0f);
    const float tropicalCore = std::max(1.0f, std::fabs(tropicalLatitudeDegrees));
    for (int y = 0; y < latitudeCells; y++)
    {
        const float latitude = 90.0f - 180.0f * static_cast<float>(y) /
            static_cast<float>(latitudeCells - 1);
        const float tropicalBlend = 1.0f - smoothstep(
            (std::fabs(latitude) - tropicalCore) / (0.5f * tropicalCore));
        double sourceRowMean = 0.0;
        double responseRowMean = 0.0;

        for (int x = 0; x < longitudeCells; x++)
        {
            double asymmetric = 0.0;
            double asymmetricWeight = 0.0;
            for (int offset = -shortZonalRadius; offset <= longZonalRadius; offset++)
            {
                const int scale = offset < 0 ? shortZonalRadius : longZonalRadius;
                const double weight = std::exp(
                    -2.0 * static_cast<double>(std::abs(offset)) /
                        static_cast<double>(scale));
                const int sourceX = wrappedColumn(x + rotationSign * offset);
                asymmetric += weight * meridional[index(sourceX, y)];
                asymmetricWeight += weight;
            }

            double symmetric = 0.0;
            double symmetricWeight = 0.0;
            for (int offset = -symmetricZonalRadius;
                 offset <= symmetricZonalRadius;
                 offset++)
            {
                const double weight = std::exp(
                    -2.0 * static_cast<double>(std::abs(offset)) /
                        static_cast<double>(symmetricZonalRadius));
                symmetric += weight * meridional[index(wrappedColumn(x + offset), y)];
                symmetricWeight += weight;
            }

            const float asymmetricValue = asymmetricWeight > 0.0
                ? static_cast<float>(asymmetric / asymmetricWeight)
                : meridional[index(x, y)];
            const float symmetricValue = symmetricWeight > 0.0
                ? static_cast<float>(symmetric / symmetricWeight)
                : meridional[index(x, y)];
            response[index(x, y)] =
                tropicalBlend * asymmetricValue + (1.0f - tropicalBlend) * symmetricValue;
            sourceRowMean += meridional[index(x, y)];
            responseRowMean += response[index(x, y)];
        }

        const float meanCorrection = static_cast<float>(
            (responseRowMean - sourceRowMean) / static_cast<double>(longitudeCells));
        for (int x = 0; x < longitudeCells; x++)
            response[index(x, y)] -= meanCorrection;
    }

    return response;
}

std::vector<float> mechanicalTopographicPressureForcingHpa(
    int longitudeCells,
    int latitudeCells,
    const std::vector<float>& terrainMetres,
    const std::vector<float>& eastWindMps,
    const std::vector<float>& southWindMps,
    float sampleDistanceCells,
    float terrainScaleMetres,
    float maximumAmplitudeHpa,
    float minimumWindMps,
    float fullStrengthWindMps,
    float minimumLatitudeDegrees,
    float fullStrengthLatitudeDegrees)
{
    const size_t cellCount = static_cast<size_t>(std::max(0, longitudeCells)) *
        static_cast<size_t>(std::max(0, latitudeCells));
    if (longitudeCells <= 0 || latitudeCells <= 1 || terrainMetres.size() != cellCount ||
        eastWindMps.size() != cellCount || southWindMps.size() != cellCount ||
        sampleDistanceCells <= 0.0f || terrainScaleMetres <= 0.0f ||
        maximumAmplitudeHpa <= 0.0f)
    {
        return std::vector<float>(cellCount, 0.0f);
    }

    const auto index = [longitudeCells](int x, int y)
    {
        return static_cast<size_t>(y) * longitudeCells + x;
    };
    const auto wrappedColumn = [longitudeCells](int x)
    {
        const int remainder = x % longitudeCells;
        return remainder < 0 ? remainder + longitudeCells : remainder;
    };
    const auto sampleTerrain = [&](float x, float y)
    {
        const int y0 = std::clamp(static_cast<int>(std::floor(y)), 0, latitudeCells - 1);
        const int y1 = std::min(latitudeCells - 1, y0 + 1);
        const int x0Unwrapped = static_cast<int>(std::floor(x));
        const int x0 = wrappedColumn(x0Unwrapped);
        const int x1 = wrappedColumn(x0Unwrapped + 1);
        const float fractionX = x - std::floor(x);
        const float fractionY = std::clamp(y - static_cast<float>(y0), 0.0f, 1.0f);
        const float north = terrainMetres[index(x0, y0)] +
            (terrainMetres[index(x1, y0)] - terrainMetres[index(x0, y0)]) * fractionX;
        const float south = terrainMetres[index(x0, y1)] +
            (terrainMetres[index(x1, y1)] - terrainMetres[index(x0, y1)]) * fractionX;
        return north + (south - north) * fractionY;
    };

    std::vector<float> forcing(cellCount, 0.0f);
    double areaWeightedTotal = 0.0;
    double areaWeightTotal = 0.0;
    for (int y = 0; y < latitudeCells; y++)
    {
        const float latitude = 90.0f - 180.0f * static_cast<float>(y) /
            static_cast<float>(latitudeCells - 1);
        const float latitudeDenominator = std::max(
            0.001f,
            fullStrengthLatitudeDegrees - minimumLatitudeDegrees);
        const float latitudeFactor = smoothstep(
            (std::fabs(latitude) - minimumLatitudeDegrees) / latitudeDenominator);
        const double areaWeight = std::max(
            0.0,
            static_cast<double>(std::cos(latitude * pi / 180.0f)));

        for (int x = 0; x < longitudeCells; x++)
        {
            const size_t cell = index(x, y);
            areaWeightTotal += areaWeight;
            const float east = eastWindMps[cell];
            const float south = southWindMps[cell];
            const float speed = std::sqrt(east * east + south * south);
            const float windDenominator = std::max(0.001f, fullStrengthWindMps - minimumWindMps);
            const float windFactor = smoothstep((speed - minimumWindMps) / windDenominator);
            if (windFactor <= 0.0f || latitudeFactor <= 0.0f)
                continue;

            const float directionX = east / speed;
            const float directionY = south / speed;
            const float terrainHere = terrainMetres[cell];
            const float terrainAhead = sampleTerrain(
                static_cast<float>(x) + directionX * sampleDistanceCells,
                static_cast<float>(y) + directionY * sampleDistanceCells);
            const float terrainUpwind = sampleTerrain(
                static_cast<float>(x) - directionX * sampleDistanceCells,
                static_cast<float>(y) - directionY * sampleDistanceCells);
            const float windwardRise = std::max(0.0f, terrainAhead - terrainHere);
            const float leeDrop = std::max(0.0f, terrainUpwind - terrainHere);
            forcing[cell] = std::clamp(
                (windwardRise - leeDrop) / terrainScaleMetres *
                    maximumAmplitudeHpa * windFactor * latitudeFactor,
                -maximumAmplitudeHpa,
                maximumAmplitudeHpa);
            areaWeightedTotal += areaWeight * forcing[cell];
        }
    }

    const float mean = areaWeightTotal > 0.0
        ? static_cast<float>(areaWeightedTotal / areaWeightTotal)
        : 0.0f;
    for (float& value : forcing)
        value -= mean;
    return forcing;
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
