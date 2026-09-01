#include "climate_atmosphere.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

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

float thermalModePressureAnomalyHpa(
    float temperatureAnomalyK,
    float airDensityKgM3,
    float lowerPressurePa,
    float upperPressurePa,
    float dryAirGasConstant)
{
    if (airDensityKgM3 <= 0.0f || lowerPressurePa <= upperPressurePa ||
        upperPressurePa <= 0.0f || dryAirGasConstant <= 0.0f)
    {
        return 0.0f;
    }

    const float geopotentialAnomaly = dryAirGasConstant *
        std::log(lowerPressurePa / upperPressurePa) * temperatureAnomalyK;
    return -airDensityKgM3 * geopotentialAnomaly / 100.0f;
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

HorizontalWind steadyQuadraticDragCoriolisWind(
    float forceEastMetresPerSecondSquared,
    float forceNorthMetresPerSecondSquared,
    float latitudeDegrees,
    float dragCoefficient,
    float boundaryLayerDepthMetres,
    float rotationRatePerSecond,
    float rotationDirection)
{
    if (dragCoefficient <= 0.0f || boundaryLayerDepthMetres <= 0.0f)
        return {};

    const float forceSquared =
        forceEastMetresPerSecondSquared * forceEastMetresPerSecondSquared +
        forceNorthMetresPerSecondSquared * forceNorthMetresPerSecondSquared;

    if (forceSquared <= 0.0f)
        return {};

    const float dragPerMetre = dragCoefficient / boundaryLayerDepthMetres;
    const float coriolis = coriolisParameterPerSecond(
        latitudeDegrees,
        rotationRatePerSecond,
        rotationDirection);
    const float coriolisSquared = coriolis * coriolis;
    const float speedSquared = 2.0f * forceSquared /
        (coriolisSquared + std::sqrt(
            coriolisSquared * coriolisSquared +
            4.0f * dragPerMetre * dragPerMetre * forceSquared));
    const float dragRate = dragPerMetre * std::sqrt(speedSquared);
    const float denominator = dragRate * dragRate + coriolisSquared;
    const float east =
        (dragRate * forceEastMetresPerSecondSquared + coriolis * forceNorthMetresPerSecondSquared) /
        denominator;
    const float north =
        (dragRate * forceNorthMetresPerSecondSquared - coriolis * forceEastMetresPerSecondSquared) /
        denominator;
    return { east, -north };
}

StationaryWaveResponse solveSteadyStationaryWavePressure(
    int longitudeCells,
    int latitudeCells,
    const std::vector<float>& equilibriumPressureAnomalyHpa,
    const std::vector<float>& dragTimeSeconds,
    float equivalentPressureDepthHpa,
    float pressureDampingTimeSeconds,
    float airDensityKgM3,
    float planetRadiusMetres,
    float rotationRatePerSecond,
    float rotationDirection,
    bool preserveZonalMean,
    int maximumIterations,
    float relativeTolerance)
{
    // Solve p' + tau_p P_e div(u(p')) = p'_eq with the steady
    // Rayleigh-Coriolis momentum balance used by the surface wind model.
    StationaryWaveResponse response;
    const size_t cellCount = static_cast<size_t>(std::max(0, longitudeCells)) *
        static_cast<size_t>(std::max(0, latitudeCells));
    response.pressureAnomalyHpa = equilibriumPressureAnomalyHpa;

    if (longitudeCells < 3 || latitudeCells < 3 ||
        equilibriumPressureAnomalyHpa.size() != cellCount ||
        dragTimeSeconds.size() != cellCount || equivalentPressureDepthHpa <= 0.0f ||
        pressureDampingTimeSeconds <= 0.0f || airDensityKgM3 <= 0.0f ||
        planetRadiusMetres <= 0.0f || maximumIterations <= 0 ||
        relativeTolerance <= 0.0f)
    {
        return response;
    }

    const auto index = [longitudeCells](int x, int y)
    {
        return static_cast<size_t>(y) * static_cast<size_t>(longitudeCells) +
            static_cast<size_t>(x);
    };
    const auto wrappedColumn = [longitudeCells](int x)
    {
        const int remainder = x % longitudeCells;
        return remainder < 0 ? remainder + longitudeCells : remainder;
    };
    const auto latitudeForRow = [latitudeCells](int y)
    {
        return 90.0 - 180.0 * static_cast<double>(y) /
            static_cast<double>(latitudeCells - 1);
    };
    const auto projectPressure = [&](std::vector<double>& values)
    {
        if (!preserveZonalMean)
        {
            for (int y : { 0, latitudeCells - 1 })
            {
                double poleMean = 0.0;
                for (int x = 0; x < longitudeCells; x++)
                    poleMean += values[index(x, y)];
                poleMean /= static_cast<double>(longitudeCells);
                for (int x = 0; x < longitudeCells; x++)
                    values[index(x, y)] = poleMean;
            }

            double weightedTotal = 0.0;
            double weightTotal = 0.0;
            for (int y = 0; y < latitudeCells; y++)
            {
                const double weight = std::max(
                    0.0,
                    std::cos(latitudeForRow(y) * static_cast<double>(pi) / 180.0));
                for (int x = 0; x < longitudeCells; x++)
                {
                    weightedTotal += values[index(x, y)] * weight;
                    weightTotal += weight;
                }
            }
            const double globalMean = weightTotal > 0.0 ?
                weightedTotal / weightTotal : 0.0;
            for (double& value : values)
                value -= globalMean;
            return;
        }

        for (int y = 0; y < latitudeCells; y++)
        {
            if (y == 0 || y == latitudeCells - 1)
            {
                for (int x = 0; x < longitudeCells; x++)
                    values[index(x, y)] = 0.0;
                continue;
            }

            double rowMean = 0.0;
            for (int x = 0; x < longitudeCells; x++)
                rowMean += values[index(x, y)];
            rowMean /= static_cast<double>(longitudeCells);
            for (int x = 0; x < longitudeCells; x++)
                values[index(x, y)] -= rowMean;
        }
    };
    const auto dot = [](const std::vector<double>& first, const std::vector<double>& second)
    {
        double total = 0.0;
        for (size_t cell = 0; cell < first.size(); cell++)
            total += first[cell] * second[cell];
        return total;
    };

    std::vector<double> eastWind(cellCount, 0.0);
    std::vector<double> southWind(cellCount, 0.0);
    std::vector<double> divergence(cellCount, 0.0);
    const auto applyOperator = [&](const std::vector<double>& pressure, std::vector<double>& result)
    {
        std::fill(eastWind.begin(), eastWind.end(), 0.0);
        std::fill(southWind.begin(), southWind.end(), 0.0);

        for (int y = 1; y < latitudeCells - 1; y++)
        {
            const double latitude = latitudeForRow(y);
            const double latitudeRadians = latitude * static_cast<double>(pi) / 180.0;
            const double cosine = std::max(0.02, std::fabs(std::cos(latitudeRadians)));
            const double zonalSpacing = 2.0 * static_cast<double>(pi) *
                static_cast<double>(planetRadiusMetres) * cosine /
                static_cast<double>(longitudeCells);
            const double meridionalSpacing = static_cast<double>(pi) *
                static_cast<double>(planetRadiusMetres) /
                static_cast<double>(latitudeCells - 1);
            const double coriolis = 2.0 * static_cast<double>(rotationRatePerSecond) *
                static_cast<double>(rotationDirection) * std::sin(latitudeRadians);

            for (int x = 0; x < longitudeCells; x++)
            {
                const size_t cell = index(x, y);
                const double dragTime = static_cast<double>(dragTimeSeconds[cell]);
                if (dragTime <= 0.0)
                    continue;

                const double pressureGradientEast =
                    (pressure[index(wrappedColumn(x + 1), y)] -
                        pressure[index(wrappedColumn(x - 1), y)]) * 100.0 /
                    (2.0 * zonalSpacing);
                const double pressureGradientNorth =
                    (pressure[index(x, y - 1)] - pressure[index(x, y + 1)]) * 100.0 /
                    (2.0 * meridionalSpacing);
                const double forceEast = -pressureGradientEast /
                    static_cast<double>(airDensityKgM3);
                const double forceNorth = -pressureGradientNorth /
                    static_cast<double>(airDensityKgM3);
                const double dragRate = 1.0 / dragTime;
                const double denominator = dragRate * dragRate + coriolis * coriolis;
                const double east =
                    (dragRate * forceEast + coriolis * forceNorth) / denominator;
                const double north =
                    (dragRate * forceNorth - coriolis * forceEast) / denominator;
                eastWind[cell] = east;
                southWind[cell] = -north;
            }
        }

        std::fill(divergence.begin(), divergence.end(), 0.0);
        for (int y = 1; y < latitudeCells - 1; y++)
        {
            const double latitude = latitudeForRow(y);
            const double latitudeRadians = latitude * static_cast<double>(pi) / 180.0;
            const double centreCosine = std::max(0.02, std::fabs(std::cos(latitudeRadians)));
            const double northCosine = std::max(
                0.0,
                std::cos(latitudeForRow(y - 1) * static_cast<double>(pi) / 180.0));
            const double southCosine = std::max(
                0.0,
                std::cos(latitudeForRow(y + 1) * static_cast<double>(pi) / 180.0));
            const double zonalSpacing = 2.0 * static_cast<double>(pi) *
                static_cast<double>(planetRadiusMetres) * centreCosine /
                static_cast<double>(longitudeCells);
            const double meridionalSpacing = static_cast<double>(pi) *
                static_cast<double>(planetRadiusMetres) /
                static_cast<double>(latitudeCells - 1);
            double rowMean = 0.0;

            for (int x = 0; x < longitudeCells; x++)
            {
                const double zonal =
                    (eastWind[index(wrappedColumn(x + 1), y)] -
                        eastWind[index(wrappedColumn(x - 1), y)]) /
                    (2.0 * zonalSpacing);
                const double meridional =
                    (southWind[index(x, y + 1)] * southCosine -
                        southWind[index(x, y - 1)] * northCosine) /
                    (2.0 * meridionalSpacing * centreCosine);
                divergence[index(x, y)] = zonal + meridional;
                rowMean += zonal + meridional;
            }

            if (preserveZonalMean)
            {
                rowMean /= static_cast<double>(longitudeCells);
                for (int x = 0; x < longitudeCells; x++)
                    divergence[index(x, y)] -= rowMean;
            }
        }

        result.resize(cellCount);
        const double coupling = static_cast<double>(pressureDampingTimeSeconds) *
            static_cast<double>(equivalentPressureDepthHpa);
        for (size_t cell = 0; cell < cellCount; cell++)
            result[cell] = pressure[cell] + coupling * divergence[cell];
        projectPressure(result);
    };

    // A local dissipative-response estimate is sufficient as a Jacobi
    // preconditioner; restarted GMRES handles the Coriolis asymmetry.
    std::vector<double> inverseDiagonal(cellCount, 1.0);
    const double coupling = static_cast<double>(pressureDampingTimeSeconds) *
        static_cast<double>(equivalentPressureDepthHpa);
    for (int y = 1; y < latitudeCells - 1; y++)
    {
        const double latitude = latitudeForRow(y);
        const double latitudeRadians = latitude * static_cast<double>(pi) / 180.0;
        const double cosine = std::max(0.02, std::fabs(std::cos(latitudeRadians)));
        const double zonalSpacing = 2.0 * static_cast<double>(pi) *
            static_cast<double>(planetRadiusMetres) * cosine /
            static_cast<double>(longitudeCells);
        const double meridionalSpacing = static_cast<double>(pi) *
            static_cast<double>(planetRadiusMetres) /
            static_cast<double>(latitudeCells - 1);
        const double coriolis = 2.0 * static_cast<double>(rotationRatePerSecond) *
            static_cast<double>(rotationDirection) * std::sin(latitudeRadians);
        for (int x = 0; x < longitudeCells; x++)
        {
            const size_t cell = index(x, y);
            const double dragTime = static_cast<double>(dragTimeSeconds[cell]);
            if (dragTime <= 0.0)
                continue;
            const double dragRate = 1.0 / dragTime;
            const double directPressureResponse =
                100.0 / static_cast<double>(airDensityKgM3) * dragRate /
                (dragRate * dragRate + coriolis * coriolis);
            const double approximateDiagonal = 1.0 + coupling * directPressureResponse *
                (0.5 / (zonalSpacing * zonalSpacing) +
                    0.5 / (meridionalSpacing * meridionalSpacing));
            inverseDiagonal[cell] = 1.0 / approximateDiagonal;
        }
    }
    const auto applyPreconditionedOperator = [&](
        const std::vector<double>& pressure,
        std::vector<double>& result)
    {
        applyOperator(pressure, result);
        for (size_t cell = 0; cell < cellCount; cell++)
            result[cell] *= inverseDiagonal[cell];
    };

    std::vector<double> rightHandSide(cellCount, 0.0);
    for (size_t cell = 0; cell < cellCount; cell++)
        rightHandSide[cell] = static_cast<double>(equilibriumPressureAnomalyHpa[cell]);
    projectPressure(rightHandSide);

    const double physicalRightHandSideNorm = std::sqrt(dot(rightHandSide, rightHandSide));
    if (physicalRightHandSideNorm <= std::numeric_limits<double>::epsilon())
    {
        response.pressureAnomalyHpa.assign(cellCount, 0.0f);
        response.converged = true;
        return response;
    }

    std::vector<double> solution = rightHandSide;
    std::vector<double> operatorSolution;
    std::vector<double> residual(cellCount, 0.0);
    constexpr int restartLength = 60;
    constexpr double breakdownTolerance = 1.0e-24;
    const int krylovColumns = std::min(restartLength, maximumIterations);
    std::vector<std::vector<double>> basis(
        krylovColumns + 1,
        std::vector<double>(cellCount, 0.0));
    std::vector<std::vector<double>> hessenberg(
        krylovColumns + 1,
        std::vector<double>(krylovColumns, 0.0));
    std::vector<double> cosine(krylovColumns, 0.0);
    std::vector<double> sine(krylovColumns, 0.0);
    std::vector<double> projectedResidual(krylovColumns + 1, 0.0);
    std::vector<double> work(cellCount, 0.0);

    while (response.iterations < maximumIterations)
    {
        applyOperator(solution, operatorSolution);
        for (size_t cell = 0; cell < cellCount; cell++)
            residual[cell] = rightHandSide[cell] - operatorSolution[cell];
        response.relativeResidual = static_cast<float>(
            std::sqrt(dot(residual, residual)) / physicalRightHandSideNorm);
        if (response.relativeResidual <= relativeTolerance)
        {
            response.converged = true;
            break;
        }

        for (size_t cell = 0; cell < cellCount; cell++)
            residual[cell] *= inverseDiagonal[cell];
        const double residualNorm = std::sqrt(dot(residual, residual));
        for (size_t cell = 0; cell < cellCount; cell++)
            basis[0][cell] = residual[cell] / residualNorm;
        for (auto& row : hessenberg)
            std::fill(row.begin(), row.end(), 0.0);
        std::fill(cosine.begin(), cosine.end(), 0.0);
        std::fill(sine.begin(), sine.end(), 0.0);
        std::fill(projectedResidual.begin(), projectedResidual.end(), 0.0);
        projectedResidual[0] = residualNorm;
        const int availableColumns = std::min(
            krylovColumns,
            maximumIterations - response.iterations);
        int usedColumns = 0;

        for (int column = 0; column < availableColumns; column++)
        {
            applyPreconditionedOperator(basis[column], work);
            for (int previous = 0; previous <= column; previous++)
            {
                hessenberg[previous][column] = dot(work, basis[previous]);
                for (size_t cell = 0; cell < cellCount; cell++)
                    work[cell] -= hessenberg[previous][column] * basis[previous][cell];
            }

            hessenberg[column + 1][column] = std::sqrt(dot(work, work));
            if (hessenberg[column + 1][column] > breakdownTolerance)
            {
                for (size_t cell = 0; cell < cellCount; cell++)
                    basis[column + 1][cell] =
                        work[cell] / hessenberg[column + 1][column];
            }

            for (int previous = 0; previous < column; previous++)
            {
                const double rotated =
                    cosine[previous] * hessenberg[previous][column] +
                    sine[previous] * hessenberg[previous + 1][column];
                hessenberg[previous + 1][column] =
                    -sine[previous] * hessenberg[previous][column] +
                    cosine[previous] * hessenberg[previous + 1][column];
                hessenberg[previous][column] = rotated;
            }

            const double rotationMagnitude = std::hypot(
                hessenberg[column][column],
                hessenberg[column + 1][column]);
            if (rotationMagnitude <= breakdownTolerance)
                break;
            cosine[column] = hessenberg[column][column] / rotationMagnitude;
            sine[column] = hessenberg[column + 1][column] / rotationMagnitude;
            hessenberg[column][column] = rotationMagnitude;
            hessenberg[column + 1][column] = 0.0;
            projectedResidual[column + 1] =
                -sine[column] * projectedResidual[column];
            projectedResidual[column] *= cosine[column];
            usedColumns = column + 1;
            response.iterations++;

        }

        if (usedColumns == 0)
            break;

        std::vector<double> coefficients(usedColumns, 0.0);
        for (int row = usedColumns - 1; row >= 0; row--)
        {
            double value = projectedResidual[row];
            for (int column = row + 1; column < usedColumns; column++)
                value -= hessenberg[row][column] * coefficients[column];
            if (std::fabs(hessenberg[row][row]) <= breakdownTolerance)
            {
                usedColumns = 0;
                break;
            }
            coefficients[row] = value / hessenberg[row][row];
        }

        if (usedColumns == 0)
            break;
        for (int column = 0; column < usedColumns; column++)
        {
            for (size_t cell = 0; cell < cellCount; cell++)
                solution[cell] += coefficients[column] * basis[column][cell];
        }
    }

    applyOperator(solution, operatorSolution);
    for (size_t cell = 0; cell < cellCount; cell++)
        residual[cell] = rightHandSide[cell] - operatorSolution[cell];
    response.relativeResidual = static_cast<float>(
        std::sqrt(dot(residual, residual)) / physicalRightHandSideNorm);
    response.converged = response.relativeResidual <= relativeTolerance;

    projectPressure(solution);
    response.pressureAnomalyHpa.resize(cellCount);
    for (size_t cell = 0; cell < cellCount; cell++)
        response.pressureAnomalyHpa[cell] = static_cast<float>(solution[cell]);
    return response;
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
