#include "thermal_response.hpp"

#include <algorithm>
#include <cmath>

namespace climateatmosphere
{
namespace
{
float smoothstep(float value)
{
    const float clamped = std::clamp(value, 0.0f, 1.0f);
    return clamped * clamped * (3.0f - 2.0f * clamped);
}

std::vector<float> smoothMeridionalAnomalies(int longitudeCells, int latitudeCells,
    const std::vector<float>& localTemperatureAnomalyK, int meridionalRadius)
{
    const size_t cellCount = localTemperatureAnomalyK.size();
    const auto index = [longitudeCells](int x, int y)
    {
        return static_cast<size_t>(y) * longitudeCells + x;
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
    return meridional;
}

std::vector<float> spreadZonalAnomalies(int longitudeCells, int latitudeCells,
    const std::vector<float>& meridional, float tropicalLatitudeDegrees,
    int longZonalRadius, int shortZonalRadius, int symmetricZonalRadius, int rotationSign)
{
    const size_t cellCount = meridional.size();
    const auto index = [longitudeCells](int x, int y)
    {
        return static_cast<size_t>(y) * longitudeCells + x;
    };
    const auto wrappedColumn = [longitudeCells](int x)
    {
        const int remainder = x % longitudeCells;
        return remainder < 0 ? remainder + longitudeCells : remainder;
    };
    std::vector<float> response(cellCount, 0.0f);
    const float tropicalCore = std::max(1.0f, std::fabs(tropicalLatitudeDegrees));
    for (int y = 0; y < latitudeCells; y++)
    {
        const float latitude = 90.0f - 180.0f *
            (static_cast<float>(y) + 0.5f) / static_cast<float>(latitudeCells);
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
    const float latitudeDegreesPerCell = 180.0f / static_cast<float>(latitudeCells);
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

    const auto meridional = smoothMeridionalAnomalies(longitudeCells, latitudeCells,
        localTemperatureAnomalyK, meridionalRadius);
    return spreadZonalAnomalies(longitudeCells, latitudeCells, meridional,
        tropicalLatitudeDegrees, longZonalRadius, shortZonalRadius,
        symmetricZonalRadius, rotationSign);
}
}
