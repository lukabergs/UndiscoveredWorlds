#include "climate_atmosphere.hpp"
#include "climate_grid.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
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

void filterStationary(int columns, int rows, int zonalModes, int meridionalModes,
    std::vector<float>& field)
{
    // Fourier in longitude, cell-centred cosine modes in latitude. Projection
    // is applied to the source, never to an already accepted solver residual.
    std::vector<float> filtered(field.size(), 0.0f);
    for (int y = 0; y < rows; ++y)
        for (int k = 1; k <= std::min(zonalModes, (columns - 1) / 2); ++k)
        {
            double a = 0.0, b = 0.0;
            for (int x = 0; x < columns; ++x)
            {
                const double angle = 2.0 * pi * k * x / columns;
                a += field[y * columns + x] * std::cos(angle) * 2.0 / columns;
                b += field[y * columns + x] * std::sin(angle) * 2.0 / columns;
            }
            for (int x = 0; x < columns; ++x)
                filtered[y * columns + x] += static_cast<float>(
                    a * std::cos(2.0 * pi * k * x / columns) +
                    b * std::sin(2.0 * pi * k * x / columns));
        }
    std::fill(field.begin(), field.end(), 0.0f);
    for (int x = 0; x < columns; ++x)
        for (int k = 0; k <= std::min(meridionalModes, rows - 1); ++k)
        {
            double coefficient = 0.0;
            for (int y = 0; y < rows; ++y)
                coefficient += filtered[y * columns + x] *
                    std::cos(pi * k * (y + 0.5) / rows) * (k == 0 ? 1.0 : 2.0) / rows;
            for (int y = 0; y < rows; ++y)
                field[y * columns + x] += static_cast<float>(coefficient *
                    std::cos(pi * k * (y + 0.5) / rows));
        }
}
}

ColumnHeating diagnoseColumnHeating(const ColumnHeatingInput& in)
{
    // Flux-divergence discretization of grey two-stream transfer; see GFDL's
    // idealized moist model (Frierson et al. 2006, doi:10.1175/JAS3753.1).
    // Optical depths are reduced-model parameters, not observed Earth fields.
    constexpr double sigma = 5.670374419e-8, latentHeat = 2.5e6;
    ColumnHeating out;
    std::array<double, 2> transmission{}, emission{}, solarTransmission{};
    for (int layer = 0; layer < 2; ++layer)
    {
        transmission[layer] = std::exp(-std::max(0.0, in.longwaveOpticalDepth[layer]));
        solarTransmission[layer] = std::exp(-std::max(0.0, in.shortwaveOpticalDepth[layer]));
        emission[layer] = (1.0 - transmission[layer]) * sigma * std::pow(std::max(1.0, in.airTemperatureK[layer]), 4);
    }
    const double surfaceEmission = sigma * std::pow(std::max(1.0, in.surfaceTemperatureK), 4);
    const double downwardMiddle = emission[1];
    const double downwardSurface = downwardMiddle * transmission[0] + emission[0];
    const double upwardMiddle = surfaceEmission * transmission[0] + emission[0];
    const double outgoingLongwave = upwardMiddle * transmission[1] + emission[1];
    out.radiativeWm2[0] = surfaceEmission + downwardMiddle - upwardMiddle - downwardSurface;
    out.radiativeWm2[1] = upwardMiddle - outgoingLongwave - downwardMiddle;
    const double solarTop = std::max(0.0, in.incomingSolarWm2);
    const double solarMiddle = solarTop * solarTransmission[1];
    const double solarSurface = solarMiddle * solarTransmission[0];
    const double reflectedSurface = solarSurface * std::clamp(in.surfaceAlbedo, 0.0, 1.0);
    const double reflectedMiddle = reflectedSurface * solarTransmission[0];
    const double reflectedTop = reflectedMiddle * solarTransmission[1];
    out.radiativeWm2[0] += solarMiddle - solarSurface + reflectedSurface - reflectedMiddle;
    out.radiativeWm2[1] += solarTop - solarMiddle + reflectedMiddle - reflectedTop;
    out.surfaceRadiativeWm2 = solarSurface - reflectedSurface + downwardSurface - surfaceEmission;
    out.topNetRadiationWm2 = solarTop - reflectedTop - outgoingLongwave;
    const double latentPerMm = latentHeat / std::max(1.0, in.accumulationSeconds);
    out.latentWm2[0] = latentPerMm * (in.condensationMm[0] - in.reevaporationMm);
    out.latentWm2[1] = latentPerMm * in.condensationMm[1];
    out.totalWm2[0] = out.radiativeWm2[0] + out.latentWm2[0] + in.sensibleHeatingWm2;
    out.totalWm2[1] = out.radiativeWm2[1] + out.latentWm2[1];
    out.surfaceNetHeatingWm2 = out.surfaceRadiativeWm2 - in.sensibleHeatingWm2 - latentPerMm * in.surfaceEvaporationMm;
    out.vapourLatentStorageWm2 = latentPerMm * (in.surfaceEvaporationMm + in.reevaporationMm -
        in.condensationMm[0] - in.condensationMm[1]);
    out.closureResidualWm2 = out.surfaceNetHeatingWm2 + out.totalWm2[0] + out.totalWm2[1] +
        out.vapourLatentStorageWm2 - out.topNetRadiationWm2;
    return out;
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
    // Reduced eddy-cell closure: the subpolar low lies poleward of the
    // midpoint between the Hadley edge and the pole.
    const float subpolarLow = hadleyEdge + 0.65f * (hemisphereSpan - hadleyEdge);

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

DiabaticHeatingBudget diagnoseDiabaticHeating(
    int longitudeCells,
    int latitudeCells,
    const std::vector<float>& absorbedShortwaveWm2,
    const std::vector<float>& outgoingLongwaveWm2,
    const std::vector<float>& sensibleHeatingWm2,
    const std::vector<float>& condensationMm,
    float accumulationSeconds,
    float verticalProjection,
    float latentHeatCouplingFraction,
    float latentHeatJkg,
    bool removeZonalMean)
{
    DiabaticHeatingBudget budget;
    const std::size_t cellCount = static_cast<std::size_t>(std::max(0, longitudeCells)) *
        static_cast<std::size_t>(std::max(0, latitudeCells));
    if (longitudeCells <= 0 || latitudeCells <= 0 ||
        absorbedShortwaveWm2.size() != cellCount ||
        outgoingLongwaveWm2.size() != cellCount ||
        sensibleHeatingWm2.size() != cellCount || condensationMm.size() != cellCount ||
        accumulationSeconds <= 0.0f || latentHeatJkg < 0.0f)
    {
        return budget;
    }

    budget.netColumnHeatingWm2.resize(cellCount, 0.0f);
    budget.stationaryProjectedHeatingWm2.resize(cellCount, 0.0f);
    const double projection = std::clamp(static_cast<double>(verticalProjection), 0.0, 1.0);
    const double latentFraction = std::clamp(
        static_cast<double>(latentHeatCouplingFraction), 0.0, 1.0);
    double weightTotal = 0.0;
    for (int y = 0; y < latitudeCells; y++)
    {
        const double latitude = pi * 0.5 -
            (static_cast<double>(y) + 0.5) * pi / static_cast<double>(latitudeCells);
        const double weight = std::max(0.0, std::cos(latitude));
        for (int x = 0; x < longitudeCells; x++)
        {
            const std::size_t cell = static_cast<std::size_t>(y) * longitudeCells + x;
            const double radiative = static_cast<double>(absorbedShortwaveWm2[cell]) -
                outgoingLongwaveWm2[cell];
            const double sensible = sensibleHeatingWm2[cell];
            const double latent = static_cast<double>(condensationMm[cell]) * latentHeatJkg /
                static_cast<double>(accumulationSeconds) * latentFraction;
            const double net = radiative + sensible + latent;
            budget.netColumnHeatingWm2[cell] = static_cast<float>(net);
            budget.stationaryProjectedHeatingWm2[cell] = static_cast<float>(net * projection);
            budget.areaWeightedRadiativeHeatingWm2 += weight * radiative;
            budget.areaWeightedSensibleHeatingWm2 += weight * sensible;
            budget.areaWeightedLatentHeatingWm2 += weight * latent;
            budget.areaWeightedProjectedHeatingWm2 += weight * net * projection;
            weightTotal += weight;
        }
    }

    if (weightTotal > 0.0)
    {
        budget.areaWeightedRadiativeHeatingWm2 /= weightTotal;
        budget.areaWeightedSensibleHeatingWm2 /= weightTotal;
        budget.areaWeightedLatentHeatingWm2 /= weightTotal;
        budget.areaWeightedProjectedHeatingWm2 /= weightTotal;
    }

    if (removeZonalMean)
    {
        for (int y = 0; y < latitudeCells; y++)
        {
            double rowMean = 0.0;
            for (int x = 0; x < longitudeCells; x++)
                rowMean += budget.stationaryProjectedHeatingWm2[
                    static_cast<std::size_t>(y) * longitudeCells + x];
            rowMean /= static_cast<double>(longitudeCells);
            double residualMean = 0.0;
            for (int x = 0; x < longitudeCells; x++)
            {
                float& value = budget.stationaryProjectedHeatingWm2[
                    static_cast<std::size_t>(y) * longitudeCells + x];
                value -= static_cast<float>(rowMean);
                residualMean += value;
            }
            budget.maximumAbsoluteRowMeanWm2 = std::max(
                budget.maximumAbsoluteRowMeanWm2,
                static_cast<float>(std::abs(residualMean / longitudeCells)));
        }
    }
    return budget;
}

StationaryParameterDiagnosis diagnoseStationaryParameters(
    float bruntVaisalaFrequencyPerSecond,
    float modeDepthMetres,
    float gravityMetresPerSecondSquared,
    float planetRadiusMetres,
    float rotationRatePerSecond,
    int longitudeCells,
    int latitudeCells,
    float nondimensionalDamping,
    float resolvedForcingScaleMetres)
{
    StationaryParameterDiagnosis diagnosis;
    if (bruntVaisalaFrequencyPerSecond <= 0.0f || modeDepthMetres <= 0.0f ||
        gravityMetresPerSecondSquared <= 0.0f || planetRadiusMetres <= 0.0f ||
        longitudeCells < 3 || latitudeCells < 3)
    {
        return diagnosis;
    }

    diagnosis.gravityWaveSpeedMps = bruntVaisalaFrequencyPerSecond * modeDepthMetres / pi;
    diagnosis.equivalentDepthMetres =
        diagnosis.gravityWaveSpeedMps * diagnosis.gravityWaveSpeedMps /
        gravityMetresPerSecondSquared;
    diagnosis.adjustmentLengthMetres = diagnosis.gravityWaveSpeedMps /
        std::max(1.0e-8f, 2.0f * std::abs(rotationRatePerSecond));
    diagnosis.dampingTimeSeconds = diagnosis.adjustmentLengthMetres /
        diagnosis.gravityWaveSpeedMps /
        std::max(0.01f, nondimensionalDamping);
    const float forcingScale = std::max(
        resolvedForcingScaleMetres,
        2.0f * pi * planetRadiusMetres / static_cast<float>(longitudeCells));
    diagnosis.maximumZonalWavenumber = std::clamp(
        static_cast<int>(std::floor(2.0f * pi * planetRadiusMetres / forcingScale)),
        1,
        longitudeCells / 2);
    diagnosis.maximumMeridionalWavenumber = std::clamp(
        static_cast<int>(std::floor(pi * planetRadiusMetres / forcingScale)),
        1,
        latitudeCells - 1);
    return diagnosis;
}

float diagnoseBruntVaisalaFrequency(
    float temperatureK, float lapseRateKPerMetre, float gravityMps2)
{
    if (temperatureK <= 0.0f || gravityMps2 <= 0.0f)
        return 0.0f;
    return std::sqrt(std::max(0.0f, gravityMps2 / temperatureK *
        (gravityMps2 / 1004.0f - lapseRateKPerMetre)));
}

std::vector<float> upperOrographicHeightForcing(
    int columns, int rows, const std::vector<float>& terrainMetres,
    const std::vector<float>& backgroundEastWindMps,
    float stabilityPerSecond, float levelHeightMetres,
    float dampingTimeSeconds, const ModeSeparatedCirculationConfig& config)
{
    const std::size_t count = static_cast<std::size_t>(std::max(0, columns)) * std::max(0, rows);
    std::vector<float> result(count, 0.0f);
    if (columns < 3 || rows < 2 || terrainMetres.size() != count ||
        backgroundEastWindMps.size() != count || stabilityPerSecond <= 0.0f ||
        dampingTimeSeconds <= 0.0f || levelHeightMetres < 0.0f ||
        config.planetRadiusMetres <= 0.0f || config.gravityMetresPerSecondSquared <= 0.0f)
        return result;
    const auto grid = climategrid::makeSphericalGrid(columns, rows, config.planetRadiusMetres);
    for (int y = 0; y < rows; ++y)
    {
        double wind = 0.0;
        for (int x = 0; x < columns; ++x)
            wind += backgroundEastWindMps[y * columns + x] / columns;
        if (std::abs(wind) < 0.1) continue;
        const double f = 2.0 * config.rotationRatePerSecond * config.rotationDirection *
            std::sin(grid.latitudeCentresRadians[y]);
        for (int mode = 1; mode <= std::min(config.upperMaximumZonalWavenumber, (columns - 1) / 2); ++mode)
        {
            const double k = mode / (config.planetRadiusMetres * std::cos(grid.latitudeCentresRadians[y]));
            // Complex intrinsic frequency damps critical levels continuously.
            const std::complex<double> omega(k * wind, -1.0 / dampingTimeSeconds);
            const auto mSquared = k * k * (static_cast<double>(stabilityPerSecond) * stabilityPerSecond - omega * omega) /
                (omega * omega - f * f);
            auto m = std::sqrt(mSquared);
            if (m.imag() < 0.0) m = -m;
            std::complex<double> terrain(0.0, 0.0);
            for (int x = 0; x < columns; ++x)
                terrain += static_cast<double>(terrainMetres[y * columns + x]) *
                    std::exp(std::complex<double>(0.0, -2.0 * pi * mode * x / columns)) *
                    (2.0 / columns);
            const auto response = terrain * (omega * omega - f * f) /
                (config.gravityMetresPerSecondSquared * k * k) *
                std::complex<double>(0.0, 1.0) * m *
                std::exp(std::complex<double>(0.0, 1.0) * m * static_cast<double>(levelHeightMetres));
            for (int x = 0; x < columns; ++x)
                result[y * columns + x] += static_cast<float>((response *
                    std::exp(std::complex<double>(0.0, 2.0 * pi * mode * x / columns))).real());
        }
    }
    return result;
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
        const float latitude = 90.0f - 180.0f *
            (static_cast<float>(y) + 0.5f) / static_cast<float>(latitudeCells);
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
    float relativeTolerance,
    int restartLength)
{
    // Solve p' + tau_p P_e div(u(p')) = p'_eq with the steady
    // Rayleigh-Coriolis momentum balance used by the surface wind model.
    StationaryWaveResponse response;
    const size_t cellCount = static_cast<size_t>(std::max(0, longitudeCells)) *
        static_cast<size_t>(std::max(0, latitudeCells));
    response.pressureAnomalyHpa = equilibriumPressureAnomalyHpa;
    response.equilibriumPressureAnomalyHpa = equilibriumPressureAnomalyHpa;

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
        return 90.0 - 180.0 * (static_cast<double>(y) + 0.5) /
            static_cast<double>(latitudeCells);
    };
    const auto projectPressure = [&](std::vector<double>& values)
    {
        if (!preserveZonalMean)
        {
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

        for (int y = 0; y < latitudeCells; y++)
        {
            const int northRow = std::max(0, y - 1);
            const int southRow = std::min(latitudeCells - 1, y + 1);
            const double latitude = latitudeForRow(y);
            const double latitudeRadians = latitude * static_cast<double>(pi) / 180.0;
            const double cosine = std::max(0.02, std::fabs(std::cos(latitudeRadians)));
            const double zonalSpacing = 2.0 * static_cast<double>(pi) *
                static_cast<double>(planetRadiusMetres) * cosine /
                static_cast<double>(longitudeCells);
            const double meridionalSpacing = static_cast<double>(pi) *
                static_cast<double>(planetRadiusMetres) /
                static_cast<double>(latitudeCells);
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
                    (pressure[index(x, northRow)] - pressure[index(x, southRow)]) * 100.0 /
                    (static_cast<double>(southRow - northRow) * meridionalSpacing);
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
        for (int y = 0; y < latitudeCells; y++)
        {
            const int northRow = std::max(0, y - 1);
            const int southRow = std::min(latitudeCells - 1, y + 1);
            const double latitude = latitudeForRow(y);
            const double latitudeRadians = latitude * static_cast<double>(pi) / 180.0;
            const double centreCosine = std::max(0.02, std::fabs(std::cos(latitudeRadians)));
            const double northCosine = std::max(
                0.0,
                std::cos(latitudeForRow(northRow) * static_cast<double>(pi) / 180.0));
            const double southCosine = std::max(
                0.0,
                std::cos(latitudeForRow(southRow) * static_cast<double>(pi) / 180.0));
            const double zonalSpacing = 2.0 * static_cast<double>(pi) *
                static_cast<double>(planetRadiusMetres) * centreCosine /
                static_cast<double>(longitudeCells);
            const double meridionalSpacing = static_cast<double>(pi) *
                static_cast<double>(planetRadiusMetres) /
                static_cast<double>(latitudeCells);
            double rowMean = 0.0;

            for (int x = 0; x < longitudeCells; x++)
            {
                const double zonal =
                    (eastWind[index(wrappedColumn(x + 1), y)] -
                        eastWind[index(wrappedColumn(x - 1), y)]) /
                    (2.0 * zonalSpacing);
                const double meridional =
                    (southWind[index(x, southRow)] * southCosine -
                        southWind[index(x, northRow)] * northCosine) /
                    (static_cast<double>(southRow - northRow) *
                        meridionalSpacing * centreCosine);
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
    for (int y = 0; y < latitudeCells; y++)
    {
        const double latitude = latitudeForRow(y);
        const double latitudeRadians = latitude * static_cast<double>(pi) / 180.0;
        const double cosine = std::max(0.02, std::fabs(std::cos(latitudeRadians)));
        const double zonalSpacing = 2.0 * static_cast<double>(pi) *
            static_cast<double>(planetRadiusMetres) * cosine /
            static_cast<double>(longitudeCells);
        const double meridionalSpacing = static_cast<double>(pi) *
            static_cast<double>(planetRadiusMetres) /
            static_cast<double>(latitudeCells);
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
    restartLength = std::clamp(restartLength, 1, maximumIterations);
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
    float previousCycleResidual = std::numeric_limits<float>::infinity();
    int stagnantCycles = 0;

    while (response.iterations < maximumIterations)
    {
        response.restartCycles++;
        applyOperator(solution, operatorSolution);
        for (size_t cell = 0; cell < cellCount; cell++)
            residual[cell] = rightHandSide[cell] - operatorSolution[cell];
        response.relativeResidual = static_cast<float>(
            std::sqrt(dot(residual, residual)) / physicalRightHandSideNorm);
        response.residualHistory.push_back(response.relativeResidual);
        if (std::isfinite(previousCycleResidual))
        {
            const float improvement = previousCycleResidual - response.relativeResidual;
            stagnantCycles = improvement <= previousCycleResidual * 1.0e-3f
                ? stagnantCycles + 1
                : 0;
            response.stagnated = stagnantCycles >= 3;
        }
        previousCycleResidual = response.relativeResidual;
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
    if (response.residualHistory.empty() ||
        response.residualHistory.back() != response.relativeResidual)
    {
        response.residualHistory.push_back(response.relativeResidual);
    }
    response.converged = response.relativeResidual <= relativeTolerance;

    projectPressure(solution);
    response.pressureAnomalyHpa.resize(cellCount);
    for (size_t cell = 0; cell < cellCount; cell++)
        response.pressureAnomalyHpa[cell] = static_cast<float>(solution[cell]);
    return response;
}

ModeSeparatedCirculation solveModeSeparatedCirculation(
    int longitudeCells,
    int latitudeCells,
    const std::vector<float>& zonalSurfacePressureHpa,
    const std::vector<float>& stationaryHeatingWm2,
    const std::vector<float>& orographicSurfaceForcingHpa,
    const ModeSeparatedCirculationConfig& config)
{
    ModeSeparatedCirculation result;
    const std::size_t cellCount = static_cast<std::size_t>(std::max(0, longitudeCells)) *
        static_cast<std::size_t>(std::max(0, latitudeCells));
    if (longitudeCells < 3 || latitudeCells < 3 ||
        zonalSurfacePressureHpa.size() != static_cast<std::size_t>(latitudeCells) ||
        stationaryHeatingWm2.size() != cellCount ||
        orographicSurfaceForcingHpa.size() != cellCount ||
        config.airDensityKgM3 <= 0.0f || config.gravityMetresPerSecondSquared <= 0.0f ||
        config.planetRadiusMetres <= 0.0f)
    {
        return result;
    }

    result.surfacePressureAnomalyHpa.assign(cellCount, 0.0f);
    result.upperHeightAnomalyMetres.assign(cellCount, 0.0f);
    result.surfaceEastWindMps.assign(cellCount, 0.0f);
    result.surfaceSouthWindMps.assign(cellCount, 0.0f);
    result.upperEastWindMps.assign(cellCount, 0.0f);
    result.upperSouthWindMps.assign(cellCount, 0.0f);
    std::vector<float> surfaceEquilibrium(cellCount, 0.0f);
    std::vector<float> upperEquilibriumPressure(cellCount, 0.0f);
    std::vector<float> surfaceDrag(cellCount, config.surfaceDragTimeSeconds);
    if (config.surfaceDragTimesSeconds.size() == cellCount)
        surfaceDrag = config.surfaceDragTimesSeconds;
    std::vector<float> upperDrag(cellCount, config.upperDragTimeSeconds);
    const float upperPressurePerMetre = config.airDensityKgM3 *
        config.gravityMetresPerSecondSquared / 100.0f;

    for (int y = 0; y < latitudeCells; y++)
    {
        for (int x = 0; x < longitudeCells; x++)
        {
            const std::size_t cell = static_cast<std::size_t>(y) * longitudeCells + x;
            const float surfaceStationary = config.enabled.stationary
                ? stationaryHeatingWm2[cell] *
                        config.surfaceHeatingPressureResponseHpaPerWm2 +
                    orographicSurfaceForcingHpa[cell]
                : 0.0f;
            surfaceEquilibrium[cell] = surfaceStationary;
            const float upperHeight = config.enabled.stationary
                ? (config.upperStationaryHeatingWm2.size() == cellCount
                    ? config.upperStationaryHeatingWm2[cell] : stationaryHeatingWm2[cell]) *
                        config.upperHeatingHeightResponseMetresPerWm2 +
                    (config.upperOrographicHeightMetres.size() == cellCount
                        ? config.upperOrographicHeightMetres[cell] : 0.0f)
                : 0.0f;
            upperEquilibriumPressure[cell] = upperHeight * upperPressurePerMetre +
                (config.enabled.surface ? surfaceStationary * config.surfaceToUpperCoupling : 0.0f);
        }
    }

    filterStationary(longitudeCells, latitudeCells, config.maximumZonalWavenumber,
        config.maximumMeridionalWavenumber, surfaceEquilibrium);
    filterStationary(longitudeCells, latitudeCells, config.upperMaximumZonalWavenumber,
        config.upperMaximumMeridionalWavenumber, upperEquilibriumPressure);
    if (config.enabled.stationary && config.enabled.surface)
    {
        result.surfaceStationarySolver = solveSteadyStationaryWavePressure(
            longitudeCells,
            latitudeCells,
            surfaceEquilibrium,
            surfaceDrag,
            config.surfaceEquivalentPressureDepthHpa,
            config.surfaceDampingTimeSeconds,
            config.airDensityKgM3,
            config.planetRadiusMetres,
            config.rotationRatePerSecond,
            config.rotationDirection,
            true,
            config.maximumIterations,
            config.relativeTolerance,
            config.solverRestartLength);
    }
    else
    {
        result.surfaceStationarySolver.equilibriumPressureAnomalyHpa.assign(cellCount, 0.0f);
        result.surfaceStationarySolver.pressureAnomalyHpa.assign(cellCount, 0.0f);
        result.surfaceStationarySolver.converged = true;
    }

    if (config.enabled.stationary && config.enabled.upper)
    {
        result.upperStationarySolver = solveSteadyStationaryWavePressure(
            longitudeCells,
            latitudeCells,
            upperEquilibriumPressure,
            upperDrag,
            config.upperEquivalentPressureDepthHpa,
            config.upperDampingTimeSeconds,
            config.airDensityKgM3,
            config.planetRadiusMetres,
            config.rotationRatePerSecond,
            config.rotationDirection,
            true,
            config.maximumIterations,
            config.relativeTolerance,
            config.solverRestartLength);
    }
    else
    {
        result.upperStationarySolver.equilibriumPressureAnomalyHpa.assign(cellCount, 0.0f);
        result.upperStationarySolver.pressureAnomalyHpa.assign(cellCount, 0.0f);
        result.upperStationarySolver.converged = true;
    }

    for (int y = 0; y < latitudeCells; y++)
    {
        for (int x = 0; x < longitudeCells; x++)
        {
            const std::size_t cell = static_cast<std::size_t>(y) * longitudeCells + x;
            if (config.enabled.surface)
            {
                result.surfacePressureAnomalyHpa[cell] =
                    (config.enabled.zonal ? zonalSurfacePressureHpa[y] : 0.0f) +
                    (result.surfaceStationarySolver.converged
                        ? result.surfaceStationarySolver.pressureAnomalyHpa[cell] : 0.0f);
            }
            if (config.enabled.upper)
            {
                result.upperHeightAnomalyMetres[cell] =
                    (config.enabled.zonal && config.zonalUpperHeightMetres.size() ==
                        static_cast<std::size_t>(latitudeCells) ? config.zonalUpperHeightMetres[y] : 0.0f) +
                    (result.upperStationarySolver.converged
                        ? result.upperStationarySolver.pressureAnomalyHpa[cell] /
                            std::max(1.0e-6f, upperPressurePerMetre) : 0.0f);
            }
        }
    }

    diagnoseModeWinds(longitudeCells, latitudeCells, config, result);
    return result;
}

void diagnoseModeWinds(int longitudeCells, int latitudeCells,
    const ModeSeparatedCirculationConfig& config, ModeSeparatedCirculation& result)
{
    const std::size_t cellCount = static_cast<std::size_t>(std::max(0, longitudeCells)) * std::max(0, latitudeCells);
    if (longitudeCells < 3 || latitudeCells < 3 || result.surfacePressureAnomalyHpa.size() != cellCount ||
        result.upperHeightAnomalyMetres.size() != cellCount)
        return;
    result.surfaceEastWindMps.assign(cellCount, 0.0f);
    result.surfaceSouthWindMps.assign(cellCount, 0.0f);
    result.upperEastWindMps.assign(cellCount, 0.0f);
    result.upperSouthWindMps.assign(cellCount, 0.0f);
    const double latitudeSpacing = pi * static_cast<double>(config.planetRadiusMetres) /
        static_cast<double>(latitudeCells);
    for (int y = 0; y < latitudeCells; y++)
    {
        const int north = std::max(0, y - 1);
        const int south = std::min(latitudeCells - 1, y + 1);
        const float latitudeDegrees = 90.0f - 180.0f *
            (static_cast<float>(y) + 0.5f) / static_cast<float>(latitudeCells);
        const double zonalSpacing = 2.0 * pi * config.planetRadiusMetres *
            std::max(0.02, std::abs(std::cos(latitudeDegrees * pi / 180.0))) /
            static_cast<double>(longitudeCells);
        for (int x = 0; x < longitudeCells; x++)
        {
            const int west = (x + longitudeCells - 1) % longitudeCells;
            const int east = (x + 1) % longitudeCells;
            const std::size_t cell = static_cast<std::size_t>(y) * longitudeCells + x;
            if (config.enabled.surface)
            {
                const float pressureGradientEast = static_cast<float>(
                    (result.surfacePressureAnomalyHpa[
                        static_cast<std::size_t>(y) * longitudeCells + east] -
                     result.surfacePressureAnomalyHpa[
                        static_cast<std::size_t>(y) * longitudeCells + west]) * 100.0 /
                    (2.0 * zonalSpacing));
                const float pressureGradientNorth = static_cast<float>(
                    (result.surfacePressureAnomalyHpa[
                        static_cast<std::size_t>(north) * longitudeCells + x] -
                     result.surfacePressureAnomalyHpa[
                        static_cast<std::size_t>(south) * longitudeCells + x]) * 100.0 /
                    ((south - north) * latitudeSpacing));
                const HorizontalWind wind = steadyQuadraticDragCoriolisWind(
                    -pressureGradientEast / config.airDensityKgM3,
                    -pressureGradientNorth / config.airDensityKgM3,
                    latitudeDegrees,
                    config.surfaceDragCoefficients.size() == cellCount
                        ? config.surfaceDragCoefficients[cell] : config.surfaceDragCoefficient,
                    config.surfaceBoundaryLayerDepthMetres,
                    config.rotationRatePerSecond,
                    config.rotationDirection);
                result.surfaceEastWindMps[cell] = wind.eastMetresPerSecond;
                result.surfaceSouthWindMps[cell] = wind.southMetresPerSecond;
            }
            if (config.enabled.upper)
            {
                const float heightGradientEast = static_cast<float>(
                    (result.upperHeightAnomalyMetres[
                        static_cast<std::size_t>(y) * longitudeCells + east] -
                     result.upperHeightAnomalyMetres[
                        static_cast<std::size_t>(y) * longitudeCells + west]) /
                    (2.0 * zonalSpacing));
                const float heightGradientNorth = static_cast<float>(
                    (result.upperHeightAnomalyMetres[
                        static_cast<std::size_t>(north) * longitudeCells + x] -
                     result.upperHeightAnomalyMetres[
                        static_cast<std::size_t>(south) * longitudeCells + x]) /
                    ((south - north) * latitudeSpacing));
                const HorizontalWind wind = steadyRayleighCoriolisWind(
                    -config.gravityMetresPerSecondSquared * heightGradientEast,
                    -config.gravityMetresPerSecondSquared * heightGradientNorth,
                    latitudeDegrees,
                    config.upperDragTimeSeconds,
                    config.rotationRatePerSecond,
                    config.rotationDirection);
                result.upperEastWindMps[cell] = wind.eastMetresPerSecond;
                result.upperSouthWindMps[cell] = wind.southMetresPerSecond;
            }
        }
    }
    // Weak, equal-and-opposite exchange of stationary momentum only. Removing
    // each transfer's row mean leaves the independent zonal closures intact.
    result.maximumMomentumExchangeResidual = 0.0;
    if (config.enabled.surface && config.enabled.upper && config.enabled.stationary)
    {
        const double lowerMass = config.surfaceEquivalentPressureDepthHpa;
        const double upperMass = config.upperEquivalentPressureDepthHpa;
        const double coupling = std::clamp(config.interlayerMomentumCoupling, 0.0f, 1.0f);
        for (int y = 0; y < latitudeCells; ++y)
            for (int component = 0; component < 2; ++component)
            {
                auto& lower = component == 0 ? result.surfaceEastWindMps : result.surfaceSouthWindMps;
                auto& upper = component == 0 ? result.upperEastWindMps : result.upperSouthWindMps;
                double rowShear = 0.0;
                for (int x = 0; x < longitudeCells; ++x)
                    rowShear += (upper[y * longitudeCells + x] - lower[y * longitudeCells + x]) / longitudeCells;
                for (int x = 0; x < longitudeCells; ++x)
                {
                    const int cell = y * longitudeCells + x;
                    const double shear = coupling * (upper[cell] - lower[cell] - rowShear);
                    const double lowerDelta = shear * upperMass / (lowerMass + upperMass);
                    const double upperDelta = -shear * lowerMass / (lowerMass + upperMass);
                    lower[cell] += static_cast<float>(lowerDelta);
                    upper[cell] += static_cast<float>(upperDelta);
                    result.maximumMomentumExchangeResidual = std::max(result.maximumMomentumExchangeResidual,
                        std::abs(lowerMass * lowerDelta + upperMass * upperDelta));
                }
            }
    }
    const auto grid = climategrid::makeSphericalGrid(longitudeCells, latitudeCells, config.planetRadiusMetres);
    result.ascentHpaPerDay.assign(cellCount, 0.0f);
    result.areaWeightedKineticEnergyJm2 = 0.0;
    result.areaWeightedDragDissipationWm2 = 0.0;
    result.areaWeightedMassAnomalyKgM2 = 0.0;
    result.maximumStationaryRowMeanHpa = 0.0f;
    double totalArea = 0.0;
    for (int y = 0; y < latitudeCells; ++y)
    {
        double rowMean = 0.0;
        for (int x = 0; x < longitudeCells; ++x)
        {
            const auto cell = grid.index(x, y);
            const double area = grid.cellAreasSquareMetres[y];
            totalArea += area;
            rowMean += result.surfaceStationarySolver.converged &&
                result.surfaceStationarySolver.pressureAnomalyHpa.size() == cellCount
                ? result.surfaceStationarySolver.pressureAnomalyHpa[cell] : 0.0;
            const double surfaceSpeed2 = std::pow(result.surfaceEastWindMps[cell], 2) +
                std::pow(result.surfaceSouthWindMps[cell], 2);
            const double upperSpeed2 = std::pow(result.upperEastWindMps[cell], 2) +
                std::pow(result.upperSouthWindMps[cell], 2);
            const double surfaceMass = config.surfaceEquivalentPressureDepthHpa * 100.0 / config.gravityMetresPerSecondSquared;
            const double upperMass = config.upperEquivalentPressureDepthHpa * 100.0 / config.gravityMetresPerSecondSquared;
            result.areaWeightedKineticEnergyJm2 += area * 0.5 * (surfaceMass * surfaceSpeed2 + upperMass * upperSpeed2);
            result.areaWeightedDragDissipationWm2 += area *
                (surfaceMass * (config.surfaceDragCoefficients.size() == cellCount
                    ? config.surfaceDragCoefficients[cell] : config.surfaceDragCoefficient) / config.surfaceBoundaryLayerDepthMetres *
                    surfaceSpeed2 * std::sqrt(surfaceSpeed2) + upperMass * upperSpeed2 / config.upperDragTimeSeconds);
            result.areaWeightedMassAnomalyKgM2 += area * result.surfacePressureAnomalyHpa[cell] *
                100.0 / config.gravityMetresPerSecondSquared;
            const auto divergence = [&](const std::vector<float>& u, const std::vector<float>& v)
            {
                const double east = 0.5 * (u[cell] + u[grid.index(x + 1, y)]) * grid.zonalFaceLengthsMetres[y];
                const double west = 0.5 * (u[cell] + u[grid.index(x - 1, y)]) * grid.zonalFaceLengthsMetres[y];
                const double southFlux = y + 1 < latitudeCells
                    ? 0.5 * (v[cell] + v[grid.index(x, y + 1)]) * grid.southFaceLengthsMetres[y] : 0.0;
                const double northFlux = y > 0
                    ? 0.5 * (v[cell] + v[grid.index(x, y - 1)]) * grid.northFaceLengthsMetres[y] : 0.0;
                return (east - west + southFlux - northFlux) / area;
            };
            // Positive is ascent. Equal-and-opposite interface mass exchange;
            // closed face fluxes make its global area integral zero.
            result.ascentHpaPerDay[cell] = static_cast<float>(43200.0 *
                (config.upperEquivalentPressureDepthHpa * divergence(result.upperEastWindMps, result.upperSouthWindMps) -
                 config.surfaceEquivalentPressureDepthHpa * divergence(result.surfaceEastWindMps, result.surfaceSouthWindMps)));
        }
        result.maximumStationaryRowMeanHpa = std::max(result.maximumStationaryRowMeanHpa,
            static_cast<float>(std::abs(rowMean / longitudeCells)));
    }
    if (totalArea > 0.0)
    {
        result.areaWeightedKineticEnergyJm2 /= totalArea;
        result.areaWeightedDragDissipationWm2 /= totalArea;
        result.areaWeightedMassAnomalyKgM2 /= totalArea;
    }
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
