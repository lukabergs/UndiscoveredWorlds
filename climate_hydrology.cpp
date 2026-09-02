#include "climate_hydrology.hpp"

#include "climate_grid.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace climatehydrology
{
namespace
{
constexpr std::array<int, monthCount> daysInMonth = {
    31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31
};

float relaxedExcess(
    float availableColumnWaterMm,
    float criticalColumnWaterMm,
    float timeStepSeconds,
    float conversionTimeSeconds)
{
    if (timeStepSeconds <= 0.0f || conversionTimeSeconds <= 0.0f)
        return 0.0f;

    const float excess = std::max(
        0.0f,
        availableColumnWaterMm - std::max(0.0f, criticalColumnWaterMm));
    const float convertedFraction = -std::expm1(-timeStepSeconds / conversionTimeSeconds);
    return std::clamp(
        excess * convertedFraction,
        0.0f,
        std::max(0.0f, availableColumnWaterMm));
}

float smoothStep(float lower, float upper, float value)
{
    if (upper <= lower)
        return value >= upper ? 1.0f : 0.0f;

    const float phase = std::clamp((value - lower) / (upper - lower), 0.0f, 1.0f);
    return phase * phase * (3.0f - 2.0f * phase);
}
}

float PrecipitationPartition::totalMm() const
{
    return stratiformMm + orographicMm + convectiveMm;
}

float FallingPrecipitation::surfaceTotalMm() const
{
    return rainMm + snowMm;
}

CalendarMonth calendarMonth(int month)
{
    const int wrappedMonth = ((month % monthCount) + monthCount) % monthCount;
    const int firstSeason = wrappedMonth / 3;
    return {
        wrappedMonth,
        firstSeason,
        (firstSeason + 1) % seasonCount,
        static_cast<float>(wrappedMonth % 3) / 3.0f,
        daysInMonth[wrappedMonth]
    };
}

ClimateGridDimensions climateGridDimensions(
    int outputColumns,
    int outputRows,
    int targetHorizontalCells)
{
    const int requested = std::min(outputColumns, targetHorizontalCells);
    const int columns = std::max(2, requested - requested % 2);
    // Climate state always uses the finite-volume 2:1 spherical contract.
    // Requested raster dimensions affect only the explicit input/output remap.
    (void)outputRows;
    return { columns, std::max(1, columns / 2) };
}

float climateCellLatitudeDegrees(int row, int rowCount)
{
    constexpr double radiansToDegrees = 57.2957795130823208768;
    return static_cast<float>(climategrid::latitudeCentreRadians(
        row, std::max(1, rowCount), climategrid::LatitudeLayout::cellCentred) * radiansToDegrees);
}

float climateCellAreaWeight(int row, int rowCount)
{
    constexpr double pi = 3.14159265358979323846;
    const int rows = std::max(1, rowCount);
    return static_cast<float>(climategrid::latitudeBandMeasure(
        row, rows, climategrid::LatitudeLayout::cellCentred) * rows / pi);
}

float polarTaperFactor(
    float latitudeDegrees,
    float taperStartDegrees,
    float taperEndDegrees)
{
    return 1.0f - smoothStep(
        std::abs(taperStartDegrees),
        std::abs(taperEndDegrees),
        std::abs(latitudeDegrees));
}

int adjacentMeridionalTransportTargetRow(
    int sourceRow,
    float displacementRows,
    int maximumRow)
{
    const int boundedMaximum = std::max(0, maximumRow);
    const int boundedSource = std::clamp(sourceRow, 0, boundedMaximum);
    const int direction = (displacementRows > 0.0f) - (displacementRows < 0.0f);
    return std::clamp(boundedSource + direction, 0, boundedMaximum);
}

SphericalTracerTransportDiagnostics advectSphericalTracer(
    int columns,
    int rows,
    const std::vector<float>& source,
    const std::vector<float>& zonalWindMps,
    const std::vector<float>& meridionalWindMps,
    float timeStepSeconds,
    float planetRadiusMetres,
    float maximumMeridionalCourantPerSubstep,
    float maximumDisplacementCells,
    std::vector<float>& destination)
{
    (void)maximumDisplacementCells;
    MpdataOptions options;
    options.maximumCourantPerSubstep = maximumMeridionalCourantPerSubstep;
    return advectSphericalTracerMpdata(
        columns,
        rows,
        source,
        zonalWindMps,
        meridionalWindMps,
        timeStepSeconds,
        planetRadiusMetres,
        options,
        destination);
}

SphericalTracerTransportDiagnostics advectSphericalTracerMpdata(
    int columns,
    int rows,
    const std::vector<float>& source,
    const std::vector<float>& zonalWindMps,
    const std::vector<float>& meridionalWindMps,
    float timeStepSeconds,
    float planetRadiusMetres,
    const MpdataOptions& options,
    std::vector<float>& destination)
{
    struct FaceFlux
    {
        std::size_t first = 0;
        std::size_t second = 0;
        double sweptArea = 0.0;
        bool meridional = false;
        double initialSweptArea = 0.0;
        double finalSweptArea = 0.0;
    };

    SphericalTracerTransportDiagnostics diagnostics;
    const std::size_t cellCount = static_cast<std::size_t>(std::max(0, columns)) *
        static_cast<std::size_t>(std::max(0, rows));
    if (columns <= 0 || rows <= 0 || source.size() != cellCount ||
        zonalWindMps.size() != cellCount || meridionalWindMps.size() != cellCount ||
        timeStepSeconds <= 0.0f || planetRadiusMetres <= 0.0f ||
        (!options.endZonalWindMps.empty() && options.endZonalWindMps.size() != cellCount) ||
        (!options.endMeridionalWindMps.empty() && options.endMeridionalWindMps.size() != cellCount))
    {
        destination = source;
        return diagnostics;
    }

    const climategrid::SphericalGrid grid = climategrid::makeSphericalGrid(
        columns, rows, planetRadiusMetres);
    diagnostics.eastIntegratedFlux.assign(cellCount, 0.0);
    diagnostics.southIntegratedFlux.assign(cellCount, 0.0);
    const auto massIntegral = [&](const std::vector<float>& field)
    {
        double mass = 0.0;
        for (int y = 0; y < rows; y++)
        {
            const double area = grid.cellAreasSquareMetres[y];
            for (int x = 0; x < columns; x++)
                mass += area * static_cast<double>(field[grid.index(x, y)]);
        }
        return mass;
    };

    std::vector<FaceFlux> fullStepFaces;
    fullStepFaces.reserve(static_cast<std::size_t>(columns) * (2 * rows - 1));
    std::vector<double> outgoingArea(cellCount, 0.0);
    for (int y = 0; y < rows; y++)
    {
        const double cellArea = grid.cellAreasSquareMetres[y];
        for (int x = 0; x < columns; x++)
        {
            const std::size_t west = grid.index(x, y);
            const std::size_t east = grid.index(x + 1, y);
            const double eastWind = 0.5 * (
                static_cast<double>(zonalWindMps[west]) + zonalWindMps[east]);
            const double sweptArea = eastWind * grid.zonalFaceLengthsMetres[y] *
                static_cast<double>(timeStepSeconds);
            fullStepFaces.push_back({ west, east, sweptArea });
            diagnostics.maximumZonalCourant = std::max(
                diagnostics.maximumZonalCourant,
                static_cast<float>(std::abs(sweptArea) / cellArea));
            outgoingArea[sweptArea >= 0.0 ? west : east] += std::abs(sweptArea);

            if (y >= rows - 1)
                continue;
            const std::size_t north = west;
            const std::size_t south = grid.index(x, y + 1);
            const double southWind = 0.5 * (
                static_cast<double>(meridionalWindMps[north]) + meridionalWindMps[south]);
            const double meridionalSweptArea = southWind *
                grid.southFaceLengthsMetres[y] * static_cast<double>(timeStepSeconds);
            fullStepFaces.push_back({ north, south, meridionalSweptArea, true });
            const double smallerArea = std::min(
                grid.cellAreasSquareMetres[y], grid.cellAreasSquareMetres[y + 1]);
            diagnostics.maximumMeridionalCourant = std::max(
                diagnostics.maximumMeridionalCourant,
                static_cast<float>(std::abs(meridionalSweptArea) / smallerArea));
            outgoingArea[meridionalSweptArea >= 0.0 ? north : south] +=
                std::abs(meridionalSweptArea);
        }
    }

    // Bound the complete time-varying route with the outgoing maxima on every
    // face. This remains safe when a face reverses direction during the step.
    std::fill(outgoingArea.begin(), outgoingArea.end(), 0.0);
    for (auto& face : fullStepFaces)
    {
        const int row = static_cast<int>(face.first / columns);
        const auto& endWind = face.meridional ? options.endMeridionalWindMps : options.endZonalWindMps;
        face.initialSweptArea = face.sweptArea;
        face.finalSweptArea = endWind.empty() ? face.sweptArea :
            0.5 * (endWind[face.first] + endWind[face.second]) * timeStepSeconds *
                (face.meridional ? grid.southFaceLengthsMetres[row] : grid.zonalFaceLengthsMetres[row]);
        outgoingArea[face.first] += std::max({0.0, face.initialSweptArea, face.finalSweptArea});
        outgoingArea[face.second] += std::max({0.0, -face.initialSweptArea, -face.finalSweptArea});
    }
    for (int y = 0; y < rows; y++)
    {
        const double inverseArea = 1.0 / grid.cellAreasSquareMetres[y];
        for (int x = 0; x < columns; x++)
        {
            diagnostics.maximumMultidimensionalCourant = std::max(
                diagnostics.maximumMultidimensionalCourant,
                static_cast<float>(outgoingArea[grid.index(x, y)] * inverseArea));
        }
    }
    const float courantLimit = std::clamp(
        options.maximumCourantPerSubstep, 0.05f, 0.95f);
    diagnostics.substeps = std::max(
        1,
        static_cast<int>(std::ceil(
            diagnostics.maximumMultidimensionalCourant / courantLimit)));
    diagnostics.correctivePasses = std::clamp(options.correctivePasses, 0, 3);
    diagnostics.initialAreaWeightedMass = massIntegral(source);

    std::vector<float> current = source;
    std::vector<double> mass(cellCount, 0.0);
    std::vector<float> donor(cellCount, 0.0f);
    std::vector<double> positiveFlux(cellCount, 0.0);
    std::vector<double> negativeFlux(cellCount, 0.0);
    std::vector<double> positiveRatio(cellCount, 1.0);
    std::vector<double> negativeRatio(cellCount, 1.0);
    std::vector<double> correction(fullStepFaces.size(), 0.0);
    std::vector<double> advector(fullStepFaces.size(), 0.0);
    std::vector<double> nextAdvector(fullStepFaces.size(), 0.0);
    std::vector<double> eastAdvector(cellCount, 0.0);
    std::vector<double> southAdvector(cellCount, 0.0);
    std::vector<double> divergence(cellCount, 0.0);
    std::vector<float> lowerBound(cellCount, 0.0f);
    std::vector<float> upperBound(cellCount, 0.0f);

    const auto localBounds = [&](const std::vector<float>& field)
    {
        for (int y = 0; y < rows; y++)
        {
            const int north = std::max(0, y - 1);
            const int south = std::min(rows - 1, y + 1);
            for (int x = 0; x < columns; x++)
            {
                const std::size_t cell = grid.index(x, y);
                float minimum = field[cell];
                float maximum = field[cell];
                for (const std::size_t neighbour : {
                         grid.index(x - 1, y), grid.index(x + 1, y),
                         grid.index(x, north), grid.index(x, south) })
                {
                    minimum = std::min(minimum, field[neighbour]);
                    maximum = std::max(maximum, field[neighbour]);
                }
                lowerBound[cell] = options.monotone ? std::max(0.0f, minimum) : 0.0f;
                upperBound[cell] = options.monotone
                    ? maximum
                    : std::numeric_limits<float>::max();
            }
        }
    };

    for (int step = 0; step < diagnostics.substeps; step++)
    {
        const double fraction = (step + 0.5) / diagnostics.substeps;
        for (auto& face : fullStepFaces)
            face.sweptArea = face.initialSweptArea + fraction * (face.finalSweptArea - face.initialSweptArea);
        for (int y = 0; y < rows; y++)
        {
            const double area = grid.cellAreasSquareMetres[y];
            for (int x = 0; x < columns; x++)
            {
                const std::size_t cell = grid.index(x, y);
                mass[cell] = area * static_cast<double>(current[cell]);
            }
        }

        for (const FaceFlux& fullFace : fullStepFaces)
        {
            const double sweptArea = fullFace.sweptArea /
                static_cast<double>(diagnostics.substeps);
            const float upwind = sweptArea >= 0.0
                ? current[fullFace.first]
                : current[fullFace.second];
            const double tracerFlux = sweptArea * static_cast<double>(upwind);
            mass[fullFace.first] -= tracerFlux;
            mass[fullFace.second] += tracerFlux;
            (fullFace.meridional ? diagnostics.southIntegratedFlux : diagnostics.eastIntegratedFlux)[fullFace.first] += tracerFlux;
        }
        for (int y = 0; y < rows; y++)
        {
            const double inverseArea = 1.0 / grid.cellAreasSquareMetres[y];
            for (int x = 0; x < columns; x++)
            {
                const std::size_t cell = grid.index(x, y);
                donor[cell] = static_cast<float>(mass[cell] * inverseArea);
            }
        }

        localBounds(current);
        current.swap(donor);
        // Allow extrema introduced by physical compression/expansion in the
        // donor pass. The correction must not undo a real divergent tendency.
        for (std::size_t cell = 0; cell < cellCount; cell++)
        {
            lowerBound[cell] = std::min(lowerBound[cell], current[cell]);
            upperBound[cell] = std::max(upperBound[cell], current[cell]);
        }
        for (std::size_t face = 0; face < fullStepFaces.size(); face++)
            advector[face] = fullStepFaces[face].sweptArea / diagnostics.substeps;
        for (int pass = 0; pass < diagnostics.correctivePasses; pass++)
        {
            std::fill(positiveFlux.begin(), positiveFlux.end(), 0.0);
            std::fill(negativeFlux.begin(), negativeFlux.end(), 0.0);
            std::fill(divergence.begin(), divergence.end(), 0.0);
            for (std::size_t face = 0; face < fullStepFaces.size(); face++)
            {
                const auto& edge = fullStepFaces[face];
                (edge.meridional ? southAdvector : eastAdvector)[edge.first] = advector[face];
                divergence[edge.first] += advector[face];
                divergence[edge.second] -= advector[face];
            }
            for (std::size_t cell = 0; cell < cellCount; cell++)
                divergence[cell] /= grid.cellAreasSquareMetres[cell / columns];
            for (std::size_t faceIndex = 0; faceIndex < fullStepFaces.size(); faceIndex++)
            {
                const FaceFlux& fullFace = fullStepFaces[faceIndex];
                const double sweptArea = advector[faceIndex];
                const double firstArea = grid.cellAreasSquareMetres[
                    static_cast<int>(fullFace.first / static_cast<std::size_t>(columns))];
                const double secondArea = grid.cellAreasSquareMetres[
                    static_cast<int>(fullFace.second / static_cast<std::size_t>(columns))];
                const double faceArea = 0.5 * (firstArea + secondArea);
                const int x = static_cast<int>(fullFace.first % columns);
                const int y = static_cast<int>(fullFace.first / columns);
                const int north = std::max(0, y - 1);
                const int south = std::min(rows - 1, y + 1);
                double transverseAdvector = 0.0;
                double transverseNumerator = 0.0;
                double transverseDenominator = 0.0;
                if (fullFace.meridional)
                {
                    transverseAdvector = 0.25 * (
                        eastAdvector[grid.index(x, y)] + eastAdvector[grid.index(x - 1, y)] +
                        eastAdvector[grid.index(x, south)] + eastAdvector[grid.index(x - 1, south)]);
                    transverseNumerator = current[grid.index(x + 1, y)] +
                        current[grid.index(x + 1, south)] - current[grid.index(x - 1, y)] -
                        current[grid.index(x - 1, south)];
                    transverseDenominator = current[grid.index(x + 1, y)] +
                        current[grid.index(x + 1, south)] + current[grid.index(x - 1, y)] +
                        current[grid.index(x - 1, south)];
                }
                else
                {
                    transverseAdvector = 0.25 * (
                        southAdvector[grid.index(x, y)] + southAdvector[grid.index(x + 1, y)] +
                        (y > 0 ? southAdvector[grid.index(x, north)] +
                            southAdvector[grid.index(x + 1, north)] : 0.0));
                    transverseNumerator = current[grid.index(x, south)] +
                        current[grid.index(x + 1, south)] - current[grid.index(x, north)] -
                        current[grid.index(x + 1, north)];
                    transverseDenominator = current[grid.index(x, south)] +
                        current[grid.index(x + 1, south)] + current[grid.index(x, north)] +
                        current[grid.index(x + 1, north)];
                }
                const double sum = current[fullFace.first] + current[fullFace.second];
                const double gradientRatio = sum > 1.0e-30
                    ? (current[fullFace.second] - current[fullFace.first]) / sum : 0.0;
                const double crossRatio = transverseDenominator > 1.0e-30
                    ? 0.5 * transverseNumerator / transverseDenominator : 0.0;
                // Generalized-coordinate MPDATA: normal, cross-dimensional,
                // and divergent-flow terms. Every corrective pass is another
                // donor-cell pass driven by the preceding pseudo-advector.
                // Jaruga et al., doi:10.5194/gmd-8-1005-2015, section 3.1.
                const double pseudoAdvector =
                    (std::abs(sweptArea) - sweptArea * sweptArea / faceArea) * gradientRatio -
                    sweptArea * transverseAdvector / faceArea * crossRatio -
                    0.25 * sweptArea * (divergence[fullFace.first] + divergence[fullFace.second]);
                nextAdvector[faceIndex] = pseudoAdvector;
                const double antidiffusiveFlux = pseudoAdvector *
                    (pseudoAdvector >= 0.0 ? current[fullFace.first] : current[fullFace.second]);
                correction[faceIndex] = antidiffusiveFlux;
                if (antidiffusiveFlux >= 0.0)
                {
                    negativeFlux[fullFace.first] += antidiffusiveFlux;
                    positiveFlux[fullFace.second] += antidiffusiveFlux;
                }
                else
                {
                    positiveFlux[fullFace.first] -= antidiffusiveFlux;
                    negativeFlux[fullFace.second] -= antidiffusiveFlux;
                }
            }

            for (int y = 0; y < rows; y++)
            {
                const double area = grid.cellAreasSquareMetres[y];
                for (int x = 0; x < columns; x++)
                {
                    const std::size_t cell = grid.index(x, y);
                    const double allowedIncrease = area * std::max(
                        0.0,
                        static_cast<double>(upperBound[cell] - current[cell]));
                    const double allowedDecrease = area * std::max(
                        0.0,
                        static_cast<double>(current[cell] - lowerBound[cell]));
                    positiveRatio[cell] = positiveFlux[cell] > 0.0
                        ? std::min(1.0, allowedIncrease / positiveFlux[cell])
                        : 1.0;
                    negativeRatio[cell] = negativeFlux[cell] > 0.0
                        ? std::min(1.0, allowedDecrease / negativeFlux[cell])
                        : 1.0;
                    mass[cell] = area * static_cast<double>(current[cell]);
                }
            }

            for (std::size_t faceIndex = 0; faceIndex < fullStepFaces.size(); faceIndex++)
            {
                const FaceFlux& face = fullStepFaces[faceIndex];
                double limitedFlux = correction[faceIndex];
                if (limitedFlux >= 0.0)
                {
                    limitedFlux *= std::min(
                        negativeRatio[face.first], positiveRatio[face.second]);
                }
                else
                {
                    limitedFlux *= std::min(
                        positiveRatio[face.first], negativeRatio[face.second]);
                }
                mass[face.first] -= limitedFlux;
                mass[face.second] += limitedFlux;
                (face.meridional ? diagnostics.southIntegratedFlux : diagnostics.eastIntegratedFlux)[face.first] += limitedFlux;
                nextAdvector[faceIndex] *= correction[faceIndex] != 0.0
                    ? limitedFlux / correction[faceIndex] : 0.0;
            }
            for (int y = 0; y < rows; y++)
            {
                const double inverseArea = 1.0 / grid.cellAreasSquareMetres[y];
                for (int x = 0; x < columns; x++)
                {
                    const std::size_t cell = grid.index(x, y);
                    current[cell] = static_cast<float>(mass[cell] * inverseArea);
                }
            }
            advector.swap(nextAdvector);
        }
    }

    destination = std::move(current);
    diagnostics.finalAreaWeightedMass = massIntegral(destination);
    diagnostics.minimumMixingRatio = *std::min_element(
        destination.begin(), destination.end());
    return diagnostics;
}

MeanMoistureTransport meanMoistureTransport(const SeasonalProcessFields& fields, double radiusMetres)
{
    MeanMoistureTransport result;
    const int cells = fields.columns * fields.rows;
    if (cells <= 0 || fields.durationSeconds <= 0.0) return result;
    const auto grid = climategrid::makeSphericalGrid(fields.columns, fields.rows, radiusMetres);
    result.convergenceMmPerDay.assign(cells, 0.0f);
    for (int layer = 0; layer < 2; ++layer)
    {
        if (fields.eastIntegratedFlux[layer].size() != static_cast<std::size_t>(cells) ||
            fields.southIntegratedFlux[layer].size() != static_cast<std::size_t>(cells)) return {};
        result.eastKgPerMetreSecond[layer].resize(cells);
        result.southKgPerMetreSecond[layer].resize(cells);
        for (int y = 0; y < fields.rows; ++y)
            for (int x = 0; x < fields.columns; ++x)
            {
                const auto cell = grid.index(x, y), west = grid.index(x - 1, y), north = grid.index(x, y - 1);
                const double e = fields.eastIntegratedFlux[layer][cell], w = fields.eastIntegratedFlux[layer][west];
                const double s = fields.southIntegratedFlux[layer][cell], n = y > 0 ? fields.southIntegratedFlux[layer][north] : 0.0;
                result.eastKgPerMetreSecond[layer][cell] = static_cast<float>((e + w) /
                    (2.0 * fields.durationSeconds * grid.zonalFaceLengthsMetres[y]));
                result.southKgPerMetreSecond[layer][cell] = static_cast<float>((s + n) /
                    (fields.durationSeconds * std::max(1.0, grid.northFaceLengthsMetres[y] + grid.southFaceLengthsMetres[y])));
                result.convergenceMmPerDay[cell] += static_cast<float>((w - e + n - s) * 86400.0 /
                    (fields.durationSeconds * grid.cellAreasSquareMetres[y]));
            }
    }
    return result;
}

WeatherPhase deterministicWeatherPhase(
    int phase,
    int phaseCount,
    float windRotationDegrees,
    float daytimeLandTemperatureAnomalyC,
    float nighttimeLandTemperatureAnomalyC,
    float daytimeSeaTemperatureAnomalyC,
    float nighttimeSeaTemperatureAnomalyC)
{
    constexpr float pi = 3.14159265358979323846f;
    const int count = std::max(1, phaseCount);
    const int wrapped = ((phase % count) + count) % count;
    const float centred = count > 1
        ? 2.0f * static_cast<float>(wrapped) / static_cast<float>(count - 1) - 1.0f
        : 0.0f;
    const float coastalDirection = count > 1
        ? (wrapped == 0 ? 1.0f : (wrapped == count - 1 ? -1.0f : 0.0f))
        : 0.0f;

    return {
        centred * windRotationDegrees * pi / 180.0f,
        2.0f * pi * static_cast<float>(wrapped) / static_cast<float>(count),
        coastalDirection,
        coastalDirection > 0.0f
            ? daytimeLandTemperatureAnomalyC
            : (coastalDirection < 0.0f ? nighttimeLandTemperatureAnomalyC : 0.0f),
        coastalDirection > 0.0f
            ? daytimeSeaTemperatureAnomalyC
            : (coastalDirection < 0.0f ? nighttimeSeaTemperatureAnomalyC : 0.0f)
    };
}

float interpolateSeasonal(float first, float second, float interpolation)
{
    const float phase = std::clamp(interpolation, 0.0f, 1.0f);
    return first + (second - first) * phase;
}

float kuoPrecipitationEfficiency(
    float relativeHumidity,
    float criticalRelativeHumidity,
    float humidityExponent)
{
    const float critical = std::clamp(criticalRelativeHumidity, 0.0f, 0.999f);
    const float humidity = std::clamp(relativeHumidity, 0.0f, 1.0f);

    if (humidity <= critical)
        return 0.0f;

    const float dryFraction = std::clamp(
        (1.0f - humidity) / (1.0f - critical),
        0.0f,
        1.0f);
    return 1.0f - std::pow(dryFraction, std::max(0.0f, humidityExponent));
}

float convectiveBuoyancyEfficiency(
    float parcelBuoyancyC,
    float activationBuoyancyC,
    float fullStrengthBuoyancyC)
{
    return smoothStep(activationBuoyancyC, fullStrengthBuoyancyC, parcelBuoyancyC);
}

float shallowConvectionExchangeFraction(
    float boundaryRelativeHumidity,
    float freeTroposphereRelativeHumidity,
    float parcelBuoyancyC,
    float verticalWindShearMps,
    float timeStepSeconds,
    float mixingTimeDays,
    float humidityOnset,
    float fullHumidity,
    float fullShearMps,
    float maximumExchangeFraction)
{
    constexpr float secondsPerDay = 86400.0f;
    const float humidity = std::max(boundaryRelativeHumidity, freeTroposphereRelativeHumidity);
    const float humidityFactor = smoothStep(humidityOnset, fullHumidity, humidity);
    const float instability = smoothStep(0.0f, 4.0f, parcelBuoyancyC);
    const float shearFactor = std::clamp(
        verticalWindShearMps / std::max(0.1f, fullShearMps),
        0.0f,
        1.0f);
    const float timeFraction = -std::expm1(
        -std::max(0.0f, timeStepSeconds) /
        (std::max(0.01f, mixingTimeDays) * secondsPerDay));
    return std::clamp(
        timeFraction * humidityFactor * instability * (0.5f + 0.5f * shearFactor),
        0.0f,
        std::max(0.0f, maximumExchangeFraction));
}

float dryConvectionExchangeFraction(
    float parcelBuoyancyC,
    float timeStepSeconds,
    float mixingTimeDays,
    float activationBuoyancyC,
    float fullStrengthBuoyancyC,
    float maximumExchangeFraction)
{
    constexpr float secondsPerDay = 86400.0f;
    const float timeFraction = -std::expm1(
        -std::max(0.0f, timeStepSeconds) /
        (std::max(0.01f, mixingTimeDays) * secondsPerDay));
    return std::clamp(
        timeFraction * convectiveBuoyancyEfficiency(
            parcelBuoyancyC,
            activationBuoyancyC,
            fullStrengthBuoyancyC),
        0.0f,
        std::max(0.0f, maximumExchangeFraction));
}

float soilMoistureStress(
    float soilMoistureMm,
    float soilMoistureCapacityMm,
    float criticalCapacityFraction,
    float exponent)
{
    const float criticalStorage = std::max(
        0.001f,
        std::max(0.0f, soilMoistureCapacityMm) *
            std::clamp(criticalCapacityFraction, 0.001f, 1.0f));
    const float availableFraction = std::clamp(
        soilMoistureMm / criticalStorage,
        0.0f,
        1.0f);
    return std::pow(availableFraction, std::max(0.0f, exponent));
}

float diagnosticCloudFraction(float relativeHumidity, float cloudOnsetRelativeHumidity)
{
    return smoothStep(
        std::clamp(cloudOnsetRelativeHumidity, 0.0f, 1.0f),
        1.0f,
        std::clamp(relativeHumidity, 0.0f, 1.0f));
}

MoistAdjustment moistSaturationAdjustment(
    float availableColumnWaterMm,
    float saturationCapacityMm,
    float temperatureC,
    float timeStepSeconds,
    float conversionTimeSeconds,
    int iterations,
    float latentHeatingCPerMillimetre,
    float capacityTemperatureSensitivityPerC)
{
    MoistAdjustment result;
    result.remainingVapourMm = std::max(0.0f, availableColumnWaterMm);
    result.adjustedTemperatureC = temperatureC;
    const float initialCapacity = std::max(0.0f, saturationCapacityMm);

    const int iterationCount = std::max(1, iterations);

    for (int iteration = 0; iteration < iterationCount; iteration++)
    {
        const float warming = result.adjustedTemperatureC - temperatureC;
        const float adjustedCapacity = initialCapacity * std::exp(
            std::clamp(capacityTemperatureSensitivityPerC * warming, -20.0f, 20.0f));
        const float condensed = relaxedExcess(
            result.remainingVapourMm,
            adjustedCapacity,
            timeStepSeconds / static_cast<float>(iterationCount),
            conversionTimeSeconds);

        result.condensedMm += condensed;
        result.remainingVapourMm -= condensed;
        result.adjustedTemperatureC += condensed *
            std::max(0.0f, latentHeatingCPerMillimetre);
    }

    result.condensedMm = std::clamp(
        result.condensedMm,
        0.0f,
        std::max(0.0f, availableColumnWaterMm));
    result.remainingVapourMm = std::max(
        0.0f,
        std::max(0.0f, availableColumnWaterMm) - result.condensedMm);
    return result;
}

MoistureLayerExchange exchangeMoistureLayers(
    float boundaryLayerMm,
    float freeTroposphereMm,
    float upwardFraction,
    float downwardFraction)
{
    const float boundary = std::max(0.0f, boundaryLayerMm);
    const float free = std::max(0.0f, freeTroposphereMm);
    const float upward = boundary * std::clamp(upwardFraction, 0.0f, 1.0f);
    const float downward = free * std::clamp(downwardFraction, 0.0f, 1.0f);

    return {
        boundary - upward + downward,
        free + upward - downward,
        upward,
        downward
    };
}

FallingPrecipitation processFallingPrecipitation(
    float condensateMm,
    float surfaceTemperatureC,
    float boundaryRelativeHumidity,
    float maximumReevaporationFraction,
    float maximumVapourUptakeMm,
    float allSnowTemperatureC,
    float allRainTemperatureC)
{
    const float condensate = std::max(0.0f, condensateMm);
    const float humidityDeficit = 1.0f - std::clamp(boundaryRelativeHumidity, 0.0f, 1.0f);
    const float reevaporation = std::min({
        condensate,
        std::max(0.0f, maximumVapourUptakeMm),
        condensate * std::clamp(maximumReevaporationFraction, 0.0f, 1.0f) *
            humidityDeficit * humidityDeficit
    });
    const float reachingSurface = condensate - reevaporation;
    const float rainFraction = smoothStep(
        allSnowTemperatureC,
        allRainTemperatureC,
        surfaceTemperatureC);

    return {
        reachingSurface * rainFraction,
        reachingSurface * (1.0f - rainFraction),
        reevaporation
    };
}

float snowMeltAmount(
    float snowWaterEquivalentMm,
    float surfaceTemperatureC,
    float timeStepSeconds,
    float degreeDayMeltMmPerDegreeC)
{
    constexpr float secondsPerDay = 86400.0f;
    const float potentialMelt = std::max(0.0f, surfaceTemperatureC) *
        std::max(0.0f, degreeDayMeltMmPerDegreeC) *
        std::max(0.0f, timeStepSeconds) / secondsPerDay;
    return std::min(std::max(0.0f, snowWaterEquivalentMm), potentialMelt);
}

SnowAccumulation accumulateSnowfall(
    float snowWaterEquivalentMm,
    float snowfallMm,
    float maximumSnowStorageMm)
{
    const float combined = std::max(0.0f, snowWaterEquivalentMm) +
        std::max(0.0f, snowfallMm);
    const float storage = std::min(combined, std::max(0.0f, maximumSnowStorageMm));
    return { storage, combined - storage };
}

PrecipitationPartition partitionPrecipitation(
    float availableColumnWaterMm,
    float nonOrographicSaturationCapacityMm,
    float terrainAdjustedSaturationCapacityMm,
    float signedMoistureFluxConvergenceMm,
    float surfaceEvaporationMm,
    float surfaceTemperatureC,
    float timeStepSeconds,
    float stratiformCriticalRelativeHumidity,
    float stratiformConversionTimeSeconds,
    float convectiveResidualRelativeHumidity,
    float convectiveConversionEfficiency,
    float convectiveActivationTemperatureC,
    float convectiveFullStrengthTemperatureC)
{
    const float available = std::max(0.0f, availableColumnWaterMm);
    const float criticalRelativeHumidity = std::clamp(
        stratiformCriticalRelativeHumidity, 0.0f, 1.0f);
    const float nonOrographicThreshold = std::max(
        0.0f,
        nonOrographicSaturationCapacityMm * criticalRelativeHumidity);
    const float terrainAdjustedThreshold = std::max(
        0.0f,
        terrainAdjustedSaturationCapacityMm * criticalRelativeHumidity);
    const float nonOrographicCondensation = relaxedExcess(
        available,
        nonOrographicThreshold,
        timeStepSeconds,
        stratiformConversionTimeSeconds);
    const float terrainAdjustedCondensation = relaxedExcess(
        available,
        terrainAdjustedThreshold,
        timeStepSeconds,
        stratiformConversionTimeSeconds);

    PrecipitationPartition result;
    result.orographicMm = std::max(
        0.0f,
        terrainAdjustedCondensation - nonOrographicCondensation);
    result.stratiformMm = std::max(
        0.0f,
        terrainAdjustedCondensation - result.orographicMm);

    const float remaining = std::max(0.0f, available - terrainAdjustedCondensation);
    const float convectiveFloor = std::max(
        0.0f,
        nonOrographicSaturationCapacityMm *
            std::clamp(convectiveResidualRelativeHumidity, 0.0f, 1.0f));
    const float convectivelyAvailable = std::max(0.0f, remaining - convectiveFloor);
    const float instability = smoothStep(
        convectiveActivationTemperatureC,
        convectiveFullStrengthTemperatureC,
        surfaceTemperatureC);
    const float convergentSupply = std::max(
        0.0f,
        signedMoistureFluxConvergenceMm + std::max(0.0f, surfaceEvaporationMm));
    result.convectiveMm = std::min(
        convectivelyAvailable,
        convergentSupply * instability *
            std::clamp(convectiveConversionEfficiency, 0.0f, 1.0f));

    const float total = result.totalMm();
    if (total > available && total > 0.0f)
    {
        const float scale = available / total;
        result.stratiformMm *= scale;
        result.orographicMm *= scale;
        result.convectiveMm *= scale;
    }

    return result;
}

PrecipitationPartition partitionTwoLayerPrecipitation(
    float boundaryLayerWaterMm,
    float freeTroposphereWaterMm,
    float boundaryLayerSaturationCapacityMm,
    float nonOrographicFreeTroposphereCapacityMm,
    float terrainAdjustedFreeTroposphereCapacityMm,
    float signedBoundaryLayerConvergenceMm,
    float surfaceEvaporationMm,
    float surfaceTemperatureC,
    float freeTroposphereTemperatureC,
    float timeStepSeconds,
    float stratiformConversionTimeSeconds,
    float kuoCriticalRelativeHumidity,
    float convectiveConversionEfficiency,
    float convectiveActivationBuoyancyC,
    float convectiveFullStrengthBuoyancyC,
    int moistAdjustmentIterations,
    float latentHeatingCPerMillimetre,
    float capacityTemperatureSensitivityPerC,
    float signedFreeTroposphereConvergenceMm,
    float elevatedMoistureAccessionFraction,
    float kuoHumidityExponent)
{
    const float boundaryWater = std::max(0.0f, boundaryLayerWaterMm);
    const float freeWater = std::max(0.0f, freeTroposphereWaterMm);
    const MoistAdjustment nonOrographicAdjustment = moistSaturationAdjustment(
        freeWater,
        nonOrographicFreeTroposphereCapacityMm,
        freeTroposphereTemperatureC,
        timeStepSeconds,
        stratiformConversionTimeSeconds,
        moistAdjustmentIterations,
        latentHeatingCPerMillimetre,
        capacityTemperatureSensitivityPerC);
    const MoistAdjustment terrainAdjustment = moistSaturationAdjustment(
        freeWater,
        terrainAdjustedFreeTroposphereCapacityMm,
        freeTroposphereTemperatureC,
        timeStepSeconds,
        stratiformConversionTimeSeconds,
        moistAdjustmentIterations,
        latentHeatingCPerMillimetre,
        capacityTemperatureSensitivityPerC);

    PrecipitationPartition result;
    result.orographicMm = std::max(
        0.0f,
        terrainAdjustment.condensedMm - nonOrographicAdjustment.condensedMm);
    result.stratiformMm = std::max(
        0.0f,
        terrainAdjustment.condensedMm - result.orographicMm);

    const float columnCapacity = std::max(
        0.05f,
        boundaryLayerSaturationCapacityMm + nonOrographicFreeTroposphereCapacityMm);
    const float columnRelativeHumidity = std::clamp(
        (boundaryWater + freeWater) / columnCapacity,
        0.0f,
        1.0f);
    const float precipitationEfficiency = kuoPrecipitationEfficiency(
        columnRelativeHumidity,
        kuoCriticalRelativeHumidity,
        kuoHumidityExponent);
    const float buoyancyEfficiency = convectiveBuoyancyEfficiency(
        surfaceTemperatureC - freeTroposphereTemperatureC,
        convectiveActivationBuoyancyC,
        convectiveFullStrengthBuoyancyC);
    const float moistureAccession = std::max(
        0.0f,
        signedBoundaryLayerConvergenceMm + std::max(0.0f, surfaceEvaporationMm) +
            std::max(0.0f, signedFreeTroposphereConvergenceMm) *
                std::clamp(elevatedMoistureAccessionFraction, 0.0f, 1.0f));
    const float convectivelyAvailable = boundaryWater;
    result.convectiveMm = std::min(
        convectivelyAvailable,
        moistureAccession * precipitationEfficiency * buoyancyEfficiency *
            std::clamp(convectiveConversionEfficiency, 0.0f, 1.0f));

    return result;
}
}
