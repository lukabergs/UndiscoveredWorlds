#include "climate_grid.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace climategrid
{
namespace
{
constexpr double pi = 3.1415926535897932384626433832795;

struct GridBounds
{
    std::vector<double> longitudeWest;
    std::vector<double> longitudeEast;
    std::vector<double> latitudeNorth;
    std::vector<double> latitudeSouth;
};

GridBounds makeBounds(int columns, int rows)
{
    GridBounds bounds;
    if (columns <= 0 || rows <= 0)
        return bounds;

    bounds.longitudeWest.resize(columns);
    bounds.longitudeEast.resize(columns);
    const double longitudeStep = 2.0 * pi / static_cast<double>(columns);
    for (int x = 0; x < columns; x++)
    {
        bounds.longitudeWest[x] = static_cast<double>(x) * longitudeStep;
        bounds.longitudeEast[x] = static_cast<double>(x + 1) * longitudeStep;
    }

    bounds.latitudeNorth.resize(rows);
    bounds.latitudeSouth.resize(rows);
    const double latitudeStep = pi / static_cast<double>(rows);
    for (int y = 0; y < rows; y++)
    {
        bounds.latitudeNorth[y] = pi * 0.5 - static_cast<double>(y) * latitudeStep;
        bounds.latitudeSouth[y] = pi * 0.5 - static_cast<double>(y + 1) * latitudeStep;
    }
    return bounds;
}

double latitudeMeasure(double north, double south)
{
    return std::max(0.0, std::sin(north) - std::sin(south));
}

std::vector<std::vector<RemapWeight>> smoothWeights(
    const std::vector<double>& source, const std::vector<double>& target, bool periodic)
{
    const int cells = static_cast<int>(source.size()) - 1;
    const double length = source.back();
    std::vector<double> halfWidth(cells + 1, 0.0);
    for (int face = 1; face < cells; ++face)
        halfWidth[face] = 0.5 * std::min(source[face] - source[face - 1], source[face + 1] - source[face]);
    if (periodic && cells > 1)
        halfWidth.front() = halfWidth.back() = 0.5 * std::min(source[1], length - source[cells - 1]);
    std::vector<std::vector<RemapWeight>> result(target.size() - 1);
    for (std::size_t cell = 0; cell < result.size(); ++cell)
    {
        const double a = target[cell], b = target[cell + 1];
        const auto add = [&](int donor, double integral) {
            if (integral <= 0.0) return;
            auto& weights = result[cell];
            const auto found = std::find_if(weights.begin(), weights.end(),
                [&](const auto& weight) { return weight.source == donor; });
            if (found == weights.end()) weights.push_back({donor, integral / (b - a), 0.0});
            else found->fraction += integral / (b - a);
        };
        const auto ramp = [&](double face, double h, int left, int right) {
            if (h <= 0.0) return;
            const double lo = std::max(a, face - h), hi = std::min(b, face + h);
            if (hi <= lo) return;
            const double rightIntegral = (hi - lo) * (0.5 + (0.5 * (lo + hi) - face) / (2.0 * h));
            add(left, (hi - lo) - rightIntegral);
            add(right, rightIntegral);
        };
        for (int donor = 0; donor < cells; ++donor)
            add(donor, std::max(0.0, std::min(b, source[donor + 1] - halfWidth[donor + 1]) -
                std::max(a, source[donor] + halfWidth[donor])));
        for (int face = 1; face < cells; ++face)
            ramp(source[face], halfWidth[face], face - 1, face);
        if (periodic)
        {
            ramp(0.0, halfWidth.front(), cells - 1, 0);
            ramp(length, halfWidth.back(), cells - 1, 0);
        }
    }
    return result;
}
}

std::size_t SphericalGrid::index(int column, int row) const
{
    return static_cast<std::size_t>(row) * static_cast<std::size_t>(columns) +
        static_cast<std::size_t>(wrapColumn(column));
}

int SphericalGrid::wrapColumn(int column) const
{
    if (columns <= 0)
        return 0;
    const int remainder = column % columns;
    return remainder < 0 ? remainder + columns : remainder;
}

SphericalGrid makeSphericalGrid(int columns, int rows, double radiusMetres)
{
    SphericalGrid grid;
    if (!validGlobalGridDimensions(columns, rows) || radiusMetres <= 0.0)
        return grid;

    grid.columns = columns;
    grid.rows = rows;
    grid.radiusMetres = radiusMetres;
    grid.longitudeSpacingRadians = 2.0 * pi / static_cast<double>(columns);
    grid.latitudeSpacingRadians = pi / static_cast<double>(rows);
    grid.latitudeCentresRadians.resize(rows);
    grid.latitudeNorthFacesRadians.resize(rows);
    grid.latitudeSouthFacesRadians.resize(rows);
    grid.cellAreasSquareMetres.resize(rows);
    grid.zonalFaceLengthsMetres.resize(rows);
    grid.northFaceLengthsMetres.resize(rows);
    grid.southFaceLengthsMetres.resize(rows);

    for (int y = 0; y < rows; y++)
    {
        const double north = pi * 0.5 - static_cast<double>(y) * grid.latitudeSpacingRadians;
        const double south = north - grid.latitudeSpacingRadians;
        const double centre = 0.5 * (north + south);
        grid.latitudeCentresRadians[y] = centre;
        grid.latitudeNorthFacesRadians[y] = north;
        grid.latitudeSouthFacesRadians[y] = south;
        grid.cellAreasSquareMetres[y] = radiusMetres * radiusMetres *
            grid.longitudeSpacingRadians * latitudeMeasure(north, south);
        grid.zonalFaceLengthsMetres[y] = radiusMetres * grid.latitudeSpacingRadians;
        grid.northFaceLengthsMetres[y] = radiusMetres * grid.longitudeSpacingRadians *
            std::max(0.0, std::cos(north));
        grid.southFaceLengthsMetres[y] = radiusMetres * grid.longitudeSpacingRadians *
            std::max(0.0, std::cos(south));
    }
    return grid;
}

bool validGlobalGridDimensions(int columns, int rows)
{
    return columns >= 2 && columns % 2 == 0 && rows == columns / 2;
}

double latitudeCentreRadians(int row, int rows)
{
    if (rows <= 0)
        return 0.0;
    const int boundedRow = std::clamp(row, 0, rows - 1);
    return pi * 0.5 - (static_cast<double>(boundedRow) + 0.5) * pi /
        static_cast<double>(rows);
}

double latitudeBandMeasure(int row, int rows)
{
    if (rows <= 0)
        return 0.0;
    const int boundedRow = std::clamp(row, 0, rows - 1);
    const double step = pi / static_cast<double>(rows);
    const double centre = latitudeCentreRadians(boundedRow, rows);
    return latitudeMeasure(std::min(pi * 0.5, centre + step * 0.5),
        std::max(-pi * 0.5, centre - step * 0.5));
}

ConservativeRemap makeConservativeRemap(
    int sourceColumns, int sourceRows,
    int destinationColumns, int destinationRows)
{
    ConservativeRemap plan;
    if (sourceColumns <= 0 || sourceRows <= 0 || destinationColumns <= 0 || destinationRows <= 0)
        return plan;
    const auto source = makeBounds(sourceColumns, sourceRows);
    const auto destination = makeBounds(destinationColumns, destinationRows);
    plan.sourceColumns = sourceColumns; plan.sourceRows = sourceRows;
    plan.reconstruct = destinationColumns > sourceColumns || destinationRows > sourceRows;
    plan.northSlopeScale.assign(sourceRows, 1.0);
    plan.southSlopeScale.assign(sourceRows, 1.0);
    for (int y = 0; y < sourceRows; ++y)
    {
        const auto centre = [&](int row) {
            return 0.5 * (std::sin(source.latitudeNorth[row]) + std::sin(source.latitudeSouth[row])); };
        const double width = latitudeMeasure(source.latitudeNorth[y], source.latitudeSouth[y]);
        if (y > 0) plan.northSlopeScale[y] = width / (centre(y - 1) - centre(y));
        if (y + 1 < sourceRows) plan.southSlopeScale[y] = width / (centre(y) - centre(y + 1));
    }
    plan.longitude.resize(destinationColumns);
    plan.latitude.resize(destinationRows);
    for (int x = 0; x < destinationColumns; x++)
    {
        const double width = destination.longitudeEast[x] - destination.longitudeWest[x];
        for (int sx = 0; sx < sourceColumns; sx++)
        {
            const double overlap = std::min(destination.longitudeEast[x], source.longitudeEast[sx]) -
                std::max(destination.longitudeWest[x], source.longitudeWest[sx]);
            if (overlap > 1.0e-14)
            {
                const double centroid = 0.5 * (std::min(destination.longitudeEast[x], source.longitudeEast[sx]) +
                    std::max(destination.longitudeWest[x], source.longitudeWest[sx]));
                const double offset = destinationColumns > sourceColumns ?
                    (centroid - 0.5 * (source.longitudeEast[sx] + source.longitudeWest[sx])) /
                        (source.longitudeEast[sx] - source.longitudeWest[sx]) : 0.0;
                plan.longitude[x].push_back({sx, overlap / width, offset});
            }
        }
    }
    for (int y = 0; y < destinationRows; y++)
    {
        const double area = latitudeMeasure(destination.latitudeNorth[y], destination.latitudeSouth[y]);
        for (int sy = 0; sy < sourceRows; sy++)
        {
            const double overlap = latitudeMeasure(
                std::min(destination.latitudeNorth[y], source.latitudeNorth[sy]),
                std::max(destination.latitudeSouth[y], source.latitudeSouth[sy]));
            if (overlap > 1.0e-14)
            {
                // Integrate the linear reconstruction in sin(latitude), the
                // equal-area coordinate. Its source-cell mean stays unchanged.
                const double centroid = 0.5 * (std::sin(std::min(destination.latitudeNorth[y], source.latitudeNorth[sy])) +
                    std::sin(std::max(destination.latitudeSouth[y], source.latitudeSouth[sy])));
                const double sourceCentre = 0.5 * (std::sin(source.latitudeNorth[sy]) + std::sin(source.latitudeSouth[sy]));
                const double offset = destinationRows > sourceRows ? (sourceCentre - centroid) /
                    latitudeMeasure(source.latitudeNorth[sy], source.latitudeSouth[sy]) : 0.0;
                plan.latitude[y].push_back({sy, overlap / area, offset});
            }
        }
    }
    return plan;
}

ConservativeRemap makeSmoothRemap(
    int sourceColumns, int sourceRows, int destinationColumns, int destinationRows)
{
    auto plan = makeConservativeRemap(sourceColumns, sourceRows, destinationColumns, destinationRows);
    if (plan.longitude.empty() || plan.latitude.empty()) return plan;
    plan.reconstruct = false;
    const auto edges = [](int count, bool latitude) {
        std::vector<double> result(count + 1);
        for (int i = 0; i <= count; ++i)
            result[i] = latitude ? 1.0 - std::cos(pi * i / count) : static_cast<double>(i) / count;
        return result;
    };
    if (destinationColumns > sourceColumns)
        plan.longitude = smoothWeights(edges(sourceColumns, false), edges(destinationColumns, false), true);
    if (destinationRows > sourceRows)
        plan.latitude = smoothWeights(edges(sourceRows, true), edges(destinationRows, true), false);
    return plan;
}

std::vector<float> remapSmoothField(
    int sourceColumns, int sourceRows, const std::vector<float>& source,
    int destinationColumns, int destinationRows)
{
    if (sourceColumns <= 0 || sourceRows <= 0 || source.size() !=
        static_cast<std::size_t>(sourceColumns) * sourceRows) return {};
    return remapField(makeSmoothRemap(sourceColumns, sourceRows, destinationColumns, destinationRows),
        source, sourceColumns);
}

std::vector<float> remapField(
    int sourceColumns,
    int sourceRows,
    const std::vector<float>& source,
    int destinationColumns,
    int destinationRows)
{
    const std::size_t sourceCount = static_cast<std::size_t>(std::max(0, sourceColumns)) *
        static_cast<std::size_t>(std::max(0, sourceRows));
    if (sourceColumns <= 0 || sourceRows <= 0 || destinationColumns <= 0 ||
        destinationRows <= 0 || source.size() != sourceCount)
    {
        return {};
    }

    const auto plan = makeConservativeRemap(
        sourceColumns, sourceRows, destinationColumns, destinationRows);
    return remapField(plan, source, sourceColumns);
}

std::vector<float> remapField(const ConservativeRemap& plan,
    const std::vector<float>& source, int sourceColumns)
{
    if (sourceColumns <= 0 || source.empty()) return {};
    for (const auto& row : plan.latitude)
        for (const auto& yw : row)
            if (yw.source < 0 || static_cast<std::size_t>(yw.source + 1) * sourceColumns > source.size()) return {};
    for (const auto& column : plan.longitude)
        for (const auto& xw : column)
            if (xw.source < 0 || xw.source >= sourceColumns) return {};
    const int columns = static_cast<int>(plan.longitude.size());
    const int rows = static_cast<int>(plan.latitude.size());
    std::vector<float> result(static_cast<std::size_t>(columns) * rows);
    for (int y = 0; y < rows; ++y)
        for (int x = 0; x < columns; ++x)
            result[y * columns + x] = sampleConservative(plan, x, y,
                [&](int xx, int yy) { return source[yy * sourceColumns + xx]; });
    return result;
}

double areaWeightedIntegral(
    int columns,
    int rows,
    const std::vector<float>& field,
    double radiusMetres)
{
    const std::size_t count = static_cast<std::size_t>(std::max(0, columns)) *
        static_cast<std::size_t>(std::max(0, rows));
    if (columns <= 0 || rows <= 0 || radiusMetres <= 0.0 || field.size() != count)
        return 0.0;

    const GridBounds bounds = makeBounds(columns, rows);
    const double longitudeStep = 2.0 * pi / static_cast<double>(columns);
    double integral = 0.0;
    for (int y = 0; y < rows; y++)
    {
        const double cellArea = radiusMetres * radiusMetres * longitudeStep *
            latitudeMeasure(bounds.latitudeNorth[y], bounds.latitudeSouth[y]);
        for (int x = 0; x < columns; x++)
            integral += cellArea * static_cast<double>(field[
                static_cast<std::size_t>(y) * columns + x]);
    }
    return integral;
}
}
