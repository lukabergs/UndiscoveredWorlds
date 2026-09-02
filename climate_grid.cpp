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

GridBounds makeBounds(int columns, int rows, LatitudeLayout layout)
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
    if (layout == LatitudeLayout::cellCentred || rows == 1)
    {
        const double latitudeStep = pi / static_cast<double>(rows);
        for (int y = 0; y < rows; y++)
        {
            bounds.latitudeNorth[y] = pi * 0.5 - static_cast<double>(y) * latitudeStep;
            bounds.latitudeSouth[y] = pi * 0.5 - static_cast<double>(y + 1) * latitudeStep;
        }
        return bounds;
    }

    const double sampleStep = pi / static_cast<double>(rows - 1);
    for (int y = 0; y < rows; y++)
    {
        const double centre = pi * 0.5 - static_cast<double>(y) * sampleStep;
        bounds.latitudeNorth[y] = y == 0 ? pi * 0.5 : centre + sampleStep * 0.5;
        bounds.latitudeSouth[y] = y == rows - 1 ? -pi * 0.5 : centre - sampleStep * 0.5;
    }
    return bounds;
}

double latitudeMeasure(double north, double south)
{
    return std::max(0.0, std::sin(north) - std::sin(south));
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
    if (columns <= 0 || rows <= 0 || radiusMetres <= 0.0)
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

double latitudeCentreRadians(int row, int rows, LatitudeLayout layout)
{
    if (rows <= 0)
        return 0.0;
    const int boundedRow = std::clamp(row, 0, rows - 1);
    if (layout == LatitudeLayout::poleInclusive && rows > 1)
    {
        return pi * 0.5 - static_cast<double>(boundedRow) * pi /
            static_cast<double>(rows - 1);
    }
    return pi * 0.5 - (static_cast<double>(boundedRow) + 0.5) * pi /
        static_cast<double>(rows);
}

double latitudeBandMeasure(int row, int rows, LatitudeLayout layout)
{
    if (rows <= 0)
        return 0.0;
    const int boundedRow = std::clamp(row, 0, rows - 1);
    const double step = pi / (layout == LatitudeLayout::poleInclusive && rows > 1
        ? rows - 1 : rows);
    const double centre = latitudeCentreRadians(boundedRow, rows, layout);
    return latitudeMeasure(std::min(pi * 0.5, centre + step * 0.5),
        std::max(-pi * 0.5, centre - step * 0.5));
}

ConservativeRemap makeConservativeRemap(
    int sourceColumns, int sourceRows, LatitudeLayout sourceLayout,
    int destinationColumns, int destinationRows, LatitudeLayout destinationLayout)
{
    ConservativeRemap plan;
    if (sourceColumns <= 0 || sourceRows <= 0 || destinationColumns <= 0 || destinationRows <= 0)
        return plan;
    const auto source = makeBounds(sourceColumns, sourceRows, sourceLayout);
    const auto destination = makeBounds(destinationColumns, destinationRows, destinationLayout);
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
                plan.longitude[x].push_back({sx, overlap / width});
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
                plan.latitude[y].push_back({sy, overlap / area});
        }
    }
    return plan;
}

std::vector<float> remapField(
    int sourceColumns,
    int sourceRows,
    LatitudeLayout sourceLayout,
    const std::vector<float>& source,
    int destinationColumns,
    int destinationRows,
    LatitudeLayout destinationLayout)
{
    const std::size_t sourceCount = static_cast<std::size_t>(std::max(0, sourceColumns)) *
        static_cast<std::size_t>(std::max(0, sourceRows));
    const std::size_t destinationCount =
        static_cast<std::size_t>(std::max(0, destinationColumns)) *
        static_cast<std::size_t>(std::max(0, destinationRows));
    if (sourceColumns <= 0 || sourceRows <= 0 || destinationColumns <= 0 ||
        destinationRows <= 0 || source.size() != sourceCount)
    {
        return {};
    }

    const auto plan = makeConservativeRemap(sourceColumns, sourceRows, sourceLayout,
        destinationColumns, destinationRows, destinationLayout);
    std::vector<float> destination(destinationCount, 0.0f);

    for (int destinationY = 0; destinationY < destinationRows; destinationY++)
    {
        for (int destinationX = 0; destinationX < destinationColumns; destinationX++)
        {
            double weightedValue = 0.0;
            for (const auto& yw : plan.latitude[destinationY])
            {
                for (const auto& xw : plan.longitude[destinationX])
                {
                    weightedValue += yw.fraction * xw.fraction * static_cast<double>(source[
                        static_cast<std::size_t>(yw.source) * sourceColumns + xw.source]);
                }
            }
            destination[static_cast<std::size_t>(destinationY) * destinationColumns +
                destinationX] = static_cast<float>(weightedValue);
        }
    }
    return destination;
}

double areaWeightedIntegral(
    int columns,
    int rows,
    LatitudeLayout layout,
    const std::vector<float>& field,
    double radiusMetres)
{
    const std::size_t count = static_cast<std::size_t>(std::max(0, columns)) *
        static_cast<std::size_t>(std::max(0, rows));
    if (columns <= 0 || rows <= 0 || radiusMetres <= 0.0 || field.size() != count)
        return 0.0;

    const GridBounds bounds = makeBounds(columns, rows, layout);
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
