#pragma once

#include <cstddef>
#include <algorithm>
#include <cmath>
#include <vector>

namespace climategrid
{
struct SphericalGrid
{
    int columns = 0;
    int rows = 0;
    double radiusMetres = 0.0;
    double longitudeSpacingRadians = 0.0;
    double latitudeSpacingRadians = 0.0;
    std::vector<double> latitudeCentresRadians;
    std::vector<double> latitudeNorthFacesRadians;
    std::vector<double> latitudeSouthFacesRadians;
    std::vector<double> cellAreasSquareMetres;
    std::vector<double> zonalFaceLengthsMetres;
    std::vector<double> northFaceLengthsMetres;
    std::vector<double> southFaceLengthsMetres;

    std::size_t index(int column, int row) const;
    int wrapColumn(int column) const;
};

SphericalGrid makeSphericalGrid(int columns, int rows, double radiusMetres);
bool validGlobalGridDimensions(int columns, int rows);
double latitudeCentreRadians(int row, int rows);
double latitudeBandMeasure(int row, int rows);

struct RemapWeight
{
    int source = 0;
    double fraction = 0.0;
    // Overlap centroid relative to source-cell centroid, in cell widths.
    double offset = 0.0;
};

// Separable, normalized overlap weights. Reuse for multiple input fields.
struct ConservativeRemap
{
    std::vector<std::vector<RemapWeight>> longitude;
    std::vector<std::vector<RemapWeight>> latitude;
    int sourceColumns = 0, sourceRows = 0;
    bool reconstruct = false;
    std::vector<double> northSlopeScale, southSlopeScale;
};
ConservativeRemap makeConservativeRemap(
    int sourceColumns, int sourceRows,
    int destinationColumns, int destinationRows);

// Continuous, bounded reconstruction for refined climate fields. Symmetric
// transitions across cell faces conserve the global area integral, but allow
// exchange between adjacent source cells. Restriction still uses exact overlaps.
ConservativeRemap makeSmoothRemap(
    int sourceColumns, int sourceRows,
    int destinationColumns, int destinationRows);
std::vector<float> remapSmoothField(
    int sourceColumns, int sourceRows, const std::vector<float>& source,
    int destinationColumns, int destinationRows);

// Apply a reusable overlap plan to any cell-centred storage layout.
template<typename Accessor>
float sampleConservative(const ConservativeRemap& plan, int x, int y, Accessor at)
{
    double value = 0.0;
    for (const auto& yw : plan.latitude[y])
        for (const auto& xw : plan.longitude[x])
        {
            const int sx = xw.source, sy = yw.source;
            const double centre = at(sx, sy);
            double reconstructed = centre;
            if (plan.reconstruct)
            {
                const double west = at((sx + plan.sourceColumns - 1) % plan.sourceColumns, sy);
                const double east = at((sx + 1) % plan.sourceColumns, sy);
                const double north = at(sx, (std::max)(0, sy - 1));
                const double south = at(sx, (std::min)(plan.sourceRows - 1, sy + 1));
                const auto minmod = [](double a, double b) {
                    return a * b <= 0.0 ? 0.0 : std::copysign((std::min)(std::abs(a), std::abs(b)), a); };
                const double dx = minmod(centre - west, east - centre);
                const double dy = minmod((centre - north) * plan.northSlopeScale[sy],
                    (south - centre) * plan.southSlopeScale[sy]);
                const double extent = 0.5 * (std::abs(dx) + std::abs(dy));
                const double bound = (std::min)(centre - (std::min)({centre, west, east, north, south}),
                    (std::max)({centre, west, east, north, south}) - centre);
                const double limiter = extent > 0.0 ? (std::min)(1.0, bound / extent) : 0.0;
                reconstructed += limiter * (dx * xw.offset + dy * yw.offset);
            }
            value += yw.fraction * xw.fraction * reconstructed;
        }
    return static_cast<float>(value);
}

std::vector<float> remapField(
    int sourceColumns,
    int sourceRows,
    const std::vector<float>& source,
    int destinationColumns,
    int destinationRows);

std::vector<float> remapField(const ConservativeRemap& plan,
    const std::vector<float>& source, int sourceColumns);

double areaWeightedIntegral(
    int columns,
    int rows,
    const std::vector<float>& field,
    double radiusMetres = 1.0);
}
