#pragma once

#include <cstddef>
#include <vector>

namespace climategrid
{
enum class LatitudeLayout
{
    cellCentred,
    poleInclusive
};

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
double latitudeCentreRadians(int row, int rows, LatitudeLayout layout);
double latitudeBandMeasure(int row, int rows, LatitudeLayout layout);

struct RemapWeight
{
    int source = 0;
    double fraction = 0.0;
};

// Separable, normalized overlap weights. Reuse for multiple input fields.
struct ConservativeRemap
{
    std::vector<std::vector<RemapWeight>> longitude;
    std::vector<std::vector<RemapWeight>> latitude;
};
ConservativeRemap makeConservativeRemap(
    int sourceColumns, int sourceRows, LatitudeLayout sourceLayout,
    int destinationColumns, int destinationRows, LatitudeLayout destinationLayout);

std::vector<float> remapField(
    int sourceColumns,
    int sourceRows,
    LatitudeLayout sourceLayout,
    const std::vector<float>& source,
    int destinationColumns,
    int destinationRows,
    LatitudeLayout destinationLayout);

double areaWeightedIntegral(
    int columns,
    int rows,
    LatitudeLayout layout,
    const std::vector<float>& field,
    double radiusMetres = 1.0);
}
