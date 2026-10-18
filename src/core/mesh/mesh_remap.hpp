#pragma once

// Remap between cell grids and point sampling (docs/architecture/sphere-mesh.md,
// section 5). makeremap builds sparse rows per destination cell, entries
// (source id, overlap area in m^2) in ascending source order:
// - full lat-lon to full lat-lon: exact overlaps, separable in longitude and
//   sin(latitude);
// - any other pair: each source cell is subdivided in its native coordinates
//   (face a, b or longitude, sin latitude) as a quadtree; a piece whose four
//   corners and centre locate in one destination cell is assigned whole with its
//   exact area, others split down to RemapOptions::depth; a leaf goes to the cell
//   of its centre. Each source cell's leaves are scaled to its exact area, so
//   column sums are exact; a row sums to its area within the band of the leaves.
//   Source cells run in parallel chunks joined in source order (independent of
//   the thread count).
// Direction-coded fields (D8 codes) are never remapped: they are recomputed on
// the destination grid.

#include "core/mesh/cell_grid.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace mesh
{
struct RemapOptions
{
    int depth = 6;  // quadtree depth limit: leaves 1 / 2^depth of the source cell's edge
};

struct RemapWeights
{
    std::size_t sourcecount = 0;
    std::size_t destinationcount = 0;
    std::vector<std::size_t> rowstart;        // destinationcount + 1 offsets into source and area_m2
    std::vector<CellId> source;               // per entry, ascending within a row
    std::vector<double> area_m2;              // per entry, the overlap area
    std::vector<double> sourcearea_m2;        // per source cell, the area its entries distribute
    std::size_t emptyrows = 0;                // destination rows without weight (filled by sampling)
};

RemapWeights makeremap(const CellGrid& source, const CellGrid& destination, RemapOptions options = {});

// Up to four cells and weights (summing to 1) of the bilinear sample at a point:
// on a cube face bilinear in (a, b) between the four surrounding centres, cells
// across a face edge from the neighbour table, at a cube corner the missing
// fourth cell's weight shared by the three in proportion (the mean of the three
// at the corner itself); on lat-lon bilinear in longitude
// (wrapped) and latitude (polar rows clamped). On a subgrid, cells outside it
// are dropped and the weights renormalised (count 0 when none is inside).
struct BilinearStencil
{
    int count = 0;
    CellId cell[4] = { 0, 0, 0, 0 };
    double weight[4] = { 0.0, 0.0, 0.0, 0.0 };
};
BilinearStencil bilinearstencil(const CellGrid& grid, const Vec3& p);

template <typename T>
T samplebilinear(const CellGrid& grid, const std::vector<T>& field, const Vec3& p)
{
    const BilinearStencil stencil = bilinearstencil(grid, p);
    if (stencil.count == 0)
        return T{};
    // Relative to the first value, so a constant field samples exactly.
    const double first = static_cast<double>(field[stencil.cell[0]]);
    double sum = 0.0;
    for (int n = 1; n < stencil.count; n++)
        sum += stencil.weight[n] * (static_cast<double>(field[stencil.cell[n]]) - first);
    return static_cast<T>(first + sum);
}

// Intensive fields: the area-weighted mean over each row (relative to the row's
// first value, so a constant field stays exact); a row without weight takes the
// bilinear sample of the source at the destination centre (RemapWeights::emptyrows).
template <typename T>
std::vector<T> remapmean(const RemapWeights& weights, const CellGrid& source, const CellGrid& destination,
    const std::vector<T>& values)
{
    if (values.size() != weights.sourcecount || destination.cellcount() != weights.destinationcount)
        throw std::invalid_argument("remapmean: field or grid does not match the weights");
    std::vector<T> result(weights.destinationcount);
    for (std::size_t d = 0; d < weights.destinationcount; d++)
    {
        const std::size_t begin = weights.rowstart[d], end = weights.rowstart[d + 1];
        if (begin == end)
        {
            result[d] = samplebilinear(source, values, destination.centre(static_cast<CellId>(d)));
            continue;
        }
        const double first = static_cast<double>(values[weights.source[begin]]);
        double weighted = 0.0, total = 0.0;
        for (std::size_t e = begin; e < end; e++)
        {
            weighted += weights.area_m2[e] * (static_cast<double>(values[weights.source[e]]) - first);
            total += weights.area_m2[e];
        }
        result[d] = static_cast<T>(first + weighted / total);
    }
    return result;
}

// Extensive fields: each source total split by its overlap shares; the global total is kept.
// The shares are summed in double; an integer type takes the sum rounded to nearest, so an
// integer total is exact only up to that rounding of each destination cell.
template <typename T>
std::vector<T> remaptotal(const RemapWeights& weights, const std::vector<T>& values)
{
    if (values.size() != weights.sourcecount)
        throw std::invalid_argument("remaptotal: field does not match the weights");
    std::vector<T> result(weights.destinationcount);
    for (std::size_t d = 0; d < weights.destinationcount; d++)
    {
        double total = 0.0;
        for (std::size_t e = weights.rowstart[d]; e < weights.rowstart[d + 1]; e++)
        {
            const CellId s = weights.source[e];
            total += static_cast<double>(values[s]) * (weights.area_m2[e] / weights.sourcearea_m2[s]);
        }
        if constexpr (std::is_integral_v<T>)
            result[d] = static_cast<T>(std::llround(total));
        else
            result[d] = static_cast<T>(total);
    }
    return result;
}

// Codes, ids and masks: the source cell at each destination centre (T{} where it
// is outside a source subgrid).
template <typename T>
std::vector<T> remapnearest(const CellGrid& source, const CellGrid& destination, const std::vector<T>& values)
{
    if (values.size() != source.cellcount())
        throw std::invalid_argument("remapnearest: field does not match the source grid");
    std::vector<T> result(destination.cellcount());
    for (std::size_t d = 0; d < result.size(); d++)
    {
        const CellId s = source.locate(destination.centre(static_cast<CellId>(d)));
        result[d] = s == nocell ? T{} : values[s];
    }
    return result;
}

// Codes over larger destination cells: the code with the largest overlap, ties
// to the smaller code; a row without weight takes the nearest source cell.
template <typename T>
std::vector<T> remapmajority(const RemapWeights& weights, const CellGrid& source, const CellGrid& destination,
    const std::vector<T>& values)
{
    if (values.size() != weights.sourcecount || destination.cellcount() != weights.destinationcount)
        throw std::invalid_argument("remapmajority: field or grid does not match the weights");
    std::vector<T> result(weights.destinationcount);
    std::vector<std::pair<T, double>> tally;
    for (std::size_t d = 0; d < weights.destinationcount; d++)
    {
        tally.clear();
        for (std::size_t e = weights.rowstart[d]; e < weights.rowstart[d + 1]; e++)
        {
            const T code = values[weights.source[e]];
            auto found = std::find_if(tally.begin(), tally.end(), [&](const auto& entry) { return entry.first == code; });
            if (found == tally.end())
                tally.emplace_back(code, weights.area_m2[e]);
            else
                found->second += weights.area_m2[e];
        }
        if (tally.empty())
        {
            const CellId s = source.locate(destination.centre(static_cast<CellId>(d)));
            result[d] = s == nocell ? T{} : values[s];
            continue;
        }
        std::size_t best = 0;
        for (std::size_t n = 1; n < tally.size(); n++)
            if (tally[n].second > tally[best].second ||
                (tally[n].second == tally[best].second && tally[n].first < tally[best].first))
                best = n;
        result[d] = tally[best].first;
    }
    return result;
}

struct EastNorth
{
    double east = 0.0, north = 0.0;
};
// A 3-D tangent vector at cell c to its east and north components, and back.
EastNorth tangenttoeastnorth(const CellGrid& grid, CellId c, const Vec3& vector);
Vec3 eastnorthtotangent(const CellGrid& grid, CellId c, const EastNorth& components);
}
