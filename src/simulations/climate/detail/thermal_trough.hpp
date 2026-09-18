#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace climateatmosphere::detail
{
// Experimental closure input, not a prognostic circulation solver. Temperature
// is simulated; there are no observed winds, named regions or Earth coordinates.
// The returned latitude is bounded around the existing global thermal trough.
inline std::vector<float> regionalThermalEquatorDegrees(
    int columns, int rows, const std::vector<float>& temperature,
    const std::vector<double>& latitudesRadians, float globalEquatorDegrees,
    float hadleyWidthDegrees, float maximumDisplacementFraction = 0.35f)
{
    if (columns < 3 || rows < 3 || temperature.size() != columns * rows ||
        latitudesRadians.size() != rows || !std::isfinite(globalEquatorDegrees) ||
        !std::isfinite(hadleyWidthDegrees) || hadleyWidthDegrees <= 0.0f ||
        !std::isfinite(maximumDisplacementFraction) || maximumDisplacementFraction < 0.0f ||
        !std::all_of(temperature.begin(), temperature.end(), [](float t) { return std::isfinite(t); }))
        throw std::invalid_argument("Invalid regional thermal trough inputs");
    constexpr double degreesPerRadian = 180.0 / 3.141592653589793;
    for (int y = 0; y < rows; ++y)
        if (!std::isfinite(latitudesRadians[y]) ||
            (y && latitudesRadians[y] >= latitudesRadians[y - 1]))
            throw std::invalid_argument("Expected descending finite latitude centres");

    auto field = temperature;
    auto scratch = field;
    // Fixed native-grid smoothing: eight longitude and two latitude binomial
    // passes suppress small warm pockets before locating a broad tropical peak.
    for (int pass = 0; pass < 10; ++pass)
    {
        for (int y = 0; y < rows; ++y)
            for (int x = 0; x < columns; ++x)
            {
                const int cell = y * columns + x;
                const int before = pass < 8 ? y * columns + (x + columns - 1) % columns : std::max(0, y - 1) * columns + x;
                const int after = pass < 8 ? y * columns + (x + 1) % columns : std::min(rows - 1, y + 1) * columns + x;
                scratch[cell] = 0.25f * (field[before] + 2.0f * field[cell] + field[after]);
            }
        field.swap(scratch);
    }
    // Continental summer maxima can lie poleward of the zonal ocean maximum.
    // Search inside the tropical branch, retaining the bounded displacement.
    const float searchLimit = std::min(30.0f, 0.90f * hadleyWidthDegrees);
    const float displacement = maximumDisplacementFraction * hadleyWidthDegrees;
    std::vector<float> target(columns, globalEquatorDegrees);
    for (int x = 0; x < columns; ++x)
    {
        int warm = -1;
        float coldest = std::numeric_limits<float>::max();
        for (int y = 0; y < rows; ++y)
            if (std::abs(latitudesRadians[y] * degreesPerRadian) <= searchLimit)
            {
                coldest = std::min(coldest, field[y * columns + x]);
                if (warm < 0 || field[y * columns + x] > field[warm * columns + x]) warm = y;
            }
        if (warm < 0 || field[warm * columns + x] - coldest < 1.0e-4f) continue;
        double latitude = latitudesRadians[warm] * degreesPerRadian;
        if (warm > 0 && warm + 1 < rows)
        {
            const double before = field[(warm - 1) * columns + x];
            const double centre = field[warm * columns + x];
            const double after = field[(warm + 1) * columns + x];
            const double curvature = before - 2.0 * centre + after;
            if (curvature < -1.0e-6)
                latitude += std::clamp(0.5 * (before - after) / curvature, -0.5, 0.5) *
                    (latitudesRadians[warm + 1] - latitudesRadians[warm]) * degreesPerRadian;
        }
        target[x] = std::clamp(static_cast<float>(latitude), globalEquatorDegrees - displacement, globalEquatorDegrees + displacement);
    }
    auto buffer = target;
    for (int pass = 0; pass < 4; ++pass)
    {
        for (int x = 0; x < columns; ++x)
            buffer[x] = 0.25f * (target[(x + columns - 1) % columns] + 2.0f * target[x] + target[(x + 1) % columns]);
        target.swap(buffer);
    }
    return target;
}
}
