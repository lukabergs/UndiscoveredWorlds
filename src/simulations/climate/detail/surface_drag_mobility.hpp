#pragma once

#include "climate_grid.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace climateatmosphere::detail
{
// Parallel surface patches share acceleration, not velocity. In the
// drag-dominated limit U = sqrt(H |F|) / sqrt(Cd), so the coarse velocity
// responds to the area mean of 1/sqrt(Cd). Arithmetic Cd instead assumes
// identical velocity over ocean and forest and overstates mixed-cell stress.
// This scalar homogenization is exact at f=0; it is not the full tensor
// mobility of a rotating, vertically resolved heterogeneous boundary layer.
inline std::vector<float> parallelSurfaceDragCoefficients(
    int sourceColumns, int sourceRows, const std::vector<float>& coefficients,
    int columns, int rows)
{
    if (!climategrid::validGlobalGridDimensions(sourceColumns, sourceRows) ||
        !climategrid::validGlobalGridDimensions(columns, rows) ||
        coefficients.size() != static_cast<std::size_t>(sourceColumns) * sourceRows ||
        std::any_of(coefficients.begin(), coefficients.end(),
            [](float cd) { return !std::isfinite(cd) || cd <= 0.0f; }))
        return {};
    std::vector<float> mobility(coefficients.size());
    for (std::size_t i = 0; i < coefficients.size(); ++i)
        mobility[i] = 1.0f / std::sqrt(coefficients[i]);
    auto result = climategrid::remapField(sourceColumns, sourceRows, mobility, columns, rows);
    for (float& value : result) value = 1.0f / (value * value);
    return result;
}
}
