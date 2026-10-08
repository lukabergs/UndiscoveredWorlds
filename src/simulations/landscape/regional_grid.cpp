#include "simulations/landscape/regional_grid.hpp"

#include <cmath>
#include <stdexcept>

#include "simulations/landscape/sphere_metric.hpp"

// The window functions are implemented in erosion plan step E2a
// (docs/architecture/erosion-hydrology-plan.md, Phase E); the patch metric below is
// part of the contracts (E0), so that E2c can use it while E2a runs.

namespace landscape
{
namespace
{
constexpr double pi = 3.14159265358979323846;
}

std::size_t PatchMetric::neighbour(std::size_t cell, int direction) const
{
    if (direction < 1 || direction > directioncount)
        return cell;
    const int x = static_cast<int>(cell % static_cast<std::size_t>(width)) + directionoffsetx(direction);
    const int y = static_cast<int>(cell / static_cast<std::size_t>(width)) + directionoffsety(direction);
    if (x < 0 || x >= width || y < 0 || y >= height)
        return cell;
    return index(x, y);
}

double PatchMetric::distance_m(std::size_t cell, int direction) const
{
    if (direction < 1 || direction > directioncount)
        return 0.0;
    const int ox = directionoffsetx(direction), oy = directionoffsety(direction);
    if (ox == 0)
        return dy_m;
    const std::size_t y = cell / static_cast<std::size_t>(width);
    if (oy == 0)
        return rowdx_m[y];
    const std::size_t other = neighbour(cell, direction) / static_cast<std::size_t>(width);  // the cell's own row beyond the edge
    return std::hypot(0.5 * (rowdx_m[y] + rowdx_m[other]), dy_m);
}

PatchMetric finepatchmetric(const RegionalWindow& window, double radius_m)
{
    if (window.scale < 1 || window.worldwidth < 1 || window.worldheight < 1 || !(radius_m > 0.0))
        throw std::invalid_argument("finepatchmetric needs a window with a positive grid and scale, and a positive radius");
    PatchMetric metric;
    metric.width = window.domainfinewidth();
    metric.height = window.domainfineheight();
    metric.radius_m = radius_m;
    const int fineheight = window.worldheight * window.scale;
    const double dlambda = 2.0 * pi / (static_cast<double>(window.worldwidth) * window.scale);
    metric.dy_m = pi * radius_m / fineheight;
    const std::size_t rows = static_cast<std::size_t>(metric.height);
    metric.rowdx_m.resize(rows);
    metric.rowarea_m2.resize(rows);
    metric.rowlatitude_rad.resize(rows);
    for (int y = 0; y < metric.height; y++)
    {
        // The same formulas as SphereMetric for row J of the planet-wide fine grid.
        const int j = window.domainrow0 * window.scale + y;
        const double north = pi / 2.0 - pi * j / fineheight;
        const double south = pi / 2.0 - pi * (j + 1) / fineheight;
        const double centre = 0.5 * (north + south);
        const std::size_t row = static_cast<std::size_t>(y);
        metric.rowarea_m2[row] = radius_m * radius_m * dlambda * (std::sin(north) - std::sin(south));
        metric.rowdx_m[row] = radius_m * dlambda * std::cos(centre);
        metric.rowlatitude_rad[row] = centre;
    }
    return metric;
}
}
