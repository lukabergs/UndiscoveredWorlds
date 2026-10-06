#include "simulations/landscape/sphere_metric.hpp"

#include <cmath>
#include <stdexcept>

namespace landscape
{
namespace
{
constexpr double pi = 3.14159265358979323846;
// Offsets by direction code: 0 none, 1 N, 2 NE, 3 E, 4 SE, 5 S, 6 SW, 7 W, 8 NW.
constexpr int offsetsx[9] = { 0, 0, 1, 1, 1, 0, -1, -1, -1 };
constexpr int offsetsy[9] = { 0, -1, -1, 0, 1, 1, 1, 0, -1 };

bool valid(int direction) { return direction >= 1 && direction <= directioncount; }
}

int directionoffsetx(int direction) { return valid(direction) ? offsetsx[direction] : 0; }
int directionoffsety(int direction) { return valid(direction) ? offsetsy[direction] : 0; }
int oppositedirection(int direction) { return valid(direction) ? (direction + 3) % directioncount + 1 : 0; }

SphereMetric::SphereMetric(int width_, int height_, double radius)
    : width(width_), height(height_), radius_m(radius)
{
    if (width < 1 || height < 1 || !(radius > 0.0))
        throw std::invalid_argument("SphereMetric needs a positive grid and radius");
    dy_m = pi * radius / height;
    const double dlambda = 2.0 * pi / width;
    rowarea_m2.resize(static_cast<std::size_t>(height));
    rowdx_m.resize(static_cast<std::size_t>(height));
    rowlatitude_rad.resize(static_cast<std::size_t>(height));
    for (int j = 0; j < height; j++)
    {
        const double north = pi / 2.0 - pi * j / height;
        const double south = pi / 2.0 - pi * (j + 1) / height;
        const double centre = 0.5 * (north + south);
        const std::size_t row = static_cast<std::size_t>(j);
        rowarea_m2[row] = radius * radius * dlambda * (std::sin(north) - std::sin(south));
        rowdx_m[row] = radius * dlambda * std::cos(centre);
        rowlatitude_rad[row] = centre;
    }
}

double SphereMetric::totalarea_m2() const
{
    double total = 0.0;
    for (double area : rowarea_m2)
        total += area * width;
    return total;
}

std::size_t SphereMetric::neighbour(std::size_t cell, int direction) const
{
    if (!valid(direction))
        return cell;
    const int y = row(cell) + offsetsy[direction];
    if (y < 0 || y >= height)
        return cell;
    const int x = (column(cell) + offsetsx[direction] + width) % width;
    return index(x, y);
}

double SphereMetric::distance_m(std::size_t cell, int direction) const
{
    if (!valid(direction))
        return 0.0;
    const int ox = offsetsx[direction], oy = offsetsy[direction];
    if (ox == 0)
        return dy_m;
    const int y = row(cell);
    if (oy == 0)
        return rowdx_m[static_cast<std::size_t>(y)];
    const int other = (y + oy < 0 || y + oy >= height) ? y : y + oy;
    return std::hypot(0.5 * (rowdx_m[static_cast<std::size_t>(y)] + rowdx_m[static_cast<std::size_t>(other)]), dy_m);
}

double SphereMetric::facelength_m(std::size_t cell, int direction) const
{
    if (!valid(direction))
        return 0.0;
    const int oy = offsetsy[direction];
    if (oy == 0)
        return dy_m;
    const int y = row(cell);
    const int other = (y + oy < 0 || y + oy >= height) ? y : y + oy;
    return 0.5 * (rowdx_m[static_cast<std::size_t>(y)] + rowdx_m[static_cast<std::size_t>(other)]);
}
}
