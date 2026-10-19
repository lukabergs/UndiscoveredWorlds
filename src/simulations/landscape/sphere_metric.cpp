#include "simulations/landscape/sphere_metric.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace landscape
{
namespace
{
constexpr double pi = 3.14159265358979323846;
// Offsets by direction code: 0 none, 1 N, 2 NE, 3 E, 4 SE, 5 S, 6 SW, 7 W, 8 NW.
constexpr int offsetsx[9] = { 0, 0, 1, 1, 1, 0, -1, -1, -1 };
constexpr int offsetsy[9] = { 0, -1, -1, 0, 1, 1, 1, 0, -1 };

bool valid(int direction) { return direction >= 1 && direction <= directioncount; }

// The row API of a cell-centred W x H grid on a sphere of radius R.
void fillrows(SphereMetric& metric, int width, int height, double radius)
{
    metric.width = width;
    metric.height = height;
    metric.radius_m = radius;
    metric.dy_m = pi * radius / height;
    const double dlambda = 2.0 * pi / width;
    metric.rowarea_m2.resize(static_cast<std::size_t>(height));
    metric.rowdx_m.resize(static_cast<std::size_t>(height));
    metric.rowlatitude_rad.resize(static_cast<std::size_t>(height));
    for (int j = 0; j < height; j++)
    {
        const double north = pi / 2.0 - pi * j / height;
        const double south = pi / 2.0 - pi * (j + 1) / height;
        const double centre = 0.5 * (north + south);
        const std::size_t row = static_cast<std::size_t>(j);
        metric.rowarea_m2[row] = radius * radius * dlambda * (std::sin(north) - std::sin(south));
        metric.rowdx_m[row] = radius * dlambda * std::cos(centre);
        metric.rowlatitude_rad[row] = centre;
    }
}
}

int directionoffsetx(int direction) { return valid(direction) ? offsetsx[direction] : 0; }
int directionoffsety(int direction) { return valid(direction) ? offsetsy[direction] : 0; }
int oppositedirection(int direction) { return valid(direction) ? (direction + 3) % directioncount + 1 : 0; }

SphereMetric::SphereMetric(int width_, int height_, double radius)
{
    if (width_ < 1 || height_ < 1 || !(radius > 0.0))
        throw std::invalid_argument("SphereMetric needs a positive grid and radius");
    fillrows(*this, width_, height_, radius);
    // The whole sphere: the latitude-longitude grid whose tables the per-cell API reads.
    if (width_ % 2 == 0 && height_ == width_ / 2)
        grid = mesh::makelatlongrid(width_, radius);
}

SphereMetric::SphereMetric(std::shared_ptr<const mesh::CellGrid> cells) : grid(std::move(cells))
{
    if (!grid)
        throw std::invalid_argument("SphereMetric needs a grid");
    radius_m = grid->radius_m();
    if (grid->kind() == mesh::GridKind::latlon)
    {
        if (!grid->issubgrid())
            fillrows(*this, grid->columns(), grid->rows(), radius_m);
        else
            dy_m = pi * radius_m / grid->rows();
    }
    else
        dy_m = pi * radius_m / (2.0 * grid->facecells());
}

double SphereMetric::totalarea_m2() const
{
    if (grid)
        return grid->totalarea_m2();
    double total = 0.0;
    for (double area : rowarea_m2)
        total += area * width;
    return total;
}

std::size_t SphereMetric::neighbour(std::size_t cell, int direction) const
{
    if (grid)
    {
        const std::int32_t next = grid->neighbour(static_cast<mesh::CellId>(cell), direction);
        return next < 0 ? cell : static_cast<std::size_t>(next);
    }
    if (!valid(direction))
        return cell;
    int y = row(cell) + offsetsy[direction];
    if (y < 0 || y >= height)
    {
        if (!wrapsrows)
            return cell;
        y = (y + height) % height;
    }
    const int x = (column(cell) + offsetsx[direction] + width) % width;
    return index(x, y);
}

std::int32_t SphereMetric::neighbourid(std::size_t cell, int direction) const
{
    if (grid)
        return grid->neighbour(static_cast<mesh::CellId>(cell), direction);
    const std::size_t next = neighbour(cell, direction);
    return next == cell ? -1 : static_cast<std::int32_t>(next);
}

int SphereMetric::reverse(std::size_t cell, int direction) const
{
    if (grid)
        return grid->reverse(static_cast<mesh::CellId>(cell), direction);
    return neighbourid(cell, direction) < 0 ? 0 : oppositedirection(direction);
}

double SphereMetric::distance_m(std::size_t cell, int direction) const
{
    if (grid)
        return grid->distance_m(static_cast<mesh::CellId>(cell), direction);
    if (!valid(direction))
        return 0.0;
    const int ox = offsetsx[direction], oy = offsetsy[direction];
    if (ox == 0)
        return dy_m;
    const int y = row(cell);
    if (oy == 0)
        return rowdx_m[static_cast<std::size_t>(y)];
    const int other = row(neighbour(cell, direction));  // y beyond a pole, wrapped on a torus
    return std::hypot(0.5 * (rowdx_m[static_cast<std::size_t>(y)] + rowdx_m[static_cast<std::size_t>(other)]), dy_m);
}

double SphereMetric::facelength_m(std::size_t cell, int direction) const
{
    if (grid)
        return grid->facelength_m(static_cast<mesh::CellId>(cell), direction);
    if (!valid(direction))
        return 0.0;
    const int oy = offsetsy[direction];
    if (oy == 0)
        return dy_m;
    const int y = row(cell);
    const int other = row(neighbour(cell, direction));  // y beyond a pole, wrapped on a torus
    return 0.5 * (rowdx_m[static_cast<std::size_t>(y)] + rowdx_m[static_cast<std::size_t>(other)]);
}

TangentOffset SphereMetric::tangentoffset(std::size_t cell, int direction) const
{
    if (!valid(direction))
        return {};
    if (grid && (grid->kind() == mesh::GridKind::cube || grid->issubgrid()))
    {
        const mesh::CellId c = static_cast<mesh::CellId>(cell);
        const double distance = grid->distance_m(c, direction), bearing = grid->bearing_rad(c, direction);
        return { distance * std::sin(bearing), distance * std::cos(bearing) };
    }
    // The row formulas (a full latitude-longitude grid fills them from its own).
    const int y = row(cell);
    const int ox = offsetsx[direction], oy = offsetsy[direction];
    const double dx = rowdx_m[static_cast<std::size_t>(y)];
    if (oy == 0)
        return { ox * dx, 0.0 };
    if (ox == 0)
        return { 0.0, -oy * dy_m };
    // The neighbour's row: y + oy on a sphere (a caller never asks beyond a pole),
    // wrapped on a torus.
    const int other = row(neighbour(cell, direction));
    const double meandx = 0.5 * (dx + rowdx_m[static_cast<std::size_t>(other)]);
    return { ox * meandx, -oy * dy_m };
}
}
