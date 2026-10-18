#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <stdexcept>

#include "simulations/society/social_generation_internal.hpp"

namespace socialgen
{
namespace
{
void latlonoffsets(const mesh::CellGrid& grid, mesh::CellId a, mesh::CellId b, int& dx, int& dy)
{
    if (grid.kind() != mesh::GridKind::latlon)
        throw std::logic_error("society: grid-step distances exist on lat-lon only");
    const int direct = std::abs(grid.column(a) - grid.column(b));
    dx = std::min(direct, grid.columns() - direct);
    dy = std::abs(grid.row(a) - grid.row(b));
}
}

int exportcolumns(const mesh::CellGrid& grid)
{
    return grid.kind() == mesh::GridKind::latlon ? grid.columns() : 4 * grid.facecells();
}

float celldistance(const mesh::CellGrid& grid, mesh::CellId a, mesh::CellId b)
{
    int dx = 0;
    int dy = 0;
    latlonoffsets(grid, a, b, dx, dy);
    return std::sqrt(static_cast<float>(dx * dx + dy * dy));
}

int cellmanhattan(const mesh::CellGrid& grid, mesh::CellId a, mesh::CellId b)
{
    int dx = 0;
    int dy = 0;
    latlonoffsets(grid, a, b, dx, dy);
    return dx + dy;
}

int celldistancesquared(const mesh::CellGrid& grid, mesh::CellId a, mesh::CellId b)
{
    int dx = 0;
    int dy = 0;
    latlonoffsets(grid, a, b, dx, dy);
    return dx * dx + dy * dy;
}

void placesettlement(Settlement& settlement, const mesh::CellGrid& grid, mesh::CellId cell)
{
    settlement.cell = cell;

    if (grid.kind() == mesh::GridKind::latlon)
    {
        settlement.x = grid.column(cell);
        settlement.y = grid.row(cell);
        return;
    }

    const double pi = 3.14159265358979323846;
    const int columns = exportcolumns(grid);
    const int rows = columns / 2;
    const double longitude = grid.longitude_rad(cell);
    const double latitude = grid.latitude_rad(cell);
    settlement.x = std::clamp(static_cast<int>(std::floor((longitude + pi) / (2.0 * pi) * columns)), 0, columns - 1);
    settlement.y = std::clamp(static_cast<int>(std::floor((0.5 * pi - latitude) / pi * rows)), 0, rows - 1);
}
}
