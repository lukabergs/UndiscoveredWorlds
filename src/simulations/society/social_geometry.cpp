#include <algorithm>
#include <cmath>
#include <cstdint>

#include "core/deterministic_random.hpp"
#include "core/mesh/mesh_geometry.hpp"
#include "simulations/society/social_generation_internal.hpp"

namespace socialgen
{
namespace
{
constexpr double pi = 3.14159265358979323846;
}

int exportcolumns(const mesh::CellGrid& grid)
{
    return grid.kind() == mesh::GridKind::latlon ? grid.columns() : 4 * grid.facecells();
}

double cellwidth_m(const mesh::CellGrid& grid)
{
    return 2.0 * pi * grid.radius_m() / static_cast<double>(exportcolumns(grid));
}

float celldistance(const mesh::CellGrid& grid, mesh::CellId a, mesh::CellId b)
{
    return static_cast<float>(mesh::greatcircle_m(grid.centre(a), grid.centre(b), grid.radius_m()) / cellwidth_m(grid));
}

double exportcellarea_m2(const mesh::CellGrid& grid)
{
    const double columns = static_cast<double>(exportcolumns(grid));
    return grid.totalarea_m2() / (columns * (columns / 2.0));
}

std::vector<double> sitedrawkeys(const mesh::CellGrid& grid, long seed)
{
    const std::uint64_t base = deterministiccontextseed(seed, 0x53495445);
    const double meanarea = grid.totalarea_m2() / static_cast<double>(grid.cellcount());
    std::vector<double> keys(grid.cellcount());

    for (mesh::CellId c = 0; c < grid.cellcount(); c++)
    {
        // Efraimidis-Spirakis: u^(1 / w) for weight w, compared as log(u) / w; u in (0, 1).
        const double u = (static_cast<double>(deterministiccombine(base, static_cast<std::int64_t>(c)) >> 11) + 0.5) * (1.0 / 9007199254740992.0);
        keys[c] = std::log(u) / (grid.area_m2(c) / meanarea);
    }

    return keys;
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

    const int columns = exportcolumns(grid);
    const int rows = columns / 2;
    const double longitude = grid.longitude_rad(cell);
    const double latitude = grid.latitude_rad(cell);
    settlement.x = std::clamp(static_cast<int>(std::floor((longitude + pi) / (2.0 * pi) * columns)), 0, columns - 1);
    settlement.y = std::clamp(static_cast<int>(std::floor((0.5 * pi - latitude) / pi * rows)), 0, rows - 1);
}
}
