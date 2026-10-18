#pragma once

// The metric the landscape kernels walk: a view of a mesh::CellGrid
// (docs/architecture/sphere-mesh.md, section 3), or, without one, the planar
// torus of the tectonics grid and the hand-built window metrics.
//
// - A grid-backed metric (`grid` set) reads every per-cell quantity from the
//   grid's tables: size, area_m2, neighbour, neighbourid, reverse, distance_m,
//   facelength_m. The grid is a cube or a latitude-longitude grid; kernels use only
//   this per-cell API, so they run on either. On a cube the slots are the cell's
//   face frame (N toward v, E toward u, also across face edges), a reverse slot is
//   reverse(cell, k), never oppositedirection(k), and the one diagonal of a corner
//   cell that would point at a cube corner has no neighbour.
// - The (W, H, R) constructor with H = W / 2 builds the cell-centred W x W/2
//   latitude-longitude grid (mesh::makelatlongrid, whose tables reproduce the
//   formulas below expression for expression) and fills the row API from the same
//   formulas. Row j spans latitudes phi_n = 90 - 180 j / H to phi_s = phi_n - 180 / H;
//   its cell area is R^2 dlambda (sin phi_n - sin phi_s) with dlambda = 2 pi / W;
//   its east-west spacing is R dlambda cos(phi) at the row centre; the north-south
//   spacing is pi R / H everywhere; a diagonal is hypot(mean dx of the two rows, dy).
//   Longitude wraps; a neighbour beyond a pole does not exist (neighbourid -1,
//   neighbour the cell itself). Any other H gives a metric without a grid (a
//   window of rows) that computes the same formulas from its row fields.
// - The row API (width, height, rowarea_m2, rowdx_m, rowlatitude_rad, index, row,
//   column) exists on latitude-longitude and gridless metrics only; a cube metric
//   has width = height = 0 and empty rows, and dy_m is its nominal spacing
//   pi R / (2N), the equatorial cell's (equal to a W = 4N lat-lon grid's dy).
//   A gridless metric whose row fields are changed after construction keeps the
//   change; a grid-backed one reads the grid and does not see it.
// - A torus metric (torusmetric, erosion plan step D2a: the tectonics grid) is a
//   planar W x H grid of square cells that wraps rows too (no pole), without a grid:
//   every cell has area cell^2, dx = dy = cell, latitude 0 and radius 0, so a
//   diagonal is cell sqrt(2) and totalarea_m2 is W H cell^2. It stays until the
//   engine drops its torus.
// Planet-free: no world class, no tuning header.

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

#include "core/mesh/cell_grid.hpp"

namespace landscape
{
// Direction codes are those of hydrosphere.river_direction: 1 north, 2 north-east,
// 3 east, 4 south-east, 5 south, 6 south-west, 7 west, 8 north-west; 0 is none
// (an outlet, or a cell that is its own receiver). They are the mesh::CellGrid
// slots. The offsets and the opposite code are latitude-longitude (and torus)
// geometry: on a cube use SphereMetric::reverse.
constexpr int directioncount = 8;
int directionoffsetx(int direction);   // -1, 0 or 1 (east positive)
int directionoffsety(int direction);   // -1, 0 or 1 (south positive; row 0 is the north pole)
int oppositedirection(int direction);  // 1 <-> 5, 2 <-> 6, 3 <-> 7, 4 <-> 8; 0 -> 0

// The offset from a cell centre to a neighbour's on the cell's tangent plane, in metres.
struct TangentOffset
{
    double east = 0.0;
    double north = 0.0;
};

struct SphereMetric
{
    int width = 0;                        // W columns (longitude); 0 on a cube
    int height = 0;                       // H = W / 2 rows (latitude); 0 on a cube
    double radius_m = 0.0;
    double dy_m = 0.0;                    // north-south spacing, pi R / H (a cube: pi R / (2N))
    std::vector<double> rowarea_m2;       // cell area per row
    std::vector<double> rowdx_m;          // east-west spacing at the row centre, per row
    std::vector<double> rowlatitude_rad;  // row-centre latitude, per row (positive north)
    bool wrapsrows = false;               // a torus: row -1 is row H - 1 and row H is row 0
    // The cell grid the per-cell API reads; null on a torus and on a window of rows.
    std::shared_ptr<const mesh::CellGrid> grid;

    SphereMetric() = default;
    SphereMetric(int width, int height, double radius_m);
    // A view of a grid (a cube, a latitude-longitude grid or a subgrid of either);
    // a full latitude-longitude grid also fills the row API. Throws
    // std::invalid_argument on a null grid.
    explicit SphereMetric(std::shared_ptr<const mesh::CellGrid> grid);

    std::size_t size() const
    {
        return grid ? grid->cellcount() : static_cast<std::size_t>(width) * static_cast<std::size_t>(height);
    }
    // Row API: latitude-longitude and gridless metrics only.
    std::size_t index(int x, int y) const
    {
        return static_cast<std::size_t>(y) * static_cast<std::size_t>(width) + static_cast<std::size_t>(x);
    }
    int column(std::size_t cell) const { return static_cast<int>(cell % static_cast<std::size_t>(width)); }
    int row(std::size_t cell) const { return static_cast<int>(cell / static_cast<std::size_t>(width)); }

    double area_m2(std::size_t cell) const
    {
        return grid ? grid->area_m2(static_cast<mesh::CellId>(cell)) : rowarea_m2[static_cast<std::size_t>(row(cell))];
    }
    double totalarea_m2() const;          // sum of the cell areas: 4 pi R^2 to relative 1e-9 on a sphere

    // The neighbour in a direction 1..8, the cell itself where none exists (beyond
    // a pole, a cube corner's missing diagonal; callers skip neighbour(cell, d) ==
    // cell); longitude wraps; on a torus (wrapsrows) rows wrap as well.
    std::size_t neighbour(std::size_t cell, int direction) const;
    // The same neighbour as a cell id, -1 where none exists: what the kernels use.
    std::int32_t neighbourid(std::size_t cell, int direction) const;
    // The direction by which neighbourid(cell, direction) reaches the cell, 0 when
    // none: oppositedirection on latitude-longitude and the torus, the grid's
    // reverse table on a cube (its frames rotate across face edges).
    int reverse(std::size_t cell, int direction) const;
    // Distance between the centres of a cell and its neighbour in a direction: on a
    // cube the great-circle distance; on latitude-longitude dy for north and south,
    // the row's dx for east and west, hypot(mean dx of the two rows, dy) for a
    // diagonal (the neighbour's row, wrapped on a torus). Undefined when the
    // neighbour does not exist.
    double distance_m(std::size_t cell, int direction) const;
    // Length of the face shared with the neighbour in a cardinal direction (a flux
    // crosses it): on a cube the arc of the shared edge (0 for a diagonal); on
    // latitude-longitude dy for east and west, the mean dx of the two rows for north
    // and south (and, kept from the row formulas, for a diagonal).
    double facelength_m(std::size_t cell, int direction) const;
    // The offset to the neighbour in a direction, east and north on the cell's
    // tangent plane: on latitude-longitude (ox dx, -oy dy) with the row's dx, the
    // diagonal with the mean dx of the two rows (its length is distance_m); on a
    // cube and on a subgrid distance_m along bearing_rad. Undefined when the
    // neighbour does not exist.
    TangentOffset tangentoffset(std::size_t cell, int direction) const;
};

// A planar torus of width x height square cells of side cell_m, wrapping in both
// directions (see the header comment). Throws std::invalid_argument unless the
// grid and the cell size are positive. Grids narrower than 3 cells in either
// direction make opposite neighbours coincide and are not supported.
SphereMetric torusmetric(int width, int height, double cell_m);
}
