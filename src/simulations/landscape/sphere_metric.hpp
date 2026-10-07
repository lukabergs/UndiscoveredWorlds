#pragma once

// Cell-centred W x W/2 latitude-longitude grid on a sphere of radius R
// (docs/architecture/erosion-hydrology-plan.md, section 3, Kernels). Row j spans
// latitudes phi_n = 90 - 180 j / H to phi_s = phi_n - 180 / H; its cell area is
// R^2 dlambda (sin phi_n - sin phi_s) with dlambda = 2 pi / W; its east-west
// spacing is R dlambda cos(phi) at the row centre; the north-south spacing is
// pi R / H everywhere; a diagonal is hypot(mean dx of the two rows, dy).
// Longitude wraps; a neighbour beyond a pole does not exist and is the cell
// itself. A torus metric (torusmetric, erosion plan step D2a: the tectonics grid)
// is a planar W x H grid of square cells that wraps rows too (no pole): every cell
// has area cell^2, dx = dy = cell, latitude 0 and radius 0, so a diagonal is
// cell sqrt(2) and totalarea_m2 is W H cell^2. Planet-free: no world class, no
// tuning header.

#include <cstddef>
#include <vector>

namespace landscape
{
// Direction codes are those of hydrosphere.river_direction: 1 north, 2 north-east,
// 3 east, 4 south-east, 5 south, 6 south-west, 7 west, 8 north-west; 0 is none
// (an outlet, or a cell that is its own receiver).
constexpr int directioncount = 8;
int directionoffsetx(int direction);   // -1, 0 or 1 (east positive)
int directionoffsety(int direction);   // -1, 0 or 1 (south positive; row 0 is the north pole)
int oppositedirection(int direction);  // 1 <-> 5, 2 <-> 6, 3 <-> 7, 4 <-> 8; 0 -> 0

struct SphereMetric
{
    int width = 0;                        // W columns (longitude)
    int height = 0;                       // H = W / 2 rows (latitude)
    double radius_m = 0.0;
    double dy_m = 0.0;                    // north-south spacing, pi R / H
    std::vector<double> rowarea_m2;       // cell area per row
    std::vector<double> rowdx_m;          // east-west spacing at the row centre, per row
    std::vector<double> rowlatitude_rad;  // row-centre latitude, per row (positive north)
    bool wrapsrows = false;               // a torus: row -1 is row H - 1 and row H is row 0

    SphereMetric() = default;
    SphereMetric(int width, int height, double radius_m);

    std::size_t size() const { return static_cast<std::size_t>(width) * static_cast<std::size_t>(height); }
    std::size_t index(int x, int y) const
    {
        return static_cast<std::size_t>(y) * static_cast<std::size_t>(width) + static_cast<std::size_t>(x);
    }
    int column(std::size_t cell) const { return static_cast<int>(cell % static_cast<std::size_t>(width)); }
    int row(std::size_t cell) const { return static_cast<int>(cell / static_cast<std::size_t>(width)); }
    double area_m2(std::size_t cell) const { return rowarea_m2[static_cast<std::size_t>(row(cell))]; }
    double totalarea_m2() const;          // sum of the cell areas: 4 pi R^2 to relative 1e-9 on a sphere

    // The neighbour in a direction 1..8: longitude wraps; a neighbour beyond a
    // pole is the cell itself (callers skip neighbour(cell, d) == cell); on a
    // torus (wrapsrows) rows wrap as well.
    std::size_t neighbour(std::size_t cell, int direction) const;
    // Distance between the centres of a cell and its neighbour in a direction:
    // dy for north and south, the row's dx for east and west, hypot(mean dx of
    // the two rows, dy) for a diagonal (the neighbour's row, wrapped on a torus).
    // Undefined when the neighbour is the cell itself.
    double distance_m(std::size_t cell, int direction) const;
    // Length of the face shared with the neighbour in a cardinal direction (a flux
    // crosses it): dy for east and west, the mean dx of the two rows for north and south.
    double facelength_m(std::size_t cell, int direction) const;
};

// A planar torus of width x height square cells of side cell_m, wrapping in both
// directions (see the header comment). Throws std::invalid_argument unless the
// grid and the cell size are positive. Grids narrower than 3 cells in either
// direction make opposite neighbours coincide and are not supported.
SphereMetric torusmetric(int width, int height, double cell_m);
}
