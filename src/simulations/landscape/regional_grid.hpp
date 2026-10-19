#pragma once

// Regional windows on the planet-wide fine grid (docs/architecture/erosion-hydrology-plan.md,
// Phase E, the contracts above its table and section 6, "Decided for Phase E").
// The world's W x H cell-centred grid, refined by a whole-number scale k, is a
// planet-wide fine grid of kW x kH cells; fine cell (I, J) is the same cell in every
// window, and global cell (i, j) holds fine columns k i .. k (i + 1) - 1 and rows
// k j .. k (j + 1) - 1. A window is a box of whole global cells (the kept part)
// widened by a halo of global cells on every side (the domain, computed and then
// discarded); longitude wraps, rows stop at the poles. Degrees follow the planet:
// global column i spans longitudes -180 + 360 i / W to -180 + 360 (i + 1) / W, row j
// spans latitudes 90 - 180 j / H to 90 - 180 (j + 1) / H. Planet-free.

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

#include "core/mesh/cell_grid.hpp"

namespace landscape
{
struct RegionalWindow
{
    int worldwidth = 0;   // W global columns
    int worldheight = 0;  // H = W / 2 global rows
    int scale = 0;        // k fine cells per global cell edge, at least 1
    int halo = 0;         // global cells computed around the kept box and discarded

    // The kept box in global cells: columns column0 .. column0 + columns - 1, each
    // taken modulo W (the box may cross the date line), rows row0 .. row0 + rows - 1.
    int column0 = 0;
    int columns = 0;
    int row0 = 0;
    int rows = 0;

    // The domain in global cells: the kept box widened by `halo` on every side,
    // rows clipped to 0 .. H - 1 (domainrow0 = max(0, row0 - halo)); columns are
    // never clipped (domaincolumn0 = column0 - halo modulo W).
    int domaincolumn0 = 0;
    int domaincolumns = 0;
    int domainrow0 = 0;
    int domainrows = 0;

    // Fine sizes: the kept grid is columns k x rows k, the domain grid
    // domaincolumns k x domainrows k; the kept grid starts at fine column
    // (column0 - domaincolumn0, modulo W) k and fine row (row0 - domainrow0) k of
    // the domain grid.
    int finewidth() const { return columns * scale; }
    int fineheight() const { return rows * scale; }
    int domainfinewidth() const { return domaincolumns * scale; }
    int domainfineheight() const { return domainrows * scale; }
    int keptfinecolumnoffset() const;
    int keptfinerowoffset() const;

    // Absolute fine indices on the kW x kH grid of a domain fine cell (x, y),
    // 0 <= x < domainfinewidth(), 0 <= y < domainfineheight(): the column wraps.
    std::int64_t absolutefinecolumn(int x) const;
    std::int64_t absolutefinerow(int y) const;
    // The fine cell centre in radians, latitude positive north, longitude measured
    // eastwards from column 0's west edge (2 pi (I + 0.5) / (k W), the convention of
    // relief_roughness.cpp for global cells).
    double finelatitude_rad(int y) const;
    double finelongitude_rad(int x) const;
    // The level-0 ancestor of domain fine cell (x, y): the id on the world's grid
    // (RegionalGlobalInputs::grid; on lat-lon J W + I, the global column I wrapped to
    // 0 .. W - 1) of the global cell holding it.
    std::size_t globalcell(int x, int y) const;
};

// The window that covers a box in degrees, snapped outwards to whole global cells.
// The box runs eastwards from west to east: west > east crosses the date line (east
// lies one turn further), so its width is east - west, plus 360 when west > east;
// -180 .. 180 is the whole circle. A box of width 0 (west == east, or west 180 and
// east -180: the same meridian) is rejected, as are a longitude with |lon| > 360, a
// non-finite value, south >= north, a latitude outside -90 .. 90, scale < 1,
// halo < 0, worldwidth < 4 or an odd worldwidth, and a domain wider than the planet
// (columns + 2 halo > W): all throw std::invalid_argument.
RegionalWindow regionalwindow(int worldwidth, int scale, int halo, double west_deg, double south_deg, double east_deg,
                              double north_deg);

// The domain's fine grid as a table-driven cell grid (sphere-mesh.md section 8, "Window
// grids"): on lat-lon mesh::latlonwindowgrid of the kW x kH planet-wide fine grid over the
// domain box, local id y domainfinewidth() + x (the domain's row-major order), nothing
// wrapping inside it (no neighbour beyond its edge: -1), latitudes, spacings, areas and
// distances those of finepatchmetric bit for bit (MeshWindow.LatLonWindowEqualsPatchMetric),
// globalkey(c) the absolute fine index I + J kW. Throws as latlonwindowgrid.
std::shared_ptr<const mesh::CellGrid> regionalwindowgrid(const RegionalWindow& window, double radius_m);

// The north-south spacing of a window grid's cells, which the units take as the level's
// cell size (on lat-lon pi R / (k H), finepatchmetric's dy_m: the north slot's distance).
double northsouthspacing_m(const mesh::CellGrid& grid);

// The metric of a domain's fine grid on a sphere of radius R, the reference that
// regionalwindowgrid reproduces (the pass reads the window grid; T8 W2-W2a): a patch, not the whole
// sphere, so nothing wraps (its edges are the domain's boundary). Row y's centre
// latitude, east-west spacing R dlambda cos(phi) with dlambda = 2 pi / (k W), the
// north-south spacing pi R / (k H) and the cell area as SphereMetric computes them
// for the planet-wide fine grid, so a cell's metric is the same in every window.
struct PatchMetric
{
    int width = 0;
    int height = 0;
    double radius_m = 0.0;
    double dy_m = 0.0;
    std::vector<double> rowdx_m;
    std::vector<double> rowarea_m2;
    std::vector<double> rowlatitude_rad;

    std::size_t size() const { return static_cast<std::size_t>(width) * static_cast<std::size_t>(height); }
    std::size_t index(int x, int y) const
    {
        return static_cast<std::size_t>(y) * static_cast<std::size_t>(width) + static_cast<std::size_t>(x);
    }
    // The neighbour in a direction 1..8 (sphere_metric.hpp's codes) or the cell
    // itself beyond the patch's edge (callers skip it); distance as SphereMetric.
    std::size_t neighbour(std::size_t cell, int direction) const;
    double distance_m(std::size_t cell, int direction) const;
};

PatchMetric finepatchmetric(const RegionalWindow& window, double radius_m);
}
