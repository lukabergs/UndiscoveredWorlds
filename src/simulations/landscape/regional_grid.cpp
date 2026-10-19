#include "simulations/landscape/regional_grid.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "simulations/landscape/sphere_metric.hpp"

// The window functions are erosion plan step E2a (docs/architecture/erosion-hydrology-plan.md,
// Phase E); the patch metric below is part of the contracts (E0).

namespace landscape
{
namespace
{
constexpr double pi = 3.14159265358979323846;

// A box edge that lies within this many global cells of a cell edge is taken to be
// on it, so that a box given by cell-edge degrees (rounded in decimal) is not
// widened by a whole cell; a box edge inside a cell snaps outwards to the cell's edge.
constexpr double edgetolerance = 1e-9;

std::int64_t wrap(std::int64_t value, std::int64_t period) { return ((value % period) + period) % period; }
}

int RegionalWindow::keptfinecolumnoffset() const
{
    return static_cast<int>(wrap(static_cast<std::int64_t>(column0) - domaincolumn0, worldwidth)) * scale;
}

int RegionalWindow::keptfinerowoffset() const { return (row0 - domainrow0) * scale; }

std::int64_t RegionalWindow::absolutefinecolumn(int x) const
{
    const std::int64_t finewidth = static_cast<std::int64_t>(worldwidth) * scale;
    return wrap(static_cast<std::int64_t>(domaincolumn0) * scale + x, finewidth);
}

std::int64_t RegionalWindow::absolutefinerow(int y) const { return static_cast<std::int64_t>(domainrow0) * scale + y; }

double RegionalWindow::finelatitude_rad(int y) const
{
    // The row-centre latitude of SphereMetric for row J of the kW x kH grid (the
    // same expressions, so that it equals PatchMetric::rowlatitude_rad bit for bit).
    const int fineheight = worldheight * scale;
    const int j = static_cast<int>(absolutefinerow(y));
    const double north = pi / 2.0 - pi * j / fineheight;
    const double south = pi / 2.0 - pi * (j + 1) / fineheight;
    return 0.5 * (north + south);
}

double RegionalWindow::finelongitude_rad(int x) const
{
    // 2 pi (I + 0.5) / (k W): relief_roughness.cpp's expression for a global cell.
    const int finewidth = worldwidth * scale;
    return 2.0 * pi * (static_cast<double>(absolutefinecolumn(x)) + 0.5) / finewidth;
}

std::size_t RegionalWindow::globalcell(int x, int y) const
{
    const std::size_t column = static_cast<std::size_t>(absolutefinecolumn(x) / scale);
    const std::size_t row = static_cast<std::size_t>(absolutefinerow(y) / scale);
    return row * static_cast<std::size_t>(worldwidth) + column;
}

RegionalWindow regionalwindow(int worldwidth, int scale, int halo, double west_deg, double south_deg, double east_deg,
                              double north_deg)
{
    if (worldwidth < 4 || worldwidth % 2 != 0)
        throw std::invalid_argument("regionalwindow: the world width must be even and at least 4");
    if (scale < 1)
        throw std::invalid_argument("regionalwindow: the scale must be at least 1");
    if (halo < 0)
        throw std::invalid_argument("regionalwindow: the halo must not be negative");
    if (!std::isfinite(west_deg) || !std::isfinite(east_deg) || !std::isfinite(south_deg) || !std::isfinite(north_deg))
        throw std::invalid_argument("regionalwindow: the box must be finite");
    if (south_deg < -90.0 || north_deg > 90.0 || south_deg >= north_deg)
        throw std::invalid_argument("regionalwindow: the latitudes must satisfy -90 <= south < north <= 90");
    // |lon| <= 360 keeps the cell indices below far from any integer overflow.
    if (std::fabs(west_deg) > 360.0 || std::fabs(east_deg) > 360.0)
        throw std::invalid_argument("regionalwindow: a longitude must lie within -360 .. 360");
    // The width eastwards from west to east; 0 is a single meridian (west == east, or
    // west 180 and east -180), not the whole circle (that is -180 .. 180).
    const double span = east_deg - west_deg + (west_deg > east_deg ? 360.0 : 0.0);
    if (span == 0.0)
        throw std::invalid_argument("regionalwindow: west and east are the same meridian");

    RegionalWindow window;
    window.worldwidth = worldwidth;
    window.worldheight = worldwidth / 2;
    window.scale = scale;
    window.halo = halo;
    const int width = worldwidth, height = window.worldheight;

    // Column i spans longitudes -180 + 360 i / W .. -180 + 360 (i + 1) / W; the box
    // snaps outwards to whole columns; west > east crosses the date line, so the
    // east edge lies one turn further east.
    const double westcells = (west_deg + 180.0) / 360.0 * width;
    double eastcells = (east_deg + 180.0) / 360.0 * width;
    if (west_deg > east_deg)
        eastcells += width;
    const std::int64_t first = static_cast<std::int64_t>(std::floor(westcells + edgetolerance));
    std::int64_t last = static_cast<std::int64_t>(std::ceil(eastcells - edgetolerance));  // one past the box
    if (last <= first)
        last = first + 1;
    if (last - first + 2 * static_cast<std::int64_t>(halo) > width)
        throw std::invalid_argument("regionalwindow: the domain (box and halo) is wider than the planet");

    // Row j spans latitudes 90 - 180 j / H .. 90 - 180 (j + 1) / H.
    const double northcells = (90.0 - north_deg) / 180.0 * height;
    const double southcells = (90.0 - south_deg) / 180.0 * height;
    int row0 = static_cast<int>(std::floor(northcells + edgetolerance));
    int rowend = static_cast<int>(std::ceil(southcells - edgetolerance));
    row0 = std::max(0, std::min(height - 1, row0));
    rowend = std::max(row0 + 1, std::min(height, rowend));

    window.column0 = static_cast<int>(wrap(first, width));
    window.columns = static_cast<int>(last - first);
    window.row0 = row0;
    window.rows = rowend - row0;
    window.domaincolumn0 = static_cast<int>(wrap(static_cast<std::int64_t>(window.column0) - halo, width));
    window.domaincolumns = window.columns + 2 * halo;
    window.domainrow0 = std::max(0, row0 - halo);
    window.domainrows = std::min(height, rowend + halo) - window.domainrow0;
    return window;
}

std::shared_ptr<const mesh::CellGrid> regionalwindowgrid(const RegionalWindow& window, double radius_m)
{
    mesh::LatLonWindowBox box;
    box.column0 = window.absolutefinecolumn(0);
    box.columns = window.domainfinewidth();
    box.row0 = static_cast<int>(window.absolutefinerow(0));
    box.rows = window.domainfineheight();
    return mesh::latlonwindowgrid(window.worldwidth * window.scale, radius_m, box);
}

double northsouthspacing_m(const mesh::CellGrid& grid)
{
    // The cube: the level's arc pi R / (2 N_L) (the equiangular cells' arc along the face centre lines;
    // T8 W2-W2d). Lat-lon: the north slot's distance, the row spacing on every cell (also beyond the edge).
    if (grid.kind() == mesh::GridKind::cube)
        return grid.facecells() > 0 ? pi * grid.radius_m() / (2.0 * static_cast<double>(grid.facecells())) : 0.0;
    return grid.cellcount() > 0 ? grid.distance_m(0, 1) : 0.0;
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
