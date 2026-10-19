#include "core/mesh/cell_grid.hpp"

#include "core/mesh/cube_geometry.hpp"

#include <algorithm>
#include <limits>
#include <stdexcept>

namespace mesh
{
namespace
{
constexpr double pi = 3.14159265358979323846;
constexpr int latloncolumnoffsets[9] = { 0, 0, 1, 1, 1, 0, -1, -1, -1 };
constexpr int latlonrowoffsets[9] = { 0, -1, -1, 0, 1, 1, 1, 0, -1 };

bool validslot(int k) { return k >= 1 && k <= slotcount; }
std::size_t slotindex(CellId c, int k) { return static_cast<std::size_t>(c) * slotcount + static_cast<std::size_t>(k - 1); }
int oppositeslot(int k) { return (k + 3) % slotcount + 1; }
double greatcircleangle(const Vec3& p, const Vec3& q) { return std::atan2(length(cross(p, q)), dot(p, q)); }

double wrapbearing(double bearing)
{
    if (bearing < 0.0)
        bearing += 2.0 * pi;
    return bearing >= 2.0 * pi ? bearing - 2.0 * pi : bearing;
}
}

Vec3 unitvector(double latitude_rad, double longitude_rad)
{
    const double c = std::cos(latitude_rad);
    return { c * std::cos(longitude_rad), c * std::sin(longitude_rad), std::sin(latitude_rad) };
}

double latitudeof(const Vec3& p) { return std::atan2(p.z, std::hypot(p.x, p.y)); }
double longitudeof(const Vec3& p) { return std::atan2(p.y, p.x); }

Vec3 eastof(double, double longitude_rad) { return { -std::sin(longitude_rad), std::cos(longitude_rad), 0.0 }; }

Vec3 northof(double latitude_rad, double longitude_rad)
{
    const double s = std::sin(latitude_rad);
    return { -s * std::cos(longitude_rad), -s * std::sin(longitude_rad), std::cos(latitude_rad) };
}

Vec3 CellGrid::centre(CellId c) const
{
    if (stored_)
        return centre_[c];
    const int x = static_cast<int>(c % static_cast<CellId>(analyticwidth()));
    const int y = static_cast<int>(c / static_cast<CellId>(analyticwidth()));
    return unitvector(rowlatitude_[static_cast<std::size_t>(y)], columnlongitude_[static_cast<std::size_t>(x)]);
}

double CellGrid::latitude_rad(CellId c) const
{
    if (compact_)
        return latitudeof(centre_[c]);
    if (stored_)
        return latitude_[c];
    return rowlatitude_[c / static_cast<CellId>(analyticwidth())];
}

double CellGrid::longitude_rad(CellId c) const
{
    if (compact_)
        return longitudeof(centre_[c]);
    if (stored_)
        return longitude_[c];
    return columnlongitude_[c % static_cast<CellId>(analyticwidth())];
}

double CellGrid::area_m2(CellId c) const
{
    if (stored_)
        return area_[c];
    return rowarea_[c / static_cast<CellId>(analyticwidth())];
}

std::int32_t CellGrid::neighbour(CellId c, int k) const
{
    if (!validslot(k))
        return -1;
    if (stored_)
        return neighbour_[slotindex(c, k)];
    if (latlonwindow_)
    {
        // A box: nothing wraps (finepatchmetric's patch).
        const int y = static_cast<int>(c / static_cast<CellId>(windowcolumns_)) + latlonrowoffsets[k];
        const int x = static_cast<int>(c % static_cast<CellId>(windowcolumns_)) + latloncolumnoffsets[k];
        if (y < 0 || y >= windowrows_ || x < 0 || x >= windowcolumns_)
            return -1;
        return static_cast<std::int32_t>(y * windowcolumns_ + x);
    }
    const int y = static_cast<int>(c / static_cast<CellId>(columns_)) + latlonrowoffsets[k];
    if (y < 0 || y >= rows_)
        return -1;
    const int x = (static_cast<int>(c % static_cast<CellId>(columns_)) + latloncolumnoffsets[k] + columns_) % columns_;
    return static_cast<std::int32_t>(y * columns_ + x);
}

int CellGrid::reverse(CellId c, int k) const
{
    if (!validslot(k))
        return 0;
    if (stored_)
        return reverse_[slotindex(c, k)];
    return neighbour(c, k) < 0 ? 0 : oppositeslot(k);
}

double CellGrid::distance_m(CellId c, int k) const
{
    if (!validslot(k))
        return 0.0;
    if (compact_)
        return cubekeydistance_m(facecells_, keys_[c], k, radius_m_);
    if (stored_)
        return distance_[slotindex(c, k)];
    // landscape::SphereMetric::distance_m, expression for expression.
    const int ox = latloncolumnoffsets[k], oy = latlonrowoffsets[k];
    if (ox == 0)
        return dy_m_;
    const int y = static_cast<int>(c / static_cast<CellId>(analyticwidth()));
    if (oy == 0)
        return rowdx_[static_cast<std::size_t>(y)];
    const std::int32_t n = neighbour(c, k);
    const int other = n < 0 ? y : n / analyticwidth();
    return std::hypot(0.5 * (rowdx_[static_cast<std::size_t>(y)] + rowdx_[static_cast<std::size_t>(other)]), dy_m_);
}

double CellGrid::facelength_m(CellId c, int k) const
{
    if (!validslot(k))
        return 0.0;
    if (compact_)
        return cubekeyfacelength_m(facecells_, keys_[c], k, radius_m_);
    if (stored_)
        return facelength_[slotindex(c, k)];
    // landscape::SphereMetric::facelength_m, expression for expression.
    const int oy = latlonrowoffsets[k];
    if (oy == 0)
        return dy_m_;
    const int y = static_cast<int>(c / static_cast<CellId>(analyticwidth()));
    const std::int32_t n = neighbour(c, k);
    const int other = n < 0 ? y : n / analyticwidth();
    return 0.5 * (rowdx_[static_cast<std::size_t>(y)] + rowdx_[static_cast<std::size_t>(other)]);
}

double CellGrid::bearing_rad(CellId c, int k) const
{
    if (!validslot(k))
        return 0.0;
    if (compact_)
        return cubekeybearing_rad(facecells_, keys_[c], k);
    if (stored_)
        return bearing_[slotindex(c, k)];
    const int ox = latloncolumnoffsets[k], oy = latlonrowoffsets[k];
    if (ox == 0)
        return oy < 0 ? 0.0 : pi;
    if (oy == 0)
        return ox > 0 ? 0.5 * pi : 1.5 * pi;
    const int y = static_cast<int>(c / static_cast<CellId>(analyticwidth()));
    const std::int32_t n = neighbour(c, k);
    const int other = n < 0 ? y : n / analyticwidth();
    const double theta = std::atan2(0.5 * (rowdx_[static_cast<std::size_t>(y)] + rowdx_[static_cast<std::size_t>(other)]), dy_m_);
    if (ox > 0)
        return oy < 0 ? theta : pi - theta;
    return oy > 0 ? pi + theta : 2.0 * pi - theta;
}

LevelKey CellGrid::locatekey(const Vec3& p) const
{
    if (kind_ == GridKind::cube)
        return cubelocatekey(facecells_, p);
    // asin(z) for the row (cheaper than atan2 with hypot; p is a unit vector).
    const double latitude = std::asin(std::clamp(p.z, -1.0, 1.0)), longitude = longitudeof(p);
    const int y = std::clamp(static_cast<int>(std::floor((0.5 * pi - latitude) / (pi / rows_))), 0, rows_ - 1);
    int x = static_cast<int>(std::floor((longitude + pi) / (2.0 * pi / columns_)));
    x = std::clamp(x, 0, columns_);
    if (x == columns_)
        x = 0;
    return static_cast<LevelKey>(y) * static_cast<LevelKey>(columns_) + static_cast<LevelKey>(x);
}

CellId CellGrid::locate(const Vec3& p) const
{
    const LevelKey key = locatekey(p);
    return subgrid_ ? localofkey(key) : static_cast<CellId>(key);
}

CellId CellGrid::globalid(CellId local) const
{
    if (!subgrid_)
        return local;
    if (!keywindow_ && !latlonwindow_)
        return globalids_[local];
    const LevelKey key = globalkey(local);
    return key < nocell ? static_cast<CellId>(key) : nocell;
}

CellId CellGrid::localid(CellId global) const
{
    if (!subgrid_)
        return global < cellcount_ ? global : nocell;
    if (keywindow_ || latlonwindow_)
        return localofkey(global);
    const auto found = std::lower_bound(globalids_.begin(), globalids_.end(), global);
    if (found == globalids_.end() || *found != global)
        return nocell;
    return static_cast<CellId>(found - globalids_.begin());
}

LevelKey CellGrid::globalkey(CellId local) const
{
    if (!subgrid_)
        return local;
    if (keywindow_)
        return keys_[local];
    if (latlonwindow_)
    {
        const std::int64_t width = columns_;
        const std::int64_t x = local % static_cast<CellId>(windowcolumns_), y = local / static_cast<CellId>(windowcolumns_);
        const std::int64_t column = (windowcolumn0_ + x) % width;
        return static_cast<LevelKey>(windowrow0_ + y) * static_cast<LevelKey>(width) + static_cast<LevelKey>(column);
    }
    return globalids_[local];
}

CellId CellGrid::localofkey(LevelKey key) const
{
    if (!subgrid_)
        return key < cellcount_ ? static_cast<CellId>(key) : nocell;
    if (keywindow_)
    {
        const auto found = std::lower_bound(keys_.begin(), keys_.end(), key);
        if (found == keys_.end() || *found != key)
            return nocell;
        return static_cast<CellId>(found - keys_.begin());
    }
    if (latlonwindow_)
    {
        const LevelKey width = static_cast<LevelKey>(columns_);
        if (key >= width * static_cast<LevelKey>(rows_))
            return nocell;
        const std::int64_t y = static_cast<std::int64_t>(key / width) - windowrow0_;
        const std::int64_t x = (static_cast<std::int64_t>(key % width) - windowcolumn0_ + columns_) % columns_;
        if (y < 0 || y >= windowrows_ || x >= windowcolumns_)
            return nocell;
        return static_cast<CellId>(y * windowcolumns_ + x);
    }
    return key < nocell ? localid(static_cast<CellId>(key)) : nocell;
}

CellId CellGrid::cellof(int x, int y) const
{
    if (kind_ != GridKind::latlon || x < 0 || x >= columns_ || y < 0 || y >= rows_)
        return nocell;
    const LevelKey key = static_cast<LevelKey>(y) * static_cast<LevelKey>(columns_) + static_cast<LevelKey>(x);
    return subgrid_ ? localofkey(key) : static_cast<CellId>(key);
}

int CellGrid::column(CellId c) const
{
    return kind_ == GridKind::latlon ? static_cast<int>(globalkey(c) % static_cast<LevelKey>(columns_)) : -1;
}

int CellGrid::row(CellId c) const
{
    return kind_ == GridKind::latlon ? static_cast<int>(globalkey(c) / static_cast<LevelKey>(columns_)) : -1;
}


std::shared_ptr<const CellGrid> makecubegrid(int facecells, double radius_m)
{
    if (facecells < 1 || facecells > 16384 || !(radius_m > 0.0))
        throw std::invalid_argument("makecubegrid needs 1 <= N <= 16384 and a positive radius");
    auto grid = std::make_shared<CellGrid>();
    CellGrid& g = *grid;
    const int n = facecells;
    g.kind_ = GridKind::cube;
    g.facecells_ = n;
    g.cellcount_ = static_cast<std::size_t>(cubefacecount) * static_cast<std::size_t>(n) * static_cast<std::size_t>(n);
    g.radius_m_ = radius_m;
    g.stored_ = true;
    const std::size_t count = g.cellcount_;
    g.centre_.resize(count);
    g.latitude_.resize(count);
    g.longitude_.resize(count);
    g.area_.resize(count);
    g.neighbour_.resize(count * slotcount);
    g.reverse_.assign(count * slotcount, 0);
    g.distance_.assign(count * slotcount, 0.0);
    g.facelength_.assign(count * slotcount, 0.0);
    g.bearing_.assign(count * slotcount, 0.0);

    const double r2 = radius_m * radius_m;
    double total = 0.0;
    for (CellId c = 0; c < count; c++)
    {
        int face = 0, i = 0, j = 0;
        cubeindices(n, c, face, i, j);
        const Vec3 p = cubepoint(face, cubecentrea(n, i), cubecentreb(n, j));
        g.centre_[c] = p;
        g.latitude_[c] = latitudeof(p);
        g.longitude_[c] = longitudeof(p);
        g.area_[c] = r2 * cubekeyunitarea(n, c);  // the one area rule of every cube level
        total += g.area_[c];
        for (int k = 1; k <= slotcount; k++)
            g.neighbour_[slotindex(c, k)] = cubeneighbour(n, c, k);
    }
    g.totalarea_m2_ = total;

    double minimum = std::numeric_limits<double>::infinity();
    for (CellId c = 0; c < count; c++)
    {
        int face = 0, i = 0, j = 0;
        cubeindices(n, c, face, i, j);
        const Vec3& p = g.centre_[c];
        const Vec3 e = eastof(g.latitude_[c], g.longitude_[c]);
        const Vec3 nn = northof(g.latitude_[c], g.longitude_[c]);
        for (int k = 1; k <= slotcount; k++)
        {
            const std::size_t s = slotindex(c, k);
            const std::int32_t other = g.neighbour_[s];
            const Vec3 q = other < 0 ? cubeslotprobe(n, c, k) : g.centre_[static_cast<std::size_t>(other)];
            g.bearing_[s] = wrapbearing(std::atan2(dot(q, e), dot(q, nn)));
            if (other < 0)
                continue;
            for (int m = 1; m <= slotcount; m++)
                if (g.neighbour_[slotindex(static_cast<CellId>(other), m)] == static_cast<std::int32_t>(c))
                {
                    g.reverse_[s] = static_cast<std::int8_t>(m);
                    break;
                }
            g.distance_[s] = radius_m * greatcircleangle(p, q);
            if (k % 2 == 1)
            {
                // The shared edge: constant b for north and south, constant a for east and west.
                Vec3 from, to;
                if (k == 1 || k == 5)
                {
                    const double b = cubeedgeb(n, k == 1 ? j : j + 1);
                    from = cubepoint(face, cubeedgea(n, i), b);
                    to = cubepoint(face, cubeedgea(n, i + 1), b);
                }
                else
                {
                    const double a = cubeedgea(n, k == 3 ? i + 1 : i);
                    from = cubepoint(face, a, cubeedgeb(n, j + 1));
                    to = cubepoint(face, a, cubeedgeb(n, j));
                }
                g.facelength_[s] = radius_m * greatcircleangle(from, to);
                minimum = std::min(minimum, g.distance_[s]);
            }
        }
    }
    g.minimumspacing_m_ = minimum;
    return grid;
}

std::shared_ptr<const CellGrid> makelatlongrid(int columns, double radius_m)
{
    if (columns < 2 || columns % 2 != 0 || columns > 65536 || !(radius_m > 0.0))
        throw std::invalid_argument("makelatlongrid needs an even W in 2..65536 and a positive radius");
    auto grid = std::make_shared<CellGrid>();
    CellGrid& g = *grid;
    const int width = columns, height = columns / 2;
    const double radius = radius_m;
    g.kind_ = GridKind::latlon;
    g.columns_ = width;
    g.rows_ = height;
    g.cellcount_ = static_cast<std::size_t>(width) * static_cast<std::size_t>(height);
    g.radius_m_ = radius_m;
    // landscape::SphereMetric's constructor, expression for expression.
    g.dy_m_ = pi * radius / height;
    const double dlambda = 2.0 * pi / width;
    g.rowarea_.resize(static_cast<std::size_t>(height));
    g.rowdx_.resize(static_cast<std::size_t>(height));
    g.rowlatitude_.resize(static_cast<std::size_t>(height));
    for (int j = 0; j < height; j++)
    {
        const double north = pi / 2.0 - pi * j / height;
        const double south = pi / 2.0 - pi * (j + 1) / height;
        const double centre = 0.5 * (north + south);
        const std::size_t row = static_cast<std::size_t>(j);
        g.rowarea_[row] = radius * radius * dlambda * (std::sin(north) - std::sin(south));
        g.rowdx_[row] = radius * dlambda * std::cos(centre);
        g.rowlatitude_[row] = centre;
    }
    g.columnlongitude_.resize(static_cast<std::size_t>(width));
    for (int i = 0; i < width; i++)
        g.columnlongitude_[static_cast<std::size_t>(i)] = -pi + (i + 0.5) * dlambda;
    double total = 0.0;
    for (double area : g.rowarea_)
        total += area * width;
    g.totalarea_m2_ = total;
    g.minimumspacing_m_ = std::min(g.dy_m_, *std::min_element(g.rowdx_.begin(), g.rowdx_.end()));
    return grid;
}

std::shared_ptr<const CellGrid> subgrid(const CellGrid& base, const std::vector<CellId>& globalids)
{
    if (base.issubgrid())
        throw std::invalid_argument("subgrid needs a full grid");
    std::vector<CellId> ids = globalids;
    std::sort(ids.begin(), ids.end());
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
    if (!ids.empty() && ids.back() >= base.cellcount())
        throw std::invalid_argument("subgrid: a cell id is out of range");
    auto grid = std::make_shared<CellGrid>();
    CellGrid& g = *grid;
    g.kind_ = base.kind_;
    g.facecells_ = base.facecells_;
    g.columns_ = base.columns_;
    g.rows_ = base.rows_;
    g.radius_m_ = base.radius_m_;
    g.minimumspacing_m_ = base.minimumspacing_m_;
    g.cellcount_ = ids.size();
    g.subgrid_ = true;
    g.stored_ = true;
    g.globalids_ = std::move(ids);
    const std::size_t count = g.cellcount_;
    g.centre_.resize(count);
    g.latitude_.resize(count);
    g.longitude_.resize(count);
    g.area_.resize(count);
    g.neighbour_.assign(count * slotcount, -1);
    g.reverse_.assign(count * slotcount, 0);
    g.distance_.resize(count * slotcount);
    g.facelength_.resize(count * slotcount);
    g.bearing_.resize(count * slotcount);
    double total = 0.0;
    for (CellId local = 0; local < count; local++)
    {
        const CellId global = g.globalids_[local];
        g.centre_[local] = base.centre(global);
        g.latitude_[local] = base.latitude_rad(global);
        g.longitude_[local] = base.longitude_rad(global);
        g.area_[local] = base.area_m2(global);
        total += g.area_[local];
        for (int k = 1; k <= slotcount; k++)
        {
            const std::size_t s = slotindex(local, k);
            const std::int32_t other = base.neighbour(global, k);
            const CellId otherlocal = other < 0 ? nocell : g.localid(static_cast<CellId>(other));
            if (otherlocal != nocell)
            {
                g.neighbour_[s] = static_cast<std::int32_t>(otherlocal);
                g.reverse_[s] = static_cast<std::int8_t>(base.reverse(global, k));
            }
            g.distance_[s] = base.distance_m(global, k);
            g.facelength_[s] = base.facelength_m(global, k);
            g.bearing_[s] = base.bearing_rad(global, k);
        }
    }
    g.totalarea_m2_ = total;
    return grid;
}

std::shared_ptr<const CellGrid> cubewindowgrid(int facecells, double radius_m, const std::vector<LevelKey>& keys, bool compact)
{
    if (facecells < 1 || facecells > maxlevelfacecells || !(radius_m > 0.0))
        throw std::invalid_argument("cubewindowgrid needs 1 <= N_L <= 2^30 and a positive radius");
    std::vector<LevelKey> sorted = keys;
    std::sort(sorted.begin(), sorted.end());
    sorted.erase(std::unique(sorted.begin(), sorted.end()), sorted.end());
    const LevelKey n = static_cast<LevelKey>(facecells);
    if (!sorted.empty() && sorted.back() >= static_cast<LevelKey>(cubefacecount) * n * n)
        throw std::invalid_argument("cubewindowgrid: a key is out of range");
    if (sorted.size() >= static_cast<std::size_t>(nocell))
        throw std::invalid_argument("cubewindowgrid: too many cells for 32-bit local ids");
    auto grid = std::make_shared<CellGrid>();
    CellGrid& g = *grid;
    g.kind_ = GridKind::cube;
    g.facecells_ = facecells;
    g.radius_m_ = radius_m;
    g.cellcount_ = sorted.size();
    g.subgrid_ = true;
    g.keywindow_ = true;
    g.stored_ = true;
    g.compact_ = compact;
    g.keys_ = std::move(sorted);
    const std::size_t count = g.cellcount_;
    g.centre_.resize(count);
    g.area_.resize(count);
    g.neighbour_.assign(count * slotcount, -1);
    g.reverse_.assign(count * slotcount, 0);
    if (!compact)
    {
        g.latitude_.resize(count);
        g.longitude_.resize(count);
        g.distance_.resize(count * slotcount);
        g.facelength_.resize(count * slotcount);
        g.bearing_.resize(count * slotcount);
    }
    const double r2 = radius_m * radius_m;
    double total = 0.0;
    for (CellId local = 0; local < count; local++)
    {
        const LevelKey key = g.keys_[local];
        const Vec3 p = cubekeycentre(facecells, key);
        g.centre_[local] = p;
        g.area_[local] = r2 * cubekeyunitarea(facecells, key);
        total += g.area_[local];
        if (!compact)
        {
            g.latitude_[local] = latitudeof(p);
            g.longitude_[local] = longitudeof(p);
        }
        for (int k = 1; k <= slotcount; k++)
        {
            const std::size_t s = slotindex(local, k);
            const LevelKey other = cubekeyneighbour(facecells, key, k);
            const CellId otherlocal = other == nokey ? nocell : g.localofkey(other);
            if (otherlocal != nocell)
            {
                g.neighbour_[s] = static_cast<std::int32_t>(otherlocal);
                g.reverse_[s] = static_cast<std::int8_t>(cubekeyreverse(facecells, key, k));
            }
            if (compact)
                continue;
            g.distance_[s] = cubekeydistance_m(facecells, key, k, radius_m);
            g.facelength_[s] = cubekeyfacelength_m(facecells, key, k, radius_m);
            g.bearing_[s] = cubekeybearing_rad(facecells, key, k);
        }
    }
    g.totalarea_m2_ = total;
    g.minimumspacing_m_ = cubeminimumspacing_m(facecells, radius_m);  // the level's own (makecubegrid's)
    return grid;
}

std::shared_ptr<const CellGrid> latlonwindowgrid(int levelcolumns, double radius_m, const LatLonWindowBox& box)
{
    if (levelcolumns < 2 || levelcolumns % 2 != 0 || !(radius_m > 0.0))
        throw std::invalid_argument("latlonwindowgrid needs an even level width of at least 2 and a positive radius");
    const int width = levelcolumns, height = levelcolumns / 2;
    if (box.columns < 1 || box.columns > width || box.rows < 1 || box.row0 < 0 || box.row0 + box.rows > height)
        throw std::invalid_argument("latlonwindowgrid: the box is empty or larger than the level grid");
    const std::size_t count = static_cast<std::size_t>(box.columns) * static_cast<std::size_t>(box.rows);
    if (count >= static_cast<std::size_t>(nocell))
        throw std::invalid_argument("latlonwindowgrid: too many cells for 32-bit local ids");
    auto grid = std::make_shared<CellGrid>();
    CellGrid& g = *grid;
    const double radius = radius_m;
    g.kind_ = GridKind::latlon;
    g.columns_ = width;
    g.rows_ = height;
    g.cellcount_ = count;
    g.radius_m_ = radius_m;
    g.subgrid_ = true;
    g.latlonwindow_ = true;
    g.windowcolumn0_ = ((box.column0 % width) + width) % width;
    g.windowrow0_ = box.row0;
    g.windowcolumns_ = box.columns;
    g.windowrows_ = box.rows;
    // makelatlongrid's (and finepatchmetric's) expressions for the level grid's rows.
    g.dy_m_ = pi * radius / height;
    const double dlambda = 2.0 * pi / width;
    g.rowarea_.resize(static_cast<std::size_t>(box.rows));
    g.rowdx_.resize(static_cast<std::size_t>(box.rows));
    g.rowlatitude_.resize(static_cast<std::size_t>(box.rows));
    for (int y = 0; y < box.rows; y++)
    {
        const int j = box.row0 + y;
        const double north = pi / 2.0 - pi * j / height;
        const double south = pi / 2.0 - pi * (j + 1) / height;
        const double centre = 0.5 * (north + south);
        const std::size_t row = static_cast<std::size_t>(y);
        g.rowarea_[row] = radius * radius * dlambda * (std::sin(north) - std::sin(south));
        g.rowdx_[row] = radius * dlambda * std::cos(centre);
        g.rowlatitude_[row] = centre;
    }
    g.columnlongitude_.resize(static_cast<std::size_t>(box.columns));
    for (int x = 0; x < box.columns; x++)
    {
        const std::int64_t i = (g.windowcolumn0_ + x) % width;
        g.columnlongitude_[static_cast<std::size_t>(x)] = -pi + (static_cast<double>(i) + 0.5) * dlambda;
    }
    double total = 0.0;
    for (double area : g.rowarea_)
        total += area * box.columns;
    g.totalarea_m2_ = total;
    // The level's own spacing (makelatlongrid's min over the rows: the polar rows').
    const auto rowdx = [&](int j) {
        const double north = pi / 2.0 - pi * j / height;
        const double south = pi / 2.0 - pi * (j + 1) / height;
        return radius * dlambda * std::cos(0.5 * (north + south));
    };
    g.minimumspacing_m_ = std::min(g.dy_m_, std::min(rowdx(0), rowdx(height - 1)));
    return grid;
}

std::shared_ptr<const CellGrid> makeplanarwindow(int columns, int rows, double cell_m)
{
    if (columns < 1 || rows < 1 || !(cell_m > 0.0) || !std::isfinite(cell_m))
        throw std::invalid_argument("makeplanarwindow needs at least one cell and a positive cell side");
    const std::size_t count = static_cast<std::size_t>(columns) * static_cast<std::size_t>(rows);
    if (count >= static_cast<std::size_t>(nocell))
        throw std::invalid_argument("makeplanarwindow: too many cells for 32-bit local ids");
    auto grid = std::make_shared<CellGrid>();
    CellGrid& g = *grid;
    g.kind_ = GridKind::latlon;
    g.columns_ = columns;
    g.rows_ = rows;
    g.cellcount_ = count;
    g.radius_m_ = 6.371e6;  // nominal: nothing reads it on a planar patch
    g.subgrid_ = true;
    g.latlonwindow_ = true;
    g.windowcolumns_ = columns;
    g.windowrows_ = rows;
    g.dy_m_ = cell_m;
    g.rowarea_.assign(static_cast<std::size_t>(rows), cell_m * cell_m);
    g.rowdx_.assign(static_cast<std::size_t>(rows), cell_m);
    g.rowlatitude_.assign(static_cast<std::size_t>(rows), 0.0);
    g.columnlongitude_.assign(static_cast<std::size_t>(columns), 0.0);
    g.totalarea_m2_ = cell_m * cell_m * static_cast<double>(count);
    g.minimumspacing_m_ = cell_m;
    return grid;
}
}
