#include "core/mesh/cube_geometry.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace mesh
{
namespace
{
constexpr double pi = 3.14159265358979323846;
constexpr double quarterpi = 0.25 * pi;
constexpr int columnoffsets[9] = { 0, 0, 1, 1, 1, 0, -1, -1, -1 };
constexpr int rowoffsets[9] = { 0, -1, -1, 0, 1, 1, 1, 0, -1 };

const CubeFace faces[cubefacecount] = {
    { { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 } },
    { { 0.0, 1.0, 0.0 }, { -1.0, 0.0, 0.0 }, { 0.0, 0.0, 1.0 } },
    { { -1.0, 0.0, 0.0 }, { 0.0, -1.0, 0.0 }, { 0.0, 0.0, 1.0 } },
    { { 0.0, -1.0, 0.0 }, { 1.0, 0.0, 0.0 }, { 0.0, 0.0, 1.0 } },
    { { 0.0, 0.0, 1.0 }, { 0.0, 1.0, 0.0 }, { -1.0, 0.0, 0.0 } },
    { { 0.0, 0.0, -1.0 }, { 0.0, 1.0, 0.0 }, { 1.0, 0.0, 0.0 } },
};

double cornerfunction(double x, double y) { return std::atan(x * y / std::sqrt(1.0 + x * x + y * y)); }
bool validslot(int slot) { return slot >= 1 && slot <= slotcount; }

// cell_grid.cpp's helpers, expression for expression (makecubegrid's metric).
double greatcircleangle(const Vec3& p, const Vec3& q) { return std::atan2(length(cross(p, q)), dot(p, q)); }
double wrapbearing(double bearing)
{
    if (bearing < 0.0)
        bearing += 2.0 * pi;
    return bearing >= 2.0 * pi ? bearing - 2.0 * pi : bearing;
}

// Gauss-Legendre, four nodes on [-1, 1].
constexpr double gaussnodes[4] = { -0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526 };
constexpr double gaussweights[4] = { 0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538 };
constexpr std::int64_t exactareafacecells = 32;  // the corner formula below this N_L

struct EdgeLinks
{
    CubeEdgeLink link[cubefacecount][4];
};

// The links read once from cubeneighbour at N 4: the cell at index 1 along each
// edge (never a corner cell) and the cell its edge slot reaches.
EdgeLinks makeedgelinks()
{
    constexpr int n = 4;
    constexpr int edgeslots[4] = { 1, 3, 5, 7 };
    EdgeLinks links;
    for (int f = 0; f < cubefacecount; f++)
        for (int e = 0; e < 4; e++)
        {
            const int column = e == 1 ? n - 1 : (e == 3 ? 0 : 1);
            const int row = e == 0 ? 0 : (e == 2 ? n - 1 : 1);
            const std::int32_t other = cubeneighbour(n, cubeid(n, f, column, row), edgeslots[e]);
            int g = 0, i = 0, j = 0;
            cubeindices(n, static_cast<CellId>(other), g, i, j);
            CubeEdgeLink& link = links.link[f][e];
            link.face = g;
            int along = 0;
            if (j == 0 || j == n - 1)
            {
                link.edge = j == 0 ? 0 : 2;
                along = i;
            }
            else
            {
                link.edge = i == n - 1 ? 1 : 3;
                along = j;
            }
            link.reversed = along != 1;
        }
    return links;
}

const EdgeLinks& edgelinks()
{
    static const EdgeLinks links = makeedgelinks();
    return links;
}
}

const CubeFace& cubeface(int face) { return faces[face]; }

// The int64 arguments convert to double as the former int ones did, so every
// value is unchanged.
double cubestep(std::int64_t facecells) { return pi / (2.0 * static_cast<double>(facecells)); }
double cubeedgea(std::int64_t facecells, std::int64_t column)
{
    return (static_cast<double>(column) - 0.5 * static_cast<double>(facecells)) * cubestep(facecells);
}
double cubeedgeb(std::int64_t facecells, std::int64_t row)
{
    return (0.5 * static_cast<double>(facecells) - static_cast<double>(row)) * cubestep(facecells);
}
double cubecentrea(std::int64_t facecells, std::int64_t column)
{
    return (static_cast<double>(column) + 0.5 - 0.5 * static_cast<double>(facecells)) * cubestep(facecells);
}
double cubecentreb(std::int64_t facecells, std::int64_t row)
{
    return (0.5 * static_cast<double>(facecells) - static_cast<double>(row) - 0.5) * cubestep(facecells);
}

Vec3 cubepoint(int face, double a, double b)
{
    const CubeFace& f = faces[face];
    return normalise(f.centre + std::tan(a) * f.right + std::tan(b) * f.up);
}

FacePosition cubepositiononface(int face, const Vec3& p)
{
    const CubeFace& f = faces[face];
    const double c = dot(p, f.centre);
    return { face, std::atan(dot(p, f.right) / c), std::atan(dot(p, f.up) / c) };
}

FacePosition cubeposition(const Vec3& p)
{
    int face = 0;
    double best = dot(p, faces[0].centre);
    for (int f = 1; f < cubefacecount; f++)
    {
        const double value = dot(p, faces[f].centre);
        if (value > best)
        {
            best = value;
            face = f;
        }
    }
    return cubepositiononface(face, p);
}

CellId cubeid(int facecells, int face, int column, int row)
{
    const CellId n = static_cast<CellId>(facecells);
    return static_cast<CellId>(face) * n * n + static_cast<CellId>(row) * n + static_cast<CellId>(column);
}

void cubeindices(int facecells, CellId c, int& face, int& column, int& row)
{
    const CellId n = static_cast<CellId>(facecells);
    face = static_cast<int>(c / (n * n));
    const CellId rest = c % (n * n);
    row = static_cast<int>(rest / n);
    column = static_cast<int>(rest % n);
}

CellId cubelocate(int facecells, const Vec3& p)
{
    const FacePosition position = cubeposition(p);
    const double step = cubestep(facecells);
    const double half = 0.5 * facecells;
    const int column = std::clamp(static_cast<int>(std::floor(position.a / step + half)), 0, facecells - 1);
    const int row = std::clamp(static_cast<int>(std::floor(half - position.b / step)), 0, facecells - 1);
    return cubeid(facecells, position.face, column, row);
}

double cubeunitarea(double a0, double a1, double b0, double b1)
{
    const double x0 = std::tan(a0), x1 = std::tan(a1), y0 = std::tan(b0), y1 = std::tan(b1);
    return cornerfunction(x1, y1) - cornerfunction(x0, y1) - cornerfunction(x1, y0) + cornerfunction(x0, y0);
}

int slotcolumnoffset(int slot) { return validslot(slot) ? columnoffsets[slot] : 0; }
int slotrowoffset(int slot) { return validslot(slot) ? rowoffsets[slot] : 0; }

std::int32_t cubeneighbour(int facecells, CellId c, int slot)
{
    if (!validslot(slot))
        return -1;
    int face = 0, column = 0, row = 0;
    cubeindices(facecells, c, face, column, row);
    const int dx = columnoffsets[slot], dy = rowoffsets[slot];
    const int targetcolumn = column + dx, targetrow = row + dy;
    const bool outcolumn = targetcolumn < 0 || targetcolumn >= facecells;
    const bool outrow = targetrow < 0 || targetrow >= facecells;
    if (outcolumn && outrow)
        return -1;
    if (!outcolumn && !outrow)
        return static_cast<std::int32_t>(cubeid(facecells, face, targetcolumn, targetrow));
    // Just beyond the crossed edge, at the in-face coordinate along it: the edge
    // is subdivided alike on both faces, so the point lies well inside the cell
    // across the edge.
    const double epsilon = 1.0e-3 * cubestep(facecells);
    const double a = outcolumn ? (dx > 0 ? quarterpi + epsilon : -quarterpi - epsilon) : cubecentrea(facecells, targetcolumn);
    const double b = outrow ? (dy < 0 ? quarterpi + epsilon : -quarterpi - epsilon) : cubecentreb(facecells, targetrow);
    return static_cast<std::int32_t>(cubelocate(facecells, cubepoint(face, a, b)));
}

Vec3 cubeslotprobe(int facecells, CellId c, int slot)
{
    int face = 0, column = 0, row = 0;
    cubeindices(facecells, c, face, column, row);
    const double step = cubestep(facecells);
    const double a = std::clamp(cubecentrea(facecells, column) + slotcolumnoffset(slot) * step, -quarterpi, quarterpi);
    const double b = std::clamp(cubecentreb(facecells, row) - slotrowoffset(slot) * step, -quarterpi, quarterpi);
    return cubepoint(face, a, b);
}

const CubeEdgeLink& cubeedgelink(int face, int edge) { return edgelinks().link[face][edge]; }

LevelKey cubekey(std::int64_t facecells, int face, std::int64_t column, std::int64_t row)
{
    const LevelKey n = static_cast<LevelKey>(facecells);
    return static_cast<LevelKey>(face) * n * n + static_cast<LevelKey>(row) * n + static_cast<LevelKey>(column);
}

CubeCell cubecell(std::int64_t facecells, LevelKey key)
{
    const LevelKey n = static_cast<LevelKey>(facecells);
    const LevelKey rest = key % (n * n);
    CubeCell cell;
    cell.face = static_cast<int>(key / (n * n));
    cell.row = static_cast<std::int64_t>(rest / n);
    cell.column = static_cast<std::int64_t>(rest % n);
    return cell;
}

LevelKey cubelocatekey(std::int64_t facecells, const Vec3& p)
{
    // cubelocate's expressions in 64 bits.
    const FacePosition position = cubeposition(p);
    const double step = cubestep(facecells);
    const double half = 0.5 * static_cast<double>(facecells);
    const std::int64_t column =
        std::clamp(static_cast<std::int64_t>(std::floor(position.a / step + half)), std::int64_t{ 0 }, facecells - 1);
    const std::int64_t row =
        std::clamp(static_cast<std::int64_t>(std::floor(half - position.b / step)), std::int64_t{ 0 }, facecells - 1);
    return cubekey(facecells, position.face, column, row);
}

Vec3 cubekeycentre(std::int64_t facecells, LevelKey key)
{
    const CubeCell cell = cubecell(facecells, key);
    return cubepoint(cell.face, cubecentrea(facecells, cell.column), cubecentreb(facecells, cell.row));
}

double cubekeyunitarea(std::int64_t facecells, LevelKey key)
{
    const CubeCell cell = cubecell(facecells, key);
    const std::int64_t i = cell.column, j = cell.row;
    if (facecells < exactareafacecells)
        return cubeunitarea(cubeedgea(facecells, i), cubeedgea(facecells, i + 1), cubeedgeb(facecells, j + 1),
                            cubeedgeb(facecells, j));
    const double half = 0.5 * cubestep(facecells);
    const double a = cubecentrea(facecells, i), b = cubecentreb(facecells, j);
    double x2[4], y2[4];
    for (int m = 0; m < 4; m++)
    {
        const double x = std::tan(a + half * gaussnodes[m]), y = std::tan(b + half * gaussnodes[m]);
        x2[m] = x * x;
        y2[m] = y * y;
    }
    double sum = 0.0;
    for (int m = 0; m < 4; m++)
        for (int l = 0; l < 4; l++)
        {
            const double s = 1.0 + x2[m] + y2[l];
            sum += gaussweights[m] * gaussweights[l] * (1.0 + x2[m]) * (1.0 + y2[l]) / (s * std::sqrt(s));
        }
    return sum * half * half;
}

LevelKey cubeacrossedge(std::int64_t facecells, int face, int edge, std::int64_t along, std::int64_t depth)
{
    const CubeEdgeLink& link = cubeedgelink(face, edge);
    const std::int64_t last = facecells - 1;
    const std::int64_t t = link.reversed ? last - along : along;
    switch (link.edge)
    {
    case 0:
        return cubekey(facecells, link.face, t, depth);
    case 1:
        return cubekey(facecells, link.face, last - depth, t);
    case 2:
        return cubekey(facecells, link.face, t, last - depth);
    default:
        return cubekey(facecells, link.face, depth, t);
    }
}

LevelKey cubekeyneighbour(std::int64_t facecells, LevelKey key, int slot)
{
    if (!validslot(slot))
        return nokey;
    const CubeCell cell = cubecell(facecells, key);
    const std::int64_t column = cell.column + columnoffsets[slot], row = cell.row + rowoffsets[slot];
    const bool outcolumn = column < 0 || column >= facecells;
    const bool outrow = row < 0 || row >= facecells;
    if (outcolumn && outrow)
        return nokey;
    if (!outcolumn && !outrow)
        return cubekey(facecells, cell.face, column, row);
    if (outcolumn)
        return cubeacrossedge(facecells, cell.face, column < 0 ? 3 : 1, row, 0);
    return cubeacrossedge(facecells, cell.face, row < 0 ? 0 : 2, column, 0);
}

int cubekeyreverse(std::int64_t facecells, LevelKey key, int slot)
{
    const LevelKey other = cubekeyneighbour(facecells, key, slot);
    if (other == nokey)
        return 0;
    for (int m = 1; m <= slotcount; m++)
        if (cubekeyneighbour(facecells, other, m) == key)
            return m;
    return 0;
}

double cubekeydistance_m(std::int64_t facecells, LevelKey key, int slot, double radius_m)
{
    const LevelKey other = cubekeyneighbour(facecells, key, slot);
    if (other == nokey)
        return 0.0;
    return radius_m * greatcircleangle(cubekeycentre(facecells, key), cubekeycentre(facecells, other));
}

double cubekeyfacelength_m(std::int64_t facecells, LevelKey key, int slot, double radius_m)
{
    if (!validslot(slot) || slot % 2 == 0 || cubekeyneighbour(facecells, key, slot) == nokey)
        return 0.0;
    const CubeCell cell = cubecell(facecells, key);
    const std::int64_t i = cell.column, j = cell.row;
    // The shared edge: constant b for north and south, constant a for east and west.
    Vec3 from, to;
    if (slot == 1 || slot == 5)
    {
        const double b = cubeedgeb(facecells, slot == 1 ? j : j + 1);
        from = cubepoint(cell.face, cubeedgea(facecells, i), b);
        to = cubepoint(cell.face, cubeedgea(facecells, i + 1), b);
    }
    else
    {
        const double a = cubeedgea(facecells, slot == 3 ? i + 1 : i);
        from = cubepoint(cell.face, a, cubeedgeb(facecells, j + 1));
        to = cubepoint(cell.face, a, cubeedgeb(facecells, j));
    }
    return radius_m * greatcircleangle(from, to);
}

double cubekeybearing_rad(std::int64_t facecells, LevelKey key, int slot)
{
    if (!validslot(slot))
        return 0.0;
    const CubeCell cell = cubecell(facecells, key);
    const Vec3 p = cubekeycentre(facecells, key);
    const double latitude = latitudeof(p), longitude = longitudeof(p);
    const Vec3 e = eastof(latitude, longitude), nn = northof(latitude, longitude);
    const LevelKey other = cubekeyneighbour(facecells, key, slot);
    Vec3 q;
    if (other != nokey)
        q = cubekeycentre(facecells, other);
    else
    {
        // cubeslotprobe's expressions in 64 bits.
        const double step = cubestep(facecells);
        const double a = std::clamp(cubecentrea(facecells, cell.column) + slotcolumnoffset(slot) * step, -quarterpi, quarterpi);
        const double b = std::clamp(cubecentreb(facecells, cell.row) - slotrowoffset(slot) * step, -quarterpi, quarterpi);
        q = cubepoint(cell.face, a, b);
    }
    return wrapbearing(std::atan2(dot(q, e), dot(q, nn)));
}

double cubeminimumspacing_m(std::int64_t facecells, double radius_m)
{
    // The metric along a face edge is smallest at the edge's middle (1 / sqrt 2 of the
    // face centre's), so the closest centres are the edge cells' around the middle row or
    // column, along the edge; every face and both sides, as makecubegrid's scan meets them.
    const std::int64_t n = facecells;
    const std::int64_t low = std::max<std::int64_t>(0, (n - 1) / 2 - 1), high = std::min<std::int64_t>(n - 1, n / 2 + 1);
    double minimum = std::numeric_limits<double>::infinity();
    for (int face = 0; face < cubefacecount; face++)
        for (std::int64_t side : { std::int64_t{ 0 }, n - 1 })
            for (std::int64_t t = low; t <= high; t++)
                for (int slot = 1; slot <= slotcount; slot += 2)
                {
                    minimum = std::min(minimum, cubekeydistance_m(n, cubekey(n, face, side, t), slot, radius_m));
                    minimum = std::min(minimum, cubekeydistance_m(n, cubekey(n, face, t, side), slot, radius_m));
                }
    return minimum;
}

LevelKey cubeparentkey(std::int64_t facecells, int factor, LevelKey key)
{
    if (factor < 1 || facecells % factor != 0)
        throw std::invalid_argument("cubeparentkey: the level's face cells must be a multiple of the factor");
    const CubeCell cell = cubecell(facecells, key);
    return cubekey(facecells / factor, cell.face, cell.column / factor, cell.row / factor);
}

std::vector<LevelKey> cubechildkeys(std::int64_t facecells, int factor, LevelKey key)
{
    if (factor < 1 || facecells > maxlevelfacecells / factor)
        throw std::invalid_argument("cubechildkeys: the children's level exceeds 2^30 face cells");
    const CubeCell cell = cubecell(facecells, key);
    const std::int64_t fine = facecells * factor;
    std::vector<LevelKey> children;
    children.reserve(static_cast<std::size_t>(factor) * static_cast<std::size_t>(factor));
    for (int dj = 0; dj < factor; dj++)
        for (int di = 0; di < factor; di++)
            children.push_back(cubekey(fine, cell.face, cell.column * factor + di, cell.row * factor + dj));
    return children;
}
}
