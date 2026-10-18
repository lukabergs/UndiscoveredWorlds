#include "core/mesh/cube_geometry.hpp"

#include <algorithm>
#include <cmath>

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
}

const CubeFace& cubeface(int face) { return faces[face]; }

double cubestep(int facecells) { return pi / (2.0 * facecells); }
double cubeedgea(int facecells, int column) { return (column - 0.5 * facecells) * cubestep(facecells); }
double cubeedgeb(int facecells, int row) { return (0.5 * facecells - row) * cubestep(facecells); }
double cubecentrea(int facecells, int column) { return (column + 0.5 - 0.5 * facecells) * cubestep(facecells); }
double cubecentreb(int facecells, int row) { return (0.5 * facecells - row - 0.5) * cubestep(facecells); }

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
}
