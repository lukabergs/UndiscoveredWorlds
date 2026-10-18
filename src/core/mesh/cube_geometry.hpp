#pragma once

// Geometry of the equiangular cubed sphere, orientation uw-cube-1
// (docs/architecture/sphere-mesh.md, section 2). Face f has centre c, right axis
// u (i grows) and up axis v (j = 0 at the top); face coordinates a, b lie in
// [-pi/4, pi/4] and the point is normalise(c + tan a u + tan b v). With
// D = pi / (2N), cell (f, i, j) spans a in [(i - N/2) D, (i + 1 - N/2) D] and
// b in [(N/2 - j - 1) D, (N/2 - j) D]; its id is f N^2 + j N + i.

#include "core/mesh/cell_grid.hpp"

#include <cstdint>

namespace mesh
{
constexpr int cubefacecount = 6;

struct CubeFace
{
    Vec3 centre, right, up;
};

struct FacePosition
{
    int face = 0;
    double a = 0.0, b = 0.0;
};

const CubeFace& cubeface(int face);  // 0 +X, 1 +Y, 2 -X, 3 -Y, 4 +Z, 5 -Z

double cubestep(int facecells);                     // D = pi / (2N)
double cubeedgea(int facecells, int column);        // a at the west edge of column i (i = N: the east face edge)
double cubeedgeb(int facecells, int row);           // b at the top edge of row j (j = N: the bottom face edge)
double cubecentrea(int facecells, int column);      // (i + 0.5 - N/2) D
double cubecentreb(int facecells, int row);         // (N/2 - j - 0.5) D

Vec3 cubepoint(int face, double a, double b);       // the unit vector at face coordinates a, b
// The face of a unit vector (the axis of the largest |coordinate| with its sign,
// ties to the lower face number) and its face coordinates.
FacePosition cubeposition(const Vec3& p);
// Face coordinates of p on a given face (p . c > 0 assumed).
FacePosition cubepositiononface(int face, const Vec3& p);

CellId cubeid(int facecells, int face, int column, int row);
void cubeindices(int facecells, CellId c, int& face, int& column, int& row);
CellId cubelocate(int facecells, const Vec3& p);

// The exact area on the unit sphere of the face-coordinate rectangle
// [a0, a1] x [b0, b1]: F(x1, y1) - F(x0, y1) - F(x1, y0) + F(x0, y0) with
// x = tan a, y = tan b and F(x, y) = atan(x y / sqrt(1 + x^2 + y^2)).
double cubeunitarea(double a0, double a1, double b0, double b1);

// Slot offsets: column (east positive) and row (south positive) per slot 0..8.
int slotcolumnoffset(int slot);
int slotrowoffset(int slot);

// The neighbour of a cell in slot 1..8 in its face frame: within the face by
// index; across one face edge the cell containing a point just beyond the edge
// (at the in-face neighbour's coordinate along it); -1 for the diagonal that
// would cross two edges at a cube corner (and for slot 0).
std::int32_t cubeneighbour(int facecells, CellId c, int slot);
// A point in the direction of slot k from the cell's centre: the neighbour's
// centre position extended in face coordinates (at most to the face edge or
// corner), for bearings where no neighbour exists.
Vec3 cubeslotprobe(int facecells, CellId c, int slot);
}
