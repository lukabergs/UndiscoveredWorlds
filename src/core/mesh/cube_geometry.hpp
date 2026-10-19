#pragma once

// Geometry of the equiangular cubed sphere, orientation uw-cube-1
// (docs/architecture/sphere-mesh.md, section 2). Face f has centre c, right axis
// u (i grows) and up axis v (j = 0 at the top); face coordinates a, b lie in
// [-pi/4, pi/4] and the point is normalise(c + tan a u + tan b v). With
// D = pi / (2N), cell (f, i, j) spans a in [(i - N/2) D, (i + 1 - N/2) D] and
// b in [(N/2 - j - 1) D, (N/2 - j) D]; its id is f N^2 + j N + i.

#include "core/mesh/cell_grid.hpp"

#include <cstdint>
#include <vector>

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

double cubestep(std::int64_t facecells);                         // D = pi / (2N)
double cubeedgea(std::int64_t facecells, std::int64_t column);    // a at the west edge of column i (i = N: the east face edge)
double cubeedgeb(std::int64_t facecells, std::int64_t row);       // b at the top edge of row j (j = N: the bottom face edge)
double cubecentrea(std::int64_t facecells, std::int64_t column);  // (i + 0.5 - N/2) D
double cubecentreb(std::int64_t facecells, std::int64_t row);     // (N/2 - j - 0.5) D

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

// ---- 64-bit level keys (sphere-mesh.md section 8): the same geometry for N_L up
// to 2^30, every function in closed form from the key. At N_L <= 16384 the
// centres, neighbours, reverses, distances, face lengths and bearings equal
// makecubegrid's bit for bit (the same expressions; makecubegrid's areas are cubekeyunitarea's).

constexpr std::int64_t maxlevelfacecells = std::int64_t{ 1 } << 30;

// The four edges of a face: 0 north (row -1 beyond, along the columns), 1 east
// (column N beyond, along the rows), 2 south (row N beyond), 3 west (column -1
// beyond). The face across an edge, the edge it arrives on, and whether the
// index along the edge runs backwards there.
struct CubeEdgeLink
{
    int face = 0;
    int edge = 0;
    bool reversed = false;
};
const CubeEdgeLink& cubeedgelink(int face, int edge);

struct CubeCell
{
    int face = 0;
    std::int64_t column = 0, row = 0;
};

LevelKey cubekey(std::int64_t facecells, int face, std::int64_t column, std::int64_t row);
CubeCell cubecell(std::int64_t facecells, LevelKey key);
LevelKey cubelocatekey(std::int64_t facecells, const Vec3& p);
Vec3 cubekeycentre(std::int64_t facecells, LevelKey key);
// The cell's area on the unit sphere: the exact corner formula (cubeunitarea)
// below N_L 32, where it is the more accurate; above, a 4 x 4 Gauss-Legendre
// rule in (a, b) on sec^2 a sec^2 b / (1 + tan^2 a + tan^2 b)^1.5, whose relative
// error stays near 1e-15 at every level (the corner formula's grows as N_L^2:
// 2e-12 at 128, 5e-10 at 4096).
double cubekeyunitarea(std::int64_t facecells, LevelKey key);
// The cell `depth` cells beyond a face edge (0: the first cell across it) at
// index `along` along the edge, on the face across: the chart of the net rotation.
LevelKey cubeacrossedge(std::int64_t facecells, int face, int edge, std::int64_t along, std::int64_t depth);
// Neighbour in slot 1..8 (nokey: slot 0, or the diagonal at a cube corner) and
// the slot by which it reaches the cell (0 when none).
LevelKey cubekeyneighbour(std::int64_t facecells, LevelKey key, int slot);
int cubekeyreverse(std::int64_t facecells, LevelKey key, int slot);
// makecubegrid's per-slot metric: the great-circle distance between the
// centres, the arc of the shared edge (cardinal slots), the initial bearing to
// the neighbour's centre or the slot's probe. Distance and face length are 0
// where no neighbour exists.
double cubekeydistance_m(std::int64_t facecells, LevelKey key, int slot, double radius_m);
double cubekeyfacelength_m(std::int64_t facecells, LevelKey key, int slot, double radius_m);
double cubekeybearing_rad(std::int64_t facecells, LevelKey key, int slot);
// The level's smallest cardinal centre spacing (makecubegrid's minimumspacing_m): the
// cells along a face edge next to its middle, along the edge.
double cubeminimumspacing_m(std::int64_t facecells, double radius_m);

// The per-face quadtree of the refinement levels: the parent at N_L / factor of
// a key at N_L (N_L a multiple of factor), and the factor^2 children at
// N_L factor of a key at N_L, ascending.
LevelKey cubeparentkey(std::int64_t facecells, int factor, LevelKey key);
std::vector<LevelKey> cubechildkeys(std::int64_t facecells, int factor, LevelKey key);
}
