#pragma once

// Cell-grid helpers of the coast, shelf, ocean-feature, arid-landform and import terrain
// functions (T8 W2-T2; docs/architecture/sphere-mesh.md sections 3, 4 and 11). They walk
// the neighbour table of any mesh::CellGrid: on the lat-lon grid a walk east or west wraps
// in longitude and a walk north or south stops at a pole, exactly as the legacy (x, y) loops
// wrapped x and clipped y; on the cube a walk continues straight across a face edge (its
// frame turns with the neighbour table's reverse slots).

#include <cstdint>
#include <type_traits>
#include <vector>

#include "core/mesh/cell_grid.hpp"
#include "core/mesh/mesh_geometry.hpp"
#include "simulations/landscape/spherical_noise.hpp"

class planet;
class boolshapetemplate;

namespace terrain::coastgrid
{
using mesh::CellId;

// The cells in the order of the legacy loops "for x { for y }": on lat-lon column by column,
// each from north to south; on the cube face by face, each face column by column (i outer, j
// inner); a subgrid in ascending id order. The order of "for y { for x }" is the id order.
std::vector<CellId> columnorder(const mesh::CellGrid& grid);

// A slot (1..8) turned clockwise by turn eighths.
inline int turned(int slot, int turn) { return ((slot - 1 + turn) % 8 + 8) % 8 + 1; }

// The neighbour of c in slot turned by turn, -1 when none. turn is updated so that the next
// step in the same slot continues straight ahead (on lat-lon it stays 0).
inline std::int32_t step(const mesh::CellGrid& grid, CellId c, int slot, int& turn)
{
    const int k = turned(slot, turn);
    const std::int32_t n = grid.neighbour(c, k);
    if (n < 0)
        return -1;
    const int ahead = (grid.reverse(c, k) + 3) % 8 + 1;  // the slot opposite the way back
    turn = ((turn + ahead - k) % 8 + 8) % 8;
    return n;
}

// The cell dx steps east (west when negative) and then dy steps south (north when negative)
// of c, the frame carried along the walk; -1 when the walk leaves the grid (beyond a lat-lon
// pole). clamped: a walk north or south stops at the last cell instead (the legacy clip of y).
std::int32_t offset(const mesh::CellGrid& grid, CellId c, int dx, int dy, bool clamped = false);

// The square box of half-width r around c in the legacy order (dx from -r to r outer, dy from
// -r to r inner): f(cell, dx, dy) for every offset that stays on the grid. A callback that
// returns bool stops the walk when it returns false.
template <typename F>
void box(const mesh::CellGrid& grid, CellId c, int r, F&& f);

// planet::outline and planet::coast on any grid: a land (sea) cell with a sea (land) cell among
// its four cardinal neighbours, false where one is missing; sea as planet::seawrap (base
// elevation in whole metres at or below sea level, no lake).
bool outline(const planet& world, CellId c);
bool coast(const planet& world, CellId c);

// The slot [x][y] of a cell in a legacy ARRAYWIDTH x ARRAYHEIGHT array (the vector<vector<>>
// arguments of functions other modules call): column and row on lat-lon; elsewhere the id folded
// into rows of ARRAYWIDTH.
void legacyslot(const mesh::CellGrid& grid, CellId c, int& x, int& y);

// disruptseacoastline on a cell: base elevations above avedepth lowered to it in the disc of
// the legacy box of half-width size / 2 around centre.
void disruptseacoastline(planet& world, CellId centre, float avedepth, bool raise, int size);

// Sphere forms of the legacy pixel geometry (tuning::terrain::sphere_noise_coast).

// The legacy pixel, 2 pi R / legacy_columns km.
double legacypixelkm(const planet& world);
// The cells whose centres lie within radius_km of p, ascending id; at least the cell containing p.
std::vector<CellId> cellswithinkm(const planet& world, const mesh::Vec3& p, double radius_km);
// Great-circle distance in km.
double distancekm(const planet& world, const mesh::Vec3& p, const mesh::Vec3& q);
// The point east_km and north_km from a cell centre in its tangent plane (gnomonic).
mesh::Vec3 movedpoint(const planet& world, CellId c, double east_km, double north_km);
// A uniform number in [0, 1) from two draws of the ambient sequence (random(0, 32767)).
double uniformdraw();
// A position drawn by area from the ambient sequence.
mesh::Vec3 randompoint();
// The chance that a cell sees at least one event of a chance per legacy pixel: 1 - (1 - p)^(A / pixel area).
double cellchance(const planet& world, CellId c, double perpixel);

// The cells under a pixel template stamped at p (pixel i east, j south of the template's
// centre, mirrored east-west or north-south), one legacy pixel apart through p's tangent frame;
// in the legacy loop order (i outer, j inner), each cell once.
std::vector<CellId> stamptemplate(const planet& world, const mesh::Vec3& p, const boolshapetemplate& shape, bool mirroreast, bool mirrorsouth);

// Ridge codes on a cell (one bit per neighbour slot, stored at both ends: at the neighbour in
// its reverse slot). deleteridge also clears the mountain height of a cell left without ridges.
bool getoceanridge(const planet& world, CellId c, int slot);
void deleteridge(planet& world, CellId c, int slot);
void deleteoceanridge(planet& world, CellId c, int slot);

// A legacy createfractal on the sphere: per cell maximum x clamp(mean + deviation z, 0, 1), at
// least 1, z the noise at the centre over its normalisation (unit deviation).
struct LegacyNoise
{
    double wavelength_km = 0.0;
    double persistence = 0.5;
    double mean = 0.5;
    double deviation = 0.2;
};
LegacyNoise grain8noise();
LegacyNoise grain16noise();
LegacyNoise grain4noise();
LegacyNoise finenoise();
// The unit-deviation noise z per cell (seed: the world seed and salt).
std::vector<double> spherenoise(const planet& world, int salt, const LegacyNoise& noise);
// The same noise as a function of any unit vector: z = sum / normalisation.
class SphereNoise
{
public:
    SphereNoise(const planet& world, int salt, const LegacyNoise& noise);
    double z(const mesh::Vec3& p) const { return field_.sum(p.x, p.y, p.z) / normalisation_; }
    int value(const mesh::Vec3& p, int maximum) const; // the legacy fractal's value there

private:
    LegacyNoise legacy_;
    landscape::FractalNoiseField field_;
    double normalisation_ = 1.0;
};
// The inverse of mesh::fromtangent: the gnomonic east and north km of p in the frame (p in the
// frame's hemisphere).
void totangent(const mesh::TangentFrame& frame, const mesh::Vec3& p, double& east_km, double& north_km);
std::vector<int> spherefractal(const planet& world, int salt, const LegacyNoise& noise, int maximum);

template <typename F>
void box(const mesh::CellGrid& grid, CellId c, int r, F&& f)
{
    std::vector<std::int32_t> columns(static_cast<std::size_t>(2 * r + 1), -1);
    std::vector<int> turns(static_cast<std::size_t>(2 * r + 1), 0);
    columns[static_cast<std::size_t>(r)] = static_cast<std::int32_t>(c);
    for (int side = -1; side <= 1; side += 2)
    {
        int turn = 0;
        std::int32_t at = static_cast<std::int32_t>(c);
        for (int n = 1; n <= r && at >= 0; n++)
        {
            at = step(grid, static_cast<CellId>(at), side < 0 ? 7 : 3, turn);
            columns[static_cast<std::size_t>(r + side * n)] = at;
            turns[static_cast<std::size_t>(r + side * n)] = turn;
        }
    }
    std::vector<std::int32_t> north(static_cast<std::size_t>(r), -1), south(static_cast<std::size_t>(r), -1);
    for (int dx = -r; dx <= r; dx++)
    {
        const std::int32_t column = columns[static_cast<std::size_t>(dx + r)];
        if (column < 0)
            continue;
        for (int side = 0; side < 2; side++)
        {
            std::vector<std::int32_t>& ray = side == 0 ? north : south;
            int turn = turns[static_cast<std::size_t>(dx + r)];
            std::int32_t at = column;
            for (int n = 0; n < r; n++)
            {
                at = at >= 0 ? step(grid, static_cast<CellId>(at), side == 0 ? 1 : 5, turn) : -1;
                ray[static_cast<std::size_t>(n)] = at;
            }
        }
        for (int dy = -r; dy <= r; dy++)
        {
            const std::int32_t cell = dy < 0 ? north[static_cast<std::size_t>(-dy - 1)]
                : dy > 0 ? south[static_cast<std::size_t>(dy - 1)] : column;
            if (cell < 0)
                continue;
            if constexpr (std::is_same_v<decltype(f(CellId{}, 0, 0)), bool>)
            {
                if (!f(static_cast<CellId>(cell), dx, dy))
                    return;
            }
            else
                f(static_cast<CellId>(cell), dx, dy);
        }
    }
}
}
