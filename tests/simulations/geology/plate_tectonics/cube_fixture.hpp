#pragma once

// The cube fixtures of the engine's tests (W2-E, E3: the engine runs on the cube alone):
// grids, the equatorial band of the four side faces, motions and cell lookups.
//
// The band is faces 3, 0, 1, 2 in the order of longitude (-135 to 225 degrees) as a periodic
// strip of 4N columns and N rows (row 0 at the north): each side face's column coordinate a
// is its longitude offset, so the band's columns are bounded by meridians and its east-west
// neighbours are the next column (across the face edges too). Rows are not parallels (b is
// atan(z / x) on face 0). Layouts periodic in x map onto it; the polar faces take a
// background value.

#include "core/mesh/cube_geometry.hpp"
#include "core/mesh/mesh_geometry.hpp"
#include "simulations/geology/plate_tectonics/tectonic_grid.hpp"
#include "simulations/geology/plate_tectonics/tectonic_sphere_transport.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <memory>
#include <vector>

namespace cubefixture {

using mesh::Vec3;
using tectonics::grid::CellId;
using tectonics::grid::TectonicGrid;
using tectonics::sphere_transport::Motion;

constexpr double earth_km = 6371.0;
constexpr double pi = 3.14159265358979323846;

inline std::shared_ptr<const TectonicGrid> cube(int n, double radius_km = earth_km) {
    return TectonicGrid::cube(n, radius_km);
}

// Geographic east and north at a cell.
inline Vec3 east(const TectonicGrid& grid, CellId c) {
    const Vec3 p = grid.position(c);
    return mesh::eastof(mesh::latitudeof(p), mesh::longitudeof(p));
}
inline Vec3 north(const TectonicGrid& grid, CellId c) {
    const Vec3 p = grid.position(c);
    return mesh::northof(mesh::latitudeof(p), mesh::longitudeof(p));
}
inline double longitude_deg(const TectonicGrid& grid, CellId c) { return mesh::longitudeof(grid.position(c)) * 180 / pi; }
inline double latitude_deg(const TectonicGrid& grid, CellId c) { return mesh::latitudeof(grid.position(c)) * 180 / pi; }
// The cell at a longitude and latitude in degrees.
inline CellId at(const TectonicGrid& grid, double longitude, double latitude) {
    return grid.mesh()->locate(mesh::unitvector(latitude * pi / 180, longitude * pi / 180));
}

struct Band {
    std::shared_ptr<const TectonicGrid> grid;
    int n = 0;

    explicit Band(int facecells, double radius_km = earth_km) : grid(cube(facecells, radius_km)), n(facecells) {}
    int columns() const { return 4 * n; }
    int rows() const { return n; }
    std::size_t cellcount() const { return grid->cellcount(); }
    double radius_m() const { return grid->radius_m(); }
    // The width of a column at the equator, R pi / 2N.
    double spacing_m() const { return grid->radius_m() * pi / (2.0 * n); }
    CellId cell(int x, int y) const {
        static constexpr int faces[4] = {3, 0, 1, 2};
        const int wrapped = ((x % columns()) + columns()) % columns();
        return mesh::cubeid(n, faces[wrapped / n], wrapped % n, y);
    }
    bool inband(CellId c) const {
        int face = 0, i = 0, j = 0;
        mesh::cubeindices(n, c, face, i, j);
        return face < 4;
    }
    // The band's column and row of a cell (-1 on the polar faces).
    int column(CellId c) const {
        static constexpr int order[4] = {1, 2, 3, 0};  // face 0 is the second block, face 3 the first
        int face = 0, i = 0, j = 0;
        mesh::cubeindices(n, c, face, i, j);
        return face < 4 ? order[face] * n + i : -1;
    }
    int row(CellId c) const {
        int face = 0, i = 0, j = 0;
        mesh::cubeindices(n, c, face, i, j);
        return face < 4 ? j : -1;
    }
    // A rotation about the z axis that carries the equator east by `cells` columns per update.
    Vec3 eastward(double cells) const { return {0, 0, cells * pi / (2.0 * n)}; }
};

// The radius at which a band of N cells per face edge has equatorial columns cell_m wide, km.
inline double radiusfor_km(int n, double cell_m) { return cell_m * 2.0 * n / pi / 1000.0; }

// One value per cell: f(column, row) on the band, the background on the polar faces.
template <class T, class F>
std::vector<T> layout(const Band& band, const T& background, F f) {
    std::vector<T> values(band.cellcount(), background);
    for (int y = 0; y < band.rows(); ++y)
        for (int x = 0; x < band.columns(); ++x) values[band.cell(x, y)] = f(x, y);
    return values;
}

// Two plates split by the plane y = 0 (longitude 0 and 180 degrees), a grid line on every face
// for even N: plate 0 at y < 0, plate 1 at y > 0. On the band plate 0 holds columns 3.5 N .. 4 N - 1
// and 0 .. 1.5 N - 1.
inline std::vector<uint32_t> hemispheres(const TectonicGrid& grid) {
    std::vector<uint32_t> owners(grid.cellcount());
    for (CellId c = 0; c < grid.cellcount(); ++c) owners[c] = grid.position(c).y < 0 ? 0 : 1;
    return owners;
}

// Nothing moves: every cell at rest (no plate).
inline Motion still(const TectonicGrid& grid) {
    Motion motion;
    motion.plate.assign(grid.cellcount(), tectonics::sphere_transport::noplate);
    return motion;
}
// One rigid rotation (axis times radians per update) of every cell, as plate 0.
inline Motion rotating(const TectonicGrid& grid, const Vec3& rotation) {
    Motion motion;
    motion.rotation = {rotation};
    motion.plate.assign(grid.cellcount(), 0);
    return motion;
}
// Plates by owner, each with its rotation per update.
inline Motion plates(const std::vector<uint32_t>& owners, std::vector<Vec3> rotation) {
    Motion motion;
    motion.rotation = std::move(rotation);
    motion.plate = owners;
    return motion;
}
// Offsets per cell (tangent vectors in metres per update) on plates at rest: every cell plate 0
// without owners, else its owner.
inline Motion offsets(const TectonicGrid& grid, std::vector<Vec3> offset, const std::vector<uint32_t>& owners = {}) {
    Motion motion;
    motion.plate = owners.empty() ? std::vector<uint32_t>(grid.cellcount(), 0) : owners;
    uint32_t count = 1;
    for (const uint32_t p : motion.plate) count = std::max(count, p + 1);
    motion.rotation.assign(count, Vec3{});
    motion.offset_m = std::move(offset);
    return motion;
}

} // namespace cubefixture
