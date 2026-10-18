#ifndef TECTONIC_GRID_HPP
#define TECTONIC_GRID_HPP

// The engine's own cell table (W2-E, step E1; docs/architecture/terrain-roadmap.md, "W2-E").
// One table, two builders: the cube over mesh::CellGrid (the engine's target grid) and a
// planar scaffold that reproduces today's torus exactly (ids y W + x, the stencils' order,
// one cell edge in metres). The scaffold is a transitional fixture: the default torus path
// and the tests use it until E3 switches the engine to the cube and deletes it.

#include "core/mesh/cell_grid.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <vector>

namespace tectonics::grid {

using CellId = mesh::CellId;
using mesh::Vec3;
constexpr CellId nocell = mesh::nocell;

// Cardinal directions in the cell's own frame, in the order of the planar stencils
// (x - 1, x + 1, y - 1, y + 1). On the cube, north is toward the face's v (row j - 1) and
// east toward u (column i + 1), also across face edges (mesh slots 7, 3, 1, 5).
enum Direction : std::uint8_t { west = 0, east = 1, north = 2, south = 3 };
constexpr int directioncount = 4;
constexpr std::uint8_t opposite(std::uint8_t k) { return static_cast<std::uint8_t>(k ^ 1U); }

// A unique cardinal edge. Planar: edge 2 i + 0 is (i, east of i), 2 i + 1 is (i, south of i),
// so edge ids are today's face ids 2 index + axis. Cube: each adjacent pair once, from the
// lower id, in the order (cell, east, south, west, north).
struct Edge {
    CellId a = 0, b = 0;
    std::uint8_t direction = east;  // b is a's neighbour in this direction (a's frame)
    std::uint8_t axis = 0;          // 0 across x or u (east, west), 1 across y or v (north, south)
    Vec3 normal;                    // unit normal from a toward b (planar: +x or +y)
    double length_m = 0;            // arc length of the shared face
    double distance_m = 0;          // centre to centre
};

// The cell one step ahead and the direction that continues straight in its own frame
// (it differs from the step's direction only across a cube face edge).
struct Step {
    CellId cell = nocell;
    std::uint8_t direction = east;
};

struct Distances {
    std::vector<double> distance_m;  // infinity where unreached
    std::vector<std::uint32_t> source;  // index into the sources, UINT32_MAX where unreached
};

class TectonicGrid {
public:
    enum class Kind { planar, cube };

    // The cube over a mesh cube grid (throws std::invalid_argument for a lat-lon or sub grid).
    static std::shared_ptr<const TectonicGrid> cube(std::shared_ptr<const mesh::CellGrid> mesh);
    static std::shared_ptr<const TectonicGrid> cube(int facecells, double radius_km);
    // Today's torus: width x height cells, each width_km / width wide (square cells).
    static std::shared_ptr<const TectonicGrid> planar(std::uint32_t width, std::uint32_t height,
                                                      double width_km);

    Kind kind() const { return kind_; }
    std::size_t cellcount() const { return cellcount_; }
    std::uint32_t width() const { return width_; }    // planar columns (0 on the cube)
    std::uint32_t height() const { return height_; }  // planar rows (0 on the cube)
    int facecells() const { return facecells_; }      // cube N (0 on the scaffold)
    const mesh::CellGrid* mesh() const { return mesh_.get(); }  // the cube's mesh (null on the scaffold)
    double cell_m() const { return cell_m_; }         // the scaffold's cell edge (0 on the cube)

    CellId neighbour(CellId c, std::uint8_t k) const { return neighbours_[4 * static_cast<std::size_t>(c) + k]; }
    // west, east, north, south
    std::array<CellId, 4> neighbours(CellId c) const {
        const std::size_t base = 4 * static_cast<std::size_t>(c);
        return {neighbours_[base], neighbours_[base + 1], neighbours_[base + 2], neighbours_[base + 3]};
    }
    Step ahead(CellId c, std::uint8_t k) const {
        const std::size_t slot = 4 * static_cast<std::size_t>(c) + k;
        return {neighbours_[slot], continuing_[slot]};
    }

    const std::vector<Edge>& edges() const { return edges_; }
    double area_m2(CellId c) const { return area_.empty() ? cell_area_m2_ : area_[c]; }
    double totalarea_m2() const { return totalarea_m2_; }
    // Planar: (x, y, 0), the cell's column and row (the mass-centre convention of the torus).
    // Cube: the centre's unit vector.
    Vec3 position(CellId c) const;

    // rings[d]: the cells d cardinal steps from c (d = 0 .. radius), ascending id.
    std::vector<std::vector<CellId>> rings(CellId c, int radius) const;
    // The cells whose centres lie within radius_km of c's centre (shortest torus offset on the
    // scaffold, great circle on the cube), ascending id.
    std::vector<CellId> withinkm(CellId c, double radius_km) const;
    // Shortest paths over the cardinal edges' centre distances from every source (distance 0);
    // a cell with passable[c] == 0 is never entered (empty: all passable); paths stop at limit_m.
    // Ties go to the smaller source index.
    Distances dijkstra(const std::vector<CellId>& sources, const std::vector<std::uint8_t>& passable = {},
                       double limit_m = std::numeric_limits<double>::infinity()) const;
    // The area-weighted quantile of one value per cell: the sorted value at which the area
    // below reaches q of the total less the last cell's area. On the scaffold (equal areas) it
    // is today's quantile_threshold: sorted[floor(q (n - 1))] in float.
    float quantile(const std::vector<float>& values, float q) const;
    // Per owner 0 .. count - 1, the area-weighted centroid of its cells. Planar: the circular
    // mean of columns and rows, in [0, width) and [0, height), z 0. Cube: the normalised mean of
    // the centres. An owner without cells gets (0, 0, 0).
    std::vector<Vec3> centroids(const std::vector<std::uint32_t>& owners, std::uint32_t count) const;
    // c's centre in the tangent plane at origin (a position as returned above), metres.
    // Planar: the shortest torus offset (x, y) times the cell edge. Cube: gnomonic (east, north);
    // NaN on the far hemisphere.
    std::array<double, 2> tangentoffset_m(const Vec3& origin, CellId c) const;

private:
    TectonicGrid() = default;
    void buildedges();

    Kind kind_ = Kind::planar;
    std::size_t cellcount_ = 0;
    std::uint32_t width_ = 0, height_ = 0;
    int facecells_ = 0;
    double cell_m_ = 0, cell_area_m2_ = 0, totalarea_m2_ = 0;
    std::shared_ptr<const mesh::CellGrid> mesh_;
    std::vector<CellId> neighbours_;            // 4 per cell
    std::vector<std::uint8_t> continuing_;      // 4 per cell
    std::vector<double> area_;                  // cube only
    std::vector<Edge> edges_;
};

} // namespace tectonics::grid
#endif
