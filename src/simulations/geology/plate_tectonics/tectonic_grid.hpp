#ifndef TECTONIC_GRID_HPP
#define TECTONIC_GRID_HPP

// The engine's cell table on the cube (W2-E; docs/architecture/terrain-roadmap.md, "W2-E"):
// the equiangular cubed sphere over mesh::CellGrid (N = sphere_face_cells), its cardinal
// steps across face edges, the table of unique edges, rings, discs in km, Dijkstra,
// area-weighted quantiles and centroids. Since E3 the engine runs on it alone.

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

// Cardinal directions in the cell's own frame: north toward the face's v (row j - 1), east
// toward u (column i + 1), also across face edges (mesh slots 7, 3, 1, 5).
enum Direction : std::uint8_t { west = 0, east = 1, north = 2, south = 3 };
constexpr int directioncount = 4;
constexpr std::uint8_t opposite(std::uint8_t k) { return static_cast<std::uint8_t>(k ^ 1U); }

// A unique cardinal edge: each adjacent pair once, from the lower id, in the order (cell,
// east, south, west, north).
struct Edge {
    CellId a = 0, b = 0;
    std::uint8_t direction = east;  // b is a's neighbour in this direction (a's frame)
    std::uint8_t axis = 0;          // 0 across u (east, west), 1 across v (north, south)
    // The shared edge's perpendicular from a toward b: the normalised cross product of the
    // edge's two vertices, tangent at its midpoint (contact forces and fluxes act across it).
    Vec3 normal;
    // The unit tangent of the line between the two centres at their midpoint, from a toward b
    // (the two-point difference's direction; up to 30 degrees off the normal near corners).
    Vec3 centreline;
    Vec3 midpoint;                  // the unit vector of the shared edge's middle
    double length_m = 0;            // arc length of the shared edge
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
    // The cube over a mesh cube grid (throws std::invalid_argument for a lat-lon or sub grid).
    static std::shared_ptr<const TectonicGrid> cube(std::shared_ptr<const mesh::CellGrid> mesh);
    static std::shared_ptr<const TectonicGrid> cube(int facecells, double radius_km);

    std::size_t cellcount() const { return cellcount_; }
    int facecells() const { return facecells_; }      // N
    const mesh::CellGrid* mesh() const { return mesh_.get(); }
    const std::shared_ptr<const mesh::CellGrid>& sharedmesh() const { return mesh_; }
    double radius_m() const { return mesh_->radius_m(); }
    // The unit of length of thresholds counted in cells: sqrt(total area / cells).
    double cellspacing_m() const { return cellspacing_m_; }
    // The shortest of a cell's four shared edges (the displacement limit's scale).
    double smallestedge_m(CellId c) const { return smallestedge_[c]; }

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
    // The edge between c and its neighbour in direction k (an index of edges()).
    std::uint32_t edge(CellId c, std::uint8_t k) const { return edgeof_[4 * static_cast<std::size_t>(c) + k]; }

    const std::vector<Edge>& edges() const { return edges_; }
    double area_m2(CellId c) const { return area_[c]; }
    double totalarea_m2() const { return totalarea_m2_; }
    Vec3 position(CellId c) const { return mesh_->centre(c); }  // the centre's unit vector

    // rings[d]: the cells d cardinal steps from c (d = 0 .. radius), ascending id.
    std::vector<std::vector<CellId>> rings(CellId c, int radius) const;
    // The cells whose centres lie within radius_km (great circle) of c's centre, ascending id.
    std::vector<CellId> withinkm(CellId c, double radius_km) const;
    // Shortest paths over the cardinal edges' centre distances from every source (distance 0);
    // a cell with passable[c] == 0 is never entered (empty: all passable); paths stop at limit_m.
    // Ties go to the smaller source index.
    Distances dijkstra(const std::vector<CellId>& sources, const std::vector<std::uint8_t>& passable = {},
                       double limit_m = std::numeric_limits<double>::infinity()) const;
    // The area-weighted quantile of one value per cell: the sorted value at which the area
    // below reaches q of the total less the last cell's area.
    float quantile(const std::vector<float>& values, float q) const;
    // Per owner 0 .. count - 1, the area-weighted mean direction of its centres
    // (boundary_cube::meandirection: a vanishing sum takes the member cell nearest it);
    // (0, 0, 0) for an owner without cells.
    std::vector<Vec3> centroids(const std::vector<std::uint32_t>& owners, std::uint32_t count) const;
    // c's centre in the tangent plane at origin (a unit vector), gnomonic (east, north) in
    // metres; NaN on the far hemisphere.
    std::array<double, 2> tangentoffset_m(const Vec3& origin, CellId c) const;

private:
    TectonicGrid() = default;
    void buildedges();

    std::size_t cellcount_ = 0;
    int facecells_ = 0;
    double totalarea_m2_ = 0, cellspacing_m_ = 0;
    std::shared_ptr<const mesh::CellGrid> mesh_;
    std::vector<CellId> neighbours_;            // 4 per cell
    std::vector<std::uint8_t> continuing_;      // 4 per cell
    std::vector<std::uint32_t> edgeof_;         // 4 per cell
    std::vector<double> area_, smallestedge_;
    std::vector<Edge> edges_;
};

} // namespace tectonics::grid

// Helpers of the boundary products on the cube (W2-E, E2A), shared by provenance, the
// boundary graph, junctions, their coupling and the deforming regions.
namespace tectonics::boundary_cube {

// The cube's unit of length for thresholds counted in cells (TectonicGrid::cellspacing_m).
double cellspacing_m(const grid::TectonicGrid& grid);
// A centroid on the cube goes into the contract's centroid_x and centroid_y as longitude and
// latitude in degrees (the unit vector's), until contract schema 8 names them (E5).
std::array<float, 2> encodecentroid(const mesh::Vec3& unit);
mesh::Vec3 decodecentroid(float x, float y);
// omega x p at the unit vector p in cells per update (radius / cellspacing), as 3-D tangent
// components (x, y, z); a null rotation table or a plate past it is at rest. rotation: per
// plate, radians per update.
std::array<double, 3> velocity3(const grid::TectonicGrid& grid, const mesh::Vec3* rotation,
                                uint32_t plate, uint32_t plate_count, const mesh::Vec3& p);
// The same velocity as (east, north) at p.
std::array<float, 2> velocity(const grid::TectonicGrid& grid, const mesh::Vec3* rotation,
                              uint32_t plate, uint32_t plate_count, const mesh::Vec3& p);
// A 3-D tangent vector at p as (east, north).
std::array<double, 2> eastnorth(const mesh::Vec3& p, const mesh::Vec3& vector);
// The mean direction of weighted unit vectors: sum / |sum|, or, when |sum| is below 1e-6 of
// the total weight (points spread around a great circle), the member nearest the sum's
// direction (ties and an exactly zero sum: the first member); (0, 0, 0) without members.
mesh::Vec3 meandirection(const mesh::Vec3& sum, double weight, const std::vector<mesh::Vec3>& members);

} // namespace tectonics::boundary_cube
#endif
