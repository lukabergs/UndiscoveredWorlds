#pragma once

// Geometric helpers on a CellGrid (docs/architecture/sphere-mesh.md, section 3):
// great-circle distances, cells within a radius, cells along an arc, a gnomonic
// tangent frame for templates in km, and positions drawn by area.

#include "core/mesh/cell_grid.hpp"

#include <cstdint>
#include <stdexcept>
#include <vector>

namespace mesh
{
double greatcircle_rad(const Vec3& p, const Vec3& q);                  // atan2(|p x q|, p . q)
double greatcircle_m(const Vec3& p, const Vec3& q, double radius_m);  // radius_m * greatcircle_rad

// The cells whose centres lie within radius_m (great-circle) of a unit vector,
// found by a flood over the neighbour table from the located cell; ascending id.
// Empty when no centre is close enough (or the point is outside a subgrid).
std::vector<CellId> cellswithin(const CellGrid& grid, const Vec3& centre, double radius_m);

// The cells along the great-circle arc from p to q (unit vectors, not antipodal):
// the arc is stepped at a quarter of the smallest centre spacing of the cell at
// the current point, each point located, consecutive duplicates (and, on a
// subgrid, points outside it) dropped. Starts with p's cell, ends with q's.
std::vector<CellId> cellsalongarc(const CellGrid& grid, const Vec3& p, const Vec3& q);

struct TangentFrame
{
    Vec3 origin, east, north;
    double radius_m = 0.0;
};

// The local frame at a unit vector (east and north as CellGrid::east, north).
TangentFrame tangentframe(const Vec3& centre, double radius_m);
// The unit vector at (east_km, north_km) in the frame's tangent plane (gnomonic).
Vec3 fromtangent(const TangentFrame& frame, double east_km, double north_km);

// The cell at a point drawn uniformly on the sphere from two uniform numbers in [0, 1):
// z = 2 u - 1, longitude = 2 pi v - pi.
Vec3 uniformpoint(double u, double v);
// A cell chosen with probability proportional to its area: two uniform draws in
// [0, 1) from the stream (a callable returning double), the point located. On a
// subgrid it draws again until the point lands inside (throws after 1e6 tries).
template <typename Stream>
CellId randomcellbyarea(const CellGrid& grid, Stream&& stream)
{
    for (int attempt = 0; attempt < 1000000; attempt++)
    {
        const double u = stream();
        const double v = stream();
        const CellId c = grid.locate(uniformpoint(u, v));
        if (c != nocell)
            return c;
    }
    throw std::runtime_error("randomcellbyarea: no draw landed in the subgrid");
}
// The uniform draw in [0, 1) of a 64-bit value (its top 53 bits), for streams
// built from deterministicmix64 (core/deterministic_random.hpp).
inline double uniformfrombits(std::uint64_t bits) { return static_cast<double>(bits >> 11) * (1.0 / 9007199254740992.0); }
}
