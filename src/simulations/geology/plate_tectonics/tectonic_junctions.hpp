#ifndef TECTONIC_JUNCTIONS_HPP
#define TECTONIC_JUNCTIONS_HPP

#include "simulations/geology/plate_tectonics/tectonic_contract.hpp"

#include <array>
#include <cstdint>
#include <vector>

namespace mesh { struct Vec3; }
namespace tectonics::grid { class TectonicGrid; }
namespace tectonics::tectonic_junctions {

struct Inputs {
    const uint32_t* plate_id = nullptr;
    const uint32_t* boundary_segment_id = nullptr;
    const contract::BoundarySegment* boundary_segments = nullptr;
    uint32_t boundary_segment_count = 0;
    uint32_t plate_count = 0;
    // The engine's cell table (required: no junctions without it).
    const grid::TectonicGrid* grid = nullptr;
    // Per plate the rotation vector in radians per update (null: every plate at rest).
    // Centroids are longitude and latitude in degrees, strikes and velocities (east, north) at
    // the junction, velocities in cells per update.
    const mesh::Vec3* plate_rotation = nullptr;
};

struct Outputs {
    std::vector<contract::TectonicJunction> junctions;
};

Outputs build(const Inputs& inputs);

// The cells around one mesh vertex of the cube: four, or three at a cube corner (the fourth
// 0xffffffff, mesh::nocell).
struct Vertex {
    std::array<uint32_t, 4> cells{0xffffffffu, 0xffffffffu, 0xffffffffu, 0xffffffffu};
    uint32_t count = 0U;
};
// The cube's 6 N^2 + 2 mesh vertices, each once (the candidate points of build): cell c's
// corner is listed by the smallest id among the cells that share it.
std::vector<Vertex> cube_vertices(const grid::TectonicGrid& grid);

} // namespace tectonics::tectonic_junctions

#endif
