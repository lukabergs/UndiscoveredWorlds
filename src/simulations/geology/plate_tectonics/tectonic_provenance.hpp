#ifndef TECTONIC_PROVENANCE_HPP
#define TECTONIC_PROVENANCE_HPP

#include "simulations/geology/plate_tectonics/tectonic_contract.hpp"

#include <array>
#include <cstdint>

namespace mesh { struct Vec3; }
namespace tectonics::grid { class TectonicGrid; }
namespace tectonics::provenance {

struct Inputs {
    uint32_t width = 0;
    uint32_t height = 0;
    const float* heightmap = nullptr;
    const uint32_t* plate_id = nullptr;
    const contract::PlateKinematics* plates = nullptr;
    uint32_t plate_count = 0;
    // Optional for synthetic inputs; the engine supplies current material-based classes.
    const uint8_t* crust_class = nullptr;
    // Zero preserves dimensionless cells/update thresholds. Positive converts native
    // velocities to cm/year, so physical boundary classes survive dt/grid changes.
    double velocity_to_cm_yr = 0.0;
    // The engine's cell table (W2-E); null for synthetic inputs, which use width x height.
    const grid::TectonicGrid* grid = nullptr;
    // Cube only: per plate the rotation vector in radians per update (null: every plate at
    // rest). The planar fields of `plates` are not read on the cube.
    const mesh::Vec3* plate_rotation = nullptr;
};

void compute_maps(const Inputs& inputs, uint8_t* convergence_score,
                  uint8_t* divergence_score, uint8_t* shear_score,
                  uint8_t* geologic_regime);

} // namespace tectonics::provenance

// The boundary products on the cube (W2-E, E2A), shared by provenance, the boundary graph,
// junctions, their coupling and the deforming regions.
namespace tectonics::boundary_cube {

// The cube's unit of length for thresholds counted in cells: sqrt(total area / cells).
double cellspacing_m(const grid::TectonicGrid& grid);
// A centroid on the cube goes into the contract's centroid_x and centroid_y as longitude and
// latitude in degrees (the unit vector's), until contract schema 8 names them (E5).
std::array<float, 2> encodecentroid(const mesh::Vec3& unit);
mesh::Vec3 decodecentroid(float x, float y);
// omega x p at the unit vector p in cells per update (radius / cellspacing), as 3-D tangent
// components (x, y, z); a null rotation table or a plate past it is at rest.
std::array<double, 3> velocity3(const grid::TectonicGrid& grid, const mesh::Vec3* rotation,
                                uint32_t plate, uint32_t plate_count, const mesh::Vec3& p);
// The same velocity as (east, north) at p.
std::array<float, 2> velocity(const grid::TectonicGrid& grid, const mesh::Vec3* rotation,
                              uint32_t plate, uint32_t plate_count, const mesh::Vec3& p);
// A 3-D tangent vector at p as (east, north).
std::array<double, 2> eastnorth(const mesh::Vec3& p, const mesh::Vec3& vector);

} // namespace tectonics::boundary_cube

#endif
