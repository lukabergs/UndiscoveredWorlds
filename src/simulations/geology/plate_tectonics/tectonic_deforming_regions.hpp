#ifndef TECTONIC_DEFORMING_REGIONS_HPP
#define TECTONIC_DEFORMING_REGIONS_HPP

#include "simulations/geology/plate_tectonics/tectonic_contract.hpp"

#include <cstdint>
#include <vector>

namespace mesh { struct Vec3; }
namespace tectonics::grid { class TectonicGrid; }
namespace tectonics::deforming_regions {

struct Inputs {
    uint32_t width = 0;
    uint32_t height = 0;
    const uint32_t* plate_id = nullptr;
    const uint8_t* crust_class = nullptr;
    const uint8_t* convergence_score = nullptr;
    const uint8_t* divergence_score = nullptr;
    const uint8_t* shear_score = nullptr;
    const uint8_t* boundary_type = nullptr;
    const uint8_t* geologic_regime = nullptr;
    const uint16_t* boundary_distance = nullptr;
    const uint32_t* nearest_boundary_id = nullptr;
    const contract::BoundarySegment* boundary_segments = nullptr;
    uint32_t boundary_segment_count = 0;
    const contract::PlateKinematics* plates = nullptr;
    uint32_t plate_count = 0;
    // The engine's cell table (W2-E); null for synthetic inputs, which use width x height.
    const grid::TectonicGrid* grid = nullptr;
    // Cube only: per plate the rotation vector in radians per update (null: every plate at
    // rest). On the cube a cell's velocity is omega x p, (east, north) in cells per update, and
    // a region's centroid longitude and latitude in degrees.
    const mesh::Vec3* plate_rotation = nullptr;
};

struct Outputs {
    std::vector<contract::DeformingRegion> regions;
    std::vector<uint32_t> deforming_region_id;
    std::vector<uint8_t> deforming_region_type;
    std::vector<float> deformation_rate;
    std::vector<float> deformation_velocity_x;
    std::vector<float> deformation_velocity_y;
};

Outputs build(const Inputs& inputs);

} // namespace tectonics::deforming_regions

#endif
