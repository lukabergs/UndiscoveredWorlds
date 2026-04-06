#ifndef TECTONIC_DEFORMING_REGIONS_HPP
#define TECTONIC_DEFORMING_REGIONS_HPP

#include "tectonic_contract.hpp"

#include <cstdint>
#include <vector>

namespace platec::deforming_regions {

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

} // namespace platec::deforming_regions

#endif
