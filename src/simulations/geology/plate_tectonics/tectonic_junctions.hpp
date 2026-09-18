#ifndef TECTONIC_JUNCTIONS_HPP
#define TECTONIC_JUNCTIONS_HPP

#include "tectonic_contract.hpp"

#include <cstdint>
#include <vector>

namespace platec::tectonic_junctions {

struct Inputs {
    uint32_t width = 0;
    uint32_t height = 0;
    const uint32_t* plate_id = nullptr;
    const uint32_t* boundary_segment_id = nullptr;
    const contract::BoundarySegment* boundary_segments = nullptr;
    uint32_t boundary_segment_count = 0;
    const contract::PlateKinematics* plates = nullptr;
    uint32_t plate_count = 0;
};

struct Outputs {
    std::vector<contract::TectonicJunction> junctions;
};

Outputs build(const Inputs& inputs);

} // namespace platec::tectonic_junctions

#endif
