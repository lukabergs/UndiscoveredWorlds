#ifndef TECTONIC_JUNCTIONS_HPP
#define TECTONIC_JUNCTIONS_HPP

#include "simulations/geology/plate_tectonics/tectonic_contract.hpp"

#include <cstdint>
#include <vector>

namespace tectonics::grid { class TectonicGrid; }
namespace tectonics::tectonic_junctions {

struct Inputs {
    uint32_t width = 0;
    uint32_t height = 0;
    const uint32_t* plate_id = nullptr;
    const uint32_t* boundary_segment_id = nullptr;
    const contract::BoundarySegment* boundary_segments = nullptr;
    uint32_t boundary_segment_count = 0;
    const contract::PlateKinematics* plates = nullptr;
    uint32_t plate_count = 0;
    // The engine's cell table (W2-E); null for synthetic inputs, which use width x height.
    const grid::TectonicGrid* grid = nullptr;
};

struct Outputs {
    std::vector<contract::TectonicJunction> junctions;
};

Outputs build(const Inputs& inputs);

} // namespace tectonics::tectonic_junctions

#endif
