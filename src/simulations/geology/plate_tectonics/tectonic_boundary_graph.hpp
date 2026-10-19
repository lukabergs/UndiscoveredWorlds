#ifndef TECTONIC_BOUNDARY_GRAPH_HPP
#define TECTONIC_BOUNDARY_GRAPH_HPP

#include "simulations/geology/plate_tectonics/tectonic_contract.hpp"

#include <cstdint>
#include <vector>

namespace tectonics::grid { class TectonicGrid; }
namespace tectonics::boundary_graph {

struct Inputs {
    const uint32_t* plate_id = nullptr;
    const uint8_t* convergence_score = nullptr;
    const uint8_t* divergence_score = nullptr;
    const uint8_t* shear_score = nullptr;
    const uint8_t* raw_geologic_regime = nullptr;
    const uint8_t* crust_class = nullptr;
    const uint32_t* previous_boundary_segment_id = nullptr;
    const uint32_t* previous_nearest_boundary_id = nullptr;
    const contract::BoundarySegment* previous_segments = nullptr;
    uint32_t previous_segment_count = 0;
    uint32_t next_boundary_id = contract::kNoBoundaryId + 1U;
    double delta_time_myr = 0.0;
    // The engine's cell table (required: empty outputs without it). Segment centroids are
    // longitude and latitude in degrees (boundary_cube::encodecentroid).
    const grid::TectonicGrid* grid = nullptr;
};

struct Outputs {
    std::vector<contract::BoundarySegment> segments;
    std::vector<uint32_t> boundary_segment_id;
    std::vector<uint8_t> boundary_type;
    std::vector<uint8_t> boundary_regime;
    uint32_t next_boundary_id = contract::kNoBoundaryId + 1U;
};

Outputs build(const Inputs& inputs);

} // namespace tectonics::boundary_graph

#endif
