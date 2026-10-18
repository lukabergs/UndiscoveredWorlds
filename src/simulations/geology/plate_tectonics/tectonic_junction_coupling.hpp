#ifndef TECTONIC_JUNCTION_COUPLING_HPP
#define TECTONIC_JUNCTION_COUPLING_HPP

#include "simulations/geology/plate_tectonics/tectonic_contract.hpp"

#include <cstdint>
#include <vector>

namespace tectonics::grid { class TectonicGrid; }
namespace tectonics::tectonic_junction_coupling {

struct Inputs {
    uint32_t width = 0;
    uint32_t height = 0;
    const contract::TectonicJunction* junctions = nullptr;
    uint32_t junction_count = 0;
    // The engine's cell table (W2-E); null for synthetic inputs, which use width x height.
    // On the cube junction centroids are longitude and latitude in degrees and the rift
    // tension (east, north) at each cell.
    const grid::TectonicGrid* grid = nullptr;
};

struct Outputs {
    std::vector<float> rift_influence;
    std::vector<float> transform_influence;
    std::vector<float> subduction_influence;
    std::vector<float> rift_tension_x;
    std::vector<float> rift_tension_y;
};

Outputs build(const Inputs& inputs);

} // namespace tectonics::tectonic_junction_coupling

#endif
