#ifndef TECTONIC_JUNCTION_COUPLING_HPP
#define TECTONIC_JUNCTION_COUPLING_HPP

#include "tectonic_contract.hpp"

#include <cstdint>
#include <vector>

namespace platec::tectonic_junction_coupling {

struct Inputs {
    uint32_t width = 0;
    uint32_t height = 0;
    const contract::TectonicJunction* junctions = nullptr;
    uint32_t junction_count = 0;
};

struct Outputs {
    std::vector<float> rift_influence;
    std::vector<float> transform_influence;
    std::vector<float> subduction_influence;
    std::vector<float> rift_tension_x;
    std::vector<float> rift_tension_y;
};

Outputs build(const Inputs& inputs);

} // namespace platec::tectonic_junction_coupling

#endif
