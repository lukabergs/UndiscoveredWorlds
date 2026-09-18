#ifndef TECTONIC_PROVENANCE_HPP
#define TECTONIC_PROVENANCE_HPP

#include "tectonic_contract.hpp"

#include <cstdint>

namespace platec::provenance {

struct Inputs {
    uint32_t width = 0;
    uint32_t height = 0;
    const float* heightmap = nullptr;
    const uint32_t* plate_id = nullptr;
    const contract::PlateKinematics* plates = nullptr;
    uint32_t plate_count = 0;
    // Optional for legacy callers; the engine supplies current material-based classes.
    const uint8_t* crust_class = nullptr;
};

void compute_maps(const Inputs& inputs, uint8_t* convergence_score,
                  uint8_t* divergence_score, uint8_t* shear_score,
                  uint8_t* geologic_regime);

} // namespace platec::provenance

#endif
