#ifndef TECTONIC_PROVENANCE_HPP
#define TECTONIC_PROVENANCE_HPP

#include "tectonic_contract.hpp"

#include <cstdint>

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
};

void compute_maps(const Inputs& inputs, uint8_t* convergence_score,
                  uint8_t* divergence_score, uint8_t* shear_score,
                  uint8_t* geologic_regime);

} // namespace tectonics::provenance

#endif
