#ifndef TECTONIC_PROVENANCE_HPP
#define TECTONIC_PROVENANCE_HPP

#include "simulations/geology/plate_tectonics/tectonic_contract.hpp"
#include "simulations/geology/plate_tectonics/tectonic_grid.hpp"

#include <cstdint>

namespace tectonics::provenance {

struct Inputs {
    const float* heightmap = nullptr;
    const uint32_t* plate_id = nullptr;
    uint32_t plate_count = 0;
    // Optional for synthetic inputs; the engine supplies current material-based classes.
    const uint8_t* crust_class = nullptr;
    // Zero preserves dimensionless cells/update thresholds. Positive converts velocities in
    // cell spacings per update to cm/year (cellspacing_m / dt), so physical boundary classes
    // survive dt/grid changes.
    double velocity_to_cm_yr = 0.0;
    // The engine's cell table (required: nothing is computed without it).
    const grid::TectonicGrid* grid = nullptr;
    // Per plate the rotation vector in radians per update (null: every plate at rest).
    const mesh::Vec3* plate_rotation = nullptr;
};

void compute_maps(const Inputs& inputs, uint8_t* convergence_score,
                  uint8_t* divergence_score, uint8_t* shear_score,
                  uint8_t* geologic_regime);

} // namespace tectonics::provenance

#endif
