#ifndef TECTONIC_INITIALIZATION_HPP
#define TECTONIC_INITIALIZATION_HPP
#include "simulations/geology/plate_tectonics/tectonic_columns.hpp"

namespace tectonics::grid { class TectonicGrid; }
namespace tectonics::initialization {
// Detailed seeded terrain conditioned on opening contacts. This is a plausible
// starting state, not a reconstruction of its geological history.
std::vector<columns::Cell> terrain(uint32_t w, uint32_t h, double cell_m, double dt,
    const std::vector<double>& elevations, const std::vector<double>& ages,
    const std::vector<uint32_t>& owners, const std::vector<transport::Displacement>& velocity,
    double thermal_strength, double margin_width_m = 0, double boundary_influence = 1, double root_kg_m3 = 0,
    double root_max_thickness_m = 0, double root_delamination_width_m = 0, const grid::TectonicGrid* grid = nullptr);
}
#endif
