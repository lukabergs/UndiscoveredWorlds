#ifndef TECTONIC_OCEAN_RELIEF_HPP
#define TECTONIC_OCEAN_RELIEF_HPP
#include "core/mesh/cell_grid.hpp"
#include <cstdint>
#include <vector>

namespace tectonics::ocean_relief {
struct Parameters {
    double crust_variation_m = 0; // Basalt thickness amplitude, not elevation amplitude.
    double wavelength_km = 600;
    double seamount_crust_m = 0;
    double seamount_radius_km = 350; // Gaussian sigma; compact support at 3 sigma.
    double seamount_spacing_km = 3500;
    double seamount_roughness = 0; // Relative coherent modulation; zero keeps elliptical profiles.
};
void validate(const Parameters& parameters);
// The newborn ocean crust thickness at each cell of a grid on the sphere, m: the ridge crust
// plus crust_variation_m x three octaves of sphere noise from wavelength_km (the value noise's
// deviation, octaves shorter than two cell spacings left out; sphere_noise::Octaves), plus
// seamount_crust_m x seamounts: one per cell of the coarse cube whose cells are about
// seamount_spacing_km wide, at a seeded offset within the middle 60 % of the cell, an ellipse
// of seeded orientation and aspect 0.75 to 1.25 by great-circle distance in its tangent
// frame. A pure function of the seed and the cell centres: no time argument or mutable state.
std::vector<double> crust(const mesh::CellGrid& grid, long seed, const Parameters& parameters);
}
#endif
