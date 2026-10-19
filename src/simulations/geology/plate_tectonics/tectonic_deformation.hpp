#ifndef TECTONIC_DEFORMATION_HPP
#define TECTONIC_DEFORMATION_HPP
#include "core/mesh/cell_grid.hpp"
#include <cstdint>
#include <vector>

namespace tectonics::columns { class Metric; }
namespace tectonics::deformation {
struct Parameters {
    double peak_distance_km = 300; // Maximum distributed shortening, inland from contact.
    double distance_variation = .2; // Smooth, seeded spatial variation; no per-step RNG.
    double subduction_coupling = .2; // Fraction of relative motion entering upper-plate deformation.
    double process_variation = 0, time_myr = 0;
    double process_wavelength_km = 1800, process_period_myr = 80;
    // The spatial variation's wavelength: the torus pattern's four cycles per world width,
    // pi R / 2 on a sphere of Earth's radius (scenario deformation_variation_wavelength_km).
    double distance_variation_wavelength_km = 10007;
};
// Cumulative compact strain curve: -d(taper)/dd peaks at peak_distance.
double taper(double distance_m, double peak_distance_m);
// Displacements per update as 3-D tangent vectors in metres (the result too): the rigid
// displacement plus the distributed share of the closing at convergent continental contacts.
// The peak distance's variation is read at each cell's centre on the sphere.
std::vector<mesh::Vec3> distribute(const columns::Metric& metric, const std::vector<uint32_t>& owners,
    const std::vector<double>& continental_fraction, const std::vector<mesh::Vec3>& rigid_m, long seed,
    const Parameters& parameters);
}
#endif
