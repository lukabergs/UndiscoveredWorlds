#ifndef TECTONIC_DEFORMATION_HPP
#define TECTONIC_DEFORMATION_HPP
#include "simulations/geology/plate_tectonics/tectonic_transport.hpp"
#include <cstdint>
#include <vector>

namespace tectonics::deformation {
struct Parameters {
    double peak_distance_km = 300; // Maximum distributed shortening, inland from contact.
    double distance_variation = .2; // Smooth, seeded spatial variation; no per-step RNG.
    double subduction_coupling = .2; // Fraction of relative motion entering upper-plate deformation.
    double process_variation = 0, time_myr = 0;
    double process_wavelength_km = 1800, process_period_myr = 80;
};
// Cumulative compact strain curve: -d(taper)/dd peaks at peak_distance.
double taper(double distance_m, double peak_distance_m);
std::vector<transport::Displacement> distribute(uint32_t width, uint32_t height, double cell_m,
    const std::vector<uint32_t>& owners, const std::vector<double>& continental_fraction,
    const std::vector<transport::Displacement>& rigid, long seed, const Parameters& parameters);
}
#endif
