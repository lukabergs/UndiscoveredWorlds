#ifndef TECTONIC_OCEAN_RELIEF_HPP
#define TECTONIC_OCEAN_RELIEF_HPP
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
// Smooth periodic value noise in [-1,1]. No time argument or mutable RNG state.
double noise(double u, double v, uint32_t periods_x, uint32_t periods_y, uint32_t seed);
std::vector<double> crust(uint32_t width, uint32_t height, double cell_m, long seed,
                          const Parameters& parameters);
}
#endif
