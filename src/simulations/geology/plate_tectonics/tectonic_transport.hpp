#ifndef TECTONIC_TRANSPORT_HPP
#define TECTONIC_TRANSPORT_HPP

#include <cstdint>
#include <vector>

namespace tectonics::transport {
struct Displacement { double x = 0, y = 0; }; // cells per update
// Positive, conservative, piecewise-parabolic finite-volume remap on a torus.
// Values are cell averages; no physical diffusion is applied here.
std::vector<double> remap(uint32_t width, uint32_t height,
    const std::vector<double>& values, const std::vector<Displacement>& displacement);
// A material moment (e.g. volume * age) uses the rock's transfer fractions.
// This conserves both integrals and bounds the transported specific quantity.
std::vector<double> remap_with_moment(uint32_t width, uint32_t height,
    const std::vector<double>& values, const std::vector<double>& moment,
    const std::vector<Displacement>& displacement, std::vector<double>& moved_moment);
}
#endif
