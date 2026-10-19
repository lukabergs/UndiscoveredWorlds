#ifndef TECTONIC_VARIATION_HPP
#define TECTONIC_VARIATION_HPP
#include "simulations/geology/plate_tectonics/tectonic_sphere_noise.hpp"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>

namespace tectonics::variation {
inline void validate(double amplitude,double wavelength_km,double period_myr) {
    if (!std::isfinite(amplitude) || amplitude<0 || amplitude>.75 ||
        !std::isfinite(wavelength_km) || wavelength_km<=0 || !std::isfinite(period_myr) || period_myr<=0)
        throw std::runtime_error("Invalid continuous process variation");
}
// Smooth, periodic environmental heterogeneity on the sphere. Absolute geological time and
// independent channel phases avoid flicker and update-order/RNG dependence. Three plane waves
// through the globe, read at the unit vector p of a sphere of radius_km: wave vectors
// (2 pi / wavelength) (e1 + e2), (2 e1 - e2) and (e1 + 2 e3), the torus pattern's lengths, in a
// frame seeded by the seed and channel (sphere_noise::frame), so the value is continuous across
// the cube's face edges and corners and periodic in time with period_myr. The third wave leaves
// the e1-e2 plane: three coplanar wave vectors would make the pattern mirror-symmetric about
// the great circle normal to e3 (E6).
inline double multiplier(const mesh::Vec3& p,double radius_km,double time_myr,
        long seed,uint32_t channel,double amplitude,double wavelength_km,double period_myr) {
    if (amplitude==0) return 1;
    constexpr double tau=6.28318530717958647692;
    const uint32_t bits=(static_cast<uint32_t>(seed)^(channel*2246822519U))*1664525U+1013904223U;
    const double phase=tau*bits/4294967296.0,t=tau*time_myr/period_myr;
    const auto axes=sphere_noise::frame(bits);
    const double k=tau*radius_km/wavelength_km;
    const double x=k*mesh::dot(p,axes[0]),y=k*mesh::dot(p,axes[1]),z=k*mesh::dot(p,axes[2]);
    const double value=.5*std::sin(x+y+phase+.6*std::sin(t))
        +.3*std::sin(2*x-y+1.7*phase-t)
        +.2*std::cos(x+2*z-.8*phase+.4*std::cos(t));
    return 1+amplitude*value;
}
}
#endif
