#ifndef TECTONIC_VARIATION_HPP
#define TECTONIC_VARIATION_HPP
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
// Smooth, periodic environmental heterogeneity. Absolute geological time and
// independent channel phases avoid flicker and update-order/RNG dependence.
inline double multiplier(double x,double y,double width_km,double height_km,double time_myr,
        long seed,uint32_t channel,double amplitude,double wavelength_km,double period_myr) {
    if (amplitude==0) return 1;
    constexpr double tau=6.28318530717958647692;
    const uint32_t bits=(static_cast<uint32_t>(seed)^(channel*2246822519U))*1664525U+1013904223U;
    const double phase=tau*bits/4294967296.0,t=tau*time_myr/period_myr;
    const double kx=std::max(1.0,std::round(width_km/wavelength_km));
    const double ky=std::max(1.0,std::round(height_km/wavelength_km));
    const double value=.5*std::sin(tau*(kx*x+ky*y)+phase+.6*std::sin(t))
        +.3*std::sin(tau*(2*kx*x-ky*y)+1.7*phase-t)
        +.2*std::cos(tau*(kx*x+2*ky*y)-.8*phase+.4*std::cos(t));
    return 1+amplitude*value;
}
}
#endif
