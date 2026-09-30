#ifndef TECTONIC_THERMAL_HPP
#define TECTONIC_THERMAL_HPP

#include <algorithm>
#include <cmath>

namespace tectonics::thermal {

// Constant-property conductive plate cooling with isostatic water loading.
// This is an idealized model (Parsons & Sclater, 1977), not a fit to refs/.
// Physical assumptions: 125 km thermal plate, 1300 K contrast, expansivity
// 3e-5 / K, diffusivity 1e-6 m^2/s, mantle/water densities 3300/1030 kg/m^3.
// Returns accumulated thermal subsidence since birth, not an absolute depth.
inline double subsidence_m(double age_myr) {
    constexpr double pi = 3.14159265358979323846;
    constexpr double thickness_m = 125000.0;
    constexpr double density_factor = 3.0e-5 * 1300.0 * 3300.0 / (3300.0 - 1030.0);
    const double diffusion_length_sq = 1.0e-6 * std::max(0.0, age_myr) * 3.15576e13;
    const double u = diffusion_length_sq / (thickness_m * thickness_m);
    // At u <= .01 the basal correction is negligible (< micrometre), while
    // the half-space solution resolves t=0 without truncating the series tail.
    if (u <= .01) {
        return density_factor * 2.0 * std::sqrt(diffusion_length_sq / pi);
    }
    double sum = 0.0;
    for (int n = 1; n < 64; n += 2) {
        const double n_sq = static_cast<double>(n * n);
        sum += std::exp(-n_sq * pi * pi * u) / n_sq;
    }
    return density_factor * thickness_m * .5 * (1.0 - 8.0 * sum / (pi * pi));
}

inline double subsidence_increment_m(double age_myr, double dt_myr) {
    return subsidence_m(age_myr + dt_myr) - subsidence_m(age_myr);
}

} // namespace tectonics::thermal
#endif
