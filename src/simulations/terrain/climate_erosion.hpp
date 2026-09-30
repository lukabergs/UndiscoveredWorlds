#pragma once
#include <cstdint>
#include <vector>

namespace climateerosion {
struct Forcing {
    std::vector<double> rain_mm_year, temperature_c, wind_m_s, snow_fraction, ice_fraction;
};
struct Parameters {
    double duration_myr=1, step_myr=.01;
    double river_coefficient=15, wave_m_myr=150, sediment_travel_km=600;
    double ice_m_myr=40, wind_m_myr=2, hillslope_m2_myr=2e8;
    double radius_m=6371000;
};
struct Result {
    std::vector<double> elevation_m, erosion_m, deposition_m;
    double eroded_m3=0, deposited_m3=0, relative_volume_residual=0;
};
// Surface rock-equivalent volume is conserved on a cell-centred spherical grid.
// This regional closure does not resolve glacier thickness, tectonic uplift or
// isostatic rebound. Four-season climate forcing is supplied by the caller.
Result evolve(uint32_t width,uint32_t height,const std::vector<double>& elevation,
              const Forcing& forcing,const Parameters& parameters);
}
