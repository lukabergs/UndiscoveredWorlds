#pragma once

#include <algorithm>
#include <cmath>

namespace climateenergy::detail
{
// Prescribed zonal cloud contribution to effective shortwave albedo. This is
// an annual-mean closure, not a prognostic cloud fraction or an Earth map.
inline double tropicalCloudAlbedo(double latitudeDegrees, double amplitude,
    double widthDegrees)
{
    if (!(widthDegrees > 0.0) || !(amplitude > 0.0)) return 0.0;
    const double latitude = latitudeDegrees / widthDegrees;
    return std::clamp(amplitude, 0.0, 0.2) * std::exp(-0.5 * latitude * latitude);
}
}
