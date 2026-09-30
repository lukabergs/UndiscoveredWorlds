#pragma once

#include <algorithm>
#include <cmath>

namespace climateocean::detail
{
struct EkmanCellMobility { double alongStressSeconds, crossStressSeconds; };

// Exact cell-area average of the steady slab response for uniform stress.
// Spherical area is proportional to d(sin latitude), hence to df. Averaging
// r/(r*r+f*f) and f/(r*r+f*f) removes the centre-sample bias near the equator
// without an arbitrary tropical band, minimum |f|, or added momentum source.
inline EkmanCellMobility ekmanCellMobility(double dampingPerSecond,
    double northCoriolisPerSecond, double southCoriolisPerSecond)
{
    const double r=dampingPerSecond, fn=northCoriolisPerSecond, fs=southCoriolisPerSecond;
    if (!(r>0.0) || !std::isfinite(r) || !std::isfinite(fn) || !std::isfinite(fs)) return {};
    const double delta=fn-fs;
    if (std::abs(delta)<1e-8*std::max({r,std::abs(fn),std::abs(fs)}))
    {
        const double f=.5*(fn+fs), inverse=1.0/(r*r+f*f);
        return {r*inverse,f*inverse};
    }
    return {(std::atan(fn/r)-std::atan(fs/r))/delta,
        .5*std::log((r*r+fn*fn)/(r*r+fs*fs))/delta};
}
}
