#include "simulations/landscape/regional_coast.hpp"

#include <algorithm>
#include <cstddef>

#include "simulations/landscape/regional_helpers.hpp"

// Erosion plan Phase R2, step R6: the contract (coastfactors' formula; sinksmallislands a
// no-op). R7 implements the sink and sets the values.

namespace landscape
{
using regionaldetail::at;

CoastFactors coastfactors(const RegionalGlobalInputs& inputs, const RegionalParameters& parameters)
{
    const CoastParameters& coast = parameters.coast;
    const std::size_t size = static_cast<std::size_t>(inputs.width) * static_cast<std::size_t>(inputs.height);
    const bool glacial = (coast.glacialamplitude_m > 0.0 || coast.glacialdropfraction > 0.0) && inputs.glaciation.size() == size;
    const bool denudation =
        coast.denudationamplitude_m > 0.0 && coast.denudationreference_m > 0.0 && inputs.denudation_m.size() == size;
    CoastFactors factors;
    if (!glacial && !denudation)
        return factors;
    factors.active = true;
    factors.amplitude_m.assign(size, 0.0f);
    factors.dropfactor.assign(size, 1.0f);
    for (std::size_t cell = 0; cell < size; cell++)
    {
        double amplitude = 0.0, factor = 1.0;
        if (glacial)
        {
            const double g = std::clamp(static_cast<double>(inputs.glaciation[cell]), 0.0, 1.0);
            amplitude += coast.glacialamplitude_m * g;
            factor = std::clamp(1.0 - coast.glacialdropfraction * g, 0.0, 1.0);
        }
        if (denudation)
        {
            const double d = std::max(0.0, static_cast<double>(inputs.denudation_m[cell]));
            amplitude += coast.denudationamplitude_m * std::min(1.0, d / coast.denudationreference_m);
        }
        factors.amplitude_m[cell] = static_cast<float>(amplitude);
        factors.dropfactor[cell] = static_cast<float>(factor);
    }
    return factors;
}

void sinksmallislands(const RegionalGlobalInputs& inputs, const RegionalWindow& box, const RegionalParameters& parameters,
                      const std::vector<uint8_t>& gate, std::vector<uint8_t>& mask)
{
    (void)inputs;
    (void)box;
    (void)gate;
    (void)mask;
    if (parameters.coast.islandmargin <= 0.0)
        return;
    // R7 implements the sink (regional_coast.hpp).
}
}
