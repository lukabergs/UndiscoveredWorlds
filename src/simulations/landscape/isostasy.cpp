#include "simulations/landscape/isostasy.hpp"

#include <cstddef>
#include <stdexcept>

namespace landscape
{
void applyisostasy(const IsostasyParameters& parameters, const std::vector<float>& erosion_m,
                   const std::vector<float>& deposition_m, std::vector<float>& elevation)
{
    if (erosion_m.size() != elevation.size() || deposition_m.size() != elevation.size())
        throw std::invalid_argument("applyisostasy: erosion, deposition and elevation differ in size");
    if (parameters.reboundfraction == 0.0 && parameters.loadingfraction == 0.0)
        return;
    for (std::size_t cell = 0; cell < elevation.size(); cell++)
    {
        const double change = parameters.reboundfraction * erosion_m[cell]
            - parameters.loadingfraction * deposition_m[cell];
        elevation[cell] = static_cast<float>(elevation[cell] + change);
    }
}
}
