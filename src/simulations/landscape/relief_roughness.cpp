#include "simulations/landscape/relief_roughness.hpp"

#include <stdexcept>

namespace landscape
{
// Erosion plan C3 implements it; amplitude 0 (the C0 default) changes nothing.
void addreliefroughness(const SphereMetric&, const ReliefRoughness& roughness, const std::vector<uint8_t>&, float,
                        std::vector<float>&)
{
    if (roughness.amplitude_m != 0.0)
        throw std::logic_error("addreliefroughness: not implemented (erosion plan C3)");
}
}
