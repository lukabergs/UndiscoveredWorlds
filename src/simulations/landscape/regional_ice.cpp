#include "simulations/landscape/regional_ice.hpp"

// Erosion plan Phase E5: the E5.0 contract (regional_ice.hpp), which step E5.3
// implements. E5.0 leaves the hook inert: it returns before reading a field.

namespace landscape
{
void iceabovetheela(const RegionalGlobalInputs& inputs, const RegionalWindow& window, const RegionalParameters& parameters,
                    const FineGlacial& glacial, const std::vector<uint8_t>& sea, const std::vector<float>& elevation_m,
                    std::vector<float>& icethickness_m)
{
    (void)inputs;
    (void)window;
    (void)glacial;
    (void)sea;
    (void)elevation_m;
    (void)icethickness_m;
    if (!(parameters.ice.thicknessfraction > 0.0) || !glacial.active)
        return;
    // E5.3 implements the contract of regional_ice.hpp here.
}
}
