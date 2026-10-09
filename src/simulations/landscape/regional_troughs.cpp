#include "simulations/landscape/regional_troughs.hpp"

// Erosion plan Phase E5: the E5.0 contract (regional_troughs.hpp), which step E5.1
// implements. E5.0 leaves the hook inert: it returns before reading a field.

namespace landscape
{
void carvetroughs(const RegionalGlobalInputs& inputs, const RegionalWindow& window, const RegionalParameters& parameters,
                  const FineGlacial& glacial, const FineChannels& channels, const std::vector<uint8_t>& sea,
                  std::vector<float>& elevation_m)
{
    (void)inputs;
    (void)window;
    (void)glacial;
    (void)channels;
    (void)sea;
    (void)elevation_m;
    if (!(parameters.troughs.width_km > 0.0) || !glacial.active)
        return;
    // E5.1 implements the contract of regional_troughs.hpp here.
}
}
