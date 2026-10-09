#include "simulations/landscape/regional_cirques.hpp"

// Erosion plan Phase E5: the E5.0 contract (regional_cirques.hpp), which step E5.2
// implements. E5.0 leaves the hook inert: it returns before reading a field.

namespace landscape
{
void carvecirques(const RegionalGlobalInputs& inputs, const RegionalWindow& window, const RegionalParameters& parameters,
                  const FineGlacial& glacial, const FineChannels& channels, const std::vector<uint8_t>& sea,
                  std::vector<float>& elevation_m)
{
    (void)inputs;
    (void)window;
    (void)glacial;
    (void)channels;
    (void)sea;
    (void)elevation_m;
    if (!(parameters.cirques.radius_km > 0.0) || !glacial.active)
        return;
    // E5.2 implements the contract of regional_cirques.hpp here.
}
}
