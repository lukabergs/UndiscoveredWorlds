#include "simulations/landscape/regional_valleys.hpp"

// Erosion plan Phase R2, step R6: the contract, inert (no incision, no valleys). R9
// implements both (regional_valleys.hpp).

namespace landscape
{
GlobalIncision globalincision(const RegionalGlobalInputs& inputs, const RegionalParameters& parameters)
{
    (void)inputs;
    GlobalIncision incision;
    if (parameters.valleys.upliftscale <= 0.0 && parameters.valleys.denudationratescale <= 0.0)
        return incision;
    // R9 implements the incision (regional_valleys.hpp).
    return incision;
}

void carvevalleys(const RegionalGlobalInputs& inputs, const RegionalWindow& window, const RegionalParameters& parameters,
                  const FineChannels& channels, const std::vector<uint8_t>& sea, std::vector<float>& elevation_m)
{
    (void)inputs;
    (void)window;
    (void)channels;
    (void)sea;
    (void)elevation_m;
    if (parameters.valleys.floorcells <= 0.0)
        return;
    // R9 implements the valley floors (regional_valleys.hpp).
}
}
