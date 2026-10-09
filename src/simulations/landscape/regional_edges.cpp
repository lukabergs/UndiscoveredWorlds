#include "simulations/landscape/regional_edges.hpp"

#include <cstddef>

#include "simulations/landscape/regional_helpers.hpp"

// Erosion plan Phase R2, step R6: the contract (nearseacells and inlandcells moved here from
// regional_surface.cpp unchanged; seagate reproduces the 3 x 3 gate; breakstraightruns is
// a no-op). R8 implements the smooth gate and the breaker.

namespace landscape
{
using regionaldetail::aroundcell;
using regionaldetail::at;

std::vector<uint8_t> nearseacells(const RegionalGlobalInputs& inputs)
{
    const int width = inputs.width, height = inputs.height;
    std::vector<uint8_t> nearsea(static_cast<std::size_t>(width) * static_cast<std::size_t>(height), 0);
    for (int j = 0; j < height; j++)
        for (int i = 0; i < width; i++)
        {
            const std::size_t global = at(i, j, width);
            const bool lake = inputs.lakesurface_m[global] > 0.0f && inputs.network.sea[global] == 0;
            nearsea[global] = !lake && aroundcell(width, height, i, j, [&](std::size_t other) {
                return inputs.network.sea[other] != 0;
            }) ? 1 : 0;
        }
    return nearsea;
}

std::vector<float> inlandcells(const RegionalGlobalInputs& inputs)
{
    const std::vector<uint8_t> nearsea = nearseacells(inputs);
    std::vector<float> inland(nearsea.size());
    for (std::size_t cell = 0; cell < nearsea.size(); cell++)
        inland[cell] = nearsea[cell] != 0 ? 0.0f : 1.0f;
    return inland;
}

std::vector<uint8_t> seagate(const RegionalGlobalInputs& inputs, const RegionalWindow& box, const FineWarp& warp,
                             const RegionalParameters& parameters)
{
    (void)warp;
    (void)parameters;  // R8: edges.smoothgate selects the smooth gate.
    const int domainwidth = box.domainfinewidth(), domainheight = box.domainfineheight();
    const std::vector<uint8_t> nearsea = nearseacells(inputs);
    std::vector<uint8_t> gate(static_cast<std::size_t>(domainwidth) * static_cast<std::size_t>(domainheight), 0);
    for (int y = 0; y < domainheight; y++)
    {
        const int row = box.globalrow(y);
        for (int x = 0; x < domainwidth; x++)
            gate[at(x, y, domainwidth)] = nearsea[at(box.globalcolumn(x), row, inputs.width)];
    }
    return gate;
}

void breakstraightruns(const RegionalGlobalInputs& inputs, const RegionalWindow& box, const RegionalParameters& parameters,
                       const std::vector<uint8_t>& gate, std::vector<uint8_t>& mask)
{
    (void)inputs;
    (void)box;
    (void)gate;
    (void)mask;
    if (parameters.edges.straightmax <= 0)
        return;
    // R8 implements the breaker (regional_edges.hpp).
}
}
