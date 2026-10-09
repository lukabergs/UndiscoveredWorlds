#include "simulations/landscape/regional_glacial.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>

#include "simulations/landscape/sphere_metric.hpp"
#include "simulations/landscape/upsampling.hpp"

// Erosion plan Phase E5, step E5.0 (regional_glacial.hpp): the fine ELA and glaciation
// fields shared by the glacial forms, and the ice catchment of the global network.

namespace landscape
{
namespace
{
constexpr std::size_t nocell = std::numeric_limits<std::size_t>::max();

// The receiver of a global land cell on the drainage tree (regional_valleys.cpp's
// receiverof): nocell for a sea cell, a sink and a direction beyond a pole.
std::size_t receiverof(const GlobalNetwork& network, std::size_t cell)
{
    const int direction = network.direction[cell];
    if (network.sea[cell] != 0 || direction < 1 || direction > directioncount)
        return nocell;
    const int i = static_cast<int>(cell % static_cast<std::size_t>(network.width));
    const int j = static_cast<int>(cell / static_cast<std::size_t>(network.width)) + directionoffsety(direction);
    if (j < 0 || j >= network.height)
        return nocell;
    return static_cast<std::size_t>(j) * static_cast<std::size_t>(network.width) +
           static_cast<std::size_t>((i + directionoffsetx(direction) + network.width) % network.width);
}
}

bool glacialformson(const RegionalParameters& parameters)
{
    return parameters.troughs.width_km > 0.0 || parameters.cirques.radius_km > 0.0 ||
           parameters.ice.thicknessfraction > 0.0;
}

FineGlacial regionalglacialfields(const RegionalGlobalInputs& inputs, const RegionalWindow& window,
                                  const RegionalParameters& parameters)
{
    FineGlacial glacial;
    const std::size_t size = static_cast<std::size_t>(inputs.width) * static_cast<std::size_t>(inputs.height);
    if (!glacialformson(parameters) || size == 0 || inputs.equilibriumline_m.size() != size)
        return glacial;
    glacial.active = true;
    glacial.runela_m = upsample(inputs.equilibriumline_m, inputs.width, inputs.height, window);
    glacial.presentela_m = glacial.runela_m;
    const double offset = std::isfinite(inputs.presentelaoffset_m) ? inputs.presentelaoffset_m : 0.0;
    for (float& value : glacial.presentela_m)
        value = static_cast<float>(static_cast<double>(value) + offset);
    if (inputs.glaciation.size() == size)
    {
        glacial.glaciation = upsample(inputs.glaciation, inputs.width, inputs.height, window);
        for (float& value : glacial.glaciation)
            value = std::clamp(value, 0.0f, 1.0f);
    }
    else
        glacial.glaciation.assign(glacial.runela_m.size(), 0.0f);
    return glacial;
}

std::vector<double> glacialcatchment_km2(const RegionalGlobalInputs& inputs)
{
    std::vector<double> catchment;
    const std::size_t size = static_cast<std::size_t>(inputs.width) * static_cast<std::size_t>(inputs.height);
    const GlobalNetwork& network = inputs.network;
    if (size == 0 || inputs.glaciation.size() != size || network.sea.size() != size || network.direction.size() != size ||
        !(inputs.radius_m > 0.0))
        return catchment;
    const SphereMetric metric(inputs.width, inputs.height, inputs.radius_m);
    catchment.assign(size, 0.0);
    for (std::size_t cell = 0; cell < size; cell++)
    {
        if (network.sea[cell] != 0)
            continue;
        const double share = std::clamp(static_cast<double>(inputs.glaciation[cell]), 0.0, 1.0);
        catchment[cell] = std::isfinite(share) ? share * metric.area_m2(cell) / 1.0e6 : 0.0;
    }
    // Kahn's order over the forest of land links: a cell's donors are added before it is
    // passed on to its receiver.
    std::vector<std::size_t> receiver(size, nocell);
    std::vector<int> donors(size, 0);
    for (std::size_t cell = 0; cell < size; cell++)
    {
        receiver[cell] = receiverof(network, cell);
        if (receiver[cell] != nocell)
            donors[receiver[cell]]++;
    }
    std::vector<std::size_t> ready;
    for (std::size_t cell = 0; cell < size; cell++)
        if (donors[cell] == 0)
            ready.push_back(cell);
    while (!ready.empty())
    {
        const std::size_t cell = ready.back();
        ready.pop_back();
        const std::size_t next = receiver[cell];
        if (next == nocell)
            continue;
        if (network.sea[next] == 0)
            catchment[next] += catchment[cell];
        if (--donors[next] == 0)
            ready.push_back(next);
    }
    // Cells on a receiver loop keep their own term (their donors' sums were added already);
    // a sea cell carries nothing.
    for (std::size_t cell = 0; cell < size; cell++)
        if (network.sea[cell] != 0)
            catchment[cell] = 0.0;
    return catchment;
}
}
