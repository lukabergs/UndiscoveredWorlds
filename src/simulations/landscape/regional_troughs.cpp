#include "simulations/landscape/regional_troughs.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>


// Erosion plan Phase E5: the E5.0 contract (regional_troughs.hpp), which step E5.1
// implements: glacial troughs, a power-law cross-section about every glaciated land channel
// cell, half-width from the ice catchment of its global cell.

namespace landscape
{
void carvetroughs(const RegionalGlobalInputs& inputs, const RegionalWindow& window, const RegionalParameters& parameters,
                  const FineGlacial& glacial, const FineChannels& channels, const std::vector<uint8_t>& sea,
                  std::vector<float>& elevation_m)
{
    const TroughParameters& troughs = parameters.troughs;
    if (!(troughs.width_km > 0.0) || !glacial.active)
        return;
    const int width = window.domainfinewidth(), height = window.domainfineheight();
    const std::size_t count = static_cast<std::size_t>(width) * static_cast<std::size_t>(height);
    if (count == 0 || channels.channel.size() != count || channels.globalcell.size() != count || sea.size() != count ||
        elevation_m.size() != count || glacial.glaciation.size() != count)
        return;
    // The ice catchment per global cell, once per call (a pure function of the global fields).
    const std::vector<double> catchment = glacialcatchment_km2(inputs);
    if (catchment.empty())
        return;
    const int k = window.scale;
    const double reference = troughs.reference_km2 > 0.0 ? troughs.reference_km2 : 1.0;
    const double cell_m = finepatchmetric(window, inputs.radius_m).dy_m;
    if (!(cell_m > 0.0))
        return;
    // The half-width in fine cells: kilometres on the ground over the fine north-south
    // spacing, capped at widthmax_km and at k (one global cell, the halo rule).
    const double widthcells = troughs.width_km * 1000.0 / cell_m;
    const double widthmax = std::min(static_cast<double>(k), std::max(0.0, troughs.widthmax_km * 1000.0 / cell_m));
    const float floorlevel = std::nextafter(inputs.sealevel, std::numeric_limits<float>::infinity());

    // The lowest target per cell over every axis that reaches it. The beds are the channel
    // cells' elevations, which no trough changes, and the targets are applied after every
    // axis is visited, so the order of the axes does not matter.
    std::vector<double> target(count, std::numeric_limits<double>::infinity());
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
        {
            const std::size_t cell = static_cast<std::size_t>(y) * static_cast<std::size_t>(width) + static_cast<std::size_t>(x);
            if (channels.channel[cell] == 0 || sea[cell] != 0)
                continue;
            if (!(glacial.glaciation[cell] >= troughs.glaciationmin))
                continue;
            const int32_t global = channels.globalcell[cell];
            if (global < 0 || static_cast<std::size_t>(global) >= catchment.size())
                continue;
            const double area = catchment[static_cast<std::size_t>(global)];
            if (!(area > 0.0))
                continue;
            const double halfwidth = std::min(widthmax, widthcells * std::pow(area / reference, troughs.widthexponent));
            if (!(halfwidth >= 1.0))
                continue;  // under one fine cell: only the axis, which is a channel
            const int extent = std::min(k, static_cast<int>(std::floor(halfwidth)));
            const double bed = static_cast<double>(elevation_m[cell]);
            if (!std::isfinite(bed))
                continue;
            for (int dy = -extent; dy <= extent; dy++)
            {
                const int ny = y + dy;
                if (ny < 0 || ny >= height)
                    continue;
                for (int dx = -extent; dx <= extent; dx++)
                {
                    const int nx = x + dx;
                    if (nx < 0 || nx >= width)
                        continue;
                    const double distance = std::sqrt(static_cast<double>(dx * dx + dy * dy));
                    if (distance > halfwidth)
                        continue;
                    const double value = bed + troughs.wall_m * std::pow(distance / halfwidth, troughs.exponent);
                    const std::size_t next = static_cast<std::size_t>(ny) * static_cast<std::size_t>(width) + static_cast<std::size_t>(nx);
                    target[next] = std::min(target[next], value);
                }
            }
        }

    for (std::size_t cell = 0; cell < count; cell++)
    {
        if (channels.channel[cell] != 0 || sea[cell] != 0 || !std::isfinite(target[cell]))
            continue;
        const float value = std::max(floorlevel, static_cast<float>(target[cell]));
        elevation_m[cell] = std::min(elevation_m[cell], value);
    }
}
}
