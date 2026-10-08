#include "simulations/landscape/regional_water.hpp"

#include <algorithm>
#include <cstddef>

// Erosion plan steps E3 and R0 (docs/architecture/erosion-hydrology-plan.md, Phase E
// and Phase R): moved from regional_pass.cpp; R1 wets the channel cells in lakes.

namespace landscape
{
namespace
{
std::size_t at(int x, int y, int width)
{
    return static_cast<std::size_t>(y) * static_cast<std::size_t>(width) + static_cast<std::size_t>(x);
}

int wrapcolumn(int column, int width) { return ((column % width) + width) % width; }

// True when global cell (i, j) or one of its 8 neighbours satisfies `test` (columns
// wrap, rows beyond the poles do not exist).
template<typename Test>
bool aroundcell(int width, int height, int i, int j, Test test)
{
    for (int dj = -1; dj <= 1; dj++)
    {
        const int row = j + dj;
        if (row < 0 || row >= height)
            continue;
        for (int di = -1; di <= 1; di++)
            if (test(at(wrapcolumn(i + di, width), row, width)))
                return true;
    }
    return false;
}
}

FineWater regionalwater(const RegionalGlobalInputs& inputs, const RegionalWindow& window, const FineWarp& warp,
                        const RegionalParameters& parameters, const std::vector<uint8_t>& sea,
                        const FineChannels& channels, std::vector<float>& elevation_m)
{
    (void)warp;
    (void)parameters;
    const int width = inputs.width, height = inputs.height;
    const int domainwidth = window.domainfinewidth(), domainheight = window.domainfineheight();
    const std::size_t globalsize = static_cast<std::size_t>(width) * static_cast<std::size_t>(height);
    const std::size_t count = static_cast<std::size_t>(domainwidth) * static_cast<std::size_t>(domainheight);

    std::vector<uint8_t> nearice(globalsize, 0);
    for (int j = 0; j < height; j++)
        for (int i = 0; i < width; i++)
            nearice[at(i, j, width)] = aroundcell(width, height, i, j, [&](std::size_t other) {
                return inputs.icethickness_m[other] > 0.0f;
            }) ? 1 : 0;
    const std::vector<float> ice = upsample(inputs.icethickness_m, width, height, window);

    FineWater water;
    water.lakesurface_m.assign(count, 0.0f);
    water.icethickness_m.assign(count, 0.0f);
    for (int y = 0; y < domainheight; y++)
    {
        const int row = window.globalrow(y);
        for (int x = 0; x < domainwidth; x++)
        {
            const std::size_t cell = at(x, y, domainwidth);
            const std::size_t global = at(window.globalcolumn(x), row, width);
            // R1: a channel cell carries the lake's water where its bed is at or below
            // the higher of the lake surface and the global cell's filled elevation,
            // and its elevation is lowered to at most the surface. Beds are pinned at
            // or below the filled elevation, which equals the surface of a lake that
            // spills and lies above that of a closed lake under its sill; the strict
            // rule left rivers dry through their lakes. A shore cell whose bed follows
            // a start surface above that level stays dry.
            const float surfacelevel = inputs.lakesurface_m[global];
            if (sea[cell] == 0 && surfacelevel > 0.0f)
            {
                if (channels.channel[cell] != 0)
                {
                    const float level = std::max(surfacelevel, inputs.network.filledelevation_m[global]);
                    if (elevation_m[cell] <= level)
                    {
                        water.lakesurface_m[cell] = surfacelevel;
                        elevation_m[cell] = std::min(elevation_m[cell], surfacelevel);
                    }
                }
                else if (elevation_m[cell] < surfacelevel)
                    water.lakesurface_m[cell] = surfacelevel;
            }
            water.icethickness_m[cell] = nearice[global] != 0 ? std::max(0.0f, ice[cell]) : 0.0f;
        }
    }
    return water;
}
}
