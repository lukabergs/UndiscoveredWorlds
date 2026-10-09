#include "simulations/landscape/regional_valleys.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>

#include "simulations/landscape/sphere_metric.hpp"

// Erosion plan Phase R2: the R6 contract, which step R9 implements (regional_valleys.hpp):
// the incision of the pinned channels by the world's uplift or denudation rate over the
// elapsed time, and the valley floors after the legacy expandrivers and pasterivercircle.

namespace landscape
{
namespace
{
constexpr std::size_t nocell = std::numeric_limits<std::size_t>::max();

// The receiver of a global land cell on the drainage tree (channel_network.cpp's
// globalreceiver): nocell for a sea cell, a sink and a direction beyond a pole.
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

// The local global relief of one cell: the largest absolute difference between its
// elevation and its 8 neighbours' (columns wrap, rows beyond the poles do not exist), the
// expression of regional_surface.cpp's localrelief.
float localrelief(const std::vector<float>& elevation, int width, int height, std::size_t cell)
{
    const int i = static_cast<int>(cell % static_cast<std::size_t>(width));
    const int j = static_cast<int>(cell / static_cast<std::size_t>(width));
    const float own = elevation[cell];
    float largest = 0.0f;
    for (int dj = -1; dj <= 1; dj++)
    {
        const int row = j + dj;
        if (row < 0 || row >= height)
            continue;
        for (int di = -1; di <= 1; di++)
        {
            const int column = ((i + di) % width + width) % width;
            const float other = elevation[static_cast<std::size_t>(row) * static_cast<std::size_t>(width) +
                                          static_cast<std::size_t>(column)];
            largest = std::max(largest, std::fabs(other - own));
        }
    }
    return largest;
}
}

GlobalIncision globalincision(const RegionalGlobalInputs& inputs, const RegionalParameters& parameters)
{
    GlobalIncision incision;
    const ValleyParameters& valleys = parameters.valleys;
    const std::size_t size = static_cast<std::size_t>(inputs.width) * static_cast<std::size_t>(inputs.height);
    const bool uplift = valleys.upliftscale > 0.0 && inputs.upliftrate_m_myr.size() == size;
    const bool denudation = valleys.denudationratescale > 0.0 && inputs.denudation_m.size() == size &&
                            inputs.erosionduration_myr > 0.0;
    const double elapsed = parameters.relief.elapsed_myr;
    const GlobalNetwork& network = inputs.network;
    if ((!uplift && !denudation) || !(elapsed > 0.0) || size == 0 || network.sea.size() != size ||
        network.direction.size() != size || network.filledelevation_m.size() != size)
        return incision;

    // The lowest incised level: the next float above the sea level, so a land channel at the
    // mouth never sits at the sea level itself (the bed law's lo, channel_network.hpp).
    const float floorlevel = std::nextafter(inputs.sealevel, std::numeric_limits<float>::infinity());
    const bool withrelief = inputs.elevation_m.size() == size;

    // The rate and the capped depth I per land cell: at most F minus the lowest level and at
    // most the cell's local global relief, where the relief saturates once the plateau
    // itself lies within the knickpoint wave's reach (regional_valleys.hpp).
    std::vector<float> rate(size, 0.0f);
    std::vector<double> depth(size, 0.0);
    for (std::size_t cell = 0; cell < size; cell++)
    {
        if (network.sea[cell] != 0)
            continue;
        double value = 0.0;
        if (uplift)
            value = std::max(value, valleys.upliftscale * std::max(0.0, static_cast<double>(inputs.upliftrate_m_myr[cell])));
        if (denudation)
            value = std::max(value, valleys.denudationratescale * static_cast<double>(inputs.denudation_m[cell]) /
                                        inputs.erosionduration_myr);
        if (!std::isfinite(value))
            value = 0.0;
        rate[cell] = static_cast<float>(value);
        const double filled = static_cast<double>(network.filledelevation_m[cell]);
        double cap = filled - static_cast<double>(floorlevel);
        if (withrelief)
            cap = std::min(cap, static_cast<double>(localrelief(inputs.elevation_m, inputs.width, inputs.height, cell)));
        depth[cell] = std::max(0.0, std::min(value * elapsed, cap));
    }

    // Kahn's order over the forest of land links: a cell after all its donors.
    std::vector<std::size_t> receiver(size, nocell);
    std::vector<int> donors(size, 0);
    for (std::size_t cell = 0; cell < size; cell++)
    {
        receiver[cell] = receiverof(network, cell);
        if (receiver[cell] != nocell)
            donors[receiver[cell]]++;
    }
    // F'(r) = max(s+, min(F(r) - I(r), min over donors F'(u))) with s+ the next float above
    // the sea level, or F(r) where F(r) < s+; `lowest` holds the donors' minimum F' as they
    // are finished.
    std::vector<float> incised(network.filledelevation_m);
    std::vector<float> lowest(size, std::numeric_limits<float>::infinity());
    std::vector<uint8_t> done(size, 0);
    std::vector<std::size_t> ready;
    for (std::size_t cell = 0; cell < size; cell++)
        if (donors[cell] == 0)
            ready.push_back(cell);
    const auto finish = [&](std::size_t cell, bool withdonors) {
        done[cell] = 1;
        if (network.sea[cell] != 0)
            return;
        const float filled = network.filledelevation_m[cell];
        if (!(filled >= floorlevel))
            return;  // at or below the sea level (or not a number): no incision
        float value = static_cast<float>(static_cast<double>(filled) - depth[cell]);
        if (withdonors)
            value = std::min(value, lowest[cell]);
        incised[cell] = std::min(filled, std::max(floorlevel, value));
    };
    while (!ready.empty())
    {
        const std::size_t cell = ready.back();
        ready.pop_back();
        finish(cell, true);
        const std::size_t next = receiver[cell];
        if (next == nocell)
            continue;
        lowest[next] = std::min(lowest[next], incised[cell]);
        if (--donors[next] == 0)
            ready.push_back(next);
    }
    // Cells on a receiver loop (bad input) and downstream of one: their own capped I.
    for (std::size_t cell = 0; cell < size; cell++)
        if (done[cell] == 0)
            finish(cell, false);

    bool any = false;
    std::vector<float> depthout(size, 0.0f);
    for (std::size_t cell = 0; cell < size; cell++)
    {
        const double value = static_cast<double>(network.filledelevation_m[cell]) - static_cast<double>(incised[cell]);
        depthout[cell] = static_cast<float>(std::max(0.0, value));
        any = any || incised[cell] < network.filledelevation_m[cell];
    }
    if (!any)
        return incision;
    incision.active = true;
    incision.incision_m = std::move(depthout);
    incision.filledelevation_m = std::move(incised);
    incision.rate_m_myr = std::move(rate);
    return incision;
}

void carvevalleys(const RegionalGlobalInputs& inputs, const RegionalWindow& window, const RegionalParameters& parameters,
                  const FineChannels& channels, const std::vector<uint8_t>& sea, std::vector<float>& elevation_m)
{
    const ValleyParameters& valleys = parameters.valleys;
    if (!(valleys.floor_km > 0.0))
        return;
    const int width = window.domainfinewidth(), height = window.domainfineheight();
    const std::size_t count = static_cast<std::size_t>(width) * static_cast<std::size_t>(height);
    if (count == 0 || channels.channel.size() != count || sea.size() != count || elevation_m.size() != count)
        return;
    const int k = window.scale;
    const double reference = valleys.reference_m3s > 0.0 ? valleys.reference_m3s : 1.0;
    const bool sides = valleys.sideslope > 0.0 && valleys.reachcells > 0;
    const int reach = sides ? std::min(k, valleys.reachcells) : 0;
    const double cell_m = finepatchmetric(window, inputs.radius_m).dy_m;
    // The floor in fine cells: kilometres on the ground over the fine north-south spacing.
    const double floorcells = valleys.floor_km * 1000.0 / cell_m;
    const double floormax = std::min(static_cast<double>(k), std::max(0.0, valleys.floormax_km * 1000.0 / cell_m));
    const double rise = valleys.sideslope * cell_m;  // per fine cell beyond the floor
    const float floorlevel = std::nextafter(inputs.sealevel, std::numeric_limits<float>::infinity());

    // The lowest target per cell over every channel that reaches it.
    std::vector<double> target(count, std::numeric_limits<double>::infinity());
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
        {
            const std::size_t cell = static_cast<std::size_t>(y) * static_cast<std::size_t>(width) + static_cast<std::size_t>(x);
            if (channels.channel[cell] == 0 || sea[cell] != 0)
                continue;
            const double discharge = std::max(0.0, static_cast<double>(channels.discharge_m3s[cell]));
            const double halfwidth = std::min(floormax, floorcells * std::sqrt(discharge / reference));
            const int extent = std::min(k, std::max(static_cast<int>(std::floor(halfwidth)), reach));
            if (extent < 1)
                continue;
            const double bed = static_cast<double>(elevation_m[cell]);
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
                    const int distance = std::max(std::abs(dx), std::abs(dy));
                    double value;
                    if (static_cast<double>(distance) <= halfwidth)
                        value = bed;
                    else if (distance <= reach)
                        value = bed + (static_cast<double>(distance) - halfwidth) * rise;
                    else
                        continue;
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
