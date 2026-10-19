#include "simulations/landscape/regional_cirques.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>

#include "core/deterministic_random.hpp"
#include "simulations/landscape/level_blocks.hpp"
#include "simulations/landscape/regional_valleys.hpp"
#include "simulations/landscape/regional_helpers.hpp"
#include "simulations/landscape/upsampling.hpp"

// Erosion plan Phase E5: the E5.0 contract (regional_cirques.hpp), which step E5.2
// implements: hashed candidate points in the glaciated global cells' blocks, the band
// above the run's ELA, and the parabolic bowls.

namespace landscape
{
namespace
{
// Salt of the cirques' candidate draws (deterministiccombine(seed, salt, i, j, n)).
constexpr std::int64_t cirquesalt = 0x43495251;  // "CIRQ"
}

void carvecirques(const RegionalGlobalInputs& inputs, const RegionalWindow& window, const RegionalParameters& parameters,
                  const FineGlacial& glacial, const FineChannels& channels, const std::vector<uint8_t>& sea,
                  std::vector<float>& elevation_m, double bandlower, double bandupper)
{
    if (!(parameters.cirques.radius_km > 0.0) || !glacial.active)
        return;
    const CirqueParameters& cirques = parameters.cirques;
    const int width = window.domainfinewidth(), height = window.domainfineheight();
    const std::size_t count = static_cast<std::size_t>(width) * static_cast<std::size_t>(height);
    const std::size_t size = static_cast<std::size_t>(inputs.width) * static_cast<std::size_t>(inputs.height);
    const int k = window.scale;
    if (count == 0 || size == 0 || k < 1 || cirques.candidates < 1 || channels.channel.size() != count ||
        sea.size() != count || elevation_m.size() != count || glacial.runela_m.size() != count ||
        inputs.network.sea.size() != size)
        return;
    const bool withglaciation = inputs.glaciation.size() == size;

    // The bowl's radius in fine cells: kilometres on the ground over the fine north-south
    // spacing, at most k (the halo rule); under one cell the centre alone.
    const std::shared_ptr<const mesh::CellGrid> grid = regionalwindowgrid(window, inputs.radius_m);
    const double cell_m = northsouthspacing_m(*grid);
    const bool loop = bandupper > 0.0;
    const double radius = std::min(static_cast<double>(k), cirques.radius_km * 1000.0 / cell_m);
    if (!(radius > 0.0))
        return;
    if (bandupper > 0.0 && (radius < bandlower || radius >= bandupper))
        return;  // carved at another level of the refinement
    const int extent = static_cast<int>(std::floor(radius));
    const double radius2 = radius * radius;
    const float floorlevel = std::nextafter(inputs.sealevel, std::numeric_limits<float>::infinity());
    // The candidate's surface for the band test and the floor: the global base elevation
    // upsampled at the fine centres (a pure function of the global fields and the fine
    // position), not the solved elevation, which is not guaranteed to agree between windows
    // on a non-channel cell of the innermost halo ring where a bowl reaching a kept cell may
    // be centred (the review of E5.0). The bowl reads each lowered cell's own elevation,
    // which agrees on the kept cells.
    const std::vector<float> base = upsample(inputs.elevation_m, inputs.width, inputs.height, window);
    if (base.size() != count)
        return;

    // The lowest target per cell over every bowl that reaches it. The elevation is read
    // only here and written after the loop, so every candidate sees the surface before any
    // cut and the order of the candidates does not matter.
    std::vector<double> target(count, std::numeric_limits<double>::infinity());
    bool any = false;
    for (int row = 0; row < window.domainrows; row++)
    {
        const int j = window.domainrow0 + row;
        for (int column = 0; column < window.domaincolumns; column++)
        {
            const int i = (window.domaincolumn0 + column) % inputs.width;
            const std::size_t global = static_cast<std::size_t>(j) * static_cast<std::size_t>(inputs.width) +
                                       static_cast<std::size_t>(i);
            if (inputs.network.sea[global] != 0)
                continue;
            const double glaciation = withglaciation ? static_cast<double>(inputs.glaciation[global]) : 0.0;
            if (!(glaciation >= cirques.glaciationmin))
                continue;
            for (int n = 0; n < cirques.candidates; n++)
            {
                const std::uint64_t hash = deterministiccombine(inputs.seed, cirquesalt, i, j, n);
                const int ox = static_cast<int>((hash & 0xffffffffull) % static_cast<std::uint64_t>(k));
                const int oy = static_cast<int>((hash >> 32) % static_cast<std::uint64_t>(k));
                const int x = column * k + ox, y = row * k + oy;
                const std::size_t centre = static_cast<std::size_t>(y) * static_cast<std::size_t>(width) + static_cast<std::size_t>(x);
                if (sea[centre] != 0 || channels.channel[centre] != 0)
                    continue;
                const double z = static_cast<double>(base[centre]);
                const double ela = static_cast<double>(glacial.runela_m[centre]);
                if (!(z >= ela) || !(z <= ela + cirques.band_m))
                    continue;
                const double floor = std::max(static_cast<double>(floorlevel), std::max(ela, z - cirques.depth_m));
                any = true;
                // S3.5: in the loop the distance on the ground, the column offsets scaled by the
                // centre row's east-west over north-south cell size; the single pass keeps R9's cells.
                const double ratio = loop ? eastwestratio(*grid, centre) : 1.0;
                const int columns = loop ? columnextent(radius, ratio, k) : extent;
                for (int dy = -extent; dy <= extent; dy++)
                {
                    const int ny = y + dy;
                    if (ny < 0 || ny >= height)
                        continue;
                    for (int dx = -columns; dx <= columns; dx++)
                    {
                        const int nx = x + dx;
                        if (nx < 0 || nx >= width)
                            continue;
                        const double east = static_cast<double>(dx) * ratio;
                        const double d2 = east * east + static_cast<double>(dy * dy);
                        if (d2 > radius2)
                            continue;
                        const std::size_t next = static_cast<std::size_t>(ny) * static_cast<std::size_t>(width) + static_cast<std::size_t>(nx);
                        const double value = floor + (static_cast<double>(elevation_m[next]) - floor) * (d2 / radius2);
                        target[next] = std::min(target[next], value);
                    }
                }
            }
        }
    }
    if (!any)
        return;

    // S3.5: a channel cell keeps its value (its bed in the single pass; its cell mean in the loop,
    // which channelcellmeans, regional_valleys.hpp, lowers with the bowl around it, never back to
    // the bed).
    for (std::size_t cell = 0; cell < count; cell++)
    {
        if (sea[cell] != 0 || !std::isfinite(target[cell]) || channels.channel[cell] != 0)
            continue;
        const float value = std::max(floorlevel, static_cast<float>(target[cell]));
        elevation_m[cell] = std::min(elevation_m[cell], value);
    }
}

void carvecirques(const RegionalGlobalInputs& inputs, const mesh::CellGrid& window, const RegionalParameters& parameters,
                  const FineGlacial& glacial, const FineChannels& channels, const std::vector<uint8_t>& sea,
                  std::vector<float>& elevation_m, double bandlower, double bandupper)
{
    if (!(parameters.cirques.radius_km > 0.0) || !glacial.active)
        return;
    const CirqueParameters& cirques = parameters.cirques;
    const std::size_t count = window.cellcount();
    if (inputs.grid == nullptr || inputs.grid->kind() != mesh::GridKind::cube || window.kind() != mesh::GridKind::cube ||
        window.facecells() % inputs.grid->facecells() != 0)
        throw std::invalid_argument("carvecirques: a cube window refining the inputs' cube level 0");
    const int k = window.facecells() / inputs.grid->facecells();
    const mesh::CellGrid& levelzero = *inputs.grid;
    const std::size_t size = levelzero.cellcount();
    if (count == 0 || cirques.candidates < 1 || channels.channel.size() != count || sea.size() != count ||
        elevation_m.size() != count || glacial.runela_m.size() != count || inputs.network.sea.size() != size)
        return;
    const bool withglaciation = inputs.glaciation.size() == size;
    const double cell_m = northsouthspacing_m(window);
    const double radius = std::min(static_cast<double>(k), cirques.radius_km * 1000.0 / cell_m);
    if (!(radius > 0.0))
        return;
    if (bandupper > 0.0 && (radius < bandlower || radius >= bandupper))
        return;
    const double radius2 = radius * radius;
    const float floorlevel = std::nextafter(inputs.sealevel, std::numeric_limits<float>::infinity());
    const std::vector<float> base = upsample(inputs.elevation_m, levelzero, window);
    const double spacing = window.minimumspacing_m() / cell_m;
    const int rings = std::min(2 * k, static_cast<int>(std::ceil(radius / spacing)));
    std::vector<double> target(count, std::numeric_limits<double>::infinity());
    std::vector<std::uint32_t> stamp;
    std::vector<std::pair<mesh::CellId, int>> around;
    std::uint32_t id = 0;
    bool any = false;
    // The candidates of every level-0 cell of the window (its blocks): the hash of the seed and the
    // level-0 key draws the child at face-index offset (ox, oy) in the block.
    const WindowBlocks blocks = windowblocks(window, k);
    for (std::size_t b = 0; b < blocks.count(); b++)
    {
        const std::size_t global = levelzero.localofkey(blocks.parents[b]);
        if (inputs.network.sea[global] != 0)
            continue;
        const double glaciation = withglaciation ? static_cast<double>(inputs.glaciation[global]) : 0.0;
        if (!(glaciation >= cirques.glaciationmin))
            continue;
        for (int n = 0; n < cirques.candidates; n++)
        {
            const std::uint64_t hash =
                deterministiccombine(inputs.seed, cirquesalt, static_cast<std::int64_t>(blocks.parents[b]), std::int64_t{ -1 }, n);
            const int ox = static_cast<int>((hash & 0xffffffffull) % static_cast<std::uint64_t>(k));
            const int oy = static_cast<int>((hash >> 32) % static_cast<std::uint64_t>(k));
            const mesh::CellId centre = blocks.childrenof(b)[oy * k + ox];
            if (sea[centre] != 0 || channels.channel[centre] != 0)
                continue;
            const double z = static_cast<double>(base[centre]);
            const double ela = static_cast<double>(glacial.runela_m[centre]);
            if (!(z >= ela) || !(z <= ela + cirques.band_m))
                continue;
            const double floor = std::max(static_cast<double>(floorlevel), std::max(ela, z - cirques.depth_m));
            any = true;
            regionaldetail::cellswithin(window, centre, rings, stamp, ++id, around);
            for (const auto& [next, ring] : around)
            {
                const double distance = regionaldetail::centredistance(window, centre, next, cell_m);
                const double d2 = distance * distance;
                if (d2 > radius2)
                    continue;
                const double value = floor + (static_cast<double>(elevation_m[next]) - floor) * (d2 / radius2);
                target[next] = std::min(target[next], value);
            }
        }
    }
    if (!any)
        return;
    for (std::size_t cell = 0; cell < count; cell++)
    {
        if (sea[cell] != 0 || !std::isfinite(target[cell]) || channels.channel[cell] != 0)
            continue;
        const float value = std::max(floorlevel, static_cast<float>(target[cell]));
        elevation_m[cell] = std::min(elevation_m[cell], value);
    }
}
}
