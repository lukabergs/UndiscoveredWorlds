#include "simulations/landscape/regional_coast.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>

#include "simulations/landscape/regional_helpers.hpp"

// Erosion plan Phase R2: R6 wrote the contract and coastfactors' formula, R7 the island
// sink and the values (configs/parameters/landscape.json, the note of the R7 keys).

namespace landscape
{
using regionaldetail::at;
using regionaldetail::wrapcolumn;

CoastFactors coastfactors(const RegionalGlobalInputs& inputs, const RegionalParameters& parameters)
{
    const CoastParameters& coast = parameters.coast;
    const std::size_t size = static_cast<std::size_t>(inputs.width) * static_cast<std::size_t>(inputs.height);
    const bool glacial = (coast.glacialamplitude_m > 0.0 || coast.glacialdropfraction > 0.0) && inputs.glaciation.size() == size;
    const bool denudation =
        coast.denudationamplitude_m > 0.0 && coast.denudationreference_m > 0.0 && inputs.denudation_m.size() == size;
    CoastFactors factors;
    if (!glacial && !denudation)
        return factors;
    factors.active = true;
    factors.amplitude_m.assign(size, 0.0f);
    factors.dropfactor.assign(size, 1.0f);
    for (std::size_t cell = 0; cell < size; cell++)
    {
        double amplitude = 0.0, factor = 1.0;
        if (glacial)
        {
            const double g = std::clamp(static_cast<double>(inputs.glaciation[cell]), 0.0, 1.0);
            amplitude += coast.glacialamplitude_m * g;
            factor = std::clamp(1.0 - coast.glacialdropfraction * g, 0.0, 1.0);
        }
        if (denudation)
        {
            const double d = std::max(0.0, static_cast<double>(inputs.denudation_m[cell]));
            amplitude += coast.denudationamplitude_m * std::min(1.0, d / coast.denudationreference_m);
        }
        factors.amplitude_m[cell] = static_cast<float>(amplitude);
        factors.dropfactor[cell] = static_cast<float>(factor);
    }
    return factors;
}

void sinksmallislands(const RegionalGlobalInputs& inputs, const RegionalWindow& box, const RegionalParameters& parameters,
                      const std::vector<uint8_t>& gate, std::vector<uint8_t>& mask)
{
    if (parameters.coast.islandmargin <= 0.0)
        return;
    const int width = inputs.width, height = inputs.height;
    const int k = box.scale;
    const int columns = box.domaincolumns, rows = box.domainrows;
    const int finewidth = box.domainfinewidth(), fineheight = box.domainfineheight();
    const bool wholeturn = columns >= width;
    const int margin = std::max(1, static_cast<int>(std::lround(parameters.coast.islandmargin * k)));
    const std::size_t size = static_cast<std::size_t>(width) * static_cast<std::size_t>(height);
    const bool islands = inputs.island.size() == size, mountainislands = inputs.mountainisland.size() == size;

    // A ring cell is land when the mask says so or when it lies beyond the box (rows end at
    // the box's edge or a pole; columns wrap only around a whole turn).
    auto ringland = [&](int x, int y) {
        if (y < 0 || y >= fineheight)
            return true;
        if (wholeturn)
            x = wrapcolumn(x, finewidth);
        else if (x < 0 || x >= finewidth)
            return true;
        return mask[at(x, y, finewidth)] == 0;
    };

    // The blocks to sink, all decided on the mask as it came in (so the result does not
    // depend on the order of the blocks), then sunk.
    std::vector<uint8_t> sink(static_cast<std::size_t>(columns) * static_cast<std::size_t>(rows), 0);
    for (int by = 0; by < rows; by++)
        for (int bx = 0; bx < columns; bx++)
        {
            const int column = wrapcolumn(box.domaincolumn0 + bx, width), row = box.domainrow0 + by;
            if (inputs.network.sea[at(column, row, width)] == 0)
                continue;
            bool land = false;
            for (int y = by * k; y < (by + 1) * k && !land; y++)
                for (int x = bx * k; x < (bx + 1) * k && !land; x++)
                    land = mask[at(x, y, finewidth)] == 0 && gate[at(x, y, finewidth)] != 0;
            if (!land)
                continue;
            // A one-cell island or a mountain island of the world in the 3 x 3 keeps the block.
            bool island = false;
            for (int dj = -1; dj <= 1 && !island; dj++)
            {
                const int j = row + dj;
                if (j < 0 || j >= height)
                    continue;
                for (int di = -1; di <= 1 && !island; di++)
                {
                    const std::size_t cell = at(wrapcolumn(column + di, width), j, width);
                    island = (islands && inputs.island[cell] != 0) || (mountainislands && inputs.mountainisland[cell] != 0);
                }
            }
            if (island)
                continue;
            // The ring at `margin` fine cells around the block.
            const int x0 = bx * k - margin, x1 = (bx + 1) * k - 1 + margin;
            const int y0 = by * k - margin, y1 = (by + 1) * k - 1 + margin;
            bool clear = true;
            for (int x = x0; x <= x1 && clear; x++)
                clear = !ringland(x, y0) && !ringland(x, y1);
            for (int y = y0 + 1; y < y1 && clear; y++)
                clear = !ringland(x0, y) && !ringland(x1, y);
            if (clear)
                sink[at(bx, by, columns)] = 1;
        }
    for (int by = 0; by < rows; by++)
        for (int bx = 0; bx < columns; bx++)
        {
            if (sink[at(bx, by, columns)] == 0)
                continue;
            for (int y = by * k; y < (by + 1) * k; y++)
                for (int x = bx * k; x < (bx + 1) * k; x++)
                {
                    const std::size_t cell = at(x, y, finewidth);
                    if (gate[cell] != 0)
                        mask[cell] = 1;
                }
        }
}
}
