#include "simulations/landscape/regional_coast.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <stdexcept>

#include "core/deterministic_random.hpp"
#include "simulations/landscape/regional_helpers.hpp"

// Erosion plan Phase R2: R6 wrote the contract and coastfactors' formula, R7 the island
// sink and the values (configs/parameters/landscape.json, the note of the R7 keys).
// Multiscale plan S3.15: the shoreline's displacement and the sea pockets.

namespace landscape
{
using regionaldetail::at;
using regionaldetail::wrapcolumn;

namespace
{
constexpr double pi = 3.14159265358979323846;
// Salts of the shoreline displacement's two components (S3.15).
constexpr std::uint64_t shorewarpxsalt = 0xc0a6;
constexpr std::uint64_t shorewarpysalt = 0xc0a7;
// The displacement's bound in standard deviations of its band (a method constant, ShorelineWarp).
constexpr double shorewarpbound = 2.0;

// The octaves o with 2^(o + 1) <= scale: those of at least two fine cells (RegionalNoise's count).
int octavesof(int scale)
{
    int octaves = 0;
    while ((std::int64_t{2} << octaves) <= scale)
        octaves++;
    return octaves;
}

// A fine cell centre's latitude and longitude from its absolute indices (RegionalWindow's
// finelatitude_rad and finelongitude_rad expressions).
double finelatitude(std::int64_t finerow, int fineheight)
{
    const int j = static_cast<int>(finerow);
    const double north = pi / 2.0 - pi * j / fineheight;
    const double south = pi / 2.0 - pi * (j + 1) / fineheight;
    return 0.5 * (north + south);
}

double finelongitude(std::int64_t finecolumn, int finewidth)
{
    const std::int64_t w = finewidth;
    return 2.0 * pi * (static_cast<double>(((finecolumn % w) + w) % w) + 0.5) / finewidth;
}

template<typename Value>
std::vector<double> displaced(const std::vector<Value>& field, const RegionalWindow& box, const ShorelineWarp& warp)
{
    const int width = box.domainfinewidth(), height = box.domainfineheight();
    const std::size_t count = static_cast<std::size_t>(width) * static_cast<std::size_t>(height);
    if (field.size() != count)
        throw std::invalid_argument("displacedsample: the field must cover the box's fine cells");
    std::vector<double> result(count);
    if (!warp.active())
    {
        for (std::size_t cell = 0; cell < count; cell++)
            result[cell] = static_cast<double>(field[cell]);
        return result;
    }
    const bool wholeturn = box.domaincolumns >= box.worldwidth;
    const int k = box.scale;
    auto column = [&](std::int64_t x) {
        return wholeturn ? static_cast<int>(((x % width) + width) % width) : static_cast<int>(std::clamp<std::int64_t>(x, 0, width - 1));
    };
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
        {
            double dx = 0.0, dy = 0.0;
            warp.at(box.absolutefinecolumn(x), box.absolutefinerow(y), dx, dy);
            const double fx = x + dx * k, fy = y + dy * k;
            const double x0 = std::floor(fx), y0 = std::floor(fy);
            const double tx = fx - x0, ty = fy - y0;
            const std::int64_t ix = static_cast<std::int64_t>(x0), iy = static_cast<std::int64_t>(y0);
            const int c0 = column(ix), c1 = column(ix + 1);
            const int r0 = static_cast<int>(std::clamp<std::int64_t>(iy, 0, height - 1));
            const int r1 = static_cast<int>(std::clamp<std::int64_t>(iy + 1, 0, height - 1));
            const double a = static_cast<double>(field[at(c0, r0, width)]), b = static_cast<double>(field[at(c1, r0, width)]);
            const double c = static_cast<double>(field[at(c0, r1, width)]), d = static_cast<double>(field[at(c1, r1, width)]);
            const double north = a + tx * (b - a), south = c + tx * (d - c);
            result[at(x, y, width)] = north + ty * (south - north);
        }
    return result;
}
}

ShorelineWarp::ShorelineWarp(int worldwidth, int worldheight, int scale, int parentscale, std::uint64_t seed,
                             const CoastParameters& coast)
    : scale_(scale), finewidth_(worldwidth * scale), fineheight_(worldheight * scale), amplitude_(coast.warpcells)
{
    const int octaves = octavesof(scale), first = parentscale > 1 ? octavesof(parentscale) : 0;
    if (!(coast.warpcells > 0.0) || worldheight < 1 || octaves <= first)
        return;
    FractalNoise noise;
    noise.radius_km = 1.0;
    noise.wavelength_km = pi / worldheight;
    noise.octaves = octaves;
    noise.persistence = coast.warppersistence;
    // The band: octave o >= first has a wavelength below that of octave first - 1 (powers of two:
    // exact), so longest_km excludes the parent's octaves.
    noise.longest_km = first > 0 ? noise.wavelength_km / std::ldexp(1.0, first - 1) : 0.0;
    double band = 0.0;
    for (int o = first; o < octaves; o++)
        band += std::pow(coast.warppersistence, 2 * o);
    noise.seed = deterministiccombine(seed, shorewarpxsalt);
    x_.emplace(noise);
    noise.seed = deterministiccombine(seed, shorewarpysalt);
    y_.emplace(noise);
    bound_ = shorewarpbound * coast.warpcells * std::sqrt(band);
    active_ = bound_ > 0.0;
}

void ShorelineWarp::at(std::int64_t finecolumn, std::int64_t finerow, double& dx, double& dy) const
{
    dx = dy = 0.0;
    if (!active_)
        return;
    double px = 0.0, py = 0.0, pz = 0.0;
    unitvector(finelatitude(finerow, fineheight_), finelongitude(finecolumn, finewidth_), px, py, pz);
    const double u = amplitude_ * x_->sum(px, py, pz) / octavedeviation;
    const double v = amplitude_ * y_->sum(px, py, pz) / octavedeviation;
    const double length = std::sqrt(u * u + v * v);
    const double factor = length > 0.0 ? bound_ * std::tanh(length / bound_) / length : 1.0;
    dx = factor * u;
    dy = factor * v;
}

int ShorelineWarp::reachcells() const
{
    return active_ ? static_cast<int>(std::ceil(bound_ * scale_)) + 1 : 0;
}

std::vector<double> displacedsample(const std::vector<double>& field, const RegionalWindow& box, const ShorelineWarp& warp)
{
    return displaced(field, box, warp);
}

std::vector<double> displacedsample(const std::vector<float>& field, const RegionalWindow& box, const ShorelineWarp& warp)
{
    return displaced(field, box, warp);
}

namespace
{
// The blocks of a box within `reach` global cells of block (bx, by): the region of seapockets and
// pocketbasins (columns wrap around a whole turn; blocks beyond the box do not exist).
struct Region
{
    int reach = 0, columns = 0, rows = 0, k = 1;
    bool wholeturn = false;
    bool within(int nbx, int nby, int bx, int by) const
    {
        if (std::abs(nby - by) > reach)
            return false;
        if (!wholeturn)
            return std::abs(nbx - bx) <= reach;
        const int d = ((nbx - bx) % columns + columns) % columns;
        return std::min(d, columns - d) <= reach;
    }
};

// The breadth-first search of seapockets and pocketbasins: from the cells in `queue` (marked `id`)
// through the 4-neighbours that `pass` admits inside the region of block (bx, by).
// The global lakes as the pocket search and the basins see them (S3.15, the merged run 47: a basin took
// a lake's bed below the sea level and lifted it out of the lake): per global cell of the world, the
// lake's surface where it is a lake cell (lake surface above zero and not sea), else the highest
// surface of a lake cell among its 8 neighbours, else 0; and whether it is a lake cell. A pure function
// of the global fields, so the rules that read it stay seamless.
struct GlobalLakes
{
    std::vector<uint8_t> own;
    std::vector<float> surface;
    bool empty() const { return own.empty(); }
};

GlobalLakes globallakes(const RegionalGlobalInputs& inputs)
{
    GlobalLakes lakes;
    const int width = inputs.width, height = inputs.height;
    const std::size_t size = static_cast<std::size_t>(width) * static_cast<std::size_t>(height);
    if (inputs.lakesurface_m.size() != size)
        return lakes;
    auto lake = [&](std::size_t cell) { return inputs.lakesurface_m[cell] > 0.0f && inputs.network.sea[cell] == 0; };
    bool any = false;
    lakes.own.assign(size, 0);
    lakes.surface.assign(size, 0.0f);
    for (int j = 0; j < height; j++)
        for (int i = 0; i < width; i++)
        {
            const std::size_t cell = at(i, j, width);
            if (lake(cell))
            {
                lakes.own[cell] = 1;
                lakes.surface[cell] = inputs.lakesurface_m[cell];
                any = true;
                continue;
            }
            float highest = 0.0f;
            regionaldetail::aroundcell(width, height, i, j, [&](std::size_t other) {
                if (lake(other))
                    highest = std::max(highest, inputs.lakesurface_m[other]);
                return false;
            });
            lakes.surface[cell] = highest;
        }
    if (!any)
        lakes = GlobalLakes{};
    return lakes;
}

template<typename Pass>
void spread(std::vector<std::size_t>& queue, std::vector<int>& mark, int id, int finewidth, int fineheight,
            const Region& region, int bx, int by, Pass pass)
{
    for (std::size_t next = 0; next < queue.size(); next++)
    {
        const int x = static_cast<int>(queue[next] % static_cast<std::size_t>(finewidth));
        const int y = static_cast<int>(queue[next] / static_cast<std::size_t>(finewidth));
        const int nx[4] = {x - 1, x + 1, x, x};
        const int ny[4] = {y, y, y - 1, y + 1};
        for (int n = 0; n < 4; n++)
        {
            int cx = nx[n];
            const int cy = ny[n];
            if (cy < 0 || cy >= fineheight)
                continue;
            if (region.wholeturn)
                cx = ((cx % finewidth) + finewidth) % finewidth;
            else if (cx < 0 || cx >= finewidth)
                continue;
            const std::size_t cell = at(cx, cy, finewidth);
            if (mark[cell] == id || !pass(cell) || !region.within(cx / region.k, cy / region.k, bx, by))
                continue;
            mark[cell] = id;
            queue.push_back(cell);
        }
    }
}
}

std::vector<uint8_t> seapockets(const RegionalGlobalInputs& inputs, const RegionalWindow& box, int reach,
                                const std::vector<uint8_t>& sea, const std::vector<uint8_t>& core)
{
    const int width = inputs.width, k = box.scale;
    const int columns = box.domaincolumns, rows = box.domainrows;
    const int finewidth = box.domainfinewidth(), fineheight = box.domainfineheight();
    const std::size_t count = static_cast<std::size_t>(finewidth) * static_cast<std::size_t>(fineheight);
    if (sea.size() != count || (!core.empty() && core.size() != count))
        throw std::invalid_argument("seapockets: the sea and the core must cover the box's fine cells");
    std::vector<uint8_t> filled(count, 0);
    if (reach <= 0)
        return filled;
    const Region region{reach, columns, rows, k, columns >= width};
    // The search never enters a global lake cell (lakes are the lake flood's, run after the coast).
    const GlobalLakes lakes = globallakes(inputs);
    auto lakecell = [&](int x, int y) {
        return !lakes.empty() && lakes.own[at(wrapcolumn(box.domaincolumn0 + x / k, width), box.domainrow0 + y / k, width)] != 0;
    };
    auto seablock = [&](int bx, int by) {
        return inputs.network.sea[at(wrapcolumn(box.domaincolumn0 + bx, width), box.domainrow0 + by, width)] != 0;
    };
    // A core cell: sea and core (without `core`: the sea of a global sea cell).
    auto iscore = [&](std::size_t cell, int bx, int by) {
        return sea[cell] != 0 && (core.empty() ? seablock(bx, by) : core[cell] != 0);
    };
    std::vector<int> mark(count, -1);  // the block whose search reached the cell
    std::vector<std::size_t> queue;
    for (int by = 0; by < rows; by++)
        for (int bx = 0; bx < columns; bx++)
        {
            // Only a block holding sea that is not core can hold a pocket; a lake cell's block holds none
            // (its fringe of sea meets the lake by S3.4's gate and is left alone).
            if (lakecell(bx * k, by * k))
                continue;
            bool any = false;
            for (int y = by * k; y < (by + 1) * k && !any; y++)
                for (int x = bx * k; x < (bx + 1) * k && !any; x++)
                    any = sea[at(x, y, finewidth)] != 0 && !iscore(at(x, y, finewidth), bx, by);
            if (!any)
                continue;
            const int id = by * columns + bx;
            // The seeds: the core sea within the reach.
            queue.clear();
            for (int dy = -reach; dy <= reach; dy++)
            {
                const int nby = by + dy;
                if (nby < 0 || nby >= rows)
                    continue;
                for (int dx = -reach; dx <= reach; dx++)
                {
                    int nbx = bx + dx;
                    if (region.wholeturn)
                        nbx = ((nbx % columns) + columns) % columns;
                    else if (nbx < 0 || nbx >= columns)
                        continue;
                    for (int y = nby * k; y < (nby + 1) * k; y++)
                        for (int x = nbx * k; x < (nbx + 1) * k; x++)
                        {
                            const std::size_t cell = at(x, y, finewidth);
                            if (mark[cell] != id && iscore(cell, nbx, nby))
                            {
                                mark[cell] = id;
                                queue.push_back(cell);
                            }
                        }
                }
            }
            // Through fine sea (4-connected), inside the blocks within the reach.
            spread(queue, mark, id, finewidth, fineheight, region, bx, by, [&](std::size_t cell) {
                return sea[cell] != 0 && !lakecell(static_cast<int>(cell % static_cast<std::size_t>(finewidth)),
                                                   static_cast<int>(cell / static_cast<std::size_t>(finewidth)));
            });
            for (int y = by * k; y < (by + 1) * k; y++)
                for (int x = bx * k; x < (bx + 1) * k; x++)
                {
                    const std::size_t cell = at(x, y, finewidth);
                    if (sea[cell] != 0 && mark[cell] != id)
                        filled[cell] = 1;
                }
        }
    return filled;
}

std::vector<uint8_t> pocketbasins(const RegionalGlobalInputs& inputs, const RegionalWindow& box, int reach,
                                  const std::vector<uint8_t>& sea, const std::vector<uint8_t>& filled,
                                  const std::vector<float>& surface, float level)
{
    const int width = inputs.width, k = box.scale;
    const int columns = box.domaincolumns, rows = box.domainrows;
    const int finewidth = box.domainfinewidth(), fineheight = box.domainfineheight();
    const std::size_t count = static_cast<std::size_t>(finewidth) * static_cast<std::size_t>(fineheight);
    if (sea.size() != count || filled.size() != count || surface.size() != count)
        throw std::invalid_argument("pocketbasins: the sea, the filled cells and the surface must cover the box's fine cells");
    std::vector<uint8_t> basin(count, 0);
    if (reach <= 0)
        return basin;
    const Region region{reach, columns, rows, k, columns >= width};
    // Never a global lake cell's fine cells, nor those of its 8 neighbours below its surface (the lake's
    // bed, which the lake flood fills after the coast).
    const GlobalLakes lakes = globallakes(inputs);
    auto lakebed = [&](std::size_t cell) {
        if (lakes.empty())
            return false;
        const int x = static_cast<int>(cell % static_cast<std::size_t>(finewidth));
        const int y = static_cast<int>(cell / static_cast<std::size_t>(finewidth));
        const std::size_t global = at(wrapcolumn(box.domaincolumn0 + x / k, width), box.domainrow0 + y / k, width);
        return lakes.own[global] != 0 || surface[cell] < lakes.surface[global];
    };
    std::vector<int> mark(count, -1);
    std::vector<std::size_t> queue;
    for (int by = 0; by < rows; by++)
        for (int bx = 0; bx < columns; bx++)
        {
            const int id = by * columns + bx;
            queue.clear();
            for (int y = by * k; y < (by + 1) * k; y++)
                for (int x = bx * k; x < (bx + 1) * k; x++)
                {
                    const std::size_t cell = at(x, y, finewidth);
                    if (filled[cell] != 0)
                    {
                        mark[cell] = id;
                        queue.push_back(cell);
                    }
                }
            if (queue.empty())
                continue;
            spread(queue, mark, id, finewidth, fineheight, region, bx, by,
                   [&](std::size_t cell) { return sea[cell] == 0 && (filled[cell] != 0 || (surface[cell] < level && !lakebed(cell))); });
            for (const std::size_t cell : queue)
                basin[cell] = 1;
        }
    return basin;
}

PocketCount countpockets(const std::vector<uint8_t>& filled, int domainwidth, int x0, int y0, int width, int height)
{
    PocketCount result;
    std::vector<uint8_t> seen(static_cast<std::size_t>(width) * static_cast<std::size_t>(height), 0);
    std::vector<std::size_t> queue;
    auto isfilled = [&](int x, int y) { return filled[at(x0 + x, y0 + y, domainwidth)] != 0; };
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
        {
            if (!isfilled(x, y))
                continue;
            result.pocketcells++;
            if (seen[at(x, y, width)] != 0)
                continue;
            result.pockets++;
            queue.assign(1, at(x, y, width));
            seen[at(x, y, width)] = 1;
            for (std::size_t next = 0; next < queue.size(); next++)
            {
                const int cx = static_cast<int>(queue[next] % static_cast<std::size_t>(width));
                const int cy = static_cast<int>(queue[next] / static_cast<std::size_t>(width));
                const int nx[4] = {cx - 1, cx + 1, cx, cx};
                const int ny[4] = {cy, cy, cy - 1, cy + 1};
                for (int n = 0; n < 4; n++)
                {
                    if (nx[n] < 0 || ny[n] < 0 || nx[n] >= width || ny[n] >= height || seen[at(nx[n], ny[n], width)] != 0 ||
                        !isfilled(nx[n], ny[n]))
                        continue;
                    seen[at(nx[n], ny[n], width)] = 1;
                    queue.push_back(at(nx[n], ny[n], width));
                }
            }
        }
    return result;
}

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
