#include "simulations/landscape/regional_edges.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>

#include "core/deterministic_random.hpp"
#include "core/mesh/window_grid.hpp"
#include "simulations/landscape/level_blocks.hpp"
#include "simulations/landscape/regional_helpers.hpp"
#include "simulations/landscape/regional_surface.hpp"

// Erosion plan Phase R2: the contract of R6 (nearseacells and inlandcells moved here from
// regional_surface.cpp unchanged; seagate reproduces the 3 x 3 gate at edges.smoothgate 0)
// and R8's smooth gate and straight-run breaker (regional_edges.hpp states the rules).

namespace landscape
{
using regionaldetail::aroundcell;
using regionaldetail::at;
using regionaldetail::wrapcolumn;

namespace
{
// Salt of the straight-run breaker's draws (deterministiccombine(seed, salt, ...)).
constexpr std::int64_t straightsalt = 0x5354524e;

// A fine cell centre in global cell units: (I + 0.5) / k (upsampling.cpp's expression).
double globalunits(std::int64_t fine, int scale) { return (static_cast<double>(fine) + 0.5) / scale; }

int wrapindex(std::int64_t index, int width)
{
    const std::int64_t w = width;
    return static_cast<int>(((index % w) + w) % w);
}

// Whether a sea cell lies within `reach` global cells of (gx, gy) in global cell units in
// each axis (half-open: the cell centre's offset d in -reach <= d < reach): columns
// i0 - reach + 1 .. i0 + reach with i0 = floor(gx - 0.5), wrapped, rows likewise, clamped to
// the grid as bicubic clamps them. reach 2 is the 4 x 4 stencil `bicubic` reads.
bool seawithin(const std::vector<uint8_t>& sea, int width, int height, double gx, double gy, int reach)
{
    const std::int64_t i0 = static_cast<std::int64_t>(std::floor(gx - 0.5));
    const std::int64_t j0 = static_cast<std::int64_t>(std::floor(gy - 0.5));
    for (int b = 0; b < 2 * reach; b++)
    {
        const int row = static_cast<int>(std::clamp<std::int64_t>(j0 - (reach - 1) + b, 0, height - 1));
        for (int a = 0; a < 2 * reach; a++)
            if (sea[at(wrapindex(i0 - (reach - 1) + a, width), row, width)] != 0)
                return true;
    }
    return false;
}

bool lakecell(const RegionalGlobalInputs& inputs, std::size_t global)
{
    return inputs.lakesurface_m[global] > 0.0f && inputs.network.sea[global] == 0;
}
}

namespace
{
// The world's level-0 grid (the gates read its table; T8 W2-W2c).
const mesh::CellGrid& levelzero(const RegionalGlobalInputs& inputs, const char* name)
{
    if (inputs.grid == nullptr || inputs.grid->cellcount() != inputs.network.sea.size())
        throw std::invalid_argument(std::string(name) + ": the inputs need the world's grid, one cell per global value");
    return *inputs.grid;
}
}

std::vector<uint8_t> nearseacells(const RegionalGlobalInputs& inputs)
{
    // The 3 x 3 of every global cell on the world's table (sphere-mesh.md section 8: gates by table
    // neighbours; on lat-lon the raster 3 x 3, columns wrapping, rows stopping at the poles).
    const mesh::CellGrid& grid = levelzero(inputs, "nearseacells");
    std::vector<uint8_t> nearsea(grid.cellcount(), 0);
    for (std::size_t global = 0; global < nearsea.size(); global++)
    {
        const bool lake = inputs.lakesurface_m[global] > 0.0f && inputs.network.sea[global] == 0;
        nearsea[global] = !lake && aroundcell(grid, global, [&](std::size_t other) {
            return inputs.network.sea[other] != 0;
        }) ? 1 : 0;
    }
    return nearsea;
}

std::vector<float> inlandcells(const RegionalGlobalInputs& inputs, const RegionalParameters& parameters)
{
    const std::vector<uint8_t> nearsea = nearseacells(inputs);
    std::vector<float> inland(nearsea.size());
    const int reach = parameters.edges.smoothgate;
    if (reach <= 2)
    {
        for (std::size_t cell = 0; cell < nearsea.size(); cell++)
            inland[cell] = nearsea[cell] != 0 ? 0.0f : 1.0f;
        return inland;
    }
    // R8: no global sea cell within reach - 1 rings (R on the world's table: on lat-lon the
    // (2 reach - 1)-wide square, columns wrapping, rows beyond the poles not existing), or a lake cell.
    const mesh::CellGrid& grid = levelzero(inputs, "inlandcells");
    std::vector<mesh::CellId> seas;
    for (std::size_t global = 0; global < nearsea.size(); global++)
        if (inputs.network.sea[global] != 0)
            seas.push_back(static_cast<mesh::CellId>(global));
    const std::vector<int> ring = seas.empty() ? std::vector<int>(nearsea.size(), -1) : mesh::ringdistances(grid, seas, reach - 1);
    for (std::size_t global = 0; global < nearsea.size(); global++)
    {
        const bool sea = ring[global] >= 0;
        inland[global] = sea && !lakecell(inputs, global) ? 0.0f : 1.0f;
    }
    return inland;
}

std::vector<float> lakecellweight(const RegionalGlobalInputs& inputs, const RegionalWindow& box, const FineWarp& warp)
{
    std::vector<float> mask(inputs.lakesurface_m.size());
    for (std::size_t global = 0; global < mask.size(); global++)
        mask[global] = lakecell(inputs, global) ? 1.0f : 0.0f;
    return upsample(mask, inputs.width, inputs.height, box, warp);
}

std::vector<uint8_t> seagate(const RegionalGlobalInputs& inputs, const RegionalWindow& box, const FineWarp& warp,
                             const RegionalParameters& parameters, const std::vector<float>* lakeweight)
{
    const int width = inputs.width, height = inputs.height, k = box.scale;
    const int domainwidth = box.domainfinewidth(), domainheight = box.domainfineheight();
    const std::size_t count = static_cast<std::size_t>(domainwidth) * static_cast<std::size_t>(domainheight);
    if (lakeweight != nullptr && lakeweight->size() != count)
        throw std::invalid_argument("seagate: the lake weight must hold one value per fine cell of the box");
    const std::vector<uint8_t> nearsea = nearseacells(inputs);
    // S3.4: with a lake weight a lake cell's block is not excluded as a whole; a cell is excluded
    // where the lake cells dominate (weight >= 0.5) and a lake cell's block is gated by the rule of
    // the other cells (the 3 x 3 rule on its own global cell, or the smooth gate's reach).
    const auto lakeside = [&](std::size_t cell) { return lakeweight != nullptr && (*lakeweight)[cell] >= 0.5f; };
    std::vector<uint8_t> gate(count, 0);
    if (parameters.edges.smoothgate <= 0)
    {
        for (int y = 0; y < domainheight; y++)
        {
            for (int x = 0; x < domainwidth; x++)
            {
                const std::size_t cell = at(x, y, domainwidth);
                const std::size_t global = box.globalcell(x, y);
                if (lakeweight == nullptr)
                    gate[cell] = nearsea[global];
                else if (!lakeside(cell))
                    // The 3 x 3 of the ancestor on the world's table.
                    gate[cell] = aroundcell(levelzero(inputs, "seagate"), global,
                                            [&](std::size_t other) { return inputs.network.sea[other] != 0; }) ? 1 : 0;
            }
        }
        return gate;
    }

    // R8: the 3 x 3 rule, or a sea cell within edges.smoothgate global cells of the cell's
    // warped position (the Surfacer's faded warp), never in a lake cell's block.
    const int reach = parameters.edges.smoothgate;
    if (!warp.empty() && (warp.dx.size() != count || warp.dy.size() != count))
        throw std::invalid_argument("seagate: the warp must hold one displacement per fine cell of the box");
    std::vector<uint8_t> sea(nearsea.size(), 0);
    for (std::size_t global = 0; global < sea.size(); global++)
        sea[global] = inputs.network.sea[global] != 0 ? 1 : 0;
    const double sealevel = static_cast<double>(inputs.sealevel);
    const double fade_m = parameters.warp.fade_m;
    for (int y = 0; y < domainheight; y++)
    {
        const double gy = globalunits(box.absolutefinerow(y), k);
        for (int x = 0; x < domainwidth; x++)
        {
            const std::size_t cell = at(x, y, domainwidth);
            const std::size_t global = box.globalcell(x, y);
            if (lakeside(cell))
                continue;
            if (nearsea[global] != 0)
            {
                gate[cell] = 1;
                continue;
            }
            if (lakeweight == nullptr && lakecell(inputs, global))
                continue;
            const double gx = globalunits(box.absolutefinecolumn(x), k);
            double wx = gx, wy = gy;
            if (!warp.empty())
            {
                double fade = 1.0;
                if (fade_m > 0.0)
                {
                    const float lapsebase = static_cast<float>(bicubic(inputs.elevation_m, width, height, gx, gy));
                    fade = std::clamp(1.0 - std::fabs(static_cast<double>(lapsebase) - sealevel) / fade_m, 0.0, 1.0);
                }
                wx = gx + fade * warp.dx[cell];
                wy = gy + fade * warp.dy[cell];
            }
            gate[cell] = seawithin(sea, width, height, wx, wy, reach) ? 1 : 0;
        }
    }
    return gate;
}

void breakstraightruns(const RegionalGlobalInputs& inputs, const RegionalWindow& box, const RegionalParameters& parameters,
                       const std::vector<uint8_t>& gate, std::vector<uint8_t>& mask)
{
    const int straightmax = parameters.edges.straightmax;
    if (straightmax <= 0)
        return;
    const int width = inputs.width, k = box.scale;
    const int finewidth = box.domainfinewidth(), fineheight = box.domainfineheight();
    const std::size_t count = static_cast<std::size_t>(finewidth) * static_cast<std::size_t>(fineheight);
    if (gate.size() != count || mask.size() != count)
        throw std::invalid_argument("breakstraightruns: the gate and the mask must cover the box's fine cells");
    const int chunkmin = std::max(1, std::min(parameters.edges.chunkmin, parameters.edges.chunkmax));
    const int chunkmax = std::max(chunkmin, parameters.edges.chunkmax);
    const bool wholeturn = box.domaincolumns >= width;
    const int blockcells = k * k;
    const double share = parameters.coast.landshare;
    // The land-share rule's bound for a global sea cell's block (regional_surface.cpp's expression).
    const int mostland = share > 0.0 ? static_cast<int>(std::floor((1.0 - share) * blockcells + 1.0e-9)) : blockcells;

    // The runs are found on the mask as it came, the chunks written to `mask`.
    const std::vector<uint8_t> source = mask;
    // 1 sea, 0 land, -1 beyond the box or a pole (unknown: never land on a run's side).
    auto kind = [&](int x, int y) -> int {
        if (y < 0 || y >= fineheight)
            return -1;
        if (x < 0 || x >= finewidth)
        {
            if (!wholeturn)
                return -1;
            x = wrapindex(x, finewidth);
        }
        return source[at(x, y, finewidth)] != 0 ? 1 : 0;
    };

    // The four sides of a run (where its land lies): north, south, west, east; a run runs
    // along x for north and south, along y for west and east.
    const int sidex[4] = {0, 0, -1, 1};
    const int sidey[4] = {-1, 1, 0, 0};
    for (int by = 0; by < box.domainrows; by++)
        for (int bx = 0; bx < box.domaincolumns; bx++)
        {
            const std::size_t global = at(wrapcolumn(box.domaincolumn0 + bx, width), box.domainrow0 + by, width);
            if (lakecell(inputs, global))
                continue;
            const bool seablock = inputs.network.sea[global] != 0;
            const int x0 = bx * k, y0 = by * k;
            int land = 0;
            for (int y = y0; y < y0 + k; y++)
                for (int x = x0; x < x0 + k; x++)
                    land += mask[at(x, y, finewidth)] == 0 ? 1 : 0;
            auto inblock = [&](int x, int y) { return x >= x0 && x < x0 + k && y >= y0 && y < y0 + k; };
            // Writes one cell of a chunk inside the block: land always, sea only on a gated cell.
            auto flip = [&](int x, int y, bool toland) {
                if (!inblock(x, y))
                    return;
                const std::size_t cell = at(x, y, finewidth);
                if (toland)
                {
                    if (mask[cell] != 0)
                    {
                        mask[cell] = 0;
                        land++;
                    }
                }
                else if (mask[cell] == 0 && gate[cell] != 0)
                {
                    mask[cell] = 1;
                    land--;
                }
            };
            for (int side = 0; side < 4; side++)
            {
                const bool alongx = side < 2;
                const int sx = sidex[side], sy = sidey[side];         // towards its land
                for (int line = 0; line < k; line++)
                {
                    // The cells of this line inside the block, position t = 0 .. k - 1.
                    auto cellx = [&](int t) { return alongx ? x0 + t : x0 + line; };
                    auto celly = [&](int t) { return alongx ? y0 + line : y0 + t; };
                    auto runcell = [&](int t) {
                        const int x = cellx(t), y = celly(t);
                        return kind(x, y) == 1 && kind(x + sx, y + sy) == 0;
                    };
                    int t = 0;
                    while (t < k)
                    {
                        if (!runcell(t))
                        {
                            t++;
                            continue;
                        }
                        int end = t;
                        while (end + 1 < k && runcell(end + 1))
                            end++;
                        // The run t .. end; pieces longer than straightmax get chunks.
                        int p = t;
                        while (end - p + 1 > straightmax)
                        {
                            const int remaining = end - p + 1;
                            // The gap before the chunk (at least one run cell, so two chunks never
                            // join into a longer edge) from the piece's first cell, the chunk's size
                            // and kind from its own start cell (absolute fine positions).
                            const std::uint64_t gaphash = deterministiccombine(
                                inputs.seed, straightsalt, box.absolutefinecolumn(cellx(p)), box.absolutefinerow(celly(p)), side, 0);
                            const int gap = 1 + static_cast<int>(gaphash % static_cast<std::uint64_t>(std::min(straightmax, remaining - 1)));
                            const int start = p + gap;
                            const std::uint64_t hash = deterministiccombine(
                                inputs.seed, straightsalt, box.absolutefinecolumn(cellx(start)), box.absolutefinerow(celly(start)), side, 1);
                            const int size = chunkmin + static_cast<int>((hash >> 20) % static_cast<std::uint64_t>(chunkmax - chunkmin + 1));
                            // The chunk: a first row of `first` cells along the run and a second
                            // of `second` cells beyond it, centred (sizes 2 .. 7: 2+0, 2+1, 3+1,
                            // 4+1, 4+2, 5+2), so its own edges stay short.
                            const int second = size / 3;
                            const int first = std::min(size - second, end - start + 1);
                            const int secondoffset = std::max(0, (first - second) / 2);
                            const int secondcount = std::min(second, first - secondoffset);
                            // Land in a land cell's block; in a sea cell's, land or sea by the draw,
                            // land only while the block keeps its share, sea only on gated cells
                            // inside the block.
                            bool toland = !seablock || ((hash >> 40) & 1u) == 0;
                            if (toland && seablock && land + first + secondcount > mostland)
                                toland = false;
                            if (!toland)
                            {
                                bool possible = false;
                                for (int o = 0; o < first && !possible; o++)
                                {
                                    const int x = cellx(start + o) + sx, y = celly(start + o) + sy;
                                    possible = inblock(x, y) && gate[at(x, y, finewidth)] != 0;
                                }
                                if (!possible)
                                    toland = !seablock || land + first + secondcount <= mostland;
                                if (!possible && !toland)
                                {
                                    p = start + first;
                                    continue;
                                }
                            }
                            // Raised: the run's own cells and the row beyond them on the sea side;
                            // lowered: the land side's row and the one beyond it.
                            const int firstrow = toland ? 0 : 1, secondrow = toland ? -1 : 2;
                            for (int o = 0; o < first; o++)
                                flip(cellx(start + o) + firstrow * sx, celly(start + o) + firstrow * sy, toland);
                            for (int o = 0; o < secondcount; o++)
                                flip(cellx(start + secondoffset + o) + secondrow * sx,
                                     celly(start + secondoffset + o) + secondrow * sy, toland);
                            p = start + first;
                        }
                        t = end + 1;
                    }
                }
            }
        }
}

namespace
{
// A cube window's fine cells per level-0 cell edge and its level 0 (the inputs' cube grid).
int cubescale(const RegionalGlobalInputs& inputs, const mesh::CellGrid& window, const char* name)
{
    if (inputs.grid == nullptr || inputs.grid->kind() != mesh::GridKind::cube || window.kind() != mesh::GridKind::cube ||
        inputs.grid->facecells() < 1 || window.facecells() % inputs.grid->facecells() != 0)
        throw std::invalid_argument(std::string(name) + ": a cube window refining the inputs' cube level 0");
    return window.facecells() / inputs.grid->facecells();
}
}

std::vector<float> lakecellweight(const RegionalGlobalInputs& inputs, const mesh::CellGrid& window, const FineWarp& warp)
{
    cubescale(inputs, window, "lakecellweight");
    const mesh::CellGrid& grid = *inputs.grid;
    if (!warp.empty() && (warp.dx.size() != window.cellcount() || warp.dy.size() != window.cellcount()))
        throw std::invalid_argument("lakecellweight: the warp must hold one displacement per window cell");
    std::vector<float> mask(inputs.lakesurface_m.size());
    for (std::size_t global = 0; global < mask.size(); global++)
        mask[global] = lakecell(inputs, global) ? 1.0f : 0.0f;
    const int n = grid.facecells();
    std::vector<float> weight(window.cellcount());
    for (mesh::CellId c = 0; c < window.cellcount(); c++)
    {
        const mesh::Vec3 p = warp.empty() ? window.centre(c)
                                          : warpedpoint(window.centre(c), warp.dx[c], warp.dy[c], grid.radius_m(), n);
        weight[c] = static_cast<float>(mesh::applystencil(mesh::cubicstencil(n, p),
                                                          [&](mesh::LevelKey key) { return mask[grid.localofkey(key)]; }));
    }
    return weight;
}

std::vector<uint8_t> seagate(const RegionalGlobalInputs& inputs, const mesh::CellGrid& window, const FineWarp& warp,
                             const RegionalParameters& parameters, const std::vector<float>* lakeweight)
{
    const int k = cubescale(inputs, window, "seagate");
    const mesh::CellGrid& grid = *inputs.grid;
    const std::size_t count = window.cellcount();
    if (lakeweight != nullptr && lakeweight->size() != count)
        throw std::invalid_argument("seagate: the lake weight must hold one value per window cell");
    if (!warp.empty() && (warp.dx.size() != count || warp.dy.size() != count))
        throw std::invalid_argument("seagate: the warp must hold one displacement per window cell");
    const std::vector<uint8_t> nearsea = nearseacells(inputs);
    const auto lakeside = [&](std::size_t cell) { return lakeweight != nullptr && (*lakeweight)[cell] >= 0.5f; };
    const auto ancestor = [&](mesh::CellId c) {
        return static_cast<std::size_t>(grid.localofkey(parentkey(window, window.globalkey(c), k)));
    };
    std::vector<uint8_t> gate(count, 0);
    const int reach = parameters.edges.smoothgate;
    if (reach <= 0)
    {
        for (mesh::CellId c = 0; c < count; c++)
        {
            const std::size_t global = ancestor(c);
            if (lakeweight == nullptr)
                gate[c] = nearsea[global];
            else if (!lakeside(c))
                gate[c] = aroundcell(grid, global, [&](std::size_t other) { return inputs.network.sea[other] != 0; }) ? 1 : 0;
        }
        return gate;
    }
    const double sealevel = static_cast<double>(inputs.sealevel);
    const double fade_m = parameters.warp.fade_m;
    const int n = grid.facecells();
    const std::vector<float> lapsebase =
        !warp.empty() && fade_m > 0.0 ? upsample(inputs.elevation_m, grid, window) : std::vector<float>{};
    for (mesh::CellId c = 0; c < count; c++)
    {
        const std::size_t global = ancestor(c);
        if (lakeside(c))
            continue;
        if (nearsea[global] != 0)
        {
            gate[c] = 1;
            continue;
        }
        if (lakeweight == nullptr && lakecell(inputs, global))
            continue;
        mesh::Vec3 p = window.centre(c);
        if (!warp.empty())
        {
            double fade = 1.0;
            if (fade_m > 0.0)
                fade = std::clamp(1.0 - std::fabs(static_cast<double>(lapsebase[c]) - sealevel) / fade_m, 0.0, 1.0);
            p = warpedpoint(p, fade * warp.dx[c], fade * warp.dy[c], grid.radius_m(), n);
        }
        bool sea = false;
        for (const mesh::LevelKey key : mesh::cubefootprint(n, p, reach))
            sea = sea || inputs.network.sea[grid.localofkey(key)] != 0;
        gate[c] = sea ? 1 : 0;
    }
    return gate;
}
}
