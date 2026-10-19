#include "simulations/landscape/upsampling.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>

#include "core/mesh/window_grid.hpp"
#include "simulations/landscape/level_blocks.hpp"

// Erosion plan step E2a (docs/architecture/erosion-hydrology-plan.md, Phase E, and
// section 6, "Faithful upsampling").

namespace landscape
{
namespace
{
int wrapcolumn(std::int64_t column, int width)
{
    return static_cast<int>(((column % width) + width) % width);
}

void checkfield(const std::vector<float>& field, int width, int height, const char* name)
{
    if (width < 1 || height < 1 || field.size() != static_cast<std::size_t>(width) * static_cast<std::size_t>(height))
        throw std::invalid_argument(std::string(name) + ": the field must hold width x height values of a positive grid");
}

void checkwindow(const RegionalWindow& window, int width, int height, const char* name)
{
    if (window.worldwidth != width || window.worldheight != height || window.scale < 1)
        throw std::invalid_argument(std::string(name) + ": the window belongs to another grid");
}

// A fine cell centre in global cell units: (I + 0.5) / k, a function of the
// absolute index only.
double globalunits(std::int64_t fine, int scale) { return (static_cast<double>(fine) + 0.5) / scale; }
}

double catmullrom(double p0, double p1, double p2, double p3, double t)
{
    return p1 + 0.5 * t * ((p2 - p0) + t * ((2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3) + t * (3.0 * (p1 - p2) + p3 - p0)));
}

double bicubic(const std::vector<float>& field, int width, int height, double gx, double gy)
{
    checkfield(field, width, height, "bicubic");
    // Cell (i, j)'s centre is at (i + 0.5, j + 0.5): the spline runs between the
    // centres i0 and i0 + 1 with i0 = floor(gx - 0.5).
    const double u = gx - 0.5, v = gy - 0.5;
    const double fu = std::floor(u), fv = std::floor(v);
    const double tx = u - fu, ty = v - fv;
    const std::int64_t i0 = static_cast<std::int64_t>(fu);
    const std::int64_t j0 = static_cast<std::int64_t>(fv);
    int columns[4];
    for (int a = 0; a < 4; a++)
        columns[a] = wrapcolumn(i0 - 1 + a, width);
    double rows[4];
    for (int b = 0; b < 4; b++)
    {
        const std::int64_t j = std::max<std::int64_t>(0, std::min<std::int64_t>(height - 1, j0 - 1 + b));
        const std::size_t base = static_cast<std::size_t>(j) * static_cast<std::size_t>(width);
        rows[b] = catmullrom(field[base + static_cast<std::size_t>(columns[0])], field[base + static_cast<std::size_t>(columns[1])],
                             field[base + static_cast<std::size_t>(columns[2])], field[base + static_cast<std::size_t>(columns[3])],
                             tx);
    }
    return catmullrom(rows[0], rows[1], rows[2], rows[3], ty);
}

std::vector<float> upsample(const std::vector<float>& field, int width, int height, const RegionalWindow& window)
{
    checkfield(field, width, height, "upsample");
    checkwindow(window, width, height, "upsample");
    const int finewidth = window.domainfinewidth(), fineheight = window.domainfineheight();
    std::vector<float> fine(static_cast<std::size_t>(finewidth) * static_cast<std::size_t>(fineheight));
    for (int y = 0; y < fineheight; y++)
    {
        const double gy = globalunits(window.absolutefinerow(y), window.scale);
        for (int x = 0; x < finewidth; x++)
        {
            const double gx = globalunits(window.absolutefinecolumn(x), window.scale);
            fine[static_cast<std::size_t>(y) * static_cast<std::size_t>(finewidth) + static_cast<std::size_t>(x)] =
                static_cast<float>(bicubic(field, width, height, gx, gy));
        }
    }
    return fine;
}

std::vector<float> upsample(const std::vector<float>& field, int width, int height, const RegionalWindow& window,
                            const FineWarp& warp)
{
    if (warp.empty())
        return upsample(field, width, height, window);
    checkfield(field, width, height, "upsample");
    checkwindow(window, width, height, "upsample");
    const int finewidth = window.domainfinewidth(), fineheight = window.domainfineheight();
    const std::size_t count = static_cast<std::size_t>(finewidth) * static_cast<std::size_t>(fineheight);
    if (warp.dx.size() != count || warp.dy.size() != count)
        throw std::invalid_argument("upsample: the warp must hold one displacement per domain fine cell");
    std::vector<float> fine(count);
    for (int y = 0; y < fineheight; y++)
    {
        const double gy = globalunits(window.absolutefinerow(y), window.scale);
        for (int x = 0; x < finewidth; x++)
        {
            const std::size_t cell = static_cast<std::size_t>(y) * static_cast<std::size_t>(finewidth) + static_cast<std::size_t>(x);
            const double gx = globalunits(window.absolutefinecolumn(x), window.scale);
            fine[cell] = static_cast<float>(bicubic(field, width, height, gx + warp.dx[cell], gy + warp.dy[cell]));
        }
    }
    return fine;
}

std::vector<float> upsample(const std::vector<float>& field, const mesh::CellGrid& levelzero, const mesh::CellGrid& window)
{
    if (field.size() != levelzero.cellcount() || levelzero.kind() != window.kind())
        throw std::invalid_argument("upsample: the field must hold one value per level-0 cell, the grids one kind");
    const std::size_t count = window.cellcount();
    std::vector<float> fine(count);
    if (window.kind() == mesh::GridKind::cube)
    {
        const int n = levelzero.facecells();
        if (n < 1 || window.facecells() % n != 0)
            throw std::invalid_argument("upsample: the window's level must refine level 0 by a whole factor");
        const auto value = [&](mesh::LevelKey key) {
            const mesh::CellId cell = levelzero.localofkey(key);
            if (cell == mesh::nocell)
                throw std::invalid_argument("upsample: a stencil node lies outside the level-0 grid");
            return field[cell];
        };
        // The stencils at the window's centres, kept with the window grid (T8 W2-W3b).
        const mesh::CentreStencils& stencils = mesh::centrestencils(window, n);
        for (mesh::CellId c = 0; c < count; c++)
            fine[c] = static_cast<float>(stencils.apply(c, value));
        return fine;
    }
    const int width = levelzero.columns(), height = levelzero.rows();
    if (levelzero.issubgrid() || width < 1 || window.columns() % width != 0)
        throw std::invalid_argument("upsample: level 0 must be the whole lat-lon grid, the window a whole refinement of it");
    const int k = window.columns() / width;
    const mesh::LevelKey columns = static_cast<mesh::LevelKey>(window.columns());
    for (mesh::CellId c = 0; c < count; c++)
    {
        const mesh::LevelKey key = window.globalkey(c);
        const double gx = globalunits(static_cast<std::int64_t>(key % columns), k);
        const double gy = globalunits(static_cast<std::int64_t>(key / columns), k);
        fine[c] = static_cast<float>(bicubic(field, width, height, gx, gy));
    }
    return fine;
}

void correctblockmeans(const std::vector<float>& field, int width, int height, const RegionalWindow& window, int passes,
                       std::vector<float>& fine)
{
    checkfield(field, width, height, "correctblockmeans");
    checkwindow(window, width, height, "correctblockmeans");
    const int k = window.scale;
    const int finewidth = window.domainfinewidth(), fineheight = window.domainfineheight();
    if (passes < 0 || fine.size() != static_cast<std::size_t>(finewidth) * static_cast<std::size_t>(fineheight))
        throw std::invalid_argument("correctblockmeans: passes must not be negative and fine must cover the domain");
    const int dcolumns = window.domaincolumns, drows = window.domainrows;
    // The domain's blocks over the per-face quadtree (level_blocks.hpp): each global cell of the domain
    // (its key, the level-0 id) with its k x k fine cells, in the domain's row-major order of global
    // cells, so block b is domain cell (b mod domaincolumns, b div domaincolumns). Keys only: the
    // window grid's radius is nominal.
    const WindowBlocks blocks = windowblocks(*regionalwindowgrid(window, 1.0), k);
    if (blocks.count() != static_cast<std::size_t>(dcolumns) * static_cast<std::size_t>(drows))
        throw std::logic_error("correctblockmeans: the domain's blocks are not its global cells");
    std::vector<double> residual(blocks.count());
    std::vector<float> residualfield(field.size());
    for (int pass = 0; pass < passes; pass++)
    {
        // The residual of every domain cell: the mean of its k x k block, summed in
        // ascending key order (the block's row-major order), minus the global value.
        for (std::size_t b = 0; b < blocks.count(); b++)
            residual[b] = blockmean(blocks, b, fine) - static_cast<double>(field[blocks.parents[b]]);

        // The residuals on the global grid for `bicubic`: the domain's cells, and the
        // two global cells beyond its edges that the stencil reaches, which take the
        // nearest domain cell's residual (a clamp, as at the poles). A cell of the
        // domain reached around the planet keeps its own residual; nothing else is read.
        const int rowfirst = std::max(0, window.domainrow0 - 2);
        const int rowlast = std::min(height - 1, window.domainrow0 + drows + 1);
        const bool wholeturn = dcolumns >= width;
        for (int row = rowfirst; row <= rowlast; row++)
        {
            const int dr = std::max(0, std::min(drows - 1, row - window.domainrow0));
            for (int rel = -2; rel < dcolumns + 2; rel++)
            {
                const bool outside = rel < 0 || rel >= dcolumns;
                if (outside && wholeturn)
                    continue;
                const int column = wrapcolumn(static_cast<std::int64_t>(window.domaincolumn0) + rel, width);
                if (outside && wrapcolumn(static_cast<std::int64_t>(column) - window.domaincolumn0, width) < dcolumns)
                    continue;
                const int dc = std::max(0, std::min(dcolumns - 1, rel));
                residualfield[static_cast<std::size_t>(row) * static_cast<std::size_t>(width) + static_cast<std::size_t>(column)] =
                    static_cast<float>(residual[static_cast<std::size_t>(dr) * static_cast<std::size_t>(dcolumns) + static_cast<std::size_t>(dc)]);
            }
        }

        // The interpolated residual subtracted at every fine centre.
        for (int y = 0; y < fineheight; y++)
        {
            const double gy = globalunits(window.absolutefinerow(y), k);
            for (int x = 0; x < finewidth; x++)
            {
                const double gx = globalunits(window.absolutefinecolumn(x), k);
                float& value = fine[static_cast<std::size_t>(y) * static_cast<std::size_t>(finewidth) + static_cast<std::size_t>(x)];
                value = static_cast<float>(static_cast<double>(value) - bicubic(residualfield, width, height, gx, gy));
            }
        }
    }
}

void correctblockmeans(const std::vector<float>& field, const mesh::CellGrid& levelzero, const mesh::CellGrid& window,
                       int passes, std::vector<float>& fine)
{
    if (levelzero.kind() != mesh::GridKind::cube || window.kind() != mesh::GridKind::cube || levelzero.issubgrid() ||
        levelzero.facecells() < 1 || window.facecells() % levelzero.facecells() != 0 || field.size() != levelzero.cellcount())
        throw std::invalid_argument("correctblockmeans: a cube level 0, its field and a cube window refining it");
    if (passes < 0 || fine.size() != window.cellcount())
        throw std::invalid_argument("correctblockmeans: passes must not be negative and fine must cover the window");
    const WindowBlocks blocks = windowblocks(window, window.facecells() / levelzero.facecells());
    std::vector<double> residual(blocks.count());
    const GatherPlan plan = passes > 0 ? gatherplan(blocks, window) : GatherPlan{};
    for (int pass = 0; pass < passes; pass++)
    {
        for (std::size_t b = 0; b < blocks.count(); b++)
            residual[b] = static_cast<double>(static_cast<float>(
                blockmean(blocks, b, fine) - static_cast<double>(field[levelzero.localofkey(blocks.parents[b])])));
        const std::vector<double> correction = gatherblocks(plan, residual);
        for (std::size_t cell = 0; cell < fine.size(); cell++)
            fine[cell] = static_cast<float>(static_cast<double>(fine[cell]) - correction[cell]);
    }
}
}
