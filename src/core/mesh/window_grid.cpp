#include "core/mesh/window_grid.hpp"

#include "core/mesh/cube_geometry.hpp"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <utility>

namespace mesh
{
namespace
{
using Terms = std::vector<std::pair<LevelKey, double>>;

int oppositeslot(int k) { return (k + 3) % slotcount + 1; }

// Catmull-Rom weights of the nodes -1, 0, 1, 2 at t in [0, 1).
void catmullrom(double t, double w[4])
{
    const double t2 = t * t, t3 = t2 * t;
    w[0] = 0.5 * (-t3 + 2.0 * t2 - t);
    w[1] = 0.5 * (3.0 * t3 - 5.0 * t2 + 2.0);
    w[2] = 0.5 * (-3.0 * t3 + 4.0 * t2 + t);
    w[3] = 0.5 * (t3 - t2);
}

// The continuous column and row indices of face coordinates (bilinearstencil's).
double columnindex(std::int64_t n, double a) { return a / cubestep(n) + 0.5 * static_cast<double>(n) - 0.5; }
double rowindex(std::int64_t n, double b) { return 0.5 * static_cast<double>(n) - b / cubestep(n) - 0.5; }

void addbilinear(Terms& terms, std::int64_t n, const Vec3& p, double scale)
{
    const KeyStencil stencil = cubebilinearstencil(n, p);
    for (std::size_t m = 0; m < stencil.key.size(); m++)
        terms.emplace_back(stencil.key[m], scale * stencil.weight[m]);
}

Vec3 extendedcentre(std::int64_t n, int face, std::int64_t column, std::int64_t row)
{
    return cubepoint(face, cubecentrea(n, column), cubecentreb(n, row));
}

// Ghosts of ghosts: a node of a ghost's 1-D stencil beyond the next face's edge
// is itself a ghost (level 1); a node beyond that (level 2) is a bilinear sample.
constexpr int ghostlevels = 2;

// A node of a stencil on face `face` at extended indices (i, j).
void addnode(Terms& terms, std::int64_t n, int face, std::int64_t i, std::int64_t j, double w, GhostRule rule, int level)
{
    const bool outi = i < 0 || i >= n, outj = j < 0 || j >= n;
    if (!outi && !outj)
    {
        terms.emplace_back(cubekey(n, face, i, j), w);
        return;
    }
    const Vec3 q = extendedcentre(n, face, i, j);
    const int edgei = i < 0 ? 3 : 1, edgej = j < 0 ? 0 : 2;
    const std::int64_t depthi = i < 0 ? -1 - i : i - n, depthj = j < 0 ? -1 - j : j - n;
    if (rule == GhostRule::indexunfolded)
    {
        // The measurement's alternative: the chart cell; beyond two edges bilinear.
        if (outi && outj)
            addbilinear(terms, n, q, w);
        else if (outi)
            terms.emplace_back(depthi < n ? cubeacrossedge(n, face, edgei, j, depthi) : nokey, w);
        else
            terms.emplace_back(depthj < n ? cubeacrossedge(n, face, edgej, i, depthj) : nokey, w);
        if (terms.back().first == nokey)
        {
            terms.pop_back();
            addbilinear(terms, n, q, w);
        }
        return;
    }
    if (level >= ghostlevels)
    {
        addbilinear(terms, n, q, w);
        return;
    }
    // The face that holds the ghost: across the crossed edge; beyond two edges
    // (a corner quadrant) the one of the two faces across that q lies on (the
    // third face at the corner).
    bool acrosscolumns = outi;
    if (outi && outj)
    {
        const int holder = cubeposition(q).face;
        if (cubeedgelink(face, edgei).face == holder)
            acrosscolumns = true;
        else if (cubeedgelink(face, edgej).face == holder)
            acrosscolumns = false;
        else
        {
            addbilinear(terms, n, q, w);
            return;
        }
    }
    const int edge = acrosscolumns ? edgei : edgej;
    const std::int64_t depth = acrosscolumns ? depthi : depthj;
    if (depth >= n)
    {
        addbilinear(terms, n, q, w);
        return;
    }
    // The ghost at its true position: a great circle of constant edge-normal
    // coordinate continues as a line of centres of the face across, parallel to
    // the edge at the same depth; the position along it is q's own coordinate
    // there, and the line's nodes are nodes of that face (ghosts again beyond it).
    const std::int64_t along = std::clamp(acrosscolumns ? j : i, std::int64_t{ 0 }, n - 1);  // only the line matters
    const CubeEdgeLink& link = cubeedgelink(face, edge);
    const CubeCell cell = cubecell(n, cubeacrossedge(n, face, edge, along, depth));
    const FacePosition position = cubepositiononface(link.face, q);
    const bool alongcolumns = link.edge == 0 || link.edge == 2;
    const double s = alongcolumns ? columnindex(n, position.a) : rowindex(n, position.b);
    const double s0 = std::floor(s);
    double ws[4];
    catmullrom(s - s0, ws);
    const std::int64_t first = static_cast<std::int64_t>(s0) - 1;
    for (int m = 0; m < 4; m++)
    {
        if (ws[m] == 0.0)
            continue;
        const std::int64_t t = first + m;
        addnode(terms, n, link.face, alongcolumns ? t : cell.column, alongcolumns ? cell.row : t, w * ws[m], rule, level + 1);
    }
}


KeyStencil merged(Terms& terms)
{
    std::stable_sort(terms.begin(), terms.end(), [](const auto& x, const auto& y) { return x.first < y.first; });
    KeyStencil stencil;
    for (const auto& term : terms)
    {
        if (!stencil.key.empty() && stencil.key.back() == term.first)
            stencil.weight.back() += term.second;
        else
        {
            stencil.key.push_back(term.first);
            stencil.weight.push_back(term.second);
        }
    }
    return stencil;
}
}

LevelKey cubechartkey(std::int64_t facecells, int face, std::int64_t column, std::int64_t row)
{
    const std::int64_t n = facecells;
    const bool outcolumn = column < 0 || column >= n, outrow = row < 0 || row >= n;
    if (!outcolumn && !outrow)
        return cubekey(n, face, column, row);
    if (outcolumn && outrow)
        return nokey;
    const std::int64_t depth = outcolumn ? (column < 0 ? -1 - column : column - n) : (row < 0 ? -1 - row : row - n);
    if (depth >= n)
        return nokey;
    if (outcolumn)
        return cubeacrossedge(n, face, column < 0 ? 3 : 1, row, depth);
    return cubeacrossedge(n, face, row < 0 ? 0 : 2, column, depth);
}

KeyStencil cubebilinearstencil(std::int64_t facecells, const Vec3& p)
{
    // bilinearstencil's cube branch on keys.
    const std::int64_t n = facecells;
    const FacePosition position = cubeposition(p);
    const double fi = columnindex(n, position.a), fj = rowindex(n, position.b);
    const std::int64_t i0 = static_cast<std::int64_t>(std::floor(fi)), j0 = static_cast<std::int64_t>(std::floor(fj));
    const double tx = fi - static_cast<double>(i0), ty = fj - static_cast<double>(j0);
    const std::int64_t is[4] = { i0, i0 + 1, i0, i0 + 1 }, js[4] = { j0, j0, j0 + 1, j0 + 1 };
    const double weight[4] = { (1.0 - tx) * (1.0 - ty), tx * (1.0 - ty), (1.0 - tx) * ty, tx * ty };
    Terms terms;
    double total = 0.0;
    for (int m = 0; m < 4; m++)
    {
        const std::int64_t ci = std::clamp(is[m], std::int64_t{ 0 }, n - 1), cj = std::clamp(js[m], std::int64_t{ 0 }, n - 1);
        const LevelKey base = cubekey(n, position.face, ci, cj);
        const bool outi = ci != is[m], outj = cj != js[m];
        int slot = 0;
        if (outi && outj)
            slot = is[m] > ci ? (js[m] < cj ? 2 : 4) : (js[m] > cj ? 6 : 8);
        else if (outi)
            slot = is[m] > ci ? 3 : 7;
        else if (outj)
            slot = js[m] < cj ? 1 : 5;
        const LevelKey key = slot == 0 ? base : cubekeyneighbour(n, base, slot);
        if (key == nokey)
            continue;
        terms.emplace_back(key, weight[m]);
        total += weight[m];
    }
    for (auto& term : terms)
        term.second = total > 0.0 ? term.second / total : 1.0 / static_cast<double>(terms.size());
    return merged(terms);
}

KeyStencil cubicstencil(std::int64_t facecells, int face, double a, double b, GhostRule rule)
{
    const std::int64_t n = facecells;
    const double u = columnindex(n, a), v = rowindex(n, b);
    const double u0 = std::floor(u), v0 = std::floor(v);
    double wx[4], wy[4];
    catmullrom(u - u0, wx);
    catmullrom(v - v0, wy);
    const std::int64_t i0 = static_cast<std::int64_t>(u0) - 1, j0 = static_cast<std::int64_t>(v0) - 1;
    Terms terms;
    for (int m = 0; m < 4; m++)
        for (int l = 0; l < 4; l++)
        {
            const double w = wx[l] * wy[m];
            if (w != 0.0)
                addnode(terms, n, face, i0 + l, j0 + m, w, rule, 0);
        }
    return merged(terms);
}

KeyStencil cubicstencil(std::int64_t facecells, const Vec3& p, GhostRule rule)
{
    const FacePosition position = cubeposition(p);
    return cubicstencil(facecells, position.face, position.a, position.b, rule);
}

WalkStep cubewalkstep(std::int64_t facecells, LevelKey key, int slot)
{
    const LevelKey next = cubekeyneighbour(facecells, key, slot);
    if (next == nokey)
        return {};
    return { next, oppositeslot(cubekeyreverse(facecells, key, slot)) };
}

std::vector<CellId> straightwalk(const CellGrid& grid, CellId start, int slot, int steps)
{
    std::vector<CellId> cells;
    CellId c = start;
    int k = slot;
    for (int s = 0; s < steps; s++)
    {
        const std::int32_t next = grid.neighbour(c, k);
        if (next < 0)
            break;
        const int r = grid.reverse(c, k);
        c = static_cast<CellId>(next);
        k = oppositeslot(r);
        cells.push_back(c);
    }
    return cells;
}

std::vector<LevelKey> dilatekeys(std::int64_t facecells, const std::vector<LevelKey>& keys, int rings)
{
    std::vector<LevelKey> set = keys;
    std::sort(set.begin(), set.end());
    set.erase(std::unique(set.begin(), set.end()), set.end());
    std::vector<LevelKey> frontier = set;
    for (int r = 0; r < rings && !frontier.empty(); r++)
    {
        std::vector<LevelKey> next;
        for (LevelKey key : frontier)
            for (int k = 1; k <= slotcount; k++)
            {
                const LevelKey other = cubekeyneighbour(facecells, key, k);
                if (other != nokey && !std::binary_search(set.begin(), set.end(), other))
                    next.push_back(other);
            }
        std::sort(next.begin(), next.end());
        next.erase(std::unique(next.begin(), next.end()), next.end());
        std::vector<LevelKey> grown;
        grown.reserve(set.size() + next.size());
        std::merge(set.begin(), set.end(), next.begin(), next.end(), std::back_inserter(grown));
        set = std::move(grown);
        frontier = std::move(next);
    }
    return set;
}

std::vector<int> ringdistances(const CellGrid& grid, const std::vector<CellId>& seeds, int rings)
{
    std::vector<int> ring(grid.cellcount(), -1);
    std::vector<CellId> frontier;
    for (CellId c : seeds)
        if (c < grid.cellcount() && ring[c] < 0)
        {
            ring[c] = 0;
            frontier.push_back(c);
        }
    for (int r = 1; r <= rings && !frontier.empty(); r++)
    {
        std::vector<CellId> next;
        for (CellId c : frontier)
            for (int k = 1; k <= slotcount; k++)
            {
                const std::int32_t other = grid.neighbour(c, k);
                if (other >= 0 && ring[static_cast<std::size_t>(other)] < 0)
                {
                    ring[static_cast<std::size_t>(other)] = r;
                    next.push_back(static_cast<CellId>(other));
                }
            }
        frontier = std::move(next);
    }
    return ring;
}
}
