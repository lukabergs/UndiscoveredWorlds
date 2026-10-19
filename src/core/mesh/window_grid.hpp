#pragma once

// Helpers of the regional units on the cube's refinement levels
// (docs/architecture/sphere-mesh.md section 8, "Cube rules of the regional
// units"), all on 64-bit level keys at N_L cells per face edge, in closed form:
// - C, the cell chart: a face's index space extended across one face edge by
//   the net rotation of the cross net;
// - G, a 4 x 4 Catmull-Rom stencil in face coordinates whose nodes beyond one
//   face edge are ghosts at their true positions (1-D Catmull-Rom values along
//   the next face's line of centres), beyond two edges (a corner quadrant) the
//   same along the third face that holds the node, the 1-D nodes beyond that
//   face's edge ghosts again, and only a ghost of a ghost's ghost bilinear;
// - S, straight walks: the next slot is the opposite of the reverse slot;
// - R, rings: dilation of a key set by 8-connected hops, and ring distances on
//   a grid's table.
// Stencils are lists of (key, weight) in ascending key order, so a value is a
// pure function of the keys and the field, in the same order in every window.

#include "core/mesh/cell_grid.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace mesh
{
// C: the cell at (column, row) of the face's index space extended across at
// most one face edge: inside the face the cell itself; beyond one edge by d
// cells (1 <= d <= N_L) the cell d - 1 cells into the face across, at the
// edge's index mapped by the net rotation; beyond two edges (a corner
// quadrant) or further nokey.
LevelKey cubechartkey(std::int64_t facecells, int face, std::int64_t column, std::int64_t row);

// A weighted sum of level cells: keys ascending and unique, weights summing to 1.
struct KeyStencil
{
    std::vector<LevelKey> key;
    std::vector<double> weight;
};

template <typename Value>
double applystencil(const KeyStencil& stencil, Value&& value)
{
    // Relative to the first value (as samplebilinear), so a constant field comes back exactly.
    if (stencil.key.empty())
        return 0.0;
    const double first = static_cast<double>(value(stencil.key[0]));
    double sum = 0.0;
    for (std::size_t n = 1; n < stencil.key.size(); n++)
        sum += stencil.weight[n] * (static_cast<double>(value(stencil.key[n])) - first);
    return first + sum;
}

// samplebilinear's cube rule on keys: bilinear in (a, b) between the four
// centres around p on p's face, cells across an edge from the neighbour slots,
// the missing corner cell's weight shared by the three in proportion.
KeyStencil cubebilinearstencil(std::int64_t facecells, const Vec3& p);
// The same terms without allocating (T8 W2-W3b; cubebilinearstencil is built from them): at most 4,
// keys ascending and unique.
struct SmallStencil
{
    int size = 0;
    LevelKey key[4] = {};
    double weight[4] = {};
};
SmallStencil cubebilinearterms(std::int64_t facecells, const Vec3& p);
template <typename Value>
double applystencil(const SmallStencil& stencil, Value&& value)
{
    // applystencil's arithmetic on a KeyStencil of the same terms.
    if (stencil.size == 0)
        return 0.0;
    const double first = static_cast<double>(value(stencil.key[0]));
    double sum = 0.0;
    for (int n = 1; n < stencil.size; n++)
        sum += stencil.weight[n] * (static_cast<double>(value(stencil.key[n])) - first);
    return first + sum;
}

// G: the 4 x 4 Catmull-Rom stencil at face coordinates (a, b) of a face, over
// the centres of the level's cells. trueposition is the rule of section 8;
// indexunfolded (a ghost takes the chart cell's value) exists to measure it.
enum class GhostRule
{
    trueposition,
    indexunfolded,
};
KeyStencil cubicstencil(std::int64_t facecells, int face, double a, double b, GhostRule rule = GhostRule::trueposition);
// The same on p's own face (cubeposition).
KeyStencil cubicstencil(std::int64_t facecells, const Vec3& p, GhostRule rule = GhostRule::trueposition);
// The keys one node of G (true-position ghosts) reads, at extended indices (column, row) of a face: the cell
// itself on the face, beyond an edge its ghost's expansion; ascending and unique (T8 W2-W3c). A point inside
// cell (i, j) reads the nodes i - 2 .. i + 2 by j - 2 .. j + 2 of its face, so their keys bound its stencil.
std::vector<LevelKey> cubicnodekeys(std::int64_t facecells, int face, std::int64_t column, std::int64_t row);

// G (true-position ghosts) of level `facecells` at every cell centre of a cube grid, without building a
// KeyStencil per call (T8 W2-W3b): a cell whose 4 x 4 nodes lie on its face is computed in place from
// its centre (cubicstencil's arithmetic); for the others the stencil's structure (its ghosts' expansions,
// their weights' factors, the merged keys and the order of their sums) depends only on the first node
// and which node weights are zero, so it is built once per such pattern and kept with the grid
// (centrestencils), and a cell's weights are its node weights times the factors, summed in that order.
// apply is applystencil(cubicstencil(facecells, grid.centre(c)), value) and visit(c, f) calls f(key,
// weight) over that stencil's terms in its order, bit for bit.
class CentreStencils
{
public:
    struct Interior
    {
        int face = 0;
        std::int64_t i0 = 0, j0 = 0;  // the first node
        double wx[4] = {}, wy[4] = {};
    };
    CentreStencils(const CellGrid& grid, std::int64_t facecells);
    std::int64_t facecells() const { return facecells_; }
    // The cell's first node and node weights; true for an on-face cell, false for a ghost pattern's.
    bool interior(CellId c, Interior& nodes) const;

    // A ghost pattern: per emitted term (cubicstencil's order before its merge) its node (row * 4 +
    // column) and the factors its node weight is multiplied by, left to right; the merged keys ascending,
    // each the sum of its terms in the stable key order.
    struct Pattern
    {
        std::vector<LevelKey> key;
        std::vector<std::uint8_t> node;
        std::vector<std::uint32_t> chain;  // terms + 1 offsets into factor
        std::vector<double> factor;
        std::vector<std::uint16_t> order;  // the terms in the stable key order
        std::vector<std::uint16_t> first;  // keys + 1 offsets into order
    };

    // cubicstencil's on-face terms in its order: f(key, weight) for the nonzero weights, rows outer.
    template <typename Visit>
    static void visitinterior(std::int64_t n, const Interior& nodes, Visit&& visitor)
    {
        const LevelKey base = static_cast<LevelKey>(nodes.face) * static_cast<LevelKey>(n) * static_cast<LevelKey>(n);
        for (int m = 0; m < 4; m++)
            for (int l = 0; l < 4; l++)
            {
                const double w = nodes.wx[l] * nodes.wy[m];
                if (w == 0.0)
                    continue;
                visitor(base + static_cast<LevelKey>(nodes.j0 + m) * static_cast<LevelKey>(n) + static_cast<LevelKey>(nodes.i0 + l), w);
            }
    }
    // applystencil's arithmetic over visited terms: relative to the first node.
    template <typename Terms, typename Value>
    static double applyterms(Terms&& terms, Value&& value)
    {
        bool started = false;
        double first = 0.0, sum = 0.0;
        terms([&](LevelKey key, double weight) {
            if (!started)
            {
                first = static_cast<double>(value(key));
                started = true;
                return;
            }
            sum += weight * (static_cast<double>(value(key)) - first);
        });
        return started ? first + sum : 0.0;
    }

    template <typename Visit>
    void visit(CellId c, Visit&& visitor) const
    {
        Interior nodes;
        if (interior(c, nodes))
        {
            visitinterior(facecells_, nodes, visitor);
            return;
        }
        const Pattern& pattern = patterns_[slot_[c]];
        const std::size_t terms = pattern.node.size();
        double stack[256];
        std::vector<double> heap;
        double* term = stack;
        if (terms > 256)
        {
            heap.resize(terms);
            term = heap.data();
        }
        for (std::size_t e = 0; e < terms; e++)
        {
            double w = nodes.wx[pattern.node[e] % 4] * nodes.wy[pattern.node[e] / 4];
            for (std::uint32_t f = pattern.chain[e]; f < pattern.chain[e + 1]; f++)
                w = w * pattern.factor[f];
            term[e] = w;
        }
        for (std::size_t s = 0; s + 1 < pattern.first.size(); s++)
        {
            double weight = term[pattern.order[pattern.first[s]]];
            for (std::uint32_t j = pattern.first[s] + 1u; j < pattern.first[s + 1]; j++)
                weight += term[pattern.order[j]];
            visitor(pattern.key[s], weight);
        }
    }

    template <typename Value>
    double apply(CellId c, Value&& value) const
    {
        return applyterms([&](auto&& visitor) { visit(c, visitor); }, value);
    }

private:
    const CellGrid* grid_ = nullptr;
    std::int64_t facecells_ = 0;
    static constexpr std::uint32_t onface = 0xffffffffu;
    std::vector<std::uint32_t> slot_;  // per cell: its pattern, or onface
    std::vector<Pattern> patterns_;
};
// The grid's CentreStencils at a level (facecells >= 1), built at the first call and kept with the grid.
const CentreStencils& centrestencils(const CellGrid& grid, std::int64_t facecells);

// G at one point (true-position ghosts), without a KeyStencil when its 4 x 4 nodes lie on p's face
// (T8 W2-W3b): apply is applystencil(cubicstencil(facecells, p), value) bit for bit.
class PointStencil
{
public:
    PointStencil(std::int64_t facecells, const Vec3& p);
    template <typename Value>
    double apply(Value&& value) const
    {
        if (onface_)
            return CentreStencils::applyterms([&](auto&& visitor) { CentreStencils::visitinterior(facecells_, nodes_, visitor); },
                                              value);
        return applystencil(stencil_, value);
    }

private:
    std::int64_t facecells_ = 0;
    bool onface_ = false;
    CentreStencils::Interior nodes_;
    KeyStencil stencil_;
};
// One node of that stencil before its Catmull-Rom weight, as a stencil of its own: column and row
// 0..3 of the 4 x 4 nodes (1 and 2 the central 2 x 2 around p), the cell itself on p's face or a
// ghost's 1-D Catmull-Rom at its true position beyond an edge (weights summing to 1); empty for an
// offset outside 0..3. Prolong's limit reads the central 2 x 2 (sphere-mesh.md section 8).
KeyStencil cubicstencilnode(std::int64_t facecells, const Vec3& p, int column, int row,
                            GhostRule rule = GhostRule::trueposition);

// The 2 r x 2 r footprint of a point (the regional pass's smooth sea gate, sphere-mesh.md section
// 8, T8 W2-W2c): the C-chart cells (cubechartkey) of columns i0 - r + 1 .. i0 + r and rows
// j0 - r + 1 .. j0 + r of p's face index space, i0 and j0 the cell centres at or before p (at r = 2
// the 4 x 4 nodes of G on the face and its charts); chart positions beyond two face edges (a corner
// quadrant, where no cell lies) or deeper than a face are left out. Ascending and unique; empty for
// reach < 1.
std::vector<LevelKey> cubefootprint(std::int64_t facecells, const Vec3& p, int reach);

// S: one straight step from a key in a slot: the neighbour and the slot that
// continues straight from it (the opposite of its reverse slot); nokey and 0
// where the slot has no neighbour (a corner diagonal).
struct WalkStep
{
    LevelKey key = nokey;
    int slot = 0;
};
WalkStep cubewalkstep(std::int64_t facecells, LevelKey key, int slot);
// The cells after start along a straight walk on a grid's table, at most steps
// of them, stopping where a slot has no neighbour (or leaves a window).
std::vector<CellId> straightwalk(const CellGrid& grid, CellId start, int slot, int steps);

// R: the keys within `rings` 8-connected hops of a key set (ascending, the set
// included), and on a grid the ring of every cell from a set of seeds (0 for a
// seed, -1 beyond `rings`).
std::vector<LevelKey> dilatekeys(std::int64_t facecells, const std::vector<LevelKey>& keys, int rings);
std::vector<int> ringdistances(const CellGrid& grid, const std::vector<CellId>& seeds, int rings);
}
