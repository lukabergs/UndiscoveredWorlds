#pragma once

// Helpers of the regional units on the cube's refinement levels
// (docs/architecture/sphere-mesh.md section 8, "Cube rules of the regional
// units"), all on 64-bit level keys at N_L cells per face edge, in closed form:
// - C, the cell chart: a face's index space extended across one face edge by
//   the net rotation of the cross net;
// - G, a 4 x 4 Catmull-Rom stencil in face coordinates whose nodes beyond one
//   face edge are ghosts at their true positions (1-D Catmull-Rom values along
//   the next face's line of centres) and beyond two edges bilinear samples;
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
    double sum = 0.0;
    for (std::size_t n = 0; n < stencil.key.size(); n++)
        sum += stencil.weight[n] * static_cast<double>(value(stencil.key[n]));
    return sum;
}

// samplebilinear's cube rule on keys: bilinear in (a, b) between the four
// centres around p on p's face, cells across an edge from the neighbour slots,
// the missing corner cell's weight shared by the three in proportion.
KeyStencil cubebilinearstencil(std::int64_t facecells, const Vec3& p);

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
