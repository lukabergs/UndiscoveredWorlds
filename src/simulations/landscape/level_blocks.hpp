#pragma once

// The per-face quadtree of a window grid (docs/architecture/sphere-mesh.md section 8, "Cube rules
// of the regional units": block means are exact over the per-face quadtree; T8 W2-W2b): a level
// key's parent `factor` cells coarser and its factor^2 children, and a window's blocks (each
// parent with its children's local ids). On lat-lon the parent of (I, J) on the level's W_L x H_L
// grid is (I div factor, J div factor) on the (W_L / factor)-column grid and the children are the
// aligned factor x factor block; on the cube mesh::cubeparentkey and mesh::cubechildkeys (no
// cross-face parentage). Children are listed in ascending key order, which on lat-lon is the
// block's row-major order (a block never straddles the date line: W_L is a multiple of factor),
// so a block's sum in that order is today's bit for bit. Planet-free.

#include <cstddef>
#include <cstdint>
#include <memory>
#include <utility>
#include <vector>

#include "core/mesh/cell_grid.hpp"

namespace landscape
{
struct LevelWindow;  // refinement_grid.hpp

// The parent key of a level key of `grid`'s level (lat-lon: grid.columns() level columns; cube:
// grid.facecells() cells per face edge), factor >= 1 dividing the level's size.
mesh::LevelKey parentkey(const mesh::CellGrid& grid, mesh::LevelKey key, int factor);
// The factor^2 keys at `grid`'s level of a parent key factor cells coarser, ascending.
std::vector<mesh::LevelKey> childkeys(const mesh::CellGrid& grid, mesh::LevelKey parent, int factor);
// True when a key is its parent's first (smallest) child.
bool firstchild(const mesh::CellGrid& grid, mesh::LevelKey key, int factor);

// A window's blocks: every parent with a child in the window, ordered by its first child's local
// id (on a lat-lon window the domain's row-major parent order), with its factor^2 children's local
// ids ascending by key. Throws std::invalid_argument when a parent's child lies outside the window
// (a window is whole parents: the lat-lon domains and the cube windows of the halo chain are).
struct WindowBlocks
{
    int factor = 0;
    std::vector<mesh::LevelKey> parents;
    std::vector<mesh::CellId> children;  // factor^2 per parent, parent by parent

    std::size_t count() const { return parents.size(); }
    std::size_t perparent() const { return static_cast<std::size_t>(factor) * static_cast<std::size_t>(factor); }
    const mesh::CellId* childrenof(std::size_t block) const { return children.data() + block * perparent(); }
};
WindowBlocks windowblocks(const mesh::CellGrid& window, int factor);

// The mean of a block's children (summed in double in the children's order, over factor^2).
double blockmean(const WindowBlocks& blocks, std::size_t block, const std::vector<float>& field);

// R for a stencil node outside a domain on the cube (sphere-mesh.md section 8, decided 2026-10-19):
// a level key that is not one of the domain's keys takes the domain key at the smallest ring
// distance (8-connected hops on the level's closed-form neighbours, mesh::cubekeyneighbour), ties
// to the smallest key; a domain key is itself. Results are cached per key.
class NearestDomainKey
{
public:
    NearestDomainKey(std::int64_t facecells, std::vector<mesh::LevelKey> domain);
    bool contains(mesh::LevelKey key) const;
    // The index of a key's nearest domain key in the domain's ascending order.
    std::size_t operator()(mesh::LevelKey key) const;
    const std::vector<mesh::LevelKey>& keys() const { return domain_; }

private:
    std::int64_t facecells_ = 0;
    std::vector<mesh::LevelKey> domain_;
    mutable std::vector<std::pair<mesh::LevelKey, std::size_t>> cache_;  // sorted by key
};

// G of per-block values (values[b] for blocks.parents[b]) at every cell of a cube window, in its
// local order: the 4 x 4 Catmull-Rom stencil of the parent level (window.facecells() / factor cells
// per face edge; mesh::cubicstencil, true-position ghosts) at the cell's centre, a node outside the
// blocks' parents taking its nearest parent's value (NearestDomainKey), applied relative to the first
// node. The cube's form of the lat-lon block and parent-mean corrections' gather (whose nodes beyond
// the domain clamp to the nearest domain row and column). Throws std::invalid_argument on a lat-lon
// window or a size mismatch.
std::vector<double> gatherblocks(const WindowBlocks& blocks, const std::vector<double>& values, const mesh::CellGrid& window);

// gatherblocks' stencils precomputed for repeated passes over one window (T8 W2-W2f): per cell its
// terms (the block whose value a node takes, the node's weight) in the stencil's key order, applied
// with gatherblocks' arithmetic, so a planned gather is bit-identical to gatherblocks.
struct GatherPlan
{
    std::size_t blocks = 0;
    std::vector<std::uint32_t> start;  // cellcount + 1 offsets into block and weight
    std::vector<std::uint32_t> block;
    std::vector<double> weight;
};
GatherPlan gatherplan(const WindowBlocks& blocks, const mesh::CellGrid& window);
std::vector<double> gatherblocks(const GatherPlan& plan, const std::vector<double>& values);

// The lat-lon key grid of a level window's domain (mesh::latlonwindowgrid of the level's
// levelwidth() columns over the domain box, local id the domain's row-major index): the blocks and
// the parents' local ids read its keys only, so its radius is nominal (1 m).
std::shared_ptr<const mesh::CellGrid> levelkeygrid(const LevelWindow& window);
}
