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
#include <memory>
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

// The lat-lon key grid of a level window's domain (mesh::latlonwindowgrid of the level's
// levelwidth() columns over the domain box, local id the domain's row-major index): the blocks and
// the parents' local ids read its keys only, so its radius is nominal (1 m).
std::shared_ptr<const mesh::CellGrid> levelkeygrid(const LevelWindow& window);
}
