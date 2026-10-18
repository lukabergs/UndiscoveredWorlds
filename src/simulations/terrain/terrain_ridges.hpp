#pragma once

// Ridge codes on the cell grid (docs/architecture/sphere-mesh.md, section 3): a cell's
// mountain_ridge and ocean_ridges hold one bit per neighbour slot (1 N .. 8 NW in the cell's
// own frame, bit slot - 1); a ridge joins a cell and its neighbour in slot k and is stored at
// both ends, at the neighbour in reverse(c, k). The (x, y) forms of terrain.hpp forward here.

#include "core/mesh/cell_grid.hpp"

class planet;

namespace terrain
{
// 1 when a ridge leaves cell c in slot (1..8).
int getridge(const planet& world, mesh::CellId c, int slot);
int getoceanridge(const planet& world, mesh::CellId c, int slot);
// Removes the ridge leaving c in slot and its other end at the neighbour (a cell left with no
// mountain ridge loses its mountain height); nothing where no ridge leaves in that slot.
void deleteridge(planet& world, mesh::CellId c, int slot);
void deleteoceanridge(planet& world, mesh::CellId c, int slot);
}
