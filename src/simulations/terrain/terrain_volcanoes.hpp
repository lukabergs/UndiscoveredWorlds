#pragma once

// Isolated volcanoes and mountainous islands on the cell grid (docs/architecture/sphere-mesh.md
// sections 3 and 12): the (x, y) forms of terrain.hpp and terrain_internal.hpp forward here.

#include <vector>

#include "core/mesh/cell_grid.hpp"
#include "core/world/world_grid.hpp"

class planet;

namespace terrain
{
// A small mountainous island from cell c: one to three steps to a random neighbour, each step
// raising both cells, marking the cells within 1.5 legacy pixels of each (the legacy 3 x 3 box)
// as mountain islands and joining them by a ridge.
void makemountainisland(planet& world, mesh::CellId c, float peakheight);

// An isolated volcano at c with a line of extinct neighbours: each next one moves by the two
// direction fields (in [1, maximum elevation], centred at half of it; east and north) times
// six legacy pixels per half range, plus up to two pixels at random, through the cell's tangent
// frame. No submarine volcano within five legacy pixels of a continental shelf.
void createisolatedvolcano(planet& world, mesh::CellId c, const worldgrid<bool>& shelves,
                           const std::vector<float>& eastdirection, const std::vector<float>& northdirection,
                           int peakheight, bool strato);
}
