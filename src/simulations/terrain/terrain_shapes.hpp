#pragma once

// Land shapes on the cell grid (docs/architecture/sphere-mesh.md sections 3, 11 and 12): the
// (x, y) form of makearchipelagos in terrain.hpp forwards here.

#include <cstdint>
#include <vector>

class planet;
class boolshapetemplate;

namespace terrain
{
// Scatters mountain islands over the sea cells where land was removed (removedland, one value
// per cell in id order), by a sphere noise (terrain.sphere_noise.archipelago_*); an island
// within two legacy pixels of other land is not marked as a mountain island.
void makearchipelagos(planet& world, const std::vector<std::uint8_t>& removedland, boolshapetemplate landshape[]);
}
