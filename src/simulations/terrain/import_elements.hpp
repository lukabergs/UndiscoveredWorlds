#pragma once

// Terrain elements a user adds to an imported map, one function per "Generate elements" button of the
// desktop application. They draw from the ambient random sequence and belong to no recipe.

#include <vector>

#include "core/classes.hpp"

class planet;

// Continental shelves and ocean trenches around the land of an imported map.
void generateimportedshelves(planet& world);
// Mid-ocean ridges on the deep sea floor.
void generateimportedoceanridges(planet& world);
// Fractal depths for the sea floor.
void generateimportedseabed(planet& world);
// Fractal elevation for the land, with canyon uplift.
void generateimportedlandelevation(planet& world);
// Mountain chains over the land (the Mountains and Hills buttons).
void generateimportedmountainchains(planet& world, boolshapetemplate landshape[], boolshapetemplate chainland[]);
// A random number of craters.
void generateimportedcraters(planet& world, std::vector<int>& squareroot);
