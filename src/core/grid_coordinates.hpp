#ifndef grid_coordinates_hpp
#define grid_coordinates_hpp

#include "classes.hpp"

// Wrap an integer into [0, max], where max is the final valid index.
int wrap(int value, int max);

// Direction 1 is north, continuing clockwise through 8; other values stay put.
// The returned coordinate is unbounded; callers choose wrapping or clipping.
twointegers getdestination(int x, int y, int dir);

#endif
