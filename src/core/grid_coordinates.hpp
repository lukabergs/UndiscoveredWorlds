#ifndef grid_coordinates_hpp
#define grid_coordinates_hpp

#include "classes.hpp"

// Wrap an integer into [0, max], where max is the final valid index.
inline int wrap(int value, int max)
{
    max++;

    value = value % max;

    if (value < 0)
        value = max + value;

    return value;
}

// Direction from (x, y) to (xx, yy): 1 is north, continuing clockwise through 8;
// 0 for the same cell.
inline int getdir(int x, int y, int xx, int yy)
{
    if (xx == x && yy < y)
        return (1);

    if (xx > x && yy < y)
        return (2);

    if (xx > x && yy == y)
        return (3);

    if (xx > x && yy > y)
        return (4);

    if (xx == x && yy > y)
        return (5);

    if (xx<x && yy>y)
        return (6);

    if (xx < x && yy == y)
        return (7);

    if (xx < x && yy < y)
        return (8);

    return (0);
}

// Direction 1 is north, continuing clockwise through 8; other values stay put.
// The returned coordinate is unbounded; callers choose wrapping or clipping.
twointegers getdestination(int x, int y, int dir);

#endif
