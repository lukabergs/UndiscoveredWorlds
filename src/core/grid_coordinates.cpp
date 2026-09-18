#include "grid_coordinates.hpp"

int wrap(int value, int max)
{
    max++;

    value = value % max;

    if (value < 0)
        value = max + value;

    return value;
}

twointegers getdestination(int x, int y, int dir)
{
    if (dir == 8 || dir == 1 || dir == 2)
        y--;

    if (dir == 4 || dir == 5 || dir == 6)
        y++;

    if (dir == 2 || dir == 3 || dir == 4)
        x++;

    if (dir == 6 || dir == 7 || dir == 8)
        x--;

    twointegers dest;

    dest.x = x;
    dest.y = y;

    return dest;
}
