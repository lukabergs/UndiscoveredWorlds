// Fractal fields shared by terrain, climate, and the editor.
// Algorithms by Jonathan Hill, extracted from globalterrain.cpp.

#include "terrain_internal.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <queue>

#include "classes.hpp"
#include "functions.hpp"
#include "planet.hpp"

using namespace std;

namespace
{
void newfractalinit(std::vector<std::vector<int>>& arr, int awidth, int aheight, int grain, int min, int max, bool extreme);
void newfractal(std::vector<std::vector<int>>& arr, int awidth, int aheight, int grain, float valuemod, float valuemod2, int min, int max, bool wrapped);
int square(std::vector<std::vector<int>>& arr, int awidth, int aheight, int s, int x, int y, int value, int min, int max, bool wrapped);
int diamond(std::vector<std::vector<int>>& arr, int awidth, int aheight, int s, int x, int y, int value, int min, int max, bool wrapped);
}

// This function creates a new fractal map.

void createfractal(vector<vector<int>>& arr, int awidth, int aheight, int grain, float valuemod, float valuemod2, int min, int max, bool extreme, bool wrapped)
{
    newfractalinit(arr, awidth, aheight, grain, min, max, extreme);

    newfractal(arr, awidth, aheight, grain, valuemod, valuemod2, min, max, wrapped);
}

namespace
{
// Initialises the random values for the fractal map generator.

void newfractalinit(vector<vector<int>>& arr, int awidth, int aheight, int grain, int min, int max, bool extreme)
{
    int s = (aheight + 1) / grain;

    int newval = 0;

    for (int i = 0; i <= awidth; i = i + s)
    {
        for (int j = 0; j <= aheight; j = j + s)
        {
            if (extreme == 1)
            {
                if (random(0, 1) == 1)
                    newval = min;

                else
                    newval = max;
            }
            else
                newval = random(min, max);

            arr[i][j] = newval;
        }
    }
}
}

namespace
{
// This is the main fractal generating routine.

void newfractal(vector<vector<int>>& arr, int awidth, int aheight, int grain, float valuemod, float valuemod2, int min, int max, bool wrapped)
{
    bool simple = 0;

    if (valuemod == valuemod2)
        simple = 1;

    int newheight = 0;

    int s = (aheight + 1) / grain;

    int focalx = awidth / 2; // Coordinates for the focal point of the value variance.
    int focaly = random(0, aheight);

    float valuediff = abs(valuemod - valuemod2); // difference between the two valuemods.
    float increment = valuediff / 100.0f;

    valuemod = valuemod * 100;
    valuemod2 = valuemod2 * 100;

    while (s != -1) // This loop will repeat, with s halving in size each time until it gets to 1.
    {
        int value = s * (int)valuemod; // This is the amount we will vary each new pixel by. The smaller the tile, the smaller this value should be.

        // First we go over the whole map doing the square step.

        for (int z = 0; z <= awidth + s; z = z + s)
        {
            int i = z;

            if (i > awidth)
                i = wrap(i, awidth);

            int ss = s / 2;
            int ii = i + ss;

            if (ii > awidth)
                ii = wrap(ii, awidth);

            for (int j = 0; j <= aheight; j = j + s)
            {
                int jj = j + ss;

                if (simple == 0)
                {
                    float focaldistance = (float)sqrt((focalx - ii) * (focalx - ii) + (focaly - jj) * (focaly - jj));
                    int perc = (int)(increment * focaldistance);

                    float newval = tilt(valuemod, valuemod2, perc);

                    newval = newval / 100.0f;
                    value = s * (int)newval;
                }

                newheight = square(arr, awidth, aheight, s, ii, jj, value, min, max, wrapped);

                if (jj <= aheight)
                    arr[ii][jj] = newheight;
            }
        }

        // Now we go over the whole map again, doing the diamond step.

        for (int z = 0; z <= awidth + s; z = z + s)
        {
            int i = z;

            if (i > awidth)
                i = wrap(i, awidth);

            int ss = s / 2;

            int ii = i + ss;

            if (ii > awidth)
                ii = wrap(ii, awidth);

            for (int j = 0; j <= aheight; j = j + s)
            {
                int jj = j + ss;

                if (simple == 0)
                {
                    float focaldistance = (float)sqrt((focalx - ii) * (focalx - ii) + (focaly - jj) * (focaly - jj));

                    int perc = (int)(increment * focaldistance);

                    float newval = tilt(valuemod, valuemod2, perc);

                    newval = newval / 100.0f;
                    value = s * (int)newval;
                }

                for (int n = 1; n <= 2; n++)
                {
                    int ii = 0;
                    int jj = 0;

                    if (n == 1)
                    {
                        ii = i + ss;
                        jj = j;
                    }

                    if (n == 2)
                    {
                        ii = i;
                        jj = j + ss;
                    }

                    if (jj >= 0 && jj <= aheight)
                    {
                        ii = wrap(ii, awidth);

                        newheight = diamond(arr, awidth, aheight, s, ii, jj, value, min, max, wrapped);

                        arr[ii][jj] = newheight;
                    }
                }
            }
        }

        if (s > 2) // Now halve the size of the tiles.
            s = s / 2;
        else
            s = -1;
    }

    if (simple == 0)
    {
        int offset = random(1, awidth);
        shift(arr, awidth, aheight, offset);
    }
}
}

namespace
{
// This does the square part of the fractal.

int square(vector<vector<int>>& arr, int awidth, int aheight, int s, int x, int y, int value, int min, int max, bool wrapped)
{
    value = value * 50;
    int dist = s / 2;
    int total = 0;
    int crount = 0;

    int coords[4][2];

    coords[0][0] = x - dist;
    coords[0][1] = y - dist;

    coords[1][0] = x + dist;
    coords[1][1] = y - dist;

    coords[2][0] = x + dist;
    coords[2][1] = y + dist;

    coords[3][0] = x - dist;
    coords[3][1] = y + dist;

    if (wrapped == 1)
    {
        int values[4];

        int nullval = 0 - max * 10;

        for (int n = 0; n < 4; n++)
            values[n] = nullval;

        for (int n = 0; n <= 3; n++)
        {
            if (coords[n][1] >= 0 && coords[n][1] <= aheight)
            {
                coords[n][0] = wrap(coords[n][0], awidth);

                int x = coords[n][0];
                int y = coords[n][1];

                values[n] = arr[x][y];
            }
        }

        int a = 0;

        if (values[0] != nullval && values[1] != nullval)
            a = wrappedaverage(values[0], values[1], max);
        else
        {
            if (values[0] == nullval)
                a = values[1];
            else
                a = values[0];
        }

        int b = 0;

        if (values[2] != nullval && values[3] != nullval)
            b = wrappedaverage(values[2], values[3], max);
        else
        {
            if (values[2] == nullval)
                b = values[3];
            else
                b = values[2];
        }

        if (a != nullval && b != nullval)
            total = wrappedaverage(a, b, max);
        else
        {
            if (a == nullval)
                total = b;
            else
                total = a;
        }
    }
    else
    {
        for (int n = 0; n <= 3; n++)
        {
            if (coords[n][1] >= 0 && coords[n][1] <= aheight)
            {
                coords[n][0] = wrap(coords[n][0], awidth);

                int x = coords[n][0];
                int y = coords[n][1];

                total = total + arr[x][y];
                crount++;
            }
        }

        total = total / crount;
    }

    int difference = randomsign(random(0, value));

    total = total + difference;

    if (wrapped == 1)
    {
        while (total > max)
            total = total - max;

        while (total < min)
            total = total + max;

    }
    else
    {
        if (total > max)
            total = max;

        if (total < min)
            total = min;
    }

    return (total);
}
}

namespace
{
// This does the diamond part of the fractal.

int diamond(vector<vector<int>>& arr, int awidth, int aheight, int s, int x, int y, int value, int min, int max, bool wrapped)
{
    value = value * 50;
    int dist = s / 2;
    int total = 0;
    int crount = 0;

    int coords[4][2];

    coords[0][0] = x;
    coords[0][1] = y - dist;

    coords[1][0] = x + dist;
    coords[1][1] = y;

    coords[2][0] = x;
    coords[2][1] = y + dist;

    coords[3][0] = x - dist;
    coords[3][1] = y;

    if (wrapped == 1)
    {
        int values[4];

        int nullval = 0 - max * 10;

        for (int n = 0; n < 4; n++)
            values[n] = nullval;

        for (int n = 0; n <= 3; n++)
        {
            if (coords[n][1] >= 0 && coords[n][1] <= aheight)
            {
                coords[n][0] = wrap(coords[n][0], awidth);

                int x = coords[n][0];
                int y = coords[n][1];

                values[n] = arr[x][y];
            }
        }

        int a = 0;

        if (values[0] != nullval && values[1] != nullval)
            a = wrappedaverage(values[0], values[1], max);
        else
        {
            if (values[0] == nullval)
                a = values[1];
            else
                a = values[0];
        }

        int b = 0;

        if (values[2] != nullval && values[3] != nullval)
            b = wrappedaverage(values[2], values[3], max);
        else
        {
            if (values[2] == nullval)
                b = values[3];
            else
                b = values[2];
        }

        if (a != nullval && b != nullval)
            total = wrappedaverage(a, b, max);
        else
        {
            if (a == nullval)
                total = b;
            else
                total = a;
        }
    }
    else
    {
        for (int n = 0; n <= 3; n++)
        {
            if (coords[n][1] >= 0 && coords[n][1] <= aheight)
            {
                coords[n][0] = wrap(coords[n][0], awidth);

                int x = coords[n][0];
                int y = coords[n][1];

                total = total + arr[x][y];
                crount++;
            }
        }

        total = total / crount;
    }

    int difference = randomsign(random(0, value));

    total = total + difference;

    if (wrapped == 1)
    {
        while (total > max)
            total = total - max;

        while (total < min)
            total = total + max;
    }
    else
    {
        if (total > max)
            total = max;

        if (total < min)
            total = min;
    }

    return (total);
}
}
