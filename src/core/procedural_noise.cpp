#include "core/procedural_noise.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <vector>

#include "core/world/field_access.hpp"

#ifndef ARRAYWIDTH
#define ARRAYWIDTH 2048
#define ARRAYHEIGHT 1024
#endif

using namespace std;

// The shared generator is a global linear congruential stream. Each stage
// reseeds it from the world seed and a salt (deterministic_random.hpp).
static long g_seed = 1;

namespace
{
// Wraps an index into [0, max], where max is the final valid index.
int wrapindex(int value, int max)
{
    max++;
    value = value % max;
    if (value < 0)
        value = max + value;
    return value;
}
}

// This function finds the average of two values, wrapped around a maximum.

int wrappedaverage(int x, int y, int max)
{
    int result = -1;

    int diffa = abs(x - y);

    int diffb = abs((x - max) - y);

    int diffc = abs(x - (y - max));

    if (diffa <= diffb && diffa <= diffc)
        result = (x + y) / 2;
    else
    {
        if (diffb <= diffa && diffb <= diffc)
            result = ((x - max) + y) / 2;
        else
            result = (x + (y - max)) / 2;
    }

    if (result < 0)
        result = result + max;

    if (result > max)
        result = result - max;

    return result;
}

// This function takes two values (a and b) and creates a new value that is tilted towards a by percentage %.

int tilt(int a, int b, int percentage)
{
    int neg = 0;
    int tilt = 0;

    float diff = (float)(a - b);

    if (diff < 1.0f)
    {
        diff = 0.0f - diff;
        neg = 1;
    }

    if (diff != 0.0f)
        tilt = (int)((diff / 100.0f) * (float)percentage);

    if (neg == 1)
        tilt = 0 - tilt;

    int newval = b + tilt;

    return newval ;
}

// This does the same thing, but for floats.

float tilt(float a, float b, int percentage)
{
    int neg = 0;
    float tilt = 0.0f;

    float diff = a - b;

    if (diff < 1.0f)
    {
        diff = 0.0f - diff;
        neg = 1;
    }

    if (diff != 0.0f)
        tilt = (diff / 100.0f) * (float)percentage;

    if (neg == 1)
        tilt = 0.0f - tilt;

    float newval = b + tilt;

    return newval;
}

// A random number from a to b inclusive.
int random(int a, int b)
{
    int range = (b - a) + 1; // This is the range of possible numbers.

    int r = a + (fast_rand() % range);

    return r;
}

// Randomises the sign of an integer.
int randomsign(int a)
{
    if (random(0, 1) == 1)
        a = 0 - a;

    return a;
}

float randomsign(float a)
{
    if (random(0, 1) == 1)
        a = 0.0f - a;

    return a;
}

// These do the same with the C library generator.
int altrandom(int a, int b)
{
    int range = (b - a) + 1; // This is the range of possible numbers.

    UW_RANDOM_DRAW();
    int r = a + (rand() % range);

    return r;
}

int altrandomsign(int a)
{
    if (altrandom(0, 1) == 1)
        a = 0 - a;

    return (a);
}

// Warps a vector by the specified amount, creating a more natural look.
void warp(vector<vector<int>>& map, int width, int height, int maxelev, int warpfactor, bool vary)
{
    vector<vector<int>> dest(width+1, vector<int>(height+1, 0));
    
    int grain = 8; // Level of detail on this fractal map.
    float valuemod = 0.2f;
    int v = random(3, 6);
    float valuemod2 = (float)v;

    vector<vector<int>> warpright(width + 1, vector<int>(height + 1, 0));
    vector<vector<int>> warpdown(width + 1, vector<int>(height + 1, 0));

    createfractal(warpright, width, height, grain, valuemod, valuemod2, 1, maxelev, 0, 0);
    createfractal(warpdown, width, height, grain, valuemod, valuemod2, 1, maxelev, 0, 0);

    // Now create another one to vary the strength of this effect (if needed).

    vector<vector<int>> broad(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    grain = 4;
    valuemod = 0.002f;
    valuemod2 = 0.002f;

    if (vary)
        createfractal(broad, width, height, grain, valuemod, valuemod2, 1, maxelev, 0, 0);

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            int x = i;
            int y = j;

            float thiswarpright = (float)warpright[i][j] / (float)maxelev;
            float thiswarpdown = (float)warpdown[i][j] / (float)maxelev;

            thiswarpright = thiswarpright - 0.5f; // So they're in the range -0.5 - 0.5.
            thiswarpdown = thiswarpdown - 0.5f;

            thiswarpright = thiswarpright * warpfactor;
            thiswarpdown = thiswarpdown * warpfactor;

            if (vary)
            {
                float varyfactor = (float)broad[i][j] / (float)maxelev;

                thiswarpright = thiswarpright * varyfactor;
                thiswarpdown = thiswarpdown * varyfactor;
            }

            x = x + (int)thiswarpright;

            if (x<0 || x>width)
                x = wrapindex(x, width);

            y = y + (int)thiswarpdown;

            if (y < 0)
                y = 0;

            if (y > height)
                y = height;

            int val = map[x][y];

            dest[i][j] = (int)val;
        }
    }

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            map[i][j] = dest[i][j];
        }
    }
}

// A "monstrous" warp - http://thingonitsown.blogspot.com/2018/11/monstrous-terrain.html
void monstrouswarp(vector<vector<int>>& map, int width, int height, int maxelev, int warpfactor)
{
    vector<vector<int>> dest(width + 1, vector<int>(height + 1, 0));

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            float point = (float)map[i][j];

            point = point / (float)maxelev; // Makes it into a value from 0.0 to 1.0

            point = point * (float)warpfactor;

            int val = map[(int)point][(int)point];
 
            dest[i][j] = val;
        }
    }

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            map[i][j] = dest[i][j];
        }
    }
}

// Shifts everything in a vector to the left by a given number of pixels.
void shift(vector<vector<int>>& map, int width, int height, int offset)
{
    vector<vector<int>> dummy(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
            dummy[i][j] = map[i][j];
    }

    for (int i = 0; i <= width; i++)
    {
        int ii = i + offset;

        if (ii<0 || ii>width)
            ii = wrapindex(ii, width);

        for (int j = 0; j <= height; j++)
            map[i][j] = dummy[ii][j];
    }

    // For some reason this leaves a vertical line of noise, so we remove that manually now.

    int i = width - offset; // This is the x-coordinate of the problematic line.
    int iminus = i - 1;
    int iplus = i + 1;

    if (iminus < 0)
        iminus = width;

    if (iplus > width)
        iplus = 0;

    for (int j = 0; j <= height; j++) // Make the line the average of its neighbours to left and right.
    {
        map[i][j] = (map[iminus][j] + map[iplus][j]) / 2;
    }
}

// Flips an array vertically and/or horizontally.
void flip(vector<vector<int>>& arr, int awidth, int aheight, bool vert, bool horiz)
{
    vector<vector<int>> dummy(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    for (int i = 0; i <= awidth; i++)
    {
        for (int j = 0; j <= aheight; j++)
            dummy[i][j] = arr[i][j];
    }

    for (int i = 0; i <= awidth; i++)
    {
        for (int j = 0; j <= aheight; j++)
        {
            int ii, jj;

            if (horiz == 0)
                ii = 1;
            else
                ii = awidth - i;

            if (vert == 0)
                jj = j;
            else
                jj = aheight - j;

            arr[i][j] = dummy[ii][jj];
        }
    }
}

// Used to seed the generator.
// Linear congruential generator, from https://stackoverflow.com/questions/26237419/faster-than-rand
void fast_srand(long seed)
{
    g_seed = seed;
}

// Compute a pseudorandom integer.
// Output value in range [0, 32767]
int fast_rand(void)
{
    UW_RANDOM_DRAW();
    g_seed = (214013 * g_seed + 2531011);
    return (g_seed >> 16) & 0x7FFF;
}

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
                i = wrapindex(i, awidth);

            int ss = s / 2;
            int ii = i + ss;

            if (ii > awidth)
                ii = wrapindex(ii, awidth);

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
                i = wrapindex(i, awidth);

            int ss = s / 2;

            int ii = i + ss;

            if (ii > awidth)
                ii = wrapindex(ii, awidth);

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
                        ii = wrapindex(ii, awidth);

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
                coords[n][0] = wrapindex(coords[n][0], awidth);

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
                coords[n][0] = wrapindex(coords[n][0], awidth);

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
                coords[n][0] = wrapindex(coords[n][0], awidth);

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
                coords[n][0] = wrapindex(coords[n][0], awidth);

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
