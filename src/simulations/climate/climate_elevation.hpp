#pragma once

// Lapse-rate adjustment of a temperature by the cell's elevation above sea level.
// World is any type with map(x, y), sealevel() and tempdecrease() (C/km):
// planet or climateworld. The integer overloads truncate as the original helpers
// did; the float overloads apply the lapse rate exactly.

#include <algorithm>

template<typename World>
int tempelevadd(const World& world, int temp, int i, int j)
{
    int elevation = world.map(i, j) - world.sealevel();

    if (elevation < 0)
        elevation = 0;

    float elevationadjust = (float)elevation;
    elevationadjust = elevationadjust / 1000.0f;
    elevationadjust = elevationadjust * world.tempdecrease();

    temp = temp - (int)elevationadjust;

    return temp;
}

// Removes the elevation element of a temperature.
template<typename World>
int tempelevremove(const World& world, int temp, int i, int j)
{
    int elevation = world.map(i, j) - world.sealevel();

    if (elevation < 0)
        elevation = 0;

    float elevationadjust = (float)elevation;
    elevationadjust = elevationadjust / 1000.0f;
    elevationadjust = elevationadjust * world.tempdecrease();

    temp = temp + (int)elevationadjust;

    return temp;
}

template<typename World>
float elevationlapse(const World& world, int i, int j)
{
    const int elevation = (std::max)(0, world.map(i, j) - world.sealevel());
    return static_cast<float>(elevation) / 1000.0f * world.tempdecrease();
}

template<typename World>
float tempelevadd(const World& world, float temp, int i, int j)
{
    return temp - elevationlapse(world, i, j);
}

template<typename World>
float tempelevremove(const World& world, float temp, int i, int j)
{
    return temp + elevationlapse(world, i, j);
}
