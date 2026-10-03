#pragma once

// Lapse-rate adjustment of a temperature by the cell's elevation above sea level.
// World is any type with map(x, y), sealevel() and tempdecrease() (C/km):
// planet or climateworld. Integer truncation matches the original helpers.

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
