#pragma once

// Lapse-rate adjustment of a temperature by the cell's elevation above sea level.
// World is any type with mapfloat(x, y), sealevel() and tempdecrease() (C/km): planet or
// climateworld. The lapse rate is applied exactly; the regional generator keeps its
// whole-number form in app/regional/region_queries.

#include <algorithm>

template<typename World>
float elevationlapse(const World& world, int i, int j)
{
    const float elevation = (std::max)(0.0f, world.mapfloat(i, j) - static_cast<float>(world.sealevel()));
    return elevation / 1000.0f * world.tempdecrease();
}

template<typename World>
float tempelevadd(const World& world, float temp, int i, int j)
{
    return temp - elevationlapse(world, i, j);
}

// Removes the elevation element of a temperature.
template<typename World>
float tempelevremove(const World& world, float temp, int i, int j)
{
    return temp + elevationlapse(world, i, j);
}
