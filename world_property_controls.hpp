#pragma once

#include "planet.hpp"

struct WorldPropertyControls
{
    int size = 0;
    float gravity = 0.0f;
    float lunar = 0.0f;
    float eccentricity = 0.0f;
    int perihelion = 0;
    int rotation = 0;
    float tilt = 0.0f;
    float tempdecrease = 0.0f;
    int northpolaradjust = 0;
    int southpolaradjust = 0;
    int averagetemp = 0;
    float waterpickup = 0.0f;
    int glacialtemp = 0;
    bool rivers = true;
    bool lakes = true;
    bool deltas = true;
};

WorldPropertyControls makeworldpropertycontrols(const planet& world);
void syncworldpropertycontrols(const planet& world, WorldPropertyControls& controls);
void applyworldpropertycontrols(planet& world, const WorldPropertyControls& controls);
void clampworldpropertycontrols(WorldPropertyControls& controls);
float getdefaultgravityforsize(int size);
