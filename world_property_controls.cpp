#include <algorithm>

#include "world_property_controls.hpp"

WorldPropertyControls makeworldpropertycontrols(const planet& world)
{
    WorldPropertyControls controls;
    syncworldpropertycontrols(world, controls);
    return controls;
}

void syncworldpropertycontrols(const planet& world, WorldPropertyControls& controls)
{
    controls.size = world.size();
    controls.gravity = world.gravity();
    controls.lunar = world.lunar();
    controls.eccentricity = world.eccentricity();
    controls.perihelion = world.perihelion();
    controls.rotation = world.rotation();
    controls.tilt = world.tilt();
    controls.tempdecrease = world.tempdecrease();
    controls.northpolaradjust = world.northpolaradjust();
    controls.southpolaradjust = world.southpolaradjust();
    controls.averagetemp = world.averagetemp();
    controls.waterpickup = world.waterpickup();
    controls.glacialtemp = world.glacialtemp();
}

void applyworldpropertycontrols(planet& world, const WorldPropertyControls& controls)
{
    world.setgravity(controls.gravity);
    world.setlunar(controls.lunar);
    world.seteccentricity(controls.eccentricity);
    world.setperihelion(controls.perihelion);
    world.setrotation((bool)controls.rotation);
    world.settilt(controls.tilt);
    world.settempdecrease(controls.tempdecrease);
    world.setnorthpolaradjust(controls.northpolaradjust);
    world.setsouthpolaradjust(controls.southpolaradjust);
    world.setaveragetemp(controls.averagetemp);
    world.setwaterpickup(controls.waterpickup);
    world.setglacialtemp(controls.glacialtemp);
}

void clampworldpropertycontrols(WorldPropertyControls& controls)
{
    controls.gravity = std::clamp(controls.gravity, 0.05f, 10.0f);
    controls.waterpickup = std::max(0.0f, controls.waterpickup);
    controls.tilt = std::clamp(controls.tilt, -90.0f, 90.0f);
    controls.eccentricity = std::clamp(controls.eccentricity, 0.0f, 0.999f);
    controls.lunar = std::clamp(controls.lunar, 0.0f, 10.0f);

    if (!controls.rivers)
    {
        controls.lakes = false;
        controls.deltas = false;
    }
}

float getdefaultgravityforsize(int size)
{
    if (size == 1)
        return 0.4f;

    if (size == 0)
        return 0.15f;

    return 1.0f;
}
