#pragma once

#include "terrain.hpp"

// Shared implementation details; callers outside terrain use terrain.hpp.
namespace terrain::detail
{
void removefloatingmountains(planet& world);
int getcode(int dir);
void deleteoceanridge(planet& world, int x, int y, int dir);
void makemountainisland(planet& world, int x, int y, int peakheight);
}
