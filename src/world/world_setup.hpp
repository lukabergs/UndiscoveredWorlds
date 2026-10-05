#pragma once

// Initial state of a world and a region, default map appearance, land and sea totals.

#include <string>
#include <vector>

#include "core/classes.hpp"

class planet;
class region;

void initialiseworld(planet& world);
void initialisemapcolours(planet& world);
void initialisegradientmapappearance(planet& world);
void setdefaultnonreliefmapappearance(planet& world);
void initialiseregion(planet& world, region& region);
void changeworldproperties(planet& world);
void getlandandseatotals(planet& world);
