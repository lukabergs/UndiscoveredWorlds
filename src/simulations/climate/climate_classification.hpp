#pragma once

// Koppen and Holdridge classification from generated climate fields; numeric codes remain stable.

#include <string>

class planet;
class region;

void createclimatemap(planet& world);
void createbiomemap(planet& world);
short getclimate(planet& world, int x, int y);
short getclimate(region& region, int x, int y);
int calculateholdridgebiome(float janTemp, float aprTemp, float julTemp, float octTemp, float janRain, float aprRain, float julRain, float octRain);
short calculateclimate(int elev, int sealevel, float wrain, float srain, float mintemp, float maxtemp);
std::string getclimatename(short climate);
std::string getclimatecode(short climate);
