#pragma once

// Koppen and Holdridge classification from generated climate fields; numeric codes remain stable.

#include <string>

class climateworld;

void createclimatemap(climateworld& world);
void createbiomemap(climateworld& world);
short getclimate(climateworld& world, int x, int y);
int calculateholdridgebiome(float janTemp, float aprTemp, float julTemp, float octTemp, float janRain, float aprRain, float julRain, float octRain);
short calculateclimate(int elev, int sealevel, float wrain, float srain, float mintemp, float maxtemp);
std::string getclimatename(short climate);
std::string getclimatecode(short climate);
