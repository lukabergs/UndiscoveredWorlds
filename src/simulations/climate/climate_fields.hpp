#pragma once

// World raster adapters, coastal temperature processing and the active all-land rainfall fallback.

#include <vector>

class climateworld;
class boolshapetemplate;

void createrainmap(climateworld& world, std::vector<std::vector<int>>& fractal,  int landtotal, int seatotal, boolshapetemplate smalllake[], boolshapetemplate shape[]);
void refreshadvectedrainfall(climateworld& world, std::vector<std::vector<int>>& fractal);
void createtemperaturemap(climateworld& world, std::vector<std::vector<int>>& fractal);
void createseaicemap(climateworld& world, std::vector<std::vector<int>>& fractal);
void createmountainprecipitation(climateworld& world);
void checkpoleclimates(climateworld& world);
