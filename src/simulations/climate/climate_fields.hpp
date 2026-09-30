#pragma once

// World raster adapters, coastal temperature processing and the active all-land rainfall fallback.

#include <vector>

class planet;
class boolshapetemplate;

void createrainmap(planet& world, std::vector<std::vector<int>>& fractal,  int landtotal, int seatotal, boolshapetemplate smalllake[], boolshapetemplate shape[]);
void refreshadvectedrainfall(planet& world, std::vector<std::vector<int>>& fractal);
void createtemperaturemap(planet& world, std::vector<std::vector<int>>& fractal);
void createseaicemap(planet& world, std::vector<std::vector<int>>& fractal);
void createmountainprecipitation(planet& world);
void checkpoleclimates(planet& world);
