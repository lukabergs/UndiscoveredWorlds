#pragma once

#include <vector>

class planet;
class boolshapetemplate;
class twointegers;

// Fractal fields shared by terrain, climate, and the editor.
void createfractal(std::vector<std::vector<int>>& arr, int awidth, int aheight, int grain, float valuemod, float valuemod2, int min, int max, bool extreme, bool wrapped);

// Template mountain chains and land shapes used by editing tools.
void createchains(planet& world, int baseheight, int conheight, std::vector<std::vector<int>>& fractal, std::vector<std::vector<int>>& plateaumap, boolshapetemplate landshape[], boolshapetemplate chainland[], twointegers focuspoints[], int focustotal, int focaldistance, int mode);
void fractaladdland(planet& world, std::vector<std::vector<int>>& fractal);
void makearchipelagos(planet& world, std::vector<std::vector<bool>>& removedland, boolshapetemplate landshape[]);

// Sea connectivity, coastal elevation, and island cleanup.
void removesmallseas(planet& world, int minseasize, int level);
void removestraights(planet& world);
void disruptseacoastline(planet& world, int centrex, int centrey, int avedepth, bool raise, int size);
void widenchannels(planet& world);
void loweroceans(planet& world);
void normalisecoasts(planet& world, int landheight, int seadepth, int severity);
void checkislands(planet& world);
void extendnoshade(planet& world);

// Ridge connectivity, mountain bases, and fjord relief.
void cleanmountainridges(planet& world);
int getridge(planet& world, int x, int y, int dir);
int getridge(std::vector<std::vector<int>>& arr, int x, int y, int dir);
int getoceanridge(planet& world, int x, int y, int dir);
void deleteridge(planet& world, int x, int y, int dir);
void deleteridge(planet& world, std::vector<std::vector<int>>& ridgesarr, std::vector<std::vector<int>>& heightsarr, int x, int y, int dir);
void raisemountainbases(planet& world, std::vector<std::vector<int>>& mountaindrainage, std::vector<std::vector<bool>>& OKmountains);
void addfjordmountains(planet& world);

// Elevation smoothing, relief noise, and depression filling.
void smoothland(planet& world, int amount);
void smoothonlysea(planet& world, int amount);
void createextraelev(planet& world);
void depressionfill(planet& world);
void addlandnoise(planet& world);

// Continental shelf masks and their connected gaps.
void makecontinentalshelves(planet& world, std::vector<std::vector<bool>>& shelves, int pointdist);

// Procedural ocean ridges, faults, and trenches.
void createoceanridges(planet& world, std::vector<std::vector<bool>>& shelves);
void createoceantrenches(planet& world, std::vector<std::vector<bool>>& shelves);

// Convert raw ridge masks from imports or FastLEM to mountain relief.
void createmountainsfromraw(planet& world, std::vector<std::vector<int>>& rawmountains, std::vector<std::vector<bool>>& OKmountains);

// Isolated volcanoes and mountainous islands.
void createisolatedvolcano(planet& world, int x, int y, std::vector<std::vector<bool>>& shelves, std::vector<std::vector<int>>& volcanodirection, int peakheight, bool strato);

// Crater placement and impact relief used by editing tools.
void createcratermap(planet& world, int cratertotal, std::vector<int>& squareroot, bool custom);
