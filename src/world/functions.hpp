//
//  functions.hpp
//  Undiscovered Worlds
//
//  Created by Jonathan Hill on 24/07/2019.
//
//  The code for this project is released under the GNU General Public Licence v3.0f - https://choosealicense.com/licenses/gpl-3.0/
//
//  Please note that this code requires the following libraries to work:
// 
//  SFML - https://www.sfml-dev.org/
//  Dear ImGui - https://github.com/ocornut/imgui
//  ImGui-SFML - https://github.com/eliasdaler/imgui-sfml
//  ImGuiFileDialog - https://github.com/aiekick/ImGuiFileDialog
//
//  main.cpp contains the functions for drawing the map images (as well as the main loop, of course), and all functions that require external libraries.
//  The former misc.cpp is divided by owner: core/text_values, core/grid_shapes, world/world_setup, world/region_queries, simulations/terrain/coast_queries, simulations/hydrology/flow_queries and app/appearance_files.
//  Simulation modules under src/simulations contain terrain, climate, and hydrology.
//  regionalmap.cpp contains the functions used to generate the regional terrain.
//
//  The code uses a number of shape templates; their data is assets/terrain/shape-templates.txt, decoded in simulations/terrain/shape_templates.cpp.
//
//  The global terrain/climate information is stored in an object of the planet class. The regional terrain/climate information is stored in an object of the region class. Both of these classes are defined in the relevant hpp/cpp files.
//
//  One oddity to be aware of: variables that hold the width/height of objects typically hold the index number of the final item, rather than the actual size. E.g. if planet.width() is set to 100, that means the width has 101 elements, with the index of the last one being 100. Hence the frequent use of loops such as for (int i=0; i<=width; i++). I'm aware that this is eccentric but I find it more intuitive.


#ifndef functions_hpp
#define functions_hpp

#include <stdio.h>
#include <stdint.h>
#include <cstdint>
#include <algorithm>
#include <array>
#include <functional>
#include <thread>
#include <vector>

#include "simulations/climate/climate_diagnostics.hpp"
#include "simulations/climate/climate_elevation.hpp"

#include "pipeline/world_generation_debug.hpp"

#include "core/classes.hpp"
#include "world/planet.hpp"
#include "world/region.hpp"

#define ARRAYWIDTH 2048
#define ARRAYHEIGHT 1024

#define REGIONALTILEWIDTH 32
#define REGIONALTILEHEIGHT 32

#define MOUNTAINTEMPLATESTOTAL 27

#define MAXCRATERRADIUS 50

using namespace std;

// Define some enums.

enum screenmodeenum { quit, createworldscreen, creatingworldscreen, globalmapscreen, regionalmapscreen, generatingregionscreen, importscreen, completingimportscreen, movingtoglobalmapscreen, exportareascreen, exportingareascreen, loadingworldscreen, savingworldscreen, generatingterrainscreen, loadfailure, settingsloadfailure };

struct ImportedClimateMaps;

constexpr std::array<const char*, BIOMEMAPCOLOURCOUNT> biomemapnames =
{ {
    "Ice",
    "Polar desert",
    "Polar dry tundra",
    "Polar moist tundra",
    "Polar wet tundra",
    "Polar rain tundra",
    "Subpolar desert",
    "Subpolar dry tundra",
    "Subpolar moist tundra",
    "Subpolar wet tundra",
    "Subpolar rain tundra",
    "Boreal desert",
    "Boreal dry bush",
    "Boreal moist forest",
    "Boreal wet forest",
    "Boreal rain forest",
    "Cool temperate desert",
    "Cool temperate desert bush",
    "Cool temperate steppe",
    "Cool temperate moist forest",
    "Cool temperate wet forest",
    "Cool temperate rain forest",
    "Warm temperate desert",
    "Warm temperate desert bush",
    "Warm temperate thorn steppe",
    "Warm temperate dry forest",
    "Warm temperate moist forest",
    "Warm temperate wet forest",
    "Warm temperate rain forest",
    "Subtropical desert",
    "Subtropical desert bush",
    "Subtropical thorn steppe",
    "Subtropical dry forest",
    "Subtropical moist forest",
    "Subtropical wet forest",
    "Subtropical rain forest",
    "Tropical desert",
    "Tropical desert bush",
    "Tropical thorn steppe",
    "Tropical very dry forest",
    "Tropical dry forest",
    "Tropical moist forest",
    "Tropical wet forest",
    "Tropical rain forest",
} };

constexpr std::array<std::array<int, 3>, BIOMEMAPCOLOURCOUNT> defaultbiomemapcolours =
{ {
    { 65, 171, 155 },
    { 173, 68, 83 },
    { 158, 51, 146 },
    { 61, 125, 179 },
    { 42, 199, 102 },
    { 103, 66, 166 },
    { 163, 60, 137 },
    { 165, 121, 109 },
    { 127, 86, 116 },
    { 39, 182, 105 },
    { 71, 111, 175 },
    { 154, 52, 191 },
    { 173, 191, 73 },
    { 194, 47, 54 },
    { 36, 166, 108 },
    { 40, 157, 184 },
    { 55, 80, 153 },
    { 199, 62, 130 },
    { 36, 41, 173 },
    { 176, 122, 56 },
    { 91, 153, 58 },
    { 66, 196, 49 },
    { 171, 36, 88 },
    { 196, 65, 152 },
    { 168, 39, 58 },
    { 59, 50, 161 },
    { 130, 59, 156 },
    { 156, 59, 125 },
    { 120, 153, 44 },
    { 194, 161, 43 },
    { 161, 68, 50 },
    { 60, 199, 155 },
    { 135, 199, 58 },
    { 201, 60, 190 },
    { 66, 197, 201 },
    { 61, 85, 191 },
    { 161, 152, 56 },
    { 56, 161, 37 },
    { 101, 67, 196 },
    { 189, 78, 42 },
    { 176, 96, 42 },
    { 101, 36, 158 },
    { 48, 201, 133 },
    { 54, 143, 162 },
} };


#include "core/parallel_rows.hpp"
#include "core/grid_coordinates.hpp"

#include "core/deterministic_random.hpp"
#include "core/procedural_noise.hpp"


#include "pipeline/generation_progress.hpp"

// The former misc.cpp: each group of functions is declared with its source.
#include "core/text_values.hpp"
#include "core/grid_shapes.hpp"
#include "world/world_setup.hpp"
#include "simulations/terrain/coast_queries.hpp"
#include "simulations/hydrology/flow_queries.hpp"
#include "world/region_queries.hpp"

// Simulation entry points live with their implementations. This header
// remains a compatibility umbrella for the desktop and regional generator.
#include "simulations/terrain/terrain.hpp"
#include "simulations/geology/world_tectonics.hpp"

#include "pipeline/climate_generation.hpp"
#include "simulations/climate/climate_fields.hpp"
#include "simulations/climate/climate_classification.hpp"
#include "simulations/hydrology/drainage.hpp"
#include "simulations/hydrology/lakes.hpp"
#include "simulations/hydrology/lake_effects.hpp"
#include "simulations/hydrology/deltas.hpp"
#include "simulations/hydrology/wetlands.hpp"
#include "simulations/hydrology/tides.hpp"
#include "simulations/terrain/climate_landforms.hpp"

// Coupled ocean/atmosphere and validation integration.
#include "simulations/climate/climate_coupling.hpp"

// Declare functions that are in regionalmap.cpp

void removesealakes(planet& world, region& region, int leftx, int lefty, int rightx, int righty);

void generateregionalmap(planet& world, region& region, boolshapetemplate smalllake[], boolshapetemplate island[], peaktemplate& peaks, vector<vector<float>>& riftblob, int riftblobsize, int partial, byteshapetemplate smudge[], byteshapetemplate smallsmudge[], vector<int>& squareroot);
void makeregionalwater(planet& world, region& region, vector<vector<bool>>& safesaltlakes, vector<vector<bool>>& disruptpoints, vector<vector<int>>& rriverscarved, vector<vector<int>>& fakesourcex, vector<vector<int>>& fakesourcey, boolshapetemplate smalllake[], boolshapetemplate island[], vector<vector<float>>& riftblob, int riftblobsize, int xleft, int xright, int ytop, int ybottom);
void makeregionalterrain(planet& world, region& region, vector<vector<bool>>& disruptpoints, vector<vector<bool>>& riverinlets, vector<vector<bool>>& globalestuaries, vector<vector<int>>& rriverscarved, vector<vector<int>>& fakesourcex, vector<vector<int>>& fakesourcey, boolshapetemplate smalllake[], peaktemplate& peaks, byteshapetemplate smallsmudge[], int xleft, int xright, int ytop, int ybottom, vector<int>& squareroot, vector<vector<int>>& warptox, vector<vector<int>>& warptoy);
void makeregionalunderseaterrain(planet& world, region& region, peaktemplate& peaks, byteshapetemplate smudge[], byteshapetemplate smallsmudge[], int xleft, int xright, int ytop, int ybottom, vector<vector<int>>& warptox, vector<vector<int>>& warptoy);
void makeregionalmiscellanies(planet& world, region& region, vector<vector<bool>>& safesaltlakes, vector<vector<bool>>& riverinlets, vector<vector<bool>>& globalestuaries, boolshapetemplate smalllake[], byteshapetemplate smallsmudge[], int xleft, int xright, int ytop, int ybottom);
void makeregionalclimates(planet& world, region& region, vector<vector<bool>>& safesaltlakes, boolshapetemplate smalllake[], int xleft, int xright, int ytop, int ybottom);
twointegers convertregionaltoglobal(planet& world, region& region, int x, int y);
void makerivertile(planet& world, region& region, int dx, int dy, int sx, int sy, vector<vector<int>>& rriverscarved, boolshapetemplate smalllake[], vector<vector<bool>>& rivercurves);
twointegers getjunctionpoint(planet& world, region& region, int dx, int dy, int sx, int sy, bool lakepresent, bool goingtolake, bool diag, int maininflow, int inx, int iny, int outx, int outy);
void makedeltatile(planet& world, region& region, int dx, int dy, int sx, int sy, vector<vector<int>>& rriverscarved);
threeintegers calculateregionalriver(planet& world, region& region, int dx, int dy, int sx, int sy, twofloats pt, twofloats mm1, twofloats mm2, twofloats mm3, int janinflow, int julinflow, int riverlength, int isittributary, int goingtolake, bool delta);
void removeweirdelevations(planet& world, region& region, int dx, int dy, int sx, int sy);
void removelowflow(region& region, int dx, int dy);
void addtoexistingregionalriver(planet& world, region& region, int dx, int dy, int sx, int sy, int x, int y, int janinflow, int julinflow, bool delta);
void carveriver(planet& world, region& region, int dx, int dy, int sx, int sy, int riverlength, int riverstartx, int riverstarty, int riverendx, int riverendy, int startlandlevel, int endlandlevel, bool goingtolake, vector<vector<int>>& rriverscarved);
void carverivertributary(planet& world, region& region, int dx, int dy, int sx, int sy, int thisriverlength, int thisriverstartx, int thisriverstarty, int thisriverendx, int thisriverendy, int thisstartlandlevel, int thisendlandlevel, int riverendx, int riverendy, int endlandlevel, vector<vector<int>>& rriverscarved);
twointegers gettributaryendpoint(region& region, int startx, int starty);
void addsprings(planet& world, region& region, int sx, int sy, int dx, int dy, int junctionpointx, int junctionpointy, int riverendx, int riverendy, int maxmidvar, vector<vector<int>>& rriverscarved);
twointegers nearestlowerriver(region& region, int dx, int dy, int i, int j, int landlevel);
bool riverdircheck(region& region, int x, int y, int x2, int y2, int delta);
int removeriverorphans(planet& world, region& region, int sx, int sy, int dx, int dy, int riverlength);
void removeriverwidows(planet& world, region& region, int sx, int sy, int dx, int dy, int junctionpointx, int junctionpointy);
void removeregionalriverstraights(region& region, int dx, int dy, bool delta, vector<vector<int>>& rriverscarved);
int findriverlength(region& region, int startx, int starty, int endx, int endy, bool delta);
void joinuprivers(region& region, int dx, int dy, int sx, int sy);
void removenegativeflow(region& region, int dx, int dy);
void createsmalllake(planet& world, region& region, int dx, int dy, int sx, int sy, int centrex, int centrey, vector<vector<int>>& mainriver, boolshapetemplate smalllake[]);
void createsmallsaltlake(planet& world, region& region, int dx, int dy, int sx, int sy, int centrex, int centrey, int surfacelevel, vector<vector<bool>>& safesaltlakes, boolshapetemplate smalllake[]);
void createsmallglaciallake(planet& world, region& region, int dx, int dy, int sx, int sy, int centrex, int centrey, vector<vector<int>>& mainriver);
void expandrivers(planet& world, region& region, int dx, int dy, int sx, int sy, bool delta, vector<vector<int>>& fakesourcex, vector<vector<int>>& fakesourcey);
int getriverwidth(region& region, int x, int y, bool delta, int season);
void pasterivercircle(region& region, int centrex, int centrey, int pixels, bool river, int dir, int janload, int julload, bool delta, vector<vector<int>>& fakesourcex, vector<vector<int>>& fakesourcey);
void turnriverstolakes(planet& world, region& region, int dx, int dy, int sx, int sy);
void finishrivers(planet& world, region& region, int leftx, int lefty, int rightx, int righty);
void drawriverline(region& region, int leftx, int lefty, int rightx, int righty, int startx, int starty, int endx, int endy, int janload, int julload, int elev);
void makeelevationtile(planet& world, region& region, int dx, int dy, int sx, int sy, int coords[4][2], boolshapetemplate smalllake[]);
int rsquare(region& region, int dx, int dy, int s, int x, int y, int value, int min, int max, int coords[4][2], bool onlyup, int negchance, int sealevel, bool nosea);
int rcoastsquare(region& region, int dx, int dy, int s, int x, int y, int value, int min, int max, int coords[4][2], int sealevel);
int rlakesquare(region& region, int dx, int dy, int s, int x, int y, int value, int min, int max, int coords[4][2], int lakesurface);
int rlakediamond(region& region, int dx, int dy, int s, int x, int y, int value, int min, int max, int coords[4][2], int lakesurface);
int rdiamond(region& region, int dx, int dy, int s, int x, int y, int value, int min, int max, int coords[4][2], bool onlyup, int negchance, int sealevel, bool nosea);
int rcoastdiamond(region& region, int dx, int dy, int s, int x, int y, int value, int min, int max, int coords[4][2], int sealevel);
//void landfill(const planet &world, region &region, int dx, int dy, int sx, int sy, int surfaceleve, boolshapetemplate smalllake[]);
void warpcoasts(planet& world, region& region, int leftx, int lefty, int rightx, int righty, vector<vector<int>>& warptox, vector<vector<int>>& warptoy);
void removestraights(planet& world, region& region, int dx, int dy, int sx, int sy, boolshapetemplate smalllake[]);
void disruptseacoastline(planet& world, region& region, int dx, int dy, int centrex, int centrey, int avedepth, bool raise, int maxsize, bool stayintile, boolshapetemplate smalllake[]);
void disruptlakecoastline(planet& world, region& region, int dx, int dy, int centrex, int centrey, int surfacelevel, int avedepth, bool raise, int size, bool stayintile, int special, boolshapetemplate smalllake[]);
void removeseatilerivers(planet& world, region& region, int dx, int dy, int sx, int sy);
void removesearivers(planet& world, region& region, int dx, int dy);
void removeextrasearivers(planet& world, region& region, int dx, int dy, int sx, int sy, vector<vector<bool>>& regionsea);
void removeriverscomingfromsea(planet& world, region& region, int dx, int dy, int sx, int sy, vector<vector<int>>& fakesourcex, vector<vector<int>>& fakesourcey, vector<vector<bool>>& regionsea);
void removeregionalriver(region& region, int dx, int dy, int sx, int sy, int startx, int starty, vector<vector<int>>& fakesourcex, vector<vector<int>>& fakesourcey);
int findsurroundingsea(region& region, int x, int y);
int findinflowinglandriver(region& region, int x, int y);
int checknearbyriver(region& region, int x, int y);
void removeriverlakeloops(planet& world, region& region, int dx, int dy, int sx, int sy, boolshapetemplate smalllake[]);
void makemountaintile(planet& world, region& region, int dx, int dy, int sx, int sy, peaktemplate& peaks, vector<vector<int>>& rmountainmap, vector<vector<int>>& ridgeids, short markgap, vector<vector<int>>& warptox, vector<vector<int>>& warptoy);
void calculateridges(planet& world, region& region, int dx, int dy, int sx, int sy, twofloats pt, twofloats mm1, twofloats mm2, twofloats mm3, int newheight, int midheight, int destheight, peaktemplate& peaks, vector<vector<int>>& rmountainmap, int ridgedir, vector<vector<int>>& ridgeids, short markgap, vector<vector<int>>& warptox, vector<vector<int>>& warptoy);
void assignridgeregions(planet& world, region& region, int dx, int dy, int sx, int sy, vector<vector<int>>& rmountainmap, vector<vector<int>>& ridgeids, vector<vector<int>>& ridgeregions, vector<vector<int>>& nearestridgepointx, vector<vector<int>>& nearestridgepointy, int smallermaxdist, int maxdist);
void findmountainedges(planet& world, region& region, int dx, int dy, int sx, int sy, vector<vector<int>>& nearestridgepointdist, vector<vector<int>>& nearestridgepointx, vector<vector<int>>& nearestridgepointy, vector<vector<bool>>& mountainedges);
void findbuttresspoints(planet& world, region& region, int dx, int dy, int sx, int sy, vector<vector<int>>& ridgeids, vector<vector<int>>& nearestridgepointdist, vector<vector<int>>& nearestridgepointx, vector<vector<int>>& nearestridgepointy, vector<vector<bool>>& mountainedges, vector<vector<bool>>& buttresspoints, int maxdist, int spacing);
void makebuttresses(planet& world, region& region, int dx, int dy, int sx, int sy, vector<vector<int>>& rmountainmap, peaktemplate& peaks, vector<vector<int>>& nearestridgepointx, vector<vector<int>>& nearestridgepointy, vector<vector<int>>& nearestridgepointdist, int maxdist, vector<vector<bool>>& buttresspoints, vector<vector<int>>& ridgeids, short markgap, bool minibuttresses);
void drawpeak(planet& world, region& region, int sx, int sy, int dx, int dy, int x, int y, int peakheight, peaktemplate& peaks, vector<vector<int>>& rmountainmap, bool buttress);
void pastepeak(planet& world, region& region, int x, int y, float peakheight, int templateno, bool leftr, bool downr, peaktemplate& peaks, vector<vector<int>>& rmountainmap);
void removepools(planet& world, region& region, int dx, int dy, int sx, int sy, vector<vector<int>>& pathchecked, int& checkno);
bool findpath(region& region, int& leftx, int& lefty, int& rightx, int& righty, int& fromx, int& fromy, int& destx, int& desty, int& checkno, vector<vector<int>>& pathchecked, int& recursion);
void turnpoolstolakes(planet& world, region& region, int dx, int dy, int sx, int sy, vector<vector<bool>>& regionsea, vector<vector<int>>& pathchecked, int& checkno);
void poolcheckrecursive(region const& region, int const currentx, int const currenty, int& tally, int maxtally, int const checkno, vector<vector<bool>> const& regionsea, vector<vector<int>>& pathchecked);
void poolcheck(region const& region, int const currentx, int const currenty, int& tally, int maxtally, int const checkno, vector<vector<bool>> const& regionsea, vector<vector<int>>& pathchecked);
void turntosea(region& region, int leftx, int rightx, int lefty, int righty, int x, int y, int newheight, int sealevel);
void removediagonalwater(region& region, int leftx, int lefty, int rightx, int righty, int sealevel);
void removelakesbysea(region& region, int leftx, int lefty, int rightx, int righty, int sealevel);
void removesinks(planet& world, region& region, int dx, int dy, int sx, int sy);
void fillregionaldepression(region& region, int dx, int dy, int x, int y, int elev);
void checkdepression(region& region, int dx, int dy, int x, int y, int elev, bool& overrun, vector<vector<bool>>& checked);
void filldepression(region& region, int dx, int dy, int x, int y, int elev);
void addinlets(planet& world, region& region, int dx, int dy, int sx, int sy, vector<vector<bool>>& riverinlets, vector<vector<bool>>& globalestuaries);
void pasteinletcircle(region& region, int centrex, int centrey, int depth, int pixels, vector<vector<bool>>& riverinlets);
//void addmountainsprings(const planet &world, region &region, int dx, int dy, int sx, int sy);
void makewetlandtile(planet& world, region& region, int dx, int dy, int sx, int sy, boolshapetemplate smalllake[]);
void pasteregionalwetlands(region& region, int centrex, int centrey, int special, int elev, int shapenumber, boolshapetemplate smalllake[]);
void convertlakestospecials(planet& world, region& region, int dx, int dy, int sx, int sy, vector<vector<bool>>& safesaltlakes);
void removedeltasea(planet& world, region& region, int dx, int dy, int sx, int sy);
void adddeltamap(planet& world, region& region, int leftx, int lefty, int rightx, int righty);
void makegenerictile(planet& world, int dx, int dy, int sx, int sy, float valuemod, int coords[4][2], vector<vector<int>>& source, vector<vector<int>>& dest, int max, int min, bool interpolate);
int genericsquare(int dx, int dy, int s, int x, int y, int value, int min, int max, int coords[4][2], vector<vector<int>>& dest);
int genericdiamond(int dx, int dy, int s, int x, int y, int value, int min, int max, int coords[4][2], vector<vector<int>>& dest);
void smoothprecipitation(region& region, int leftx, int lefty, int rightx, int righty, int amount);
int getsurroundingice(planet& world, int x, int y);
void makesmallsaltpans(planet& world, region& region, int dx, int dy, int sx, int sy, vector<vector<bool>>& safesaltlakes, boolshapetemplate smalllake[]);
void makesaltpantile(planet& world, region& region, int dx, int dy, int sx, int sy, boolshapetemplate smalllake[]);
void pastesaltpan(region& region, int centrex, int centrey, int surfacelevel, boolshapetemplate smalllake[]);
void removewetsaltpans(region& region, int leftx, int lefty, int rightx, int righty);
void addbarrierislands(planet& world, region& region, int dx, int dy, int sx, int sy, vector<vector<bool>>& riverinlets);
void addregionalglaciers(planet& world, region& region, int dx, int dy, int sx, int sy);
void makeriftlaketemplates(planet& world, region& region, int dx, int dy, int sx, int sy, int extra, vector<vector<bool>>& riftlakemap);
void makeriftlaketile(planet& world, region& region, int dx, int dy, int sx, int sy, int extra, vector<vector<bool>>& riftlakemap, vector<vector<float>>& riftblob, int riftblobsize);
void makelaketile(planet& world, region& region, int dx, int dy, int sx, int sy, vector<vector<bool>>& lakemap, int surfacelevel, int coords[4][2], vector<vector<int>>& source, vector<vector<int>>& dest, vector<vector<bool>>& safesaltlakes, boolshapetemplate smalllake[]);
void makelakeislands(planet& world, region& region, int dx, int dy, int sx, int sy, int surfacelevel, boolshapetemplate island[], vector<vector<bool>>& lakeislands);
void createlakeisland(planet& world, region& region, int centrex, int centrey, int surfacelevel, boolshapetemplate island[], vector<vector<bool>>& lakeislands, bool nooverlap, int special);
void complicatecoastlines(planet& world, region& region, int dx, int dy, int sx, int sy, boolshapetemplate smalllake[], int chance);
void addregionalmountainprecipitation(planet& world, region& region, int dx, int dy, int sx, int sy, bool summer);
void removeregionalstraightrivers(planet& world, region& region, int dx, int dy, int sx, int sy, vector<vector<bool>>& rivercurves);
bool checkrivercurve(region& region, int dx, int dy, twofloats pt, twofloats mm1, twofloats mm2, twofloats mm3);
void makenewrivercurve(region& region, twofloats pt, twofloats mm1, twofloats mm2, twofloats mm3, int riverjan, int riverjul);
void addterraces(planet& world, region& region, int dx, int dy, int sx, int sy, boolshapetemplate smalllake[], byteshapetemplate smudge[]);
void addlaketerraces(planet& world, region& region, int dx, int dy, int sx, int sy, boolshapetemplate smalllake[], byteshapetemplate smudge[]);
void disruptcliff(planet& world, region& region, int dx, int dy, int sx, int sy, int startx, int starty, int endx, int endy, boolshapetemplate smalllake[], byteshapetemplate smudge[]);
void disruptlakecliff(planet& world, region& region, int dx, int dy, int sx, int sy, int startx, int starty, int endx, int endy, boolshapetemplate smalllake[], byteshapetemplate smudge[]);
void disruptland(region& region, int centrex, int centrey, int newheight, boolshapetemplate smalllake[]);
void disruptlakebed(region& region, int centrex, int centrey, int newheight, boolshapetemplate smalllake[]);
//void checkgrid(const planet &world, region &region, int dx, int dy, int sx, int sy, vector<vector<int>> &elevs, vector<vector<int>> &severities);
void makesubmarineelevationtile(planet& world, region& region, int dx, int dy, int sx, int sy, vector<vector<int>>& underseamap, int coords[4][2], int extra);
int submarinesquare(vector<vector<int>>& underseamap, int dx, int dy, int s, int x, int y, int value, int min, int max, int coords[4][2], bool onlyup);
int submarinediamond(vector<vector<int>>& underseamap, int dx, int dy, int s, int x, int y, int value, int min, int max, int coords[4][2], bool onlyup);
int disruptsubmarineelevationtile(planet& world, region& region, int dx, int dy, int sx, int sy, vector<vector<int>>& underseamap, byteshapetemplate smudge[], int extra);
void smudgesubmarineterrain(planet& world, region& region, int centrex, int centrey, int searchdist, vector<vector<int>>& underseamap, byteshapetemplate smudge[]);
void smudgeterrain(planet& world, region& region, int centrex, int centrey, int searchdist, byteshapetemplate smudge[]);
void makesubmarineridgelines(planet& world, region& region, int dx, int dy, int s, int sy, vector<vector<bool>>& undersearidgelines, vector<vector<int>>& warptox, vector<vector<int>>& warptoy);
void marksubmarineridgeline(region& region, vector<vector<bool>>& undersearidgelines, int fromx, int fromy, int tox, int toy, vector<vector<int>>& warptox, vector<vector<int>>& warptoy);
void drawsubmarineridges(planet& world, region& region, int dx, int dy, int sx, int sy, vector<vector<bool>>& undersearidgelines, peaktemplate& peaks, vector<vector<int>>& undersearidges);
void pastesubmarinepeak(planet& world, region& region, int x, int y, float peakheight, int templateno, peaktemplate& peaks, vector<vector<int>>& undersearidges, int rwidth, int rheight);
void makesubmarineriftradiations(planet& world, region& region, int dx, int dy, int sx, int sy, vector<vector<int>>& undersearidges, peaktemplate& peaks, int extra);
void drawriftspine(planet& world, region& region, int dx, int dy, int sx, int sy, int x1, int y1, int x2, int y2, float peakheight, float heightstep, vector<vector<int>>& underseaspikes, peaktemplate& peaks, bool lower);
void makesubmarineriftmountains(planet& world, region& region, int dx, int dy, int sx, int sy, vector<vector<int>>& undersearidges, peaktemplate& peaks, vector<vector<int>>& warptox, vector<vector<int>>& warptoy);
void makeoceanicriftmountains(planet& world, region& region, int sx, int sy, int fromx, int fromy, int tox, int toy, vector<vector<int>>& undersearidges, peaktemplate& peaks, vector<vector<int>>& warptox, vector<vector<int>>& warptoy);
void makesubmarinerift(planet& world, region& region, int dx, int dy, int sx, int sy, vector<vector<int>>& undersearidges, byteshapetemplate smudge[], vector<vector<int>>& warptox, vector<vector<int>>& warptoy);
void makesubmarineriftvalley(planet& world, region& region, int sx, int sy, int fromx, int fromy, int tox, int toy, int tileriftheight, vector<vector<int>>& undersearidges, byteshapetemplate smudge[], vector<vector<int>>& warptox, vector<vector<int>>& warptoy);
void trimmountainislands(planet& world, region& region, int dx, int dy, int sx, int sy, vector<vector<int>>& rmountainmap, vector<vector<bool>>& riverinlets, vector<vector<bool>>& globalestuaries, boolshapetemplate smalllake[]);
void removetoohighelevations(planet& world, region& region, int dx, int dy, int sx, int sy);
void makevolcano(planet& world, region& region, int dx, int dy, int sx, int sy, peaktemplate& peaks, vector<vector<int>>& rmountainmap, vector<vector<int>>& ridgeids, int templateno);
void makesubmarinevolcano(planet& world, region& region, int dx, int dy, int sx, int sy, peaktemplate& peaks, vector<vector<int>>& undersearidges, vector<vector<bool>>& volcanomap);
void fixrivers(planet& world, region& region, int leftx, int lefty, int rightx, int righty);
void findcoastdiagonals(planet& world, region& region, int dx, int dy, int sx, int sy, vector<vector<bool>>& disruptpoints);
void removecoastdiagonals(planet& world, region& region, int dx, int dy, int sx, int sy, vector<vector<bool>>& disruptpoints, boolshapetemplate smalllake[]);
void findlakecoastdiagonals(planet& world, region& region, int dx, int dy, int sx, int sy, vector<vector<bool>>& disruptpoints);
void removelakecoastdiagonals(planet& world, region& region, int dx, int dy, int sx, int sy, vector<vector<bool>>& disruptpoints, boolshapetemplate smalllake[]);
void rotatetileedges(planet& world, region& region, int dx, int dy, int sx, int sy, vector<vector<int>>& rotatearray, bool lakes);
void rotateland(planet& world, region& region, int centrex, int centrey, int maxradius, float angle, vector<vector<int>>& rotatearray);
void rotatelakes(planet& world, region& region, int centrex, int centrey, int maxradius, float angle, vector<vector<int>>& rotatearray);
void rotatetileedgesarray(planet& world, region& region, int dx, int dy, int sx, int sy, vector<vector<int>>& destarray, vector<vector<int>>& rotatearray, int min);
void rotatelandarray(planet& world, region& region, int centrex, int centrey, int maxradius, float angle, vector<vector<int>>& destarray, vector<vector<int>>& rotatearray, int min);
void removesmallislands(planet& world, region& region, int dx, int dy, int sx, int sy);
void removeparallelislands(region& region, int leftx, int lefty, int rightx, int righty, int sealevel);
void removeseaiceglitches(region& region);
void smoothlakebeds(region& region);
void removetoolow(planet& world, region& region, int dx, int dy, int sx, int sy);
void removelakeseas(planet& world, region& region, int dx, int dy, int sx, int sy);
void checklakebeds(region& region, int leftx, int lefty, int rightx, int righty);
void checkanomalies(planet& world, region& region, vector<vector<int>>& dest, int leftx, int lefty, int rightx, int righty);
void makeregionalcraterrim(planet& world, region& region, int dx, int dy, int sx, int sy, int leftx, int lefty, int rightx, int righty, int craterno, vector<vector<int>>& rcratermap, peaktemplate& peaks, vector<int>squareroot);
void makemudflatseeds(planet& world, region& region, int dx, int dy, int sx, int sy, vector<vector<bool>>& riverinlets, vector<vector<float>>& siltstrength);
void spreadsilt(planet& world, region& region, int leftx, int lefty, int rightx, int righty, vector<vector<float>>& siltstrength, boolshapetemplate smalllake[]);
void makesandseeds(planet& world, region& region, int dx, int dy, int sx, int sy, vector<vector<bool>>& riverinlets, vector<vector<float>>& sandstrength, vector<vector<float>>& shinglestrength);
bool spreadsand(planet& world, region& region, int dx, int dy, int sx, int sy, float reduce, float diagreduce, float landreduce, float landdiagreduce, float coastreduce, float coastdiagreduce, int crount, bool shingle, vector<vector<bool>>& checked, vector<vector<float>>& sandstrength);
void createbeaches(planet& world, region& region, int leftx, int lefty, int rightx, int righty, bool shingle, vector<vector<float>>& sandstrength);
void putmudonsand(planet& world, region& region, int leftx, int lefty, int rightx, int righty);
void checkbeachcoasts(planet& world, region& region, int leftx, int lefty, int rightx, int righty);
void makemudbarriers(planet& world, region& region, int dx, int dy, int sx, int sy);
void checklakesurfaces(planet& world, region& region);


// Declare functions that are in shape_templates.cpp

void loadpeaktemplates(peaktemplate& peaks);
int createmountainrangetemplate(vector<vector<unsigned char>>& dirtemplate, vector<vector<unsigned char>>& heighttemplate, int tempno, int dir);
int createsmalllaketemplates(boolshapetemplate smalllake[]);
int createislandtemplates(boolshapetemplate island[]);
int createsmudgetemplates(byteshapetemplate smudge[]);
int createsmallsmudgetemplates(byteshapetemplate smallsmudge[]);
int createlargelaketemplates(boolshapetemplate largelake[]);
int createlandshapetemplates(boolshapetemplate landshape[]);
int createchainlandtemplates(boolshapetemplate chainland[]);

#endif /* functions_hpp */
