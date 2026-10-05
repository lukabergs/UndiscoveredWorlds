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
//  The former misc.cpp is divided by owner: core/text_values, core/grid_shapes, world/world_setup, app/regional/region_queries, simulations/terrain/coast_queries, simulations/hydrology/flow_queries and app/appearance_files.
//  Simulation modules under src/simulations contain terrain, climate, and hydrology.
//  app/regional/regional_map.cpp contains the functions used to generate the regional terrain.
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

#define ARRAYWIDTH 2048
#define ARRAYHEIGHT 1024


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

// Simulation entry points live with their implementations. This header
// remains a compatibility umbrella for the desktop and the simulation modules.
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
