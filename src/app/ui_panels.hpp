#pragma once

#include <string>

#include "imgui.h"

#include "world/planet.hpp"
#include "simulations/regional/regional_map.hpp"

void drawglobalinfowindow(const ImGuiViewport* main_viewport, ImGuiWindowFlags window_flags, planet& world, bool focused, int poix, int poiy, const std::string& degree, const std::string& cube, bool& showtemperaturechart, bool& showrainfallchart, bool newworld, bool& showabout);
void drawregionalinfowindow(const ImGuiViewport* main_viewport, ImGuiWindowFlags window_flags, planet& world, region& region, bool focused, int poix, int poiy, const std::string& degree, const std::string& cube, bool& showtemperaturechart, bool& showrainfallchart, bool& showabout);
