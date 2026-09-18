#pragma once

#include <string>

#include "imgui.h"

#include "world_property_controls.hpp"

bool drawcustomworldsizewindow(const ImGuiViewport* main_viewport, ImGuiWindowFlags window_flags, bool& show, bool brandnew, WorldPropertyControls& controls);
bool drawworldgenerationoptionswindow(const ImGuiViewport* main_viewport, ImGuiWindowFlags window_flags, bool& show, WorldPropertyControls& controls);
bool drawterrainchooserwindow(const ImGuiViewport* main_viewport, ImGuiWindowFlags window_flags, bool& show);
bool drawareawarningwindow(const ImGuiViewport* main_viewport, ImGuiWindowFlags window_flags, bool& show);
void drawaboutwindow(const ImGuiViewport* main_viewport, ImGuiWindowFlags window_flags, bool& show, float currentversion);
void drawworldeditpropertieswindow(const ImGuiViewport* main_viewport, ImGuiWindowFlags window_flags, WorldPropertyControls& controls, const std::string& degree, bool& show);
