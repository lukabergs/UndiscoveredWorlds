#pragma once

#include <string>

#include "imgui.h"

#include "world_generation_debug.hpp"
#include "world_property_controls.hpp"

bool drawcustomworldsizewindow(const ImGuiViewport* main_viewport, ImGuiWindowFlags window_flags, bool& show, bool brandnew, int& currentsize);
bool drawworldgenerationoptionswindow(const ImGuiViewport* main_viewport, ImGuiWindowFlags window_flags, bool& show, WorldGenerationDebugOptions& options);
bool drawtectonicchooserwindow(const ImGuiViewport* main_viewport, ImGuiWindowFlags window_flags, bool& show, int& landmass, int& mergefactor);
bool drawnontectonicchooserwindow(const ImGuiViewport* main_viewport, ImGuiWindowFlags window_flags, bool& show, int& sealeveleditable, int& iterations);
bool drawareawarningwindow(const ImGuiViewport* main_viewport, ImGuiWindowFlags window_flags, bool& show);
void drawaboutwindow(const ImGuiViewport* main_viewport, ImGuiWindowFlags window_flags, bool& show, float currentversion);
void drawworldeditpropertieswindow(const ImGuiViewport* main_viewport, ImGuiWindowFlags window_flags, WorldPropertyControls& controls, const std::string& degree, bool& show);
