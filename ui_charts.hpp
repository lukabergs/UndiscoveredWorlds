#pragma once

#include "imgui.h"

#include "planet.hpp"
#include "region.hpp"

void drawglobaltemperaturechartwindow(const ImGuiViewport* main_viewport, ImGuiWindowFlags window_flags, planet& world, int poix, int poiy);
void drawglobalrainfallchartwindow(const ImGuiViewport* main_viewport, ImGuiWindowFlags window_flags, planet& world, int poix, int poiy);
void drawregionaltemperaturechartwindow(const ImGuiViewport* main_viewport, ImGuiWindowFlags window_flags, planet& world, region& region, int poix, int poiy);
void drawregionalrainfallchartwindow(const ImGuiViewport* main_viewport, ImGuiWindowFlags window_flags, region& region, int poix, int poiy);
