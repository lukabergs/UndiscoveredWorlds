#pragma once

#include <vector>

#include "imgui.h"

#include "appearance_settings.hpp"
#include "planet.hpp"
#include "region.hpp"

bool drawheighthistogramwidget(const char* label, const std::vector<float>& values, float minimumvalue, float maximumvalue, float& sealevelvalue, const MapGradientSettings* gradient, const char* unittext, const ImVec2& size = ImVec2(0.0f, 0.0f));
bool drawhypsometriccurvewidget(const char* label, const std::vector<ImVec2>& worldcurve, std::vector<ImVec2>& controlpoints, std::vector<float>& weights, int& selectedcontrolpoint, float& sealevelvalue, const ImVec2& size = ImVec2(0.0f, 0.0f), const std::vector<ImVec2>* previewcurve = nullptr, bool editable = true);

void drawglobaltemperaturechartwindow(const ImGuiViewport* main_viewport, ImGuiWindowFlags window_flags, planet& world, int poix, int poiy);
void drawglobalrainfallchartwindow(const ImGuiViewport* main_viewport, ImGuiWindowFlags window_flags, planet& world, int poix, int poiy);
void drawregionaltemperaturechartwindow(const ImGuiViewport* main_viewport, ImGuiWindowFlags window_flags, planet& world, region& region, int poix, int poiy);
void drawregionalrainfallchartwindow(const ImGuiViewport* main_viewport, ImGuiWindowFlags window_flags, region& region, int poix, int poiy);
