#pragma once

#include <array>
#include <functional>

#include "appearance_settings.hpp"
#include "functions.hpp"

void drawmapviewbuttons(const std::function<void(mapviewenum)>& onSelect);
void drawmapviewappearancetab(const mapviewdefinition& definition, planet& world, AppearanceSettings& appearance, std::array<int, MAPGRADIENTTYPECOUNT>& selectedgradientstops, int colouralign, int otheralign);
void normalizegradientstops(MapGradientSettings& gradient);
ImVec4 samplegradientcolour(const MapGradientSettings& gradient, int value);
bool drawgradienteditor(const char* id, MapGradientSettings& gradient, int minimum, int maximum, const char* unittext, int& selectedstop);
