#pragma once

#include "imgui.h"

#include "planet.hpp"

struct AppearanceSettings
{
    ImVec4 oceancolour;
    ImVec4 deepoceancolour;
    ImVec4 basecolour;
    ImVec4 grasscolour;
    ImVec4 basetempcolour;
    ImVec4 highbasecolour;
    ImVec4 desertcolour;
    ImVec4 highdesertcolour;
    ImVec4 colddesertcolour;
    ImVec4 eqtundracolour;
    ImVec4 tundracolour;
    ImVec4 coldcolour;
    ImVec4 seaicecolour;
    ImVec4 glaciercolour;
    ImVec4 saltpancolour;
    ImVec4 ergcolour;
    ImVec4 wetlandscolour;
    ImVec4 lakecolour;
    ImVec4 rivercolour;
    ImVec4 sandcolour;
    ImVec4 mudcolour;
    ImVec4 shinglecolour;
    ImVec4 mangrovecolour;
    ImVec4 highlightcolour;

    float shadingland = 0.0f;
    float shadinglake = 0.0f;
    float shadingsea = 0.0f;
    float marblingland = 0.0f;
    float marblinglake = 0.0f;
    float marblingsea = 0.0f;

    int globalriversentry = 0;
    int regionalriversentry = 0;
    int shadingdir = 0;
    int snowchange = 0;
    int seaiceappearance = 0;
    bool colourcliffs = false;
    bool mangroves = false;
};

struct AppearanceChangeFlags
{
    bool mapAppearanceChanged = false;
    bool highlightChanged = false;

    bool any() const
    {
        return mapAppearanceChanged || highlightChanged;
    }
};

AppearanceSettings makeappearancesettings(const planet& world);
void syncappearancesettings(const planet& world, AppearanceSettings& settings);
AppearanceChangeFlags getappearancechanges(const planet& world, const AppearanceSettings& settings);
void applyappearancesettings(planet& world, AppearanceSettings& settings);
