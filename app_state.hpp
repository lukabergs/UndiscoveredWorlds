#pragma once

#include "functions.hpp"

struct NavigationState
{
    bool brandnew = true;
    bool newworld = false;
    screenmodeenum screenmode = createworldscreen;
    screenmodeenum oldscreenmode = createworldscreen;
    mapviewenum mapview = relief;
};

struct FocusState
{
    int poix = 0;
    int poiy = 0;
    bool focused = false;
    bool regionalFocused = false;
};

struct AreaSelectionState
{
    int northwestX = -1;
    int northwestY = -1;
    int northeastX = -1;
    int northeastY = -1;
    int southeastX = -1;
    int southeastY = -1;
    int southwestX = -1;
    int southwestY = -1;
    bool fromRegional = false;
};

struct PanelState
{
    bool showColourOptions = false;
    bool showWorldProperties = false;
    bool showGlobalTemperatureChart = false;
    bool showGlobalRainfallChart = false;
    bool showRegionalTemperatureChart = false;
    bool showRegionalRainfallChart = false;
    bool showSetSize = false;
    bool showTectonicChooser = false;
    bool showNonTectonicChooser = false;
    bool showWorldGenerationOptions = false;
    bool showWorldEditProperties = false;
    bool showAreaWarning = false;
    bool showAbout = false;
    bool coloursChanged = false;
};

struct CustomWorldUiState
{
    int seedentry = 0;
    int newx = -1;
    int newy = -1;
    int landmass = 5;
    int mergefactor = 15;
    int iterations = 4;
    int sealeveleditable = 5;
    bool compareClimateWorkbook = true;
};

struct ProgressPassState
{
    short creatingWorld = 0;
    short completingImport = 0;
    short loadingWorld = 0;
    short savingWorld = 0;
    short exportingArea = 0;
    short generatingRegion = 0;
    short generatingTectonic = 0;
    short generatingNonTectonic = 0;
};

struct WorldGenerationDebugState
{
    WorldGenerationDebugOptions options;
};
