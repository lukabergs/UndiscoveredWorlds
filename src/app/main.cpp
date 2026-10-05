
#include "imgui.h"
#include "imgui-SFML.h"
#include "ImGuiFileDialog.h"
#include "app/create_world_screen.hpp"

#include <cstdio>

#include <windows.h>
#include <shellapi.h>

#include <algorithm>
#include <array>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <filesystem>
#include <cmath>
#include <fstream>
#include <memory>
#include <stdio.h>
#include <string>
#include <type_traits>
#include <utility>
#include <chrono>
#include <thread>
#include <queue>
#include <stdint.h>
#include <SFML/Graphics.hpp>
#include <SFML/Graphics/CircleShape.hpp>
#include <SFML/Graphics/RenderWindow.hpp>
#include <SFML/System/Clock.hpp>
#include <SFML/Window/Event.hpp>

#include "core/classes.hpp"
#include "world/planet.hpp"
#include "pipeline/physical_layers.hpp"
#include "simulations/society/social_generation.hpp"
#include "world/region.hpp"
#include "app/app_dialogs.hpp"
#include "core/app_environment.hpp"
#include "app/app_state.hpp"
#include "app/app_windows.hpp"
#include "app/rendering/appearance_settings.hpp"
#include "world/functions.hpp"
#include "app/appearance_files.hpp"
#include "app/rendering/map_views.hpp"
#include "app/generation_workbench.hpp"
#include "app/rendering/map_appearance.hpp"
#include "io/map_imports.hpp"
#include "app/ui_charts.hpp"
#include "app/ui_panels.hpp"
#include "app/world_property_controls.hpp"
#include "cli/command_line_generation.hpp"

// nodiscard attribute
#pragma warning (disable: 4834)

#define REGIONALCREATIONSTEPS 84

#define REGIONALTILEWIDTH 32
#define REGIONALTILEHEIGHT 32

using namespace std;

namespace
{
constexpr const char* kRasterImportFilter = ".png,.tif,.tiff";
constexpr const char* kGradientStripFilter = ".png";
constexpr const char* kVolcanoImportFilter = ".png";

const char* crustclasslabel(CrustClass crustclass)
{
    switch (crustclass)
    {
    case CrustClass::oceanic:
        return "oceanic";
    case CrustClass::transitional:
        return "transitional";
    case CrustClass::continental:
        return "continental";
    case CrustClass::none:
    default:
        return "none";
    }
}

const char* boundarytypelabel(BoundaryType boundarytype)
{
    switch (boundarytype)
    {
    case BoundaryType::convergent:
        return "convergent";
    case BoundaryType::divergent:
        return "divergent";
    case BoundaryType::transform:
        return "transform";
    case BoundaryType::passive_margin:
        return "passive_margin";
    case BoundaryType::none:
    default:
        return "none";
    }
}

const char* deformingregiontypelabel(DeformingRegionType regiontype)
{
    switch (regiontype)
    {
    case DeformingRegionType::continental_rift:
        return "continental_rift";
    case DeformingRegionType::diffuse_collision:
        return "diffuse_collision";
    case DeformingRegionType::none:
    default:
        return "none";
    }
}

const char* geologicregimelabel(GeologicRegime regime)
{
    switch (regime)
    {
    case GeologicRegime::convergent_arc:
        return "convergent_arc";
    case GeologicRegime::continent_collision:
        return "continent_collision";
    case GeologicRegime::divergent_rift:
        return "divergent_rift";
    case GeologicRegime::transform:
        return "transform";
    case GeologicRegime::passive_margin:
        return "passive_margin";
    case GeologicRegime::mid_ocean_ridge:
        return "mid_ocean_ridge";
    case GeologicRegime::trench_adjacent:
        return "trench_adjacent";
    case GeologicRegime::stable:
    default:
        return "stable";
    }
}

void drawboundarysegmentsummary(const char* label, int requestedid, const TectonicBoundarySegment* segment)
{
    if (requestedid <= 0)
        return;

    if (segment == nullptr)
    {
        ImGui::Text("%s: id %d (not found)", label, requestedid);
        return;
    }

    ImGui::Text("%s: id %d | type %s | plates %d/%d | cells %d | age %.2f Myr",
        label, segment->id, boundarytypelabel(segment->boundaryType), segment->leftPlateId, segment->rightPlateId,
        segment->cellCount, static_cast<float>(segment->ageMyr));
    ImGui::Text("  centroid(%.1f, %.1f) | normal %.3f | shear %.3f | len %.1f | regime %s",
        segment->centroidX, segment->centroidY, segment->averageNormalMotion, segment->averageShearMotion,
        segment->lengthCells, geologicregimelabel(segment->geologicRegime));
}

void drawdeformingregionsummary(const char* label, int requestedid, const TectonicDeformingRegion* region)
{
    if (requestedid <= 0)
        return;

    if (region == nullptr)
    {
        ImGui::Text("%s: id %d (not found)", label, requestedid);
        return;
    }

    ImGui::Text("%s: id %d | type %s | plates %d/%d | cells %d | age %.2f Myr",
        label, region->id, deformingregiontypelabel(region->type), region->primaryPlateId, region->secondaryPlateId,
        region->cellCount, static_cast<float>(region->ageMyr));
    ImGui::Text("  centroid(%.1f, %.1f) | rate %.3f | vel(%.3f, %.3f) | normal %.3f | shear %.3f",
        region->centroidX, region->centroidY, region->averageDeformationRate,
        region->averageInterpolatedVelocityX, region->averageInterpolatedVelocityY,
        region->averageNormalMotion, region->averageShearMotion);
}

// Deposits in the inspected cell, with the history that explains each one.
void drawmineraldepositsummary(const planet& world, int x, int y)
{
    const int columns = world.width() + 1;
    const int rows = world.height() + 1;
    int shown = 0;
    for (const mineralsystems::Deposit& deposit : world.mineraldeposits())
    {
        const int depositx = std::clamp(static_cast<int>(deposit.u * columns), 0, columns - 1);
        const int deposity = std::clamp(static_cast<int>(deposit.v * rows), 0, rows - 1);
        if (depositx != x || deposity != y)
            continue;
        if (shown++ == 0)
        {
            ImGui::Spacing();
            ImGui::Text("Mineral deposits");
            ImGui::Separator();
        }
        ImGui::TextWrapped("%s", mineralsystems::describedeposit(deposit, &world.mineraldeposits()).c_str());
    }
}

bool gettectonicinspectionpoint(const planet& world, const region& region, screenmodeenum screenmode, bool focused, int poix, int poiy, int& x, int& y)
{
    if (focused == false)
        return false;

    if (screenmode == globalmapscreen || screenmode == exportareascreen || screenmode == importscreen)
    {
        x = wrap(poix, world.width());
        y = std::clamp(poiy, 0, world.height());
        return true;
    }

    if (screenmode == regionalmapscreen || screenmode == generatingregionscreen)
    {
        x = wrap(region.leftx() + poix / 16, world.width());
        y = std::clamp(region.lefty() + poiy / 16, 0, world.height());
        return true;
    }

    return false;
}
}

int createrandomseednumber()
{
    int seed = random(0, 9);

    for (int n = 1; n <= 7; n++)
    {
        if (n == 7)
            seed = seed + (random(1, 9) * (int)pow(10, n));
        else
            seed = seed + (random(0, 9) * (int)pow(10, n));
    }

    return seed;
}

struct GlobalMapDisplay
{
    sf::Vector2i textureSize;
    sf::Vector2i displayTextureSize;
    mapcache maps;
    sf::Texture texture;
    sf::Sprite sprite;
    sf::Sprite minimap;
    float x = 180.f;
    float y = 20.f;
    float minimapX = 190.f + ARRAYWIDTH / 4;
    float minimapY = 20.f;
};

struct HighlightDisplay
{
    int size = 3;
    sf::Image image;
    sf::Texture texture;
    sf::Sprite sprite;
    int miniSize = 8;
    sf::Image miniImage;
    sf::Texture miniTexture;
    sf::Sprite miniSprite;
};

struct RegionalMapDisplay
{
    sf::Vector2i textureSize = { 514, 514 };
    int imageWidth = textureSize.x;
    int imageHeight = textureSize.y;
    mapcache maps;
    sf::Texture texture;
    sf::Sprite sprite;
    float x = 180.f;
    float y = 20.f;
};

void initglobalmapdisplay(planet& world, GlobalMapDisplay& display)
{
    display.textureSize.x = world.width() + 1;
    display.textureSize.y = world.height() + 1;

    for (mapviewenum thismapview : allmapviews)
        getmapimage(display.maps, thismapview).create(display.textureSize.x, display.textureSize.y, sf::Color::Black);

    display.displayTextureSize.x = DISPLAYMAPSIZEX;
    display.displayTextureSize.y = DISPLAYMAPSIZEY;

    for (mapviewenum thismapview : allmapviews)
        getdisplaymapimage(display.maps, thismapview).create(display.displayTextureSize.x, display.displayTextureSize.y, sf::Color::Black);

    updateTextureFromImage(display.texture, getdisplaymapimage(display.maps, relief));

    display.sprite.setTexture(display.texture);
    display.sprite.setPosition(sf::Vector2f(display.x, display.y));

    display.minimap.setTexture(display.texture);
    display.minimap.setScale(0.5f, 0.5f);
    display.minimap.setPosition(sf::Vector2f(display.minimapX, display.minimapY));
}

void inithighlightdisplay(planet& world, HighlightDisplay& display)
{
    display.image.create(display.size, display.size, sf::Color::Transparent);
    display.miniImage.create(display.miniSize, display.miniSize, sf::Color::Transparent);

    drawhighlightobjects(world, display.image, display.size, display.miniImage, display.miniSize);

    updateTextureFromImage(display.texture, display.image);
    display.sprite.setTexture(display.texture);
    display.sprite.setOrigin((float)(display.size / 2), (float)(display.size / 2));

    updateTextureFromImage(display.miniTexture, display.miniImage);
    display.miniSprite.setTexture(display.miniTexture);
    display.miniSprite.setOrigin((float)(display.miniSize / 2), (float)(display.miniSize / 2));
}

void initregionalmapdisplay(RegionalMapDisplay& display)
{
    for (mapviewenum thismapview : allmapviews)
        getmapimage(display.maps, thismapview).create(display.imageWidth, display.imageHeight, sf::Color::Black);

    updateTextureFromImage(display.texture, getmapimage(display.maps, relief));

    display.sprite.setTexture(display.texture);
    display.sprite.setPosition(sf::Vector2f(display.x, display.y));
}

void refreshminihighlightdisplay(HighlightDisplay& display)
{
    updateTextureFromImage(display.miniTexture, display.miniImage);
    display.miniSprite.setTexture(display.miniTexture);
    display.miniSprite.setOrigin((float)display.miniSize / 2.0f, (float)display.miniSize / 2.0f);
}

void refreshhighlightdisplay(planet& world, HighlightDisplay& display)
{
    drawhighlightobjects(world, display.image, display.size, display.miniImage, display.miniSize);

    updateTextureFromImage(display.texture, display.image);
    display.sprite.setTexture(display.texture);

    refreshminihighlightdisplay(display);
}

void resizeglobaldisplayforworld(planet& world, GlobalMapDisplay& globaldisplay, HighlightDisplay& highlightdisplay)
{
    adjustforsize(world, globaldisplay.textureSize, globaldisplay.maps, highlightdisplay.image, highlightdisplay.size, highlightdisplay.miniImage, highlightdisplay.miniSize);
    refreshminihighlightdisplay(highlightdisplay);
}

void configureimguistyle()
{
    ImGuiStyle& style = ImGui::GetStyle();
    const float highlight1 = 0.40f;
    const float highlight2 = 0.40f;
    const float highlight3 = 0.40f;

    style.WindowMinSize = ImVec2(160, 20);
    style.FramePadding = ImVec2(4, 2);
    style.ItemSpacing = ImVec2(6, 4);
    style.ItemInnerSpacing = ImVec2(2, 4);
    style.Alpha = 0.95f;
    style.WindowRounding = 4.0f;
    style.FrameRounding = 2.0f;
    style.IndentSpacing = 6.0f;
    style.ColumnsMinSpacing = 50.0f;
    style.GrabMinSize = 14.0f;
    style.GrabRounding = 16.0f;
    style.ScrollbarSize = 12.0f;
    style.ScrollbarRounding = 16.0f;

    style.Colors[ImGuiCol_Text] = ImVec4(0.86f, 0.93f, 0.89f, 0.78f);
    style.Colors[ImGuiCol_TextDisabled] = ImVec4(0.86f, 0.93f, 0.89f, 0.28f);
    style.Colors[ImGuiCol_WindowBg] = ImVec4(0.13f, 0.14f, 0.17f, 1.00f);
    style.Colors[ImGuiCol_ChildBg] = ImVec4(highlight1, highlight2, highlight3, 1.00f);
    style.Colors[ImGuiCol_PopupBg] = ImVec4(0.20f, 0.22f, 0.27f, 0.9f);
    style.Colors[ImGuiCol_Border] = ImVec4(0.31f, 0.31f, 1.00f, 0.00f);
    style.Colors[ImGuiCol_BorderShadow] = ImVec4(0.00f, 0.00f, 0.00f, 0.00f);
    style.Colors[ImGuiCol_FrameBg] = ImVec4(0.20f, 0.22f, 0.27f, 1.00f);
    style.Colors[ImGuiCol_FrameBgHovered] = ImVec4(0.20f, 0.22f, 0.27f, 1.00f);
    style.Colors[ImGuiCol_FrameBgActive] = ImVec4(0.20f, 0.22f, 0.27f, 1.00f);
    style.Colors[ImGuiCol_TitleBg] = ImVec4(0.20f, 0.22f, 0.27f, 1.00f);
    style.Colors[ImGuiCol_TitleBgCollapsed] = ImVec4(0.20f, 0.22f, 0.27f, 0.75f);
    style.Colors[ImGuiCol_TitleBgActive] = ImVec4(0.20f, 0.22f, 0.27f, 1.00f);
    style.Colors[ImGuiCol_MenuBarBg] = ImVec4(0.20f, 0.22f, 0.27f, 0.47f);
    style.Colors[ImGuiCol_ScrollbarBg] = ImVec4(0.20f, 0.22f, 0.27f, 1.00f);
    style.Colors[ImGuiCol_ScrollbarGrab] = ImVec4(0.09f, 0.15f, 0.16f, 1.00f);
    style.Colors[ImGuiCol_ScrollbarGrabHovered] = ImVec4(0.29f, 0.18f, 0.92f, 0.78f);
    style.Colors[ImGuiCol_ScrollbarGrabActive] = ImVec4(highlight1, highlight2, highlight3, 1.00f);
    style.Colors[ImGuiCol_CheckMark] = ImVec4(0.86f, 0.93f, 0.89f, 0.78f);
    style.Colors[ImGuiCol_SliderGrab] = ImVec4(0.47f, 0.77f, 0.83f, 0.14f);
    style.Colors[ImGuiCol_SliderGrabActive] = ImVec4(0.28f, 0.18f, 0.92f, 1.00f);
    style.Colors[ImGuiCol_Button] = ImVec4(0.47f, 0.77f, 0.83f, 0.14f);
    style.Colors[ImGuiCol_ButtonHovered] = ImVec4(highlight1, highlight2, highlight3, 0.86f);
    style.Colors[ImGuiCol_ButtonActive] = ImVec4(highlight1, highlight2, highlight3, 1.00f);
    style.Colors[ImGuiCol_Header] = ImVec4(highlight1, highlight2, highlight3, 0.76f);
    style.Colors[ImGuiCol_HeaderHovered] = ImVec4(highlight1, highlight2, highlight3, 0.86f);
    style.Colors[ImGuiCol_HeaderActive] = ImVec4(highlight1, highlight2, highlight3, 1.00f);
    style.Colors[ImGuiCol_Separator] = ImVec4(0.14f, 0.16f, 0.19f, 1.00f);
    style.Colors[ImGuiCol_SeparatorHovered] = ImVec4(highlight1, highlight2, highlight3, 0.78f);
    style.Colors[ImGuiCol_SeparatorActive] = ImVec4(highlight1, highlight2, highlight3, 1.00f);
    style.Colors[ImGuiCol_ResizeGrip] = ImVec4(0.47f, 0.77f, 0.83f, 0.04f);
    style.Colors[ImGuiCol_ResizeGripHovered] = ImVec4(highlight1, highlight2, highlight3, 0.78f);
    style.Colors[ImGuiCol_ResizeGripActive] = ImVec4(highlight1, highlight2, highlight3, 1.00f);
    style.Colors[ImGuiCol_PlotLines] = ImVec4(0.86f, 0.93f, 0.89f, 0.63f);
    style.Colors[ImGuiCol_PlotLinesHovered] = ImVec4(highlight1, highlight2, highlight3, 1.00f);
    style.Colors[ImGuiCol_PlotHistogram] = ImVec4(0.86f, 0.93f, 0.89f, 0.63f);
    style.Colors[ImGuiCol_PlotHistogramHovered] = ImVec4(highlight1, highlight2, highlight3, 1.00f);
    style.Colors[ImGuiCol_TextSelectedBg] = ImVec4(highlight1, highlight2, highlight3, 0.43f);
    style.Colors[ImGuiCol_Tab] = ImVec4(0.20f, 0.22f, 0.27f, 1.00f);
    style.Colors[ImGuiCol_TabHovered] = ImVec4(highlight1, highlight2, highlight3, 0.86f);
    style.Colors[ImGuiCol_TabActive] = ImVec4(highlight1, highlight2, highlight3, 0.86f);
    style.Colors[ImGuiCol_TabUnfocused] = ImVec4(0.20f, 0.22f, 0.27f, 1.00f);
    style.Colors[ImGuiCol_TabUnfocusedActive] = ImVec4(highlight1, highlight2, highlight3, 0.86f);
}

int main()
{
    ensureworkingdirectory();
    reloadappenvironment();

    int argc = 0;
    if (LPWSTR* argv = CommandLineToArgvW(GetCommandLineW(), &argc))
        LocalFree(argv);
    if (argc > 1)
    {
        cerr << "UndiscoveredWorlds is the desktop application and takes no arguments.\n"
                "Use world_pipeline.exe for headless generation and climate benchmarks.\n";
        return 1;
    }

    float currentversion = 1.0f;

    string currentversionstring = "Current version: " + to_string(currentversion);

    updatereport(currentversionstring.c_str());

    // Set up the window.

    int scwidth = 1224;
    int scheight = 768;

    sf::RenderWindow window(sf::VideoMode(scwidth, scheight), "Undiscovered Worlds");
    window.setFramerateLimit(60);
    ImGui::SFML::Init(window);
    ImGui::GetIO().IniFilename = "configs/imgui.ini";

    const ImGuiViewport* main_viewport = ImGui::GetMainViewport();

    // Setting up Dear ImGUI style
    // Adapted from Cinder-ImGui by Simon Geilfus - https://github.com/simongeilfus/Cinder-ImGui

    configureimguistyle();

    ImGuiIO& io = ImGui::GetIO();
    ImFont* font1 = io.Fonts->AddFontDefault();
    io.Fonts->Build();
    IM_ASSERT(font1 != NULL);
    ImGui::SFML::UpdateFontTexture();

    auto openFileDialog = [](const char* filter)
    {
        IGFD::FileDialogConfig dialogConfig;
        const AppEnvironmentConfig& appenv = getappenvironment();
        string dialogpath = appenv.defaultWorldDirectory.string();

        if (strcmp(filter, ".uws") == 0)
            dialogpath = appenv.defaultAppearanceDirectory.string();
        else if (strstr(filter, ".png") != nullptr || strstr(filter, ".tif") != nullptr)
            dialogpath = appenv.defaultImageDirectory.string();

        if (dialogpath.empty())
            dialogpath = ".";

        dialogConfig.path = dialogpath;
        ImGuiFileDialog::Instance()->OpenDialog("ChooseFileDlgKey", "Choose File", filter, dialogConfig);
    };

    ImGuiWindowFlags window_flags = 0;
    window_flags |= ImGuiWindowFlags_NoCollapse;
    window_flags |= ImGuiWindowFlags_NoResize;

    string degree = "\xC2\xB0"; // 0xC2 0xB0 (c2b0) ï¿½
    string cube = "\xC2\xB3"; // 0xC2 0xB3 (c2b3) ï¿½

    // Now create the main world object and initialise its basic variables.

    auto world = std::make_unique<planet>();

    initialiseworld(*world);
    initialisemapcolours(*world);

    int width = world->width();
    int height = world->height();

    // Now do the same thing for the region object.

    auto region = std::make_unique<class region>();

    initialiseregion(*world, *region);

    int regionalmapimagewidth = region->regwidthend() - region->regwidthbegin() + 1;
    int regionalmapimageheight = region->regheightend() - region->regheightbegin() + 1;

    // Create display/state bundles for the global map, highlights, and regional map.

    GlobalMapDisplay globaldisplay;
    initglobalmapdisplay(*world, globaldisplay);
    sf::Vector2i& globaltexturesize = globaldisplay.textureSize;
    sf::Vector2i& displayglobaltexturesize = globaldisplay.displayTextureSize;
    mapcache& globalmaps = globaldisplay.maps;
    sf::Texture& globalmaptexture = globaldisplay.texture;
    sf::Sprite& globalmap = globaldisplay.sprite;
    sf::Sprite& minimap = globaldisplay.minimap;
    float globalmapxpos = globaldisplay.x;
    float globalmapypos = globaldisplay.y;
    float minimapxpos = globaldisplay.minimapX;
    float minimapypos = globaldisplay.minimapY;

    sf::Image& globalelevationimage = getmapimage(globalmaps, elevation);
    sf::Image& globaltemperatureimage = getmapimage(globalmaps, temperature);
    sf::Image& globalprecipitationimage = getmapimage(globalmaps, precipitation);
    sf::Image& globalclimateimage = getmapimage(globalmaps, climate);
    sf::Image& globalbiomeimage = getmapimage(globalmaps, biomes);
    sf::Image& globalriversimage = getmapimage(globalmaps, rivers);
    sf::Image& globalreliefimage = getmapimage(globalmaps, relief);

    sf::Image& displayglobalelevationimage = getdisplaymapimage(globalmaps, elevation);
    sf::Image& displayglobaltemperatureimage = getdisplaymapimage(globalmaps, temperature);
    sf::Image& displayglobalprecipitationimage = getdisplaymapimage(globalmaps, precipitation);
    sf::Image& displayglobalclimateimage = getdisplaymapimage(globalmaps, climate);
    sf::Image& displayglobalbiomeimage = getdisplaymapimage(globalmaps, biomes);
    sf::Image& displayglobalriversimage = getdisplaymapimage(globalmaps, rivers);
    sf::Image& displayglobalreliefimage = getdisplaymapimage(globalmaps, relief);

    HighlightDisplay highlightdisplay;
    inithighlightdisplay(*world, highlightdisplay);
    int& highlightsize = highlightdisplay.size;
    sf::Image& highlightimage = highlightdisplay.image;
    sf::Texture& highlighttexture = highlightdisplay.texture;
    sf::Sprite& highlight = highlightdisplay.sprite;
    int& minihighlightsize = highlightdisplay.miniSize;
    sf::Image& minihighlightimage = highlightdisplay.miniImage;
    sf::Texture& minihighlighttexture = highlightdisplay.miniTexture;
    sf::Sprite& minihighlight = highlightdisplay.miniSprite;

    RegionalMapDisplay regionaldisplay;
    initregionalmapdisplay(regionaldisplay);
    sf::Vector2i& regionaltexturesize = regionaldisplay.textureSize;
    int& regionalimagewidth = regionaldisplay.imageWidth;
    int& regionalimageheight = regionaldisplay.imageHeight;
    mapcache& regionalmaps = regionaldisplay.maps;
    sf::Texture& regionalmaptexture = regionaldisplay.texture;
    sf::Sprite& regionalmap = regionaldisplay.sprite;
    float regionalmapxpos = regionaldisplay.x;
    float regionalmapypos = regionaldisplay.y;

    sf::Image& regionalelevationimage = getmapimage(regionalmaps, elevation);
    sf::Image& regionaltemperatureimage = getmapimage(regionalmaps, temperature);
    sf::Image& regionalprecipitationimage = getmapimage(regionalmaps, precipitation);
    sf::Image& regionalclimateimage = getmapimage(regionalmaps, climate);
    sf::Image& regionalbiomeimage = getmapimage(regionalmaps, biomes);
    sf::Image& regionalriversimage = getmapimage(regionalmaps, rivers);
    sf::Image& regionalreliefimage = getmapimage(regionalmaps, relief);

    // A rectangle for the area select screen

    sf::RectangleShape arearectangle;
    arearectangle.setOutlineThickness(2);
    arearectangle.setFillColor(sf::Color::Transparent);

    // Now we need to load the template images for various kinds of terrain creation.

    boolshapetemplate landshape[12]; // Basic land shape
    createlandshapetemplates(landshape);

    boolshapetemplate chainland[2]; // Mountain chain land shapes.
    createchainlandtemplates(chainland);

    boolshapetemplate smalllake[12]; // Small lake shapes
    createsmalllaketemplates(smalllake);

    boolshapetemplate largelake[10]; // Large lake shapes
    createlargelaketemplates(largelake);

    boolshapetemplate island[12]; // Small island shapes
    createislandtemplates(island);

    byteshapetemplate smudge[6]; // Smudge shapes
    createsmudgetemplates(smudge);

    byteshapetemplate smallsmudge[6]; // Small smudge shapes
    createsmallsmudgetemplates(smallsmudge);

    peaktemplate peaks; // Peak shapes.

    loadpeaktemplates(peaks);

    // Template for sections of rift lakes (I don't think these are used any more, but still)

    int riftblobsize = 40;

    vector<vector<float>> riftblob(riftblobsize + 2, vector<float>(riftblobsize + 2, 0.0f));

    createriftblob(riftblob, riftblobsize / 2);

    // Other bits

    fast_srand((long)time(0));

    NavigationState navigation;
    bool& brandnew = navigation.brandnew; // This means the program has just started and it's the first world to be generated.
    bool& newworld = navigation.newworld; // Whether to show the new world message.
    screenmodeenum& screenmode = navigation.screenmode; // This is used to keep track of which screen mode we're in.
    screenmodeenum& oldscreenmode = navigation.oldscreenmode;
    mapviewenum& mapview = navigation.mapview; // This is used to keep track of which kind of information we want the map to show.

    FileDialogState filedialogs;
    bool& loadingworld = filedialogs.loadingWorld; // This would mean that we're trying to load in a new world.
    bool& savingworld = filedialogs.savingWorld; // This would mean that we're trying to save a world.
    bool& loadingsettings = filedialogs.loadingSettings; // This would mean that we're trying to load appearance settings.
    bool& savingsettings = filedialogs.savingSettings; // This would mean that we're trying to save appearance settings.
    bool& exportingworldmaps = filedialogs.exportingWorldMaps;
    bool& exportingregionalmaps = filedialogs.exportingRegionalMaps;
    bool& exportingareamaps = filedialogs.exportingAreaMaps;
    bool& importinglandmap = filedialogs.importingLandMap;
    bool& importingseamap = filedialogs.importingSeaMap;
    bool& importingmountainsmap = filedialogs.importingMountainsMap;
    bool& importingvolcanoesmap = filedialogs.importingVolcanoesMap;
    bool& importingtemperaturemap = filedialogs.importingTemperatureMap;
    bool& importingprecipitationmap = filedialogs.importingPrecipitationMap;
    bool& importinggradientstrip = filedialogs.importingGradientStrip;

    vector<int> squareroot((MAXCRATERRADIUS* MAXCRATERRADIUS + MAXCRATERRADIUS + 1) * 24);

    for (int n = 0; n <= ((MAXCRATERRADIUS * MAXCRATERRADIUS + MAXCRATERRADIUS) * 24); n++)
        squareroot[n] = (int)sqrt(n);

    FocusState focusstate;
    int& poix = focusstate.poix; // Coordinates of the point of interest, if there is one.
    int& poiy = focusstate.poiy;
    bool& focused = focusstate.focused; // This will track whether we're focusing on one point or not.
    bool& rfocused = focusstate.regionalFocused; // Same thing, for the regional map.

    AreaSelectionState areaselection;
    int& areanwx = areaselection.northwestX;
    int& areanwy = areaselection.northwestY;
    int& areanex = areaselection.northeastX;
    int& areaney = areaselection.northeastY;
    int& areasex = areaselection.southeastX;
    int& areasey = areaselection.southeastY;
    int& areaswx = areaselection.southwestX;
    int& areaswy = areaselection.southwestY; // Corners of the export area, if there is one.
    bool& areafromregional = areaselection.fromRegional; // If this is 1 then the area screen was opened from the regional map, not the global map.

    CustomWorldUiState customworldui;
    ImportedClimateMaps importedclimatemaps;
    MapImportSettings mapimportsettings;
    int& seedentry = customworldui.seedentry; // The value currently entered into the seed box in the create world screen.
    int& newx = customworldui.newx;
    int& newy = customworldui.newy; // These are used to locate the new region.

    short regionmargin = 17; // The centre of the regional map can't be closer than this to the northern/southern edges of the global map.

    resetmapcache(globalmaps);
    resetmapcache(regionalmaps);

    PanelState panels;
    bool& showcolouroptions = panels.showColourOptions; // If this is 1 then we display the appearance preferences options.
    bool& showworldproperties = panels.showWorldProperties; // If this is 1 then we display the properties of the current world.
    bool& showglobaltemperaturechart = panels.showGlobalTemperatureChart; // If thjis is 1 then we show monthly temperatures for the selected point.
    bool& showglobalrainfallchart = panels.showGlobalRainfallChart; // If thjis is 1 then we show monthly rainfall for the selected point.
    bool& showregionaltemperaturechart = panels.showRegionalTemperatureChart; // If thjis is 1 then we show monthly temperatures for the selected point.
    bool& showregionalrainfallchart = panels.showRegionalRainfallChart; // If thjis is 1 then we show monthly rainfall for the selected point.
    bool& showsetsize = panels.showSetSize; // If this is 1 then a window is show to set the size for the custom world.
    bool& showterrainchooser = panels.showTerrainChooser; // If this is 1 then we show the panel for creating custom world terrain.
    bool& showworldgenerationoptions = panels.showWorldGenerationOptions;
    bool& showworldeditproperties = panels.showWorldEditProperties; // If this is 1 then we show the panel for editing custom world properties.
    bool& showareawarning = panels.showAreaWarning; // If this is 1 then we show a warning about too-large areas.
    bool& showabout = panels.showAbout; // If this is 1 then we display information about the program.
    bool& colourschanged = panels.coloursChanged; // If this is 1 then the colours have been changed and the maps need to be redrawn.

    ProgressPassState progresspasses;
    short& creatingworldpass = progresspasses.creatingWorld; // Tracks which pass we're on for this section, to ensure that widgets are correctly displayed.
    short& completingimportpass = progresspasses.completingImport; // Tracks which pass we're on for this section, to ensure that widgets are correctly displayed.
    short& loadingworldpass = progresspasses.loadingWorld; // Tracks which pass we're on for this section, to ensure that widgets are correctly displayed.
    short& savingworldpass = progresspasses.savingWorld; // Tracks which pass we're on for this section, to ensure that widgets are correctly displayed.
    short& exportingareapass = progresspasses.exportingArea; // Tracks which pass we're on for this section, to ensure that widgets are correctly displayed.
    short& generatingregionpass = progresspasses.generatingRegion; // Tracks which pass we're on for this section, to ensure that widgets are correctly displayed.
    short& generatingterrainpass = progresspasses.generatingTerrain; // Tracks which pass we're on for this section, to ensure that widgets are correctly displayed.

    WorldGenerationDebugState worldgenerationdebug;
    auto generationworkbenchstate = std::make_unique<GenerationSessionState>();
    GenerationSessionState& generationworkbench = *generationworkbenchstate;

    auto invalidateGlobalMaps = [&]()
    {
        resetmapcache(globalmaps);
    };

    auto invalidateRegionalMaps = [&]()
    {
        resetmapcache(regionalmaps);
    };

    auto currentglobalworld = [&]() -> planet&
    {
        if (screenmode == importscreen && generationworkbench.active && generationworkbench.ui.previewEnabled && generationworkbench.previewAvailable)
            return generationworkbench.previewWorld;

        return *world;
    };

    auto blendcolour = [](const sf::Color& low, const sf::Color& high, float factor) -> sf::Color
    {
        factor = std::clamp(factor, 0.0f, 1.0f);

        auto lerpchannel = [&](sf::Uint8 left, sf::Uint8 right)
        {
            return static_cast<sf::Uint8>(std::clamp(static_cast<int>(std::round(static_cast<float>(left) + (static_cast<float>(right) - static_cast<float>(left)) * factor)), 0, 255));
        };

        return sf::Color(lerpchannel(low.r, high.r), lerpchannel(low.g, high.g), lerpchannel(low.b, high.b));
    };

    auto tosfcolour = [](const ImVec4& colour, float alpha = 1.0f) -> sf::Color
    {
        return sf::Color(
            static_cast<sf::Uint8>(std::clamp(static_cast<int>(std::round(colour.x * 255.0f)), 0, 255)),
            static_cast<sf::Uint8>(std::clamp(static_cast<int>(std::round(colour.y * 255.0f)), 0, 255)),
            static_cast<sf::Uint8>(std::clamp(static_cast<int>(std::round(colour.z * 255.0f)), 0, 255)),
            static_cast<sf::Uint8>(std::clamp(static_cast<int>(std::round(alpha * 255.0f)), 0, 255)));
    };

    auto terrainrendercolour = [&](const planet& renderworld, int worldx, int worldy) -> sf::Color
    {
        const int altitude = renderworld.map(worldx, worldy) - renderworld.sealevel();
        return tosfcolour(samplegradientcolour(generationworkbench.ui.terrainMapGradient, altitude));
    };

    auto tectonicboundaryoverlaycolour = [&](const planet& renderworld, int worldx, int worldy, bool& visible) -> sf::Color
    {
        visible = false;
        const BoundaryType boundarytype = renderworld.tectonicboundarytype(worldx, worldy);
        const DeformingRegionType regiontype = renderworld.tectonicdeformingregiontype(worldx, worldy);

        if (boundarytype == BoundaryType::none && regiontype == DeformingRegionType::none)
            return sf::Color::Transparent;

        visible = true;
        sf::Color target(194, 177, 127);

        switch (boundarytype)
        {
        case BoundaryType::convergent:
            target = sf::Color(186, 86, 66);
            break;
        case BoundaryType::divergent:
            target = sf::Color(216, 166, 78);
            break;
        case BoundaryType::transform:
            target = sf::Color(57, 133, 142);
            break;
        case BoundaryType::passive_margin:
            target = sf::Color(194, 177, 127);
            break;
        case BoundaryType::none:
        default:
            break;
        }

        if (regiontype == DeformingRegionType::continental_rift)
            target = blendcolour(target, sf::Color(235, 193, 88), 0.45f);
        else if (regiontype == DeformingRegionType::diffuse_collision)
            target = blendcolour(target, sf::Color(160, 92, 72), 0.55f);

        return target;
    };

    auto styleworkbenchgloballayer = [&](mapviewenum view, planet& renderworld)
    {
        if (screenmode != importscreen || generationworkbench.active == false || workbenchfinished(generationworkbench))
            return;

        const GenerationStageId stageid = getcurrentgenerationstage(generationworkbench).id;
        const GenerationScratch& visiblescratch = (generationworkbench.previewAvailable && generationworkbench.ui.previewEnabled)
            ? generationworkbench.previewScratch
            : generationworkbench.committedScratch;
        const bool terrainrenderactive = isterraingenerationstage(stageid);
        const bool showplateboundaries = generationworkbench.ui.showPlateBoundaries;
        const bool showfeatureoverlay = generationworkbench.ui.showFeatureOverlay && visiblescratch.stageFeatureCellCount > 0;
        const bool inlandseaselection = stageid == GenerationStageId::inland_seas;

        if (!terrainrenderactive && !showplateboundaries && !showfeatureoverlay && !inlandseaselection)
            return;

        auto applystyle = [&](sf::Image& image)
        {
            const int imagewidth = static_cast<int>(image.getSize().x);
            const int imageheight = static_cast<int>(image.getSize().y);
            const int worldwidth = renderworld.width() + 1;
            const int worldheight = renderworld.height() + 1;

            for (int imagey = 0; imagey < imageheight; imagey++)
            {
                const int starty = (imagey * worldheight) / imageheight;
                const int endy = std::max(starty, ((imagey + 1) * worldheight - 1) / imageheight);

                for (int imagex = 0; imagex < imagewidth; imagex++)
                {
                    const int startx = (imagex * worldwidth) / imagewidth;
                    const int endx = std::max(startx, ((imagex + 1) * worldwidth - 1) / imagewidth);
                    const int samplex = std::clamp((startx + endx) / 2, 0, renderworld.width());
                    const int sampley = std::clamp((starty + endy) / 2, 0, renderworld.height());

                    if (terrainrenderactive)
                        image.setPixel(imagex, imagey, terrainrendercolour(renderworld, samplex, sampley));

                    if (inlandseaselection)
                    {
                        bool overlayvisible = false;
                        sf::Color overlaycolour = sf::Color::Transparent;

                        for (int worldy = starty; worldy <= endy && !overlayvisible; worldy++)
                        {
                            for (int worldx = startx; worldx <= endx && !overlayvisible; worldx++)
                            {
                                const int wrappedx = wrap(worldx, renderworld.width());
                                const int clampedy = std::clamp(worldy, 0, renderworld.height());
                                const int componentid = visiblescratch.inlandSeaComponentIds.empty() ? 0 : visiblescratch.inlandSeaComponentIds[wrappedx][clampedy];
                                if (componentid <= 0)
                                    continue;

                                for (const GenerationWaterComponent& component : visiblescratch.inlandSeaComponents)
                                {
                                    if (component.id != componentid)
                                        continue;

                                    overlayvisible = true;
                                    switch (component.policy)
                                    {
                                    case GenerationComponentPolicy::keep:
                                        overlaycolour = sf::Color(255, 0, 0);
                                        break;
                                    case GenerationComponentPolicy::drain:
                                        overlaycolour = sf::Color(0, 255, 0);
                                        break;
                                    case GenerationComponentPolicy::fill:
                                        overlaycolour = sf::Color(0, 0, 255);
                                        break;
                                    }
                                    break;
                                }
                            }
                        }

                        if (overlayvisible)
                            image.setPixel(imagex, imagey, overlaycolour);
                    }

                    if (showfeatureoverlay)
                    {
                        bool overlayvisible = false;

                        for (int worldy = starty; worldy <= endy && !overlayvisible; worldy++)
                        {
                            for (int worldx = startx; worldx <= endx && !overlayvisible; worldx++)
                            {
                                const int wrappedx = wrap(worldx, renderworld.width());
                                const int clampedy = std::clamp(worldy, 0, renderworld.height());
                                if (!visiblescratch.stageFeatureMask.empty() && visiblescratch.stageFeatureMask[wrappedx][clampedy] != 0)
                                {
                                    overlayvisible = true;
                                    break;
                                }
                            }
                        }

                        if (overlayvisible)
                            image.setPixel(imagex, imagey, blendcolour(image.getPixel(imagex, imagey), sf::Color(255, 0, 0), generationworkbench.ui.featureOverlayOpacity));
                    }

                    if (showplateboundaries)
                    {
                        bool overlayvisible = false;
                        sf::Color overlaycolour = sf::Color::Transparent;

                        for (int worldy = starty; worldy <= endy && overlayvisible == false; worldy++)
                        {
                            for (int worldx = startx; worldx <= endx; worldx++)
                            {
                                overlaycolour = tectonicboundaryoverlaycolour(renderworld, wrap(worldx, renderworld.width()), std::clamp(worldy, 0, renderworld.height()), overlayvisible);
                                if (overlayvisible)
                                    break;
                            }
                        }

                        if (overlayvisible)
                            image.setPixel(imagex, imagey, blendcolour(image.getPixel(imagex, imagey), overlaycolour, 0.72f));
                    }
                }
            }
        };

        maplayer& layer = getmaplayer(globalmaps, view);
        applystyle(layer.image);
        applystyle(layer.displayimage);
    };

    auto withworkbenchglobalrenderoverrides = [&](planet& renderworld, auto&& action)
    {
        const bool useworkbenchoverrides = screenmode == importscreen && generationworkbench.active;
        const bool originalshowmapoutline = renderworld.showmapoutline();
        const float originallandshading = renderworld.landshading();
        const float originallakeshading = renderworld.lakeshading();
        const float originalseashading = renderworld.seashading();
        const float originallandmarbling = renderworld.landmarbling();
        const float originallakemarbling = renderworld.lakemarbling();
        const float originalseamarbling = renderworld.seamarbling();

        if (useworkbenchoverrides)
        {
            renderworld.setshowmapoutline(generationworkbench.ui.showOceanLandContour);
            renderworld.setlandshading(0.0f);
            renderworld.setlakeshading(0.0f);
            renderworld.setseashading(0.0f);
            renderworld.setlandmarbling(0.0f);
            renderworld.setlakemarbling(0.0f);
            renderworld.setseamarbling(0.0f);
        }

        action();

        if (useworkbenchoverrides)
        {
            renderworld.setshowmapoutline(originalshowmapoutline);
            renderworld.setlandshading(originallandshading);
            renderworld.setlakeshading(originallakeshading);
            renderworld.setseashading(originalseashading);
            renderworld.setlandmarbling(originallandmarbling);
            renderworld.setlakemarbling(originallakemarbling);
            renderworld.setseamarbling(originalseamarbling);
        }
    };

    auto showGlobalMapView = [&](mapviewenum view, bool updateMiniMap = false)
    {
        mapview = view;
        planet& renderworld = currentglobalworld();
        withworkbenchglobalrenderoverrides(renderworld, [&]()
        {
            drawglobalmapimage(view, renderworld, globalmaps);
            styleworkbenchgloballayer(view, renderworld);
            updateTextureFromImage(globalmaptexture, getdisplaymapimage(globalmaps, view));
            globalmap.setTexture(globalmaptexture);

            if (updateMiniMap)
                minimap.setTexture(globalmaptexture);
        });
    };

    auto redrawGlobalRelief = [&](bool updateMiniMap = false)
    {
        showGlobalMapView(relief, updateMiniMap);
    };

    auto showRegionalMapView = [&](mapviewenum view)
    {
        mapview = view;
        applyregionalmapview(view, *world, *region, regionalmaps, regionalmaptexture, regionalmap);
    };

    auto redrawRegionalRelief = [&]()
    {
        showRegionalMapView(relief);
    };

    auto rebuildAllGlobalMaps = [&]()
    {
        planet& renderworld = currentglobalworld();
        withworkbenchglobalrenderoverrides(renderworld, [&]()
        {
            drawallglobalmapimages(renderworld, globalmaps);
        });
    };

    auto rebuildAllRegionalMaps = [&]()
    {
        drawallregionalmapimages(*world, *region, regionalmaps);
    };

    auto socialoptionsfromdebug = [&]() -> SocialGenerationOptions
    {
        return makesocialoptions(
            worldgenerationdebug.options.socialEnabled,
            worldgenerationdebug.options.socialMode,
            worldgenerationdebug.options.usePrehistory,
            worldgenerationdebug.options.historyYears);
    };

    auto applytectonicsoptionsfromdebug = [&](GenerationWorkbenchUiState& ui)
    {
        ui.tectonicDurationMyr = worldgenerationdebug.options.plateTectonicsDurationMyr;
        ui.tectonicPlateCount = worldgenerationdebug.options.plateTectonicsPlateCount;
        ui.tectonicDeltaTimeMyr = worldgenerationdebug.options.plateTectonicsDeltaTimeMyr;
    };

    auto runworkbenchaction = [&](auto&& action)
    {
        action();
    };

    float linespace = 8.0f; // Gap between groups of buttons.

    string& filepathname = filedialogs.filepathname;
    string& filepath = filedialogs.filepath;

    vector<vector<bool>> OKmountains(ARRAYWIDTH, vector<bool>(ARRAYHEIGHT, 0)); // This will track mountains that have been imported by the user and which therefore should not be altered to take account of gravity.

    // Now world property variables that can be directly manipulated in the world settings edit panel.

    WorldPropertyControls worldpropertycontrols = makeworldpropertycontrols(*world);
    int& currentsize = worldpropertycontrols.size;
    float& currentgravity = worldpropertycontrols.gravity;
    float& currentlunar = worldpropertycontrols.lunar;
    float& currenteccentricity = worldpropertycontrols.eccentricity;
    int& currentperihelion = worldpropertycontrols.perihelion;
    int& currentrotation = worldpropertycontrols.rotation;
    float& currenttilt = worldpropertycontrols.tilt;
    float& currenttempdecrease = worldpropertycontrols.tempdecrease;
    int& currentnorthpolaradjust = worldpropertycontrols.northpolaradjust;
    int& currentsouthpolaradjust = worldpropertycontrols.southpolaradjust;
    int& currentaveragetemp = worldpropertycontrols.averagetemp;
    float& currentwaterpickup = worldpropertycontrols.waterpickup;
    int& currentglacialtemp = worldpropertycontrols.glacialtemp;
    bool& currentrivers = worldpropertycontrols.rivers; // This one controls whether or not rivers will be calculated for the custom world.
    bool& currentlakes = worldpropertycontrols.lakes; // This one controls whether or not lakes will be calculated for the custom world.
    bool& currentdeltas = worldpropertycontrols.deltas; // This one controls whether or not river deltas will be calculated for the custom world.
    bool& currentsocialenabled = worldpropertycontrols.socialEnabled;
    SocialGenerationOptions::Mode& currentsocialmode = worldpropertycontrols.socialMode;
    bool& currentuseprehistory = worldpropertycontrols.usePrehistory;
    int& currenthistoryyears = worldpropertycontrols.historyYears;

    auto socialoptionsfromworldcontrols = [&]() -> SocialGenerationOptions
    {
        return makesocialoptions(
            currentsocialenabled,
            currentsocialmode,
            currentuseprehistory,
            currenthistoryyears);
    };

    auto startproceduralworkbenchsession = [&]()
    {
        initialiseworld(*world);
        world->clear();
        resetworkbenchsession(generationworkbench);
        clearimportedclimatemaps(importedclimatemaps);
        mapimportsettings = MapImportSettings{};

        changeworldproperties(*world);
        applyworlddimensions(*world, worldpropertycontrols);
        world->setseed(seedentry);
        syncworldpropertycontrols(*world, worldpropertycontrols);
        currentrivers = true;
        currentlakes = true;
        currentdeltas = true;

        resizeglobaldisplayforworld(*world, globaldisplay, highlightdisplay);
        fast_srand(world->seed());

        initializeproceduralworkbenchsession(generationworkbench, *world, worldgenerationdebug.options.plateTectonicsDurationMyr, worldgenerationdebug.options.plateTectonicsPlateCount);
        applytectonicsoptionsfromdebug(generationworkbench.ui);

        screenmode = importscreen;
        brandnew = 0;

        GenerationExecutionContext context;
        context.dorivers = currentrivers;
        context.dolakes = currentlakes;
        context.dodeltas = currentdeltas;
        context.socialoptions = socialoptionsfromworldcontrols();
        context.landshape = landshape;
        context.chainland = chainland;
        context.smalllake = smalllake;
        context.largelake = largelake;
        context.okmountains = &OKmountains;

        runworkbenchaction([&]()
        {
            previewcurrentgenerationstage(generationworkbench, *world, context, nullptr);
        });

        invalidateGlobalMaps();
        showGlobalMapView(preferredworkbenchmapview(getcurrentgenerationstage(generationworkbench).id));
    };

    AppearanceSettings appearance = makeappearancesettings(*world);
    ImVec4& oceancolour = appearance.oceancolour;
    ImVec4& deepoceancolour = appearance.deepoceancolour;
    ImVec4& basecolour = appearance.basecolour;
    ImVec4& grasscolour = appearance.grasscolour;
    ImVec4& basetempcolour = appearance.basetempcolour;
    ImVec4& highbasecolour = appearance.highbasecolour;
    ImVec4& desertcolour = appearance.desertcolour;
    ImVec4& highdesertcolour = appearance.highdesertcolour;
    ImVec4& colddesertcolour = appearance.colddesertcolour;
    ImVec4& eqtundracolour = appearance.eqtundracolour;
    ImVec4& tundracolour = appearance.tundracolour;
    ImVec4& coldcolour = appearance.coldcolour;
    ImVec4& seaicecolour = appearance.seaicecolour;
    ImVec4& glaciercolour = appearance.glaciercolour;
    ImVec4& saltpancolour = appearance.saltpancolour;
    ImVec4& ergcolour = appearance.ergcolour;
    ImVec4& wetlandscolour = appearance.wetlandscolour;
    ImVec4& lakecolour = appearance.lakecolour;
    ImVec4& rivercolour = appearance.rivercolour;
    ImVec4& sandcolour = appearance.sandcolour;
    ImVec4& mudcolour = appearance.mudcolour;
    ImVec4& shinglecolour = appearance.shinglecolour;
    ImVec4& mangrovecolour = appearance.mangrovecolour;
    ImVec4& highlightcolour = appearance.highlightcolour;
    ImVec4& outlinecolour = appearance.outlinecolour;
    ImVec4& elevationlowcolour = appearance.elevationlowcolour;
    ImVec4& elevationhighcolour = appearance.elevationhighcolour;
    ImVec4& temperaturecoldcolour = appearance.temperaturecoldcolour;
    ImVec4& temperaturetemperatecolour = appearance.temperaturetemperatecolour;
    ImVec4& temperaturehotcolour = appearance.temperaturehotcolour;
    ImVec4& precipitationdrycolour = appearance.precipitationdrycolour;
    ImVec4& precipitationwetcolour = appearance.precipitationwetcolour;
    float& shadingland = appearance.shadingland;
    float& shadinglake = appearance.shadinglake;
    float& shadingsea = appearance.shadingsea;
    float& marblingland = appearance.marblingland;
    float& marblinglake = appearance.marblinglake;
    float& marblingsea = appearance.marblingsea;
    int& globalriversentry = appearance.globalriversentry;
    int& regionalriversentry = appearance.regionalriversentry;
    int& shadingdir = appearance.shadingdir;
    int& snowchange = appearance.snowchange;
    int& seaiceappearance = appearance.seaiceappearance;
    bool& showmapoutline = appearance.showmapoutline;
    bool& colourcliffs = appearance.colourcliffs;
    bool& mangroves = appearance.mangroves;
    array<int, MAPGRADIENTTYPECOUNT> selectedgradientstops = {};

    auto showdeferredworkwindow = [&](short& pass, const char* windowTitle, const char* message, ImVec2 position, ImVec2 size) -> bool
    {
        if (pass >= 10)
            return true;

        ImGui::SetNextWindowPos(position, ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(size, ImGuiCond_FirstUseEver);
        ImGui::Begin(windowTitle, NULL, window_flags);
        ImGui::Text(message);
        ImGui::End();

        pass++;
        return false;
    };

    auto preparecustomworldgeneration = [&]()
    {
        initialiseworld(*world);
        world->clear();
        world->setsize(currentsize);
        applyworlddimensions(*world, worldpropertycontrols);

        resizeglobaldisplayforworld(*world, globaldisplay, highlightdisplay);

        world->setseed(0 - createrandomseednumber());

        updatereport("Generating custom world terrain from seed: " + to_string(world->seed()) + ":");
        updatereport("");
        begintimedreporting();

        invalidateGlobalMaps();
        fast_srand(world->seed());
    };

    auto finishcustomworldgeneration = [&](short& pass)
    {
        endtimedreporting();
        resetmapcache(globalmaps);

        mapview = relief;
        showGlobalMapView(mapview);

        pass = 0;
        screenmode = importscreen;
    };

    sf::Clock deltaClock;
    while (window.isOpen())
    {
        sf::Event event;
        while (window.pollEvent(event))
        {
            ImGui::SFML::ProcessEvent(window, event);

            if (event.type == sf::Event::Closed)
                window.close();
        }

        ImGuiIO& io = ImGui::GetIO();

        ImGui::SFML::Update(window, deltaClock.restart());

        if (screenmode == movingtoglobalmapscreen)
            screenmode = globalmapscreen;

        ImGui::PushFont(font1);

        //ImGui::ShowDemoWindow();

        // First, draw the GUI.

        // Create world screen

#include "app/screens/create_world_screen.inc"

        // Creating world screen

#include "app/screens/creating_world_screen.inc"

        // Global map screen

#include "app/screens/global_map_screen.inc"

        // Regional map screen

#include "app/screens/regional_map_screen.inc"

        // Area export screen (this is where the user selects the area to export)

#include "app/screens/export_area_screen.inc"

        // Custom world screen

#include "app/screens/import_screen.inc"

        // These screens all display a "Please wait" message ten times (for some reason doing it once or twice doesn't actually display it) and then do something time-consuming.
#include "app/screens/completing_import_screen.inc"

#include "app/screens/loading_world_screen.inc"

#include "app/screens/saving_world_screen.inc"

#include "app/screens/exporting_area_screen.inc"

#include "app/screens/generating_region_screen.inc"

#include "app/screens/load_failure_screen.inc"

#include "app/screens/settings_load_failure_screen.inc"

        window.clear();

        // Colour options, if being shown

#include "app/screens/appearance_panel.inc"

        // World properties window, if being shown

#include "app/screens/world_properties_panel.inc"

        // World edit properties screen, if being shown

        drawworldeditpropertieswindow(main_viewport, window_flags, worldpropertycontrols, degree, showworldeditproperties);

        // Global temperature chart, if being shown

        if (showglobaltemperaturechart)
            drawglobaltemperaturechartwindow(main_viewport, window_flags, *world, poix, poiy);

        // Global precipitation chart, if being shown

        if (showglobalrainfallchart)
            drawglobalrainfallchartwindow(main_viewport, window_flags, *world, poix, poiy);

        // Regional temperature chart, if being shown

        if (showregionaltemperaturechart)
            drawregionaltemperaturechartwindow(main_viewport, window_flags, *world, *region, poix, poiy);

        // Regional precipitation chart, if being shown

        if (showregionalrainfallchart)
            drawregionalrainfallchartwindow(main_viewport, window_flags, *region, poix, poiy);

        // Now draw the graphical elements.

        if (screenmode == globalmapscreen || screenmode == exportareascreen || screenmode == importscreen)
            window.draw(globalmap);

        if (screenmode == exportareascreen) // Area selection rectangle.
        {
            if (areaswx != -1)
            {
                float mult = ((float)world->width() + 1.0f) / (float)DISPLAYMAPSIZEX;

                float sizex = ((float)areanex - (float)areanwx) / mult;
                float sizey = ((float)areaswy - (float)areanwy) / mult;

                arearectangle.setSize(sf::Vector2f(sizex, sizey));
                arearectangle.setOutlineColor(sf::Color(world->highlight1(), world->highlight2(), world->highlight3()));

                sf::Vector2f position;

                position.x = (float)areanwx / mult + (float)globalmapxpos;
                position.y = (float)areanwy / mult + (float)globalmapypos;

                arearectangle.setPosition(position);

                window.draw(arearectangle);
            }
        }

        if (screenmode == regionalmapscreen || screenmode == generatingregionscreen)
        {
            window.draw(regionalmap);

            window.draw(minimap);

            float mult = (float)DISPLAYMAPSIZEX / ((float)world->width() + 1.0f);
            mult = mult * 0.5f;

            float posx = (float)minimapxpos + (float)region->centrex() * mult;
            float posy = (float)minimapypos + (float)region->centrey() * mult;

            minihighlight.setPosition(sf::Vector2f(posx, posy));

            if (screenmode != generatingregionscreen)
            {
                window.draw(minihighlight);
            }
        }

        if (screenmode == globalmapscreen || screenmode == regionalmapscreen)
        {
            if (focused == 1)
                window.draw(highlight);
        }

        // Now draw file dialogues if needed

        handlefiledialog(loadingworld, filepathname, filepath, [&](const std::string&, const std::string&)
        {
            resetworkbenchsession(generationworkbench);
            screenmode = loadingworldscreen;
        });

        handlefiledialog(savingworld, filepathname, filepath, [&](const std::string&, const std::string&)
        {
            screenmode = savingworldscreen;
        });

        handlefiledialog(loadingsettings, filepathname, filepath, [&](const std::string& selectedpath, const std::string&)
        {
            bool found = loadsettings(*world, selectedpath);

            if (found == 1)
            {
                syncappearancesettings(*world, appearance);
                refreshhighlightdisplay(*world, highlightdisplay);

                colourschanged = 1;
            }
            else
            {
                oldscreenmode = screenmode;
                screenmode = settingsloadfailure;
            }
        });

        handlefiledialog(savingsettings, filepathname, filepath, [&](const std::string& selectedpath, const std::string&)
        {
            savesettings(*world, selectedpath);
        });

        handlefiledialog(exportingworldmaps, filepathname, filepath, [&](const std::string& selectedpath, const std::string&)
        {
            rebuildAllGlobalMaps();

            std::string exportpath = selectedpath;
            exportpath.resize(exportpath.size() - 4);

            for (const mapviewdefinition& definition : allmapviewdefinitions)
                getmapimage(globalmaps, definition.view).saveToFile(exportpath + " " + definition.exportstem + ".png");
        });

        handlefiledialog(exportingregionalmaps, filepathname, filepath, [&](const std::string& selectedpath, const std::string&)
        {
            rebuildAllRegionalMaps();

            std::string exportpath = selectedpath;
            exportpath.resize(exportpath.size() - 4);

            for (const mapviewdefinition& definition : allmapviewdefinitions)
                getmapimage(regionalmaps, definition.view).saveToFile(exportpath + " " + definition.exportstem + ".png");
        });

        handlefiledialog(exportingareamaps, filepathname, filepath, [&](const std::string&, const std::string&)
        {
            screenmode = exportingareascreen;
        });

        // These sections are for loading in maps in the import screen.

        handlefiledialog(importinglandmap, filepathname, filepath, [&](const std::string& selectedpath, const std::string&)
        {
            if (importlandheightmap(*world, selectedpath, &mapimportsettings))
            {
                resetworkbenchsession(generationworkbench);
                invalidateGlobalMaps();
                showGlobalMapView(relief);
            }
        });

        handlefiledialog(importingseamap, filepathname, filepath, [&](const std::string& selectedpath, const std::string&)
        {
            if (importseadepthmap(*world, selectedpath, &mapimportsettings))
            {
                resetworkbenchsession(generationworkbench);
                invalidateGlobalMaps();
                showGlobalMapView(relief);
            }
        });

        handlefiledialog(importingmountainsmap, filepathname, filepath, [&](const std::string& selectedpath, const std::string&)
        {
            if (importmountainmap(*world, selectedpath, OKmountains, &mapimportsettings))
            {
                resetworkbenchsession(generationworkbench);
                invalidateGlobalMaps();
                showGlobalMapView(relief);
            }
        });

        handlefiledialog(importingvolcanoesmap, filepathname, filepath, [&](const std::string& selectedpath, const std::string&)
        {
            if (importvolcanomap(*world, selectedpath, &mapimportsettings))
            {
                resetworkbenchsession(generationworkbench);
                invalidateGlobalMaps();
                showGlobalMapView(relief);
            }
        });

        handlefiledialog(importingtemperaturemap, filepathname, filepath, [&](const std::string& selectedpath, const std::string&)
        {
            if (importtemperaturemap(*world, selectedpath, importedclimatemaps, &mapimportsettings))
            {
                resetworkbenchsession(generationworkbench);
                invalidateGlobalMaps();
                showGlobalMapView(temperature);
            }
        });

        handlefiledialog(importingprecipitationmap, filepathname, filepath, [&](const std::string& selectedpath, const std::string&)
        {
            if (importprecipitationmap(*world, selectedpath, importedclimatemaps, &mapimportsettings))
            {
                resetworkbenchsession(generationworkbench);
                invalidateGlobalMaps();
                showGlobalMapView(precipitation);
            }
        });

        handlefiledialog(importinggradientstrip, filepathname, filepath, [&](const std::string& selectedpath, const std::string&)
        {
            mapimportsettings.gradientStripPath = selectedpath;
        });

        // Set size window, for custom worlds.

        if (drawcustomworldsizewindow(main_viewport, window_flags, showsetsize, brandnew, worldpropertycontrols))
        {
            initialiseworld(*world);
            world->clear();
            resetworkbenchsession(generationworkbench);
            clearimportedclimatemaps(importedclimatemaps);
            mapimportsettings = MapImportSettings{};

            world->setsize(currentsize);
            applyworlddimensions(*world, worldpropertycontrols);

            resizeglobaldisplayforworld(*world, globaldisplay, highlightdisplay);
            world->setgravity(getdefaultgravityforsize(currentsize));

            syncworldpropertycontrols(*world, worldpropertycontrols);
            currentrivers = true;
            currentlakes = true;
            currentdeltas = true;

            int width = world->width();
            int height = world->height();
            int val = world->sealevel() - 5000;

            for (int i = 0; i <= width; i++)
            {
                for (int j = 0; j <= height; j++)
                    world->setnom(i, j, val);
            }

            mapview = relief;

            resetmapcache(globalmaps);

            showGlobalMapView(mapview);

            world->setseed(0 - createrandomseednumber());

            screenmode = importscreen;
        }

        // Window for creating custom worlds.

        if (drawterrainchooserwindow(main_viewport, window_flags, showterrainchooser))
        {
            resetworkbenchsession(generationworkbench);
            initializeproceduralworkbenchsession(generationworkbench, *world, worldgenerationdebug.options.plateTectonicsDurationMyr, worldgenerationdebug.options.plateTectonicsPlateCount);
            applytectonicsoptionsfromdebug(generationworkbench.ui);
            previewcurrentgenerationstage(generationworkbench, *world, GenerationExecutionContext{}, nullptr);
            invalidateGlobalMaps();
            showGlobalMapView(preferredworkbenchmapview(getcurrentgenerationstage(generationworkbench).id));
        }

        // Warning window for over-large area maps.

        if (drawareawarningwindow(main_viewport, window_flags, showareawarning))
        {
            openFileDialog(".png");
            exportingareamaps = 1;
        }

        // Window for displaying information about the program.

        drawaboutwindow(main_viewport, window_flags, showabout, currentversion);

        ImGui::PopFont();

        ImGui::SFML::Render(window);
        window.display();

        // Now update the colours if necessary.

        const AppearanceChangeFlags appearanceChanges = getappearancechanges(*world, appearance);
        const bool shouldApplyAppearance = showcolouroptions == 1 && (colourschanged == 1 || appearanceChanges.any());
        const bool shouldRedrawAppearanceMaps = colourschanged == 1 || appearanceChanges.mapAppearanceChanged;

        if (shouldApplyAppearance)
        {
            applyappearancesettings(*world, appearance);

            if (appearanceChanges.highlightChanged)
                refreshhighlightdisplay(*world, highlightdisplay);

            if (shouldRedrawAppearanceMaps)
            {
                invalidateGlobalMaps();
                showGlobalMapView(mapview, true);

                if (screenmode == regionalmapscreen)
                {
                    invalidateRegionalMaps();
                    showRegionalMapView(mapview);
                }
            }

            colourschanged = 0;
        }
    }

    ImGui::SFML::Shutdown();

    return 0;
}
