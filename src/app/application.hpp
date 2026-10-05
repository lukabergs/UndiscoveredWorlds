#pragma once

// The desktop application's state and its screens. One Application is built by
// main(); run() is the frame loop, and each screen or panel is a member function
// defined in app/screens/<name>.cpp.

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
#include "world/planet_bundle.hpp"
#include "pipeline/world_processes.hpp"
#include "pipeline/world_stages.hpp"
#include "simulations/terrain/import_elements.hpp"
#include "pipeline/physical_layers.hpp"
#include "simulations/society/social_generation.hpp"
#include "app/regional/regional_map.hpp"
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

// A saved world, or the manifest.json of a run bundle.
inline constexpr const char* kWorldLoadFilter = ".uww,.json";
inline constexpr const char* kRasterImportFilter = ".png,.tif,.tiff";
inline constexpr const char* kGradientStripFilter = ".png";
inline constexpr const char* kVolcanoImportFilter = ".png";

const char* crustclasslabel(CrustClass crustclass);
const char* boundarytypelabel(BoundaryType boundarytype);
const char* deformingregiontypelabel(DeformingRegionType regiontype);
const char* geologicregimelabel(GeologicRegime regime);
void drawboundarysegmentsummary(const char* label, int requestedid, const TectonicBoundarySegment* segment);
void drawdeformingregionsummary(const char* label, int requestedid, const TectonicDeformingRegion* region);
void drawmineraldepositsummary(const planet& world, int x, int y);
bool gettectonicinspectionpoint(const planet& world, const region& region, screenmodeenum screenmode, bool focused, int poix, int poiy, int& x, int& y);
int createrandomseednumber();

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

void initglobalmapdisplay(planet& world, GlobalMapDisplay& display);
void inithighlightdisplay(planet& world, HighlightDisplay& display);
void initregionalmapdisplay(RegionalMapDisplay& display);
void refreshminihighlightdisplay(HighlightDisplay& display);
void refreshhighlightdisplay(planet& world, HighlightDisplay& display);
void resizeglobaldisplayforworld(planet& world, GlobalMapDisplay& globaldisplay, HighlightDisplay& highlightdisplay);
void configureimguistyle();

struct Application
{
    Application();
    void run();

    // Screens: each draws its windows and acts on its input while its screen mode is current.
    void updatecreateworldscreen();
    void updatecreatingworldscreen();
    void updateglobalmapscreen();
    void updateregionalmapscreen();
    void updateexportareascreen();
    void updateimportscreen();
    void updatecompletingimportscreen();
    void updateloadingworldscreen();
    bool loadworldfile(const std::string& filename);
    bool runprocess(const char* id);
    void updatesavingworldscreen();
    void updateexportingareascreen();
    void updategeneratingregionscreen();
    void updateloadfailurescreen();
    void updatesettingsloadfailurescreen();
    void updateappearancepanel();
    void updateworldpropertiespanel();

    // The rest of a frame, in the order run() calls them.
    void drawchartwindows();
    void drawmaps();
    void handlefiledialogs();
    void drawdialogwindows();
    void applyappearancechanges();

    static void openFileDialog(const char* filter);
    void invalidateGlobalMaps();
    void invalidateRegionalMaps();
    planet& currentglobalworld();
    static sf::Color blendcolour(const sf::Color& low, const sf::Color& high, float factor);
    static sf::Color tosfcolour(const ImVec4& colour, float alpha = 1.0f);
    sf::Color terrainrendercolour(const planet& renderworld, int worldx, int worldy);
    sf::Color tectonicboundaryoverlaycolour(const planet& renderworld, int worldx, int worldy, bool& visible);
    void styleworkbenchgloballayer(mapviewenum view, planet& renderworld);
    template <typename Action> void withworkbenchglobalrenderoverrides(planet& renderworld, Action&& action);
    void showGlobalMapView(mapviewenum view, bool updateMiniMap = false);
    void redrawGlobalRelief(bool updateMiniMap = false);
    void showRegionalMapView(mapviewenum view);
    void redrawRegionalRelief();
    void rebuildAllGlobalMaps();
    void rebuildAllRegionalMaps();
    SocialGenerationOptions socialoptionsfromdebug();
    void applytectonicsoptionsfromdebug(GenerationWorkbenchUiState& ui);
    template <typename Action> void runworkbenchaction(Action&& action) { action(); }
    SocialGenerationOptions socialoptionsfromworldcontrols();
    void startproceduralworkbenchsession();
    bool showdeferredworkwindow(short& pass, const char* windowTitle, const char* message, ImVec2 position, ImVec2 size);
    void preparecustomworldgeneration();
    void finishcustomworldgeneration(short& pass);

    float currentversion = 1.0f;

    // The window.

    int scwidth = 1224;
    int scheight = 768;

    sf::RenderWindow window{ sf::VideoMode(scwidth, scheight), "Undiscovered Worlds" };

    const ImGuiViewport* main_viewport = nullptr;
    ImFont* font1 = nullptr;

    ImGuiWindowFlags window_flags = ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoResize;

    string degree = "\xC2\xB0"; // 0xC2 0xB0 (c2b0) ï¿½
    string cube = "\xC2\xB3"; // 0xC2 0xB3 (c2b3) ï¿½

    // The main world object.

    std::unique_ptr<planet> world = std::make_unique<planet>();

    int width = 0;
    int height = 0;

    // The region object.

    std::unique_ptr<class region> region = std::make_unique<class region>();

    int regionalmapimagewidth = 0;
    int regionalmapimageheight = 0;

    // Create display/state bundles for the global map, highlights, and regional map.

    GlobalMapDisplay globaldisplay;
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
    int& highlightsize = highlightdisplay.size;
    sf::Image& highlightimage = highlightdisplay.image;
    sf::Texture& highlighttexture = highlightdisplay.texture;
    sf::Sprite& highlight = highlightdisplay.sprite;
    int& minihighlightsize = highlightdisplay.miniSize;
    sf::Image& minihighlightimage = highlightdisplay.miniImage;
    sf::Texture& minihighlighttexture = highlightdisplay.miniTexture;
    sf::Sprite& minihighlight = highlightdisplay.miniSprite;

    RegionalMapDisplay regionaldisplay;
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

    // Now we need to load the template images for various kinds of terrain creation.

    boolshapetemplate landshape[12]; // Basic land shape

    boolshapetemplate chainland[2]; // Mountain chain land shapes.

    boolshapetemplate smalllake[12]; // Small lake shapes

    boolshapetemplate largelake[10]; // Large lake shapes

    boolshapetemplate island[12]; // Small island shapes

    byteshapetemplate smudge[6]; // Smudge shapes

    byteshapetemplate smallsmudge[6]; // Small smudge shapes

    peaktemplate peaks; // Peak shapes.

    // Template for sections of rift lakes (I don't think these are used any more, but still)

    int riftblobsize = 40;

    vector<vector<float>> riftblob = vector<vector<float>>(riftblobsize + 2, vector<float>(riftblobsize + 2, 0.0f));

    // Other bits

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

    vector<int> squareroot = vector<int>((MAXCRATERRADIUS * MAXCRATERRADIUS + MAXCRATERRADIUS + 1) * 24);

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
    std::unique_ptr<GenerationSessionState> generationworkbenchstate = std::make_unique<GenerationSessionState>();
    GenerationSessionState& generationworkbench = *generationworkbenchstate;

    float linespace = 8.0f; // Gap between groups of buttons.

    string& filepathname = filedialogs.filepathname;
    string& filepath = filedialogs.filepath;

    vector<vector<bool>> OKmountains = vector<vector<bool>>(ARRAYWIDTH, vector<bool>(ARRAYHEIGHT, 0)); // This will track mountains that have been imported by the user and which therefore should not be altered to take account of gravity.

    // Now world property variables that can be directly manipulated in the world settings edit panel.

    WorldPropertyControls worldpropertycontrols;
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

    AppearanceSettings appearance;
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
};

template <typename Action>
void Application::withworkbenchglobalrenderoverrides(planet& renderworld, Action&& action)
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
}
