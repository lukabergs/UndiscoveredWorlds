
#include "imgui.h"
#include "imgui-SFML.h"
#include "ImGuiFileDialog.h"
#include "create_world_screen.hpp"

#pragma comment(lib, "urlmon.lib")

#include <urlmon.h>
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
#include <chrono>
#include <thread>
#include <queue>
#include <stdint.h>
#include <SFML/Graphics.hpp>
#include <SFML/Graphics/CircleShape.hpp>
#include <SFML/Graphics/RenderWindow.hpp>
#include <SFML/System/Clock.hpp>
#include <SFML/Window/Event.hpp>

#include "classes.hpp"
#include "planet.hpp"
#include "physical_layers.hpp"
#include "region.hpp"
#include "app_dialogs.hpp"
#include "app_environment.hpp"
#include "app_state.hpp"
#include "app_windows.hpp"
#include "appearance_settings.hpp"
#include "functions.hpp"
#include "map_appearance.hpp"
#include "map_imports.hpp"
#include "ui_charts.hpp"
#include "ui_panels.hpp"
#include "world_property_controls.hpp"

#define getURL URLOpenBlockingStreamA

// nodiscard attribute
#pragma warning (disable: 4834)

#define REGIONALCREATIONSTEPS 84

#define GLOBALTERRAINCREATIONSTEPS1 26
#define GLOBALTERRAINCREATIONSTEPS2 31
#define GLOBALCLIMATECREATIONSTEPS 29

#define REGIONALTILEWIDTH 32
#define REGIONALTILEHEIGHT 32

using namespace std;

// Random number generator. From https://stackoverflow.com/questions/26237419/faster-than-rand
// It uses a global variable. But this is the only one in the program, honest!

static long g_seed = 1;

namespace
{
struct CommandLineGenerationOptions
{
    bool run = false;
    bool showHelp = false;
    bool earthClimateBenchmark = false;
    bool printClimateRelativeError = false;
    bool useFastLEMMountains = false;
    bool usePlateTectonicsSimulation = false;
    bool rivers = true;
    bool lakes = true;
    bool deltas = true;
    bool logToProfilingWorkbook = true;
    bool appendClimateWorkbook = false;
    bool hasSeed = false;
    bool hasSavePath = false;
    bool hasReferencePath = false;
    bool hasImportLandPath = false;
    bool hasImportSeaPath = false;
    long seed = 0;
    int plateTectonicsCycleCount = 4;
    string savePath;
    string referencePath;
    string importLandPath;
    string importSeaPath;
};

string narrowargument(const wchar_t* value)
{
    return filesystem::path(value).string();
}

long defaultcommandlineseed()
{
    const long long ticks = chrono::high_resolution_clock::now().time_since_epoch().count();
    const long long normalized = llabs(ticks % 90000000ll) + 10000000ll;
    return static_cast<long>(0 - normalized);
}

void printcommandlineusage()
{
    cout << "Undiscovered Worlds command line usage:\n";
    cout << "  --generate-world [options]\n";
    cout << "Options:\n";
    cout << "  --seed <number>         Use a fixed world seed.\n";
    cout << "  --save <path>           Save the generated world as a .uww file.\n";
    cout << "  --fastlem               Enable FastLEM mountains.\n";
    cout << "  --tectonics             Enable plate tectonics simulation.\n";
    cout << "  --plate-cycles <count>  Set tectonic cycle count.\n";
    cout << "  --reference-precip <csv> Compare against a precipitation reference grid.\n";
    cout << "  --import-land <png>     Import a land height map instead of generating terrain.\n";
    cout << "  --import-sea <png>      Import a sea depth map instead of generating terrain.\n";
    cout << "  --earth-climate-benchmark Run the imported Earth benchmark workflow.\n";
    cout << "  --print-climate-relative-error Print per-climate relative error against the Earth benchmark counts.\n";
    cout << "  --no-climate-workbook   Skip appending benchmark counts to climate.xlsx.\n";
    cout << "  --no-rivers             Disable river generation.\n";
    cout << "  --no-lakes              Disable lake generation.\n";
    cout << "  --no-deltas             Disable delta generation.\n";
    cout << "  --no-profile-log        Skip profiling workbook logging.\n";
    cout << "  --help                  Show this message.\n";
}

bool parsecommandlineoptions(CommandLineGenerationOptions& options)
{
    int argc = 0;
    LPWSTR* argv = CommandLineToArgvW(GetCommandLineW(), &argc);

    if (argv == nullptr)
        return true;

    const auto cleanup = [&]()
    {
        LocalFree(argv);
    };

    for (int index = 1; index < argc; index++)
    {
        const string argument = narrowargument(argv[index]);

        if (argument == "--help" || argument == "-h" || argument == "/?")
        {
            options.showHelp = true;
            continue;
        }

        if (argument == "--generate-world")
        {
            options.run = true;
            continue;
        }

        if (argument == "--earth-climate-benchmark")
        {
            options.run = true;
            options.earthClimateBenchmark = true;
            options.appendClimateWorkbook = true;
            options.logToProfilingWorkbook = false;
            continue;
        }

        if (argument == "--print-climate-relative-error")
        {
            options.printClimateRelativeError = true;
            continue;
        }

        if (argument == "--fastlem")
        {
            options.useFastLEMMountains = true;
            continue;
        }

        if (argument == "--tectonics")
        {
            options.usePlateTectonicsSimulation = true;
            continue;
        }

        if (argument == "--no-rivers")
        {
            options.rivers = false;
            continue;
        }

        if (argument == "--no-lakes")
        {
            options.lakes = false;
            continue;
        }

        if (argument == "--no-deltas")
        {
            options.deltas = false;
            continue;
        }

        if (argument == "--no-profile-log")
        {
            options.logToProfilingWorkbook = false;
            continue;
        }

        if (argument == "--no-climate-workbook")
        {
            options.appendClimateWorkbook = false;
            continue;
        }

        if (argument == "--seed" || argument == "--save" || argument == "--plate-cycles" || argument == "--reference-precip" || argument == "--import-land" || argument == "--import-sea")
        {
            if (index + 1 >= argc)
            {
                cerr << "Missing value for " << argument << '\n';
                cleanup();
                return false;
            }

            const string value = narrowargument(argv[++index]);

            try
            {
                if (argument == "--seed")
                {
                    options.seed = stol(value);
                    options.hasSeed = true;
                }
                else if (argument == "--plate-cycles")
                {
                    options.plateTectonicsCycleCount = max(1, stoi(value));
                }
                else if (argument == "--reference-precip")
                {
                    options.referencePath = value;
                    options.hasReferencePath = true;
                }
                else if (argument == "--import-land")
                {
                    options.importLandPath = value;
                    options.hasImportLandPath = true;
                }
                else if (argument == "--import-sea")
                {
                    options.importSeaPath = value;
                    options.hasImportSeaPath = true;
                }
                else
                {
                    options.savePath = value;
                    options.hasSavePath = true;
                }
            }
            catch (const exception&)
            {
                cerr << "Invalid value for " << argument << ": " << value << '\n';
                cleanup();
                return false;
            }

            continue;
        }

        cerr << "Unknown argument: " << argument << '\n';
        cleanup();
        return false;
    }

    cleanup();
    return true;
}

filesystem::path commandlinevalidationdirectory(long seed)
{
    filesystem::path outputroot = getappenvironment().profilingWorkbookPath.parent_path();

    if (outputroot.empty())
        outputroot = filesystem::current_path();

    return outputroot / "validation" / ("seed_" + to_string(seed));
}

void completeimportedworldgeneration(planet& world, bool dorivers, bool dolakes, bool dodeltas, bool appendclimateworkbook, boolshapetemplate smalllake[], boolshapetemplate largelake[], boolshapetemplate landshape[], vector<vector<bool>>& okmountains, const ImportedClimateMaps* importedclimate = nullptr)
{
    const int width = world.width();
    const int height = world.height();

    vector<vector<int>> mountaindrainage(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));
    vector<vector<bool>> shelves(ARRAYWIDTH, vector<bool>(ARRAYHEIGHT, 0));

    world.cleartectonicprovenance();

    world.setmaxelevation(200000);

    updatereport("Raising mountain bases");

    raisemountainbases(world, mountaindrainage, okmountains);

    getlandandseatotals(world);

    if (world.seatotal() > 10)
    {
        updatereport("Filling depressions");

        depressionfill(world);
        addlandnoise(world);
        depressionfill(world);

        updatereport("Adjusting coastlines");

        for (int n = 1; n <= 2; n++)
            normalisecoasts(world, 13, 11, 4);

        clamp(world);

        updatereport("Checking islands");

        checkislands(world);
    }

    updatereport("Creating roughness map");

    vector<vector<int>> roughness(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    createfractal(roughness, width, height, 8, 0.2f, 0.6f, 1, world.maxelevation(), 0, 0);

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
            world.setroughness(i, j, static_cast<float>(roughness[i][j]));
    }

    generateglobalclimate(world, dorivers, dolakes, dodeltas, smalllake, largelake, landshape, mountaindrainage, shelves, importedclimate);
    generatephysicalworldlayers(world, shelves);

    if (appendclimateworkbook && appendclimatebenchmarkworkbook(world) == false)
        updatereport("Climate workbook benchmark update failed");
}

int runcommandlineworldgeneration(const CommandLineGenerationOptions& options)
{
    auto world = std::make_unique<planet>();
    initialiseworld(*world);
    initialisemapcolours(*world);
    world->clear();

    const bool importedworld = options.earthClimateBenchmark || options.hasImportLandPath || options.hasImportSeaPath;
    const long seed = options.hasSeed ? options.seed : (importedworld ? 1 : defaultcommandlineseed());
    world->setseed(seed);

    fast_srand(world->seed());

    if (options.hasReferencePath)
        setreferenceprecipitationgridpath(options.referencePath);

    if (options.earthClimateBenchmark == false && (options.hasImportLandPath != options.hasImportSeaPath))
    {
        cerr << "Imported runs require both --import-land and --import-sea unless --earth-climate-benchmark is used.\n";
        return 1;
    }

    boolshapetemplate landshape[12];
    createlandshapetemplates(landshape);
    boolshapetemplate chainland[2];
    createchainlandtemplates(chainland);
    boolshapetemplate smalllake[12];
    createsmalllaketemplates(smalllake);
    boolshapetemplate largelake[10];
    createlargelaketemplates(largelake);

    vector<int> squareroot((MAXCRATERRADIUS * MAXCRATERRADIUS + MAXCRATERRADIUS + 1) * 24);

    for (int n = 1; n < squareroot.size(); n++)
        squareroot[n] = (int)sqrt(n);

    WorldGenerationDebugOptions debugoptions;
    debugoptions.visualizeEachStep = false;
    debugoptions.logToProfilingWorkbook = options.logToProfilingWorkbook;
    debugoptions.useFastLEMMountains = options.useFastLEMMountains;
    debugoptions.usePlateTectonicsSimulation = options.usePlateTectonicsSimulation;
    debugoptions.plateTectonicsCycleCount = options.plateTectonicsCycleCount;

    if (importedworld == false)
        changeworldproperties(*world);

    updatereport((importedworld ? "Generating imported world from seed: " : "Generating world from seed: ") + to_string(world->seed()) + ":");
    updatereport("");

    beginworldgendebugrun(world->seed(), &debugoptions);
    begintimedreporting();

    if (importedworld)
    {
        const AppEnvironmentConfig& appenv = getappenvironment();
        const string landpath = options.hasImportLandPath ? options.importLandPath : appenv.earthBenchmarkLandPath.string();
        const string seapath = options.hasImportSeaPath ? options.importSeaPath : appenv.earthBenchmarkSeaPath.string();
        vector<vector<bool>> okmountains(ARRAYWIDTH, vector<bool>(ARRAYHEIGHT, false));

        if (importlandheightmap(*world, landpath) == false)
        {
            endtimedreporting();
            endworldgendebugrun();
            cerr << "Failed to import land map: " << landpath << '\n';
            return 1;
        }

        if (importseadepthmap(*world, seapath) == false)
        {
            endtimedreporting();
            endworldgendebugrun();
            cerr << "Failed to import sea map: " << seapath << '\n';
            return 1;
        }

        completeimportedworldgeneration(*world, options.rivers, options.lakes, options.deltas, options.appendClimateWorkbook, smalllake, largelake, landshape, okmountains);
    }
    else
    {
        int iterations = 4;

        for (int n = 0; n < 3; n++)
        {
            if (random(1, 2) == 1)
                iterations++;
            else
                iterations--;
        }

        if (iterations < 1)
            iterations = 1;

        if (iterations > 7)
            iterations = 7;

        int mergefactor = random(1, 15);

        if (random(1, 12) == 1)
            mergefactor = random(1, 25);

        vector<vector<int>> mountaindrainage(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));
        vector<vector<bool>> shelves(ARRAYWIDTH, vector<bool>(ARRAYHEIGHT, 0));

        generateglobalterrain(*world, 0, iterations, mergefactor, -1, -1, landshape, chainland, mountaindrainage, shelves, squareroot);
        generateglobalclimate(*world, options.rivers, options.lakes, options.deltas, smalllake, largelake, landshape, mountaindrainage, shelves);
        generatephysicalworldlayers(*world, shelves);
    }

    endtimedreporting();
    endworldgendebugrun();

    if (options.printClimateRelativeError)
        printclimaterelativeerrorreport(*world);

    updatereport("");
    updatereport("World generation completed.");
    updatereport("Validation output: " + commandlinevalidationdirectory(world->seed()).string());

    if (options.hasSavePath)
    {
        world->saveworld(options.savePath);
        updatereport("Saved world: " + options.savePath);
    }

    return 0;
}
}

// Used to seed the generator.
void fast_srand(long seed)
{
    g_seed = seed;
}

// Compute a pseudorandom integer.
// Output value in range [0, 32767]
int fast_rand(void)
{
    g_seed = (214013 * g_seed + 2531011);
    return (g_seed >> 16) & 0x7FFF;
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
    display.textureSize.y = world.height() + 2;

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

void ensureworkingdirectory()
{
    if (GetFileAttributesA("version.txt") != INVALID_FILE_ATTRIBUTES)
        return;

    char exepath[MAX_PATH] = { 0 };
    GetModuleFileNameA(NULL, exepath, MAX_PATH);

    string executablepath = exepath;
    size_t lastslash = executablepath.find_last_of("\\/");

    if (lastslash == string::npos)
        return;

    executablepath = executablepath.substr(0, lastslash);
    SetCurrentDirectoryA(executablepath.c_str());

    if (GetFileAttributesA("version.txt") == INVALID_FILE_ATTRIBUTES && GetFileAttributesA("..\\..\\version.txt") != INVALID_FILE_ATTRIBUTES)
        SetCurrentDirectoryA("..\\..");
    else if (GetFileAttributesA("version.txt") == INVALID_FILE_ATTRIBUTES && GetFileAttributesA("..\\version.txt") != INVALID_FILE_ATTRIBUTES)
        SetCurrentDirectoryA("..");
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

    CommandLineGenerationOptions commandlineoptions;

    if (parsecommandlineoptions(commandlineoptions) == false)
    {
        printcommandlineusage();
        return 1;
    }

    if (commandlineoptions.showHelp)
    {
        printcommandlineusage();

        if (commandlineoptions.run == false)
            return 0;
    }

    if (commandlineoptions.run)
        return runcommandlineworldgeneration(commandlineoptions);

    float currentversion = 1.0f;
    float latestversion = getlatestversion();

    string currentversionstring = "Current version: " + to_string(currentversion);
    string latestversionstring = "Latest version: " + to_string(latestversion);

    updatereport(currentversionstring.c_str());
    updatereport(latestversionstring.c_str());

    // Set up the window.

    int scwidth = 1224;
    int scheight = 768;

    sf::RenderWindow window(sf::VideoMode(scwidth, scheight), "Undiscovered Worlds");
    window.setFramerateLimit(60);
    ImGui::SFML::Init(window);

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
        else if (strcmp(filter, ".png") == 0)
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
    int& landmass = customworldui.landmass; // Rough amount of land coverage, for custom worlds.
    int& mergefactor = customworldui.mergefactor; // Amount continents will be removed by merging with the fractal map, for custom worlds.-----------------
    int& iterations = customworldui.iterations; // Number of terrain iterations for worlds of terrain type 4.
    int& sealeveleditable = customworldui.sealeveleditable; // Sea level (0-10) for worlds of terrain type 4.
    bool& compareclimateworkbook = customworldui.compareClimateWorkbook;

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
    bool& showtectonicchooser = panels.showTectonicChooser; // If this is 1 then we show the panel for creating tectonic-based custom world terrain.
    bool& shownontectonicchooser = panels.showNonTectonicChooser; // If this is 1 then we show the panel for creating non-tectonic-based custom world terrain.
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
    short& generatingtectonicpass = progresspasses.generatingTectonic; // Tracks which pass we're on for this section, to ensure that widgets are correctly displayed.
    short& generatingnontectonicpass = progresspasses.generatingNonTectonic; // Tracks which pass we're on for this section, to ensure that widgets are correctly displayed.

    WorldGenerationDebugState worldgenerationdebug;

    auto invalidateGlobalMaps = [&]()
    {
        resetmapcache(globalmaps);
    };

    auto invalidateRegionalMaps = [&]()
    {
        resetmapcache(regionalmaps);
    };

    auto showGlobalMapView = [&](mapviewenum view, bool updateMiniMap = false)
    {
        mapview = view;
        applyglobalmapview(view, *world, globalmaps, globalmaptexture, globalmap, updateMiniMap ? &minimap : nullptr);
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
        drawallglobalmapimages(*world, globalmaps);
    };

    auto rebuildAllRegionalMaps = [&]()
    {
        drawallregionalmapimages(*world, *region, regionalmaps);
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

        if (screenmode == createworldscreen)
        {
            showglobaltemperaturechart = 0;
            showglobalrainfallchart = 0;
            showregionaltemperaturechart = 0;
            showregionalrainfallchart = 0;
            showworldproperties = 0;
            showworldeditproperties = 0;
            newworld = 0;

            const CreateWorldScreenActions actions = drawcreateworldscreen(main_viewport, window_flags, currentversion, latestversion, brandnew == 1, seedentry, createrandomseednumber);

            if (actions.openLoadDialog)
            {
                openFileDialog(".uww");
                loadingworld = 1;
            }

            if (actions.returnToGlobalMap)
                screenmode = globalmapscreen;

            if (actions.openCustomWorld)
                showsetsize = 1;

            if (actions.startWorldGeneration)
                showworldgenerationoptions = 1;

            if (drawworldgenerationoptionswindow(main_viewport, window_flags, showworldgenerationoptions, worldgenerationdebug.options))
            {
                world->setseed(seedentry);
                screenmode = creatingworldscreen;
            }
        }

        // Creating world screen

        if (screenmode == creatingworldscreen)
        {
            if (showdeferredworkwindow(creatingworldpass, "Please wait...##createworld ", "Generating world...", ImVec2(main_viewport->WorkPos.x + 506, main_viewport->WorkPos.y + 171), ImVec2(173, 68)))
            {
                updatereport("Generating world from seed: " + to_string(world->seed()) + ":");
                updatereport("");
                beginworldgendebugrun(world->seed(), &worldgenerationdebug.options);

                if (worldgenerationdebug.options.visualizeEachStep)
                {
                    setworldgenvisualizationcallback([&]()
                    {
                        invalidateGlobalMaps();
                        showGlobalMapView(relief);
                        window.clear();
                        window.draw(globalmap);
                        window.display();
                    });
                }
                else
                    clearworldgenvisualizationcallback();

                begintimedreporting();

                initialiseworld(*world);
                world->clear();

                changeworldproperties(*world);

                fast_srand(world->seed());

                resizeglobaldisplayforworld(*world, globaldisplay, highlightdisplay);

                invalidateGlobalMaps();

                int contno = random(1, 9);

                int thismergefactor = random(1, 15);

                if (random(1, 12) == 1) // Fairly rarely, have more fragmented continents
                    thismergefactor = random(1, 25);

                iterations = 4; // This is for worlds of terrain type 4.

                for (int n = 0; n < 3; n++)
                {
                    if (random(1, 2) == 1)
                        iterations++;
                    else
                        iterations--;
                }

                if (iterations < 1)
                    iterations = 1;

                if (iterations > 7)
                    iterations = 7;

                vector<vector<int>> mountaindrainage(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));
                vector<vector<bool>> shelves(ARRAYWIDTH, vector<bool>(ARRAYHEIGHT, 0));

                // Actually generate the world

                generateglobalterrain(*world, 0, iterations, thismergefactor, -1, -1, landshape, chainland, mountaindrainage, shelves,squareroot);
                generateglobalclimate(*world, 1, 1, 1, smalllake, largelake, landshape, mountaindrainage, shelves);
                generatephysicalworldlayers(*world, shelves);

                // Now draw a new map

                mapview = relief;

                showGlobalMapView(mapview);

                endtimedreporting();
                endworldgendebugrun();
                updatereport("");
                updatereport("World generation completed.");
                updatereport("");

                focused = 0;

                creatingworldpass = 0;
                screenmode = movingtoglobalmapscreen;
                newworld = 1;
            }
        }

        // Global map screen

        if (screenmode == globalmapscreen)
        {
            showworldeditproperties = 0;
            showregionaltemperaturechart = 0;
            showregionalrainfallchart = 0;
            
            areafromregional = 0;
            brandnew = 0;
            
            // Main controls.

            string title;

            if (world->seed() >= 0)
                title = "Seed: " + to_string(world->seed());
           else
                title = "Custom";

            ImGui::SetNextWindowPos(ImVec2(main_viewport->WorkPos.x + 10, main_viewport->WorkPos.y + 20), ImGuiCond_FirstUseEver);
            ImGui::SetNextWindowSize(ImVec2(161, 549), ImGuiCond_FirstUseEver);

            ImGui::Begin(title.c_str(), NULL, window_flags);

            ImVec2 pos = ImGui::GetWindowPos();

            ImGui::Text("World controls:");

            if (standardbutton("New world"))
            {
                brandnew = 0;
                seedentry = 0;

                screenmode = createworldscreen;
            }

            if (standardbutton("Load world"))
            {
                openFileDialog(".uww");

                loadingworld = 1;
            }

            if (standardbutton("Save world"))
            {
                openFileDialog(".uww");

                savingworld = 1;
            }

            if (standardbutton("Custom world"))
            {
                showsetsize = 1;
            }

            ImGui::Dummy(ImVec2(0.0f, linespace));

            ImGui::SetNextItemWidth(0);

            ImGui::Text("Export options:");

            ImGui::PushItemWidth(100.0f);

            if (standardbutton("World maps"))
            {
                openFileDialog(".png");

                exportingworldmaps = 1;
            }

            if (standardbutton("Area maps"))
            {
                redrawGlobalRelief();
                
                screenmode = exportareascreen;
                areafromregional = 0;
            }

            ImGui::Dummy(ImVec2(0.0f, linespace));

            ImGui::SetNextItemWidth(0);

            ImGui::Text("Display map type:");

            ImGui::PushItemWidth(100.0f);

            drawmapviewbuttons(showGlobalMapView);

            ImGui::Dummy(ImVec2(0.0f, linespace));

            ImGui::SetNextItemWidth(0);

            ImGui::Text("Other controls:");

            ImGui::PushItemWidth(100.0f);

            if (standardbutton("Properties"))
                toggle(showworldproperties);

            if (standardbutton("Appearance"))
                toggle(showcolouroptions);

            if (standardbutton("Zoom"))
            {
                if (focused == 1)
                {
                    redrawGlobalRelief(true);
                    newx = poix;
                    newy = poiy;

                    screenmode = generatingregionscreen;
                }
            }

            ImGui::End();

            drawglobalinfowindow(main_viewport, window_flags, *world, focused == 1, poix, poiy, degree, cube, showglobaltemperaturechart, showglobalrainfallchart, newworld == 1, showabout);

            // Now check to see if the map has been clicked on.

            if (window.hasFocus() && sf::Mouse::isButtonPressed(sf::Mouse::Left) && io.WantCaptureMouse == 0)
            {
                sf::Vector2i mousepos = sf::Mouse::getPosition(window);

                float mult = ((float)world->width()+1.0f) / (float)DISPLAYMAPSIZEX;

                float fpoix = (float)(mousepos.x - globalmapxpos) * mult;
                float fpoiy = (float)(mousepos.y - globalmapypos) * mult;

                poix = (int)fpoix;
                poiy = (int)fpoiy;

                if (poix >= 0 && poix <= world->width() && poiy >= 0 && poiy <= world->height())
                {
                    focused = 1;
                    newworld = 0;

                    highlight.setPosition(sf::Vector2f((float)mousepos.x, (float)mousepos.y));
                }
                else
                {
                    if (focused == 1)
                    {
                        focused = 0;
                        poix = -1;
                        poiy = -1;
                    }
                }
            }
        }

        // Regional map screen

        if (screenmode == regionalmapscreen)
        {
            showglobaltemperaturechart = 0;
            showglobalrainfallchart = 0;
            showworldeditproperties = 0;
            
            areafromregional = 1;
            
            showsetsize = 0;
            newworld = 0;

            // Main controls.

            string title;

            if (world->seed() >= 0)
                title = "Seed: " + to_string(world->seed());
            else
                title = "Custom";

            title = title + "##regional";

            ImGui::SetNextWindowPos(ImVec2(main_viewport->WorkPos.x + 10, main_viewport->WorkPos.y + 20), ImGuiCond_FirstUseEver);
            ImGui::SetNextWindowSize(ImVec2(161, 449), ImGuiCond_FirstUseEver);

            ImGui::Begin(title.c_str(), NULL, window_flags);

            ImGui::Text("World controls:");

            ImGui::PushItemWidth(100.0f);

            if (standardbutton("World map"))
            {
                focused = 0;
                poix = -1;
                poiy = -1;

                regionalreliefimage.create(regionalimagewidth, regionalimageheight, sf::Color::Black);

                redrawRegionalRelief();

                screenmode = globalmapscreen;

                showcolouroptions = 0;
                showworldproperties = 0;
                showregionalrainfallchart = 0;
                showregionaltemperaturechart = 0;
            }

            ImGui::Dummy(ImVec2(0.0f, linespace));

            ImGui::SetNextItemWidth(0);

            ImGui::Text("Export options:");

            ImGui::PushItemWidth(100.0f);

            if (standardbutton("Regional maps"))
            {
                openFileDialog(".png");

                exportingregionalmaps = 1;
            }

            if (standardbutton("Area maps"))
            {
                screenmode = exportareascreen;
                areafromregional = 1;
            }

            ImGui::Dummy(ImVec2(0.0f, linespace));

            ImGui::SetNextItemWidth(0);

            ImGui::Text("Display map type:");

            ImGui::PushItemWidth(100.0f);

            drawmapviewbuttons(showRegionalMapView);

            ImGui::Dummy(ImVec2(0.0f, linespace));

            ImGui::SetNextItemWidth(0);

            ImGui::Text("Other controls:");

            ImGui::PushItemWidth(100.0f);

            if (standardbutton("Properties"))
                toggle(showworldproperties);

            if (standardbutton("Appearance"))
                toggle(showcolouroptions);

            ImGui::End();

            drawregionalinfowindow(main_viewport, window_flags, *world, *region, focused == 1, poix, poiy, degree, cube, showregionaltemperaturechart, showregionalrainfallchart, showabout);

            if (window.hasFocus())
            {
                // Now check to see if the map has been clicked on.
                
                if (sf::Mouse::isButtonPressed(sf::Mouse::Left) && io.WantCaptureMouse == 0)
                {
                    sf::Vector2i mousepos = sf::Mouse::getPosition(window);

                    poix = mousepos.x - (int)regionalmapxpos;
                    poiy = mousepos.y - (int)regionalmapypos;

                    if (poix >= 0 && poix < regionalmapimagewidth && poiy >= 0 && poiy < regionalmapimageheight)
                    {
                        poix = poix + region->regwidthbegin();
                        poiy = poiy + region->regheightbegin();

                        focused = 1;

                        highlight.setPosition(sf::Vector2f((float)mousepos.x, (float)mousepos.y));
                    }
                    else
                    {
                        if (focused == 1)
                        {
                            focused = 0;
                            poix = -1;
                            poiy = -1;
                        }
                    }

                    float mult = ((float)world->width() + 1.0f) / (float)DISPLAYMAPSIZEX;

                    float fpoix = (float)(mousepos.x - minimapxpos) * mult * 2.0f;
                    float fpoiy = (float)(mousepos.y - minimapypos) * mult * 2.0f;

                    int minipoix = (int)fpoix;
                    int minipoiy = (int)fpoiy;

                    if (minipoix >= 0 && minipoiy <= world->width() && minipoiy >= 0 && minipoiy <= world->height())// If the minimap has been clicked on.
                    {
                        screenmode = generatingregionscreen;
                        newx = minipoix;
                        newy = minipoiy;
                    }
                }

                // Check to see if the cursor keys have been pressed.

                if (sf::Keyboard::isKeyPressed(sf::Keyboard::Left))
                {
                    int regionalmapmove = REGIONALTILEWIDTH;

                    newx = region->centrex() - regionalmapmove;

                    if (newx > 0)
                    {
                        newx = newx / 32;
                        newx = newx * 32;
                        newx = newx + 16;
                    }

                    if (newx < regionalmapmove * 2)
                        newx = wrap(newx, world->width());

                    newy = region->centrey();

                    screenmode = generatingregionscreen;
                }

                if (sf::Keyboard::isKeyPressed(sf::Keyboard::Right))
                {
                    int regionalmapmove = REGIONALTILEWIDTH;

                    newx = region->centrex() + regionalmapmove;

                    newx = newx / 32;
                    newx = newx * 32;
                    newx = newx + 16;

                    if (newx > world->width())
                        newx = wrap(newx, world->width());

                    newy = region->centrey();

                    screenmode = generatingregionscreen;
                }

                if (sf::Keyboard::isKeyPressed(sf::Keyboard::Up))
                {
                    int regionalmapmove = REGIONALTILEHEIGHT;

                    newy = region->centrey() - regionalmapmove;

                    newy = newy / 32;
                    newy = newy * 32;
                    newy = newy + 16;

                    if (newy > regionalmapmove)
                    {
                        newx = region->centrex();
                        screenmode = generatingregionscreen;
                    }
                }

                if (sf::Keyboard::isKeyPressed(sf::Keyboard::Down))
                {
                    int regionalmapmove = REGIONALTILEHEIGHT;

                    newy = region->centrey() + regionalmapmove;

                    newy = newy / 32;
                    newy = newy * 32;
                    newy = newy + 16;

                    if (newy < world->height())
                    {
                        newx = region->centrex();
                        screenmode = generatingregionscreen;
                    }
                }
            }
        }

        // Area export screen (this is where the user selects the area to export)

        if (screenmode == exportareascreen)
        {
            showglobaltemperaturechart = 0;
            showglobalrainfallchart = 0;
            showregionaltemperaturechart = 0;
            showregionalrainfallchart = 0;
            showworldproperties = 0;
            showworldeditproperties = 0;
            showcolouroptions = 0;
            newworld = 0;
            
            showsetsize = 0;

            // Work out the size of the currently selected area.

            int totalregions = 0;
            int maxtotalregions = 100; // Query area maps larger than this.
            
            if (areanex != -1)
            {
                float regiontilewidth = (float)REGIONALTILEWIDTH; //30;
                float regiontileheight = (float)REGIONALTILEHEIGHT; //30; // The width and height of the visible regional map, in tiles.

                int regionwidth = (int)regiontilewidth * 16;
                int regionheight = (int)regiontileheight * 16; // The width and height of the visible regional map, in pixels.

                int newareanwx = areanwx; // This is because the regions we'll be making will start to the north and west of the defined area.
                int newareanwy = areanwy;

                int newareanex = areanex;
                int newareaney = areaney;

                int newareaswx = areaswx;
                int newareaswy = areaswy;

                int newareasey = areasey;

                newareanwx = newareanwx / (int)regiontilewidth;
                newareanwx = newareanwx * (int)regiontilewidth;

                newareanwy = newareanwy / (int)regiontileheight;
                newareanwy = newareanwy * (int)regiontileheight;

                newareaswx = newareanwx;
                newareaney = newareanwy;

                float areatilewidth = (float)(newareanex - newareanwx);
                float areatileheight = (float)(newareasey - newareaney);

                int areawidth = (int)(areatilewidth * 16.0f);
                int areaheight = (int)(areatileheight * 16.0f);

                float fregionswide = areatilewidth / regiontilewidth;
                float fregionshigh = areatileheight / regiontileheight;

                int regionswide = (int)fregionswide;
                int regionshigh = (int)fregionshigh;

                if (regionswide != fregionswide)
                    regionswide++;

                if (regionshigh != fregionshigh)
                    regionshigh++;

                totalregions = regionswide * regionshigh; // This is how many regional maps we're going to have to do.
            }

            // Main controls.

            ImGui::SetNextWindowPos(ImVec2(main_viewport->WorkPos.x + 10, main_viewport->WorkPos.y + 20), ImGuiCond_FirstUseEver);
            ImGui::SetNextWindowSize(ImVec2(161, 138), ImGuiCond_FirstUseEver);

            ImGui::Begin("Export custom area", NULL, window_flags);

            ImGui::PushItemWidth(100.0f);

            if (standardbutton("Clear selection"))
            {
                areanex = -1;
                areaney = -1;
                areasex = -1;
                areasey = -1;
                areaswx = -1;
                areaswy = -1;
                areanwx = -1;
                areanwy = -1;
            }

            if (standardbutton("Export maps"))
            {
                if (areanex != -1)
                {
                    if (totalregions <= maxtotalregions)
                    {
                        openFileDialog(".png");

                        exportingareamaps = 1;
                    }
                    else
                        showareawarning = 1;
                }
            }

            ImGui::Text("  ");

            if (standardbutton("Cancel"))
            {
                areanex = -1;
                areaney = -1;
                areasex = -1;
                areasey = -1;
                areaswx = -1;
                areaswy = -1;
                areanwx = -1;
                areanwy = -1;

                if (areafromregional == 1)
                    screenmode = regionalmapscreen;
                else
                    screenmode = globalmapscreen;
            }

            ImGui::End();

            // Now the text box.

            ImGui::SetNextWindowPos(ImVec2(main_viewport->WorkPos.x + 180, main_viewport->WorkPos.y + 542), ImGuiCond_FirstUseEver);
            ImGui::SetNextWindowSize(ImVec2(1023, 151), ImGuiCond_FirstUseEver);
            string title = "            ";
            string areatext = "Total regions: " + to_string(totalregions); // 117 OK; 143 not OK; 132 not OK; 121 not OK; 112 OK; 116 OK; 110

            ImGui::Begin(title.c_str(), NULL, window_flags);
            ImGui::PushItemWidth((float)world->width() / 2.0f);
            areatext="This screen allows you to export maps at the same scale as the regional map, but of larger areas.\n\nClick on the map to pick the corners of the area you want to export. You can re-select corners to fine-tune the area.\n\nWhen you are done, click on 'export maps'. The program will ask you to specify the filename under which to save the maps, and then create them.";
            
            ImGui::Text(areatext.c_str(), world->width() / 2);
            ImGui::End();

            // Now check to see if the map has been clicked on.

            if (window.hasFocus() && sf::Mouse::isButtonPressed(sf::Mouse::Left) && io.WantCaptureMouse == 0)
            {
                sf::Vector2i mousepos = sf::Mouse::getPosition(window);

                float mult = ((float)world->width() + 1.0f) / (float)DISPLAYMAPSIZEX;

                float xpos = ((float)mousepos.x - (float)globalmapxpos) * mult;
                float ypos = ((float)mousepos.y - (float)globalmapypos) * mult;

                poix = (int)xpos;
                poiy = (int)ypos;

                if (poix >= 0 && poix <= world->width() && poiy >= 0 && poiy <= world->height())
                {
                    if (areaswx == -1) // If we don't have any corners yet
                    {
                        areanex = poix + 2;
                        areaney = poiy;
                        areasex = poix + 2;
                        areasey = poiy + 2;
                        areaswx = poix;
                        areaswy = poiy + 2;
                        areanwx = poix;
                        areanwy = poiy;
                    }
                    else // If we do have all corners
                    {
                        int nedistx = areanex - poix;
                        int nedisty = areaney - poiy;

                        int sedistx = areasex - poix;
                        int sedisty = areasey - poiy;

                        int swdistx = areaswx - poix;
                        int swdisty = areaswy - poiy;

                        int nwdistx = areanwx - poix;
                        int nwdisty = areanwy - poiy;

                        float nedist = (float)sqrt(nedistx * nedistx + nedisty * nedisty);
                        float sedist = (float)sqrt(sedistx * sedistx + sedisty * sedisty);
                        float swdist = (float)sqrt(swdistx * swdistx + swdisty * swdisty);
                        float nwdist = (float)sqrt(nwdistx * nwdistx + nwdisty * nwdisty);

                        short id = 2; // This will identify which corner this is going to be.
                        float mindist = nedist;

                        if (sedist < mindist)
                        {
                            id = 4;
                            mindist = sedist;
                        }

                        if (swdist < mindist)
                        {
                            id = 6;
                            mindist = swdist;
                        }

                        if (nwdist < mindist)
                            id = 8;

                        if (id == 2)
                        {
                            areanex = poix;
                            areaney = poiy;
                            areasex = poix;
                            areanwy = poiy;
                        }

                        if (id == 4)
                        {
                            areasex = poix;
                            areasey = poiy;
                            areanex = poix;
                            areaswy = poiy;
                        }

                        if (id == 6)
                        {
                            areaswx = poix;
                            areaswy = poiy;
                            areanwx = poix;
                            areasey = poiy;
                        }

                        if (id == 8)
                        {
                            areanwx = poix;
                            areanwy = poiy;
                            areaswx = poix;
                            areaney = poiy;
                        }

                        if (areanwx > areanex)
                            swap(areanwx, areanex);

                        if (areaswx > areasex)
                            swap(areaswx, areasex);

                        if (areanwy > areaswy)
                            swap(areanwy, areaswy);

                        if (areaney > areasey)
                            swap(areaney, areasey);
                    }
                }
            }
        }

        // Custom world screen

        if (screenmode == importscreen)
        {
            showglobaltemperaturechart = 0;
            showglobalrainfallchart = 0;
            showregionaltemperaturechart = 0;
            showregionalrainfallchart = 0;
            showworldproperties = 0;
            
            showcolouroptions = 0;
            showworldproperties = 0;
            
            // Main controls.

            ImGui::SetNextWindowPos(ImVec2(main_viewport->WorkPos.x + 10, main_viewport->WorkPos.y + 20), ImGuiCond_FirstUseEver);
            ImGui::SetNextWindowSize(ImVec2(161, 735), ImGuiCond_FirstUseEver);

            ImGui::Begin("Custom##custom", NULL, window_flags);

            ImGui::SetNextItemWidth(0);
            ImGui::Text("Value mode:");

            int importvaluemode = static_cast<int>(mapimportsettings.valueMode);
            const char* importmodeitems[] = { "Red", "Strip" };

            ImGui::PushItemWidth(100.0f);

            if (ImGui::Combo("##importvaluemode", &importvaluemode, importmodeitems, IM_ARRAYSIZE(importmodeitems)))
                mapimportsettings.valueMode = static_cast<MapImportValueMode>(importvaluemode);

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Choose whether imports read the red channel directly or decode colours through a 1-pixel-tall gradient strip.");

            if (mapimportsettings.valueMode == MapImportValueMode::gradientStrip)
            {
                if (standardbutton("Scale strip"))
                {
                    openFileDialog(".png");
                    importinggradientstrip = 1;
                }

                if (ImGui::IsItemHovered())
                    ImGui::SetTooltip("Load a 1-pixel-tall PNG strip. The leftmost colour maps to the minimum value, and each pixel to the right adds the chosen increment.");

                ImGui::InputFloat("Minimum", &mapimportsettings.gradientMinimum, 1.0f, 10.0f, "%.3f");

                if (ImGui::IsItemHovered())
                    ImGui::SetTooltip("Value assigned to the leftmost colour in the gradient strip.");

                ImGui::InputFloat("Step", &mapimportsettings.gradientIncrement, 0.1f, 1.0f, "%.3f");

                if (ImGui::IsItemHovered())
                    ImGui::SetTooltip("Amount added for each colour step to the right in the gradient strip.");

                const char* stripstatus = mapimportsettings.gradientStripPath.empty() ? "Strip: none" : "Strip: loaded";
                ImGui::TextUnformatted(stripstatus);

                if (!mapimportsettings.gradientStripPath.empty() && ImGui::IsItemHovered())
                    ImGui::SetTooltip(mapimportsettings.gradientStripPath.c_str());
            }

            ImGui::Dummy(ImVec2(0.0f, linespace));

            ImGui::Text("Import:");

            if (standardbutton("Land map"))
            {
                openFileDialog(".png");

                importinglandmap = 1;
            }

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Red mode: red 0 is sea and higher values are elevation above sea level in 10-metre steps. Strip mode: decoded values are treated as elevation above sea level.");

            if (standardbutton("Sea map"))
            {
                openFileDialog(".png");

                importingseamap = 1;
            }

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Red mode: red 0 is land and higher values are depth below sea level in 50-metre steps. Strip mode: decoded values are treated as depth below sea level.");

            if (standardbutton("Mountains"))
            {
                openFileDialog(".png");

                importingmountainsmap = 1;
            }

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Red mode: red values scale the imported mountain height. Strip mode: decoded values are used directly as imported mountain height.");

            if (standardbutton("Volcanoes"))
            {
                openFileDialog(".png");

                importingvolcanoesmap = 1;
            }

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Red or strip decoding sets volcano height. Green=0 for shield volcano, or higher for stratovolcano. Blue=0 for extinct, or higher for active.");

            if (standardbutton("Temperature"))
            {
                openFileDialog(".png");

                importingtemperaturemap = 1;
            }

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Red mode maps 0..255 linearly to -60 C..+60 C mean annual temperature. Strip mode uses the decoded strip value directly.");

            if (standardbutton("Precipitation"))
            {
                openFileDialog(".png");

                importingprecipitationmap = 1;
            }

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Red mode maps 0..255 linearly to 0..1020 mean precipitation. Strip mode uses the decoded strip value directly.");

            ImGui::Dummy(ImVec2(0.0f, linespace));

            ImGui::SetNextItemWidth(0);

            ImGui::Text("Generate terrain:");

            ImGui::PushItemWidth(100.0f);

            if (standardbutton("Tectonic"))
                toggle(showtectonicchooser);

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Create tectonic-based global terrain, with continents, mountains, coastal shelves, and oceanic ridges.");

            if (standardbutton("Non-tectonic"))
                toggle(shownontectonicchooser);

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Create non-tectonic-based global terrain.");

            ImGui::Dummy(ImVec2(0.0f, linespace));

            ImGui::SetNextItemWidth(0);

            ImGui::Text("Generate elements:");

            ImGui::PushItemWidth(100.0f);

            if (standardbutton("Coastlines"))
            {
                removestraights(*world);

                // Now redraw the map.

                resetmapcache(globalmaps);

                mapview = relief;
                showGlobalMapView(mapview);
            }

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Disrupt straight edges on coastlines.");

            if (standardbutton("Shelves"))
            {
                int width = world->width();
                int height = world->height();
                int maxelev = world->maxelevation();
                int sealevel = world->sealevel();

                vector<vector<bool>> shelves(ARRAYWIDTH, vector<bool>(ARRAYHEIGHT, 0));

                makecontinentalshelves(*world, shelves, 4);
                createoceantrenches(*world, shelves);

                int grain = 8; // Level of detail on this fractal map.
                float valuemod = 0.2f;
                int v = random(3, 6);
                float valuemod2 = (float)v;

                vector<vector<int>> seafractal(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

                createfractal(seafractal, width, height, grain, valuemod, valuemod2, 1, maxelev, 0, 0);

                float coastalvarreduce = (float)maxelev / 3000.0f;
                float oceanvarreduce = (float)maxelev / 100.0f;

                for (int i = 0; i <= width; i++)
                {
                    for (int j = 0; j <= height; j++)
                    {
                        if (world->sea(i, j) == 1 && shelves[i][j] == 1)
                            world->setnom(i, j, sealevel - 100);
                    }
                }

                // Now redraw the map.

                resetmapcache(globalmaps);

                mapview = relief;
                showGlobalMapView(mapview);
            }

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Generate continental shelves.");

            if (standardbutton("Oceanic ridges"))
            {
                int width = world->width();
                int height = world->height();
                int maxelev = world->maxelevation();
                int sealevel = world->sealevel();

                vector<vector<bool>> shelves(ARRAYWIDTH, vector<bool>(ARRAYHEIGHT, 0));

                for (int i = 0; i <= width; i++)
                {
                    for (int j = 0; j <= height; j++)
                    {
                        if (world->nom(i, j) > sealevel - 400)
                            shelves[i][j] = 1;
                    }
                }

                createoceanridges(*world, shelves);

                // Now redraw the map.

                resetmapcache(globalmaps);

                mapview = relief;
                showGlobalMapView(mapview);
            }

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Generate mid-ocean ridges.");

            if (standardbutton("Sea bed"))
            {
                int width = world->width();
                int height = world->height();
                int maxelev = world->maxelevation();
                int sealevel = world->sealevel();

                int grain = 8; // Level of detail on this fractal map.
                float valuemod = 0.2f;
                int v = random(3, 6);
                float valuemod2 = (float)v;

                vector<vector<int>> seafractal(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

                createfractal(seafractal, width, height, grain, valuemod, valuemod2, 1, maxelev, 0, 0);

                int warpfactor = random(20, 80);
                warp(seafractal, width, height, maxelev, warpfactor, 1);

                float coastalvarreduce = (float)maxelev / 3000.0f;
                float oceanvarreduce = (float)maxelev / 1000.0f;

                for (int i = 0; i <= width; i++)
                {
                    for (int j = 0; j <= height; j++)
                    {
                        if (world->sea(i, j) == 1)
                        {
                            float var = (float)(seafractal[i][j] - maxelev / 2);
                            var = var / coastalvarreduce;
                            
                            int newval = world->nom(i, j) + (int)var;

                            if (newval > sealevel - 10)
                                newval = sealevel - 10;

                            if (newval < 1)
                                newval = 1;

                            world->setnom(i, j, newval);
                        }
                    }
                }

                // Smooth the seabed.

                smoothonlysea(*world, 2);

                // Now redraw the map.

                resetmapcache(globalmaps);

                mapview = relief;
                showGlobalMapView(mapview);
            }

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Create random depth variation across the oceans.");

            if (standardbutton("Land elevation"))
            {
                int width = world->width();
                int height = world->height();
                int maxelev = world->maxelevation();
                int sealevel = world->sealevel();

                // First, make a fractal map.

                int grain = 8; // Level of detail on this fractal map.
                float valuemod = 0.2f;
                int v = random(3, 6);
                float valuemod2 = (float)v;

                vector<vector<int>> fractal(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

                createfractal(fractal, width, height, grain, valuemod, valuemod2, 1, 12750, 0, 0);

                int warpfactor = random(20, 80);
                warp(fractal, width, height, maxelev, warpfactor, 1);

                int fractaladd = sealevel - 2500;

                for (int i = 0; i <= width; i++)
                {
                    for (int j = 0; j <= height; j++)
                        fractal[i][j] = fractal[i][j] + fractaladd;
                }

                // Now use it to change the land heights.

                fractaladdland(*world, fractal);

                // Smooth the land.

                //smoothonlyland(*world, 2);

                // Also, create extra elevation.

                createextraelev(*world);

                // Now redraw the map.

                resetmapcache(globalmaps);

                mapview = relief;
                showGlobalMapView(mapview);
            }

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Generate random elevation variation across the land.");

            if (standardbutton("Mountains##2"))
            {
                int width = world->width();
                int height = world->height();
                int maxelev = world->maxelevation();
                int sealevel = world->sealevel();

                int baseheight = sealevel - 4500;
                if (baseheight < 1)
                    baseheight = 1;
                int conheight = sealevel + 50;

                vector<vector<int>> plateaumap(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

                // First, make a fractal map.

                int grain = 8; // Level of detail on this fractal map.
                float valuemod = 0.2f;
                int v = random(3, 6);
                float valuemod2 = (float)v;

                vector<vector<int>> fractal(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

                createfractal(fractal, width, height, grain, valuemod, valuemod2, 1, 12750, 0, 0);

                int fractaladd = sealevel - 2500;

                for (int i = 0; i <= width; i++)
                {
                    for (int j = 0; j <= height; j++)
                        fractal[i][j] = fractal[i][j] + fractaladd;
                }

                twointegers dummy[1];

                createchains(*world, baseheight, conheight, fractal, plateaumap, landshape, chainland, dummy, 0, 0, 5);

                // Now redraw the map.

                resetmapcache(globalmaps);

                mapview = relief;
                showGlobalMapView(mapview);
            }

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Generate random mountain ranges.");

            if (standardbutton("Hills"))
            {
                int width = world->width();
                int height = world->height();
                int maxelev = world->maxelevation();
                int sealevel = world->sealevel();

                int baseheight = sealevel - 4500;
                if (baseheight < 1)
                    baseheight = 1;
                int conheight = sealevel + 50;

                vector<vector<int>> plateaumap(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

                // First, make a fractal map.

                int grain = 8; // Level of detail on this fractal map.
                float valuemod = 0.2f;
                int v = random(3, 6);
                float valuemod2 = (float)v;

                vector<vector<int>> fractal(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

                createfractal(fractal, width, height, grain, valuemod, valuemod2, 1, 12750, 0, 0);

                int fractaladd = sealevel - 2500;

                for (int i = 0; i <= width; i++)
                {
                    for (int j = 0; j <= height; j++)
                        fractal[i][j] = fractal[i][j] + fractaladd;
                }

                twointegers dummy[1];

                createchains(*world, baseheight, conheight, fractal, plateaumap, landshape, chainland, dummy, 0, 0, 5);

                // Now redraw the map.

                resetmapcache(globalmaps);

                mapview = relief;
                showGlobalMapView(mapview);
            }

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Generate random ranges of hills.");

            if (standardbutton("Craters"))
            {
                int width = world->width();
                int height = world->height();
                int maxelev = world->maxelevation();
                int sealevel = world->sealevel();

                vector<vector<int>> oldterrain(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0)); // Copy the terrain as it is before adding craters. This is so that we can add some variation to the craters afterwards, if there's sea on this world, so the depression filling doesn't fill them completely.

                for (int i = 0; i <= width; i++)
                {
                    for (int j = 0; j <= height; j++)
                        oldterrain[i][j] = world->nom(i, j);
                }

                int cratertotal = random(500, 10000);

                createcratermap(*world, cratertotal, squareroot, 1);

                float minseaproportion = (float)random(1, 500);
                minseaproportion = minseaproportion / 1000.0f;

                int minseasize = (int)(((float)width * (float)height) * minseaproportion);

                if (random(1, 4) != 1) // That may have produced seas inside craters, so probably remove those now.
                    removesmallseas(*world, minseasize, sealevel + 1);
                else // Definitely remove any really little bits of sea.
                    removesmallseas(*world, 20, sealevel + 1);

                int totalsea = 0;

                for (int i = 0; i <= -width; i++)
                {
                    for (int j = 0; j <= height; j++)
                    {
                        if (world->sea(i, j))
                            totalsea++;
                    }
                }

                if (totalsea > 40) // If there's sea, add some variation to how prominent the craters are. This will create some gaps in the craters so they don't get entirely filled up by the depression filling.
                {
                    vector<vector<int>> fractal(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

                    int grain = 8; // Level of detail on this fractal map.
                    float valuemod = 0.2f;
                    int v = random(3, 6);
                    float valuemod2 = (float)v;

                    createfractal(fractal, width, height, grain, valuemod, valuemod2, 1, maxelev, 1, 0);

                    for (int i = 0; i <= width; i++)
                    {
                        for (int j = 0; j <= height; j++)
                        {
                            float oldmult = (float)fractal[i][j] / (float)maxelev;
                            float newmult = 1.0f - oldmult;

                            float thiselev = (float)oldterrain[i][j] * oldmult + (float)world->nom(i, j) * newmult;

                            world->setnom(i, j, (int)thiselev);
                        }
                    }
                }

                // Now redraw the map.

                resetmapcache(globalmaps);

                mapview = relief;
                showGlobalMapView(mapview);
            }

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Generate random craters.");

            ImGui::Dummy(ImVec2(0.0f, linespace));

            ImGui::SetNextItemWidth(0);

            ImGui::Text("Other controls:");

            ImGui::PushItemWidth(100.0f);

            if (standardbutton("Properties"))
                toggle(showworldeditproperties);

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Adjust the world's properties.");

            if (standardbutton("Finish"))
            {
                screenmode = completingimportscreen;
            }

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Calculate climates, lakes, and rivers, and finish the world.");

            ImGui::Checkbox("Append climate.xlsx benchmark row", &compareclimateworkbook);

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Appends this world's land-climate counts to the RAW_PIXELS sheet in climate.xlsx.");

            if (standardbutton("Cancel"))
            {
                brandnew = 1;
                seedentry = 0;

                screenmode = createworldscreen;
            }

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Return to the world creation screen.");

            ImGui::End();

            // Now the text box.

            ImGui::SetNextWindowPos(ImVec2(main_viewport->WorkPos.x + 180, main_viewport->WorkPos.y + 542), ImGuiCond_FirstUseEver);
            ImGui::SetNextWindowSize(ImVec2(1023, 155), ImGuiCond_FirstUseEver);

            string title = "            ";

            string importtext = "You can use the 'import' buttons to load in your own maps. These must be " + formatnumber(world->width() + 1) + " x " + formatnumber(world->height() + 1) + " pixels, in .png format.\nTerrain maps define land, sea, mountains, and volcanoes. Temperature and precipitation imports set the finished climate targets. Imports use the current value mode shown on the left.\nAlternatively, you can use the 'World terrain' button to generate a map from scratch.\nAfter you have imported or generated the map, you can use the other 'generate' buttons to tweak it or to add extra features.\nYou can use the 'Properties' panel to change settings such as global temperatures or rainfall.\nWhen you are done, click 'Finish' to finish the world.";

            ImGui::Begin(title.c_str(), NULL, window_flags);
            ImGui::PushItemWidth((float)(world->width() / 2));
            ImGui::Text(importtext.c_str(), world->width() / 2);
            ImGui::End();
        }

        // These screens all display a "Please wait" message ten times (for some reason doing it once or twice doesn't actually display it) and then do something time-consuming.

        if (screenmode == generatingtectonicscreen)
        {
            if (showdeferredworkwindow(generatingtectonicpass, "Please wait...##generatetectonic ", "Generating terrain...", ImVec2(main_viewport->WorkPos.x + 507, main_viewport->WorkPos.y + 173), ImVec2(173, 68)))
            {
                preparecustomworldgeneration();
                world->settype(2); // This terrain type gives large continents.

                int clusterno = -1;
                int clustersize = -1;

                switch (landmass)
                {
                case 0:
                    world->settype(3);
                    break;

                case 1:
                    world->settype(1);
                    break;

                case 2:
                    clusterno = 1;
                    clustersize = 1;
                    break;

                case 3:
                    clusterno = 1;
                    clustersize = 6;
                    break;

                case 4:
                    clusterno = 2;
                    clustersize = 3;
                    break;

                case 5:
                    clusterno = 2;
                    clustersize = 7;
                    break;

                case 6:
                    clusterno = 3;
                    clustersize = 3;
                    break;

                case 7:
                    clusterno = 3;
                    clustersize = 8;
                    break;

                case 8:
                    clusterno = 3;
                    clustersize = 9;
                    break;

                case 9:
                    clusterno = 4;
                    clustersize = 5;
                    break;

                case 10:
                    clusterno = 4;
                    clustersize = 9;
                    break;
                }

                vector<vector<int>> mountaindrainage(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));
                vector<vector<bool>> shelves(ARRAYWIDTH, vector<bool>(ARRAYHEIGHT, 0));

                // Now generate the terrain.

                generateglobalterrain(*world, 1, iterations, mergefactor - 5, clusterno, clustersize, landshape, chainland, mountaindrainage, shelves, squareroot);
                finishcustomworldgeneration(generatingtectonicpass);
            }
        }

        if (screenmode == generatingnontectonicscreen)
        {
            if (showdeferredworkwindow(generatingnontectonicpass, "Please wait...##generatenontectonic ", "Generating terrain...", ImVec2(main_viewport->WorkPos.x + 507, main_viewport->WorkPos.y + 173), ImVec2(173, 68)))
            {
                preparecustomworldgeneration();
                world->settype(4); // This terrain type gives alien-type terrain.

                float sealevel = (float)sealeveleditable / 10.0f;

                sealevel = (float)world->maxelevation() * sealevel;

                if (sealevel < 1.0f)
                    sealevel = 1.0f;

                world->setsealevel((int)sealevel);

                vector<vector<int>> mountaindrainage(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));
                vector<vector<bool>> shelves(ARRAYWIDTH, vector<bool>(ARRAYHEIGHT, 0));

                // Now generate the terrain.

                generateglobalterrain(*world, 1, iterations, mergefactor - 5, 1, 1, landshape, chainland, mountaindrainage, shelves, squareroot);
                finishcustomworldgeneration(generatingnontectonicpass);
            }
        }

        if (screenmode == completingimportscreen)
        {
            if (showdeferredworkwindow(completingimportpass, "Please wait...##completeimport ", "Finishing world...", ImVec2(main_viewport->WorkPos.x + 507, main_viewport->WorkPos.y + 173), ImVec2(173, 68)))
            {
                updatereport("Generating custom world:");
                updatereport("");
                begintimedreporting();

                // Plug in the world settings.

                applyworldpropertycontrols(*world, worldpropertycontrols);
                completeimportedworldgeneration(*world, currentrivers, currentlakes, currentdeltas, compareclimateworkbook, smalllake, largelake, landshape, OKmountains, &importedclimatemaps);

                // Now draw a new map

                invalidateGlobalMaps();

                mapview = relief;

                showGlobalMapView(mapview);

                endtimedreporting();
                updatereport("");
                updatereport("World generation completed.");
                updatereport("");

                focused = 0;

                completingimportpass = 0;
                screenmode = movingtoglobalmapscreen;
                newworld = 1;
            }
        }

        if (screenmode == loadingworldscreen)
        {
            if (showdeferredworkwindow(loadingworldpass, "Please wait...##loadworld ", "Loading world...", ImVec2(main_viewport->WorkPos.x + 507, main_viewport->WorkPos.y + 173), ImVec2(173, 68)))
            {
                bool success=world->loadworld(filepathname);

                if (success == 0) // Failed to load
                {
                    screenmode = loadfailure;
                    loadingworld = 0;
                    loadingworldpass = 0;
                }
                else
                {
                    syncworldpropertycontrols(*world, worldpropertycontrols);
                    syncappearancesettings(*world, appearance);
                    
                    // Now draw the new world.
                    
                    resizeglobaldisplayforworld(*world, globaldisplay, highlightdisplay);

                    resetmapcache(globalmaps);

                    mapview = relief;
                    showGlobalMapView(mapview);

                    focused = 0;
                    newworld = 1;

                    filepathname = "";
                    filepath = "";

                    screenmode = globalmapscreen;
                    loadingworld = 0;
                    loadingworldpass = 0;
                }
            }
        }

        if (screenmode == savingworldscreen)
        {
            if (showdeferredworkwindow(savingworldpass, "Please wait...##saveworld ", "Saving world...", ImVec2(main_viewport->WorkPos.x + 507, main_viewport->WorkPos.y + 173), ImVec2(173, 68)))
            {
                world->saveworld(filepathname);

                filepathname = "";
                filepath = "";

                savingworld = 0;

                savingworldpass = 0;

                screenmode = movingtoglobalmapscreen;
            }
        }

        if (screenmode == exportingareascreen)
        {
            if (showdeferredworkwindow(exportingareapass, "Please wait...##exportarea ", "Generating area maps...", ImVec2(main_viewport->WorkPos.x + 300, main_viewport->WorkPos.y + 200), ImVec2(173, 68)))
            {
                int oldregionalcentrex = region->centrex();
                int oldregionalcentrey = region->centrey();

                mapviewenum oldmapview = mapview;

                float regiontilewidth = (float)REGIONALTILEWIDTH; //30;
                float regiontileheight = (float)REGIONALTILEHEIGHT; //30; // The width and height of the visible regional map, in tiles.

                int regionwidth = (int)(regiontilewidth * 16.0f);
                int regionheight = (int)(regiontileheight * 16.0f); // The width and height of the visible regional map, in pixels.

                int origareanwx = areanwx; // This is because the regions we'll be making will start to the north and west of the defined area.
                int origareanwy = areanwy;

                int origareanex = areanex;
                int origareaney = areaney;

                int origareaswx = areaswx;
                int origareaswy = areaswy;

                areanwx = areanwx / (int)regiontilewidth;
                areanwx = areanwx * (int)regiontilewidth;

                areanwy = areanwy / (int)regiontileheight;
                areanwy = areanwy * (int)regiontileheight;

                areaswx = areanwx;
                areaney = areanwy;

                int woffset = (origareanwx - areanwx) * 16;
                int noffset = (origareanwy - areanwy) * 16;

                float areatilewidth = (float)(areanex - areanwx);
                float areatileheight = (float)(areasey - areaney);

                int areawidth = (int)areatilewidth * 16;
                int areaheight = (int)areatileheight * 16;

                float imageareatilewidth = (float)(origareanex - origareanwx);
                float imageareatileheight = (float)(areasey - origareaney);

                int imageareawidth = (int)imageareatilewidth * 16;
                int imageareaheight = (int)imageareatileheight * 16;

                float fregionswide = areatilewidth / regiontilewidth;
                float fregionshigh = areatileheight / regiontileheight;

                int regionswide = (int)fregionswide;
                int regionshigh = (int)fregionshigh;

                if (regionswide != fregionswide)
                    regionswide++;

                if (regionshigh != fregionshigh)
                    regionshigh++;

                int totalregions = regionswide * regionshigh; // This is how many regional maps we're going to have to do.

                if (areafromregional == 1)
                    totalregions++; // Because we'll have to redo the regional map we came from.

                initialiseregion(*world, *region); // We'll do all this using the same region object as usual. We could create a new region object for it, but that seems to lead to inexplicable crashes, so we won't.

                // Now we need to prepare the images that we're going to copy the regional maps onto.

                sf::Image areareliefimage;
                sf::Image areaelevationimage;
                sf::Image areatemperatureimage;
                sf::Image areaprecipitationimage;
                sf::Image areaclimateimage;
                sf::Image areabiomeimage;
                sf::Image areariversimage;
                sf::Image areageologyimage;
                sf::Image areabasinsimage;
                sf::Image areaerosionimage;
                sf::Image areadepositionimage;
                sf::Image areafertilityimage;
                sf::Image arearesourcesimage;

                auto getAreaExportImage = [&](mapviewenum view) -> sf::Image&
                {
                    switch (view)
                    {
                    case relief: return areareliefimage;
                    case elevation: return areaelevationimage;
                    case temperature: return areatemperatureimage;
                    case precipitation: return areaprecipitationimage;
                    case climate: return areaclimateimage;
                    case biomes: return areabiomeimage;
                    case rivers: return areariversimage;
                    case geology: return areageologyimage;
                    case basins: return areabasinsimage;
                    case erosion: return areaerosionimage;
                    case deposition: return areadepositionimage;
                    case fertility: return areafertilityimage;
                    case resources: return arearesourcesimage;
                    }

                    return areareliefimage;
                };

                for (const mapviewdefinition& definition : allmapviewdefinitions)
                    getAreaExportImage(definition.view).create(imageareawidth + 1, imageareaheight + 1, sf::Color::Black);

                // Now it's time to make the regions, one by one, generate their maps, and copy those maps over onto the area images.

                for (int i = 0; i < regionswide; i++)
                {
                    int centrex = i * (int)regiontilewidth + ((int)regiontilewidth / 2) + areanwx;

                    for (int j = 0; j < regionshigh; j++)
                    {
                        int centrey = j * (int)regiontileheight + ((int)regiontileheight / 2) + areanwy;

                        // First, create the new region.

                        region->setcentrex(centrex);
                        region->setcentrey(centrey);

                        generateregionalmap(*world, *region, smalllake, island, peaks, riftblob, riftblobsize, 0, smudge, smallsmudge,squareroot);

                        // Now generate the maps.

                        resetmapcache(regionalmaps);

                        rebuildAllRegionalMaps();

                        // Now copy those maps into the images that will be exported.

                        const int destx = i * regionwidth - woffset;
                        const int desty = j * regionheight - noffset;

                        auto copyRegionalLayer = [&](sf::Image& source, sf::Image& destination)
                        {
                            const int copyx0 = max(1, destx);
                            const int copyy0 = max(1, desty);
                            const int copyx1 = min(imageareawidth - 1, destx + regionalimagewidth - 1);
                            const int copyy1 = min(imageareaheight - 1, desty + regionalimageheight - 1);

                            if (copyx0 > copyx1 || copyy0 > copyy1)
                                return;

                            destination.copy(source, copyx0, copyy0, sf::IntRect(copyx0 - destx, copyy0 - desty, copyx1 - copyx0 + 1, copyy1 - copyy0 + 1), false);
                        };

                        for (const mapviewdefinition& definition : allmapviewdefinitions)
                            copyRegionalLayer(getmapimage(regionalmaps, definition.view), getAreaExportImage(definition.view));
                    }
                }

                region->setcentrex(oldregionalcentrex); // Move the region back to where it started.
                region->setcentrey(oldregionalcentrey);

                if (areafromregional == 1) // If we're going to go back to the regional map, we need to redo it.
                    generateregionalmap(*world, *region, smalllake, island, peaks, riftblob, riftblobsize, 0, smudge, smallsmudge, squareroot);

                // Now just save the images.

                filepathname.resize(filepathname.size() - 4);

                for (const mapviewdefinition& definition : allmapviewdefinitions)
                    getAreaExportImage(definition.view).saveToFile(filepathname + " " + definition.exportstem + ".png");

                // Clean up.

                filepathname = "";
                filepath = "";

                areanex = -1;
                areaney = -1;
                areasex = -1;
                areasey = -1;
                areaswx = -1;
                areaswy = -1;
                areanwx = -1;
                areanwy = -1;

                mapview = oldmapview;

                resetmapcache(regionalmaps);

                exportingareapass = 0;

                if (areafromregional == 1)
                    screenmode = regionalmapscreen;
                else
                    screenmode = globalmapscreen;
            }
        }

        if (screenmode == generatingregionscreen)
        {
            showcolouroptions = 0;
            showworldproperties = 0;
            showglobalrainfallchart = 0;
            showglobaltemperaturechart = 0;
            showregionalrainfallchart = 0;
            showregionaltemperaturechart = 0;
            
            if (generatingregionpass < 10) // This one has a non-functioning copy of the regional map screen controls, plus the "Please wait" message.
            {
                string title;

                if (world->seed() >= 0)
                    title = "Seed: " + to_string(world->seed());
                else
                    title = "Custom";

                title = title + "##regional";

                ImGui::SetNextWindowPos(ImVec2(main_viewport->WorkPos.x + 10, main_viewport->WorkPos.y + 20), ImGuiCond_FirstUseEver);
                ImGui::SetNextWindowSize(ImVec2(161, 449), ImGuiCond_FirstUseEver);

                ImGui::Begin(title.c_str(), NULL, window_flags);

                ImGui::Text("World controls:");

                ImGui::PushItemWidth(100.0f);

                standardbutton("World map");

                ImGui::Dummy(ImVec2(0.0f, linespace));

                ImGui::SetNextItemWidth(0);

                ImGui::Text("Export options:");

                ImGui::PushItemWidth(100.0f);

                standardbutton("Regional maps");

                standardbutton("Area maps");

                ImGui::Dummy(ImVec2(0.0f, linespace));

                ImGui::SetNextItemWidth(0);

                ImGui::Text("Display map type:");

                ImGui::PushItemWidth(100.0f);

                for (const mapviewdefinition& definition : allmapviewdefinitions)
                    standardbutton(definition.label);

                ImGui::Dummy(ImVec2(0.0f, linespace));

                ImGui::SetNextItemWidth(0);

                ImGui::Text("Other controls:");

                ImGui::PushItemWidth(100.0f);

                standardbutton("Properties");

                standardbutton("Appearance");

                ImGui::End();

                // Now the text box.

                ImGui::SetNextWindowPos(ImVec2(main_viewport->WorkPos.x + 180, main_viewport->WorkPos.y + 542), ImGuiCond_FirstUseEver);
                ImGui::SetNextWindowSize(ImVec2(1023, 140), ImGuiCond_FirstUseEver);

                title = "               ";

                ImGui::Begin(title.c_str(), NULL, window_flags);

                ImGui::End();

                // Now the additional element.

                showdeferredworkwindow(generatingregionpass, "Please wait...##generateregion ", "Generating region", ImVec2(main_viewport->WorkPos.x + 849, main_viewport->WorkPos.y + 364), ImVec2(200, 60));
            }
            else
            {
                newx = newx / 32;
                newy = newy / 32;

                newx = newx * 32;
                newy = newy * 32;

                newx = newx + 16;
                newy = newy + 16;

                region->setcentrex(newx);
                region->setcentrey(newy);

                mapview = relief;

                resetmapcache(regionalmaps);

                float progressstep = 1.0f / REGIONALCREATIONSTEPS;

                // Blank the regional map image first

                regionalreliefimage.create(regionalimagewidth, regionalimageheight, sf::Color::Black);

                updateTextureFromImage(regionalmaptexture, regionalreliefimage);
                regionalmap.setTexture(regionalmaptexture);

                // Now generate the regional map

                generateregionalmap(*world, *region, smalllake, island, peaks, riftblob, riftblobsize, 0, smudge, smallsmudge, squareroot);

                // Now draw the regional map image

                showRegionalMapView(mapview);

                // Sort out the minimap

                updateTextureFromImage(globalmaptexture, displayglobalreliefimage);
                minimap.setTexture(globalmaptexture);

                focused = 0;
                generatingregionpass = 0;

                screenmode = regionalmapscreen;
            }
        }

        if (screenmode == loadfailure)
        {
            ImGui::SetNextWindowPos(ImVec2(main_viewport->WorkPos.x + 446, main_viewport->WorkPos.y + 174), ImGuiCond_FirstUseEver);
            ImGui::SetNextWindowSize(ImVec2(296, 129), ImGuiCond_FirstUseEver);

            ImGui::Begin("Loading unsuccessful!", NULL, window_flags);

            ImGui::Text("This world file is not compatible with the");
            ImGui::Text("current version of Undiscovered Worlds.");

            ImGui::Text(" ");
            ImGui::Text(" ");

            ImGui::SameLine((float)135);

            if (ImGui::Button("OK"))
            {
                if (brandnew)
                    screenmode = createworldscreen;
                else
                    screenmode = globalmapscreen;
            }

            ImGui::End();
        }

        if (screenmode == settingsloadfailure)
        {
            ImGui::SetNextWindowPos(ImVec2(main_viewport->WorkPos.x + 446, main_viewport->WorkPos.y + 174), ImGuiCond_FirstUseEver);
            ImGui::SetNextWindowSize(ImVec2(310, 129), ImGuiCond_FirstUseEver);

            ImGui::Begin("Loading unsuccessful!##settings ", NULL, window_flags);

            ImGui::Text("This settings file is not compatible with the");
            ImGui::Text("current version of Undiscovered Worlds.");

            ImGui::Text(" ");
            ImGui::Text(" ");

            ImGui::SameLine((float)135);

            if (ImGui::Button("OK"))
                screenmode = oldscreenmode;

            ImGui::End();
        }

        window.clear();

        // Colour options, if being shown

        if (showcolouroptions && screenmode != settingsloadfailure)
        {
            const int colouralign = 390;
            const int otheralign = 360;

            ImGui::SetNextWindowPos(ImVec2(main_viewport->WorkPos.x + 300, main_viewport->WorkPos.y + 20), ImGuiCond_FirstUseEver);
            ImGui::SetNextWindowSize(ImVec2(780, 640), ImGuiCond_FirstUseEver);

            ImGui::Begin("Map appearance", NULL, window_flags);

            ImGui::Checkbox("Show landmass outline", &showmapoutline);
            ImGui::SameLine(260.0f);
            ImGui::ColorEdit3("Outline", (float*)&outlinecolour);
            ImGui::SameLine(520.0f);
            ImGui::TextUnformatted("Applies to every map view.");

            if (ImGui::BeginTabBar("MapAppearanceTabs"))
            {
                for (const mapviewdefinition& definition : allmapviewdefinitions)
                {
                    if (!ImGui::BeginTabItem(definition.label))
                        continue;

                    drawmapviewappearancetab(definition, *world, appearance, selectedgradientstops, colouralign, otheralign);

                    ImGui::EndTabItem();
                }

                ImGui::EndTabBar();
            }

            ImGui::SetCursorPos(ImVec2(140.0f, 600.0f));

            if (ImGui::Button("Save", ImVec2(120.0f, 0.0f)))
            {
                openFileDialog(".uws");

                savingsettings = 1;
            }
            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Save map appearance settings.");

            ImGui::SameLine();
            if (ImGui::Button("Load", ImVec2(120.0f, 0.0f)))
            {
                openFileDialog(".uws");

                loadingsettings = 1;
            }
            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Load map appearance settings.");

            ImGui::SameLine();
            if (ImGui::Button("Default", ImVec2(120.0f, 0.0f)))
            {
                initialisemapcolours(*world);
                syncappearancesettings(*world, appearance);
                refreshhighlightdisplay(*world, highlightdisplay);

                colourschanged = 1;
            }
            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Restore the default map appearance settings.");

            ImGui::SameLine();
            if (ImGui::Button("Close", ImVec2(120.0f, 0.0f)))
            {
                showcolouroptions = 0;
            }
            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Close the map appearance panel.");

            ImGui::End();
        }

        // World properties window, if being shown

        if (showworldproperties)
        {
            int rightalign = 290;

            int topleftfigures = 105;
            int bottomleftfigures = 225;
            int toprightfigures = 410;
            int bottomrightfigures = 460;

            ImGui::SetNextWindowPos(ImVec2(main_viewport->WorkPos.x + 265, main_viewport->WorkPos.y + 60), ImGuiCond_FirstUseEver);
            ImGui::SetNextWindowSize(ImVec2(640, 500), ImGuiCond_FirstUseEver);

            ImGui::Begin("World properties", NULL, window_flags & ~ImGuiWindowFlags_NoResize);

            std::array<int, 32> climateareas = {};
            int totallandarea = 0;

            for (int i = 0; i <= world->width(); i++)
            {
                for (int j = 0; j <= world->height(); j++)
                {
                    short thisclimate = world->climate(i, j);

                    if (thisclimate > 0 && thisclimate < static_cast<short>(climateareas.size()))
                    {
                        climateareas[thisclimate]++;
                        totallandarea++;
                    }
                }
            }

            if (ImGui::BeginTabBar("WorldPropertiesTabs"))
            {
                if (ImGui::BeginTabItem("Overview"))
                {

            string sizeinfo = "Size:";

            string sizevalue = "";

            if (world->size() == 0)
                sizevalue = "small";

            if (world->size() == 1)
                sizevalue = "medium";

            if (world->size() == 2)
                sizevalue = "large";

            stringstream ss3;
            ss3 << fixed << setprecision(5) << world->eccentricity();

            string eccentricityinfo = "Eccentricity:";
            string eccentricityvalue = ss3.str();

            stringstream ss4;
            ss4 << fixed << setprecision(2) << world->gravity();

            string gravityinfo = "Surface gravity: ";
            string gravityvalue = ss4.str() + "g";

            string perihelioninfo = "Perihelion:";
            string perihelionvalue = "";

            if (world->perihelion() == 0)
                perihelionvalue = "January";

            if (world->perihelion() == 1)
                perihelionvalue = "July";

            stringstream ss5;
            ss5 << fixed << setprecision(2) << world->lunar();

            string lunarinfo = "Lunar pull:";
            string lunarvalue = ss5.str();

            string typeinfo = "Category: ";
            string typevalue = to_string(world->type());

            if (typevalue == "1")
                typevalue = "tectonic (small)";

            if (typevalue == "2")
                typevalue = "tectonic (large)";

            if (typevalue == "3")
                typevalue = "oceanic";

            if (typevalue == "4")
                typevalue = "non-tectonic";

            string rotationinfo = "Rotation:";
            string rotationvalue = "";

            if (world->rotation())
                rotationvalue = "west to east";
            else
                rotationvalue = "east to west";

            stringstream ss;
            ss << fixed << setprecision(2) << world->tilt();

            string tiltinfo = "Obliquity:";
            string tiltvalue = ss.str() + degree;

            stringstream ss2;
            ss2 << fixed << setprecision(2) << world->tempdecrease();

            string tempdecreaseinfo = "Heat decrease per vertical km:";
            string tempdecreasevalue = ss2.str() + degree;

            string northpoleinfo = "North pole adjustment:";
            string northpolevalue = to_string(world->northpolaradjust()) + degree;

            string southpoleinfo = "South pole adjustment:";
            string southpolevalue = to_string(world->southpolaradjust()) + degree;

            string averageinfo = "Average global temperature:";
            string averagevalue = to_string(world->averagetemp()) + degree;

            stringstream ss6;
            ss6 << fixed << setprecision(2) << world->waterpickup();

            string moistureinfo = "Moisture pickup rate:";
            string moisturevalue = ss6.str();

            string glacialinfo = "Glaciation temperature:";
            string glacialvalue=to_string(world->glacialtemp())+degree;

            ImGui::Text(sizeinfo.c_str());

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Size of planet. (Earth: large; Mars: medium; Moon: small)");

            ImGui::SameLine((float)topleftfigures);
            ImGui::Text(sizevalue.c_str());

            ImGui::SameLine((float)rightalign);

            ImGui::Text(gravityinfo.c_str());

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Affects mountain and valley sizes. (Earth: 1.00g)");

            ImGui::SameLine((float)toprightfigures);
            ImGui::Text(gravityvalue.c_str());

            ImGui::Text(typeinfo.c_str());

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Terrain category. (Earth: large tectonic)");

            ImGui::SameLine((float)topleftfigures);
            ImGui::Text(typevalue.c_str());

            ImGui::SameLine((float)rightalign);

            ImGui::Text(lunarinfo.c_str());

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Affects tides and coastal regions. (Earth: 1.00)");

            ImGui::SameLine((float)toprightfigures);
            ImGui::Text(lunarvalue.c_str());

            ImGui::Text(rotationinfo.c_str());

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Affects weather patterns. (Earth: west to east)");

            ImGui::SameLine((float)topleftfigures);
            ImGui::Text(rotationvalue.c_str());

            ImGui::SameLine((float)rightalign);

            ImGui::Text(tiltinfo.c_str());

            string tilttip = "Affects seasonal variation in temperature. (Earth: 22.5" + degree + ")";

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip(tilttip.c_str());

            ImGui::SameLine((float)toprightfigures);
            ImGui::Text(tiltvalue.c_str());

            ImGui::Text(eccentricityinfo.c_str());

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("How elliptical the orbit is. (Earth: 0.0167)");

            ImGui::SameLine((float)topleftfigures);
            ImGui::Text(eccentricityvalue.c_str());

            ImGui::SameLine((float)rightalign);

            ImGui::Text(perihelioninfo.c_str());

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("When the planet is closest to the sun. (Earth: January)");

            ImGui::SameLine((float)toprightfigures);
            ImGui::Text(perihelionvalue.c_str());

            ImGui::Text("   ");

            ImGui::Text(tempdecreaseinfo.c_str());

            string tempdecreasetip = "Affects how much colder it gets higher up. (Earth: 6.5" + degree + ")";

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip(tempdecreasetip.c_str());

            ImGui::SameLine((float)bottomleftfigures);
            ImGui::Text(tempdecreasevalue.c_str());

            ImGui::SameLine((float)rightalign);

            ImGui::Text(moistureinfo.c_str());

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Affects how much moisture wind picks up from the ocean. (Earth: 1.0)");

            ImGui::SameLine((float)bottomrightfigures);
            ImGui::Text(moisturevalue.c_str());

            ImGui::Text("   ");

            ImGui::Text(averageinfo.c_str());

            string avetip = "Earth: 14" + degree;

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip(avetip.c_str());

            ImGui::SameLine((float)bottomleftfigures);
            ImGui::Text(averagevalue.c_str());

            ImGui::SameLine((float)rightalign);

            ImGui::Text(glacialinfo.c_str());

            string glacialtip = "Areas below this average temperature may show signs of past glaciation. (Earth: 4" + degree + ")";

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip(glacialtip.c_str());

            ImGui::SameLine((float)bottomrightfigures);
            ImGui::Text(glacialvalue.c_str());

            ImGui::Text(northpoleinfo.c_str());

            string northtip = "Adjustment to north pole temperature. (Earth: +3" + degree + ")";

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip(northtip.c_str());

            ImGui::SameLine((float)bottomleftfigures);
            ImGui::Text(northpolevalue.c_str());

            ImGui::Text(southpoleinfo.c_str());

            string southtip = "Adjustment to south pole temperature. (Earth: -3" + degree + ")";

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip(southtip.c_str());

            ImGui::SameLine((float)bottomleftfigures);
            ImGui::Text(southpolevalue.c_str());
                    ImGui::EndTabItem();
                }

                if (ImGui::BeginTabItem("Climate"))
                {
                    ImGui::Text("Percent of land area");
                    ImGui::Text(" ");

                    if (ImGui::BeginTable("ClimateAreas", 3, ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersInnerV | ImGuiTableFlags_Resizable | ImGuiTableFlags_ScrollY, ImVec2(0.0f, 360.0f)))
                    {
                        ImGui::TableSetupColumn("Code", ImGuiTableColumnFlags_WidthFixed, 55.0f);
                        ImGui::TableSetupColumn("Climate");
                        ImGui::TableSetupColumn("Area %", ImGuiTableColumnFlags_WidthFixed, 80.0f);
                        ImGui::TableHeadersRow();

                        for (int climate = 1; climate <= world->climatenumber() && climate < static_cast<int>(climateareas.size()); climate++)
                        {
                            if (climateareas[climate] == 0)
                                continue;

                            const float percentage = totallandarea > 0 ? (static_cast<float>(climateareas[climate]) * 100.0f) / static_cast<float>(totallandarea) : 0.0f;

                            ImGui::TableNextRow();
                            ImGui::TableSetColumnIndex(0);
                            ImGui::Text("%s", getclimatecode((short)climate).c_str());
                            ImGui::TableSetColumnIndex(1);
                            ImGui::Text("%s", getclimatename((short)climate).c_str());
                            ImGui::TableSetColumnIndex(2);
                            ImGui::Text("%.2f", percentage);
                        }

                        ImGui::EndTable();
                    }

                    ImGui::EndTabItem();
                }

                ImGui::EndTabBar();
            }

            if (ImGui::Button("Close", ImVec2(120.0f, 0.0f)))
                showworldproperties = 0;

            ImGui::End();
        }

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
                invalidateGlobalMaps();
                showGlobalMapView(relief);
            }
        });

        handlefiledialog(importingseamap, filepathname, filepath, [&](const std::string& selectedpath, const std::string&)
        {
            if (importseadepthmap(*world, selectedpath, &mapimportsettings))
            {
                invalidateGlobalMaps();
                showGlobalMapView(relief);
            }
        });

        handlefiledialog(importingmountainsmap, filepathname, filepath, [&](const std::string& selectedpath, const std::string&)
        {
            if (importmountainmap(*world, selectedpath, OKmountains, &mapimportsettings))
            {
                invalidateGlobalMaps();
                showGlobalMapView(relief);
            }
        });

        handlefiledialog(importingvolcanoesmap, filepathname, filepath, [&](const std::string& selectedpath, const std::string&)
        {
            if (importvolcanomap(*world, selectedpath, &mapimportsettings))
            {
                invalidateGlobalMaps();
                showGlobalMapView(relief);
            }
        });

        handlefiledialog(importingtemperaturemap, filepathname, filepath, [&](const std::string& selectedpath, const std::string&)
        {
            if (importtemperaturemap(*world, selectedpath, importedclimatemaps, &mapimportsettings))
            {
                invalidateGlobalMaps();
                showGlobalMapView(temperature);
            }
        });

        handlefiledialog(importingprecipitationmap, filepathname, filepath, [&](const std::string& selectedpath, const std::string&)
        {
            if (importprecipitationmap(*world, selectedpath, importedclimatemaps, &mapimportsettings))
            {
                invalidateGlobalMaps();
                showGlobalMapView(precipitation);
            }
        });

        handlefiledialog(importinggradientstrip, filepathname, filepath, [&](const std::string& selectedpath, const std::string&)
        {
            mapimportsettings.gradientStripPath = selectedpath;
        });

        // Set size window, for custom worlds.

        if (drawcustomworldsizewindow(main_viewport, window_flags, showsetsize, brandnew, currentsize))
        {
            initialiseworld(*world);
            world->clear();
            clearimportedclimatemaps(importedclimatemaps);
            mapimportsettings = MapImportSettings{};

            world->setsize(currentsize);

            resizeglobaldisplayforworld(*world, globaldisplay, highlightdisplay);
            world->setgravity(getdefaultgravityforsize(currentsize));

            landmass = 5;
            mergefactor = 15;
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

        // Window for creating custom tectonic-based worlds.

        if (drawtectonicchooserwindow(main_viewport, window_flags, showtectonicchooser, landmass, mergefactor))
        {
            screenmode = generatingtectonicscreen;
            generatingtectonicpass = 0;
        }

        // Window for creating custom non-tectonic-based worlds.

        if (drawnontectonicchooserwindow(main_viewport, window_flags, shownontectonicchooser, sealeveleditable, iterations))
        {
            screenmode = generatingnontectonicscreen;
            generatingnontectonicpass = 0;
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

// This looks up the latest version of the program. Based on code by Parveen: https://hoven.in/cpp-network/c-program-download-file-from-url.html

float getlatestversion()
{
    IStream* stream;

    const char* URL = "https://raw.githubusercontent.com/JonathanCRH/Undiscovered_Worlds/main/version.txt";

    if (getURL(0, URL, &stream, 0, 0)) // Didn't work for some reason.
        return 0.0f;

    // this char array will be cyclically filled with bytes from URL
    char buff[100];

    // we shall keep appending the bytes to this string
    string s;

    unsigned long bytesRead;

    while (true)
    {
        // Reads a specified number of bytes from the stream object into char array and stores the actual bytes read to "bytesRead"
        stream->Read(buff, 100, &bytesRead);

        if (0U == bytesRead)
            break;

        // append and collect to the string
        s.append(buff, bytesRead);
    };

    // release the interface
    stream->Release();

    float val = stof(s);

    return val;
}
