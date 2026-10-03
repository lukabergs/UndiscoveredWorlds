#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include "app_environment.hpp"
#include "generation_progress.hpp"
#include "profiling_workbook.hpp"
#include "world_generation_debug.hpp"

using namespace std;

namespace
{
const WorldGenerationDebugOptions* activeworldgenerationoptions = nullptr;
long activeworldgenerationseed = 0;
vector<pair<string, double>> activeworldgenerationtimings;
function<void()> worldgenerationvisualizationcallback;
}

bool beginworldgenstep(const char* label)
{
    updatereport(label);
    return true;
}

void beginworldgendebugrun(long seed, const WorldGenerationDebugOptions* options)
{
    activeworldgenerationseed = seed;
    activeworldgenerationoptions = options;
    activeworldgenerationtimings.clear();
}

void endworldgendebugrun()
{
    if (activeworldgenerationoptions != nullptr && activeworldgenerationoptions->logToProfilingWorkbook)
        appendprofilingworkbookrow(getappenvironment().profilingWorkbookPath, activeworldgenerationseed, activeworldgenerationtimings);

    activeworldgenerationtimings.clear();
    activeworldgenerationoptions = nullptr;
    activeworldgenerationseed = 0;
    worldgenerationvisualizationcallback = nullptr;
}

void onworldgenstepcompleted(const string& label, double elapsedms)
{
    if (activeworldgenerationoptions == nullptr)
        return;

    activeworldgenerationtimings.emplace_back(label, elapsedms);
}

void setworldgenvisualizationcallback(function<void()> callback)
{
    worldgenerationvisualizationcallback = std::move(callback);
}

void clearworldgenvisualizationcallback()
{
    worldgenerationvisualizationcallback = nullptr;
}

bool hasworldgenvisualizationcallback()
{
    return static_cast<bool>(worldgenerationvisualizationcallback);
}

void requestworldgenvisualization()
{
    if (worldgenerationvisualizationcallback)
        worldgenerationvisualizationcallback();
}

bool isworldgendebugrunactive()
{
    return activeworldgenerationoptions != nullptr;
}

long worldgenerationdebugseed()
{
    return activeworldgenerationseed;
}

double tectonicdurationmyr()
{
    return activeworldgenerationoptions ? activeworldgenerationoptions->plateTectonicsDurationMyr : 0;
}

int platetectonicsplatecount()
{
    return activeworldgenerationoptions ? activeworldgenerationoptions->plateTectonicsPlateCount : 0;
}

