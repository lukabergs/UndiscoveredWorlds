#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include "core/app_environment.hpp"
#include "pipeline/generation_progress.hpp"
#include "pipeline/world_generation_debug.hpp"

using namespace std;

namespace
{
const WorldGenerationDebugOptions* activeworldgenerationoptions = nullptr;
long activeworldgenerationseed = 0;
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
}

void endworldgendebugrun()
{
    activeworldgenerationoptions = nullptr;
    activeworldgenerationseed = 0;
    worldgenerationvisualizationcallback = nullptr;
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

