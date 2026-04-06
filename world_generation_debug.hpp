#pragma once

#include <functional>
#include <string>
#include <vector>

#include "social_generation.hpp"

struct WorldGenerationDebugOptions
{
    bool logToProfilingWorkbook = true;
    bool visualizeEachStep = false;
    bool useFastLEMMountains = false;
    bool usePlateTectonicsSimulation = false;
    int plateTectonicsCycleCount = 4;
    bool socialEnabled = false;
    SocialGenerationOptions::Mode socialMode = SocialGenerationOptions::Mode::static_ex_nihilo;
    bool usePrehistory = true;
    int historyYears = 1200;
    std::vector<bool> enabledSteps;

    WorldGenerationDebugOptions();
};

const std::vector<std::string>& getworldgenerationstepoptions();
bool beginworldgenstep(const char* label);
void beginworldgendebugrun(long seed, const WorldGenerationDebugOptions* options);
void endworldgendebugrun();
void onworldgenstepcompleted(const std::string& label, double elapsedms);
void setworldgenvisualizationcallback(std::function<void()> callback);
void clearworldgenvisualizationcallback();
bool isworldgendebugrunactive();
long worldgenerationdebugseed();
bool usefastlemmountains();
bool useplatetectonicssimulation();
int platetectonicscyclecount();
