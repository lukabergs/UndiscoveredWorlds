#pragma once

#include <functional>
#include <string>
#include <vector>

struct WorldGenerationDebugOptions
{
    bool logToProfilingWorkbook = true;
    bool visualizeEachStep = false;
    bool useFastLEMMountains = false;
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
bool usefastlemmountains();
