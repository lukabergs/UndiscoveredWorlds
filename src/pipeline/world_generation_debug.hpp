#pragma once

// Generation session hooks: progress labels, step timings for the profiling
// workbook, tectonics visualization callbacks and session options.

#include <functional>
#include <string>

#include "social_generation.hpp"

struct WorldGenerationDebugOptions
{
    bool logToProfilingWorkbook = true;
    bool visualizePlateTectonicsRealtime = false;
    double plateTectonicsDurationMyr = 0;
    int plateTectonicsPlateCount = 0;
    float plateTectonicsDeltaTimeMyr = 0.0f;
    bool socialEnabled = false;
    SocialGenerationOptions::Mode socialMode = SocialGenerationOptions::Mode::static_ex_nihilo;
    bool usePrehistory = true;
    int historyYears = 1200;
};

bool beginworldgenstep(const char* label);
void beginworldgendebugrun(long seed, const WorldGenerationDebugOptions* options);
void endworldgendebugrun();
void onworldgenstepcompleted(const std::string& label, double elapsedms);
void setworldgenvisualizationcallback(std::function<void()> callback);
void clearworldgenvisualizationcallback();
bool hasworldgenvisualizationcallback();
void requestworldgenvisualization();
bool isworldgendebugrunactive();
long worldgenerationdebugseed();
double tectonicdurationmyr();
int platetectonicsplatecount();
