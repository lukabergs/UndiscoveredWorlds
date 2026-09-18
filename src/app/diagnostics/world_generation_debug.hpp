#pragma once

#include <functional>
#include <string>

#include "social_generation.hpp"

struct WorldGenerationDebugOptions
{
    bool logToProfilingWorkbook = true;
    bool visualizePlateTectonicsRealtime = false;
    int plateTectonicsCycleCount = 2;
    int plateTectonicsCycleStepLimit = 600;
    int plateTectonicsPlateCount = 10;
    bool plateTectonicsUseSeaLevelMeters = false;
    int plateTectonicsSeaLevelMeters = 31043;
    int plateTectonicsAggregationOverlapAbsolute = -1;
    float plateTectonicsAggregationOverlapRelative = 0.20f;
    float plateTectonicsFoldingRatio = 0.08f;
    int plateTectonicsErosionPeriod = 60;
    float plateTectonicsErosionStrength = 1.0f;
    float plateTectonicsLandmassRotation = 0.20f;
    float plateTectonicsRotationStrength = 1.0f;
    float plateTectonicsSubductionStrength = 1.0f;
    float plateTectonicsDivergentCarveStrength = 0.015f;
    float plateTectonicsDeltaTimeMyr = 1.0f;
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
int platetectonicscyclecount();
int platetectonicsplatecount();
