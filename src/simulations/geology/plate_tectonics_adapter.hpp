#pragma once

#include <vector>

class planet;

// Application-facing configuration for the bundled plate tectonics simulation.
// Defaults and overloads match the existing world-generation contract.
struct PlateTectonicsSimulationOptions
{
    int cycleCount = 2;
    int cycleStepLimit = 600;
    int plateCount = 10;
    bool useSeaLevelMeters = false;
    int seaLevelMeters = 31043;
    int aggregationOverlapAbsolute = -1;
    float aggregationOverlapRelative = 0.20f;
    float foldingRatio = 0.08f;
    int erosionPeriod = 60;
    float erosionStrength = 1.0f;
    float landmassRotation = 0.20f;
    float rotationStrength = 1.0f;
    float subductionStrength = 1.0f;
    float divergentCarveStrength = 0.015f;
    float deltaTimeMyr = 1.0f;
};

int defaultplatetectonicsaggregationoverlapabs(int width, int height);
void applyplatetectonicssimulation(planet& world, std::vector<std::vector<bool>>& shelves, int cyclecount, int platecount);
void applyplatetectonicssimulation(planet& world, std::vector<std::vector<bool>>& shelves, const PlateTectonicsSimulationOptions& options);
