#pragma once

#include <filesystem>
#include <string>

struct AppEnvironmentConfig
{
    std::filesystem::path defaultWorldDirectory = ".";
    std::filesystem::path defaultAppearanceDirectory = LR"(C:\dev\UndiscoveredWorlds\extra\appearance)";
    std::filesystem::path defaultImageDirectory = LR"(C:\dev\UndiscoveredWorlds\extra\img)";
    std::filesystem::path profilingWorkbookPath = LR"(C:\dev\UndiscoveredWorlds\extra\profiling.xlsx)";
};

const AppEnvironmentConfig& getappenvironment();
void reloadappenvironment();
