#pragma once

#include <filesystem>
#include <string>

struct AppEnvironmentConfig
{
    std::filesystem::path defaultWorldDirectory = ".";
    std::filesystem::path defaultAppearanceDirectory = LR"(C:\dev\UndiscoveredWorlds\extra\appearance)";
    std::filesystem::path defaultImageDirectory = LR"(C:\dev\UndiscoveredWorlds\extra\img)";
    std::filesystem::path profilingWorkbookPath = LR"(C:\dev\UndiscoveredWorlds\extra\profiling.xlsx)";
    std::filesystem::path referencePrecipitationGridPath = LR"(C:\dev\UndiscoveredWorlds\extra\reference\earth_precipitation_grid.csv)";
    std::filesystem::path climateWorkbookPath = LR"(C:\dev\UndiscoveredWorlds\extra\climate.xlsx)";
    std::filesystem::path earthKoppenImagePath = LR"(C:\dev\UndiscoveredWorlds\extra\img\earth\in\earth_l_koppen.png)";
    std::filesystem::path earthBenchmarkLandPath = LR"(C:\dev\UndiscoveredWorlds\extra\img\earth\in\earth_land_l_3.png)";
    std::filesystem::path earthBenchmarkSeaPath = LR"(C:\dev\UndiscoveredWorlds\extra\img\earth\in\earth_sea_l_1.png)";
};

const AppEnvironmentConfig& getappenvironment();
void reloadappenvironment();
void setreferenceprecipitationgridpath(const std::filesystem::path& path);
