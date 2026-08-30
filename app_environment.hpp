#pragma once

#include <filesystem>
#include <string>

struct AppEnvironmentConfig
{
    std::filesystem::path defaultWorldDirectory = ".";
    std::filesystem::path defaultAppearanceDirectory = LR"(extra\appearance)";
    std::filesystem::path defaultImageDirectory = LR"(extra\img)";
    std::filesystem::path profilingWorkbookPath = LR"(extra\profiling.xlsx)";
    std::filesystem::path referencePrecipitationGridPath = LR"(extra\reference\earth_precipitation_grid.csv)";
    std::filesystem::path referenceClimateDirectory = LR"(extra\reference)";
    std::filesystem::path climateWorkbookPath = LR"(extra\climate.xlsx)";
    std::filesystem::path climateBenchmarkRunLogPath = "climate_benchmark_runs.json";
    std::filesystem::path climateBenchmarkImageDirectory = LR"(extra\img\earth\benchmark)";
    std::filesystem::path earthKoppenImagePath = LR"(extra\img\earth\in\earth_l_koppen.png)";
    std::filesystem::path earthBenchmarkLandPath = LR"(extra\img\earth\in\earth_land_l_3.png)";
    std::filesystem::path earthBenchmarkSeaPath = LR"(extra\img\earth\in\earth_sea_l_1.png)";
};

const AppEnvironmentConfig& getappenvironment();
void reloadappenvironment();
void setreferenceprecipitationgridpath(const std::filesystem::path& path);
