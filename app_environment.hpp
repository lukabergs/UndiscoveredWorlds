#pragma once

#include <filesystem>
#include <string>

struct AppEnvironmentConfig
{
    std::filesystem::path defaultWorldDirectory = ".";
    std::filesystem::path defaultAppearanceDirectory = LR"(extra\appearance)";
    std::filesystem::path defaultImageDirectory = LR"(extra\img)";
    std::filesystem::path profilingWorkbookPath = LR"(extra\climate\workbooks\profiling.xlsx)";
    std::filesystem::path referencePrecipitationGridPath = LR"(extra\reference\climate\processed\earth_precipitation_grid.csv)";
    std::filesystem::path referenceClimateDirectory = LR"(extra\reference\climate\processed)";
    std::filesystem::path referenceClimatePreviewDirectory = LR"(extra\reference\climate\previews)";
    std::filesystem::path climateWorkbookPath = LR"(extra\climate\workbooks\climate.xlsx)";
    std::filesystem::path climateBenchmarkRunLogPath = LR"(data\climate\benchmark_runs.json)";
    std::filesystem::path climateBenchmarkImageDirectory = LR"(extra\climate\benchmarks\maps)";
    std::filesystem::path climateBenchmarkRunDirectory = LR"(extra\climate\benchmarks\runs)";
    std::filesystem::path climateValidationDirectory = LR"(extra\climate\benchmarks\work)";
    std::filesystem::path earthKoppenImagePath = LR"(extra\reference\earth\base-maps\earth_l_koppen.png)";
    std::filesystem::path earthBenchmarkLandPath = LR"(extra\reference\earth\base-maps\earth_land_l_3.png)";
    std::filesystem::path earthBenchmarkSeaPath = LR"(extra\reference\earth\base-maps\earth_sea_l_1.png)";
};

const AppEnvironmentConfig& getappenvironment();
void reloadappenvironment();
void setreferenceprecipitationgridpath(const std::filesystem::path& path);
