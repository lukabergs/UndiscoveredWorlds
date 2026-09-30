#pragma once

#include <filesystem>
#include <string>

struct AppEnvironmentConfig
{
    std::filesystem::path defaultWorldDirectory = ".";
    std::filesystem::path defaultAppearanceDirectory = "assets/appearance";
    std::filesystem::path defaultImageDirectory = "runs/maps";
    std::filesystem::path profilingWorkbookPath = "runs/metrics/profiling.xlsx";
    std::filesystem::path referencePrecipitationGridPath = "refs/processed/climate/earth_precipitation_grid.csv";
    std::filesystem::path referenceClimateDirectory = "refs/processed/climate";
    std::filesystem::path referenceClimatePreviewDirectory = "refs/processed/climate/previews";
    std::filesystem::path climateWorkbookPath = "runs/metrics/climate.xlsx";
    std::filesystem::path climateBenchmarkRunLogPath = "runs/registry/climate.json";
    // Organized previews and durable run diagnostics; seed-scoped working files are separate.
    std::filesystem::path climateBenchmarkImageDirectory = "runs/maps";
    std::filesystem::path climateBenchmarkRunDirectory = "runs/diagnostics/climate";
    std::filesystem::path climateValidationDirectory = "runs/work/climate";
    // The old unsigned benchmark terrain and categorical reference were removed.
    // Supply compatible inputs through configs/app.env or environment overrides.
    std::filesystem::path earthKoppenImagePath;
    std::filesystem::path earthBenchmarkLandPath;
    std::filesystem::path earthBenchmarkSeaPath;
};

const AppEnvironmentConfig& getappenvironment();
void reloadappenvironment();
void setreferenceprecipitationgridpath(const std::filesystem::path& path);
