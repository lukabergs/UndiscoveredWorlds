#pragma once

#include <filesystem>
#include <string>

struct AppEnvironmentConfig
{
    std::filesystem::path defaultWorldDirectory = ".";
    std::filesystem::path defaultAppearanceDirectory = "assets/appearance";
    std::filesystem::path defaultImageDirectory = "out/images";
    // Shape templates of the terrain module (simulations/terrain/shape_templates.cpp).
    std::filesystem::path terrainTemplatePath = "assets/terrain/shape-templates.txt";
    // Tables of the climate solvers written by a run (simulations/climate/climate_diagnostics.cpp).
    std::filesystem::path climateDiagnosticsDirectory = "out/climate-diagnostics";
    // Earth terrain of the imported recipes; resolution variants are resolved beside these files.
    std::filesystem::path earthBenchmarkLandPath;
    std::filesystem::path earthBenchmarkSeaPath;
};

const AppEnvironmentConfig& getappenvironment();
void reloadappenvironment();

// Moves the working directory to the repository root (the folder with version.txt)
// when started elsewhere, searching upwards from the executable.
void ensureworkingdirectory();
// The _WIDTHxHEIGHT variant of a reference file when it exists, else the file itself.
std::filesystem::path resolutionvariant(const std::filesystem::path& source, int columns, int rows);
