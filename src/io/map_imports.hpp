#pragma once

#include <string>
#include <vector>

#include "planet.hpp"

enum class MapImportValueMode
{
    redChannel = 0,
    gradientStrip = 1
};

struct MapImportSettings
{
    MapImportValueMode valueMode = MapImportValueMode::redChannel;
    std::string gradientStripPath;
    float gradientMinimum = 0.0f;
    float gradientIncrement = 1.0f;
};

struct ImportedClimateMaps
{
    int width = 0;
    int height = 0;
    bool hasTemperature = false;
    bool hasPrecipitation = false;
    std::vector<float> annualTemperature;
    std::vector<float> annualPrecipitation;
};

void clearimportedclimatemaps(ImportedClimateMaps& maps);
bool importlandheightmap(planet& world, const std::string& filepathname, const MapImportSettings* settings = nullptr);
bool importseadepthmap(planet& world, const std::string& filepathname, const MapImportSettings* settings = nullptr);
bool importmountainmap(planet& world, const std::string& filepathname, std::vector<std::vector<bool>>& okmountains, const MapImportSettings* settings = nullptr);
bool importvolcanomap(planet& world, const std::string& filepathname, const MapImportSettings* settings = nullptr);
bool importtemperaturemap(planet& world, const std::string& filepathname, ImportedClimateMaps& maps, const MapImportSettings* settings = nullptr);
bool importprecipitationmap(planet& world, const std::string& filepathname, ImportedClimateMaps& maps, const MapImportSettings* settings = nullptr);
