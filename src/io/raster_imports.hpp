#pragma once

// World-independent raster decoding for map imports: grayscale 16-bit TIFF or
// PNG (red channel or gradient strip), checked against the world grid size.
// map_imports.hpp wraps these for planet.

#include <cstdint>
#include <functional>
#include <string>
#include <unordered_map>
#include <vector>

#include "io/raster_image.hpp"

struct MapImportSettings;

namespace rasterimports
{
struct GradientStripDecoder
{
    bool enabled = false;
    float minimum = 0.0f;
    float increment = 1.0f;
    std::vector<rasterimage::Color> colours;
    std::unordered_map<unsigned int, int> exactmatches;
};

struct ImportedRaster
{
    bool isTiff16 = false;
    rasterimage::Image image;
    std::vector<uint16_t> grayscale16;
};

// width and height are the largest indices of the world grid.
bool hastiffextension(const std::string& filepathname);
bool loadworldsizedimage(int width, int height, const std::string& filepathname, rasterimage::Image& importimage);
bool loadworldsizedraster(int width, int height, const std::string& filepathname, ImportedRaster& raster);
uint16_t rastervalueat(const ImportedRaster& raster, int width, int x, int y);
int roundimportvalue(float value);
bool loadgradientstripdecoder(const MapImportSettings* settings, GradientStripDecoder& decoder);
float decodegradientvalue(const rasterimage::Color& colour, const GradientStripDecoder& decoder);
int importtemperaturevalue(rasterimage::Uint8 value);
int importprecipitationvalue(rasterimage::Uint8 value);
}

// Land heights (metres above sea level) and sea depths (metres below) on the
// datum used by planet: setnom receives sealevel +/- the imported value.
bool importlandheightmapcells(const std::string& filepathname, int width, int height, int sealevel,
    const MapImportSettings* settings, const std::function<void(int, int, int)>& setnom);
bool importseadepthmapcells(const std::string& filepathname, int width, int height, int sealevel,
    const MapImportSettings* settings, const std::function<void(int, int, int)>& setnom);
