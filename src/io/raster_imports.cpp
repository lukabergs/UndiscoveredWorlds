#include "raster_imports.hpp"

#include <cctype>
#include <cmath>
#include <filesystem>
#include <limits>

#include "map_imports.hpp"
#include "topography_tiff.hpp"

using namespace std;
namespace rasterimports
{
    namespace fs = std::filesystem;

    constexpr float importedtemperatureminimum = -60.0f;
    constexpr float importedtemperaturemaximum = 60.0f;
    constexpr float importedprecipitationmaximum = 1020.0f;





    namespace
    {
    std::string lowercasestring(std::string value)
    {
        for (char& ch : value)
            ch = static_cast<char>(tolower(static_cast<unsigned char>(ch)));

        return value;
    }

    }

    bool hastiffextension(const string& filepathname)
    {
        const std::string extension = lowercasestring(fs::path(filepathname).extension().string());
        return extension == ".tif" || extension == ".tiff";
    }

    bool loadworldsizedimage(int width, int height, const string& filepathname, rasterimage::Image& importimage)
    {
        if (!importimage.loadFromFile(filepathname))
            return false;

        const rasterimage::Size imagesize = importimage.getSize();

        return imagesize.x == static_cast<unsigned int>(width + 1) && imagesize.y == static_cast<unsigned int>(height + 1);
    }

    bool loadworldsizedtiff16(int worldwidth, int worldheight, const string& filepathname, std::vector<uint16_t>& importvalues)
    {
        uint32_t width = 0;
        uint32_t height = 0;

        if (!topographyio::readGrayscaleTiff16(filepathname.c_str(), importvalues, width, height))
            return false;

        return width == static_cast<uint32_t>(worldwidth + 1)
            && height == static_cast<uint32_t>(worldheight + 1);
    }

    bool loadworldsizedraster(int width, int height, const string& filepathname, ImportedRaster& raster)
    {
        raster = {};

        if (hastiffextension(filepathname))
        {
            raster.isTiff16 = true;
            return loadworldsizedtiff16(width, height, filepathname, raster.grayscale16);
        }

        return loadworldsizedimage(width, height, filepathname, raster.image);
    }

    uint16_t rastervalueat(const ImportedRaster& raster, int width, int x, int y)
    {
        return raster.grayscale16[static_cast<size_t>(y) * static_cast<size_t>(width + 1) + static_cast<size_t>(x)];
    }

    

    unsigned int packrgb(const rasterimage::Color& colour)
    {
        return static_cast<unsigned int>(colour.r)
            | (static_cast<unsigned int>(colour.g) << 8)
            | (static_cast<unsigned int>(colour.b) << 16);
    }

    int roundimportvalue(float value)
    {
        return static_cast<int>(roundf(value));
    }

    bool loadgradientstripdecoder(const MapImportSettings* settings, GradientStripDecoder& decoder)
    {
        decoder = {};

        if (settings == nullptr || settings->valueMode == MapImportValueMode::redChannel)
            return true;

        if (settings->gradientStripPath.empty())
            return false;

        rasterimage::Image gradientimage;

        if (!gradientimage.loadFromFile(settings->gradientStripPath))
            return false;

        const rasterimage::Size gradientsize = gradientimage.getSize();

        if (gradientsize.y != 1 || gradientsize.x == 0)
            return false;

        decoder.enabled = true;
        decoder.minimum = settings->gradientMinimum;
        decoder.increment = settings->gradientIncrement;
        decoder.colours.reserve(gradientsize.x);

        for (unsigned int x = 0; x < gradientsize.x; x++)
        {
            const rasterimage::Color colour = gradientimage.getPixel(x, 0);
            decoder.colours.push_back(colour);
            decoder.exactmatches.emplace(packrgb(colour), static_cast<int>(x));
        }

        return true;
    }

    float decodegradientvalue(const rasterimage::Color& colour, const GradientStripDecoder& decoder)
    {
        if (!decoder.enabled)
            return static_cast<float>(colour.r);

        const auto exactmatch = decoder.exactmatches.find(packrgb(colour));

        if (exactmatch != decoder.exactmatches.end())
            return decoder.minimum + static_cast<float>(exactmatch->second) * decoder.increment;

        int nearestindex = 0;
        long bestdistance = std::numeric_limits<long>::max();

        for (int index = 0; index < static_cast<int>(decoder.colours.size()); index++)
        {
            const rasterimage::Color candidate = decoder.colours[index];
            const int dr = static_cast<int>(colour.r) - static_cast<int>(candidate.r);
            const int dg = static_cast<int>(colour.g) - static_cast<int>(candidate.g);
            const int db = static_cast<int>(colour.b) - static_cast<int>(candidate.b);
            const long distance = static_cast<long>(dr * dr + dg * dg + db * db);

            if (distance < bestdistance)
            {
                bestdistance = distance;
                nearestindex = index;
            }
        }

        return decoder.minimum + static_cast<float>(nearestindex) * decoder.increment;
    }

    int importtemperaturevalue(rasterimage::Uint8 value)
    {
        const float factor = static_cast<float>(value) / 255.0f;
        return static_cast<int>(roundf(importedtemperatureminimum + factor * (importedtemperaturemaximum - importedtemperatureminimum)));
    }

    int importprecipitationvalue(rasterimage::Uint8 value)
    {
        const float factor = static_cast<float>(value) / 255.0f;
        return static_cast<int>(roundf(factor * importedprecipitationmaximum));
    }
}

using namespace rasterimports;

bool importlandheightmapcells(const string& filepathname, int width, int height, int sealevel,
    const MapImportSettings* settings, const function<void(int, int, int)>& setnom)
{
    ImportedRaster raster;

    if (!loadworldsizedraster(width, height, filepathname, raster))
        return false;

    GradientStripDecoder decoder;

    if (!loadgradientstripdecoder(settings, decoder))
        return false;


    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            float elevation = 0.0f;

            if (raster.isTiff16)
                elevation = static_cast<float>(rastervalueat(raster, width, i, j));
            else
            {
                const rasterimage::Color colour = raster.image.getPixel(i, j);
                elevation = decoder.enabled ? decodegradientvalue(colour, decoder) : static_cast<float>(colour.r * 10);
            }

            if (elevation > 0.0f)
                setnom(i, j, sealevel + roundimportvalue(elevation));
            else
                setnom(i, j, sealevel - 5000);
        }
    }

    return true;
}

bool importseadepthmapcells(const string& filepathname, int width, int height, int sealevel,
    const MapImportSettings* settings, const function<void(int, int, int)>& setnom)
{
    ImportedRaster raster;

    if (!loadworldsizedraster(width, height, filepathname, raster))
        return false;

    GradientStripDecoder decoder;

    if (!loadgradientstripdecoder(settings, decoder))
        return false;


    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            float depth = 0.0f;

            if (raster.isTiff16)
                depth = static_cast<float>(rastervalueat(raster, width, i, j));
            else
            {
                const rasterimage::Color colour = raster.image.getPixel(i, j);
                depth = decoder.enabled ? decodegradientvalue(colour, decoder) : static_cast<float>(colour.r * 50);
            }

            if (depth > 0.0f)
            {
                int elev = sealevel - roundimportvalue(depth);

                if (elev < 1)
                    elev = 1;

                setnom(i, j, elev);
            }
        }
    }

    return true;
}
