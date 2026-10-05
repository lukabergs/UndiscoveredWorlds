#include <cmath>
#include <cctype>
#include <filesystem>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

#include "io/raster_image.hpp"
#include "io/raster_imports.hpp"

#include "io/map_imports.hpp"
#include "world/functions.hpp"
#include "io/topography_tiff.hpp"

using namespace std;

namespace
{
    using namespace rasterimports;

    void ensureimportedclimatemapsize(planet& world, ImportedClimateMaps& maps)
    {
        const int width = world.width();
        const int height = world.height();
        const int cellcount = (width + 1) * (height + 1);

        if (maps.width != width || maps.height != height)
        {
            maps.width = width;
            maps.height = height;
            maps.hasTemperature = false;
            maps.hasPrecipitation = false;
            maps.annualTemperature.assign(cellcount, 0);
            maps.annualPrecipitation.assign(cellcount, 0);
            return;
        }

        if (static_cast<int>(maps.annualTemperature.size()) != cellcount)
            maps.annualTemperature.assign(cellcount, 0);

        if (static_cast<int>(maps.annualPrecipitation.size()) != cellcount)
            maps.annualPrecipitation.assign(cellcount, 0);
    }
}

void clearimportedclimatemaps(ImportedClimateMaps& maps)
{
    maps.width = 0;
    maps.height = 0;
    maps.hasTemperature = false;
    maps.hasPrecipitation = false;
    maps.annualTemperature.clear();
    maps.annualPrecipitation.clear();
}

bool importlandheightmap(planet& world, const string& filepathname, const MapImportSettings* settings)
{
    return importlandheightmapcells(filepathname, world.width(), world.height(), world.sealevel(), settings,
        [&world](int x, int y, int value) { world.setnom(x, y, value); });
}

bool importseadepthmap(planet& world, const string& filepathname, const MapImportSettings* settings)
{
    return importseadepthmapcells(filepathname, world.width(), world.height(), world.sealevel(), settings,
        [&world](int x, int y, int value) { world.setnom(x, y, value); });
}

bool importmountainmap(planet& world, const string& filepathname, vector<vector<bool>>& okmountains, const MapImportSettings* settings)
{
    ImportedRaster raster;

    if (!loadworldsizedraster(world.width(), world.height(), filepathname, raster))
        return false;

    GradientStripDecoder decoder;

    if (!raster.isTiff16 && !loadgradientstripdecoder(settings, decoder))
        return false;

    const int width = world.width();
    const int height = world.height();
    vector<vector<int>> rawmountains(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            float mountainheight = 0.0f;

            if (raster.isTiff16)
                mountainheight = static_cast<float>(rastervalueat(raster, width, i, j));
            else
            {
                const rasterimage::Color colour = raster.image.getPixel(i, j);
                mountainheight = decoder.enabled ? decodegradientvalue(colour, decoder) : static_cast<float>(colour.r * 65);
            }

            rawmountains[i][j] = mountainheight > 0.0f ? roundimportvalue(mountainheight) : 0;
        }
    }

    createmountainsfromraw(world, rawmountains, okmountains);
    return true;
}

bool importvolcanomap(planet& world, const string& filepathname, const MapImportSettings* settings)
{
    if (hastiffextension(filepathname))
        return false;

    rasterimage::Image importimage;

    if (!loadworldsizedimage(world.width(), world.height(), filepathname, importimage))
        return false;

    GradientStripDecoder decoder;

    if (!loadgradientstripdecoder(settings, decoder))
        return false;

    const int width = world.width();
    const int height = world.height();

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            const rasterimage::Color colour = importimage.getPixel(i, j);
            const float magnitude = decoder.enabled ? decodegradientvalue(colour, decoder) : static_cast<float>(colour.r * 45);

            if (magnitude > 0.0f)
            {
                int elev = roundimportvalue(magnitude);
                bool strato = colour.g > 0;

                if (!strato)
                    elev = elev / 2;

                if (colour.b == 0)
                    elev = 0 - elev;

                world.setvolcano(i, j, elev);
                world.setstrato(i, j, strato);
            }
        }
    }

    return true;
}

bool importtemperaturemap(planet& world, const string& filepathname, ImportedClimateMaps& maps, const MapImportSettings* settings)
{
    ImportedRaster raster;

    if (!loadworldsizedraster(world.width(), world.height(), filepathname, raster))
        return false;

    GradientStripDecoder decoder;

    if (!raster.isTiff16 && !loadgradientstripdecoder(settings, decoder))
        return false;

    ensureimportedclimatemapsize(world, maps);

    const int width = world.width();
    const int height = world.height();

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            float temperature = 0.0f;

            if (raster.isTiff16)
                temperature = static_cast<float>(rastervalueat(raster, width, i, j));
            else
            {
                const rasterimage::Color colour = raster.image.getPixel(i, j);
                temperature = decoder.enabled ? decodegradientvalue(colour, decoder) : static_cast<float>(importtemperaturevalue(colour.r));
            }

            const int index = j * (width + 1) + i;

            maps.annualTemperature[index] = temperature;

            const int previewtemperature = roundimportvalue(temperature);
            world.setjantemp(i, j, previewtemperature);
            world.setjultemp(i, j, previewtemperature);

            for (int season = 0; season < CLIMATESEASONCOUNT; season++)
                world.setseasonaltemp(season, i, j, previewtemperature);
        }
    }

    maps.hasTemperature = true;
    return true;
}

bool importprecipitationmap(planet& world, const string& filepathname, ImportedClimateMaps& maps, const MapImportSettings* settings)
{
    ImportedRaster raster;

    if (!loadworldsizedraster(world.width(), world.height(), filepathname, raster))
        return false;

    GradientStripDecoder decoder;

    if (!loadgradientstripdecoder(settings, decoder))
        return false;

    ensureimportedclimatemapsize(world, maps);

    const int width = world.width();
    const int height = world.height();

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            float precipitation = 0.0f;

            if (raster.isTiff16)
                precipitation = static_cast<float>(rastervalueat(raster, width, i, j));
            else
            {
                const rasterimage::Color colour = raster.image.getPixel(i, j);
                precipitation = decoder.enabled ? decodegradientvalue(colour, decoder) : static_cast<float>(importprecipitationvalue(colour.r));
            }

            const int index = j * (width + 1) + i;

            maps.annualPrecipitation[index] = precipitation;

            const int previewprecipitation = std::max(0, roundimportvalue(precipitation));
            world.setjanrain(i, j, previewprecipitation);
            world.setjulrain(i, j, previewprecipitation);
            world.setseasonalrain(seasonjanuary, i, j, previewprecipitation);
            world.setseasonalrain(seasonjuly, i, j, previewprecipitation);
            world.setseasonalrain(seasonapril, i, j, 0);
            world.setseasonalrain(seasonoctober, i, j, 0);
        }
    }

    maps.hasPrecipitation = true;
    return true;
}
