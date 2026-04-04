#include <cmath>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

#include <SFML/Graphics.hpp>

#include "map_imports.hpp"
#include "functions.hpp"

using namespace std;

namespace
{
    constexpr float importedtemperatureminimum = -60.0f;
    constexpr float importedtemperaturemaximum = 60.0f;
    constexpr float importedprecipitationmaximum = 1020.0f;

    struct GradientStripDecoder
    {
        bool enabled = false;
        float minimum = 0.0f;
        float increment = 1.0f;
        std::vector<sf::Color> colours;
        std::unordered_map<unsigned int, int> exactmatches;
    };

    bool loadworldsizedimage(planet& world, const string& filepathname, sf::Image& importimage)
    {
        if (!importimage.loadFromFile(filepathname))
            return false;

        const sf::Vector2u imagesize = importimage.getSize();

        return imagesize.x == world.width() + 1 && imagesize.y == world.height() + 1;
    }

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

    unsigned int packrgb(const sf::Color& colour)
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

        sf::Image gradientimage;

        if (!gradientimage.loadFromFile(settings->gradientStripPath))
            return false;

        const sf::Vector2u gradientsize = gradientimage.getSize();

        if (gradientsize.y != 1 || gradientsize.x == 0)
            return false;

        decoder.enabled = true;
        decoder.minimum = settings->gradientMinimum;
        decoder.increment = settings->gradientIncrement;
        decoder.colours.reserve(gradientsize.x);

        for (unsigned int x = 0; x < gradientsize.x; x++)
        {
            const sf::Color colour = gradientimage.getPixel(x, 0);
            decoder.colours.push_back(colour);
            decoder.exactmatches.emplace(packrgb(colour), static_cast<int>(x));
        }

        return true;
    }

    float decodegradientvalue(const sf::Color& colour, const GradientStripDecoder& decoder)
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
            const sf::Color candidate = decoder.colours[index];
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

    int importtemperaturevalue(sf::Uint8 value)
    {
        const float factor = static_cast<float>(value) / 255.0f;
        return static_cast<int>(roundf(importedtemperatureminimum + factor * (importedtemperaturemaximum - importedtemperatureminimum)));
    }

    int importprecipitationvalue(sf::Uint8 value)
    {
        const float factor = static_cast<float>(value) / 255.0f;
        return static_cast<int>(roundf(factor * importedprecipitationmaximum));
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
    sf::Image importimage;

    if (!loadworldsizedimage(world, filepathname, importimage))
        return false;

    GradientStripDecoder decoder;

    if (!loadgradientstripdecoder(settings, decoder))
        return false;

    const int width = world.width();
    const int height = world.height();
    const int sealevel = world.sealevel();

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            const sf::Color colour = importimage.getPixel(i, j);
            const float elevation = decoder.enabled ? decodegradientvalue(colour, decoder) : static_cast<float>(colour.r * 10);

            if (elevation > 0.0f)
                world.setnom(i, j, sealevel + roundimportvalue(elevation));
            else
                world.setnom(i, j, sealevel - 5000);
        }
    }

    return true;
}

bool importseadepthmap(planet& world, const string& filepathname, const MapImportSettings* settings)
{
    sf::Image importimage;

    if (!loadworldsizedimage(world, filepathname, importimage))
        return false;

    GradientStripDecoder decoder;

    if (!loadgradientstripdecoder(settings, decoder))
        return false;

    const int width = world.width();
    const int height = world.height();
    const int sealevel = world.sealevel();

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            const sf::Color colour = importimage.getPixel(i, j);
            const float depth = decoder.enabled ? decodegradientvalue(colour, decoder) : static_cast<float>(colour.r * 50);

            if (depth > 0.0f)
            {
                int elev = sealevel - roundimportvalue(depth);

                if (elev < 1)
                    elev = 1;

                world.setnom(i, j, elev);
            }
        }
    }

    return true;
}

bool importmountainmap(planet& world, const string& filepathname, vector<vector<bool>>& okmountains, const MapImportSettings* settings)
{
    sf::Image importimage;

    if (!loadworldsizedimage(world, filepathname, importimage))
        return false;

    GradientStripDecoder decoder;

    if (!loadgradientstripdecoder(settings, decoder))
        return false;

    const int width = world.width();
    const int height = world.height();
    vector<vector<int>> rawmountains(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            const sf::Color colour = importimage.getPixel(i, j);
            const float mountainheight = decoder.enabled ? decodegradientvalue(colour, decoder) : static_cast<float>(colour.r * 65);
            rawmountains[i][j] = mountainheight > 0.0f ? roundimportvalue(mountainheight) : 0;
        }
    }

    createmountainsfromraw(world, rawmountains, okmountains);
    return true;
}

bool importvolcanomap(planet& world, const string& filepathname, const MapImportSettings* settings)
{
    sf::Image importimage;

    if (!loadworldsizedimage(world, filepathname, importimage))
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
            const sf::Color colour = importimage.getPixel(i, j);
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
    sf::Image importimage;

    if (!loadworldsizedimage(world, filepathname, importimage))
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
            const sf::Color colour = importimage.getPixel(i, j);
            const float temperature = decoder.enabled ? decodegradientvalue(colour, decoder) : static_cast<float>(importtemperaturevalue(colour.r));
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
    sf::Image importimage;

    if (!loadworldsizedimage(world, filepathname, importimage))
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
            const sf::Color colour = importimage.getPixel(i, j);
            const float precipitation = decoder.enabled ? decodegradientvalue(colour, decoder) : static_cast<float>(importprecipitationvalue(colour.r));
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
