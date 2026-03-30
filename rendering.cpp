#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>

#include "imgui.h"
#include "cuda_renderers.hpp"
#include "functions.hpp"

using namespace std;

namespace
{
    constexpr int rgbachannels = 4;
    bool timedreportingenabled = false;
    bool hastimedreport = false;
    string timedreporttext;
    chrono::steady_clock::time_point timedreportstart;

    template <typename T>
    void fillgrid(vector<vector<T>>& grid, const T& value)
    {
        for (vector<T>& row : grid)
            fill(row.begin(), row.end(), value);
    }

    size_t pixeloffset(int width, int x, int y)
    {
        return (static_cast<size_t>(y) * static_cast<size_t>(width) + static_cast<size_t>(x)) * rgbachannels;
    }

    void resizepixelbuffer(vector<sf::Uint8>& pixels, int width, int height)
    {
        pixels.resize(static_cast<size_t>(width) * static_cast<size_t>(height) * rgbachannels);
    }

    sf::Uint8 clampchannel(int value)
    {
        return static_cast<sf::Uint8>(clamp(value, 0, 255));
    }

    float clampunit(float value)
    {
        return clamp(value, 0.0f, 1.0f);
    }

    sf::Color makecolour(int r, int g, int b)
    {
        return sf::Color(clampchannel(r), clampchannel(g), clampchannel(b));
    }

    struct RuntimeGradientStop
    {
        int position = 0;
        sf::Color colour = sf::Color::Black;
    };

    struct MapGradientRuntime
    {
        int stopcount = 0;
        bool discrete = false;
        std::array<RuntimeGradientStop, MAPGRADIENTMAXSTOPS> stops{};
    };

    sf::Color lerpcolour(const sf::Color& low, const sf::Color& high, float factor)
    {
        factor = clampunit(factor);

        return sf::Color(
            clampchannel(static_cast<int>(roundf(static_cast<float>(low.r) + (static_cast<float>(high.r) - static_cast<float>(low.r)) * factor))),
            clampchannel(static_cast<int>(roundf(static_cast<float>(low.g) + (static_cast<float>(high.g) - static_cast<float>(low.g)) * factor))),
            clampchannel(static_cast<int>(roundf(static_cast<float>(low.b) + (static_cast<float>(high.b) - static_cast<float>(low.b)) * factor))));
    }

    sf::Color lerpcolour(const sf::Color& low, const sf::Color& mid, const sf::Color& high, float factor)
    {
        factor = clampunit(factor);

        if (factor <= 0.5f)
            return lerpcolour(low, mid, factor * 2.0f);

        return lerpcolour(mid, high, (factor - 0.5f) * 2.0f);
    }

    sf::Color getoutlinecolour(const planet& world)
    {
        return makecolour(world.outline1(), world.outline2(), world.outline3());
    }

    sf::Color getworldcolour(int r, int g, int b)
    {
        return makecolour(r, g, b);
    }

    MapGradientRuntime buildmapgradient(const planet& world, int gradient)
    {
        MapGradientRuntime result;
        result.stopcount = clamp(world.mapgradientstopcount(gradient), 2, MAPGRADIENTMAXSTOPS);
        result.discrete = world.mapgradientdiscrete(gradient);

        for (int i = 0; i < result.stopcount; i++)
        {
            result.stops[i].position = world.mapgradientposition(gradient, i);
            result.stops[i].colour = getworldcolour(
                world.mapgradientcolour(gradient, i, 0),
                world.mapgradientcolour(gradient, i, 1),
                world.mapgradientcolour(gradient, i, 2));
        }

        sort(result.stops.begin(), result.stops.begin() + result.stopcount, [](const RuntimeGradientStop& left, const RuntimeGradientStop& right)
        {
            return left.position < right.position;
        });

        return result;
    }

    sf::Color samplemapgradient(const MapGradientRuntime& gradient, int value)
    {
        if (gradient.stopcount <= 0)
            return sf::Color::Black;

        if (value <= gradient.stops[0].position)
            return gradient.stops[0].colour;

        for (int i = 0; i < gradient.stopcount - 1; i++)
        {
            const RuntimeGradientStop& low = gradient.stops[i];
            const RuntimeGradientStop& high = gradient.stops[i + 1];

            if (value <= high.position)
            {
                if (gradient.discrete || high.position <= low.position)
                    return low.colour;

                const float factor = static_cast<float>(value - low.position) / static_cast<float>(high.position - low.position);
                return lerpcolour(low.colour, high.colour, factor);
            }
        }

        return gradient.stops[gradient.stopcount - 1].colour;
    }

    sf::Color getelevationmapcolour(const MapGradientRuntime& gradient, int altitude)
    {
        return samplemapgradient(gradient, altitude);
    }

    sf::Color gettemperaturemapcolour(const MapGradientRuntime& gradient, int temperature)
    {
        return samplemapgradient(gradient, temperature);
    }

    sf::Color getprecipitationmapcolour(const MapGradientRuntime& gradient, int rainfall)
    {
        return samplemapgradient(gradient, rainfall);
    }

    sf::Color getclimateseacolour(const planet& world, int seaice)
    {
        int slot = climateopensea;

        if (seaice == 2)
            slot = climatepermanentseaice;
        else if (seaice == 1)
            slot = climateseasonalseaice;

        return getworldcolour(
            world.climatemapseacolour(slot, 0),
            world.climatemapseacolour(slot, 1),
            world.climatemapseacolour(slot, 2));
    }

    sf::Color getclimatecoloursforworld(const planet& world, short climate)
    {
        return getworldcolour(
            world.climatemapcolour(climate, 0),
            world.climatemapcolour(climate, 1),
            world.climatemapcolour(climate, 2));
    }

    sf::Color averagecolour(const sf::Color& left, const sf::Color& right)
    {
        return sf::Color(
            clampchannel((static_cast<int>(left.r) + static_cast<int>(right.r)) / 2),
            clampchannel((static_cast<int>(left.g) + static_cast<int>(right.g)) / 2),
            clampchannel((static_cast<int>(left.b) + static_cast<int>(right.b)) / 2));
    }

    enum class holdridgethermalband
    {
        ice,
        polar,
        subpolar,
        boreal,
        cooltemperate,
        warmtemperate,
        subtropical,
        tropical
    };

    enum class holdridgehumidityband
    {
        superarid,
        perarid,
        arid,
        semiarid,
        subhumid,
        humid,
        perhumid,
        superhumid
    };

    float clampbiotemperature(float temperature)
    {
        return clamp(temperature, 0.0f, 30.0f);
    }

    float calculateholdridgebiotemperature(float jan, float apr, float jul, float oct)
    {
        return (clampbiotemperature(jan) + clampbiotemperature(apr) + clampbiotemperature(jul) + clampbiotemperature(oct)) / 4.0f;
    }

    float calculateholdridgeannualprecipitation(float jan, float apr, float jul, float oct)
    {
        return max(1.0f, (jan + apr + jul + oct) * 3.0f);
    }

    holdridgethermalband classifyholdridgethermalband(float biotemperature)
    {
        if (biotemperature <= 0.0f)
            return holdridgethermalband::ice;
        if (biotemperature < 1.5f)
            return holdridgethermalband::polar;
        if (biotemperature < 3.0f)
            return holdridgethermalband::subpolar;
        if (biotemperature < 6.0f)
            return holdridgethermalband::boreal;
        if (biotemperature < 12.0f)
            return holdridgethermalband::cooltemperate;
        if (biotemperature < 18.0f)
            return holdridgethermalband::warmtemperate;
        if (biotemperature < 24.0f)
            return holdridgethermalband::subtropical;

        return holdridgethermalband::tropical;
    }

    holdridgehumidityband classifyholdridgehumidityband(float petratio)
    {
        if (petratio > 16.0f)
            return holdridgehumidityband::superarid;
        if (petratio > 8.0f)
            return holdridgehumidityband::perarid;
        if (petratio > 4.0f)
            return holdridgehumidityband::arid;
        if (petratio > 2.0f)
            return holdridgehumidityband::semiarid;
        if (petratio > 1.0f)
            return holdridgehumidityband::subhumid;
        if (petratio > 0.5f)
            return holdridgehumidityband::humid;
        if (petratio > 0.25f)
            return holdridgehumidityband::perhumid;

        return holdridgehumidityband::superhumid;
    }

    sf::Color getbiomemapcolour(const planet& world, int slot)
    {
        return makecolour(
            world.biomemapcolour(slot, 0),
            world.biomemapcolour(slot, 1),
            world.biomemapcolour(slot, 2));
    }

    int getholdridgebiomeslot(holdridgethermalband thermalband, holdridgehumidityband humidityband)
    {
        switch (thermalband)
        {
        case holdridgethermalband::ice:
            return biomeice;

        case holdridgethermalband::polar:
            switch (humidityband)
            {
            case holdridgehumidityband::superarid:
                return biomepolardesert;
            case holdridgehumidityband::perarid:
            case holdridgehumidityband::arid:
                return biomepolardrytundra;
            case holdridgehumidityband::semiarid:
                return biomepolarmoisttundra;
            case holdridgehumidityband::subhumid:
            case holdridgehumidityband::humid:
                return biomepolarwettundra;
            case holdridgehumidityband::perhumid:
            case holdridgehumidityband::superhumid:
                return biomepolarraindtundra;
            }
            break;

        case holdridgethermalband::subpolar:
            switch (humidityband)
            {
            case holdridgehumidityband::superarid:
                return biomesubpolardesert;
            case holdridgehumidityband::perarid:
            case holdridgehumidityband::arid:
                return biomesubpolardrytundra;
            case holdridgehumidityband::semiarid:
                return biomesubpolarmoisttundra;
            case holdridgehumidityband::subhumid:
            case holdridgehumidityband::humid:
                return biomesubpolarwettundra;
            case holdridgehumidityband::perhumid:
            case holdridgehumidityband::superhumid:
                return biomesubpolarraindtundra;
            }
            break;

        case holdridgethermalband::boreal:
            switch (humidityband)
            {
            case holdridgehumidityband::superarid:
                return biomeborealdesert;
            case holdridgehumidityband::perarid:
            case holdridgehumidityband::arid:
                return biomeborealdrybush;
            case holdridgehumidityband::semiarid:
            case holdridgehumidityband::subhumid:
                return biomeborealmoistforest;
            case holdridgehumidityband::humid:
            case holdridgehumidityband::perhumid:
                return biomeborealwetforest;
            case holdridgehumidityband::superhumid:
                return biomeborealrainforest;
            }
            break;

        case holdridgethermalband::cooltemperate:
            switch (humidityband)
            {
            case holdridgehumidityband::superarid:
                return biomecooltemperatedesert;
            case holdridgehumidityband::perarid:
                return biomecooltemperatedesertbush;
            case holdridgehumidityband::arid:
            case holdridgehumidityband::semiarid:
                return biomecooltemperatesteppe;
            case holdridgehumidityband::subhumid:
                return biomecooltemperatemoistforest;
            case holdridgehumidityband::humid:
            case holdridgehumidityband::perhumid:
                return biomecooltemperatewetforest;
            case holdridgehumidityband::superhumid:
                return biomecooltemperaterainforest;
            }
            break;

        case holdridgethermalband::warmtemperate:
            switch (humidityband)
            {
            case holdridgehumidityband::superarid:
                return biomewarmtemperatedesert;
            case holdridgehumidityband::perarid:
                return biomewarmtemperatedesertbush;
            case holdridgehumidityband::arid:
                return biomewarmtemperatethornsteppe;
            case holdridgehumidityband::semiarid:
                return biomewarmtemperatedryforest;
            case holdridgehumidityband::subhumid:
                return biomewarmtemperatemoistforest;
            case holdridgehumidityband::humid:
            case holdridgehumidityband::perhumid:
                return biomewarmtemperatewetforest;
            case holdridgehumidityband::superhumid:
                return biomewarmtemperaterainforest;
            }
            break;

        case holdridgethermalband::subtropical:
            switch (humidityband)
            {
            case holdridgehumidityband::superarid:
                return biomesubtropicaldesert;
            case holdridgehumidityband::perarid:
                return biomesubtropicaldesertbush;
            case holdridgehumidityband::arid:
                return biomesubtropicalthornsteppe;
            case holdridgehumidityband::semiarid:
                return biomesubtropicaldryforest;
            case holdridgehumidityband::subhumid:
                return biomesubtropicalmoistforest;
            case holdridgehumidityband::humid:
            case holdridgehumidityband::perhumid:
                return biomesubtropicalwetforest;
            case holdridgehumidityband::superhumid:
                return biomesubtropicalrainforest;
            }
            break;

        case holdridgethermalband::tropical:
            switch (humidityband)
            {
            case holdridgehumidityband::superarid:
                return biometropicaldesert;
            case holdridgehumidityband::perarid:
                return biometropicaldesertbush;
            case holdridgehumidityband::arid:
                return biometropicalthornsteppe;
            case holdridgehumidityband::semiarid:
                return biometropicalverydryforest;
            case holdridgehumidityband::subhumid:
                return biometropicaldryforest;
            case holdridgehumidityband::humid:
                return biometropicalmoistforest;
            case holdridgehumidityband::perhumid:
                return biometropicalwetforest;
            case holdridgehumidityband::superhumid:
                return biometropicalrainforest;
            }
            break;
        }

        return biomeice;
    }

    int getholdridgebiomeslot(float janTemp, float aprTemp, float julTemp, float octTemp, float janRain, float aprRain, float julRain, float octRain)
    {
        const float biotemperature = calculateholdridgebiotemperature(janTemp, aprTemp, julTemp, octTemp);
        const float annualprecipitation = calculateholdridgeannualprecipitation(janRain, aprRain, julRain, octRain);
        const float petratio = biotemperature <= 0.0f ? 1000.0f : (biotemperature * 58.93f) / annualprecipitation;

        return getholdridgebiomeslot(classifyholdridgethermalband(biotemperature), classifyholdridgehumidityband(petratio));
    }

    void applyglobaloutlineoverlay(planet& world, sf::Image& image)
    {
        if (world.showmapoutline() == false)
            return;

        const sf::Color outlinecolour = getoutlinecolour(world);

        for (int i = 0; i <= world.width(); i++)
        {
            for (int j = 0; j <= world.height(); j++)
            {
                if (world.outline(i, j))
                    image.setPixel(i, j, outlinecolour);
            }
        }
    }

    void applyglobaldisplayoutlineoverlay(planet& world, sf::Image& displayimage)
    {
        if (world.showmapoutline() == false)
            return;

        const sf::Color outlinecolour = getoutlinecolour(world);
        const int displaywidth = static_cast<int>(displayimage.getSize().x);
        const int displayheight = static_cast<int>(displayimage.getSize().y);
        const float mapdiv = (static_cast<float>(world.width()) + 1.0f) / static_cast<float>(displaywidth);

        for (int i = 0; i < displaywidth; i++)
        {
            const int sourcei = static_cast<int>(static_cast<float>(i) * mapdiv);

            for (int j = 0; j < displayheight; j++)
            {
                const int sourcej = static_cast<int>(static_cast<float>(j) * mapdiv);
                bool outlined = false;

                for (int k = max(0, sourcei - 1); k <= min(world.width(), sourcei); k++)
                {
                    for (int l = max(0, sourcej - 1); l <= min(world.height(), sourcej); l++)
                    {
                        if (world.outline(k, l))
                        {
                            outlined = true;
                            break;
                        }
                    }

                    if (outlined)
                        break;
                }

                if (outlined)
                    displayimage.setPixel(i, j, outlinecolour);
            }
        }
    }

    void applyregionaloutlineoverlay(planet& world, region& region, sf::Image& image)
    {
        if (world.showmapoutline() == false)
            return;

        const sf::Color outlinecolour = getoutlinecolour(world);

        for (int i = region.regwidthbegin(); i <= region.regwidthend(); i++)
        {
            for (int j = region.regheightbegin(); j <= region.regheightend(); j++)
            {
                if (region.outline(i, j))
                    image.setPixel(i - region.regwidthbegin(), j - region.regheightbegin(), outlinecolour);
            }
        }
    }

    sf::Color getriverflowcolour(const MapGradientRuntime& gradient, int flow)
    {
        return samplemapgradient(gradient, flow);
    }

    void flushtimedreport()
    {
        if (hastimedreport == false)
            return;

        const double elapsedms = chrono::duration<double, milli>(chrono::steady_clock::now() - timedreportstart).count();
        const streamsize oldprecision = cout.precision();
        const ios::fmtflags oldflags = cout.flags();

        cout << '\r' << left << setw(48) << timedreporttext << " [" << fixed << setprecision(2) << elapsedms << "ms]" << endl;
        cout.flags(oldflags);
        cout.precision(oldprecision);

        onworldgenstepcompleted(timedreporttext, elapsedms);

        hastimedreport = false;
        timedreporttext.clear();
    }

    void setpixel(vector<sf::Uint8>& pixels, int width, int x, int y, sf::Uint8 r, sf::Uint8 g, sf::Uint8 b, sf::Uint8 a = 255)
    {
        const size_t offset = pixeloffset(width, x, y);
        pixels[offset] = r;
        pixels[offset + 1] = g;
        pixels[offset + 2] = b;
        pixels[offset + 3] = a;
    }

    const sf::Uint8* getpixel(const sf::Uint8* pixels, int width, int x, int y)
    {
        return pixels + pixeloffset(width, x, y);
    }

    bool isblackpixel(const sf::Uint8* pixel)
    {
        return pixel[0] == 0 && pixel[1] == 0 && pixel[2] == 0;
    }

    void createnearestdisplayimage(const sf::Uint8* sourcepixels, int sourcewidth, int worldwidth, sf::Image& displayimage)
    {
        static vector<sf::Uint8> displaypixels;

        const int displaywidth = static_cast<int>(displayimage.getSize().x);
        const int displayheight = static_cast<int>(displayimage.getSize().y);

        resizepixelbuffer(displaypixels, displaywidth, displayheight);

        const float mapdiv = (static_cast<float>(worldwidth) + 1.0f) / static_cast<float>(displaywidth);

        for (int i = 0; i < displaywidth; i++)
        {
            const int sourcei = static_cast<int>(static_cast<float>(i) * mapdiv);

            for (int j = 0; j < displayheight; j++)
            {
                const int sourcej = static_cast<int>(static_cast<float>(j) * mapdiv);
                const sf::Uint8* sourcepixel = getpixel(sourcepixels, sourcewidth, sourcei, sourcej);
                setpixel(displaypixels, displaywidth, i, j, sourcepixel[0], sourcepixel[1], sourcepixel[2], sourcepixel[3]);
            }
        }

        displayimage.create(displaywidth, displayheight, displaypixels.data());
    }

    void createoutlineddisplayimage(const sf::Uint8* sourcepixels, int sourcewidth, int worldwidth, sf::Image& displayimage)
    {
        static vector<sf::Uint8> displaypixels;

        const int displaywidth = static_cast<int>(displayimage.getSize().x);
        const int displayheight = static_cast<int>(displayimage.getSize().y);

        resizepixelbuffer(displaypixels, displaywidth, displayheight);

        const float mapdiv = (static_cast<float>(worldwidth) + 1.0f) / static_cast<float>(displaywidth);
        const bool shouldpreserveoutline = worldwidth > displaywidth;

        for (int i = 0; i < displaywidth; i++)
        {
            const int sourcei = static_cast<int>(static_cast<float>(i) * mapdiv);

            for (int j = 0; j < displayheight; j++)
            {
                const int sourcej = static_cast<int>(static_cast<float>(j) * mapdiv);
                const sf::Uint8* sourcepixel = getpixel(sourcepixels, sourcewidth, sourcei, sourcej);
                const size_t offset = pixeloffset(displaywidth, i, j);

                displaypixels[offset] = sourcepixel[0];
                displaypixels[offset + 1] = sourcepixel[1];
                displaypixels[offset + 2] = sourcepixel[2];
                displaypixels[offset + 3] = sourcepixel[3];

                if (shouldpreserveoutline && sourcei > 0 && sourcej > 0)
                {
                    for (int k = sourcei - 1; k <= sourcei; k++)
                    {
                        for (int l = sourcej - 1; l <= sourcej; l++)
                        {
                            const sf::Uint8* outlinepixel = getpixel(sourcepixels, sourcewidth, k, l);

                            if (isblackpixel(outlinepixel))
                            {
                                displaypixels[offset] = outlinepixel[0];
                                displaypixels[offset + 1] = outlinepixel[1];
                                displaypixels[offset + 2] = outlinepixel[2];
                                displaypixels[offset + 3] = outlinepixel[3];
                            }
                        }
                    }
                }
            }
        }

        displayimage.create(displaywidth, displayheight, displaypixels.data());
    }

    void createriversdisplayimage(const sf::Uint8* sourcepixels, int sourcewidth, int worldwidth, int worldsize, sf::Image& displayimage)
    {
        static vector<sf::Uint8> displaypixels;

        const int displaywidth = static_cast<int>(displayimage.getSize().x);
        const int displayheight = static_cast<int>(displayimage.getSize().y);

        resizepixelbuffer(displaypixels, displaywidth, displayheight);

        if (worldsize == 0)
        {
            for (int i = 0; i < displaywidth; i++)
            {
                for (int j = 0; j < displayheight; j++)
                {
                    const sf::Uint8* sourcepixel = getpixel(sourcepixels, sourcewidth, i / 2, j / 2);
                    setpixel(displaypixels, displaywidth, i, j, sourcepixel[0], sourcepixel[1], sourcepixel[2], sourcepixel[3]);
                }
            }

            displayimage.create(displaywidth, displayheight, displaypixels.data());
            return;
        }

        if (worldsize == 1)
        {
            for (int i = 0; i < displaywidth; i++)
            {
                for (int j = 0; j < displayheight; j++)
                {
                    const sf::Uint8* sourcepixel = getpixel(sourcepixels, sourcewidth, i, j);
                    setpixel(displaypixels, displaywidth, i, j, sourcepixel[0], sourcepixel[1], sourcepixel[2], sourcepixel[3]);
                }
            }

            displayimage.create(displaywidth, displayheight, displaypixels.data());
            return;
        }

        const float mapdiv = (static_cast<float>(worldwidth) + 1.0f) / static_cast<float>(displaywidth);
        const bool shouldpreserveoutline = worldwidth > displaywidth;

        for (int i = 0; i < displaywidth; i++)
        {
            const int sourcei = static_cast<int>(static_cast<float>(i) * mapdiv);

            for (int j = 0; j < displayheight; j++)
            {
                const int sourcej = static_cast<int>(static_cast<float>(j) * mapdiv);
                const size_t offset = pixeloffset(displaywidth, i, j);
                const sf::Uint8* sourcepixel = getpixel(sourcepixels, sourcewidth, sourcei, sourcej);

                displaypixels[offset] = sourcepixel[0];
                displaypixels[offset + 1] = sourcepixel[1];
                displaypixels[offset + 2] = sourcepixel[2];
                displaypixels[offset + 3] = sourcepixel[3];

                if (sourcei > 0 && sourcej > 0)
                {
                    int lowestred = 255;
                    int lowestgreen = 255;
                    int lowestblue = 255;

                    for (int k = sourcei - 1; k <= sourcei; k++)
                    {
                        for (int l = sourcej - 1; l <= sourcej; l++)
                        {
                            const sf::Uint8* nearbypixel = getpixel(sourcepixels, sourcewidth, k, l);

                            if (nearbypixel[0] < lowestred && nearbypixel[0] == nearbypixel[1])
                            {
                                lowestred = nearbypixel[0];
                                lowestgreen = nearbypixel[1];
                                lowestblue = nearbypixel[2];
                            }
                        }
                    }

                    if (lowestred < 255)
                    {
                        displaypixels[offset] = static_cast<sf::Uint8>(lowestred);
                        displaypixels[offset + 1] = static_cast<sf::Uint8>(lowestgreen);
                        displaypixels[offset + 2] = static_cast<sf::Uint8>(lowestblue);
                    }

                    if (shouldpreserveoutline)
                    {
                        for (int k = sourcei - 1; k <= sourcei; k++)
                        {
                            for (int l = sourcej - 1; l <= sourcej; l++)
                            {
                                const sf::Uint8* outlinepixel = getpixel(sourcepixels, sourcewidth, k, l);

                                if (isblackpixel(outlinepixel))
                                {
                                    displaypixels[offset] = outlinepixel[0];
                                    displaypixels[offset + 1] = outlinepixel[1];
                                    displaypixels[offset + 2] = outlinepixel[2];
                                    displaypixels[offset + 3] = outlinepixel[3];
                                }
                            }
                        }
                    }
                }
            }
        }

        displayimage.create(displaywidth, displayheight, displaypixels.data());
    }

    vector<sf::Color> buildtundracolours(planet& world)
    {
        vector<sf::Color> colours(world.height() + 1);

        for (int row = 0; row <= world.height(); row++)
        {
            int lat = 0;
            int latminutes = 0;
            int latseconds = 0;
            bool latneg = false;

            world.latitude(row, lat, latminutes, latseconds, latneg);

            const int lat2 = 90 - lat;

            colours[row] = sf::Color(
                (lat * world.tundra1() + lat2 * world.eqtundra1()) / 90,
                (lat * world.tundra2() + lat2 * world.eqtundra2()) / 90,
                (lat * world.tundra3() + lat2 * world.eqtundra3()) / 90);
        }

        return colours;
    }

    bool isglobalnonreliefmapview(mapviewenum mapview)
    {
        return mapview != relief;
    }

    RgbTriplet makergb(int r, int g, int b)
    {
        RgbTriplet colour;
        colour.r = r;
        colour.g = g;
        colour.b = b;
        return colour;
    }

    void createimagefromrgba(const std::vector<std::uint8_t>& pixels, int width, int height, sf::Image& image)
    {
        if (pixels.empty())
            return;

        image.create(width, height, reinterpret_cast<const sf::Uint8*>(pixels.data()));
    }

    GlobalCudaRendererInputs makeglobalcudarendererinputs(planet& world)
    {
        GlobalCudaRendererInputs inputs;

        inputs.width = world.width();
        inputs.height = world.height();
        inputs.sealevel = world.sealevel();
        inputs.maxelevation = world.maxelevation();
        inputs.maxriverflow = world.maxriverflow();
        inputs.worldsize = world.size();
        inputs.perihelion = world.perihelion();
        inputs.eccentricity = world.eccentricity();
        inputs.tilt = world.tilt();

        inputs.mapnom = world.rawmapnom();
        inputs.lakemap = world.rawlakemap();
        inputs.oceanridgeheightmap = world.rawoceanridgeheightmap();
        inputs.mountainheights = world.rawmountainheights();
        inputs.volcanomap = world.rawvolcanomap();
        inputs.extraelevmap = world.rawextraelevmap();
        inputs.craterrims = world.rawcraterrims();
        inputs.cratercentres = world.rawcratercentres();
        inputs.oceanridgemap = world.rawoceanridgemap();
        inputs.oceanriftmap = world.rawoceanriftmap();
        inputs.jantempmap = world.rawjantempmap();
        inputs.jultempmap = world.rawjultempmap();
        inputs.janrainmap = world.rawjanrainmap();
        inputs.julrainmap = world.rawjulrainmap();
        inputs.seaicemap = world.rawseaicemap();
        inputs.climatemap = world.rawclimatemap();
        inputs.rivermapjan = world.rawrivermapjan();
        inputs.rivermapjul = world.rawrivermapjul();
        inputs.specials = world.rawspecials();
        inputs.deltamapdir = world.rawdeltamapdir();
        inputs.deltamapjan = world.rawdeltamapjan();
        inputs.deltamapjul = world.rawdeltamapjul();
        inputs.riftlakemapsurface = world.rawriftlakemapsurface();
        inputs.noshademap = world.rawnoshademap();

        inputs.relief.seed = static_cast<std::uint64_t>(world.seed());
        inputs.relief.minriverflow = world.minriverflowglobal();
        inputs.relief.shadingdir = world.shadingdir();
        inputs.relief.seaiceappearance = world.seaiceappearance();
        inputs.relief.snowchange = world.snowchange();
        inputs.relief.colourcliffs = world.colourcliffs();
        inputs.relief.landshading = world.landshading();
        inputs.relief.lakeshading = world.lakeshading();
        inputs.relief.seashading = world.seashading();
        inputs.relief.landmarbling = world.landmarbling() * 2.0f;
        inputs.relief.lakemarbling = world.lakemarbling() * 2.0f;
        inputs.relief.seamarbling = world.seamarbling() * 2.0f;
        inputs.relief.showmangroves = world.showmangroves();
        inputs.relief.seaice = makergb(world.seaice1(), world.seaice2(), world.seaice3());
        inputs.relief.ocean = makergb(world.ocean1(), world.ocean2(), world.ocean3());
        inputs.relief.deepocean = makergb(world.deepocean1(), world.deepocean2(), world.deepocean3());
        inputs.relief.river = makergb(world.river1(), world.river2(), world.river3());
        inputs.relief.saltpan = makergb(world.saltpan1(), world.saltpan2(), world.saltpan3());
        inputs.relief.base = makergb(world.base1(), world.base2(), world.base3());
        inputs.relief.basetemp = makergb(world.basetemp1(), world.basetemp2(), world.basetemp3());
        inputs.relief.highbase = makergb(world.highbase1(), world.highbase2(), world.highbase3());
        inputs.relief.desert = makergb(world.desert1(), world.desert2(), world.desert3());
        inputs.relief.highdesert = makergb(world.highdesert1(), world.highdesert2(), world.highdesert3());
        inputs.relief.colddesert = makergb(world.colddesert1(), world.colddesert2(), world.colddesert3());
        inputs.relief.grass = makergb(world.grass1(), world.grass2(), world.grass3());
        inputs.relief.cold = makergb(world.cold1(), world.cold2(), world.cold3());
        inputs.relief.eqtundra = makergb(world.eqtundra1(), world.eqtundra2(), world.eqtundra3());
        inputs.relief.tundra = makergb(world.tundra1(), world.tundra2(), world.tundra3());
        inputs.relief.erg = makergb(world.erg1(), world.erg2(), world.erg3());
        inputs.relief.wetlands = makergb(world.wetlands1(), world.wetlands2(), world.wetlands3());
        inputs.relief.lake = makergb(world.lake1(), world.lake2(), world.lake3());
        inputs.relief.sand = makergb(world.sand1(), world.sand2(), world.sand3());
        inputs.relief.shingle = makergb(world.shingle1(), world.shingle2(), world.shingle3());
        inputs.relief.mud = makergb(world.mud1(), world.mud2(), world.mud3());
        inputs.relief.mangrove = makergb(world.mangrove1(), world.mangrove2(), world.mangrove3());
        inputs.relief.glacier = makergb(world.glacier1(), world.glacier2(), world.glacier3());

        return inputs;
    }

    RegionalCudaRendererInputs makeregionalcudarendererinputs(planet& world, region& region)
    {
        RegionalCudaRendererInputs inputs;

        inputs.sealevel = region.sealevel();
        inputs.maxelevation = world.maxelevation();
        inputs.maxriverflow = world.maxriverflow();
        inputs.regwidthbegin = region.regwidthbegin();
        inputs.regwidthend = region.regwidthend();
        inputs.regheightbegin = region.regheightbegin();
        inputs.regheightend = region.regheightend();
        inputs.regionwidth = region.rwidth();
        inputs.regionheight = region.rheight();
        inputs.lefty = region.lefty();
        inputs.globalheight = world.height();
        inputs.map = region.rawmap();
        inputs.lakemap = region.rawlakemap();
        inputs.jantempmap = region.rawjantempmap();
        inputs.jultempmap = region.rawjultempmap();
        inputs.extrajantempmap = region.rawextrajantempmap();
        inputs.extrajultempmap = region.rawextrajultempmap();
        inputs.janrainmap = region.rawjanrainmap();
        inputs.julrainmap = region.rawjulrainmap();
        inputs.climatemap = region.rawclimatemap();
        inputs.seaicemap = region.rawseaicemap();
        inputs.rivermapdir = region.rawrivermapdir();
        inputs.rivermapjan = region.rawrivermapjan();
        inputs.rivermapjul = region.rawrivermapjul();
        inputs.fakeriversdir = region.rawfakeriversdir();
        inputs.fakeriversjan = region.rawfakeriversjan();
        inputs.fakeriversjul = region.rawfakeriversjul();
        inputs.specials = region.rawspecials();
        inputs.deltamapdir = region.rawdeltamapdir();
        inputs.deltamapjan = region.rawdeltamapjan();
        inputs.deltamapjul = region.rawdeltamapjul();
        inputs.volcanomap = region.rawvolcanomap();
        inputs.mudmap = region.rawmudmap();
        inputs.sandmap = region.rawsandmap();
        inputs.shinglemap = region.rawshinglemap();
        inputs.testmap = region.rawtestmap();
        inputs.relief = makeglobalcudarendererinputs(world).relief;
        inputs.relief.minriverflow = world.minriverflowregional();

        return inputs;
    }

    void assignglobalcudamaplayer(maplayer& layer, const std::vector<std::uint8_t>& fullpixels, const std::vector<std::uint8_t>& displaypixels, int imagewidth, int imageheight)
    {
        createimagefromrgba(fullpixels, imagewidth, imageheight, layer.image);
        createimagefromrgba(displaypixels, DISPLAYMAPSIZEX, DISPLAYMAPSIZEY, layer.displayimage);
        layer.created = true;
    }

    bool drawglobalnonreliefmapimagescuda(planet& world, mapcache& maps)
    {
        GlobalCudaRendererOutputs outputs;

        if (!renderglobalnonreliefmapscuda(makeglobalcudarendererinputs(world), outputs))
            return false;

        const int imagewidth = world.width() + 1;
        const int imageheight = world.height() + 2;

        assignglobalcudamaplayer(getmaplayer(maps, elevation), outputs.elevation, outputs.elevationDisplay, imagewidth, imageheight);
        assignglobalcudamaplayer(getmaplayer(maps, temperature), outputs.temperature, outputs.temperatureDisplay, imagewidth, imageheight);
        assignglobalcudamaplayer(getmaplayer(maps, precipitation), outputs.precipitation, outputs.precipitationDisplay, imagewidth, imageheight);
        assignglobalcudamaplayer(getmaplayer(maps, climate), outputs.climate, outputs.climateDisplay, imagewidth, imageheight);
        assignglobalcudamaplayer(getmaplayer(maps, rivers), outputs.rivers, outputs.riversDisplay, imagewidth, imageheight);
        return true;
    }

    bool drawglobalreliefmapimagecuda(planet& world, maplayer& layer)
    {
        std::vector<std::uint8_t> reliefpixels;
        std::vector<std::uint8_t> displaypixels;

        if (!renderglobalreliefmapcuda(makeglobalcudarendererinputs(world), reliefpixels, displaypixels))
            return false;

        assignglobalcudamaplayer(layer, reliefpixels, displaypixels, world.width() + 1, world.height() + 2);
        return true;
    }

    void assignregionalcudamaplayer(maplayer& layer, const std::vector<std::uint8_t>& pixels, int imagewidth, int imageheight)
    {
        createimagefromrgba(pixels, imagewidth, imageheight, layer.image);
        layer.created = true;
    }

    bool drawregionalnonreliefmapimagescuda(planet& world, region& region, mapcache& maps)
    {
        RegionalCudaRendererOutputs outputs;

        if (!renderregionalnonreliefmapscuda(makeregionalcudarendererinputs(world, region), outputs))
            return false;

        const int imagewidth = region.regwidthend() - region.regwidthbegin() + 1;
        const int imageheight = region.regheightend() - region.regheightbegin() + 1;

        assignregionalcudamaplayer(getmaplayer(maps, elevation), outputs.elevation, imagewidth, imageheight);
        assignregionalcudamaplayer(getmaplayer(maps, temperature), outputs.temperature, imagewidth, imageheight);
        assignregionalcudamaplayer(getmaplayer(maps, precipitation), outputs.precipitation, imagewidth, imageheight);
        assignregionalcudamaplayer(getmaplayer(maps, climate), outputs.climate, imagewidth, imageheight);
        assignregionalcudamaplayer(getmaplayer(maps, rivers), outputs.rivers, imagewidth, imageheight);
        return true;
    }

    bool drawregionalreliefmapimagecuda(planet& world, region& region, maplayer& layer)
    {
        std::vector<std::uint8_t> reliefpixels;

        if (!renderregionalreliefmapcuda(makeregionalcudarendererinputs(world, region), reliefpixels))
            return false;

        const int imagewidth = region.regwidthend() - region.regwidthbegin() + 1;
        const int imageheight = region.regheightend() - region.regheightbegin() + 1;
        assignregionalcudamaplayer(layer, reliefpixels, imagewidth, imageheight);
        return true;
    }
}

void resetmapcache(mapcache& maps)
{
    for (maplayer& layer : maps.layers)
        layer.created = false;
}

maplayer& getmaplayer(mapcache& maps, mapviewenum mapview)
{
    return maps.layers[mapviewindex(mapview)];
}

const maplayer& getmaplayer(const mapcache& maps, mapviewenum mapview)
{
    return maps.layers[mapviewindex(mapview)];
}

sf::Image& getmapimage(mapcache& maps, mapviewenum mapview)
{
    return getmaplayer(maps, mapview).image;
}

const sf::Image& getmapimage(const mapcache& maps, mapviewenum mapview)
{
    return getmaplayer(maps, mapview).image;
}

sf::Image& getdisplaymapimage(mapcache& maps, mapviewenum mapview)
{
    return getmaplayer(maps, mapview).displayimage;
}

const sf::Image& getdisplaymapimage(const mapcache& maps, mapviewenum mapview)
{
    return getmaplayer(maps, mapview).displayimage;
}

void updateTextureFromImage(sf::Texture& texture, const sf::Image& image)
{
    if (texture.getSize() != image.getSize())
        texture.loadFromImage(image);
    else
        texture.update(image);
}


// This makes adjustments to accommodate different world sizes. (Basically resize the map images.)

void adjustforsize(planet& world, sf::Vector2i& globaltexturesize, mapcache& globalmaps, sf::Image& highlightimage, int highlightsize, sf::Image& minihighlightimage, int& minihighlightsize)
{
    int size = world.size();
    int width = 0;
    int height = 0;

    if (size == 0) // Small
    {
        width = 511;
        height = 256;
    }

    if (size == 1) // Medium
    {
        width = 1023;
        height = 512;
    }

    if (size == 2) // Large
    {
        width = 2047;
        height = 1024;
    }

    world.setwidth(width);
    world.setheight(height);

    globaltexturesize.x = world.width() + 1;
    globaltexturesize.y = world.height() + 2;

    for (mapviewenum mapview : allmapviews)
        getmapimage(globalmaps, mapview).create(globaltexturesize.x, globaltexturesize.y, sf::Color::Black);

    resetmapcache(globalmaps);

    drawhighlightobjects(world, highlightimage, highlightsize, minihighlightimage, minihighlightsize);
}

// This draws the highlight objects.

void drawhighlightobjects(planet& world, sf::Image& highlightimage, int highlightsize, sf::Image& minihighlightimage, int& minihighlightsize)
{
    sf::Color currenthighlightcolour;

    currenthighlightcolour.r = world.highlight1();
    currenthighlightcolour.g = world.highlight2();
    currenthighlightcolour.b = world.highlight3();

    // Do the highlight point

    for (int i = 0; i < highlightsize; i++)
    {
        for (int j = 0; j < highlightsize; j++)
            highlightimage.setPixel(i, j, currenthighlightcolour);
    }

    // And the minimap highlight

    int size = world.size();

    if (size == 0)
        minihighlightsize = 32;

    if (size == 1)
        minihighlightsize = 16;

    if (size == 2)
        minihighlightsize = 8;

    minihighlightimage.create(32, 32, sf::Color::Transparent);

    for (int i = 0; i < minihighlightsize; i++)
    {
        minihighlightimage.setPixel(i, 0, currenthighlightcolour);
        minihighlightimage.setPixel(i, minihighlightsize - 1, currenthighlightcolour);
    }

    for (int j = 0; j < minihighlightsize; j++)
    {
        minihighlightimage.setPixel(0, j, currenthighlightcolour);
        minihighlightimage.setPixel(minihighlightsize - 1, j, currenthighlightcolour);
    }
}

// This prints out an update text.

void updatereport(string text)
{
    if (timedreportingenabled == false)
    {
        cout << text << endl;
        return;
    }

    flushtimedreport();

    if (text.empty())
    {
        cout << endl;
        return;
    }

    timedreporttext = text;
    timedreportstart = chrono::steady_clock::now();
    hastimedreport = true;

    cout << text << "..." << flush;
}

void begintimedreporting()
{
    timedreportingenabled = true;
    hastimedreport = false;
    timedreporttext.clear();
}

void endtimedreporting()
{
    flushtimedreport();
    timedreportingenabled = false;
}

// This makes a button in a standard size and indentation.

bool standardbutton(const char* label)
{
    ImGui::SetCursorPosX(20);

    return ImGui::Button(label, ImVec2(120.0f, 0.0f));
}

// These functions draw a global map image (ready to be applied to a texture).

void drawglobalmapimage(mapviewenum mapview, planet& world, mapcache& maps)
{
    maplayer& layer = getmaplayer(maps, mapview);

    if (layer.created)
        return;

    getmapviewdefinition(mapview).drawglobal(world, layer);

    layer.created = true;
}

void drawallglobalmapimages(planet& world, mapcache& maps)
{
    for (mapviewenum mapview : allmapviews)
        drawglobalmapimage(mapview, world, maps);
}

void applyglobalmapview(mapviewenum mapview, planet& world, mapcache& maps, sf::Texture& texture, sf::Sprite& sprite, sf::Sprite* minimap)
{
    drawglobalmapimage(mapview, world, maps);
    updateTextureFromImage(texture, getdisplaymapimage(maps, mapview));
    sprite.setTexture(texture);

    if (minimap != nullptr)
        minimap->setTexture(texture);
}

void drawglobalelevationmapimage(planet& world, maplayer& layer)
{
    sf::Image& globalelevationimage = layer.image;
    sf::Image& displayglobalelevationimage = layer.displayimage;

    int width = world.width();
    int height = world.height();
    const MapGradientRuntime gradient = buildmapgradient(world, mapgradientelevation);

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            const int altitude = world.map(i, j) - world.sealevel();
            globalelevationimage.setPixel(i, j, getelevationmapcolour(gradient, altitude));
        }
    }

    createnearestdisplayimage(globalelevationimage.getPixelsPtr(), static_cast<int>(globalelevationimage.getSize().x), width, displayglobalelevationimage);
    applyglobaloutlineoverlay(world, globalelevationimage);
    applyglobaldisplayoutlineoverlay(world, displayglobalelevationimage);
}

void drawglobaltemperaturemapimage(planet& world, maplayer& layer)
{
    sf::Image& globaltemperatureimage = layer.image;
    sf::Image& displayglobaltemperatureimage = layer.displayimage;

    int width = world.width();
    int height = world.height();
    const MapGradientRuntime gradient = buildmapgradient(world, mapgradienttemperature);

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            const int jantemp = world.jantemp(i, j);
            const int jultemp = world.jultemp(i, j);
            const int aprtemp = world.aprtemp(i, j);
            const int temperature = ((jantemp + aprtemp + jultemp + aprtemp) / 4) + 10;

            globaltemperatureimage.setPixel(i, j, gettemperaturemapcolour(gradient, temperature));
        }
    }

    createnearestdisplayimage(globaltemperatureimage.getPixelsPtr(), static_cast<int>(globaltemperatureimage.getSize().x), width, displayglobaltemperatureimage);
    applyglobaloutlineoverlay(world, globaltemperatureimage);
    applyglobaldisplayoutlineoverlay(world, displayglobaltemperatureimage);
}

void drawglobalprecipitationmapimage(planet& world, maplayer& layer)
{
    sf::Image& globalprecipitationimage = layer.image;
    sf::Image& displayglobalprecipitationimage = layer.displayimage;

    int width = world.width();
    int height = world.height();
    const MapGradientRuntime gradient = buildmapgradient(world, mapgradientprecipitation);

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            const int rainfall = (world.summerrain(i, j) + world.winterrain(i, j)) / 2;
            globalprecipitationimage.setPixel(i, j, getprecipitationmapcolour(gradient, rainfall));
        }
    }

    createnearestdisplayimage(globalprecipitationimage.getPixelsPtr(), static_cast<int>(globalprecipitationimage.getSize().x), width, displayglobalprecipitationimage);
    applyglobaloutlineoverlay(world, globalprecipitationimage);
    applyglobaldisplayoutlineoverlay(world, displayglobalprecipitationimage);
}

void drawglobalclimatemapimage(planet& world, maplayer& layer)
{
    sf::Image& globalclimateimage = layer.image;
    sf::Image& displayglobalclimateimage = layer.displayimage;

    int width = world.width();
    int height = world.height();

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (world.sea(i, j))
                globalclimateimage.setPixel(i, j, getclimateseacolour(world, world.seaice(i, j)));
            else
                globalclimateimage.setPixel(i, j, getclimatecolours(world, world.climate(i, j)));
        }
    }

    createnearestdisplayimage(globalclimateimage.getPixelsPtr(), static_cast<int>(globalclimateimage.getSize().x), width, displayglobalclimateimage);
    applyglobaloutlineoverlay(world, globalclimateimage);
    applyglobaldisplayoutlineoverlay(world, displayglobalclimateimage);
}

void drawglobalbiomemapimage(planet& world, maplayer& layer)
{
    sf::Image& globalbiomeimage = layer.image;
    sf::Image& displayglobalbiomeimage = layer.displayimage;

    int width = world.width();
    int height = world.height();

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            const bool water = world.sea(i, j) || world.truelake(i, j) != 0;

            if (water)
            {
                const int seaice = world.sea(i, j) ? world.seaice(i, j) : 0;
                globalbiomeimage.setPixel(i, j, getclimateseacolour(world, seaice));
                continue;
            }

            globalbiomeimage.setPixel(i, j, getbiomemapcolour(world, getholdridgebiomeslot(
                static_cast<float>(world.jantemp(i, j)),
                static_cast<float>(world.aprtemp(i, j)),
                static_cast<float>(world.jultemp(i, j)),
                static_cast<float>(world.octtemp(i, j)),
                static_cast<float>(world.janrain(i, j)),
                static_cast<float>(world.aprrain(i, j)),
                static_cast<float>(world.julrain(i, j)),
                static_cast<float>(world.octrain(i, j)))));
        }
    }

    createnearestdisplayimage(globalbiomeimage.getPixelsPtr(), static_cast<int>(globalbiomeimage.getSize().x), width, displayglobalbiomeimage);
    applyglobaloutlineoverlay(world, globalbiomeimage);
    applyglobaldisplayoutlineoverlay(world, displayglobalbiomeimage);
}

void drawglobalriversmapimage(planet& world, maplayer& layer)
{
    sf::Image& globalriversimage = layer.image;
    sf::Image& displayglobalriversimage = layer.displayimage;

    int width = world.width();
    int height = world.height();
    const MapGradientRuntime gradient = buildmapgradient(world, mapgradientriverflow);

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            const bool hasriver = world.riveraveflow(i, j) > 0 && world.sea(i, j) == 0;
            const bool hasdelta = world.deltadir(i, j) != 0 && world.sea(i, j) == 0;
            const bool haslake = world.truelake(i, j) != 0;
            const bool hasvolcano = world.volcano(i, j) > 0;
            sf::Color colour = getworldcolour(
                world.rivermapcolour(rivermapbackground, 0),
                world.rivermapcolour(rivermapbackground, 1),
                world.rivermapcolour(rivermapbackground, 2));

            if (hasriver)
            {
                colour = getriverflowcolour(gradient, world.riveraveflow(i, j));
            }
            else if (hasdelta)
            {
                const int flow = (world.deltajan(i, j) + world.deltajul(i, j)) / 2;
                colour = getriverflowcolour(gradient, flow);
            }
            else if (haslake)
            {
                if (world.showrivermapfeature(rivermapshowlakes))
                {
                    colour = getworldcolour(
                        world.rivermapcolour(rivermaplake, 0),
                        world.rivermapcolour(rivermaplake, 1),
                        world.rivermapcolour(rivermaplake, 2));
                }
            }
            else
            {
                const int special = world.special(i, j);

                if (special == 110 && world.showrivermapfeature(rivermapshowsaltpans))
                {
                    colour = getworldcolour(
                        world.rivermapcolour(rivermapsaltpan, 0),
                        world.rivermapcolour(rivermapsaltpan, 1),
                        world.rivermapcolour(rivermapsaltpan, 2));
                }
                else if (special == 120 && world.showrivermapfeature(rivermapshowsand))
                {
                    colour = getworldcolour(
                        world.rivermapcolour(rivermapsand, 0),
                        world.rivermapcolour(rivermapsand, 1),
                        world.rivermapcolour(rivermapsand, 2));
                }
                else if (special >= 130 && world.showrivermapfeature(rivermapshowwetlands))
                {
                    colour = getworldcolour(
                        world.rivermapcolour(rivermapwetlands, 0),
                        world.rivermapcolour(rivermapwetlands, 1),
                        world.rivermapcolour(rivermapwetlands, 2));
                }
            }

            if (hasvolcano && world.showrivermapfeature(rivermapshowvolcanoes))
            {
                colour = getworldcolour(
                    world.rivermapcolour(rivermapvolcano, 0),
                    world.rivermapcolour(rivermapvolcano, 1),
                    world.rivermapcolour(rivermapvolcano, 2));
            }

            globalriversimage.setPixel(i, j, colour);
        }
    }

    createriversdisplayimage(globalriversimage.getPixelsPtr(), static_cast<int>(globalriversimage.getSize().x), width, world.size(), displayglobalriversimage);
    applyglobaloutlineoverlay(world, globalriversimage);
    applyglobaldisplayoutlineoverlay(world, displayglobalriversimage);
}

void drawglobalreliefmapimage(planet& world, maplayer& layer)
{
    if (cudaglobalrenderersavailable() && drawglobalreliefmapimagecuda(world, layer))
    {
        applyglobaloutlineoverlay(world, layer.image);
        applyglobaldisplayoutlineoverlay(world, layer.displayimage);
        return;
    }

    sf::Image& globalreliefimage = layer.image;
    sf::Image& displayglobalreliefimage = layer.displayimage;

    int width = world.width();
    int height = world.height();
    int sealevel = world.sealevel();
    int type = world.type();

    int colour1, colour2, colour3;

    int landdiv = ((world.maxelevation() - sealevel) / 2) / 255;
    int seadiv = sealevel / 255;

    int div = world.maxelevation() / 255;
    int base = world.maxelevation() / 4;

    int mult = world.maxriverflow() / 255;

    int minriverflow = world.minriverflowglobal(); // Rivers of this size or larger will be shown on the map.

    float landshading = world.landshading();
    float lakeshading = world.lakeshading();
    float seashading = world.seashading();

    int shadingdir = world.shadingdir();
    bool colourcliffs = world.colourcliffs();

    static vector<vector<short>> reliefmap1(ARRAYWIDTH, vector<short>(ARRAYHEIGHT, 0));
    static vector<vector<short>> reliefmap2(ARRAYWIDTH, vector<short>(ARRAYHEIGHT, 0));
    static vector<vector<short>> reliefmap3(ARRAYWIDTH, vector<short>(ARRAYHEIGHT, 0));
    static vector<sf::Uint8> reliefpixels;
    const vector<sf::Color> tundracolours = buildtundracolours(world);

    const int reliefimagewidth = static_cast<int>(globalreliefimage.getSize().x);
    const int reliefimageheight = static_cast<int>(globalreliefimage.getSize().y);
    const int seaiceappearance = world.seaiceappearance();
    const int snowchange = world.snowchange();

    resizepixelbuffer(reliefpixels, reliefimagewidth, reliefimageheight);

    fillgrid(reliefmap1, static_cast<short>(0));
    fillgrid(reliefmap2, static_cast<short>(0));
    fillgrid(reliefmap3, static_cast<short>(0));

    int var = 0; // Amount colours may be varied to make the map seem more speckledy.

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            var = 10;

            const int mapvalue = world.map(i, j);
            const int sea = world.sea(i, j);
            const int special = world.special(i, j);
            const int riverflow = (world.riverjan(i, j) + world.riverjul(i, j)) / 2;
            const bool truelake = world.truelake(i, j) != 0;
            const bool riftlake = world.riftlakesurface(i, j) != 0;

            if (sea == 1)
            {
                if ((world.seaice(i, j) == 2 && (seaiceappearance == 1 || seaiceappearance == 3)) || (world.seaice(i, j) == 1 && seaiceappearance == 3)) // Sea ice
                {
                    colour1 = world.seaice1();
                    colour2 = world.seaice2();
                    colour3 = world.seaice3();

                    var = 0;
                }
                else
                {
                    colour1 = (world.ocean1() * mapvalue + world.deepocean1() * (sealevel - mapvalue)) / sealevel;
                    colour2 = (world.ocean2() * mapvalue + world.deepocean2() * (sealevel - mapvalue)) / sealevel;
                    colour3 = (world.ocean3() * mapvalue + world.deepocean3() * (sealevel - mapvalue)) / sealevel;

                    var = 5;
                }
            }
            else
            {
                if (riverflow >= minriverflow)
                {
                    colour1 = world.river1();
                    colour2 = world.river2();
                    colour3 = world.river3();
                }
                else
                {
                    if (special == 110) // Salt pan
                    {
                        colour1 = world.saltpan1();
                        colour2 = world.saltpan2();
                        colour3 = world.saltpan3();

                        var = 20;
                    }
                    else
                    {
                        int avetemp = world.avetemp(i, j) + 10;

                        // First, adjust the base colours depending on temperature.

                        int thisbase1, thisbase2, thisbase3, newdesert1, newdesert2, newdesert3;

                        if (avetemp > 30)
                        {
                            thisbase1 = world.base1();
                            thisbase2 = world.base2();
                            thisbase3 = world.base3();

                            newdesert1 = world.desert1();
                            newdesert2 = world.desert2();
                            newdesert3 = world.desert3();
                        }
                        else
                        {
                            int hotno = avetemp / 3;
                            int coldno = 10 - hotno;

                            thisbase1 = (hotno * world.base1() + coldno * world.basetemp1()) / 10;
                            thisbase2 = (hotno * world.base2() + coldno * world.basetemp2()) / 10;
                            thisbase3 = (hotno * world.base3() + coldno * world.basetemp3()) / 10;
                        }

                        if (avetemp > 30)
                        {
                            newdesert1 = world.desert1();
                            newdesert2 = world.desert2();
                            newdesert3 = world.desert3();
                        }
                        else
                        {
                            if (avetemp <= 10)
                            {
                                newdesert1 = world.colddesert1();
                                newdesert2 = world.colddesert2();
                                newdesert3 = world.colddesert3();
                            }
                            else
                            {
                                int hotno = avetemp - 10;
                                int coldno = 20 - hotno;

                                newdesert1 = (hotno * world.desert1() + coldno * world.colddesert1()) / 20;
                                newdesert2 = (hotno * world.desert2() + coldno * world.colddesert2()) / 20;
                                newdesert3 = (hotno * world.desert3() + coldno * world.colddesert3()) / 20;
                            }
                        }

                        // Now, adjust for the presence of monsoon.

                        float winterrain = (float)world.winterrain(i, j);
                        float summerrain = (float)world.summerrain(i, j);

                        float totalrain = winterrain + summerrain;

                        float monsoon = 0.0f;

                        if (winterrain < 1.0f)
                            winterrain = 1.0f;

                        if (winterrain < summerrain)
                        {
                            monsoon = summerrain - winterrain;

                            monsoon = monsoon / 1000.0f; // 410

                            if (monsoon > 0.99f)
                                monsoon = 0.99f;
                        }

                        // The closer it is to tropical rainforest, the more we intensify the rain effect.

                        float rainforestmult = (float)world.mintemp(i, j) / 18.0f; //9.0f;

                        rainforestmult = rainforestmult * (float)world.winterrain(i, j) / 80.0f;

                        if (rainforestmult < 1.0f)
                            rainforestmult = 1.0f;

                        totalrain = totalrain * rainforestmult;

                        // Now adjust the colours for height.

                        int mapelev = mapvalue - sealevel;
                        int desertmapelev = mapelev; // We won't mess about with this one.

                        // If this setting is chosen, pretend that the elevation is much lower for flat areas.

                        if (colourcliffs == 1)
                        {
                            int biggestslope = 0;

                            for (int k = i - 1; k <= i + 1; k++)
                            {
                                int kk = k;

                                if (kk<0 || kk>width)
                                    kk = wrap(kk, width);

                                for (int l = j - 1; l <= j + 1; l++)
                                {
                                    if (l >= 0 && l <= height)
                                    {
                                        int thisslope = mapelev + sealevel - world.map(kk, l);

                                        if (thisslope > biggestslope)
                                            biggestslope = thisslope;
                                    }
                                }
                            }

                            biggestslope = biggestslope - 240; // 180

                            if (biggestslope < 0)
                                biggestslope = 0;

                            float adjustedelev = (float)mapelev;

                            adjustedelev = adjustedelev * (biggestslope / 240.0f);

                            mapelev = (int)adjustedelev;
                        }

                        int newbase1, newbase2, newbase3, newgrass1, newgrass2, newgrass3;

                        if (desertmapelev > 2000)
                        {
                            newdesert1 = world.highdesert1();
                            newdesert2 = world.highdesert2();
                            newdesert3 = world.highdesert3();
                        }
                        else
                        {
                            int highno = desertmapelev / 50;
                            int lowno = 40 - highno;

                            newdesert1 = (highno * world.highdesert1() + lowno * newdesert1) / 40;
                            newdesert2 = (highno * world.highdesert2() + lowno * newdesert2) / 40;
                            newdesert3 = (highno * world.highdesert3() + lowno * newdesert3) / 40;
                        }

                        if (mapelev > 3000)
                        {
                            newbase1 = world.highbase1();
                            newbase2 = world.highbase2();
                            newbase3 = world.highbase3();

                            newgrass1 = world.highbase1();
                            newgrass2 = world.highbase2();
                            newgrass3 = world.highbase3();
                        }
                        else
                        {
                            int highno = mapelev / 75;
                            int lowno = 40 - highno;

                            newbase1 = (highno * world.highbase1() + lowno * thisbase1) / 40;
                            newbase2 = (highno * world.highbase2() + lowno * thisbase2) / 40;
                            newbase3 = (highno * world.highbase3() + lowno * thisbase3) / 40;

                            newgrass1 = (highno * world.highbase1() + lowno * world.grass1()) / 40;
                            newgrass2 = (highno * world.highbase2() + lowno * world.grass2()) / 40;
                            newgrass3 = (highno * world.highbase3() + lowno * world.grass3()) / 40;
                        }

                        // Now we need to mix these according to how dry the location is.

                        if (totalrain > 800.0f) // 800
                        {
                            colour1 = newbase1;
                            colour2 = newbase2;
                            colour3 = newbase3;
                        }
                        else
                        {
                            if (totalrain > 200.0f) //400
                            {
                                int wetno = ((int)totalrain - 200) / 40; //400 20
                                if (wetno > 20)
                                    wetno = 20;
                                int dryno = 20 - wetno;

                                colour1 = (wetno * newbase1 + dryno * newgrass1) / 20;
                                colour2 = (wetno * newbase2 + dryno * newgrass2) / 20;
                                colour3 = (wetno * newbase3 + dryno * newgrass3) / 20;
                            }
                            else
                            {
                                float ftotalrain = 200.0f - totalrain; // 400

                                ftotalrain = ftotalrain / 200.0f; // 400

                                int powamount = (int)totalrain - 150; // 350 This is to make a smoother transition.

                                if (powamount < 3)
                                    powamount = 3;

                                ftotalrain = (float)pow(ftotalrain, powamount);

                                ftotalrain = ftotalrain * 200.0f; // 400

                                totalrain = 200.0f - ftotalrain; // 400

                                int wetno = (int)totalrain;
                                int dryno = 200 - wetno;

                                colour1 = (wetno * newgrass1 + dryno * newdesert1) / 200;
                                colour2 = (wetno * newgrass2 + dryno * newdesert2) / 200;
                                colour3 = (wetno * newgrass3 + dryno * newdesert3) / 200;
                            }
                        }

                        // Now we need to alter that according to how cold the location is.

                        if (avetemp <= 0)
                        {
                            colour1 = world.cold1();
                            colour2 = world.cold2();
                            colour3 = world.cold3();
                        }
                        else
                        {
                            // Get the right tundra colour, depending on latitude.

                            const sf::Color thistundra = tundracolours[j];
                            const int thistundra1 = thistundra.r;
                            const int thistundra2 = thistundra.g;
                            const int thistundra3 = thistundra.b;

                            if (snowchange == 1) // Abrupt transition
                            {
                                if (avetemp < 20)
                                {
                                    if (avetemp < 6)
                                    {
                                        colour1 = world.cold1();
                                        colour2 = world.cold2();
                                        colour3 = world.cold3();
                                    }
                                    else
                                    {
                                        if (avetemp < 10)
                                        {
                                            colour1 = thistundra1;
                                            colour2 = thistundra2;
                                            colour3 = thistundra3;
                                        }
                                        else
                                        {
                                            int hotno = avetemp - 10;
                                            int coldno = 10 - hotno;

                                            colour1 = (hotno * colour1 + coldno * thistundra1) / 10;
                                            colour2 = (hotno * colour2 + coldno * thistundra2) / 10;
                                            colour3 = (hotno * colour3 + coldno * thistundra3) / 10;
                                        }
                                    }
                                }
                            }

                            if (snowchange == 2) // Speckled transition
                            {
                                if (avetemp < 20)
                                {
                                    if (avetemp < 6)
                                    {
                                        colour1 = world.cold1();
                                        colour2 = world.cold2();
                                        colour3 = world.cold3();
                                    }
                                    else
                                    {
                                        if (avetemp < 10)
                                        {
                                            if (random(6, 10) < avetemp)
                                            {
                                                colour1 = thistundra1;
                                                colour2 = thistundra2;
                                                colour3 = thistundra3;
                                            }
                                            else
                                            {
                                                colour1 = world.cold1();
                                                colour2 = world.cold2();
                                                colour3 = world.cold3();
                                            }
                                        }
                                        else
                                        {
                                            int hotno = avetemp - 10;
                                            int coldno = 10 - hotno;

                                            colour1 = (hotno * colour1 + coldno * thistundra1) / 10;
                                            colour2 = (hotno * colour2 + coldno * thistundra2) / 10;
                                            colour3 = (hotno * colour3 + coldno * thistundra3) / 10;
                                        }
                                    }
                                }
                            }

                            if (snowchange == 3) // Gradual transition
                            {
                                if (avetemp < 20)
                                {
                                    if (avetemp < 10)
                                    {
                                        int hotno = avetemp;
                                        int coldno = 10 - hotno;

                                        colour1 = (hotno * thistundra1 + coldno * world.cold1()) / 10;
                                        colour2 = (hotno * thistundra2 + coldno * world.cold2()) / 10;
                                        colour3 = (hotno * thistundra3 + coldno * world.cold3()) / 10;
                                    }
                                    else
                                    {
                                        int hotno = avetemp - 10;
                                        int coldno = 10 - hotno;

                                        colour1 = (hotno * colour1 + coldno * thistundra1) / 10;
                                        colour2 = (hotno * colour2 + coldno * thistundra2) / 10;
                                        colour3 = (hotno * colour3 + coldno * thistundra3) / 10;
                                    }
                                }
                            }
                        }

                        // Now add sand, if need be.

                        if (special == 120)
                        {
                            colour1 = (colour1 * 2 + world.erg1()) / 3;
                            colour2 = (colour2 * 2 + world.erg2()) / 3;
                            colour3 = (colour3 * 2 + world.erg3()) / 3;
                        }

                        // Now wetlands.

                        if (special >= 130 && special < 140)
                        {
                            colour1 = (colour1 * 2 + world.wetlands1()) / 3;
                            colour2 = (colour2 * 2 + world.wetlands2()) / 3;
                            colour3 = (colour3 * 2 + world.wetlands3()) / 3;
                        }
                    }
                }
            }

            if (world.sea(i, j) == 1)
            {
                int amount = randomsign(random(0, var));

                colour1 = colour1 + amount;
                colour2 = colour2 + amount;
                colour3 = colour3 + amount;
            }
            else
            {
                colour1 = colour1 + randomsign(random(0, var));
                colour2 = colour2 + randomsign(random(0, var));
                colour3 = colour3 + randomsign(random(0, var));

                if (truelake)
                {
                    colour1 = world.lake1();
                    colour2 = world.lake2();
                    colour3 = world.lake3();
                }
                else
                {
                    if (riftlake)
                    {
                        colour1 = world.lake1();
                        colour2 = world.lake2();
                        colour3 = world.lake3();
                    }
                }
            }

            if (colour1 > 255)
                colour1 = 255;
            if (colour2 > 255)
                colour2 = 255;
            if (colour3 > 255)
                colour3 = 255;

            if (colour1 < 0)
                colour1 = 0;
            if (colour2 < 0)
                colour2 = 0;
            if (colour3 < 0)
                colour3 = 0;

            reliefmap1[i][j] = colour1;
            reliefmap2[i][j] = colour2;
            reliefmap3[i][j] = colour3;
        }
    }

    // Now apply that to the image, adding shading for slopes where appropriate.

    short r, g, b;

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            r = reliefmap1[i][j];
            g = reliefmap2[i][j];
            b = reliefmap3[i][j];
            const int mapvalue = world.map(i, j);
            const int sea = world.sea(i, j);
            const bool truelake = world.truelake(i, j) == 1;

            if (world.noshade(i, j) == 0)
            {
                bool goahead = 1;

                if ((world.seaice(i, j) == 2 && (seaiceappearance == 1 || seaiceappearance == 3)) || (world.seaice(i, j) == 1 && seaiceappearance == 3)) // Sea ice
                    goahead = 0;

                if (goahead == 1 || sea == 0)
                {
                    auto getslopecached = [&](int xx, int yy, int currentelevation, int& slope) -> bool
                    {
                        if (yy < 0 || yy > height)
                            return false;

                        if (xx < 0 || xx > width)
                            xx = wrap(xx, width);

                        slope = world.map(xx, yy) - currentelevation;
                        return true;
                    };

                    int slope1 = 0;
                    int slope2 = 0;
                    bool hasslope1 = false;
                    bool hasslope2 = false;

                    if (shadingdir == 2)
                    {
                        hasslope1 = getslopecached(i - 1, j, mapvalue, slope1);
                        hasslope2 = getslopecached(i, j + 1, mapvalue, slope2);
                    }

                    if (shadingdir == 4)
                    {
                        hasslope1 = getslopecached(i - 1, j, mapvalue, slope1);
                        hasslope2 = getslopecached(i, j - 1, mapvalue, slope2);
                    }

                    if (shadingdir == 6)
                    {
                        hasslope1 = getslopecached(i + 1, j, mapvalue, slope1);
                        hasslope2 = getslopecached(i, j - 1, mapvalue, slope2);
                    }

                    if (shadingdir == 8)
                    {
                        hasslope1 = getslopecached(i + 1, j, mapvalue, slope1);
                        hasslope2 = getslopecached(i, j + 1, mapvalue, slope2);
                    }

                    if (hasslope1 && hasslope2)
                    {
                        int totalslope = (slope1 + slope2) / 10;

                        if (totalslope > 40)
                            totalslope = 40;

                        if (totalslope < -40)
                            totalslope = -40;

                        if (sea == 1)
                            totalslope = (int)((float)totalslope * (seashading * 2.0f));
                        else
                        {
                            if (truelake)
                                totalslope = (int)((float)totalslope * (lakeshading * 2.0f));
                            else
                                totalslope = (int)((float)totalslope * (landshading * 2.0f));
                        }

                        if (mapvalue <= sealevel && world.oceanrifts(i, j) == 0) // Reduce the shading effect around ocean ridges.
                        {
                            int amount = 1;
                            int amount2 = 3;
                            bool found = 0;
                            bool ignore = 0;

                            for (int k = i - amount2; k <= i + amount2; k++) // don't do this around the rift itself
                            {
                                int kk = k;

                                if (kk<0 || kk>width)
                                    kk = wrap(kk, width);

                                for (int l = j - amount2; l <= j + amount2; l++)
                                {
                                    if (l >= 0 && l <= height)
                                    {
                                        if (world.oceanrifts(kk, l) != 0)
                                        {
                                            ignore = 1;
                                            k = i + amount2;
                                            l = j + amount2;
                                        }
                                    }
                                }
                            }

                            if (ignore == 0)
                            {
                                for (int k = i - amount; k <= i + amount; k++)
                                {
                                    int kk = k;

                                    if (kk<0 || kk>width)
                                        kk = wrap(kk, width);

                                    for (int l = j - amount; l <= j + amount; l++)
                                    {
                                        if (l >= 0 && l <= height)
                                        {
                                            if (world.oceanridges(kk, l) != 0)
                                            {
                                                found = 1;
                                                k = i + amount;
                                                l = j + amount;
                                            }
                                        }
                                    }
                                }

                                if (found == 1)
                                    totalslope = totalslope / 4;
                            }
                        }


                        r = r + totalslope;
                        g = g + totalslope;
                        b = b + totalslope;
                    }

                    if (r < 0)
                        r = 0;
                    if (g < 0)
                        g = 0;
                    if (b < 0)
                        b = 0;

                    if (r > 255)
                        r = 255;
                    if (g > 255)
                        g = 255;
                    if (b > 255)
                        b = 255;
                }
            }

            setpixel(reliefpixels, reliefimagewidth, i, j, static_cast<sf::Uint8>(r), static_cast<sf::Uint8>(g), static_cast<sf::Uint8>(b));
        }
    }

    globalreliefimage.create(reliefimagewidth, reliefimageheight, reliefpixels.data());
    createnearestdisplayimage(reliefpixels.data(), reliefimagewidth, width, displayglobalreliefimage);
    applyglobaloutlineoverlay(world, globalreliefimage);
    applyglobaldisplayoutlineoverlay(world, displayglobalreliefimage);
}

// This function gets colours for drawing climate maps.

sf::Color getclimatecolours(const planet& world, short climate)
{
    return getclimatecoloursforworld(world, climate);
}

// These functions draw a regional map image (ready to be applied to a texture).

void drawregionalmapimage(mapviewenum mapview, planet& world, region& region, mapcache& maps)
{
    maplayer& layer = getmaplayer(maps, mapview);

    if (layer.created)
        return;

    getmapviewdefinition(mapview).drawregional(world, region, layer);

    layer.created = true;
}

void drawallregionalmapimages(planet& world, region& region, mapcache& maps)
{
    for (mapviewenum mapview : allmapviews)
        drawregionalmapimage(mapview, world, region, maps);
}

void applyregionalmapview(mapviewenum mapview, planet& world, region& region, mapcache& maps, sf::Texture& texture, sf::Sprite& sprite)
{
    drawregionalmapimage(mapview, world, region, maps);
    updateTextureFromImage(texture, getmapimage(maps, mapview));
    sprite.setTexture(texture);
}

void drawregionalelevationmapimage(planet& world, region& region, maplayer& layer)
{
    sf::Image& regionalelevationimage = layer.image;

    int origregwidthbegin = region.regwidthbegin();
    int origregwidthend = region.regwidthend();
    int origregheightbegin = region.regheightbegin();
    int origregheightend = region.regheightend();

    int regwidthbegin = origregwidthbegin;
    int regwidthend = origregwidthend;
    int regheightbegin = origregheightbegin;
    int regheightend = origregheightend;
    const MapGradientRuntime gradient = buildmapgradient(world, mapgradientelevation);

    for (int i = regwidthbegin; i <= regwidthend; i++)
    {
        for (int j = regheightbegin; j <= regheightend; j++)
        {
            int heightpoint = region.map(i, j);

            if (region.special(i, j) > 100 && region.special(i, j) < 130)
                heightpoint = region.lakesurface(i, j);

            regionalelevationimage.setPixel(i - origregwidthbegin, j - origregheightbegin, getelevationmapcolour(gradient, heightpoint - world.sealevel()));
        }
    }

    applyregionaloutlineoverlay(world, region, regionalelevationimage);
}

void drawregionaltemperaturemapimage(planet& world, region& region, maplayer& layer)
{
    sf::Image& regionaltemperatureimage = layer.image;

    int origregwidthbegin = region.regwidthbegin();
    int origregwidthend = region.regwidthend();
    int origregheightbegin = region.regheightbegin();
    int origregheightend = region.regheightend();

    int regwidthbegin = origregwidthbegin;
    int regwidthend = origregwidthend;
    int regheightbegin = origregheightbegin;
    int regheightend = origregheightend;
    const MapGradientRuntime gradient = buildmapgradient(world, mapgradienttemperature);

    for (int i = regwidthbegin; i <= regwidthend; i++)
    {
        for (int j = regheightbegin; j <= regheightend; j++)
        {
            int temperature = (region.mintemp(i, j) + region.maxtemp(i, j)) / 2;
            temperature = temperature + 10;
            regionaltemperatureimage.setPixel(i - origregwidthbegin, j - origregheightbegin, gettemperaturemapcolour(gradient, temperature));
        }
    }

    applyregionaloutlineoverlay(world, region, regionaltemperatureimage);
}

void drawregionalprecipitationmapimage(planet& world, region& region, maplayer& layer)
{
    sf::Image& regionalprecipitationimage = layer.image;

    int origregwidthbegin = region.regwidthbegin();
    int origregwidthend = region.regwidthend();
    int origregheightbegin = region.regheightbegin();
    int origregheightend = region.regheightend();

    int regwidthbegin = origregwidthbegin;
    int regwidthend = origregwidthend;
    int regheightbegin = origregheightbegin;
    int regheightend = origregheightend;
    const MapGradientRuntime gradient = buildmapgradient(world, mapgradientprecipitation);

    for (int i = regwidthbegin; i <= regwidthend; i++)
    {
        for (int j = regheightbegin; j <= regheightend; j++)
        {
            int rainfall = (region.summerrain(i, j) + region.winterrain(i, j)) / 2;

            sf::Color colour = getprecipitationmapcolour(gradient, rainfall);

            if (region.test(i, j) != 0)
                colour = sf::Color(255, 0, 255);

            regionalprecipitationimage.setPixel(i - origregwidthbegin, j - origregheightbegin, colour);
        }
    }

    applyregionaloutlineoverlay(world, region, regionalprecipitationimage);
}

void drawregionalclimatemapimage(planet& world, region& region, maplayer& layer)
{
    sf::Image& regionalclimateimage = layer.image;

    int origregwidthbegin = region.regwidthbegin();
    int origregwidthend = region.regwidthend();
    int origregheightbegin = region.regheightbegin();
    int origregheightend = region.regheightend();

    int regwidthbegin = origregwidthbegin;
    int regwidthend = origregwidthend;
    int regheightbegin = origregheightbegin;
    int regheightend = origregheightend;

    for (int i = regwidthbegin; i <= regwidthend; i++)
    {
        for (int j = regheightbegin; j <= regheightend; j++)
        {
            if (region.sea(i, j))
                regionalclimateimage.setPixel(i - origregwidthbegin, j - origregheightbegin, getclimateseacolour(world, region.seaice(i, j)));
            else
                regionalclimateimage.setPixel(i - origregwidthbegin, j - origregheightbegin, getclimatecolours(world, region.climate(i, j)));
        }
    }

    applyregionaloutlineoverlay(world, region, regionalclimateimage);
}

void drawregionalbiomemapimage(planet& world, region& region, maplayer& layer)
{
    sf::Image& regionalbiomeimage = layer.image;

    int origregwidthbegin = region.regwidthbegin();
    int origregwidthend = region.regwidthend();
    int origregheightbegin = region.regheightbegin();
    int origregheightend = region.regheightend();

    int regwidthbegin = origregwidthbegin;
    int regwidthend = origregwidthend;
    int regheightbegin = origregheightbegin;
    int regheightend = origregheightend;

    for (int i = regwidthbegin; i <= regwidthend; i++)
    {
        for (int j = regheightbegin; j <= regheightend; j++)
        {
            const bool water = region.sea(i, j) || region.truelake(i, j) != 0;

            if (water)
            {
                const int seaice = region.sea(i, j) ? region.seaice(i, j) : 0;
                regionalbiomeimage.setPixel(i - origregwidthbegin, j - origregheightbegin, getclimateseacolour(world, seaice));
                continue;
            }

            const int globaly = clamp(region.lefty() + j / 16, 0, world.height());

            regionalbiomeimage.setPixel(i - origregwidthbegin, j - origregheightbegin, getbiomemapcolour(world, getholdridgebiomeslot(
                static_cast<float>(region.jantemp(i, j)),
                static_cast<float>(region.aprtemp(i, j, globaly, world.height(), world.tilt(), world.eccentricity(), world.perihelion())),
                static_cast<float>(region.jultemp(i, j)),
                static_cast<float>(region.octtemp(i, j, globaly, world.height(), world.tilt(), world.eccentricity(), world.perihelion())),
                static_cast<float>(region.janrain(i, j)),
                static_cast<float>(region.aprrain(i, j)),
                static_cast<float>(region.julrain(i, j)),
                static_cast<float>(region.octrain(i, j)))));
        }
    }

    applyregionaloutlineoverlay(world, region, regionalbiomeimage);
}

void drawregionalriversmapimage(planet& world, region& region, maplayer& layer)
{
    sf::Image& regionalriversimage = layer.image;

    int origregwidthbegin = region.regwidthbegin();
    int origregwidthend = region.regwidthend();
    int origregheightbegin = region.regheightbegin();
    int origregheightend = region.regheightend();

    int regwidthbegin = origregwidthbegin;
    int regwidthend = origregwidthend;
    int regheightbegin = origregheightbegin;
    int regheightend = origregheightend;
    const MapGradientRuntime gradient = buildmapgradient(world, mapgradientriverflow);

    for (int i = regwidthbegin; i <= regwidthend; i++)
    {
        for (int j = regheightbegin; j <= regheightend; j++)
        {
            const int flow = region.riveraveflow(i, j);
            const bool hasriver = flow > 0;
            const bool haslake = region.truelake(i, j) != 0;
            const bool hasvolcano = region.volcano(i, j) == 1;
            sf::Color colour = getworldcolour(
                world.rivermapcolour(rivermapbackground, 0),
                world.rivermapcolour(rivermapbackground, 1),
                world.rivermapcolour(rivermapbackground, 2));

            if (hasriver)
            {
                colour = getriverflowcolour(gradient, flow);
            }
            else if (haslake)
            {
                if (world.showrivermapfeature(rivermapshowlakes))
                {
                    colour = getworldcolour(
                        world.rivermapcolour(rivermaplake, 0),
                        world.rivermapcolour(rivermaplake, 1),
                        world.rivermapcolour(rivermaplake, 2));
                }
            }
            else
            {
                const int special = region.special(i, j);

                if (special == 110 && world.showrivermapfeature(rivermapshowsaltpans))
                {
                    colour = getworldcolour(
                        world.rivermapcolour(rivermapsaltpan, 0),
                        world.rivermapcolour(rivermapsaltpan, 1),
                        world.rivermapcolour(rivermapsaltpan, 2));
                }
                else if (special >= 130 && special < 140 && world.showrivermapfeature(rivermapshowwetlands))
                {
                    colour = getworldcolour(
                        world.rivermapcolour(rivermapwetlands, 0),
                        world.rivermapcolour(rivermapwetlands, 1),
                        world.rivermapcolour(rivermapwetlands, 2));
                }

                if (region.mud(i, j) && world.showrivermapfeature(rivermapshowmud))
                {
                    colour = getworldcolour(
                        world.rivermapcolour(rivermapmud, 0),
                        world.rivermapcolour(rivermapmud, 1),
                        world.rivermapcolour(rivermapmud, 2));
                }

                if (region.sand(i, j) && world.showrivermapfeature(rivermapshowsand))
                {
                    colour = getworldcolour(
                        world.rivermapcolour(rivermapsand, 0),
                        world.rivermapcolour(rivermapsand, 1),
                        world.rivermapcolour(rivermapsand, 2));
                }

                if (region.shingle(i, j) && world.showrivermapfeature(rivermapshowshingle))
                {
                    colour = getworldcolour(
                        world.rivermapcolour(rivermapshingle, 0),
                        world.rivermapcolour(rivermapshingle, 1),
                        world.rivermapcolour(rivermapshingle, 2));
                }
            }

            if (hasvolcano && world.showrivermapfeature(rivermapshowvolcanoes))
            {
                colour = getworldcolour(
                    world.rivermapcolour(rivermapvolcano, 0),
                    world.rivermapcolour(rivermapvolcano, 1),
                    world.rivermapcolour(rivermapvolcano, 2));
            }

            regionalriversimage.setPixel(i - origregwidthbegin, j - origregheightbegin, colour);
        }
    }

    applyregionaloutlineoverlay(world, region, regionalriversimage);
}

void drawregionalreliefmapimage(planet& world, region& region, maplayer& layer)
{
    sf::Image& regionalreliefimage = layer.image;

    int origregwidthbegin = region.regwidthbegin();
    int origregwidthend = region.regwidthend();
    int origregheightbegin = region.regheightbegin();
    int origregheightend = region.regheightend();

    int regwidthbegin = origregwidthbegin;
    int regwidthend = origregwidthend;
    int regheightbegin = origregheightbegin;
    int regheightend = origregheightend;

    int colour1, colour2, colour3;

    int div = world.maxelevation() / 255;
    int base = world.maxelevation() / 4;

    int mult = world.maxriverflow() / 400;

    int leftx = region.leftx();
    int lefty = region.lefty();

    int rwidth = region.rwidth();
    int rheight = region.rheight();

    int type = world.type();
    int sealevel = world.sealevel();
    int minriverflow = world.minriverflowregional(); // Rivers of this size or larger will be shown on the map.
    int shadingdir = world.shadingdir();

    float landshading = world.landshading();
    float lakeshading = world.lakeshading();
    float seashading = world.seashading();

    float landmarbling = world.landmarbling() * 2;
    float lakemarbling = world.lakemarbling() * 2;
    float seamarbling = world.seamarbling() * 2;

    bool colourcliffs = world.colourcliffs();

    int width = world.width();
    int height = world.height();
    int xleft = 0;
    int xright = 35;
    int ytop = 0;
    int ybottom = 35;

    static vector<vector<short>> reliefmap1(RARRAYWIDTH, vector<short>(RARRAYHEIGHT, 0));
    static vector<vector<short>> reliefmap2(RARRAYWIDTH, vector<short>(RARRAYHEIGHT, 0));
    static vector<vector<short>> reliefmap3(RARRAYWIDTH, vector<short>(RARRAYHEIGHT, 0));
    static vector<sf::Uint8> reliefpixels;

    // Make a fractal based on rainfall, which will be used to add stripes to vary the colours.

    static vector<vector<int>> source(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));
    static vector<vector<int>> stripefractal(RARRAYWIDTH, vector<int>(RARRAYHEIGHT, -5000));
    const vector<sf::Color> tundracolours = buildtundracolours(world);
    const int regionalreliefimagewidth = static_cast<int>(regionalreliefimage.getSize().x);
    const int regionalreliefimageheight = static_cast<int>(regionalreliefimage.getSize().y);
    const int seaiceappearance = world.seaiceappearance();
    const int snowchange = world.snowchange();
    const auto hasvisibleflow = [&](int x, int y)
    {
        return (region.riverjan(x, y) + region.riverjul(x, y)) / 2 >= minriverflow
            || (region.fakejan(x, y) + region.fakejul(x, y)) / 2 >= minriverflow
            || region.deltajan(x, y) != 0
            || region.deltajul(x, y) != 0;
    };

    fillgrid(reliefmap1, static_cast<short>(0));
    fillgrid(reliefmap2, static_cast<short>(0));
    fillgrid(reliefmap3, static_cast<short>(0));
    fillgrid(source, 0);
    fillgrid(stripefractal, -5000);
    resizepixelbuffer(reliefpixels, regionalreliefimagewidth, regionalreliefimageheight);

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (world.wintermountainraindir(i, j) == 0)
                source[i][j] = world.winterrain(i, j);
            else
                source[i][j] = world.wintermountainrain(i, j);
        }
    }

    int var = 0; // Amount colours may be varied to make the map seem more speckledy.

    int lat = 0;
    int lat2 = 0;
    int latminutes = 0;
    int latseconds = 0;
    bool latneg = 0;

    for (int i = regwidthbegin; i <= regwidthend; i++)
    {
        for (int j = regheightbegin; j <= regheightend; j++)
        {
            bool stripe = 1; // Indicates whether we're going to do a stripe here.

            int xx = leftx + (i / 16);
            int yy = lefty + (j / 16); // Coordinates of the relevant global cell.
            const int mapvalue = region.map(i, j);
            const int special = region.special(i, j);
            const int riverflow = (region.riverjan(i, j) + region.riverjul(i, j)) / 2;
            const int fakeflow = (region.fakejan(i, j) + region.fakejul(i, j)) / 2;
            const int sea = region.sea(i, j);

            int avetemp = (region.extramaxtemp(i, j) + region.extramintemp(i, j)) / 2;

            avetemp = avetemp + 1000;

            int totalrain = region.summerrain(i, j) + region.winterrain(i, j);

            if (region.special(i, j) == 120)
                totalrain = totalrain / 4;

            var = 10;

            if (sea == 1)
            {
                if ((region.seaice(i, j) == 2 && (seaiceappearance == 1 || seaiceappearance == 3)) || (region.seaice(i, j) == 1 && seaiceappearance == 3)) // Sea ice
                {
                    colour1 = world.seaice1();
                    colour2 = world.seaice2();
                    colour3 = world.seaice3();

                    stripe = 0;

                    var = 0;
                }
                else
                {
                    colour1 = (world.ocean1() * mapvalue + world.deepocean1() * (sealevel - mapvalue)) / sealevel;
                    colour2 = (world.ocean2() * mapvalue + world.deepocean2() * (sealevel - mapvalue)) / sealevel;
                    colour3 = (world.ocean3() * mapvalue + world.deepocean3() * (sealevel - mapvalue)) / sealevel;

                    var = 5;
                }
            }
            else
            {
                if (riverflow >= minriverflow || fakeflow >= minriverflow)
                {
                    colour1 = world.river1();
                    colour2 = world.river2();
                    colour3 = world.river3();

                    stripe = 0;
                }
                else
                {
                    if (special == 110) // Salt pan
                    {
                        colour1 = world.saltpan1();
                        colour2 = world.saltpan2();
                        colour3 = world.saltpan3();

                        var = 20;

                        stripe = 0;
                    }
                    else
                    {
                        if (region.sand(i, j) || region.shingle(i,j))
                        {
                            if (region.sand(i, j) && region.shingle(i, j) == 0)
                            {
                                colour1 = world.sand1();
                                colour2 = world.sand2();
                                colour3 = world.sand3();

                                var = 5;
                            }

                            if (region.shingle(i, j) && region.sand(i, j) == 0)
                            {
                                colour1 = world.shingle1();
                                colour2 = world.shingle2();
                                colour3 = world.shingle3();

                                var = 40;
                            }

                            if (region.sand(i, j) && region.shingle(i, j))
                            {
                                colour1 = (world.sand1() + world.shingle1()) / 2;
                                colour2 = (world.sand2() + world.shingle2()) / 2;
                                colour3 = (world.sand3() + world.shingle3()) / 2;

                                var = 20;
                            }

                            if (region.mud(i, j))
                            {
                                colour1 = (colour1 + world.mud1()) / 2;
                                colour2 = (colour2 + world.mud2()) / 2;
                                colour3 = (colour3 + world.mud3()) / 2;
                            }
                        }
                        else
                        {
                            if (region.mangrove(i, j) && world.showmangroves())
                            {
                                colour1 = world.mangrove1();
                                colour2 = world.mangrove2();
                                colour3 = world.mangrove3();

                                var = 20;
                            }
                            else
                            {
                                // First, adjust the base colours depending on temperature.

                                int thisbase1, thisbase2, thisbase3, newdesert1, newdesert2, newdesert3;

                                if (avetemp > 3000)
                                {
                                    thisbase1 = world.base1();
                                    thisbase2 = world.base2();
                                    thisbase3 = world.base3();
                                }
                                else
                                {
                                    int hotno = avetemp / 3;
                                    int coldno = 1000 - hotno;

                                    thisbase1 = (hotno * world.base1() + coldno * world.basetemp1()) / 1000;
                                    thisbase2 = (hotno * world.base2() + coldno * world.basetemp2()) / 1000;
                                    thisbase3 = (hotno * world.base3() + coldno * world.basetemp3()) / 1000;
                                }

                                if (avetemp > 30)
                                {
                                    newdesert1 = world.desert1();
                                    newdesert2 = world.desert2();
                                    newdesert3 = world.desert3();
                                }
                                else
                                {
                                    if (avetemp <= 10)
                                    {
                                        newdesert1 = world.colddesert1();
                                        newdesert2 = world.colddesert2();
                                        newdesert3 = world.colddesert3();
                                    }
                                    else
                                    {
                                        int hotno = avetemp - 10;
                                        int coldno = 20 - hotno;

                                        newdesert1 = (hotno * world.desert1() + coldno * world.colddesert1()) / 20;
                                        newdesert2 = (hotno * world.desert2() + coldno * world.colddesert2()) / 20;
                                        newdesert3 = (hotno * world.desert3() + coldno * world.colddesert3()) / 20;
                                    }
                                }

                                // The closer it is to tropical rainforest, the more we intensify the rain effect.

                                float rainforestmult = (float)region.mintemp(i, j) / 18.0f;

                                rainforestmult = rainforestmult * (float)region.winterrain(i, j) / 80.0f; //30.0f;

                                if (rainforestmult < 1.0f)
                                    rainforestmult = 1.0f;

                                totalrain = totalrain * (int)rainforestmult;

                                // Now adjust the colours for height.

                                int mapelev = mapvalue - sealevel;

                                if (special == 110 || special == 120) // If this is going to be an erg or salt pan
                                    mapelev = region.lakesurface(i, j) - sealevel; // Use the surface elevation of the erg or salt pan

                                int desertmapelev = mapelev; // We won't mess about with this one.

                                // If it's flat here and this setting is enabled, pretend the elevation is lower than it really is.

                                if (colourcliffs==1)
                                {
                                    int biggestslope = 0;

                                    for (int k = i - 1; k <= i + 1; k++)
                                    {
                                        if (k >= 0 && k <= rwidth)
                                        {
                                            for (int l = j - 1; l <= j + 1; l++)
                                            {
                                                if (l >= 0 && l <= rheight)
                                                {
                                                    int thisslope = mapelev + sealevel - region.map(k, l);

                                                    if (thisslope > biggestslope)
                                                        biggestslope = thisslope;
                                                }
                                            }
                                        }
                                    }

                                    biggestslope = biggestslope - 240; //180

                                    if (biggestslope < 0)
                                        biggestslope = 0;

                                    float adjustedelev = (float)mapelev;

                                    adjustedelev = adjustedelev * (biggestslope / 240.0f);

                                    if (adjustedelev > (float)mapelev)
                                        adjustedelev = (float)mapelev;

                                    mapelev = (int)adjustedelev;
                                }

                                int newbase1, newbase2, newbase3, newgrass1, newgrass2, newgrass3;

                                if (desertmapelev > 2000)
                                {
                                    newdesert1 = world.highdesert1();
                                    newdesert2 = world.highdesert2();
                                    newdesert3 = world.highdesert3();
                                }
                                else
                                {
                                    int highno = desertmapelev / 50;
                                    int lowno = 40 - highno;

                                    newdesert1 = (highno * world.highdesert1() + lowno * newdesert1) / 40;
                                    newdesert2 = (highno * world.highdesert2() + lowno * newdesert2) / 40;
                                    newdesert3 = (highno * world.highdesert3() + lowno * newdesert3) / 40;
                                }

                                if (mapelev > 3000)
                                {
                                    newbase1 = world.highbase1();
                                    newbase2 = world.highbase2();
                                    newbase3 = world.highbase3();

                                    newgrass1 = world.highbase1();
                                    newgrass2 = world.highbase2();
                                    newgrass3 = world.highbase3();
                                }
                                else
                                {
                                    int highno = mapelev / 75;
                                    int lowno = 40 - highno;

                                    newbase1 = (highno * world.highbase1() + lowno * thisbase1) / 40;
                                    newbase2 = (highno * world.highbase2() + lowno * thisbase2) / 40;
                                    newbase3 = (highno * world.highbase3() + lowno * thisbase3) / 40;

                                    newgrass1 = (highno * world.highbase1() + lowno * world.grass1()) / 40;
                                    newgrass2 = (highno * world.highbase2() + lowno * world.grass2()) / 40;
                                    newgrass3 = (highno * world.highbase3() + lowno * world.grass3()) / 40;
                                }

                                // Now we need to mix these according to how dry the location is.

                                if (region.rivervalley(i, j) == 1 || region.riverdir(i,j)!=0 || region.fakedir(i,j)!=0) // If this is a river valley, it's wetter than it would otherwise be.
                                {
                                    float biggestflow = 0.0f;

                                    for (int k = i - 20; k <= i + 20; k++)
                                    {
                                        for (int l = j - 20; l <= j + 20; l++)
                                        {
                                            if (k >= 0 && k <= rwidth && l >= 0 && l <= rheight)
                                            {
                                                if (region.riverjan(k, l) + region.riverjul(k, l) > biggestflow)
                                                    biggestflow = (float)(region.riverjan(k, l) + region.riverjul(k, l));
                                            }
                                        }
                                    }

                                    if (biggestflow == 0.0f)
                                    {
                                        twointegers nearest = findclosestriverquickly(region, i, j);

                                        if (nearest.x != -1)
                                            biggestflow = (float)(region.riverjan(nearest.x, nearest.y) + region.riverjul(nearest.x, nearest.y));
                                    }

                                    if (biggestflow > 12000.0f)
                                        biggestflow = 1200.0f;

                                    float mult = (float)totalrain;

                                    if (mult < 1.0f)
                                        mult = 1.0f;

                                    if (mult > 1000.0f)
                                        mult = 1000.0f;


                                    biggestflow = biggestflow / mult;
                                    totalrain = totalrain + (int)biggestflow;
                                }

                                if (totalrain > 800) // 600
                                {
                                    colour1 = newbase1;
                                    colour2 = newbase2;
                                    colour3 = newbase3;
                                }
                                else
                                {
                                    if (totalrain > 200) //300
                                    {
                                        int wetno = (totalrain - 200) / 40; //20
                                        if (wetno > 20) // 40
                                            wetno = 20;
                                        int dryno = 20 - wetno;


                                        colour1 = (wetno * newbase1 + dryno * newgrass1) / 20;
                                        colour2 = (wetno * newbase2 + dryno * newgrass2) / 20;
                                        colour3 = (wetno * newbase3 + dryno * newgrass3) / 20;
                                    }
                                    else
                                    {
                                        float ftotalrain = 200.0f - (float)totalrain;

                                        ftotalrain = ftotalrain / 200.0f;

                                        int powamount = totalrain - 150; // This is to make a smoother transition.

                                        if (powamount < 3)
                                            powamount = 3;

                                        ftotalrain = (float)pow((double)ftotalrain, (double)powamount);

                                        ftotalrain = ftotalrain * 200.0f;

                                        totalrain = 200 - (int)ftotalrain;

                                        int wetno = totalrain;
                                        int dryno = 200 - wetno;

                                        colour1 = (wetno * newgrass1 + dryno * newdesert1) / 200;
                                        colour2 = (wetno * newgrass2 + dryno * newdesert2) / 200;
                                        colour3 = (wetno * newgrass3 + dryno * newdesert3) / 200;
                                    }
                                }

                                // Now we need to alter that according to how cold the location is.

                                if (avetemp <= 0 || yy > height - 3) // This is because it has an odd tendency to show the very southernmost tiles in non-cold colours.
                                {
                                    colour1 = world.cold1();
                                    colour2 = world.cold2();
                                    colour3 = world.cold3();
                                }
                                else
                                {
                                    // Get the right tundra colour, depending on latitude.

                                    const sf::Color thistundra = tundracolours[yy];
                                    const int thistundra1 = thistundra.r;
                                    const int thistundra2 = thistundra.g;
                                    const int thistundra3 = thistundra.b;

                                    if (snowchange == 1) // Abrupt transition
                                    {
                                        if (avetemp < 2000)
                                        {
                                            if (avetemp < 600)
                                            {
                                                colour1 = world.cold1();
                                                colour2 = world.cold2();
                                                colour3 = world.cold3();
                                            }
                                            else
                                            {
                                                if (avetemp < 1000)
                                                {
                                                    colour1 = thistundra1;
                                                    colour2 = thistundra2;
                                                    colour3 = thistundra3;
                                                }
                                                else
                                                {
                                                    int hotno = avetemp - 1000;
                                                    int coldno = 1000 - hotno;

                                                    colour1 = (hotno * colour1 + coldno * thistundra1) / 1000;
                                                    colour2 = (hotno * colour2 + coldno * thistundra2) / 1000;
                                                    colour3 = (hotno * colour3 + coldno * thistundra3) / 1000;
                                                }
                                            }
                                        }
                                    }

                                    if (snowchange == 2) // Speckled transition
                                    {
                                        if (avetemp < 2000)
                                        {
                                            if (avetemp < 600)
                                            {
                                                colour1 = world.cold1();
                                                colour2 = world.cold2();
                                                colour3 = world.cold3();
                                            }
                                            else
                                            {
                                                if (avetemp < 1000)
                                                {
                                                    if (random(600, 1000) < avetemp)
                                                    {
                                                        colour1 = thistundra1;
                                                        colour2 = thistundra2;
                                                        colour3 = thistundra3;
                                                    }
                                                    else
                                                    {
                                                        colour1 = world.cold1();
                                                        colour2 = world.cold2();
                                                        colour3 = world.cold3();
                                                    }
                                                }
                                                else
                                                {
                                                    int hotno = avetemp - 1000;
                                                    int coldno = 1000 - hotno;

                                                    colour1 = (hotno * colour1 + coldno * thistundra1) / 1000;
                                                    colour2 = (hotno * colour2 + coldno * thistundra2) / 1000;
                                                    colour3 = (hotno * colour3 + coldno * thistundra3) / 1000;
                                                }
                                            }
                                        }
                                    }

                                    if (snowchange == 3) // Gradual transition
                                    {
                                        if (avetemp < 2000)
                                        {
                                            if (avetemp < 1000)
                                            {
                                                int hotno = avetemp;
                                                int coldno = 1000 - hotno;

                                                colour1 = (hotno * thistundra1 + coldno * world.cold1()) / 1000;
                                                colour2 = (hotno * thistundra2 + coldno * world.cold2()) / 1000;
                                                colour3 = (hotno * thistundra3 + coldno * world.cold3()) / 1000;
                                            }
                                            else
                                            {
                                                int hotno = avetemp - 1000;
                                                int coldno = 1000 - hotno;

                                                colour1 = (hotno * colour1 + coldno * thistundra1) / 1000;
                                                colour2 = (hotno * colour2 + coldno * thistundra2) / 1000;
                                                colour3 = (hotno * colour3 + coldno * thistundra3) / 1000;
                                            }
                                        }
                                    }
                                }

                                // Now add dunes, if need be.

                                if (special == 120)
                                {                                    
                                    colour1 = (colour1 * 6 + world.erg1()) / 7;
                                    colour2 = (colour2 * 6 + world.erg2()) / 7;
                                    colour3 = (colour3 * 6 + world.erg3()) / 7;

                                    var = 10; //4;
                                }

                                // Same thing for mud flats.

                                if (region.mud(i, j))
                                {
                                    colour1 = (colour1 + world.mud1() * 2) / 3;
                                    colour2 = (colour2 + world.mud2() * 2) / 3;
                                    colour3 = (colour3 + world.mud3() * 2) / 3;

                                    var = 10;

                                }

                                // Now wetlands.

                                if (special >= 130 && special < 140)
                                {
                                    colour1 = (colour1 * 2 + world.wetlands1()) / 3;
                                    colour2 = (colour2 * 2 + world.wetlands2()) / 3;
                                    colour3 = (colour3 * 2 + world.wetlands3()) / 3;
                                }
                            }
                        }
                    }
                }
            }

            if (region.sea(i, j) == 1)
            {
                int amount = altrandomsign(altrandom(0, var));

                colour1 = colour1 + amount;
                colour2 = colour2 + amount;
                colour3 = colour3 + amount;
            }
            else
            {
                colour1 = colour1 + altrandomsign(altrandom(0, var));
                colour2 = colour2 + altrandomsign(altrandom(0, var));
                colour3 = colour3 + altrandomsign(altrandom(0, var));

                if (region.truelake(i, j) != 0)
                {
                    colour1 = world.lake1();
                    colour2 = world.lake2();
                    colour3 = world.lake3();

                    for (int k = i - 1; k <= i + 1; k++) // If a river is flowing into/out of the lake, make it slightly river coloured.
                    {
                        for (int l = j - 1; l <= j + 1; l++)
                        {
                            if (k >= 0 && k <= rwidth && l >= 0 && l <= rheight)
                            {
                                if (region.riverdir(k, l) != 0 || region.fakedir(k, l) > 0)
                                {
                                    if (region.lakesurface(k, l) == 0)
                                    {
                                        colour1 = (world.lake1() + world.river1()) / 2;
                                        colour2 = (world.lake2() + world.river2()) / 2;
                                        colour3 = (world.lake3() + world.river3()) / 2;

                                        k = i + 1;
                                        l = j + 1;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            if (region.special(i, j) == 140)
            {
                colour1 = world.glacier1();
                colour2 = world.glacier2();
                colour3 = world.glacier3();

                stripe = 0;
            }

            if (stripe == 1) // If we're doing a stripe here.
            {
                //fast_srand(region.maxtemp(i,j)+region.winterrain(i,j));
                fast_srand(region.maxtemp(i, j) + stripefractal[i][j]);

                float stripevar = (float)avetemp;

                if (stripevar > 5.0f)
                    stripevar = 5.0f;

                if (stripevar < 0.0f)
                    stripevar = 0.0f;

                if (region.special(i, j) > 100)
                    stripevar = 1.0f;

                if (region.sea(i, j) == 1)
                    stripevar = stripevar * seamarbling;
                else
                {
                    if (region.truelake(i, j) == 1)
                        stripevar = stripevar * lakemarbling;
                    else
                        stripevar = stripevar * landmarbling;
                }

                colour1 = colour1 + randomsign(random(0, (int)stripevar));
                colour2 = colour2 + randomsign(random(0, (int)stripevar));
                colour3 = colour3 + randomsign(random(0, (int)stripevar));
            }

            if (colour1 > 255)
                colour1 = 255;
            if (colour2 > 255)
                colour2 = 255;
            if (colour3 > 255)
                colour3 = 255;

            if (colour1 < 0)
                colour1 = 0;
            if (colour2 < 0)
                colour2 = 0;
            if (colour3 < 0)
                colour3 = 0;

            reliefmap1[i][j] = colour1;
            reliefmap2[i][j] = colour2;
            reliefmap3[i][j] = colour3;
        }
    }

    for (int i = regwidthbegin + 1; i <= regwidthend - 1; i++) // Blur ergs, mud flats, and wetlands, and also river valleys..
    {
        for (int j = regheightbegin + 1; j <= regheightend - 1; j++)
        {
            int special = region.special(i, j);

            if ((special == 130 || special == 131 || special == 132 || region.mud(i, j)) && (region.mangrove(i, j) == 0 || world.showmangroves() == 0))
            {
                if (hasvisibleflow(i, j) == false) // Don't do it to the rivers themselves.
                {
                    bool touchesvisibleflow = false;

                    for (int k = i - 1; k <= i + 1 && touchesvisibleflow == false; k++)
                    {
                        for (int l = j - 1; l <= j + 1; l++)
                        {
                            if ((k != i || l != j) && hasvisibleflow(k, l))
                            {
                                touchesvisibleflow = true;
                                break;
                            }
                        }
                    }

                    if (touchesvisibleflow)
                        continue;

                    float colred = 0.0f;
                    float colgreen = 0.0f;
                    float colblue = 0.0f;
                    float crount = 0.0f;

                    for (int k = i - 1; k <= i + 1; k++)
                    {
                        for (int l = j - 1; l <= j + 1; l++)
                        {
                            if (hasvisibleflow(k, l))
                                continue;

                            colred = colred + (float)reliefmap1[k][l];
                            colgreen = colgreen + (float)reliefmap2[k][l];
                            colblue = colblue + (float)reliefmap3[k][l];
                            crount++;
                        }
                    }

                    if (crount <= 0.0f)
                        continue;

                    colred = colred / crount;
                    colgreen = colgreen / crount;
                    colblue = colblue / crount;

                    reliefmap1[i][j] = (short)colred;
                    reliefmap2[i][j] = (short)colgreen;
                    reliefmap3[i][j] = (short)colblue;
                }
            }

            if (region.rivervalley(i, j) == 1 && region.special(i, j) < 130)
            {
                if (hasvisibleflow(i, j) == false)
                {
                    float colred = 0.0f;
                    float colgreen = 0.0f;
                    float colblue = 0.0f;

                    float crount = .0f;

                    for (int k = i - 1; k <= i + 1; k++)
                    {
                        for (int l = j - 1; l <= j + 1; l++)
                        {
                            if (hasvisibleflow(k, l) == false)
                            {
                                colred = colred + (float)reliefmap1[k][l];
                                colgreen = colgreen + (float)reliefmap2[k][l];
                                colblue = colblue +(float)reliefmap3[k][l];

                                crount++;
                            }
                        }
                    }

                    if (crount <= 0.0f)
                        continue;

                    colred = colred / crount;
                    colgreen = colgreen / crount;
                    colblue = colblue / crount;

                    reliefmap1[i][j] = (short)colred;
                    reliefmap2[i][j] = (short)colgreen;
                    reliefmap3[i][j] = (short)colblue;
                }
            }
        }
    }

    // Do the rivers again, as they might have got messed up by the blurring.

    for (int i = regwidthbegin; i <= regwidthend; i++)
    {
        for (int j = regheightbegin; j <= regheightend; j++)
        {
            if (region.sea(i, j) == 0 && region.truelake(i, j) == 0)
            {
                if ((region.riverjan(i, j) + region.riverjul(i, j)) / 2 >= minriverflow || (region.fakejan(i, j) + region.fakejul(i, j)) / 2 >= minriverflow)
                {
                    colour1 = world.river1();
                    colour2 = world.river2();
                    colour3 = world.river3();

                    fast_srand(region.maxtemp(i, j) + region.winterrain(i, j));

                    int stripevar = 5;

                    colour1 = colour1 + randomsign(random(0, stripevar));
                    colour2 = colour2 + randomsign(random(0, stripevar));
                    colour3 = colour3 + randomsign(random(0, stripevar));

                    colour1 = clamp(colour1, 0, 255);
                    colour2 = clamp(colour2, 0, 255);
                    colour3 = clamp(colour3, 0, 255);

                    reliefmap1[i][j] = colour1;
                    reliefmap2[i][j] = colour2;
                    reliefmap3[i][j] = colour3;
                }
            }
        }
    }

    // Now apply that to the image, adding shading for slopes where appropriate.

    short r, g, b;

    for (int i = regwidthbegin; i <= regwidthend; i++)
    {
        for (int j = regheightbegin; j <= regheightend; j++)
        {
            r = reliefmap1[i][j];
            g = reliefmap2[i][j];
            b = reliefmap3[i][j];
            const int mapvalue = region.map(i, j);
            const int sea = region.sea(i, j);
            const bool truelake = region.truelake(i, j) == 1;

            if (region.special(i, j) == 0 && region.riverdir(i, j) == 0 && region.fakedir(i, j) == 0)
            {
                bool goahead = 1;

                if ((region.seaice(i, j) == 2 && (seaiceappearance == 1 || seaiceappearance == 3)) || (region.seaice(i, j) == 1 && seaiceappearance == 3)) // Sea ice
                    goahead = 0;

                if (goahead == 1 || sea == 0)
                {
                    auto getslopecached = [&](int xx, int yy, int currentelevation, int& slope) -> bool
                    {
                        if (xx < 0 || xx > rwidth || yy < 0 || yy > rheight)
                            return false;

                        slope = region.map(xx, yy) - currentelevation;
                        return true;
                    };

                    int slope1 = 0;
                    int slope2 = 0;
                    bool hasslope1 = false;
                    bool hasslope2 = false;

                    if (shadingdir == 2)
                    {
                        hasslope1 = getslopecached(i - 1, j, mapvalue, slope1);
                        hasslope2 = getslopecached(i, j + 1, mapvalue, slope2);
                    }

                    if (shadingdir == 4)
                    {
                        hasslope1 = getslopecached(i - 1, j, mapvalue, slope1);
                        hasslope2 = getslopecached(i, j - 1, mapvalue, slope2);
                    }

                    if (shadingdir == 6)
                    {
                        hasslope1 = getslopecached(i + 1, j, mapvalue, slope1);
                        hasslope2 = getslopecached(i, j - 1, mapvalue, slope2);
                    }

                    if (shadingdir == 8)
                    {
                        hasslope1 = getslopecached(i + 1, j, mapvalue, slope1);
                        hasslope2 = getslopecached(i, j + 1, mapvalue, slope2);
                    }

                    if (hasslope1 && hasslope2)
                    {
                        int totalslope = (slope1 + slope2) / 10;

                        if (totalslope > 40)
                            totalslope = 40;

                        if (totalslope < -40)
                            totalslope = -40;

                        float thisshading = landshading;

                        if (truelake)
                            thisshading = lakeshading;

                        if (sea == 1)
                            thisshading = seashading;

                        totalslope = (int)((float)totalslope * (thisshading * 2.0f));

                        r = r + totalslope;
                        g = g + totalslope;
                        b = b + totalslope;
                    }

                    if (r < 0)
                        r = 0;
                    if (g < 0)
                        g = 0;
                    if (b < 0)
                        b = 0;

                    if (r > 255)
                        r = 255;
                    if (g > 255)
                        g = 255;
                    if (b > 255)
                        b = 255;
                }
            }

            setpixel(reliefpixels, regionalreliefimagewidth, i - origregwidthbegin, j - origregheightbegin, static_cast<sf::Uint8>(r), static_cast<sf::Uint8>(g), static_cast<sf::Uint8>(b));
        }
    }

    regionalreliefimage.create(regionalreliefimagewidth, regionalreliefimageheight, reliefpixels.data());
    applyregionaloutlineoverlay(world, region, regionalreliefimage);
}
