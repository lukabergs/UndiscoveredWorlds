#include "climate_circulation_diagnostics.hpp"

#include "climate_atmosphere.hpp"
#include "climate_flow.hpp"
#include "climate_grid.hpp"
#include "climate_tiff.hpp"
#include "generation_tuning.hpp"
#include "planet.hpp"

#include <SFML/Graphics/Image.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <limits>
#include <random>
#include <string>
#include <utility>
#include <vector>

namespace climatevalidation
{
namespace
{
constexpr float pi = 3.14159265358979323846f;
constexpr float secondsperday = 86400.0f;
constexpr float particletracetimestepseconds = 10800.0f;
constexpr float particlespeeddisplaymaximum = 25.0f;
constexpr float vectorerrordisplaymaximum = 25.0f;
constexpr float surfacedivergencedisplaymaximum = 1.0f;
constexpr float moistureconvergencedisplaymaximum = 20.0f;
constexpr int particlesteps = 36;
constexpr int particlerenderscale = 1;
constexpr int licstepsperdirection = 12;
constexpr float licstepcells = 0.65f;
constexpr int licrenderscale = 1;
constexpr int maxflowvisualizationcolumns = 1024;
constexpr int maxflowvisualizationcells = 1024 * 513;
constexpr int referencevisualizationversion = 2; // adaptive RK2, output-grid scale

constexpr std::array<const char*, CLIMATESEASONCOUNT> seasonnames = {
    "jan", "apr", "jul", "oct"
};

std::string diagnosticmapfilename(climatebenchmarkmapkind kind, int season)
{
    const std::string filename = climatebenchmarkmapfilename({ kind, season }, 0);
    return filename.substr(2);
}

struct seasonfields
{
    std::vector<float> surfaceu;
    std::vector<float> surfacev;
    std::vector<float> upperu;
    std::vector<float> upperv;
    std::vector<float> moisture;
    std::vector<float> surfacespeed;
    std::vector<float> upperspeed;
    std::vector<float> surfacedivergence;
    std::vector<float> upperdivergence;
    std::vector<float> moisturefluxconvergence;
    std::array<std::vector<float>, 2> moisturefluxu, moisturefluxv;
    std::vector<float> ascent, heating;
};

struct capturedwindcache
{
    const planet* source = nullptr;
    int columns = 0;
    int rows = 0;
    std::array<std::vector<float>, CLIMATESEASONCOUNT> surfaceu;
    std::array<std::vector<float>, CLIMATESEASONCOUNT> surfacev;
    std::array<std::vector<float>, CLIMATESEASONCOUNT> upperu;
    std::array<std::vector<float>, CLIMATESEASONCOUNT> upperv;
    std::array<bool, CLIMATESEASONCOUNT> populated{};

    void reset(const planet& world)
    {
        source = &world;
        columns = world.width() + 1;
        rows = world.height() + 1;
        populated.fill(false);
    }
};

capturedwindcache windcache;
const planet* weatherworld = nullptr;
int weathercolumns = 0, weatherrows = 0;
std::array<std::array<climateweather::WeatherStatistics, CLIMATESEASONCOUNT>, 2> weatherstatistics;
const planet* processworld = nullptr;
std::array<climatehydrology::SeasonalProcessFields, CLIMATESEASONCOUNT> processstatistics;
const planet* oceanworld = nullptr;
int oceancolumns = 0, oceanrows = 0;
std::array<climateocean::OceanState, CLIMATESEASONCOUNT> oceanstatistics;

enum class scalarpalette
{
    speed,
    vectorerror,
    divergence,
    convergence
};

std::size_t fieldindex(int x, int y, int columns)
{
    return static_cast<std::size_t>(y) * static_cast<std::size_t>(columns) +
        static_cast<std::size_t>(x);
}

float latitudeforrow(int y, int rows)
{
    constexpr double radiansToDegrees = 57.2957795130823208768;
    return static_cast<float>(climategrid::latitudeCentreRadians(
        y, rows, climategrid::LatitudeLayout::poleInclusive) *
        radiansToDegrees);
}

int wrapcolumn(int x, int columns)
{
    if (columns <= 0)
        return 0;

    x %= columns;
    return x < 0 ? x + columns : x;
}

bool island(planet& world, int x, int y)
{
    return world.sea(x, y) == 0 && world.truelake(x, y) == 0 &&
        world.riftlakesurface(x, y) == 0;
}

sf::Color blendcolour(const sf::Color& first, const sf::Color& second, float phase)
{
    const float amount = std::clamp(phase, 0.0f, 1.0f);
    auto blendchannel = [&](sf::Uint8 left, sf::Uint8 right)
    {
        return static_cast<sf::Uint8>(std::clamp(
            static_cast<int>(std::lround(
                static_cast<float>(left) +
                (static_cast<float>(right) - static_cast<float>(left)) * amount)),
            0,
            255));
    };

    return sf::Color(
        blendchannel(first.r, second.r),
        blendchannel(first.g, second.g),
        blendchannel(first.b, second.b));
}

sf::Color scalarcolour(float value, float displaymaximum, scalarpalette palette)
{
    if (palette == scalarpalette::speed)
    {
        const float phase = std::clamp(value / std::max(0.001f, displaymaximum), 0.0f, 1.0f);
        if (phase < 0.55f)
            return blendcolour(sf::Color(5, 12, 24), sf::Color(18, 125, 175), phase / 0.55f);

        return blendcolour(sf::Color(18, 125, 175), sf::Color(250, 235, 125),
            (phase - 0.55f) / 0.45f);
    }

    if (palette == scalarpalette::vectorerror)
    {
        const float phase = std::clamp(value /
            std::max(0.001f, displaymaximum), 0.0f, 1.0f);
        if (phase < 0.5f)
        {
            return blendcolour(
                sf::Color(12, 14, 20),
                sf::Color(130, 55, 205),
                phase * 2.0f);
        }

        return blendcolour(
            sf::Color(130, 55, 205),
            sf::Color(255, 225, 70),
            (phase - 0.5f) * 2.0f);
    }

    const float signedphase = std::clamp(
        value / std::max(0.001f, displaymaximum), -1.0f, 1.0f);
    const sf::Color neutral(18, 20, 25);

    if (palette == scalarpalette::divergence)
    {
        return signedphase < 0.0f
            ? blendcolour(neutral, sf::Color(25, 175, 245), -signedphase)
            : blendcolour(neutral, sf::Color(245, 120, 35), signedphase);
    }

    return signedphase < 0.0f
        ? blendcolour(neutral, sf::Color(225, 75, 45), -signedphase)
        : blendcolour(neutral, sf::Color(45, 220, 185), signedphase);
}

void blendpixel(sf::Image& image, int x, int y, const sf::Color& colour, float alpha)
{
    const sf::Vector2u size = image.getSize();
    if (x < 0 || y < 0 || x >= static_cast<int>(size.x) || y >= static_cast<int>(size.y))
        return;

    const sf::Color destination = image.getPixel(static_cast<unsigned int>(x), static_cast<unsigned int>(y));
    image.setPixel(
        static_cast<unsigned int>(x),
        static_cast<unsigned int>(y),
        blendcolour(destination, colour, alpha));
}

void drawline(
    sf::Image& image,
    float x0,
    float y0,
    float x1,
    float y1,
    const sf::Color& colour,
    float alpha)
{
    const float dx = x1 - x0;
    const float dy = y1 - y0;
    const int steps = std::max(1, static_cast<int>(std::ceil(std::max(std::abs(dx), std::abs(dy)))));

    for (int step = 0; step <= steps; step++)
    {
        const float phase = static_cast<float>(step) / static_cast<float>(steps);
        blendpixel(
            image,
            static_cast<int>(std::lround(x0 + dx * phase)),
            static_cast<int>(std::lround(y0 + dy * phase)),
            colour,
            alpha);
    }
}

std::pair<float, float> samplevector(
    const std::vector<float>& u,
    const std::vector<float>& v,
    float x,
    float y,
    int columns,
    int rows)
{
    const float wrappedx = std::fmod(std::fmod(x, static_cast<float>(columns)) +
        static_cast<float>(columns), static_cast<float>(columns));
    const float boundedy = std::clamp(y, 0.0f, static_cast<float>(rows - 1));
    const int x0 = wrapcolumn(static_cast<int>(std::floor(wrappedx)), columns);
    const int x1 = wrapcolumn(x0 + 1, columns);
    const int y0 = std::clamp(static_cast<int>(std::floor(boundedy)), 0, rows - 1);
    const int y1 = std::min(rows - 1, y0 + 1);
    const float xf = wrappedx - std::floor(wrappedx);
    const float yf = boundedy - std::floor(boundedy);

    auto sample = [&](const std::vector<float>& field)
    {
        const float north = field[fieldindex(x0, y0, columns)] * (1.0f - xf) +
            field[fieldindex(x1, y0, columns)] * xf;
        const float south = field[fieldindex(x0, y1, columns)] * (1.0f - xf) +
            field[fieldindex(x1, y1, columns)] * xf;
        return north * (1.0f - yf) + south * yf;
    };

    return { sample(u), sample(v) };
}

bool coastline(planet& world, int x, int y)
{
    const int columns = world.width() + 1;
    return
        island(world, x, y) != island(world, wrapcolumn(x + 1, columns), y) ||
        island(world, x, y) != island(world, wrapcolumn(x - 1, columns), y) ||
        (y > 0 && island(world, x, y) != island(world, x, y - 1)) ||
        (y < world.height() && island(world, x, y) != island(world, x, y + 1));
}

sf::Color windbackgroundcolour(planet& world, int x, int y, float speed)
{
    sf::Color background = scalarcolour(
        speed,
        particlespeeddisplaymaximum,
        scalarpalette::speed);
    background = blendcolour(
        background,
        island(world, x, y) ? sf::Color(32, 27, 23) : sf::Color(4, 13, 25),
        0.55f);

    return coastline(world, x, y)
        ? blendcolour(background, sf::Color(155, 155, 145), 0.45f)
        : background;
}

void drawdirectionarrows(
    sf::Image& image,
    const std::vector<float>& u,
    const std::vector<float>& v,
    int columns,
    int rows,
    int renderscale)
{
    constexpr int arrowcolumns = 16;
    constexpr int arrowrows = 8;
    for (int arrowy = 0; arrowy < arrowrows; arrowy++)
    {
        const float y = (static_cast<float>(arrowy) + 0.5f) *
            static_cast<float>(rows - 1) / static_cast<float>(arrowrows);

        for (int arrowx = 0; arrowx < arrowcolumns; arrowx++)
        {
            const float x = (static_cast<float>(arrowx) + 0.5f) *
                static_cast<float>(columns) / static_cast<float>(arrowcolumns);
            const auto wind = samplevector(u, v, x, y, columns, rows);
            const auto spacing = climateatmosphere::cellSpacingMetres(
                latitudeforrow(static_cast<int>(std::lround(y)), rows),
                columns,
                rows,
                tuning::climate::atmosphere::referencePlanetRadiusMetres);
            float dx = wind.first / spacing.zonalMetres;
            float dy = wind.second / spacing.meridionalMetres;
            const float magnitude = std::sqrt(dx * dx + dy * dy);

            if (magnitude <= 1.0e-9f)
                continue;

            dx /= magnitude;
            dy /= magnitude;
            const float length = 7.0f * renderscale;
            const float centrex = x * renderscale;
            const float centrey = y * renderscale;
            const float startx = centrex - dx * length * 0.5f;
            const float starty = centrey - dy * length * 0.5f;
            const float endx = centrex + dx * length * 0.5f;
            const float endy = centrey + dy * length * 0.5f;
            const sf::Color arrowcolour(255, 210, 55);
            drawline(image, startx, starty, endx, endy, arrowcolour, 0.95f);
            const float headlength = 2.5f * renderscale;
            const float perpendicularx = -dy;
            const float perpendiculary = dx;
            drawline(
                image,
                endx,
                endy,
                endx - dx * headlength + perpendicularx * headlength * 0.55f,
                endy - dy * headlength + perpendiculary * headlength * 0.55f,
                arrowcolour,
                0.95f);
            drawline(
                image,
                endx,
                endy,
                endx - dx * headlength - perpendicularx * headlength * 0.55f,
                endy - dy * headlength - perpendiculary * headlength * 0.55f,
                arrowcolour,
                0.95f);
        }
    }
}

bool renderscalarimage(
    const std::filesystem::path& path,
    const std::vector<float>& field,
    int columns,
    int rows,
    float displaymaximum,
    scalarpalette palette)
{
    sf::Image image;
    image.create(static_cast<unsigned int>(columns), static_cast<unsigned int>(rows), sf::Color::Black);

    for (int y = 0; y < rows; y++)
    {
        for (int x = 0; x < columns; x++)
        {
            image.setPixel(
                static_cast<unsigned int>(x),
                static_cast<unsigned int>(y),
                scalarcolour(field[fieldindex(x, y, columns)], displaymaximum, palette));
        }
    }

    return image.saveToFile(path.string());
}

bool rendervectorerrorimage(
    const std::filesystem::path& path,
    planet& world,
    const std::vector<float>& erroru,
    const std::vector<float>& errorv,
    const std::vector<float>& errormagnitude)
{
    const int columns = world.width() + 1;
    const int rows = world.height() + 1;
    sf::Image image;
    image.create(
        static_cast<unsigned int>(columns),
        static_cast<unsigned int>(rows),
        sf::Color::Black);

    for (int y = 0; y < rows; y++)
    {
        for (int x = 0; x < columns; x++)
        {
            sf::Color colour = scalarcolour(
                errormagnitude[fieldindex(x, y, columns)],
                vectorerrordisplaymaximum,
                scalarpalette::vectorerror);
            if (coastline(world, x, y))
                colour = blendcolour(colour, sf::Color(210, 210, 205), 0.55f);
            image.setPixel(
                static_cast<unsigned int>(x),
                static_cast<unsigned int>(y),
                colour);
        }
    }

    drawdirectionarrows(image, erroru, errorv, columns, rows, 1);
    return image.saveToFile(path.string());
}

bool renderlicimage(
    const std::filesystem::path& path,
    planet& world,
    const std::vector<float>& u,
    const std::vector<float>& v,
    std::uint32_t seed)
{
    const int columns = world.width() + 1;
    const int rows = world.height() + 1;
    const int imagecolumns = columns * licrenderscale;
    const int imagerows = rows * licrenderscale;
    const std::size_t noisecount = static_cast<std::size_t>(imagecolumns) * imagerows;
    std::mt19937 generator(seed);
    std::uniform_real_distribution<float> noisedistribution(0.0f, 1.0f);
    std::vector<float> noise(noisecount, 0.0f);
    for (float& value : noise)
        value = noisedistribution(generator);

    auto samplenoise = [&](float x, float y)
    {
        const float wrappedx = std::fmod(std::fmod(x, static_cast<float>(imagecolumns)) +
            static_cast<float>(imagecolumns), static_cast<float>(imagecolumns));
        const float boundedy = std::clamp(y, 0.0f, static_cast<float>(imagerows - 1));
        const int x0 = wrapcolumn(static_cast<int>(std::floor(wrappedx)), imagecolumns);
        const int x1 = wrapcolumn(x0 + 1, imagecolumns);
        const int y0 = std::clamp(static_cast<int>(std::floor(boundedy)), 0, imagerows - 1);
        const int y1 = std::min(imagerows - 1, y0 + 1);
        const float xf = wrappedx - std::floor(wrappedx);
        const float yf = boundedy - std::floor(boundedy);
        const auto noiseindex = [imagecolumns](int xx, int yy)
        {
            return static_cast<std::size_t>(yy) * imagecolumns + xx;
        };
        const float north = noise[noiseindex(x0, y0)] * (1.0f - xf) +
            noise[noiseindex(x1, y0)] * xf;
        const float south = noise[noiseindex(x0, y1)] * (1.0f - xf) +
            noise[noiseindex(x1, y1)] * xf;
        return north * (1.0f - yf) + south * yf;
    };

    auto streamdirection = [&](float x, float y)
    {
        const auto wind = samplevector(u, v, x, y, columns, rows);
        const float latitude = 90.0f - 180.0f * std::clamp(
            y / static_cast<float>(std::max(1, rows - 1)),
            0.0f,
            1.0f);
        const auto spacing = climateatmosphere::cellSpacingMetres(
            latitude,
            columns,
            rows,
            tuning::climate::atmosphere::referencePlanetRadiusMetres);
        float dx = wind.first / spacing.zonalMetres;
        float dy = wind.second / spacing.meridionalMetres;
        const float magnitude = std::sqrt(dx * dx + dy * dy);
        if (magnitude <= 1.0e-12f)
            return std::pair<float, float>(0.0f, 0.0f);
        return std::pair<float, float>(dx / magnitude, dy / magnitude);
    };

    sf::Image image;
    image.create(
        static_cast<unsigned int>(imagecolumns),
        static_cast<unsigned int>(imagerows),
        sf::Color::Black);
    for (int imagey = 0; imagey < imagerows; imagey++)
    {
        const float sourcey = std::clamp(
            (static_cast<float>(imagey) + 0.5f) /
                static_cast<float>(licrenderscale) - 0.5f,
            0.0f,
            static_cast<float>(rows - 1));
        for (int imagex = 0; imagex < imagecolumns; imagex++)
        {
            const float sourcex = (static_cast<float>(imagex) + 0.5f) /
                static_cast<float>(licrenderscale) - 0.5f;
            float convolution = noise[static_cast<std::size_t>(imagey) * imagecolumns + imagex];
            float weighttotal = 1.0f;

            for (int direction : { -1, 1 })
            {
                float x = sourcex;
                float y = sourcey;
                for (int step = 1; step <= licstepsperdirection; step++)
                {
                    const auto initialdirection = streamdirection(x, y);
                    if (initialdirection.first == 0.0f && initialdirection.second == 0.0f)
                        break;
                    float midpointx = x + static_cast<float>(direction) *
                        initialdirection.first * licstepcells * 0.5f;
                    const float midpointy = y + static_cast<float>(direction) *
                        initialdirection.second * licstepcells * 0.5f;
                    midpointx = std::fmod(std::fmod(midpointx, static_cast<float>(columns)) +
                        static_cast<float>(columns), static_cast<float>(columns));
                    if (midpointy < 0.0f || midpointy > static_cast<float>(rows - 1))
                        break;

                    const auto midpointdirection = streamdirection(midpointx, midpointy);
                    x += static_cast<float>(direction) * midpointdirection.first * licstepcells;
                    y += static_cast<float>(direction) * midpointdirection.second * licstepcells;
                    x = std::fmod(std::fmod(x, static_cast<float>(columns)) +
                        static_cast<float>(columns), static_cast<float>(columns));
                    if (y < 0.0f || y > static_cast<float>(rows - 1))
                        break;

                    const float weight = 0.5f + 0.5f * std::cos(
                        pi * static_cast<float>(step) /
                        static_cast<float>(licstepsperdirection + 1));
                    convolution += weight * samplenoise(
                        x * licrenderscale,
                        y * licrenderscale);
                    weighttotal += weight;
                }
            }

            const float luminance = std::clamp(
                ((convolution / weighttotal) - 0.5f) * 2.2f + 0.5f,
                0.0f,
                1.0f);
            const int cellx = wrapcolumn(
                static_cast<int>(sourcex),
                columns);
            const int celly = std::clamp(static_cast<int>(sourcey), 0, rows - 1);
            const auto wind = samplevector(u, v, sourcex, sourcey, columns, rows);
            const float speed = std::sqrt(wind.first * wind.first + wind.second * wind.second);
            sf::Color colour = windbackgroundcolour(world, cellx, celly, speed);
            colour = luminance < 0.5f
                ? blendcolour(colour, sf::Color::Black, (0.5f - luminance) * 1.4f)
                : blendcolour(colour, sf::Color::White, (luminance - 0.5f) * 1.1f);
            image.setPixel(
                static_cast<unsigned int>(imagex),
                static_cast<unsigned int>(imagey),
                colour);
        }
    }

    drawdirectionarrows(
        image,
        u,
        v,
        columns,
        rows,
        licrenderscale);
    return image.saveToFile(path.string());
}

bool renderparticleimage(
    const std::filesystem::path& path,
    planet& world,
    const std::vector<float>& u,
    const std::vector<float>& v,
    std::uint32_t seed)
{
    const int columns = world.width() + 1;
    const int rows = world.height() + 1;
    const int imagecolumns = columns * particlerenderscale;
    const int imagerows = rows * particlerenderscale;
    sf::Image image;
    image.create(
        static_cast<unsigned int>(imagecolumns),
        static_cast<unsigned int>(imagerows),
        sf::Color::Black);

    for (int y = 0; y < rows; y++)
    {
        for (int x = 0; x < columns; x++)
        {
            const std::size_t index = fieldindex(x, y, columns);
            const float speed = std::sqrt(u[index] * u[index] + v[index] * v[index]);
            const sf::Color background = windbackgroundcolour(world, x, y, speed);

            for (int py = 0; py < particlerenderscale; py++)
            {
                for (int px = 0; px < particlerenderscale; px++)
                {
                    image.setPixel(
                        static_cast<unsigned int>(x * particlerenderscale + px),
                        static_cast<unsigned int>(y * particlerenderscale + py),
                        background);
                }
            }
        }
    }

    std::mt19937 generator(seed);
    std::uniform_real_distribution<float> longitudedistribution(0.0f, static_cast<float>(columns));
    std::uniform_real_distribution<float> sinlatitudedistribution(-0.995f, 0.995f);
    const int particlecount = std::clamp(columns * rows / 10, 5000, 20000);

    for (int particle = 0; particle < particlecount; particle++)
    {
        float x = longitudedistribution(generator);
        const float latitude = std::asin(sinlatitudedistribution(generator)) * 180.0f / pi;
        float y = (90.0f - latitude) / 180.0f * static_cast<float>(rows - 1);

        for (int step = 0; step < particlesteps; step++)
        {
            climateflow::ParticleTraceConfig traceconfig;
            traceconfig.planetRadiusMetres =
                tuning::climate::atmosphere::referencePlanetRadiusMetres;
            const auto trace = climateflow::advectParticleRk2(
                columns,
                rows,
                u,
                v,
                x,
                y,
                particletracetimestepseconds,
                traceconfig);
            if (trace.points.size() < 2)
                break;
            for (std::size_t point = 1; point < trace.points.size(); point++)
            {
                const auto& previous = trace.points[point - 1];
                const auto& next = trace.points[point];
                if (std::abs(next.x - previous.x) < static_cast<float>(columns) * 0.5f)
                {
                    drawline(
                        image,
                        previous.x * particlerenderscale,
                        previous.y * particlerenderscale,
                        next.x * particlerenderscale,
                        next.y * particlerenderscale,
                        sf::Color(245, 247, 245),
                        0.16f);
                }
            }
            x = trace.points.back().x;
            y = trace.points.back().y;
            if (!trace.remainedInDomain)
                break;
        }
    }

    drawdirectionarrows(
        image,
        u,
        v,
        columns,
        rows,
        particlerenderscale);

    return image.saveToFile(path.string());
}

std::pair<float, float> finitefieldrange(const std::vector<float>& field)
{
    float minimum = std::numeric_limits<float>::infinity();
    float maximum = -std::numeric_limits<float>::infinity();

    for (float value : field)
    {
        if (!std::isfinite(value))
            continue;

        minimum = std::min(minimum, value);
        maximum = std::max(maximum, value);
    }

    if (!std::isfinite(minimum) || !std::isfinite(maximum))
        return { 0.0f, 0.0f };

    return { minimum, maximum };
}

bool writefield(
    const std::filesystem::path& outputdirectory,
    std::ofstream& manifest,
    const char* season,
    const char* fieldname,
    const char* units,
    const std::vector<float>& field,
    int columns,
    int rows,
    float displaylimit = 0.0f)
{
    const std::string filename = std::string(season) + "_" + fieldname + ".tif";
    const std::filesystem::path path = outputdirectory / filename;
    const bool written = climateio::writefloat32geotiff(
        path.string().c_str(),
        static_cast<std::uint32_t>(columns),
        static_cast<std::uint32_t>(rows),
        field.data());
    const auto range = finitefieldrange(field);
    manifest << season << ',' << fieldname << ',' << units << ','
        << range.first << ',' << range.second << ',' << displaylimit << ','
        << filename << ',' << (written ? "ok" : "error") << '\n';
    return written;
}

seasonfields buildseasonfields(planet& world, int season)
{
    const int columns = world.width() + 1;
    const int rows = world.height() + 1;
    const std::size_t cellcount = static_cast<std::size_t>(columns) * static_cast<std::size_t>(rows);
    seasonfields fields;
    fields.surfaceu.resize(cellcount);
    fields.surfacev.resize(cellcount);
    fields.upperu.resize(cellcount);
    fields.upperv.resize(cellcount);
    fields.moisture.resize(cellcount);
    fields.surfacespeed.resize(cellcount);
    fields.upperspeed.resize(cellcount);
    fields.surfacedivergence.resize(cellcount);
    fields.upperdivergence.resize(cellcount);
    fields.moisturefluxconvergence.resize(cellcount);

    for (int y = 0; y < rows; y++)
    {
        for (int x = 0; x < columns; x++)
        {
            const std::size_t index = fieldindex(x, y, columns);
            if (!capturedcirculationwind(
                    world,
                    season,
                    x,
                    y,
                    false,
                    fields.surfaceu[index],
                    fields.surfacev[index]))
            {
                fields.surfaceu[index] = static_cast<float>(world.seasonaluwind(season, x, y));
                fields.surfacev[index] = static_cast<float>(world.seasonalvwind(season, x, y));
            }

            if (!capturedcirculationwind(
                    world,
                    season,
                    x,
                    y,
                    true,
                    fields.upperu[index],
                    fields.upperv[index]))
            {
                fields.upperu[index] = static_cast<float>(world.seasonalupperuwind(season, x, y));
                fields.upperv[index] = static_cast<float>(world.seasonaluppervwind(season, x, y));
            }
            fields.moisture[index] = world.seasonalmoisture(season, x, y);
            fields.surfacespeed[index] = std::sqrt(
                fields.surfaceu[index] * fields.surfaceu[index] +
                fields.surfacev[index] * fields.surfacev[index]);
            fields.upperspeed[index] = std::sqrt(
                fields.upperu[index] * fields.upperu[index] +
                fields.upperv[index] * fields.upperv[index]);
        }
    }

    auto divergence = [&](const std::vector<float>& u,
                          const std::vector<float>& v,
                          const std::vector<float>* scalar,
                          int x,
                          int y)
    {
        const int xwest = wrapcolumn(x - 1, columns);
        const int xeast = wrapcolumn(x + 1, columns);
        const int ynorth = std::max(0, y - 1);
        const int ysouth = std::min(rows - 1, y + 1);
        const float latitude = latitudeforrow(y, rows);
        const float centrecosine = std::max(0.02f, std::abs(std::cos(latitude * pi / 180.0f)));
        const float northcosine = std::max(0.0f, std::cos(latitudeforrow(ynorth, rows) * pi / 180.0f));
        const float southcosine = std::max(0.0f, std::cos(latitudeforrow(ysouth, rows) * pi / 180.0f));
        const auto spacing = climateatmosphere::cellSpacingMetres(
            latitude,
            columns,
            rows,
            tuning::climate::atmosphere::referencePlanetRadiusMetres);

        auto weighted = [&](const std::vector<float>& field, int xx, int yy)
        {
            const std::size_t index = fieldindex(xx, yy, columns);
            return field[index] * (scalar == nullptr ? 1.0f : (*scalar)[index]);
        };

        return ((weighted(u, xeast, y) - weighted(u, xwest, y)) /
                (2.0f * spacing.zonalMetres) +
            (weighted(v, x, ysouth) * southcosine -
                weighted(v, x, ynorth) * northcosine) /
                (2.0f * spacing.meridionalMetres * centrecosine)) * secondsperday;
    };

    for (int y = 0; y < rows; y++)
    {
        for (int x = 0; x < columns; x++)
        {
            const std::size_t index = fieldindex(x, y, columns);
            fields.surfacedivergence[index] = divergence(
                fields.surfaceu, fields.surfacev, nullptr, x, y);
            fields.upperdivergence[index] = divergence(
                fields.upperu, fields.upperv, nullptr, x, y);
            fields.moisturefluxconvergence[index] = std::numeric_limits<float>::quiet_NaN();
        }
    }

    if (processworld == &world && processstatistics[season].durationSeconds > 0.0)
    {
        const auto& processes = processstatistics[season];
        const auto flux = climatehydrology::meanMoistureTransport(processes, tuning::climate::atmosphere::referencePlanetRadiusMetres);
        const auto remap = [&](const std::vector<float>& values) {
            return climategrid::remapField(processes.columns, processes.rows, climategrid::LatitudeLayout::cellCentred,
                values, columns, rows, climategrid::LatitudeLayout::poleInclusive); };
        fields.moisturefluxconvergence = remap(flux.convergenceMmPerDay);
        for (int layer = 0; layer < 2; ++layer)
        {
            fields.moisturefluxu[layer] = remap(flux.eastKgPerMetreSecond[layer]);
            fields.moisturefluxv[layer] = remap(flux.southKgPerMetreSecond[layer]);
        }
        fields.ascent = remap(processes.ascentHpaPerDay);
        std::vector<float> heat(processes.columns * processes.rows);
        for (std::size_t cell = 0; cell < heat.size(); ++cell)
            heat[cell] = processes.radiativeHeatingWm2[0][cell] + processes.radiativeHeatingWm2[1][cell] +
                processes.latentHeatingWm2[0][cell] + processes.latentHeatingWm2[1][cell] + processes.sensibleHeatingWm2[cell];
        fields.heating = remap(heat);
    }
    return fields;
}
}

void captureoceanfields(const planet& world, int season, int columns, int rows,
    const climateocean::OceanState& state)
{
    if (season < 0 || season >= CLIMATESEASONCOUNT) return;
    if (oceanworld != &world || season == 0) oceanstatistics = {};
    oceanworld = &world; oceancolumns = columns; oceanrows = rows;
    oceanstatistics[season] = state;
}

void captureclimateprocessfields(const planet& world, int season,
    const climatehydrology::SeasonalProcessFields& fields)
{
    if (season < 0 || season >= CLIMATESEASONCOUNT) return;
    if (processworld != &world || season == 0) processstatistics = {};
    processworld = &world;
    processstatistics[season] = fields;
}

void captureweatherstatistics(const planet& world, int season, int columns, int rows,
    const climateweather::WeatherStatistics& statistics, int layer)
{
    if (season < 0 || season >= CLIMATESEASONCOUNT || layer < 0 || layer > 1 || columns < 3 || rows < 2) return;
    if (weatherworld != &world || weathercolumns != columns || weatherrows != rows || (season == 0 && layer == 0))
        weatherstatistics = {};
    weatherworld = &world;
    weathercolumns = columns;
    weatherrows = rows;
    weatherstatistics[layer][season] = statistics;
}

void capturecirculationwindfields(
    planet& world,
    int season,
    const std::vector<std::vector<float>>& surfaceu,
    const std::vector<std::vector<float>>& surfacev,
    const std::vector<std::vector<float>>& upperu,
    const std::vector<std::vector<float>>& upperv)
{
    if (season < 0 || season >= CLIMATESEASONCOUNT)
        return;

    const int columns = world.width() + 1;
    const int rows = world.height() + 1;
    auto dimensionsmatch = [&](const std::vector<std::vector<float>>& field)
    {
        if (static_cast<int>(field.size()) != columns)
            return false;

        return std::all_of(field.begin(), field.end(), [&](const std::vector<float>& column)
        {
            return static_cast<int>(column.size()) == rows;
        });
    };

    if (!dimensionsmatch(surfaceu) || !dimensionsmatch(surfacev) ||
        !dimensionsmatch(upperu) || !dimensionsmatch(upperv))
    {
        return;
    }

    if (season == 0 || windcache.source != &world ||
        windcache.columns != columns || windcache.rows != rows)
    {
        windcache.reset(world);
    }

    const std::size_t cellcount = static_cast<std::size_t>(columns) * static_cast<std::size_t>(rows);
    windcache.surfaceu[season].resize(cellcount);
    windcache.surfacev[season].resize(cellcount);
    windcache.upperu[season].resize(cellcount);
    windcache.upperv[season].resize(cellcount);

    for (int y = 0; y < rows; y++)
    {
        for (int x = 0; x < columns; x++)
        {
            const std::size_t index = fieldindex(x, y, columns);
            windcache.surfaceu[season][index] = surfaceu[x][y];
            windcache.surfacev[season][index] = surfacev[x][y];
            windcache.upperu[season][index] = upperu[x][y];
            windcache.upperv[season][index] = upperv[x][y];
        }
    }

    windcache.populated[season] = true;
}

bool capturedcirculationwind(
    const planet& world,
    int season,
    int x,
    int y,
    bool upper,
    float& u,
    float& v)
{
    if (windcache.source != &world || season < 0 || season >= CLIMATESEASONCOUNT ||
        !windcache.populated[season] || x < 0 || x >= windcache.columns ||
        y < 0 || y >= windcache.rows)
    {
        return false;
    }

    const std::size_t index = fieldindex(x, y, windcache.columns);
    u = upper ? windcache.upperu[season][index] : windcache.surfaceu[season][index];
    v = upper ? windcache.upperv[season][index] : windcache.surfacev[season][index];
    return true;
}

bool exportcirculationdiagnostics(
    const std::filesystem::path& outputdirectory,
    planet& world,
    const climatebenchmarkmapselection& selection)
{
    std::error_code error;
    std::filesystem::create_directories(outputdirectory, error);
    if (error)
        return false;

    std::ofstream manifest(outputdirectory / "circulation_maps.csv");
    if (!manifest.is_open())
        return false;

    manifest << "map_id,file,status\n";
    const int columns = world.width() + 1;
    const int rows = world.height() + 1;
    bool success = true;

    if (processworld == &world)
    {
        std::ofstream processes(outputdirectory / "climate_process_fields.csv");
        processes << "season,latitude,longitude,duration_s,lower_flux_east_kg_m_s,lower_flux_south_kg_m_s,upper_flux_east_kg_m_s,upper_flux_south_kg_m_s,transport_convergence_mm_day,mean_column_water_mm,ascent_hpa_day,lower_radiative_wm2,upper_radiative_wm2,lower_latent_wm2,upper_latent_wm2,sensible_wm2,surface_net_heating_wm2,column_energy_residual_wm2\n";
        for (int season = 0; season < CLIMATESEASONCOUNT; ++season)
        {
            const auto& p = processstatistics[season];
            const auto transport = climatehydrology::meanMoistureTransport(p, tuning::climate::atmosphere::referencePlanetRadiusMetres);
            for (int y = 0; y < p.rows; ++y)
                for (int x = 0; x < p.columns; ++x)
                {
                    const int cell = y * p.columns + x;
                    processes << season << ',' << 90.0 - 180.0 * (y + 0.5) / p.rows << ',' << -180.0 + 360.0 * (x + 0.5) / p.columns << ',' << p.durationSeconds;
                    for (int layer = 0; layer < 2; ++layer)
                        processes << ',' << transport.eastKgPerMetreSecond[layer][cell] << ',' << transport.southKgPerMetreSecond[layer][cell];
                    processes << ',' << transport.convergenceMmPerDay[cell] << ',' << p.columnWaterMm[cell] << ',' << p.ascentHpaPerDay[cell];
                    for (int layer = 0; layer < 2; ++layer) processes << ',' << p.radiativeHeatingWm2[layer][cell];
                    for (int layer = 0; layer < 2; ++layer) processes << ',' << p.latentHeatingWm2[layer][cell];
                    processes << ',' << p.sensibleHeatingWm2[cell] << ',' << p.surfaceNetHeatingWm2[cell] << ',' << p.columnEnergyResidualWm2[cell] << '\n';
                }
        }
        success = processes.good() && success;
    }
    if (oceanworld == &world)
    {
        std::ofstream ocean(outputdirectory / "ocean_fields.csv");
        ocean << "season,latitude,longitude,east_current_mps,south_current_mps,sst_c,ekman_upwelling_mps,solver_converged,relative_heat_budget_residual,ice_thickness_m,surface_skin_temperature_c\n";
        for (int season = 0; season < CLIMATESEASONCOUNT; ++season)
        {
            const auto& state = oceanstatistics[season];
            for (int y = 0; y < oceanrows; ++y)
                for (int x = 0; x < oceancolumns; ++x)
                {
                    const int cell = y * oceancolumns + x;
                    ocean << season << ',' << 90.0 - 180.0 * (y + 0.5) / oceanrows << ',' << -180.0 + 360.0 * (x + 0.5) / oceancolumns << ','
                        << state.eastCurrentMps[cell] << ',' << state.southCurrentMps[cell] << ',' << state.sstC[cell] << ','
                        << (state.ekmanUpwellingMps.size() > cell ? state.ekmanUpwellingMps[cell] : 0.0f) << ','
                        << state.converged << ',' << state.relativeHeatBudgetResidual << ','
                        << state.iceThicknessMetres[cell] << ',' << state.surfaceSkinTemperatureC[cell] << '\n';
                }
        }
        success = ocean.good() && success;
    }

    if (weatherworld == &world)
    {
        std::ofstream statistics(outputdirectory / "weather_statistics.csv");
        statistics << "season,layer,latitude,longitude,samples,duration_s,effective_samples,mean_u_mps,mean_v_south_mps,mean_speed_mps,directional_consistency,speed_stddev_mps,independent_sample_stderr_mps,decorrelated_samples,correlated_stderr_mps\n";
        for (int layer = 0; layer < 2; ++layer)
        for (int season = 0; season < CLIMATESEASONCOUNT; ++season)
        {
            const auto& data = weatherstatistics[layer][season];
            if (data.meanSpeedMps.size() != static_cast<std::size_t>(weathercolumns * weatherrows)) continue;
            for (int y = 0; y < weatherrows; ++y)
                for (int x = 0; x < weathercolumns; ++x)
                {
                    const int cell = y * weathercolumns + x;
                    statistics << season << ',' << layer << ',' << 90.0 - 180.0 * (y + 0.5) / weatherrows << ','
                        << -180.0 + 360.0 * (x + 0.5) / weathercolumns << ',' << data.sampleCount << ','
                        << data.durationSeconds << ',' << data.effectiveSampleCount << ',' << data.meanEastWindMps[cell] << ','
                        << data.meanSouthWindMps[cell] << ',' << data.meanSpeedMps[cell] << ',' << data.directionalConsistency[cell] << ','
                        << data.speedStandardDeviationMps[cell] << ',' << data.speedStandardErrorMps[cell] << ','
                        << data.decorrelatedSampleCount[cell] << ',' << data.correlatedSpeedStandardErrorMps[cell] << '\n';
                }
        }
        success = statistics.good() && success;
    }

    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
    {
        const bool surfaceparticles = selection.includes(
            climatebenchmarkmapkind::surfacewindparticles, season);
        const bool upperparticles = selection.includes(
            climatebenchmarkmapkind::upperwindparticles, season);
        const bool surfacelic = selection.includes(
            climatebenchmarkmapkind::surfacewindlic, season);
        const bool upperlic = selection.includes(
            climatebenchmarkmapkind::upperwindlic, season);
        const bool surfacespeed = selection.includes(
            climatebenchmarkmapkind::surfacewindspeed, season);
        const bool surfacedivergence = selection.includes(
            climatebenchmarkmapkind::surfacedivergence, season);
        const bool moistureconvergence = selection.includes(
            climatebenchmarkmapkind::moisturefluxconvergence, season);
        const bool anyrequested = std::any_of(selection.requests().begin(), selection.requests().end(), [&](const auto& request) {
            return request.season == season && !climatebenchmarkmapisreference(request.kind); });
        if (!anyrequested) continue;
        const seasonfields fields = buildseasonfields(world, season);
        const bool flowenabled = circulationflowvisualizationenabled(world);
        const auto renderstatus = [&](climatebenchmarkmapkind kind, bool requested, bool rendered)
        {
            if (!requested)
                return;
            const climatebenchmarkmaprequest request{ kind, season };
            manifest << climatebenchmarkmapid(request) << ','
                << diagnosticmapfilename(kind, season) << ','
                << (rendered ? "written" : "skipped-resolution") << '\n';
        };

        if (surfaceparticles && flowenabled)
        {
            success = renderparticleimage(
                outputdirectory / diagnosticmapfilename(
                    climatebenchmarkmapkind::surfacewindparticles, season),
                world,
                fields.surfaceu,
                fields.surfacev,
                0x1185a11u + static_cast<std::uint32_t>(season) * 17u) && success;
        }
        renderstatus(climatebenchmarkmapkind::surfacewindparticles, surfaceparticles,
            surfaceparticles && flowenabled);
        if (upperparticles && flowenabled)
        {
            success = renderparticleimage(
                outputdirectory / diagnosticmapfilename(
                    climatebenchmarkmapkind::upperwindparticles, season),
                world,
                fields.upperu,
                fields.upperv,
                0x1185a91u + static_cast<std::uint32_t>(season) * 17u) && success;
        }
        renderstatus(climatebenchmarkmapkind::upperwindparticles, upperparticles,
            upperparticles && flowenabled);
        if (surfacelic && flowenabled)
        {
            success = renderlicimage(
                outputdirectory / diagnosticmapfilename(
                    climatebenchmarkmapkind::surfacewindlic, season),
                world,
                fields.surfaceu,
                fields.surfacev,
                0x1185b11u + static_cast<std::uint32_t>(season) * 17u) && success;
        }
        renderstatus(climatebenchmarkmapkind::surfacewindlic, surfacelic,
            surfacelic && flowenabled);
        if (upperlic && flowenabled)
        {
            success = renderlicimage(
                outputdirectory / diagnosticmapfilename(
                    climatebenchmarkmapkind::upperwindlic, season),
                world,
                fields.upperu,
                fields.upperv,
                0x1185b91u + static_cast<std::uint32_t>(season) * 17u) && success;
        }
        renderstatus(climatebenchmarkmapkind::upperwindlic, upperlic,
            upperlic && flowenabled);
        if (surfacespeed)
        {
            success = renderscalarimage(
            outputdirectory / diagnosticmapfilename(
                climatebenchmarkmapkind::surfacewindspeed, season),
            fields.surfacespeed,
            columns,
            rows,
            particlespeeddisplaymaximum,
            scalarpalette::speed) && success;
            renderstatus(climatebenchmarkmapkind::surfacewindspeed, true, true);
        }
        if (surfacedivergence)
        {
            success = renderscalarimage(
            outputdirectory / diagnosticmapfilename(
                climatebenchmarkmapkind::surfacedivergence, season),
            fields.surfacedivergence,
            columns,
            rows,
            surfacedivergencedisplaymaximum,
            scalarpalette::divergence) && success;
            renderstatus(climatebenchmarkmapkind::surfacedivergence, true, true);
        }
        if (moistureconvergence)
        {
            success = renderscalarimage(
            outputdirectory / diagnosticmapfilename(
                climatebenchmarkmapkind::moisturefluxconvergence, season),
            fields.moisturefluxconvergence,
            columns,
            rows,
            moistureconvergencedisplaymaximum,
            scalarpalette::convergence) && success;
            renderstatus(climatebenchmarkmapkind::moisturefluxconvergence, true, true);
        }
        for (const auto& request : selection.requests())
        {
            const auto kind = request.kind;
            if (request.season != season || kind < climatebenchmarkmapkind::boundarymoistureflux ||
                kind > climatebenchmarkmapkind::columnheating) continue;
            std::vector<float> values, u, v, consistency;
            float limit = 500.0f;
            auto palette = scalarpalette::speed;
            bool oceanonly = false;
            const auto remap = [&](const std::vector<float>& input, int sourcecolumns, int sourcerows) {
                return climategrid::remapField(sourcecolumns, sourcerows, climategrid::LatitudeLayout::cellCentred,
                    input, columns, rows, climategrid::LatitudeLayout::poleInclusive); };
            if (kind == climatebenchmarkmapkind::boundarymoistureflux || kind == climatebenchmarkmapkind::freemoistureflux ||
                kind == climatebenchmarkmapkind::columnmoistureflux)
            {
                const int layer = kind == climatebenchmarkmapkind::freemoistureflux ? 1 : 0;
                u = fields.moisturefluxu[layer]; v = fields.moisturefluxv[layer];
                if (kind == climatebenchmarkmapkind::columnmoistureflux && !u.empty())
                    for (std::size_t cell = 0; cell < u.size(); ++cell)
                    { u[cell] += fields.moisturefluxu[1][cell]; v[cell] += fields.moisturefluxv[1][cell]; }
            }
            else if (kind == climatebenchmarkmapkind::verticalascent || kind == climatebenchmarkmapkind::columnheating)
            {
                values = kind == climatebenchmarkmapkind::verticalascent ? fields.ascent : fields.heating;
                limit = kind == climatebenchmarkmapkind::verticalascent ? 100.0f : 300.0f;
                palette = scalarpalette::convergence;
            }
            else if (kind == climatebenchmarkmapkind::surfacewindconsistency || kind == climatebenchmarkmapkind::upperwindconsistency)
            {
                const auto& data = weatherstatistics[kind == climatebenchmarkmapkind::upperwindconsistency ? 1 : 0][season];
                if (weatherworld == &world && !data.meanSpeedMps.empty())
                {
                    values = remap(data.meanSpeedMps, weathercolumns, weatherrows);
                    consistency = remap(data.directionalConsistency, weathercolumns, weatherrows);
                    u = remap(data.meanEastWindMps, weathercolumns, weatherrows);
                    v = remap(data.meanSouthWindMps, weathercolumns, weatherrows);
                }
                limit = 25.0f;
            }
            else if (oceanworld == &world)
            {
                oceanonly = true;
                const auto& data = oceanstatistics[season];
                if (kind == climatebenchmarkmapkind::seasurfacetemperature)
                { values = remap(data.sstC, oceancolumns, oceanrows); limit = 35.0f; palette = scalarpalette::divergence; }
                else
                { u = remap(data.eastCurrentMps, oceancolumns, oceanrows); v = remap(data.southCurrentMps, oceancolumns, oceanrows); limit = 2.0f; }
            }
            if (values.empty() && !u.empty())
            {
                values.resize(u.size());
                for (std::size_t cell = 0; cell < u.size(); ++cell) values[cell] = std::hypot(u[cell], v[cell]);
            }
            if (values.size() != static_cast<std::size_t>(columns * rows))
            {
                manifest << climatebenchmarkmapid(request) << ",,unavailable-process-fields\n";
                success = false;
                continue;
            }
            sf::Image image;
            image.create(columns, rows, sf::Color::Black);
            for (int y = 0; y < rows; ++y)
                for (int x = 0; x < columns; ++x)
                {
                    const auto cell = fieldindex(x, y, columns);
                    if (oceanonly && !world.sea(x, y))
                    { if (!u.empty()) { u[cell] = 0.0f; v[cell] = 0.0f; } continue; }
                    auto colour = scalarcolour(values[cell], limit, palette);
                    if (!consistency.empty())
                    {
                        colour = blendcolour(sf::Color(30, 210, 240), sf::Color(240, 40, 25), consistency[cell]);
                        colour = blendcolour(sf::Color(8, 10, 15), colour, 0.2f + 0.8f * std::clamp(values[cell] / limit, 0.0f, 1.0f));
                    }
                    image.setPixel(x, y, colour);
                }
            if (!u.empty()) drawdirectionarrows(image, u, v, columns, rows, 1);
            const bool written = image.saveToFile((outputdirectory / diagnosticmapfilename(kind, season)).string());
            success = written && success;
            renderstatus(kind, true, written);
        }
    }

    return success;
}

bool exportcirculationreferencecomparisons(
    const std::filesystem::path& diagnosticoutputdirectory,
    const std::filesystem::path& referencecacheoutputdirectory,
    planet& world,
    const circulationreferencewindfields& referencefields,
    const climatebenchmarkmapselection& selection)
{
    if (!circulationflowvisualizationenabled(world))
        return true;

    const int columns = world.width() + 1;
    const int rows = world.height() + 1;
    const std::size_t cellcount = static_cast<std::size_t>(columns) * rows;
    if (referencefields.columns != columns || referencefields.rows != rows)
        return false;

    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
    {
        if (referencefields.surfaceu[season].size() != cellcount ||
            referencefields.surfacev[season].size() != cellcount ||
            referencefields.upperu[season].size() != cellcount ||
            referencefields.upperv[season].size() != cellcount)
        {
            return false;
        }
    }

    std::error_code error;
    std::filesystem::create_directories(diagnosticoutputdirectory, error);
    if (error)
        return false;
    std::filesystem::create_directories(referencecacheoutputdirectory, error);
    if (error)
        return false;

    const std::string referenceprefix =
        "era5_v" + std::to_string(referencevisualizationversion) + "_" +
        std::to_string(columns) + "x" + std::to_string(rows);
    std::ofstream comparison(
        diagnosticoutputdirectory / "wind_vector_comparison.csv");
    if (!comparison.is_open())
        return false;
    comparison
        << "season,layer,compared_cells,area_weighted_mean_vector_error_m_s,"
        << "area_weighted_vector_rmse_m_s,area_weighted_speed_bias_m_s,"
        << "area_weighted_mean_direction_error_degrees\n"
        << std::fixed << std::setprecision(7);

    bool success = true;
    bool renderedreference = false;
    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
    {
        const std::string seasonname = seasonnames[season];
        const seasonfields simulated = buildseasonfields(world, season);
        for (int layer = 0; layer < 2; layer++)
        {
            const bool upper = layer == 1;
            const char* layername = upper ? "upper" : "surface";
            const climatebenchmarkmapkind referenceparticlekind = upper
                ? climatebenchmarkmapkind::era5upperwindparticles
                : climatebenchmarkmapkind::era5surfacewindparticles;
            const climatebenchmarkmapkind referencelickind = upper
                ? climatebenchmarkmapkind::era5upperwindlic
                : climatebenchmarkmapkind::era5surfacewindlic;
            const climatebenchmarkmapkind errorkind = upper
                ? climatebenchmarkmapkind::upperwindvectorerror
                : climatebenchmarkmapkind::surfacewindvectorerror;
            const bool referenceparticlesrequested =
                selection.includes(referenceparticlekind, season);
            const bool referencelicrequested = selection.includes(referencelickind, season);
            const bool errorrequested = selection.includes(errorkind, season);
            const std::vector<float>& simulatedu = upper
                ? simulated.upperu : simulated.surfaceu;
            const std::vector<float>& simulatedv = upper
                ? simulated.upperv : simulated.surfacev;
            const std::vector<float>& referenceu = upper
                ? referencefields.upperu[season] : referencefields.surfaceu[season];
            const std::vector<float>& referencev = upper
                ? referencefields.upperv[season] : referencefields.surfacev[season];
            const std::filesystem::path referenceparticlepath =
                referencecacheoutputdirectory /
                (referenceprefix + "_" + seasonname + "_" +
                    (upper ? "u" : "s") + "_wind_part.png");
            const std::filesystem::path referencelicpath =
                referencecacheoutputdirectory /
                (referenceprefix + "_" + seasonname + "_" +
                    (upper ? "u" : "s") + "_wind_lic.png");
            const std::uint32_t particlebase = upper ? 0x1185a91u : 0x1185a11u;
            const std::uint32_t licbase = upper ? 0x1185b91u : 0x1185b11u;

            if (referenceparticlesrequested &&
                !std::filesystem::exists(referenceparticlepath))
            {
                success = renderparticleimage(
                    referenceparticlepath,
                    world,
                    referenceu,
                    referencev,
                    particlebase + static_cast<std::uint32_t>(season) * 17u) && success;
                renderedreference = true;
            }
            if (referencelicrequested && !std::filesystem::exists(referencelicpath))
            {
                success = renderlicimage(
                    referencelicpath,
                    world,
                    referenceu,
                    referencev,
                    licbase + static_cast<std::uint32_t>(season) * 17u) && success;
                renderedreference = true;
            }

            std::vector<float> erroru(cellcount, 0.0f);
            std::vector<float> errorv(cellcount, 0.0f);
            std::vector<float> errormagnitude(cellcount, 0.0f);
            double weighttotal = 0.0;
            double directionweighttotal = 0.0;
            double weightedabsoluteerror = 0.0;
            double weightedsquarederror = 0.0;
            double weightedspeedbias = 0.0;
            double weighteddirectionerror = 0.0;
            std::size_t comparedcells = 0;
            for (int y = 0; y < rows; y++)
            {
                const double areaweight = std::max(
                    0.0,
                    std::cos(static_cast<double>(latitudeforrow(y, rows)) * pi / 180.0));
                for (int x = 0; x < columns; x++)
                {
                    const std::size_t index = fieldindex(x, y, columns);
                    const float su = simulatedu[index];
                    const float sv = simulatedv[index];
                    const float ru = referenceu[index];
                    const float rv = referencev[index];
                    if (!std::isfinite(su) || !std::isfinite(sv) ||
                        !std::isfinite(ru) || !std::isfinite(rv))
                    {
                        continue;
                    }

                    erroru[index] = su - ru;
                    errorv[index] = sv - rv;
                    const double magnitude = std::hypot(
                        static_cast<double>(erroru[index]),
                        static_cast<double>(errorv[index]));
                    errormagnitude[index] = static_cast<float>(magnitude);
                    const double simulatedspeed = std::hypot(
                        static_cast<double>(su), static_cast<double>(sv));
                    const double referencespeed = std::hypot(
                        static_cast<double>(ru), static_cast<double>(rv));
                    weighttotal += areaweight;
                    weightedabsoluteerror += areaweight * magnitude;
                    weightedsquarederror += areaweight * magnitude * magnitude;
                    weightedspeedbias += areaweight * (simulatedspeed - referencespeed);
                    comparedcells++;

                    if (simulatedspeed > 1.0e-6 && referencespeed > 1.0e-6)
                    {
                        const double cosine = std::clamp(
                            (static_cast<double>(su) * ru + static_cast<double>(sv) * rv) /
                                (simulatedspeed * referencespeed),
                            -1.0,
                            1.0);
                        weighteddirectionerror += areaweight * std::acos(cosine) *
                            180.0 / static_cast<double>(pi);
                        directionweighttotal += areaweight;
                    }
                }
            }

            if (errorrequested)
            {
                success = rendervectorerrorimage(
                    diagnosticoutputdirectory / diagnosticmapfilename(errorkind, season),
                    world,
                    erroru,
                    errorv,
                    errormagnitude) && success;
            }

            comparison
                << seasonname << ',' << layername << ',' << comparedcells << ','
                << (weighttotal > 0.0 ? weightedabsoluteerror / weighttotal : 0.0) << ','
                << (weighttotal > 0.0
                    ? std::sqrt(weightedsquarederror / weighttotal) : 0.0) << ','
                << (weighttotal > 0.0 ? weightedspeedbias / weighttotal : 0.0) << ','
                << (directionweighttotal > 0.0
                    ? weighteddirectionerror / directionweighttotal : 0.0) << '\n';
        }
    }

    std::ofstream cachestatus(
        diagnosticoutputdirectory / "wind_reference_visualization_cache.txt");
    if (!cachestatus.is_open())
        return false;
    cachestatus << "reference_prefix=" << referenceprefix << '\n';
    cachestatus << "renderer_version=" << referencevisualizationversion << '\n';
    cachestatus << "reference_images_rendered=" << (renderedreference ? 1 : 0) << '\n';
    cachestatus << "reference_images_reused=" << (renderedreference ? 0 : 1) << '\n';
    cachestatus << "vector_error=simulated_minus_era5\n";
    cachestatus << "vector_error_display_range_m_s=0_to_" << vectorerrordisplaymaximum << '\n';
    return success;
}

bool circulationflowvisualizationenabled(const planet& world)
{
    const int columns = world.width() + 1;
    const int rows = world.height() + 1;
    return columns <= maxflowvisualizationcolumns &&
        static_cast<long long>(columns) * rows <= maxflowvisualizationcells;
}

int circulationflowvisualizationmaxcolumns()
{
    return maxflowvisualizationcolumns;
}

int circulationflowvisualizationmaxcells()
{
    return maxflowvisualizationcells;
}
}
