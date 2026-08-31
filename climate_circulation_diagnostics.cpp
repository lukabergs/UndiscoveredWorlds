#include "climate_circulation_diagnostics.hpp"

#include "climate_atmosphere.hpp"
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
constexpr float surfacedivergencedisplaymaximum = 1.0f;
constexpr float moistureconvergencedisplaymaximum = 20.0f;
constexpr int particlesteps = 36;
constexpr int particlerenderscale = 2;

constexpr std::array<const char*, CLIMATESEASONCOUNT> seasonnames = {
    "january", "april", "july", "october"
};

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

enum class scalarpalette
{
    speed,
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
    if (rows <= 1)
        return 0.0f;

    return 90.0f - 180.0f * static_cast<float>(y) /
        static_cast<float>(rows - 1);
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
            sf::Color background = scalarcolour(
                speed,
                particlespeeddisplaymaximum,
                scalarpalette::speed);
            background = blendcolour(
                background,
                island(world, x, y) ? sf::Color(32, 27, 23) : sf::Color(4, 13, 25),
                0.55f);

            const bool coastline = island(world, x, y) != island(world, wrapcolumn(x + 1, columns), y) ||
                island(world, x, y) != island(world, wrapcolumn(x - 1, columns), y) ||
                (y > 0 && island(world, x, y) != island(world, x, y - 1)) ||
                (y + 1 < rows && island(world, x, y) != island(world, x, y + 1));
            if (coastline)
                background = blendcolour(background, sf::Color(155, 155, 145), 0.45f);

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
            const auto wind = samplevector(u, v, x, y, columns, rows);
            const float localLatitude = latitudeforrow(
                std::clamp(static_cast<int>(std::lround(y)), 0, rows - 1),
                rows);
            const auto spacing = climateatmosphere::cellSpacingMetres(
                localLatitude,
                columns,
                rows,
                tuning::climate::atmosphere::referencePlanetRadiusMetres);
            float dx = wind.first * particletracetimestepseconds / spacing.zonalMetres;
            float dy = wind.second * particletracetimestepseconds / spacing.meridionalMetres;
            const float maximumstep = std::max(std::abs(dx), std::abs(dy));

            if (maximumstep < 0.01f)
                break;

            if (maximumstep > 1.75f)
            {
                dx *= 1.75f / maximumstep;
                dy *= 1.75f / maximumstep;
            }

            float nextx = x + dx;
            const float nexty = y + dy;
            nextx = std::fmod(std::fmod(nextx, static_cast<float>(columns)) +
                static_cast<float>(columns), static_cast<float>(columns));

            if (nexty < 0.0f || nexty > static_cast<float>(rows - 1))
                break;

            if (std::abs(nextx - x) < static_cast<float>(columns) * 0.5f)
            {
                drawline(
                    image,
                    x * particlerenderscale,
                    y * particlerenderscale,
                    nextx * particlerenderscale,
                    nexty * particlerenderscale,
                    sf::Color(245, 247, 245),
                    0.16f);
            }

            x = nextx;
            y = nexty;
        }
    }

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
            const float length = 7.0f * particlerenderscale;
            const float centrex = x * particlerenderscale;
            const float centrey = y * particlerenderscale;
            const float startx = centrex - dx * length * 0.5f;
            const float starty = centrey - dy * length * 0.5f;
            const float endx = centrex + dx * length * 0.5f;
            const float endy = centrey + dy * length * 0.5f;
            const sf::Color arrowcolour(255, 210, 55);
            drawline(image, startx, starty, endx, endy, arrowcolour, 0.95f);
            const float headlength = 2.5f * particlerenderscale;
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
            fields.moisturefluxconvergence[index] = -divergence(
                fields.surfaceu, fields.surfacev, &fields.moisture, x, y);
        }
    }

    return fields;
}
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
    planet& world)
{
    std::error_code error;
    std::filesystem::create_directories(outputdirectory, error);
    if (error)
        return false;

    std::ofstream manifest(outputdirectory / "circulation_diagnostics.csv");
    if (!manifest.is_open())
        return false;

    manifest << "season,field,units,minimum,maximum,display_limit,file,status\n";
    manifest << std::fixed << std::setprecision(7);
    const int columns = world.width() + 1;
    const int rows = world.height() + 1;
    bool success = true;

    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
    {
        const char* seasonname = seasonnames[season];
        const seasonfields fields = buildseasonfields(world, season);
        std::vector<float> surfacenorthward = fields.surfacev;
        std::vector<float> uppernorthward = fields.upperv;
        for (float& value : surfacenorthward)
            value = -value;
        for (float& value : uppernorthward)
            value = -value;

        success = writefield(outputdirectory, manifest, seasonname,
            "surface_u_eastward_m_s", "m s-1", fields.surfaceu, columns, rows) && success;
        success = writefield(outputdirectory, manifest, seasonname,
            "surface_v_northward_m_s", "m s-1", surfacenorthward, columns, rows) && success;
        success = writefield(outputdirectory, manifest, seasonname,
            "upper_u_eastward_m_s", "m s-1", fields.upperu, columns, rows) && success;
        success = writefield(outputdirectory, manifest, seasonname,
            "upper_v_northward_m_s", "m s-1", uppernorthward, columns, rows) && success;
        success = writefield(outputdirectory, manifest, seasonname,
            "column_water_mm", "kg m-2", fields.moisture, columns, rows) && success;
        success = writefield(outputdirectory, manifest, seasonname,
            "surface_wind_speed_m_s", "m s-1", fields.surfacespeed,
            columns, rows, particlespeeddisplaymaximum) && success;
        success = writefield(outputdirectory, manifest, seasonname,
            "upper_wind_speed_m_s", "m s-1", fields.upperspeed,
            columns, rows, particlespeeddisplaymaximum) && success;
        success = writefield(outputdirectory, manifest, seasonname,
            "surface_divergence_per_day", "day-1", fields.surfacedivergence,
            columns, rows, surfacedivergencedisplaymaximum) && success;
        success = writefield(outputdirectory, manifest, seasonname,
            "upper_divergence_per_day", "day-1", fields.upperdivergence,
            columns, rows, surfacedivergencedisplaymaximum) && success;
        success = writefield(outputdirectory, manifest, seasonname,
            "surface_moisture_flux_convergence_mm_day", "kg m-2 day-1",
            fields.moisturefluxconvergence, columns, rows,
            moistureconvergencedisplaymaximum) && success;

        success = renderparticleimage(
            outputdirectory / (std::string(seasonname) + "_surface_wind_particles.png"),
            world,
            fields.surfaceu,
            fields.surfacev,
            0x1185a11u + static_cast<std::uint32_t>(season) * 17u) && success;
        success = renderparticleimage(
            outputdirectory / (std::string(seasonname) + "_upper_wind_particles.png"),
            world,
            fields.upperu,
            fields.upperv,
            0x1185a91u + static_cast<std::uint32_t>(season) * 17u) && success;
        success = renderscalarimage(
            outputdirectory / (std::string(seasonname) + "_surface_wind_speed.png"),
            fields.surfacespeed,
            columns,
            rows,
            particlespeeddisplaymaximum,
            scalarpalette::speed) && success;
        success = renderscalarimage(
            outputdirectory / (std::string(seasonname) + "_surface_divergence.png"),
            fields.surfacedivergence,
            columns,
            rows,
            surfacedivergencedisplaymaximum,
            scalarpalette::divergence) && success;
        success = renderscalarimage(
            outputdirectory / (std::string(seasonname) + "_moisture_flux_convergence.png"),
            fields.moisturefluxconvergence,
            columns,
            rows,
            moistureconvergencedisplaymaximum,
            scalarpalette::convergence) && success;
    }

    return success;
}
}
