#include "app_environment.hpp"
#include "climate_atmosphere.hpp"
#include "climate_energy.hpp"
#include "climate_koppen.hpp"
#include "climate_physics.hpp"
#include "climate_reference.hpp"
#include "climate_tiff.hpp"
#include "functions.hpp"
#include "generation_tuning.hpp"

#include <windows.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

namespace
{
struct zonalstats
{
    int cells = 0;
    int landcells = 0;
    int oceancells = 0;
    double annualrain = 0.0;
    double januaryrain = 0.0;
    double julyrain = 0.0;
    double landannualrain = 0.0;
    double oceanannualrain = 0.0;
    double januarypressure = 0.0;
    double julypressure = 0.0;
    double januaryuwind = 0.0;
    double januaryvwind = 0.0;
    double julyuwind = 0.0;
    double julyvwind = 0.0;
    double januarysst = 0.0;
    double julysst = 0.0;
    double januarycurrentu = 0.0;
    double januarycurrentv = 0.0;
    double julycurrentu = 0.0;
    double julycurrentv = 0.0;
    double januaryevaporation = 0.0;
    double julyevaporation = 0.0;
    double januarymoisture = 0.0;
    double julymoisture = 0.0;
};

struct comparisonmetrics
{
    bool referencefound = false;
    bool dimensionsmatch = false;
    int comparedcells = 0;
    double simulatedmean = 0.0;
    double referencemean = 0.0;
    double meanbias = 0.0;
    double meanabsoluteerror = 0.0;
    double rmse = 0.0;
    double correlation = 0.0;
    double tropicalmeanbias = 0.0;
    double areaweight = 0.0;
    double areaweightedsimulatedmean = 0.0;
    double areaweightedreferencemean = 0.0;
    double areaweightedmeanbias = 0.0;
    double areaweightedmeanabsoluteerror = 0.0;
    double areaweightedrmse = 0.0;
    double areaweightedcorrelation = 0.0;
};

struct climatespatialmetrics
{
    bool referencefound = false;
    bool dimensionsmatch = false;
    long long comparedcells = 0;
    long long exactmatches = 0;
    long long groupmatches = 0;
    double exactaccuracy = 0.0;
    double groupaccuracy = 0.0;
    double kappa = 0.0;
    double areaweight = 0.0;
    double areaweightedexactaccuracy = 0.0;
    double areaweightedgroupaccuracy = 0.0;
    double areaweightedkappa = 0.0;
};

struct benchmarkphysicsmetrics
{
    double areaweightedlandrawzeroprecipitationfraction = 0.0;
    double areaweightedlandstoredzeroprecipitationfraction = 0.0;
    double areaweightedlandbelowonemillimetrefraction = 0.0;
    double areaweightedlandevaporationoceanevaporationratio = 0.0;
    double areaweightedrunofflandprecipitationratio = 0.0;
    double areaweightedlandprecipitationrecyclingratio = 0.0;
    double areaweightedlandprecipitationcorrelation = 0.0;
    double areaweightedlandprecipitationcorrelation5degree = 0.0;
    double areaweightedlandrainfallstorageclippingfraction = 0.0;
    double areaweightedatmosphericwaterresidencetimedays = 0.0;
    double areaweightedgriddedatmosphericwaterresidencetimedays = 0.0;
    double areaweightedactivecondensationcellstepfraction = 0.0;
    double areaweightedactivecondensationwaterfraction = 0.0;
    double areaweightedcondensableexcesswaterfraction = 0.0;
    double areaweightedtemperaturegrouplockedfraction = 0.0;
    double areaweightedthermalsignaturelockedfraction = 0.0;
    double areaweightednorthern5070temperaturegrouplockedfraction = 0.0;
    double areaweightednorthern5070thermalsignaturelockedfraction = 0.0;
    double areaweightednorthern5070simulatedwarmesttemperaturec = 0.0;
    double areaweightednorthern5070referencewarmesttemperaturec = 0.0;
    double areaweightednorthern5070precipitationratio = 0.0;
    double areaweightedsouthern6090simulatedeffraction = 0.0;
    double areaweightedsouthern6090referenceeffraction = 0.0;
    double areaweightedsouthern6090simulateddfraction = 0.0;
    double areaweightedsouthern6090referencedfraction = 0.0;
    double areaweightedsouthern6090simulatedwarmesttemperaturec = 0.0;
    double areaweightedsouthern6090referencewarmesttemperaturec = 0.0;
    double era5globalpressurecorrelation = 0.0;
    double era5globalpressurermsehpa = 0.0;
    double era5surfaceeastwardwindcorrelation = 0.0;
    double era5surfacenorthwardwindcorrelation = 0.0;
    double era5transporteastwardwindcorrelation = 0.0;
    double era5transportnorthwardwindcorrelation = 0.0;
    double era5columnwatercorrelation = 0.0;
    double era5columnwatermeanbiaskgm2 = 0.0;
    double era5verticalascentcorrelation = 0.0;
    double era5verticalascentrmsehpaday = 0.0;
    double era5landprecipitationcorrelation = 0.0;
};

struct climatedriverstats
{
    long long cells = 0;
    long long bwcells = 0;
    long long aridcells = 0;
    long long dthermalcells = 0;
    long long polarcells = 0;
    long long winterdrycells = 0;
    double meanannualtemp = 0.0;
    double mintemp = 0.0;
    double maxtemp = 0.0;
    double annualrain = 0.0;
    double drythreshold = 0.0;
    double raintothreshold = 0.0;
    double warmhalffraction = 0.0;
    double driestcoldrain = 0.0;
    double wettestwarmrain = 0.0;
};

constexpr int benchmarkcolourdistancelimit = 65;
constexpr int benchmarkcolourdistancelimitsquared = benchmarkcolourdistancelimit * benchmarkcolourdistancelimit;

constexpr array<array<int, 3>, 31> benchmarkclimatecolours =
{
    array<int, 3>{ 0, 0, 254 }, { 1, 119, 255 }, { 70, 169, 250 }, { 70, 169, 250 },
    { 249, 15, 0 }, { 251, 150, 149 }, { 245, 163, 1 }, { 254, 219, 99 },
    { 255, 255, 0 }, { 198, 199, 1 }, { 184, 184, 114 }, { 138, 255, 162 },
    { 86, 199, 112 }, { 30, 150, 66 }, { 192, 254, 109 }, { 76, 255, 93 },
    { 19, 203, 74 }, { 255, 8, 245 }, { 204, 3, 192 }, { 154, 51, 144 },
    { 153, 100, 146 }, { 172, 178, 249 }, { 91, 121, 213 }, { 78, 83, 175 },
    { 54, 3, 130 }, { 0, 255, 245 }, { 32, 200, 250 }, { 0, 126, 125 },
    { 0, 69, 92 }, { 178, 178, 178 }, { 104, 104, 104 }
};

sf::Color benchmarkclimatecolour(short climate)
{
    if (climate < 1 || climate > static_cast<short>(benchmarkclimatecolours.size()))
        return sf::Color::Black;

    const auto& colour = benchmarkclimatecolours[climate - 1];
    return sf::Color(colour[0], colour[1], colour[2]);
}

struct benchmarkcolouranchor
{
    float value = 0.0f;
    sf::Color colour;
};

template <size_t AnchorCount>
sf::Color benchmarkscalarcolour(
    float value,
    const array<benchmarkcolouranchor, AnchorCount>& anchors)
{
    if (value <= anchors.front().value)
        return anchors.front().colour;

    if (value >= anchors.back().value)
        return anchors.back().colour;

    for (size_t index = 1; index < anchors.size(); index++)
    {
        if (value > anchors[index].value)
            continue;

        const benchmarkcolouranchor& lower = anchors[index - 1];
        const benchmarkcolouranchor& upper = anchors[index];
        const float fraction = (value - lower.value) / (upper.value - lower.value);
        const auto interpolate = [&](sf::Uint8 first, sf::Uint8 second)
        {
            return static_cast<sf::Uint8>(std::clamp(
                static_cast<int>(std::round(
                    static_cast<float>(first) + fraction *
                    (static_cast<float>(second) - static_cast<float>(first)))),
                0,
                255));
        };

        return sf::Color(
            interpolate(lower.colour.r, upper.colour.r),
            interpolate(lower.colour.g, upper.colour.g),
            interpolate(lower.colour.b, upper.colour.b));
    }

    return anchors.back().colour;
}

sf::Color benchmarktemperaturecolour(float temperaturec)
{
    static const array<benchmarkcolouranchor, 8> anchors =
    {
        benchmarkcolouranchor{ -54.0f, sf::Color(225, 225, 225) },
        { -40.0f, sf::Color(185, 0, 255) },
        { -25.0f, sf::Color(45, 0, 255) },
        { -10.0f, sf::Color(0, 105, 255) },
        { 0.0f, sf::Color(0, 220, 230) },
        { 10.0f, sf::Color(0, 205, 70) },
        { 20.0f, sf::Color(255, 230, 0) },
        { 30.0f, sf::Color(220, 45, 0) }
    };

    return benchmarkscalarcolour(temperaturec, anchors);
}

sf::Color benchmarkprecipitationcolour(float precipitationmmyear)
{
    static const array<benchmarkcolouranchor, 9> anchors =
    {
        benchmarkcolouranchor{ 0.0f, sf::Color(4, 10, 35) },
        { 250.0f, sf::Color(12, 35, 95) },
        { 500.0f, sf::Color(15, 80, 160) },
        { 1000.0f, sf::Color(10, 150, 205) },
        { 1500.0f, sf::Color(15, 175, 90) },
        { 2500.0f, sf::Color(135, 205, 35) },
        { 3500.0f, sf::Color(245, 225, 25) },
        { 4500.0f, sf::Color(250, 115, 20) },
        { 6000.0f, sf::Color(220, 25, 20) }
    };

    return benchmarkscalarcolour(precipitationmmyear, anchors);
}

short nearestbenchmarkclimate(const sf::Color& pixel, int& distancesquared)
{
    short nearest = 0;
    distancesquared = 3 * 255 * 255;

    for (short climate = 1; climate <= static_cast<short>(benchmarkclimatecolours.size()); climate++)
    {
        const auto& colour = benchmarkclimatecolours[climate - 1];
        const int red = static_cast<int>(pixel.r) - colour[0];
        const int green = static_cast<int>(pixel.g) - colour[1];
        const int blue = static_cast<int>(pixel.b) - colour[2];
        const int candidate = red * red + green * green + blue * blue;

        if (candidate < distancesquared)
        {
            nearest = climate;
            distancesquared = candidate;
        }
    }

    return nearest;
}

short comparableclimate(short climate)
{
    // Aw and As share a reference-map colour, so the image cannot distinguish them.
    return climate == 4 ? 3 : climate;
}

int climatemajorgroup(short climate)
{
    if (climate >= 1 && climate <= 4) return 0;
    if (climate >= 5 && climate <= 8) return 1;
    if (climate >= 9 && climate <= 17) return 2;
    if (climate >= 18 && climate <= 29) return 3;
    if (climate >= 30 && climate <= 31) return 4;
    return -1;
}

double gridcellareaweight(int row, int height)
{
    if (height <= 0)
        return 0.0;

    constexpr double pi = 3.14159265358979323846;
    const double latitude = 90.0 - 180.0 * static_cast<double>(row) / static_cast<double>(height);
    return max(0.0, cos(latitude * pi / 180.0));
}

bool isvalidationland(planet& world, int x, int y);

climatespatialmetrics compareclimatespatially(planet& world)
{
    climatespatialmetrics metrics;
    sf::Image reference;
    const filesystem::path referencepath = getappenvironment().earthKoppenImagePath;

    if (reference.loadFromFile(referencepath.string()) == false)
        return metrics;

    metrics.referencefound = true;
    const int width = world.width();
    const int height = world.height();
    const sf::Vector2u referencesize = reference.getSize();

    if (referencesize.x != static_cast<unsigned int>(width + 1) || referencesize.y != static_cast<unsigned int>(height + 1))
        return metrics;

    metrics.dimensionsmatch = true;
    array<array<long long, 31>, 31> confusion{};
    array<array<double, 31>, 31> areaweightedconfusion{};
    double areaweightedexactmatches = 0.0;
    double areaweightedgroupmatches = 0.0;

    for (int y = 0; y <= height; y++)
    {
        const double areaweight = gridcellareaweight(y, height);

        for (int x = 0; x <= width; x++)
        {
            if (!isvalidationland(world, x, y))
                continue;

            short simulated = comparableclimate(static_cast<short>(world.climate(x, y)));

            if (simulated < 1 || simulated > 31)
                continue;

            int distancesquared = 0;
            short expected = nearestbenchmarkclimate(reference.getPixel(x, y), distancesquared);

            if (distancesquared > benchmarkcolourdistancelimitsquared)
                continue;

            expected = comparableclimate(expected);
            metrics.comparedcells++;
            confusion[simulated - 1][expected - 1]++;
            metrics.areaweight += areaweight;
            areaweightedconfusion[simulated - 1][expected - 1] += areaweight;

            if (simulated == expected)
            {
                metrics.exactmatches++;
                areaweightedexactmatches += areaweight;
            }

            if (climatemajorgroup(simulated) == climatemajorgroup(expected))
            {
                metrics.groupmatches++;
                areaweightedgroupmatches += areaweight;
            }
        }
    }

    if (metrics.comparedcells <= 0)
        return metrics;

    metrics.exactaccuracy = static_cast<double>(metrics.exactmatches) / static_cast<double>(metrics.comparedcells);
    metrics.groupaccuracy = static_cast<double>(metrics.groupmatches) / static_cast<double>(metrics.comparedcells);
    metrics.areaweightedexactaccuracy = areaweightedexactmatches / metrics.areaweight;
    metrics.areaweightedgroupaccuracy = areaweightedgroupmatches / metrics.areaweight;

    array<long long, 31> simulatedtotals{};
    array<long long, 31> referencetotals{};

    for (int simulated = 0; simulated < 31; simulated++)
    {
        for (int expected = 0; expected < 31; expected++)
        {
            simulatedtotals[simulated] += confusion[simulated][expected];
            referencetotals[expected] += confusion[simulated][expected];
        }
    }

    double chanceagreement = 0.0;
    const double comparedsquared = static_cast<double>(metrics.comparedcells) * static_cast<double>(metrics.comparedcells);

    for (int climate = 0; climate < 31; climate++)
        chanceagreement += static_cast<double>(simulatedtotals[climate]) * static_cast<double>(referencetotals[climate]) / comparedsquared;

    if (chanceagreement < 1.0)
        metrics.kappa = (metrics.exactaccuracy - chanceagreement) / (1.0 - chanceagreement);

    array<double, 31> areaweightedsimulatedtotals{};
    array<double, 31> areaweightedreferencetotals{};

    for (int simulated = 0; simulated < 31; simulated++)
    {
        for (int expected = 0; expected < 31; expected++)
        {
            areaweightedsimulatedtotals[simulated] += areaweightedconfusion[simulated][expected];
            areaweightedreferencetotals[expected] += areaweightedconfusion[simulated][expected];
        }
    }

    double areaweightedchanceagreement = 0.0;

    for (int climate = 0; climate < 31; climate++)
    {
        areaweightedchanceagreement +=
            areaweightedsimulatedtotals[climate] * areaweightedreferencetotals[climate] /
            (metrics.areaweight * metrics.areaweight);
    }

    if (areaweightedchanceagreement < 1.0)
    {
        metrics.areaweightedkappa =
            (metrics.areaweightedexactaccuracy - areaweightedchanceagreement) /
            (1.0 - areaweightedchanceagreement);
    }

    return metrics;
}

double safeaverage(double total, int count)
{
    if (count <= 0)
        return 0.0;

    return total / static_cast<double>(count);
}

string csvescape(const string& value);

enum class monthlyreferencefield
{
    temperature,
    precipitation,
    surfacewindspeed,
    surfacepressure,
    surfaceuwind,
    surfacevwind,
    transportuwind,
    transportvwind,
    columnwater,
    verticalascent
};

enum class referencescope
{
    global,
    land,
    ocean
};

bool isvalidationland(planet& world, int x, int y)
{
    return world.sea(x, y) == 0 &&
        world.truelake(x, y) == 0 &&
        world.riftlakesurface(x, y) == 0;
}

struct weightedcorrelationmoments
{
    double weight = 0.0;
    double simulated = 0.0;
    double reference = 0.0;
    double simulatedsquared = 0.0;
    double referencesquared = 0.0;
    double cross = 0.0;

    void add(double simulatedvalue, double referencevalue, double sampleweight)
    {
        weight += sampleweight;
        simulated += sampleweight * simulatedvalue;
        reference += sampleweight * referencevalue;
        simulatedsquared += sampleweight * simulatedvalue * simulatedvalue;
        referencesquared += sampleweight * referencevalue * referencevalue;
        cross += sampleweight * simulatedvalue * referencevalue;
    }

    double correlation() const
    {
        if (weight <= 0.0)
            return 0.0;

        const double covariance = cross - simulated * reference / weight;
        const double simulatedvariance = simulatedsquared - simulated * simulated / weight;
        const double referencevariance = referencesquared - reference * reference / weight;
        const double denominator = std::sqrt(
            max(0.0, simulatedvariance) * max(0.0, referencevariance));
        return denominator > 0.0 ? covariance / denominator : 0.0;
    }
};

struct coarseprecipitationcell
{
    double areaweight = 0.0;
    double simulatedsum = 0.0;
    double referencesum = 0.0;
};

double fivedegreelandprecipitationcorrelation(
    planet& world,
    const climatereference::MonthlyGrid& reference)
{
    constexpr int longitudebins = 72;
    constexpr int latitudebins = 36;
    array<coarseprecipitationcell, longitudebins * latitudebins> cells{};

    for (int y = 0; y <= world.height(); y++)
    {
        const double areaweight = gridcellareaweight(y, world.height());
        const int biny = min(latitudebins - 1, y * latitudebins / (world.height() + 1));

        for (int x = 0; x <= world.width(); x++)
        {
            if (!isvalidationland(world, x, y))
                continue;

            const double observed = reference.value(0, x, y);

            if (!std::isfinite(observed))
                continue;

            double simulated = 0.0;

            for (int season = 0; season < CLIMATESEASONCOUNT; season++)
                simulated += world.seasonalrainfloat(season, x, y);

            simulated = simulated * 12.0 / static_cast<double>(CLIMATESEASONCOUNT);
            const int binx = min(longitudebins - 1, x * longitudebins / (world.width() + 1));
            coarseprecipitationcell& cell = cells[biny * longitudebins + binx];
            cell.areaweight += areaweight;
            cell.simulatedsum += areaweight * simulated;
            cell.referencesum += areaweight * observed;
        }
    }

    weightedcorrelationmoments moments;

    for (const coarseprecipitationcell& cell : cells)
    {
        if (cell.areaweight <= 0.0)
            continue;

        moments.add(
            cell.simulatedsum / cell.areaweight,
            cell.referencesum / cell.areaweight,
            cell.areaweight);
    }

    return moments.correlation();
}

bool temperatureallowsreferenceclimate(
    short expected,
    const array<float, CLIMATESEASONCOUNT>& temperatures)
{
    const string expectedcode = getclimatecode(expected);

    if (expectedcode.empty())
        return false;

    array<float, CLIMATESEASONCOUNT> ordered = temperatures;
    sort(ordered.begin(), ordered.end());
    const float mintemp = ordered.front();
    const float maxtemp = ordered.back();
    const float secondwarmesttemp = ordered[CLIMATESEASONCOUNT - 2];
    double meanannualtemp = 0.0;

    for (float temperature : temperatures)
        meanannualtemp += temperature;

    meanannualtemp /= static_cast<double>(CLIMATESEASONCOUNT);

    if (maxtemp <= 10.0f)
        return expectedcode == (maxtemp < 0.0f ? "EF" : "ET");

    if (expectedcode[0] == 'E')
        return false;

    if (expectedcode[0] == 'B')
        return expectedcode.back() == (meanannualtemp >= 18.0 ? 'h' : 'k');

    char simulatedgroup = 'D';

    if (mintemp >= 18.0f)
        simulatedgroup = 'A';
    else if (mintemp > -3.0f)
        simulatedgroup = 'C';

    if (expectedcode[0] != simulatedgroup)
        return false;

    if (simulatedgroup == 'A')
        return true;

    char simulatedheat = 'c';

    if (maxtemp >= 22.0f)
        simulatedheat = 'a';
    else if (secondwarmesttemp >= 8.0f)
        simulatedheat = 'b';
    else if (mintemp <= -38.0f)
        simulatedheat = 'd';

    return expectedcode.back() == simulatedheat;
}

bool temperatureallowsreferencemajorgroup(
    short expected,
    const array<float, CLIMATESEASONCOUNT>& temperatures)
{
    const string expectedcode = getclimatecode(expected);

    if (expectedcode.empty())
        return false;

    const float mintemp = *min_element(temperatures.begin(), temperatures.end());
    const float maxtemp = *max_element(temperatures.begin(), temperatures.end());

    if (maxtemp <= 10.0f)
        return expectedcode[0] == 'E';

    if (expectedcode[0] == 'E')
        return false;

    if (expectedcode[0] == 'B')
    {
        double meanannualtemp = 0.0;

        for (float temperature : temperatures)
            meanannualtemp += temperature;

        meanannualtemp /= static_cast<double>(CLIMATESEASONCOUNT);
        return meanannualtemp * 20.0 + 280.0 >= 0.0;
    }

    if (mintemp >= 18.0f)
        return expectedcode[0] == 'A';

    if (mintemp > -3.0f)
        return expectedcode[0] == 'C';

    return expectedcode[0] == 'D';
}

double simulatedmonthlyvalue(planet& world, monthlyreferencefield field, int season, int x, int y)
{
    if (field == monthlyreferencefield::temperature)
        return world.seasonaltemp(season, x, y);

    if (field == monthlyreferencefield::precipitation)
        return world.seasonalrainfloat(season, x, y);

    if (field == monthlyreferencefield::surfacewindspeed)
    {
        const double u = world.seasonaluwind(season, x, y);
        const double v = world.seasonalvwind(season, x, y);
        return std::sqrt(u * u + v * v);
    }

    if (field == monthlyreferencefield::surfacepressure)
        return world.seasonalpressure(season, x, y);

    if (field == monthlyreferencefield::surfaceuwind)
        return world.seasonaluwind(season, x, y);

    if (field == monthlyreferencefield::surfacevwind)
        return -world.seasonalvwind(season, x, y);

    if (field == monthlyreferencefield::transportuwind)
    {
        const double surface = world.seasonaluwind(season, x, y);
        const double upper = world.seasonalupperuwind(season, x, y);
        return surface + (upper - surface) * tuning::climate::moistureadvection::upperWindTransportFraction;
    }

    if (field == monthlyreferencefield::transportvwind)
    {
        const double surface = -world.seasonalvwind(season, x, y);
        const double upper = -world.seasonaluppervwind(season, x, y);
        return surface + (upper - surface) * tuning::climate::moistureadvection::upperWindTransportFraction;
    }

    if (field == monthlyreferencefield::columnwater)
        return world.seasonalmoisture(season, x, y);

    return static_cast<double>(world.seasonalverticalvelocity(season, x, y)) /
        tuning::climate::circulation::verticalVelocityStorageScale;
}

comparisonmetrics comparemonthlyfield(
    planet& world,
    const climatereference::MonthlyGrid& reference,
    int season,
    monthlyreferencefield field,
    referencescope scope = referencescope::land,
    double simulatedscale = 1.0)
{
    comparisonmetrics metrics;
    metrics.referencefound = true;
    metrics.dimensionsmatch = reference.width == world.width() + 1 && reference.height == world.height() + 1;

    if (!metrics.dimensionsmatch)
        return metrics;

    constexpr array<int, CLIMATESEASONCOUNT> referenceMonths = { 0, 3, 6, 9 };
    double simulatedsum = 0.0;
    double referencesum = 0.0;
    double biassum = 0.0;
    double absoluteerrorsum = 0.0;
    double squarederrorsum = 0.0;
    double sumsimulatedsquared = 0.0;
    double sumreferencesquared = 0.0;
    double sumcross = 0.0;
    double areaweightedsimulatedsum = 0.0;
    double areaweightedreferencesum = 0.0;
    double areaweightedbiassum = 0.0;
    double areaweightedabsoluteerrorsum = 0.0;
    double areaweightedsquarederrorsum = 0.0;
    double areaweightedsumsimulatedsquared = 0.0;
    double areaweightedsumreferencesquared = 0.0;
    double areaweightedsumcross = 0.0;

    for (int y = 0; y <= world.height(); y++)
    {
        const double areaweight = gridcellareaweight(y, world.height());

        for (int x = 0; x <= world.width(); x++)
        {
            if (scope == referencescope::land && !isvalidationland(world, x, y))
                continue;

            if (scope == referencescope::ocean && world.sea(x, y) == 0)
                continue;

            double simulated = 0.0;
            double observed = 0.0;

            if (season >= 0)
            {
                simulated = simulatedmonthlyvalue(world, field, season, x, y);

                if (field == monthlyreferencefield::precipitation &&
                    reference.monthCount >= 12)
                {
                    observed = 0.0;

                    for (int month = season * 3; month < season * 3 + 3; month++)
                        observed += reference.value(month, x, y);

                    observed /= 3.0;
                }
                else
                {
                    observed = reference.value(referenceMonths[season], x, y);
                }
            }
            else
            {
                for (int modelseason = 0; modelseason < CLIMATESEASONCOUNT; modelseason++)
                    simulated += simulatedmonthlyvalue(world, field, modelseason, x, y);

                simulated /= static_cast<double>(CLIMATESEASONCOUNT);

                int validmonths = 0;

                for (int month = 0; month < reference.monthCount; month++)
                {
                    const float value = reference.value(month, x, y);

                    if (std::isfinite(value))
                    {
                        observed += value;
                        validmonths++;
                    }
                }

                if (validmonths != reference.monthCount)
                    continue;

                observed /= static_cast<double>(validmonths);
            }

            simulated *= simulatedscale;

            if (!std::isfinite(observed))
                continue;

            const double difference = simulated - observed;
            metrics.comparedcells++;
            simulatedsum += simulated;
            referencesum += observed;
            biassum += difference;
            absoluteerrorsum += std::abs(difference);
            squarederrorsum += difference * difference;
            sumsimulatedsquared += simulated * simulated;
            sumreferencesquared += observed * observed;
            sumcross += simulated * observed;
            metrics.areaweight += areaweight;
            areaweightedsimulatedsum += areaweight * simulated;
            areaweightedreferencesum += areaweight * observed;
            areaweightedbiassum += areaweight * difference;
            areaweightedabsoluteerrorsum += areaweight * std::abs(difference);
            areaweightedsquarederrorsum += areaweight * difference * difference;
            areaweightedsumsimulatedsquared += areaweight * simulated * simulated;
            areaweightedsumreferencesquared += areaweight * observed * observed;
            areaweightedsumcross += areaweight * simulated * observed;
        }
    }

    metrics.simulatedmean = safeaverage(simulatedsum, metrics.comparedcells);
    metrics.referencemean = safeaverage(referencesum, metrics.comparedcells);
    metrics.meanbias = safeaverage(biassum, metrics.comparedcells);
    metrics.meanabsoluteerror = safeaverage(absoluteerrorsum, metrics.comparedcells);
    metrics.rmse = std::sqrt(safeaverage(squarederrorsum, metrics.comparedcells));

    if (metrics.comparedcells > 0)
    {
        const double count = static_cast<double>(metrics.comparedcells);
        const double numerator = sumcross - simulatedsum * referencesum / count;
        const double simulatedvariance = sumsimulatedsquared - simulatedsum * simulatedsum / count;
        const double referencevariance = sumreferencesquared - referencesum * referencesum / count;
        const double denominator = std::sqrt(max(0.0, simulatedvariance) * max(0.0, referencevariance));

        if (denominator > 0.0)
            metrics.correlation = numerator / denominator;
    }

    if (metrics.areaweight > 0.0)
    {
        metrics.areaweightedsimulatedmean = areaweightedsimulatedsum / metrics.areaweight;
        metrics.areaweightedreferencemean = areaweightedreferencesum / metrics.areaweight;
        metrics.areaweightedmeanbias = areaweightedbiassum / metrics.areaweight;
        metrics.areaweightedmeanabsoluteerror = areaweightedabsoluteerrorsum / metrics.areaweight;
        metrics.areaweightedrmse = std::sqrt(areaweightedsquarederrorsum / metrics.areaweight);
        const double numerator = areaweightedsumcross -
            areaweightedsimulatedsum * areaweightedreferencesum / metrics.areaweight;
        const double simulatedvariance = areaweightedsumsimulatedsquared -
            areaweightedsimulatedsum * areaweightedsimulatedsum / metrics.areaweight;
        const double referencevariance = areaweightedsumreferencesquared -
            areaweightedreferencesum * areaweightedreferencesum / metrics.areaweight;
        const double denominator = std::sqrt(max(0.0, simulatedvariance) * max(0.0, referencevariance));

        if (denominator > 0.0)
            metrics.areaweightedcorrelation = numerator / denominator;
    }

    return metrics;
}

void writemonthlyreferencecomparison(const filesystem::path& outputdir, planet& world)
{
    const filesystem::path referencedirectory = getappenvironment().referenceClimateDirectory;
    climatereference::MonthlyGrid temperature;
    climatereference::MonthlyGrid precipitation;
    climatereference::MonthlyGrid surfacewind;
    string temperatureerror;
    string precipitationerror;
    string surfacewinderror;
    const bool hastemperature = climatereference::loadMonthlyGrid(
        referencedirectory / "worldclim_tavg_monthly.uwclim", "tavg", temperature, &temperatureerror);
    const bool hasprecipitation = climatereference::loadMonthlyGrid(
        referencedirectory / "worldclim_prec_monthly.uwclim", "prec", precipitation, &precipitationerror);
    const bool hassurfacewind = climatereference::loadMonthlyGrid(
        referencedirectory / "worldclim_wind_monthly.uwclim", "wind", surfacewind, &surfacewinderror);
    ofstream output(outputdir / "monthly_climate_reference_comparison.csv");

    if (!output.is_open())
        return;

    output << "variable,period,reference_month,compared_cells,simulated_mean,reference_mean,mean_bias,mae,rmse,correlation,area_weight,area_weighted_simulated_mean,area_weighted_reference_mean,area_weighted_mean_bias,area_weighted_mae,area_weighted_rmse,area_weighted_correlation\n";
    output << fixed << setprecision(6);

    if (!hastemperature)
        output << "temperature,error," << csvescape(temperatureerror) << "\n";

    if (!hasprecipitation)
        output << "precipitation,error," << csvescape(precipitationerror) << "\n";

    if (!hassurfacewind)
        output << "surface_wind_speed,error," << csvescape(surfacewinderror) << "\n";

    constexpr array<const char*, CLIMATESEASONCOUNT> periodNames = { "january", "april", "july", "october" };
    constexpr array<const char*, CLIMATESEASONCOUNT> precipitationPeriodNames = {
        "january-march", "april-june", "july-september", "october-december"
    };
    constexpr array<int, CLIMATESEASONCOUNT> referenceMonths = { 1, 4, 7, 10 };

    auto writecomparisons = [&](const char* variable, const climatereference::MonthlyGrid& grid, monthlyreferencefield field)
    {
        for (int season = -1; season < CLIMATESEASONCOUNT; season++)
        {
            const comparisonmetrics metrics = comparemonthlyfield(world, grid, season, field);
            output
                << variable << ','
                << (season < 0
                    ? "annual_mean"
                    : (field == monthlyreferencefield::precipitation
                        ? precipitationPeriodNames[season]
                        : periodNames[season])) << ','
                << (season < 0 || field == monthlyreferencefield::precipitation
                    ? 0
                    : referenceMonths[season]) << ','
                << metrics.comparedcells << ','
                << metrics.simulatedmean << ','
                << metrics.referencemean << ','
                << metrics.meanbias << ','
                << metrics.meanabsoluteerror << ','
                << metrics.rmse << ','
                << metrics.correlation << ','
                << metrics.areaweight << ','
                << metrics.areaweightedsimulatedmean << ','
                << metrics.areaweightedreferencemean << ','
                << metrics.areaweightedmeanbias << ','
                << metrics.areaweightedmeanabsoluteerror << ','
                << metrics.areaweightedrmse << ','
                << metrics.areaweightedcorrelation << '\n';
        }
    };

    if (hastemperature)
        writecomparisons("temperature_c", temperature, monthlyreferencefield::temperature);

    if (hasprecipitation)
        writecomparisons("precipitation_mm_month", precipitation, monthlyreferencefield::precipitation);

    if (hassurfacewind)
        writecomparisons("surface_wind_speed_m_s", surfacewind, monthlyreferencefield::surfacewindspeed);
}

void writeimergreferencecomparison(const filesystem::path& outputdir, planet& world)
{
    climatereference::MonthlyGrid precipitation;
    string error;
    const filesystem::path referencepath =
        getappenvironment().referenceClimateDirectory / "imerg_prec_annual.uwclim";
    const bool hasprecipitation = climatereference::loadMonthlyGrid(
        referencepath, "prec_annual", precipitation, &error);
    ofstream output(outputdir / "annual_imerg_precipitation_comparison.csv");

    if (!output.is_open())
        return;

    output << "scope,compared_cells,simulated_mean_mm_year,reference_mean_mm_year,mean_bias_mm_year,mae_mm_year,rmse_mm_year,correlation,area_weight,area_weighted_simulated_mean_mm_year,area_weighted_reference_mean_mm_year,area_weighted_mean_bias_mm_year,area_weighted_mae_mm_year,area_weighted_rmse_mm_year,area_weighted_correlation\n";

    if (!hasprecipitation)
    {
        output << "error," << csvescape(error) << '\n';
        return;
    }

    constexpr array<const char*, 3> names = { "global", "land", "ocean" };
    constexpr array<referencescope, 3> scopes = {
        referencescope::global,
        referencescope::land,
        referencescope::ocean
    };
    output << fixed << setprecision(6);

    for (size_t index = 0; index < scopes.size(); index++)
    {
        const comparisonmetrics metrics = comparemonthlyfield(
            world,
            precipitation,
            -1,
            monthlyreferencefield::precipitation,
            scopes[index],
            12.0);
        output
            << names[index] << ','
            << metrics.comparedcells << ','
            << metrics.simulatedmean << ','
            << metrics.referencemean << ','
            << metrics.meanbias << ','
            << metrics.meanabsoluteerror << ','
            << metrics.rmse << ','
            << metrics.correlation << ','
            << metrics.areaweight << ','
            << metrics.areaweightedsimulatedmean << ','
            << metrics.areaweightedreferencemean << ','
            << metrics.areaweightedmeanbias << ','
            << metrics.areaweightedmeanabsoluteerror << ','
            << metrics.areaweightedrmse << ','
            << metrics.areaweightedcorrelation << '\n';
    }
}

struct physicalreferencespecification
{
    const char* filename;
    const char* referencevariable;
    const char* outputvariable;
    monthlyreferencefield simulatedfield;
};

void writephysicalreferencecomparison(const filesystem::path& outputdir, planet& world)
{
    constexpr array<physicalreferencespecification, 8> specifications = {
        physicalreferencespecification{ "era5_slp_anom_monthly.uwclim", "slp_anom", "sea_level_pressure_anomaly_hpa", monthlyreferencefield::surfacepressure },
        { "era5_u10m_monthly.uwclim", "u10m", "surface_eastward_wind_m_s", monthlyreferencefield::surfaceuwind },
        { "era5_v10m_monthly.uwclim", "v10m", "surface_northward_wind_m_s", monthlyreferencefield::surfacevwind },
        { "era5_u850_monthly.uwclim", "u850", "transport_eastward_wind_m_s", monthlyreferencefield::transportuwind },
        { "era5_v850_monthly.uwclim", "v850", "transport_northward_wind_m_s", monthlyreferencefield::transportvwind },
        { "era5_tcwv_monthly.uwclim", "tcwv", "column_water_kg_m2", monthlyreferencefield::columnwater },
        { "era5_w500_ascent_monthly.uwclim", "w500_ascent", "midlevel_ascent_hpa_day", monthlyreferencefield::verticalascent },
        { "era5_pr_monthly.uwclim", "pr", "precipitation_mm_month", monthlyreferencefield::precipitation }
    };
    constexpr array<const char*, 3> scopeNames = { "global", "land", "ocean" };
    constexpr array<referencescope, 3> scopes = {
        referencescope::global,
        referencescope::land,
        referencescope::ocean
    };
    constexpr array<const char*, CLIMATESEASONCOUNT> periodNames = {
        "january", "april", "july", "october"
    };
    constexpr array<const char*, CLIMATESEASONCOUNT> precipitationPeriodNames = {
        "january-march", "april-june", "july-september", "october-december"
    };
    constexpr array<int, CLIMATESEASONCOUNT> referenceMonths = { 1, 4, 7, 10 };
    const filesystem::path referencedirectory = getappenvironment().referenceClimateDirectory;
    ofstream output(outputdir / "monthly_physical_reference_comparison.csv");

    if (!output.is_open())
        return;

    output << "variable,scope,period,reference_month,compared_cells,simulated_mean,reference_mean,mean_bias,mae,rmse,correlation,area_weight,area_weighted_simulated_mean,area_weighted_reference_mean,area_weighted_mean_bias,area_weighted_mae,area_weighted_rmse,area_weighted_correlation\n";
    output << fixed << setprecision(6);

    for (const physicalreferencespecification& specification : specifications)
    {
        climatereference::MonthlyGrid reference;
        string error;

        if (!climatereference::loadMonthlyGrid(
            referencedirectory / specification.filename,
            specification.referencevariable,
            reference,
            &error))
        {
            output << specification.outputvariable << ",error," << csvescape(error) << '\n';
            continue;
        }

        for (size_t scopeindex = 0; scopeindex < scopes.size(); scopeindex++)
        {
            for (int season = -1; season < CLIMATESEASONCOUNT; season++)
            {
                const comparisonmetrics metrics = comparemonthlyfield(
                    world,
                    reference,
                    season,
                    specification.simulatedfield,
                    scopes[scopeindex]);
                output
                    << specification.outputvariable << ','
                    << scopeNames[scopeindex] << ','
                    << (season < 0
                        ? "annual_mean"
                        : (specification.simulatedfield == monthlyreferencefield::precipitation
                            ? precipitationPeriodNames[season]
                            : periodNames[season])) << ','
                    << (season < 0 ||
                        specification.simulatedfield == monthlyreferencefield::precipitation
                            ? 0
                            : referenceMonths[season]) << ','
                    << metrics.comparedcells << ','
                    << metrics.simulatedmean << ','
                    << metrics.referencemean << ','
                    << metrics.meanbias << ','
                    << metrics.meanabsoluteerror << ','
                    << metrics.rmse << ','
                    << metrics.correlation << ','
                    << metrics.areaweight << ','
                    << metrics.areaweightedsimulatedmean << ','
                    << metrics.areaweightedreferencemean << ','
                    << metrics.areaweightedmeanbias << ','
                    << metrics.areaweightedmeanabsoluteerror << ','
                    << metrics.areaweightedrmse << ','
                    << metrics.areaweightedcorrelation << '\n';
            }
        }
    }
}

struct ipccregionlabel
{
    int id = -1;
    string acronym;
    string name;
    string continent;
    string type;
};

struct ipccregiondiagnostics
{
    long long landcells = 0;
    double landweight = 0.0;
    double climateweight = 0.0;
    double exactweight = 0.0;
    double groupweight = 0.0;
    double temperatureweight = 0.0;
    double simulatedtemperature = 0.0;
    double referencetemperature = 0.0;
    double simulatedwarmesttemperature = 0.0;
    double referencewarmesttemperature = 0.0;
    double precipitationweight = 0.0;
    double simulatedprecipitation = 0.0;
    double referenceprecipitation = 0.0;
    double storedzeroprecipitationweight = 0.0;
    weightedcorrelationmoments temperaturecorrelation;
    weightedcorrelationmoments precipitationcorrelation;
    array<double, 5> simulatedgroupweights{};
    array<double, 5> referencegroupweights{};
    double simulatedetweight = 0.0;
    double referenceetweight = 0.0;
    double simulatedefweight = 0.0;
    double referenceefweight = 0.0;
};

vector<ipccregionlabel> loadipccregionlabels(const filesystem::path& path)
{
    ifstream input(path);
    vector<ipccregionlabel> labels;
    string line;
    getline(input, line);

    while (getline(input, line))
    {
        array<string, 5> fields{};
        size_t start = 0;

        for (size_t field = 0; field < fields.size(); field++)
        {
            const size_t separator = line.find('\t', start);
            fields[field] = line.substr(start, separator == string::npos ? string::npos : separator - start);
            start = separator == string::npos ? line.size() : separator + 1;
        }

        try
        {
            ipccregionlabel label;
            label.id = stoi(fields[0]);
            label.acronym = fields[1];
            label.name = fields[2];
            label.continent = fields[3];
            label.type = fields[4];

            if (label.id >= static_cast<int>(labels.size()))
                labels.resize(label.id + 1);

            labels[label.id] = label;
        }
        catch (const exception&)
        {
            return {};
        }
    }

    return labels;
}

void writeipccregioncomparison(const filesystem::path& outputdir, planet& world)
{
    const filesystem::path referencedirectory = getappenvironment().referenceClimateDirectory;
    const vector<ipccregionlabel> labels = loadipccregionlabels(
        referencedirectory / "ipcc_ar6_regions.tsv");
    sf::Image regionmask;
    sf::Image koppenreference;
    climatereference::MonthlyGrid temperaturereference;
    climatereference::MonthlyGrid precipitationreference;
    string temperatureerror;
    string precipitationerror;
    ofstream output(outputdir / "ipcc_region_climate_comparison.csv");

    if (!output.is_open())
        return;

    if (labels.empty() ||
        !regionmask.loadFromFile((referencedirectory / "ipcc_ar6_regions.png").string()) ||
        !koppenreference.loadFromFile(getappenvironment().earthKoppenImagePath.string()) ||
        !climatereference::loadMonthlyGrid(
            referencedirectory / "worldclim_tavg_monthly.uwclim",
            "tavg",
            temperaturereference,
            &temperatureerror) ||
        !climatereference::loadMonthlyGrid(
            referencedirectory / "worldclim_prec_monthly.uwclim",
            "prec",
            precipitationreference,
            &precipitationerror))
    {
        output << "status,error\nreference_unavailable,"
            << csvescape(temperatureerror.empty() ? precipitationerror : temperatureerror) << '\n';
        return;
    }

    const sf::Vector2u regionsize = regionmask.getSize();
    const sf::Vector2u koppensize = koppenreference.getSize();
    const bool dimensionsmatch =
        regionsize.x == static_cast<unsigned int>(world.width() + 1) &&
        regionsize.y == static_cast<unsigned int>(world.height() + 1) &&
        koppensize == regionsize &&
        temperaturereference.width == world.width() + 1 &&
        temperaturereference.height == world.height() + 1 &&
        precipitationreference.width == world.width() + 1 &&
        precipitationreference.height == world.height() + 1;

    if (!dimensionsmatch)
    {
        output << "status,error\ndimension_mismatch,reference grids do not match benchmark world\n";
        return;
    }

    vector<ipccregiondiagnostics> regions(labels.size());

    for (int y = 0; y <= world.height(); y++)
    {
        const double areaweight = gridcellareaweight(y, world.height());

        for (int x = 0; x <= world.width(); x++)
        {
            if (!isvalidationland(world, x, y))
                continue;

            const int regionid = static_cast<int>(regionmask.getPixel(x, y).r) - 1;

            if (regionid < 0 || regionid >= static_cast<int>(regions.size()) || labels[regionid].id < 0)
                continue;

            ipccregiondiagnostics& region = regions[regionid];
            region.landcells++;
            region.landweight += areaweight;
            bool storedzero = true;
            double simulatedtemperature = 0.0;
            double simulatedwarmest = -numeric_limits<double>::infinity();
            double simulatedprecipitation = 0.0;

            for (int season = 0; season < CLIMATESEASONCOUNT; season++)
            {
                const double temperature = world.seasonaltemp(season, x, y);
                const double precipitation = world.seasonalrainfloat(season, x, y);
                simulatedtemperature += temperature;
                simulatedwarmest = max(simulatedwarmest, temperature);
                simulatedprecipitation += precipitation;
                storedzero = storedzero && precipitation <= 0.0;
            }

            simulatedtemperature /= static_cast<double>(CLIMATESEASONCOUNT);
            simulatedprecipitation *= 12.0 / static_cast<double>(CLIMATESEASONCOUNT);

            if (storedzero)
                region.storedzeroprecipitationweight += areaweight;

            double referencetemperature = 0.0;
            double referencewarmest = -numeric_limits<double>::infinity();
            double referenceprecipitation = 0.0;
            bool validtemperature = true;
            bool validprecipitation = true;

            for (int month = 0; month < 12; month++)
            {
                const double temperature = temperaturereference.value(month, x, y);
                const double precipitation = precipitationreference.value(month, x, y);
                validtemperature = validtemperature && std::isfinite(temperature);
                validprecipitation = validprecipitation && std::isfinite(precipitation);

                if (std::isfinite(temperature))
                {
                    referencetemperature += temperature;
                    referencewarmest = max(referencewarmest, temperature);
                }

                if (std::isfinite(precipitation))
                    referenceprecipitation += precipitation;
            }

            if (validtemperature)
            {
                referencetemperature /= 12.0;
                region.temperatureweight += areaweight;
                region.simulatedtemperature += areaweight * simulatedtemperature;
                region.referencetemperature += areaweight * referencetemperature;
                region.simulatedwarmesttemperature += areaweight * simulatedwarmest;
                region.referencewarmesttemperature += areaweight * referencewarmest;
                region.temperaturecorrelation.add(
                    simulatedtemperature, referencetemperature, areaweight);
            }

            if (validprecipitation)
            {
                region.precipitationweight += areaweight;
                region.simulatedprecipitation += areaweight * simulatedprecipitation;
                region.referenceprecipitation += areaweight * referenceprecipitation;
                region.precipitationcorrelation.add(
                    simulatedprecipitation, referenceprecipitation, areaweight);
            }

            const short simulated = comparableclimate(static_cast<short>(world.climate(x, y)));
            int colourdistance = 0;
            const short expected = comparableclimate(
                nearestbenchmarkclimate(koppenreference.getPixel(x, y), colourdistance));

            if (simulated < 1 || simulated > 31 || expected < 1 || expected > 31 ||
                colourdistance > benchmarkcolourdistancelimitsquared)
            {
                continue;
            }

            region.climateweight += areaweight;
            const int simulatedgroup = climatemajorgroup(simulated);
            const int expectedgroup = climatemajorgroup(expected);
            region.simulatedgroupweights[simulatedgroup] += areaweight;
            region.referencegroupweights[expectedgroup] += areaweight;
            region.simulatedetweight += simulated == 30 ? areaweight : 0.0;
            region.referenceetweight += expected == 30 ? areaweight : 0.0;
            region.simulatedefweight += simulated == 31 ? areaweight : 0.0;
            region.referenceefweight += expected == 31 ? areaweight : 0.0;

            if (simulated == expected)
                region.exactweight += areaweight;
            if (simulatedgroup == expectedgroup)
                region.groupweight += areaweight;
        }
    }

    output << "region_id,acronym,name,continent,type,land_cells,area_weight,climate_compared_area,exact_accuracy,group_accuracy,simulated_annual_temperature_c,reference_annual_temperature_c,temperature_bias_c,temperature_correlation,simulated_warmest_temperature_c,reference_warmest_temperature_c,simulated_annual_precipitation_mm,reference_annual_precipitation_mm,precipitation_ratio,precipitation_correlation,stored_zero_precipitation_fraction,simulated_a_fraction,reference_a_fraction,simulated_b_fraction,reference_b_fraction,simulated_c_fraction,reference_c_fraction,simulated_d_fraction,reference_d_fraction,simulated_e_fraction,reference_e_fraction,simulated_et_fraction,reference_et_fraction,simulated_ef_fraction,reference_ef_fraction\n";
    output << fixed << setprecision(6);

    for (size_t regionid = 0; regionid < labels.size(); regionid++)
    {
        const ipccregionlabel& label = labels[regionid];
        const ipccregiondiagnostics& region = regions[regionid];

        if (label.id < 0 || region.landcells == 0)
            continue;

        const double simulatedprecipitation = region.precipitationweight > 0.0 ?
            region.simulatedprecipitation / region.precipitationweight : 0.0;
        const double referenceprecipitation = region.precipitationweight > 0.0 ?
            region.referenceprecipitation / region.precipitationweight : 0.0;
        const double simulatedtemperature = region.temperatureweight > 0.0 ?
            region.simulatedtemperature / region.temperatureweight : 0.0;
        const double referencetemperature = region.temperatureweight > 0.0 ?
            region.referencetemperature / region.temperatureweight : 0.0;
        output
            << label.id << ','
            << label.acronym << ','
            << csvescape(label.name) << ','
            << label.continent << ','
            << label.type << ','
            << region.landcells << ','
            << region.landweight << ','
            << region.climateweight << ','
            << (region.climateweight > 0.0 ? region.exactweight / region.climateweight : 0.0) << ','
            << (region.climateweight > 0.0 ? region.groupweight / region.climateweight : 0.0) << ','
            << simulatedtemperature << ','
            << referencetemperature << ','
            << simulatedtemperature - referencetemperature << ','
            << region.temperaturecorrelation.correlation() << ','
            << (region.temperatureweight > 0.0 ? region.simulatedwarmesttemperature / region.temperatureweight : 0.0) << ','
            << (region.temperatureweight > 0.0 ? region.referencewarmesttemperature / region.temperatureweight : 0.0) << ','
            << simulatedprecipitation << ','
            << referenceprecipitation << ','
            << (referenceprecipitation > 0.0 ? simulatedprecipitation / referenceprecipitation : 0.0) << ','
            << region.precipitationcorrelation.correlation() << ','
            << (region.landweight > 0.0 ? region.storedzeroprecipitationweight / region.landweight : 0.0);

        for (int group = 0; group < 5; group++)
        {
            output
                << ',' << (region.climateweight > 0.0 ? region.simulatedgroupweights[group] / region.climateweight : 0.0)
                << ',' << (region.climateweight > 0.0 ? region.referencegroupweights[group] / region.climateweight : 0.0);
        }

        output
            << ',' << (region.climateweight > 0.0 ? region.simulatedetweight / region.climateweight : 0.0)
            << ',' << (region.climateweight > 0.0 ? region.referenceetweight / region.climateweight : 0.0)
            << ',' << (region.climateweight > 0.0 ? region.simulatedefweight / region.climateweight : 0.0)
            << ',' << (region.climateweight > 0.0 ? region.referenceefweight / region.climateweight : 0.0);

        output << '\n';
    }
}

short climatefromcode(const string& code)
{
    for (short candidate = 1; candidate <= 31; candidate++)
    {
        if (getclimatecode(candidate) == code)
            return candidate;
    }

    return 0;
}

string csvescape(const string& value)
{
    if (value.find(',') == string::npos && value.find('"') == string::npos)
        return value;

    string escaped = "\"";

    for (char ch : value)
    {
        if (ch == '"')
            escaped += "\"\"";
        else
            escaped += ch;
    }

    escaped += '"';
    return escaped;
}

float latitudeforrow(int row, int height)
{
    if (height <= 0)
        return 0.0f;

    return 90.0f - (180.0f * static_cast<float>(row) / static_cast<float>(height));
}

filesystem::path climatevalidationoutputdirectory(long seed)
{
    filesystem::path outputroot = getappenvironment().profilingWorkbookPath.parent_path();

    if (outputroot.empty())
        outputroot = filesystem::current_path();

    outputroot /= "validation";
    outputroot /= "seed_" + to_string(seed);
    filesystem::create_directories(outputroot);

    return outputroot;
}

filesystem::path climatevalidationoutputdirectory()
{
    return climatevalidationoutputdirectory(worldgenerationdebugseed());
}

vector<string> orderedclimatecodes()
{
    vector<string> codes;
    codes.reserve(31);

    for (short climate = 1; climate <= 31; climate++)
        codes.push_back(getclimatecode(climate));

    return codes;
}

vector<long long> collectsimulatedclimatecounts(planet& world)
{
    vector<long long> counts(31, 0);
    const int width = world.width();
    const int height = world.height();

    for (int y = 0; y <= height; y++)
    {
        for (int x = 0; x <= width; x++)
        {
            if (!isvalidationland(world, x, y))
                continue;

            const short climate = static_cast<short>(world.climate(x, y));

            if (climate >= 1 && climate <= 31)
                counts[climate - 1]++;
        }
    }

    return counts;
}

void printreferencedfcdriverreport(planet& world)
{
    sf::Image reference;
    const filesystem::path referencepath = getappenvironment().earthKoppenImagePath;

    if (reference.loadFromFile(referencepath.string()) == false)
        return;

    const int width = world.width();
    const int height = world.height();
    const sf::Vector2u referencesize = reference.getSize();

    if (referencesize.x != static_cast<unsigned int>(width + 1) || referencesize.y != static_cast<unsigned int>(height + 1))
        return;

    array<climatedriverstats, 31> bysimulated{};

    for (int y = 0; y <= height; y++)
    {
        for (int x = 0; x <= width; x++)
        {
            if (!isvalidationland(world, x, y))
                continue;

            int distancesquared = 0;
            const short expected = nearestbenchmarkclimate(reference.getPixel(x, y), distancesquared);

            if (distancesquared > benchmarkcolourdistancelimitsquared || expected != 28)
                continue;

            const short simulated = comparableclimate(static_cast<short>(world.climate(x, y)));

            if (simulated < 1 || simulated > 31)
                continue;

            array<float, CLIMATESEASONCOUNT> temps{};
            array<float, CLIMATESEASONCOUNT> rains{};
            array<int, CLIMATESEASONCOUNT> order = { 0, 1, 2, 3 };
            float meanannualtemp = 0.0f;
            float annualrain = 0.0f;

            for (int season = 0; season < CLIMATESEASONCOUNT; season++)
            {
                temps[season] = static_cast<float>(world.seasonaltemp(season, x, y));
                rains[season] = world.seasonalrainfloat(season, x, y);
                meanannualtemp += temps[season];
                annualrain += rains[season];
            }

            meanannualtemp /= static_cast<float>(CLIMATESEASONCOUNT);
            annualrain *= 3.0f;
            sort(order.begin(), order.end(), [&](int left, int right) { return temps[left] < temps[right]; });

            const float mintemp = temps[order[0]];
            const float maxtemp = temps[order[3]];
            const float driestcoldrain = min(rains[order[0]], rains[order[1]]);
            const float wettestwarmrain = max(rains[order[2]], rains[order[3]]);
            const float totalseasonalrain = annualrain / 3.0f;
            const float warmhalffraction = totalseasonalrain > 0.0f
                ? (rains[order[2]] + rains[order[3]]) / totalseasonalrain
                : 0.0f;
            float drythreshold = meanannualtemp * 20.0f;

            if (warmhalffraction >= 0.7f)
                drythreshold += 280.0f;
            else if (warmhalffraction >= 0.3f)
                drythreshold += 140.0f;

            climatedriverstats& stats = bysimulated[simulated - 1];
            stats.cells++;
            stats.meanannualtemp += meanannualtemp;
            stats.mintemp += mintemp;
            stats.maxtemp += maxtemp;
            stats.annualrain += annualrain;
            stats.drythreshold += drythreshold;
            stats.warmhalffraction += warmhalffraction;
            stats.driestcoldrain += driestcoldrain;
            stats.wettestwarmrain += wettestwarmrain;

            if (drythreshold > 0.0f)
                stats.raintothreshold += annualrain / drythreshold;

            if (annualrain < drythreshold * 0.5f)
                stats.bwcells++;

            if (annualrain <= drythreshold)
                stats.aridcells++;

            if (maxtemp > 10.0f && mintemp <= -3.0f)
                stats.dthermalcells++;

            if (maxtemp <= 10.0f)
                stats.polarcells++;

            if (climatekoppen::isWinterDry(driestcoldrain, wettestwarmrain))
                stats.winterdrycells++;
        }
    }

    vector<int> order(31);

    for (int climate = 0; climate < 31; climate++)
        order[climate] = climate;

    sort(order.begin(), order.end(), [&](int left, int right)
    {
        return bysimulated[left].cells > bysimulated[right].cells;
    });

    cout << "Reference Dfc climate-driver summary by simulated class:" << '\n';

    for (int climate : order)
    {
        const climatedriverstats& stats = bysimulated[climate];

        if (stats.cells == 0)
            continue;

        const double cells = static_cast<double>(stats.cells);
        cout
            << "reference_Dfc simulated=" << getclimatecode(static_cast<short>(climate + 1))
            << " cells=" << stats.cells
            << " mean_annual_temp=" << stats.meanannualtemp / cells
            << " mean_min_temp=" << stats.mintemp / cells
            << " mean_max_temp=" << stats.maxtemp / cells
            << " mean_annual_rain=" << stats.annualrain / cells
            << " mean_dry_threshold=" << stats.drythreshold / cells
            << " mean_rain_to_threshold=" << stats.raintothreshold / cells
            << " mean_warm_half_fraction=" << stats.warmhalffraction / cells
            << " mean_driest_cold_rain=" << stats.driestcoldrain / cells
            << " mean_wettest_warm_rain=" << stats.wettestwarmrain / cells
            << " bw_fraction=" << static_cast<double>(stats.bwcells) / cells
            << " arid_fraction=" << static_cast<double>(stats.aridcells) / cells
            << " d_thermal_fraction=" << static_cast<double>(stats.dthermalcells) / cells
            << " polar_fraction=" << static_cast<double>(stats.polarcells) / cells
            << " winter_dry_fraction=" << static_cast<double>(stats.winterdrycells) / cells
            << '\n';
    }
}

void printsimulatedclassconfusionreport(planet& world, short targetclimate)
{
    sf::Image reference;
    const filesystem::path referencepath = getappenvironment().earthKoppenImagePath;

    if (reference.loadFromFile(referencepath.string()) == false)
        return;

    const int width = world.width();
    const int height = world.height();
    const sf::Vector2u referencesize = reference.getSize();

    if (referencesize.x != static_cast<unsigned int>(width + 1) || referencesize.y != static_cast<unsigned int>(height + 1))
        return;

    array<long long, 31> expectedcounts{};
    long long total = 0;

    for (int y = 0; y <= height; y++)
    {
        for (int x = 0; x <= width; x++)
        {
            if (world.sea(x, y) == 1)
                continue;

            const short simulated = comparableclimate(static_cast<short>(world.climate(x, y)));

            if (simulated != targetclimate)
                continue;

            int distancesquared = 0;
            const short expected = comparableclimate(nearestbenchmarkclimate(reference.getPixel(x, y), distancesquared));

            if (distancesquared > benchmarkcolourdistancelimitsquared || expected < 1 || expected > 31)
                continue;

            expectedcounts[expected - 1]++;
            total++;
        }
    }

    vector<int> order(31);

    for (int climate = 0; climate < 31; climate++)
        order[climate] = climate;

    sort(order.begin(), order.end(), [&](int left, int right)
    {
        return expectedcounts[left] > expectedcounts[right];
    });

    cout << "Simulated " << getclimatecode(targetclimate) << " reference-class summary:" << '\n';

    for (int climate : order)
    {
        if (expectedcounts[climate] == 0)
            continue;

        cout
            << "simulated_" << getclimatecode(targetclimate)
            << " expected=" << getclimatecode(static_cast<short>(climate + 1))
            << " cells=" << expectedcounts[climate]
            << " fraction=" << static_cast<double>(expectedcounts[climate]) / static_cast<double>(total)
            << '\n';
    }
}

vector<long long> referenceclimatecounts()
{
    return
    {
        20086, 13831, 24368, 24368, 64970, 24332, 24110, 29277,
        5124, 3226, 6, 11798, 4443, 6, 19304, 12052, 17,
        758, 2188, 9114, 65, 3831, 6466, 12105, 1494,
        6413, 40215, 86356, 757, 48284, 215627
    };
}

string jsonescape(const string& value)
{
    ostringstream escaped;
    escaped << hex << setfill('0');

    for (const unsigned char ch : value)
    {
        switch (ch)
        {
        case '"': escaped << "\\\""; break;
        case '\\': escaped << "\\\\"; break;
        case '\b': escaped << "\\b"; break;
        case '\f': escaped << "\\f"; break;
        case '\n': escaped << "\\n"; break;
        case '\r': escaped << "\\r"; break;
        case '\t': escaped << "\\t"; break;
        default:
            if (ch < 0x20)
                escaped << "\\u" << setw(4) << static_cast<int>(ch);
            else
                escaped << static_cast<char>(ch);
            break;
        }
    }

    return escaped.str();
}

string climatebenchmarktimestamp()
{
    const time_t now = chrono::system_clock::to_time_t(chrono::system_clock::now());
    tm localtime{};
    localtime_s(&localtime, &now);

    ostringstream timestamp;
    timestamp << put_time(&localtime, "%Y%m%d%H%M%S");
    return timestamp.str();
}

int nextclimatebenchmarkrunid()
{
    const filesystem::path logpath = getappenvironment().climateBenchmarkRunLogPath;
    ifstream logfile(logpath);
    int maximumid = 1;

    if (!logfile.is_open())
        return maximumid + 1;

    const regex idpattern(R"(^\s*"id"\s*:\s*([0-9]+)\s*,?\s*$)");
    string line;

    while (getline(logfile, line))
    {
        smatch found;

        if (!regex_match(line, found, idpattern))
            continue;

        try
        {
            maximumid = max(maximumid, stoi(found[1].str()));
        }
        catch (const exception&)
        {
            return -1;
        }
    }

    return maximumid + 1;
}

void collectbenchmarktemperaturemetrics(planet& world, benchmarkphysicsmetrics& metrics)
{
    double landweight = 0.0;
    double clippedweight = 0.0;

    for (int y = 0; y <= world.height(); y++)
    {
        const double areaweight = gridcellareaweight(y, world.height());

        for (int x = 0; x <= world.width(); x++)
        {
            if (!isvalidationland(world, x, y))
                continue;

            landweight += areaweight;

            for (int season = 0; season < CLIMATESEASONCOUNT; season++)
            {
                if (world.seasonalrain(season, x, y) >= (numeric_limits<short>::max)())
                {
                    clippedweight += areaweight;
                    break;
                }
            }
        }
    }

    if (landweight > 0.0)
        metrics.areaweightedlandrainfallstorageclippingfraction = clippedweight / landweight;

    sf::Image koppenreference;

    if (!koppenreference.loadFromFile(getappenvironment().earthKoppenImagePath.string()))
        return;

    if (koppenreference.getSize().x != static_cast<unsigned int>(world.width() + 1) ||
        koppenreference.getSize().y != static_cast<unsigned int>(world.height() + 1))
    {
        return;
    }

    climatereference::MonthlyGrid temperaturereference;
    string temperatureerror;
    const bool hastemperaturereference = climatereference::loadMonthlyGrid(
        getappenvironment().referenceClimateDirectory / "worldclim_tavg_monthly.uwclim",
        "tavg",
        temperaturereference,
        &temperatureerror);
    double comparedweight = 0.0;
    double grouplockedweight = 0.0;
    double signaturelockedweight = 0.0;
    double northernweight = 0.0;
    double northerngrouplockedweight = 0.0;
    double northernsignaturelockedweight = 0.0;
    double northerntemperatureweight = 0.0;
    double northernsimulatedwarmest = 0.0;
    double northernreferencewarmest = 0.0;
    double southernweight = 0.0;
    double southernsimulatedefweight = 0.0;
    double southernreferenceefweight = 0.0;
    double southernsimulateddweight = 0.0;
    double southernreferencedweight = 0.0;
    double southerntemperatureweight = 0.0;
    double southernsimulatedwarmest = 0.0;
    double southernreferencewarmest = 0.0;

    for (int y = 0; y <= world.height(); y++)
    {
        const double areaweight = gridcellareaweight(y, world.height());
        const float latitude = latitudeforrow(y, world.height());
        const bool northern5070 = latitude >= 50.0f && latitude <= 70.0f;
        const bool southern6090 = latitude <= -60.0f;

        for (int x = 0; x <= world.width(); x++)
        {
            if (!isvalidationland(world, x, y))
                continue;

            int distancesquared = 0;
            const short expected = comparableclimate(
                nearestbenchmarkclimate(koppenreference.getPixel(x, y), distancesquared));

            if (distancesquared > benchmarkcolourdistancelimitsquared || expected < 1 || expected > 31)
                continue;

            array<float, CLIMATESEASONCOUNT> temperatures{};

            for (int season = 0; season < CLIMATESEASONCOUNT; season++)
                temperatures[season] = static_cast<float>(world.seasonaltemp(season, x, y));

            const float maxtemp = *max_element(temperatures.begin(), temperatures.end());
            const bool grouplocked = !temperatureallowsreferencemajorgroup(expected, temperatures);
            const bool signaturelocked = !temperatureallowsreferenceclimate(expected, temperatures);
            comparedweight += areaweight;

            if (grouplocked)
                grouplockedweight += areaweight;

            if (signaturelocked)
                signaturelockedweight += areaweight;

            if (northern5070)
            {
                northernweight += areaweight;

                if (grouplocked)
                    northerngrouplockedweight += areaweight;

                if (signaturelocked)
                    northernsignaturelockedweight += areaweight;
            }

            if (southern6090)
            {
                const short simulated = comparableclimate(static_cast<short>(world.climate(x, y)));
                southernweight += areaweight;

                if (simulated == 31)
                    southernsimulatedefweight += areaweight;

                if (expected == 31)
                    southernreferenceefweight += areaweight;

                if (simulated >= 18 && simulated <= 29)
                    southernsimulateddweight += areaweight;

                if (expected >= 18 && expected <= 29)
                    southernreferencedweight += areaweight;
            }

            if (!hastemperaturereference || (!northern5070 && !southern6090))
                continue;

            double observedwarmest = -numeric_limits<double>::infinity();
            bool validreference = true;

            for (int month = 0; month < temperaturereference.monthCount; month++)
            {
                const double observed = temperaturereference.value(month, x, y);

                if (!std::isfinite(observed))
                {
                    validreference = false;
                    break;
                }

                observedwarmest = max(observedwarmest, observed);
            }

            if (!validreference)
                continue;

            if (northern5070)
            {
                northerntemperatureweight += areaweight;
                northernsimulatedwarmest += areaweight * maxtemp;
                northernreferencewarmest += areaweight * observedwarmest;
            }

            if (southern6090)
            {
                southerntemperatureweight += areaweight;
                southernsimulatedwarmest += areaweight * maxtemp;
                southernreferencewarmest += areaweight * observedwarmest;
            }
        }
    }

    if (comparedweight > 0.0)
    {
        metrics.areaweightedtemperaturegrouplockedfraction = grouplockedweight / comparedweight;
        metrics.areaweightedthermalsignaturelockedfraction = signaturelockedweight / comparedweight;
    }

    if (northernweight > 0.0)
    {
        metrics.areaweightednorthern5070temperaturegrouplockedfraction =
            northerngrouplockedweight / northernweight;
        metrics.areaweightednorthern5070thermalsignaturelockedfraction =
            northernsignaturelockedweight / northernweight;
    }

    if (northerntemperatureweight > 0.0)
    {
        metrics.areaweightednorthern5070simulatedwarmesttemperaturec =
            northernsimulatedwarmest / northerntemperatureweight;
        metrics.areaweightednorthern5070referencewarmesttemperaturec =
            northernreferencewarmest / northerntemperatureweight;
    }

    if (southernweight > 0.0)
    {
        metrics.areaweightedsouthern6090simulatedeffraction =
            southernsimulatedefweight / southernweight;
        metrics.areaweightedsouthern6090referenceeffraction =
            southernreferenceefweight / southernweight;
        metrics.areaweightedsouthern6090simulateddfraction =
            southernsimulateddweight / southernweight;
        metrics.areaweightedsouthern6090referencedfraction =
            southernreferencedweight / southernweight;
    }

    if (southerntemperatureweight > 0.0)
    {
        metrics.areaweightedsouthern6090simulatedwarmesttemperaturec =
            southernsimulatedwarmest / southerntemperatureweight;
        metrics.areaweightedsouthern6090referencewarmesttemperaturec =
            southernreferencewarmest / southerntemperatureweight;
    }
}

benchmarkphysicsmetrics collectbenchmarkphysicsmetrics(planet& world)
{
    benchmarkphysicsmetrics metrics;
    collectbenchmarktemperaturemetrics(world, metrics);
    const auto& precipitation = climatephysics::lastPrecipitationDistributionDiagnostics();
    const auto& budgets = climatephysics::lastAreaWeightedWaterBudgets();
    double oceanevaporation = 0.0;
    double landevaporation = 0.0;
    double oceanprecipitation = 0.0;
    double landprecipitation = 0.0;
    double runoff = 0.0;
    double atmosphericstorage = 0.0;

    for (const climatephysics::WaterBudget& budget : budgets)
    {
        oceanevaporation += budget.oceanEvaporation;
        landevaporation += budget.landEvaporation;
        oceanprecipitation += budget.oceanPrecipitation;
        landprecipitation += budget.landPrecipitation;
        runoff += budget.runoff;
        atmosphericstorage += budget.atmosphericStorage;
    }

    const double oceanimport = oceanevaporation - oceanprecipitation;
    metrics.areaweightedlandrawzeroprecipitationfraction =
        precipitation.land.areaWeightedRawZeroFraction;
    metrics.areaweightedlandstoredzeroprecipitationfraction =
        precipitation.land.areaWeightedStoredZeroFraction;
    metrics.areaweightedlandbelowonemillimetrefraction =
        precipitation.land.areaWeightedBelowOneMillimetreFraction;

    if (oceanevaporation > 0.0)
        metrics.areaweightedlandevaporationoceanevaporationratio = landevaporation / oceanevaporation;

    if (landprecipitation > 0.0)
        metrics.areaweightedrunofflandprecipitationratio = runoff / landprecipitation;

    if (oceanimport > 0.0)
        metrics.areaweightedlandprecipitationrecyclingratio = landprecipitation / oceanimport;

    const double totalprecipitation = oceanprecipitation + landprecipitation;
    constexpr double representativeperioddays = 365.0 / CLIMATESEASONCOUNT;

    if (totalprecipitation > 0.0)
    {
        metrics.areaweightedatmosphericwaterresidencetimedays =
            representativeperioddays * atmosphericstorage / totalprecipitation;
    }

    double griddedatmosphericstorage = 0.0;
    double griddedprecipitation = 0.0;

    for (int y = 0; y <= world.height(); y++)
    {
        const double areaweight = gridcellareaweight(y, world.height());

        for (int x = 0; x <= world.width(); x++)
        {
            for (int season = 0; season < CLIMATESEASONCOUNT; season++)
            {
                griddedatmosphericstorage +=
                    areaweight * world.seasonalmoisture(season, x, y);
                griddedprecipitation +=
                    areaweight * world.seasonalrainfloat(season, x, y);
            }
        }
    }

    if (griddedprecipitation > 0.0)
    {
        constexpr double griddedrepresentativeperioddays =
            365.0 / (CLIMATESEASONCOUNT * 3.0);
        metrics.areaweightedgriddedatmosphericwaterresidencetimedays =
            griddedrepresentativeperioddays * griddedatmosphericstorage /
                griddedprecipitation;
    }

    climatephysics::CondensationActivityDiagnostics condensation;

    for (const climatephysics::CondensationActivityDiagnostics& season :
        climatephysics::lastCondensationActivityDiagnostics())
    {
        condensation.cellStepAreaWeight += season.cellStepAreaWeight;
        condensation.activeCellStepAreaWeight += season.activeCellStepAreaWeight;
        condensation.atmosphericWaterAreaWeighted += season.atmosphericWaterAreaWeighted;
        condensation.activeAtmosphericWaterAreaWeighted +=
            season.activeAtmosphericWaterAreaWeighted;
        condensation.excessWaterAreaWeighted += season.excessWaterAreaWeighted;
        condensation.precipitationAreaWeighted += season.precipitationAreaWeighted;
    }

    if (condensation.cellStepAreaWeight > 0.0)
    {
        metrics.areaweightedactivecondensationcellstepfraction =
            condensation.activeCellStepAreaWeight / condensation.cellStepAreaWeight;
    }

    if (condensation.atmosphericWaterAreaWeighted > 0.0)
    {
        metrics.areaweightedactivecondensationwaterfraction =
            condensation.activeAtmosphericWaterAreaWeighted /
            condensation.atmosphericWaterAreaWeighted;
        metrics.areaweightedcondensableexcesswaterfraction =
            condensation.excessWaterAreaWeighted /
            condensation.atmosphericWaterAreaWeighted;
    }

    climatereference::MonthlyGrid imergprecipitation;
    string error;
    const filesystem::path referencepath =
        getappenvironment().referenceClimateDirectory / "imerg_prec_annual.uwclim";

    if (climatereference::loadMonthlyGrid(
        referencepath, "prec_annual", imergprecipitation, &error))
    {
        const comparisonmetrics comparison = comparemonthlyfield(
            world,
            imergprecipitation,
            -1,
            monthlyreferencefield::precipitation,
            referencescope::land,
            12.0);
        metrics.areaweightedlandprecipitationcorrelation =
            comparison.areaweightedcorrelation;
        metrics.areaweightedlandprecipitationcorrelation5degree =
            fivedegreelandprecipitationcorrelation(world, imergprecipitation);

        double simulatednorthernprecipitation = 0.0;
        double referencenorthernprecipitation = 0.0;

        for (int y = 0; y <= world.height(); y++)
        {
            const float latitude = latitudeforrow(y, world.height());

            if (latitude < 50.0f || latitude > 70.0f)
                continue;

            const double areaweight = gridcellareaweight(y, world.height());

            for (int x = 0; x <= world.width(); x++)
            {
                if (!isvalidationland(world, x, y))
                    continue;

                const double observed = imergprecipitation.value(0, x, y);

                if (!std::isfinite(observed))
                    continue;

                double simulated = 0.0;

                for (int season = 0; season < CLIMATESEASONCOUNT; season++)
                    simulated += world.seasonalrainfloat(season, x, y);

                simulated *= 12.0 / static_cast<double>(CLIMATESEASONCOUNT);
                simulatednorthernprecipitation += areaweight * simulated;
                referencenorthernprecipitation += areaweight * observed;
            }
        }

        if (referencenorthernprecipitation > 0.0)
        {
            metrics.areaweightednorthern5070precipitationratio =
                simulatednorthernprecipitation / referencenorthernprecipitation;
        }
    }

    auto compareera5 = [&](const char* filename,
                           const char* variable,
                           monthlyreferencefield field,
                           referencescope scope)
    {
        climatereference::MonthlyGrid reference;
        string loaderror;

        if (!climatereference::loadMonthlyGrid(
            getappenvironment().referenceClimateDirectory / filename,
            variable,
            reference,
            &loaderror))
        {
            return comparisonmetrics{};
        }

        return comparemonthlyfield(world, reference, -1, field, scope);
    };

    const comparisonmetrics pressurecomparison = compareera5(
        "era5_slp_anom_monthly.uwclim",
        "slp_anom",
        monthlyreferencefield::surfacepressure,
        referencescope::global);
    metrics.era5globalpressurecorrelation = pressurecomparison.areaweightedcorrelation;
    metrics.era5globalpressurermsehpa = pressurecomparison.areaweightedrmse;
    metrics.era5surfaceeastwardwindcorrelation = compareera5(
        "era5_u10m_monthly.uwclim",
        "u10m",
        monthlyreferencefield::surfaceuwind,
        referencescope::global).areaweightedcorrelation;
    metrics.era5surfacenorthwardwindcorrelation = compareera5(
        "era5_v10m_monthly.uwclim",
        "v10m",
        monthlyreferencefield::surfacevwind,
        referencescope::global).areaweightedcorrelation;
    metrics.era5transporteastwardwindcorrelation = compareera5(
        "era5_u850_monthly.uwclim",
        "u850",
        monthlyreferencefield::transportuwind,
        referencescope::global).areaweightedcorrelation;
    metrics.era5transportnorthwardwindcorrelation = compareera5(
        "era5_v850_monthly.uwclim",
        "v850",
        monthlyreferencefield::transportvwind,
        referencescope::global).areaweightedcorrelation;
    const comparisonmetrics columnwatercomparison = compareera5(
        "era5_tcwv_monthly.uwclim",
        "tcwv",
        monthlyreferencefield::columnwater,
        referencescope::global);
    metrics.era5columnwatercorrelation = columnwatercomparison.areaweightedcorrelation;
    metrics.era5columnwatermeanbiaskgm2 = columnwatercomparison.areaweightedmeanbias;
    const comparisonmetrics ascentcomparison = compareera5(
        "era5_w500_ascent_monthly.uwclim",
        "w500_ascent",
        monthlyreferencefield::verticalascent,
        referencescope::global);
    metrics.era5verticalascentcorrelation = ascentcomparison.areaweightedcorrelation;
    metrics.era5verticalascentrmsehpaday = ascentcomparison.areaweightedrmse;
    metrics.era5landprecipitationcorrelation = compareera5(
        "era5_pr_monthly.uwclim",
        "pr",
        monthlyreferencefield::precipitation,
        referencescope::land).areaweightedcorrelation;

    return metrics;
}

bool appendclimatebenchmarkrunlog(
    int runid,
    const string& timestamp,
    const string& information,
    double weightedrelativeerror,
    const climatespatialmetrics& spatial,
    const benchmarkphysicsmetrics& physics)
{
    const filesystem::path logpath = getappenvironment().climateBenchmarkRunLogPath;
    string content;

    {
        ifstream logfile(logpath);

        if (logfile.is_open())
            content.assign(istreambuf_iterator<char>(logfile), istreambuf_iterator<char>());
    }

    if (content.find_first_not_of(" \t\r\n") == string::npos)
        content = "{\n  \"runs\": []\n}\n";

    const size_t runskey = content.find("\"runs\"");
    const size_t arraybegin = runskey == string::npos ? string::npos : content.find('[', runskey);
    const size_t arrayend = arraybegin == string::npos ? string::npos : content.rfind(']');

    if (arraybegin == string::npos || arrayend == string::npos || arrayend < arraybegin)
        return false;

    const bool hasentries = content.find_first_not_of(" \t\r\n", arraybegin + 1) < arrayend;
    ostringstream entry;

    if (hasentries)
        entry << ',';

    entry
        << "\n    {\n"
        << "      \"id\": " << runid << ",\n"
        << "      \"datetime\": \"" << timestamp << "\",\n"
        << "      \"information\": \"" << jsonescape(information) << "\",\n"
        << "      \"diagnostics_directory\": \""
        << jsonescape((filesystem::path("extra") / "validation" / "runs" / to_string(runid)).generic_string())
        << "\",\n"
        << "      \"metrics\": {\n"
        << fixed << setprecision(10)
        << "        \"weighted_relative_error\": " << weightedrelativeerror << ",\n"
        << "        \"spatial_compared_cells\": " << spatial.comparedcells << ",\n"
        << "        \"spatial_exact_accuracy\": " << spatial.exactaccuracy << ",\n"
        << "        \"spatial_group_accuracy\": " << spatial.groupaccuracy << ",\n"
        << "        \"spatial_kappa\": " << spatial.kappa << ",\n"
        << "        \"area_weighted_spatial_exact_accuracy\": " << spatial.areaweightedexactaccuracy << ",\n"
        << "        \"area_weighted_spatial_group_accuracy\": " << spatial.areaweightedgroupaccuracy << ",\n"
        << "        \"area_weighted_spatial_kappa\": " << spatial.areaweightedkappa << ",\n"
        << "        \"area_weighted_land_raw_zero_precipitation_fraction\": "
        << physics.areaweightedlandrawzeroprecipitationfraction << ",\n"
        << "        \"area_weighted_land_stored_zero_precipitation_fraction\": "
        << physics.areaweightedlandstoredzeroprecipitationfraction << ",\n"
        << "        \"area_weighted_land_below_1mm_fraction\": "
        << physics.areaweightedlandbelowonemillimetrefraction << ",\n"
        << "        \"area_weighted_land_evaporation_ocean_evaporation_ratio\": "
        << physics.areaweightedlandevaporationoceanevaporationratio << ",\n"
        << "        \"area_weighted_runoff_land_precipitation_ratio\": "
        << physics.areaweightedrunofflandprecipitationratio << ",\n"
        << "        \"area_weighted_land_precipitation_recycling_ratio\": "
        << physics.areaweightedlandprecipitationrecyclingratio << ",\n"
        << "        \"area_weighted_land_precipitation_correlation\": "
        << physics.areaweightedlandprecipitationcorrelation << ",\n"
        << "        \"area_weighted_land_precipitation_correlation_5degree\": "
        << physics.areaweightedlandprecipitationcorrelation5degree << ",\n"
        << "        \"area_weighted_land_rainfall_storage_clipping_fraction\": "
        << physics.areaweightedlandrainfallstorageclippingfraction << ",\n"
        << "        \"area_weighted_atmospheric_water_residence_time_days\": "
        << physics.areaweightedatmosphericwaterresidencetimedays << ",\n"
        << "        \"area_weighted_budget_atmospheric_water_residence_time_days\": "
        << physics.areaweightedatmosphericwaterresidencetimedays << ",\n"
        << "        \"area_weighted_gridded_atmospheric_water_residence_time_days\": "
        << physics.areaweightedgriddedatmosphericwaterresidencetimedays << ",\n"
        << "        \"area_weighted_active_condensation_cell_step_fraction\": "
        << physics.areaweightedactivecondensationcellstepfraction << ",\n"
        << "        \"area_weighted_active_condensation_water_fraction\": "
        << physics.areaweightedactivecondensationwaterfraction << ",\n"
        << "        \"area_weighted_condensable_excess_water_fraction\": "
        << physics.areaweightedcondensableexcesswaterfraction << ",\n"
        << "        \"area_weighted_temperature_group_locked_fraction\": "
        << physics.areaweightedtemperaturegrouplockedfraction << ",\n"
        << "        \"area_weighted_thermal_signature_locked_fraction\": "
        << physics.areaweightedthermalsignaturelockedfraction << ",\n"
        << "        \"area_weighted_northern_50_70_temperature_group_locked_fraction\": "
        << physics.areaweightednorthern5070temperaturegrouplockedfraction << ",\n"
        << "        \"area_weighted_northern_50_70_thermal_signature_locked_fraction\": "
        << physics.areaweightednorthern5070thermalsignaturelockedfraction << ",\n"
        << "        \"area_weighted_northern_50_70_simulated_warmest_temperature_c\": "
        << physics.areaweightednorthern5070simulatedwarmesttemperaturec << ",\n"
        << "        \"area_weighted_northern_50_70_reference_warmest_temperature_c\": "
        << physics.areaweightednorthern5070referencewarmesttemperaturec << ",\n"
        << "        \"area_weighted_northern_50_70_precipitation_ratio\": "
        << physics.areaweightednorthern5070precipitationratio << ",\n"
        << "        \"area_weighted_southern_60_90_simulated_ef_fraction\": "
        << physics.areaweightedsouthern6090simulatedeffraction << ",\n"
        << "        \"area_weighted_southern_60_90_reference_ef_fraction\": "
        << physics.areaweightedsouthern6090referenceeffraction << ",\n"
        << "        \"area_weighted_southern_60_90_simulated_d_fraction\": "
        << physics.areaweightedsouthern6090simulateddfraction << ",\n"
        << "        \"area_weighted_southern_60_90_reference_d_fraction\": "
        << physics.areaweightedsouthern6090referencedfraction << ",\n"
        << "        \"area_weighted_southern_60_90_simulated_warmest_temperature_c\": "
        << physics.areaweightedsouthern6090simulatedwarmesttemperaturec << ",\n"
        << "        \"area_weighted_southern_60_90_reference_warmest_temperature_c\": "
        << physics.areaweightedsouthern6090referencewarmesttemperaturec << ",\n"
        << "        \"era5_global_pressure_correlation\": "
        << physics.era5globalpressurecorrelation << ",\n"
        << "        \"era5_global_pressure_rmse_hpa\": "
        << physics.era5globalpressurermsehpa << ",\n"
        << "        \"era5_surface_eastward_wind_correlation\": "
        << physics.era5surfaceeastwardwindcorrelation << ",\n"
        << "        \"era5_surface_northward_wind_correlation\": "
        << physics.era5surfacenorthwardwindcorrelation << ",\n"
        << "        \"era5_transport_eastward_wind_correlation\": "
        << physics.era5transporteastwardwindcorrelation << ",\n"
        << "        \"era5_transport_northward_wind_correlation\": "
        << physics.era5transportnorthwardwindcorrelation << ",\n"
        << "        \"era5_column_water_correlation\": "
        << physics.era5columnwatercorrelation << ",\n"
        << "        \"era5_column_water_mean_bias_kg_m2\": "
        << physics.era5columnwatermeanbiaskgm2 << ",\n"
        << "        \"era5_vertical_ascent_correlation\": "
        << physics.era5verticalascentcorrelation << ",\n"
        << "        \"era5_vertical_ascent_rmse_hpa_day\": "
        << physics.era5verticalascentrmsehpaday << ",\n"
        << "        \"era5_land_precipitation_correlation\": "
        << physics.era5landprecipitationcorrelation << "\n"
        << "      }\n"
        << "    }\n  ";

    content.insert(arrayend, entry.str());

    if (logpath.has_parent_path())
        filesystem::create_directories(logpath.parent_path());

    filesystem::path temppath = logpath;
    temppath += ".tmp";

    {
        ofstream tempfile(temppath, ios::trunc);

        if (!tempfile.is_open())
            return false;

        tempfile << content;
    }

    if (MoveFileExW(temppath.wstring().c_str(), logpath.wstring().c_str(), MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH) == FALSE)
    {
        error_code ignored;
        filesystem::remove(temppath, ignored);
        return false;
    }

    return true;
}

double saferelativeerror(long long simulated, long long reference)
{
    if (reference == 0)
        return simulated == 0 ? 0.0 : 1.0;

    return fabs(static_cast<double>(simulated - reference)) / fabs(static_cast<double>(reference));
}

bool runhiddenprocessandwait(const wstring& commandline)
{
    STARTUPINFOW startupinfo{};
    PROCESS_INFORMATION processinfo{};
    startupinfo.cb = sizeof(startupinfo);
    startupinfo.dwFlags = STARTF_USESHOWWINDOW;
    startupinfo.wShowWindow = SW_HIDE;

    vector<wchar_t> mutablecommand(commandline.begin(), commandline.end());
    mutablecommand.push_back(L'\0');

    if (CreateProcessW(nullptr, mutablecommand.data(), nullptr, nullptr, FALSE, CREATE_NO_WINDOW, nullptr, nullptr, &startupinfo, &processinfo) == FALSE)
        return false;

    WaitForSingleObject(processinfo.hProcess, INFINITE);

    DWORD exitcode = 1;
    GetExitCodeProcess(processinfo.hProcess, &exitcode);

    CloseHandle(processinfo.hThread);
    CloseHandle(processinfo.hProcess);

    return exitcode == 0;
}

bool updateclimatebenchmarkworkbook(
    int runid,
    const vector<string>& codes,
    const vector<long long>& simulationcounts,
    string* failuremessage)
{
    const filesystem::path workbookpath = filesystem::absolute(getappenvironment().climateWorkbookPath).lexically_normal();

    if (filesystem::exists(workbookpath) == false || simulationcounts.size() != codes.size())
    {
        if (failuremessage != nullptr)
            *failuremessage = "Climate benchmark workbook is missing or the climate columns do not match";
        return false;
    }

    const filesystem::path temproot = filesystem::temp_directory_path();
    const filesystem::path datapath = temproot / "uw_climate_benchmark_input.txt";
    const filesystem::path scriptpath = temproot / "uw_climate_benchmark_excel.ps1";
    const filesystem::path errorpath = temproot / "uw_climate_benchmark_excel_error.txt";

    {
        ofstream datafile(datapath);

        if (!datafile.is_open())
            return false;

        datafile << "ID";

        for (const string& code : codes)
            datafile << ',' << code;

        datafile << '\n';
        datafile << runid;

        for (const long long count : simulationcounts)
            datafile << ',' << count;

        datafile << '\n';
    }

    {
        ofstream scriptfile(scriptpath);

        if (!scriptfile.is_open())
            return false;

        scriptfile << "param([string]$WorkbookPath, [string]$DataPath, [string]$ErrorPath)\n";
        scriptfile << "$ErrorActionPreference = 'Stop'\n";
        scriptfile << "Remove-Item -LiteralPath $ErrorPath -ErrorAction SilentlyContinue\n";
        scriptfile << "$excel = $null\n";
        scriptfile << "$workbook = $null\n";
        scriptfile << "$sheet = $null\n";
        scriptfile << "try {\n";
        scriptfile << "$lines = Get-Content -Path $DataPath\n";
        scriptfile << "if ($lines.Count -lt 2) { throw 'Benchmark input is incomplete' }\n";
        scriptfile << "$headers = $lines[0].Split(',')\n";
        scriptfile << "$simulation = $lines[1].Split(',')\n";
        scriptfile << "$excel = New-Object -ComObject Excel.Application\n";
        scriptfile << "$excel.Visible = $false\n";
        scriptfile << "$excel.DisplayAlerts = $false\n";
        scriptfile << "$workbook = $excel.Workbooks.Open($WorkbookPath)\n";
        scriptfile << "if ($workbook.ReadOnly) { throw 'Climate benchmark workbook is open read-only; close it in Excel and retry' }\n";
        scriptfile << "$sheet = $workbook.Worksheets.Item('RAW_PIXELS')\n";
        scriptfile << "for ($index = 0; $index -lt $headers.Count; $index++) {\n";
        scriptfile << "    $actual = [string]$sheet.Cells.Item(1, $index + 1).Text\n";
        scriptfile << "    if ($actual -ne $headers[$index]) { throw \"Workbook header mismatch at column $($index + 1): expected '$($headers[$index])', found '$actual'\" }\n";
        scriptfile << "}\n";
        scriptfile << "$runId = [int]$simulation[0]\n";
        scriptfile << "if ($runId -lt 2) { throw 'Benchmark run ID must be at least 2' }\n";
        scriptfile << "$row = $runId + 2\n";
        scriptfile << "$existingId = $sheet.Cells.Item($row, 1).Value2\n";
        scriptfile << "if ($null -ne $existingId -and -not [string]::IsNullOrWhiteSpace([string]$existingId) -and [int]$existingId -ne $runId) { throw \"Workbook row $row already belongs to run $existingId\" }\n";
        scriptfile << "for ($column = $simulation.Count + 1; $column -le 34; $column++) {\n";
        scriptfile << "    if (-not $sheet.Cells.Item($row, $column).HasFormula) { throw \"RAW_PIXELS formula missing at row $row, column $column\" }\n";
        scriptfile << "}\n";
        scriptfile << "for ($index = 0; $index -lt $simulation.Count; $index++) {\n";
        scriptfile << "    $sheet.Cells.Item($row, $index + 1).Value2 = [double]$simulation[$index]\n";
        scriptfile << "}\n";
        scriptfile << "$workbook.Save()\n";
        scriptfile << "}\n";
        scriptfile << "catch {\n";
        scriptfile << "    Set-Content -LiteralPath $ErrorPath -Value $_.Exception.Message\n";
        scriptfile << "    exit 1\n";
        scriptfile << "}\n";
        scriptfile << "finally {\n";
        scriptfile << "    if ($null -ne $sheet) { [void][System.Runtime.Interopservices.Marshal]::ReleaseComObject($sheet) }\n";
        scriptfile << "    if ($null -ne $workbook) { try { $workbook.Close($false) } catch {}; [void][System.Runtime.Interopservices.Marshal]::ReleaseComObject($workbook) }\n";
        scriptfile << "    if ($null -ne $excel) { try { $excel.Quit() } catch {}; [void][System.Runtime.Interopservices.Marshal]::ReleaseComObject($excel) }\n";
        scriptfile << "    [GC]::Collect()\n";
        scriptfile << "    [GC]::WaitForPendingFinalizers()\n";
        scriptfile << "}\n";
    }

    wstring commandline = L"powershell.exe -NoProfile -ExecutionPolicy Bypass -File \"";
    commandline += scriptpath.wstring();
    commandline += L"\" -WorkbookPath \"";
    commandline += workbookpath.wstring();
    commandline += L"\" -DataPath \"";
    commandline += datapath.wstring();
    commandline += L"\" -ErrorPath \"";
    commandline += errorpath.wstring();
    commandline += L"\"";

    const bool succeeded = runhiddenprocessandwait(commandline);

    if (!succeeded)
    {
        ifstream errorfile(errorpath);
        string errormessage;

        if (getline(errorfile, errormessage) && !errormessage.empty())
        {
            cerr << "Climate benchmark workbook error: " << errormessage << '\n';
            if (failuremessage != nullptr)
                *failuremessage = errormessage;
        }
        else if (failuremessage != nullptr)
        {
            *failuremessage = "Excel benchmark update failed without a detailed error";
        }
    }

    return succeeded;
}

bool exportclimatebenchmarkimages(planet& world, int runid)
{
    const AppEnvironmentConfig& appenv = getappenvironment();
    const filesystem::path outputdir = appenv.climateBenchmarkImageDirectory;
    const filesystem::path referencepath = appenv.earthKoppenImagePath;
    const int width = world.width();
    const int height = world.height();

    if (outputdir.empty())
    {
        cerr << "Climate benchmark image directory is not configured.\n";
        return false;
    }

    error_code filesystemerror;
    filesystem::create_directories(outputdir, filesystemerror);

    if (filesystemerror)
    {
        cerr << "Failed to create climate benchmark image directory: " << filesystemerror.message() << '\n';
        return false;
    }

    const filesystem::path benchmarkreferencepath = outputdir / "0.png";
    bool writereference = filesystem::exists(benchmarkreferencepath) == false;

    if (!writereference)
    {
        sf::Image existingreference;
        writereference = existingreference.loadFromFile(benchmarkreferencepath.string()) == false ||
            existingreference.getSize().x != static_cast<unsigned int>(width + 1) ||
            existingreference.getSize().y != static_cast<unsigned int>(height + 1);

        if (!writereference)
        {
            const sf::Color efcolour = benchmarkclimatecolour(31);
            bool containsEf = false;

            for (unsigned int y = 0; y < existingreference.getSize().y && !containsEf; y++)
            {
                for (unsigned int x = 0; x < existingreference.getSize().x; x++)
                {
                    if (existingreference.getPixel(x, y) == efcolour)
                    {
                        containsEf = true;
                        break;
                    }
                }
            }

            writereference = !containsEf;
        }
    }

    if (writereference)
    {
        if (filesystem::exists(referencepath) == false)
        {
            cerr << "Climate benchmark reference image not found: " << referencepath.string() << '\n';
            return false;
        }

        sf::Image referenceimage;

        if (referenceimage.loadFromFile(referencepath.string()) == false)
        {
            cerr << "Failed to load climate benchmark reference image: " << referencepath.string() << '\n';
            return false;
        }

        const sf::Vector2u referencesize = referenceimage.getSize();

        for (unsigned int y = 0; y < referencesize.y; y++)
        {
            for (unsigned int x = 0; x < referencesize.x; x++)
            {
                int distancesquared = 0;
                const short climate = nearestbenchmarkclimate(referenceimage.getPixel(x, y), distancesquared);

                if (distancesquared <= benchmarkcolourdistancelimitsquared)
                    referenceimage.setPixel(x, y, benchmarkclimatecolour(climate));
                else
                    referenceimage.setPixel(x, y, sf::Color::Black);
            }
        }

        if (referenceimage.saveToFile(benchmarkreferencepath.string()) == false)
        {
            cerr << "Failed to save climate benchmark reference image: " << benchmarkreferencepath.string() << '\n';
            return false;
        }
    }

    sf::Image simulatedimage;
    simulatedimage.create(width + 1, height + 1, sf::Color::Black);

    for (int y = 0; y <= height; y++)
    {
        for (int x = 0; x <= width; x++)
        {
            const short climate = static_cast<short>(world.climate(x, y));
            const bool water = world.sea(x, y) == 1 || world.truelake(x, y) != 0 || world.riftlakesurface(x, y) != 0;

            if (water == false && climate >= 1 && climate <= 31)
                simulatedimage.setPixel(x, y, benchmarkclimatecolour(climate));
        }
    }

    const filesystem::path simulatedpath = outputdir / (to_string(runid) + ".png");

    if (simulatedimage.saveToFile(simulatedpath.string()) == false)
    {
        cerr << "Failed to save simulated climate image: " << simulatedpath.string() << '\n';
        return false;
    }

    sf::Image temperatureimage;
    sf::Image precipitationimage;
    temperatureimage.create(width + 1, height + 1, sf::Color::Black);
    precipitationimage.create(width + 1, height + 1, sf::Color::Black);
    vector<float> annualprecipitationvalues(
        static_cast<size_t>(width + 1) * static_cast<size_t>(height + 1));

    for (int y = 0; y <= height; y++)
    {
        for (int x = 0; x <= width; x++)
        {
            const bool water = world.sea(x, y) == 1 || world.truelake(x, y) != 0 ||
                world.riftlakesurface(x, y) != 0;

            if (!water)
            {
                temperatureimage.setPixel(
                    x,
                    y,
                    benchmarktemperaturecolour(static_cast<float>(world.avetemp(x, y))));
            }

            const float annualprecipitation =
                world.averainfloat(x, y) * 12.0f;
            annualprecipitationvalues[
                static_cast<size_t>(y) * static_cast<size_t>(width + 1) +
                static_cast<size_t>(x)] = annualprecipitation;
            precipitationimage.setPixel(
                x,
                y,
                benchmarkprecipitationcolour(annualprecipitation));
        }
    }

    const filesystem::path temperaturepath =
        outputdir / (to_string(runid) + "_temperature_c.png");
    const filesystem::path precipitationpreviewpath =
        outputdir / (to_string(runid) + "_precipitation_mm_year_preview.png");
    const filesystem::path precipitationtiffpath =
        outputdir / (to_string(runid) + "_precipitation_mm_year.tif");

    if (temperatureimage.saveToFile(temperaturepath.string()) == false)
    {
        cerr << "Failed to save simulated temperature image: " << temperaturepath.string() << '\n';
        return false;
    }

    if (precipitationimage.saveToFile(precipitationpreviewpath.string()) == false)
    {
        cerr << "Failed to save simulated precipitation preview: " << precipitationpreviewpath.string() << '\n';
        return false;
    }

    if (climateio::writefloat32geotiff(
            precipitationtiffpath.string().c_str(),
            static_cast<uint32_t>(width + 1),
            static_cast<uint32_t>(height + 1),
            annualprecipitationvalues.data()) == false)
    {
        cerr << "Failed to save simulated precipitation GeoTIFF: " << precipitationtiffpath.string() << '\n';
        return false;
    }

    ofstream scalefile(outputdir / "physical_map_scales.txt");

    if (!scalefile.is_open())
    {
        cerr << "Failed to save benchmark physical-map scale description.\n";
        return false;
    }

    scalefile << "temperature_file={run_id}_temperature_c.png\n";
    scalefile << "temperature_quantity=four-season mean surface air temperature\n";
    scalefile << "temperature_units=degrees C\n";
    scalefile << "temperature_range=-54 to 30 (values outside the range are clamped)\n";
    scalefile << "temperature_water=black\n";
    scalefile << "temperature_reference=extra/img/earth/in/earth_temp_l.png\n";
    scalefile << "precipitation_file={run_id}_precipitation_mm_year.tif\n";
    scalefile << "precipitation_preview={run_id}_precipitation_mm_year_preview.png\n";
    scalefile << "precipitation_quantity=four-season mean monthly precipitation multiplied by 12\n";
    scalefile << "precipitation_units=mm/year\n";
    scalefile << "precipitation_sample=single-band IEEE float32 (unclamped physical values)\n";
    scalefile << "precipitation_nodata=-9999.9\n";
    scalefile << "precipitation_georeference=WGS84 EPSG:4326, global equirectangular, north-to-south\n";
    scalefile << "precipitation_preview_range=0 to 6000 mm/year (values above the range are clamped)\n";
    scalefile << "precipitation_water=included\n";
    scalefile << "precipitation_reference=extra/reference/imerg_precipitation_mm_year.tif\n";

    return true;
}

bool loadprecipitationgrid(const filesystem::path& filepath, vector<vector<double>>& grid)
{
    ifstream infile(filepath);

    if (!infile.is_open())
        return false;

    string line;

    if (!getline(infile, line))
        return false;

    while (getline(infile, line))
    {
        if (line.empty())
            continue;

        vector<double> rowvalues;
        string token;
        stringstream linestream(line);
        int column = 0;

        while (getline(linestream, token, ','))
        {
            if (column >= 2)
            {
                rowvalues.push_back(token.empty() ?
                    numeric_limits<double>::quiet_NaN() :
                    stod(token));
            }

            column++;
        }

        if (rowvalues.empty() == false)
            grid.push_back(rowvalues);
    }

    return grid.empty() == false;
}

comparisonmetrics compareannualprecipitation(const filesystem::path& outputdir, planet& world, const vector<zonalstats>& rows)
{
    comparisonmetrics metrics;
    const filesystem::path referencepath = getappenvironment().referencePrecipitationGridPath;
    const filesystem::path comparisonpath = outputdir / "annual_precipitation_comparison.txt";

    if (filesystem::exists(referencepath) == false)
    {
        ofstream comparisonfile(comparisonpath);

        if (comparisonfile.is_open())
            comparisonfile << "status=reference_not_found\nreference_grid_path=" << referencepath.string() << '\n';

        return metrics;
    }

    metrics.referencefound = true;

    vector<vector<double>> referencegrid;

    if (loadprecipitationgrid(referencepath, referencegrid) == false)
    {
        ofstream comparisonfile(comparisonpath);

        if (comparisonfile.is_open())
            comparisonfile << "status=reference_unreadable\nreference_grid_path=" << referencepath.string() << '\n';

        return metrics;
    }

    const int width = world.width();
    const int height = world.height();

    if (static_cast<int>(referencegrid.size()) != height + 1)
    {
        ofstream comparisonfile(comparisonpath);

        if (comparisonfile.is_open())
        {
            comparisonfile << "status=dimension_mismatch\n";
            comparisonfile << "reference_grid_path=" << referencepath.string() << '\n';
            comparisonfile << "expected_height=" << height + 1 << '\n';
            comparisonfile << "actual_height=" << referencegrid.size() << '\n';
        }

        return metrics;
    }

    for (const auto& row : referencegrid)
    {
        if (static_cast<int>(row.size()) != width + 1)
        {
            ofstream comparisonfile(comparisonpath);

            if (comparisonfile.is_open())
            {
                comparisonfile << "status=dimension_mismatch\n";
                comparisonfile << "reference_grid_path=" << referencepath.string() << '\n';
                comparisonfile << "expected_width=" << width + 1 << '\n';
                comparisonfile << "actual_width=" << row.size() << '\n';
            }

            return metrics;
        }
    }

    metrics.dimensionsmatch = true;

    double simulatedsum = 0.0;
    double referencesum = 0.0;
    double biassum = 0.0;
    double absoluteerrorsum = 0.0;
    double squarederrorsum = 0.0;
    double sumsim2 = 0.0;
    double sumref2 = 0.0;
    double sumcross = 0.0;
    double tropicalbiassum = 0.0;
    int tropicalcells = 0;
    double areaweightedsimulatedsum = 0.0;
    double areaweightedreferencesum = 0.0;
    double areaweightedbiassum = 0.0;
    double areaweightedabsoluteerrorsum = 0.0;
    double areaweightedsquarederrorsum = 0.0;
    double areaweightedsumsimulatedsquared = 0.0;
    double areaweightedsumreferencesquared = 0.0;
    double areaweightedsumcross = 0.0;

    ofstream zonalcomparisonfile(outputdir / "annual_precipitation_zonal_comparison.csv");

    if (zonalcomparisonfile.is_open())
    {
        zonalcomparisonfile << "y,latitude,simulated_mean_annual_rain,reference_mean_annual_rain,bias\n";
        zonalcomparisonfile << fixed << setprecision(4);
    }

    for (int y = 0; y <= height; y++)
    {
        const double areaweight = gridcellareaweight(y, height);
        double simulatedsumrow = 0.0;
        double referencesumrow = 0.0;
        int referencecellsrow = 0;

        for (int x = 0; x <= width; x++)
        {
            if (!isvalidationland(world, x, y))
                continue;

            const double simulated = static_cast<double>(world.averainfloat(x, y));
            const double reference = referencegrid[y][x];

            if (!std::isfinite(reference))
                continue;

            const double diff = simulated - reference;

            metrics.comparedcells++;
            simulatedsum = simulatedsum + simulated;
            referencesum = referencesum + reference;
            biassum = biassum + diff;
            absoluteerrorsum = absoluteerrorsum + fabs(diff);
            squarederrorsum = squarederrorsum + diff * diff;
            sumsim2 = sumsim2 + simulated * simulated;
            sumref2 = sumref2 + reference * reference;
            sumcross = sumcross + simulated * reference;
            simulatedsumrow = simulatedsumrow + simulated;
            referencesumrow = referencesumrow + reference;
            referencecellsrow++;
            metrics.areaweight += areaweight;
            areaweightedsimulatedsum += areaweight * simulated;
            areaweightedreferencesum += areaweight * reference;
            areaweightedbiassum += areaweight * diff;
            areaweightedabsoluteerrorsum += areaweight * fabs(diff);
            areaweightedsquarederrorsum += areaweight * diff * diff;
            areaweightedsumsimulatedsquared += areaweight * simulated * simulated;
            areaweightedsumreferencesquared += areaweight * reference * reference;
            areaweightedsumcross += areaweight * simulated * reference;

            if (fabs(latitudeforrow(y, height)) <= 30.0f)
            {
                tropicalbiassum = tropicalbiassum + diff;
                tropicalcells++;
            }
        }

        if (zonalcomparisonfile.is_open())
        {
            const double simulatedmeanrow = safeaverage(simulatedsumrow, referencecellsrow);
            const double referencemeanrow = safeaverage(referencesumrow, referencecellsrow);

            zonalcomparisonfile
                << y << ','
                << latitudeforrow(y, height) << ','
                << simulatedmeanrow << ','
                << referencemeanrow << ','
                << simulatedmeanrow - referencemeanrow << '\n';
        }
    }

    metrics.simulatedmean = safeaverage(simulatedsum, metrics.comparedcells);
    metrics.referencemean = safeaverage(referencesum, metrics.comparedcells);
    metrics.meanbias = safeaverage(biassum, metrics.comparedcells);
    metrics.meanabsoluteerror = safeaverage(absoluteerrorsum, metrics.comparedcells);
    metrics.rmse = sqrt(safeaverage(squarederrorsum, metrics.comparedcells));
    metrics.tropicalmeanbias = safeaverage(tropicalbiassum, tropicalcells);

    const double numerator = sumcross - (simulatedsum * referencesum / static_cast<double>(metrics.comparedcells));
    const double simulatedvariance = sumsim2 - (simulatedsum * simulatedsum / static_cast<double>(metrics.comparedcells));
    const double referencevariance = sumref2 - (referencesum * referencesum / static_cast<double>(metrics.comparedcells));
    const double denominator = sqrt(max(0.0, simulatedvariance) * max(0.0, referencevariance));

    if (denominator > 0.0)
        metrics.correlation = numerator / denominator;

    if (metrics.areaweight > 0.0)
    {
        metrics.areaweightedsimulatedmean = areaweightedsimulatedsum / metrics.areaweight;
        metrics.areaweightedreferencemean = areaweightedreferencesum / metrics.areaweight;
        metrics.areaweightedmeanbias = areaweightedbiassum / metrics.areaweight;
        metrics.areaweightedmeanabsoluteerror = areaweightedabsoluteerrorsum / metrics.areaweight;
        metrics.areaweightedrmse = sqrt(areaweightedsquarederrorsum / metrics.areaweight);
        const double areaweightednumerator = areaweightedsumcross -
            areaweightedsimulatedsum * areaweightedreferencesum / metrics.areaweight;
        const double areaweightedsimulatedvariance = areaweightedsumsimulatedsquared -
            areaweightedsimulatedsum * areaweightedsimulatedsum / metrics.areaweight;
        const double areaweightedreferencevariance = areaweightedsumreferencesquared -
            areaweightedreferencesum * areaweightedreferencesum / metrics.areaweight;
        const double areaweighteddenominator = sqrt(
            max(0.0, areaweightedsimulatedvariance) * max(0.0, areaweightedreferencevariance));

        if (areaweighteddenominator > 0.0)
            metrics.areaweightedcorrelation = areaweightednumerator / areaweighteddenominator;
    }

    ofstream comparisonfile(comparisonpath);

    if (comparisonfile.is_open())
    {
        comparisonfile << "status=ok\n";
        comparisonfile << fixed << setprecision(6);
        comparisonfile << "reference_grid_path=" << referencepath.string() << '\n';
        comparisonfile << "compared_cells=" << metrics.comparedcells << '\n';
        comparisonfile << "simulated_mean=" << metrics.simulatedmean << '\n';
        comparisonfile << "reference_mean=" << metrics.referencemean << '\n';
        comparisonfile << "mean_bias=" << metrics.meanbias << '\n';
        comparisonfile << "mean_absolute_error=" << metrics.meanabsoluteerror << '\n';
        comparisonfile << "rmse=" << metrics.rmse << '\n';
        comparisonfile << "correlation=" << metrics.correlation << '\n';
        comparisonfile << "tropical_mean_bias=" << metrics.tropicalmeanbias << '\n';
        comparisonfile << "area_weight=" << metrics.areaweight << '\n';
        comparisonfile << "area_weighted_simulated_mean=" << metrics.areaweightedsimulatedmean << '\n';
        comparisonfile << "area_weighted_reference_mean=" << metrics.areaweightedreferencemean << '\n';
        comparisonfile << "area_weighted_mean_bias=" << metrics.areaweightedmeanbias << '\n';
        comparisonfile << "area_weighted_mae=" << metrics.areaweightedmeanabsoluteerror << '\n';
        comparisonfile << "area_weighted_rmse=" << metrics.areaweightedrmse << '\n';
        comparisonfile << "area_weighted_correlation=" << metrics.areaweightedcorrelation << '\n';
    }

    return metrics;
}

void writeatmosphericbudget(const filesystem::path& outputdir, planet& world)
{
    constexpr array<const char*, CLIMATESEASONCOUNT> seasonnames = { "january", "april", "july", "october" };
    constexpr double verticalstoragescale = tuning::climate::circulation::verticalVelocityStorageScale;
    const int width = world.width();
    const int height = world.height();
    ofstream outfile(outputdir / "climate_atmosphere_budget.csv");

    if (!outfile.is_open())
        return;

    outfile << "season,area_weighted_mean_pressure_anomaly_hpa,rms_pressure_anomaly_hpa,rms_zonal_mean_pressure_anomaly_hpa,rms_stationary_pressure_anomaly_hpa,pressure_limiter_fraction,area_weighted_mean_surface_divergence_per_day,area_weighted_mean_upper_divergence_per_day,area_weighted_mean_vertical_velocity_hpa_per_day,rms_surface_wind_m_s,surface_wind_limiter_fraction,rms_upper_wind_m_s,rms_vertical_velocity_hpa_per_day,vertical_limiter_fraction,rms_column_divergence_per_day,mean_orographic_uplift_m,maximum_orographic_uplift_m,uplift_above_500m_fraction,mean_orographic_descent_m,maximum_orographic_descent_m,descent_above_500m_fraction\n";
    outfile << fixed << setprecision(9);

    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
    {
        double weighttotal = 0.0;
        double pressuretotal = 0.0;
        double pressuresquared = 0.0;
        double zonalpressuresquared = 0.0;
        double stationarypressuresquared = 0.0;
        double pressurelimitedweight = 0.0;
        double surfacedivergencetotal = 0.0;
        double upperdivergencetotal = 0.0;
        double verticaltotal = 0.0;
        double surfacewindsquared = 0.0;
        double upperwindsquared = 0.0;
        double verticalsquared = 0.0;
        double columndivergencesquared = 0.0;
        double surfacewindlimitedweight = 0.0;
        double verticallimitedweight = 0.0;
        double upliftmetrestotal = 0.0;
        double descentmetrestotal = 0.0;
        double upliftabove500weight = 0.0;
        double descentabove500weight = 0.0;
        double maximumupliftmetres = 0.0;
        double maximumdescentmetres = 0.0;
        vector<double> zonalpressure(height + 1, 0.0);

        for (int y = 0; y <= height; y++)
        {
            for (int x = 0; x <= width; x++)
                zonalpressure[y] += world.seasonalpressure(season, x, y);

            zonalpressure[y] /= static_cast<double>(width + 1);
        }

        for (int y = 0; y <= height; y++)
        {
            const int ynorth = y > 0 ? y - 1 : y;
            const int ysouth = y < height ? y + 1 : y;
            const double latitude = latitudeforrow(y, height);
            const double centrecosine = (std::max)(
                0.02,
                std::abs(std::cos(latitude * 3.14159265358979323846 / 180.0)));
            const double northcosine = (std::max)(
                0.0,
                std::cos(latitudeforrow(ynorth, height) * 3.14159265358979323846 / 180.0));
            const double southcosine = (std::max)(
                0.0,
                std::cos(latitudeforrow(ysouth, height) * 3.14159265358979323846 / 180.0));
            const double weight = gridcellareaweight(y, height);
            const auto spacing = climateatmosphere::cellSpacingMetres(
                latitude,
                width + 1,
                height + 1,
                tuning::climate::atmosphere::referencePlanetRadiusMetres);

            for (int x = 0; x <= width; x++)
            {
                const int xwest = x > 0 ? x - 1 : width - 1;
                const int xeast = x < width ? x + 1 : 1;
                const double surfaceu = world.seasonaluwind(season, x, y);
                const double surfacev = world.seasonalvwind(season, x, y);
                const double upperu = world.seasonalupperuwind(season, x, y);
                const double upperv = world.seasonaluppervwind(season, x, y);
                const double surfacedivergence =
                    ((world.seasonaluwind(season, xeast, y) - world.seasonaluwind(season, xwest, y)) /
                        (2.0 * spacing.zonalMetres) +
                        (world.seasonalvwind(season, x, ysouth) * southcosine -
                            world.seasonalvwind(season, x, ynorth) * northcosine) /
                        (2.0 * spacing.meridionalMetres * centrecosine)) *
                    tuning::climate::circulation::secondsPerDay;
                const double upperdivergence =
                    ((world.seasonalupperuwind(season, xeast, y) - world.seasonalupperuwind(season, xwest, y)) /
                        (2.0 * spacing.zonalMetres) +
                        (world.seasonaluppervwind(season, x, ysouth) * southcosine -
                            world.seasonaluppervwind(season, x, ynorth) * northcosine) /
                        (2.0 * spacing.meridionalMetres * centrecosine)) *
                    tuning::climate::circulation::secondsPerDay;
                const double verticalvelocity =
                    static_cast<double>(world.seasonalverticalvelocity(season, x, y)) / verticalstoragescale;
                const double columndivergence = surfacedivergence + upperdivergence;
                const double upliftmetres = static_cast<double>(
                    world.seasonaluplift(season, x, y)) /
                    tuning::climate::atmosphere::topographyVerticalMotionStorageScale;
                const double descentmetres = static_cast<double>(
                    world.seasonalsubsidence(season, x, y)) /
                    tuning::climate::atmosphere::topographyVerticalMotionStorageScale;

                weighttotal += weight;
                pressuretotal += weight * world.seasonalpressure(season, x, y);
                pressuresquared += weight * world.seasonalpressure(season, x, y) *
                    world.seasonalpressure(season, x, y);
                zonalpressuresquared += weight * zonalpressure[y] * zonalpressure[y];
                const double stationarypressure =
                    world.seasonalpressure(season, x, y) - zonalpressure[y];
                stationarypressuresquared +=
                    weight * stationarypressure * stationarypressure;
                surfacedivergencetotal += weight * surfacedivergence;
                upperdivergencetotal += weight * upperdivergence;
                verticaltotal += weight * verticalvelocity;
                surfacewindsquared += weight * (surfaceu * surfaceu + surfacev * surfacev);
                upperwindsquared += weight * (upperu * upperu + upperv * upperv);
                verticalsquared += weight * verticalvelocity * verticalvelocity;
                columndivergencesquared += weight * columndivergence * columndivergence;
                upliftmetrestotal += weight * upliftmetres;
                descentmetrestotal += weight * descentmetres;
                maximumupliftmetres = (std::max)(maximumupliftmetres, upliftmetres);
                maximumdescentmetres = (std::max)(maximumdescentmetres, descentmetres);
                upliftabove500weight += upliftmetres >= 500.0 ? weight : 0.0;
                descentabove500weight += descentmetres >= 500.0 ? weight : 0.0;

                if (std::abs(world.seasonalpressure(season, x, y)) >=
                    tuning::climate::circulation::maximumSurfacePressureAnomalyHpa - 1.0f)
                {
                    pressurelimitedweight += weight;
                }

                if (std::sqrt(surfaceu * surfaceu + surfacev * surfacev) >=
                    tuning::climate::atmosphere::maxVectorWind - 1.0f)
                {
                    surfacewindlimitedweight += weight;
                }

                if (std::abs(verticalvelocity) >=
                    tuning::climate::circulation::maximumVerticalVelocity - 1.0f)
                {
                    verticallimitedweight += weight;
                }
            }
        }

        const double inverseweight = weighttotal > 0.0 ? 1.0 / weighttotal : 0.0;
        outfile
            << seasonnames[season] << ','
            << pressuretotal * inverseweight << ','
            << sqrt(pressuresquared * inverseweight) << ','
            << sqrt(zonalpressuresquared * inverseweight) << ','
            << sqrt(stationarypressuresquared * inverseweight) << ','
            << pressurelimitedweight * inverseweight << ','
            << surfacedivergencetotal * inverseweight << ','
            << upperdivergencetotal * inverseweight << ','
            << verticaltotal * inverseweight << ','
            << sqrt(surfacewindsquared * inverseweight) << ','
            << surfacewindlimitedweight * inverseweight << ','
            << sqrt(upperwindsquared * inverseweight) << ','
            << sqrt(verticalsquared * inverseweight) << ','
            << verticallimitedweight * inverseweight << ','
            << sqrt(columndivergencesquared * inverseweight) << ','
            << upliftmetrestotal * inverseweight << ','
            << maximumupliftmetres << ','
            << upliftabove500weight * inverseweight << ','
            << descentmetrestotal * inverseweight << ','
            << maximumdescentmetres << ','
            << descentabove500weight * inverseweight << '\n';
    }
}

void writecirculationprecision(const filesystem::path& outputdir)
{
    constexpr array<const char*, CLIMATESEASONCOUNT> seasonnames = {
        "january", "april", "july", "october"
    };
    ofstream output(outputdir / "climate_circulation_precision.csv");

    if (!output.is_open())
        return;

    output << "season,stage,float_pressure_gradient_rms_pa_m,rounded_pressure_gradient_rms_pa_m,pressure_gradient_difference_rms_pa_m,float_surface_wind_rms_m_s,rounded_surface_wind_rms_m_s,surface_wind_difference_rms_m_s,float_upper_height_gradient_rms,rounded_upper_height_gradient_rms,upper_height_gradient_difference_rms,float_upper_wind_rms_m_s,rounded_upper_wind_rms_m_s,upper_wind_difference_rms_m_s\n";
    output << fixed << setprecision(9);
    const auto& diagnostics = climateatmosphere::lastCirculationPrecisionDiagnostics();

    auto writerow = [&](int season, const char* stage,
        const climateatmosphere::CirculationPrecisionStageDiagnostics& values)
    {
        output
            << seasonnames[season] << ',' << stage << ','
            << values.areaWeightedFloatPressureGradientRmsPaPerMetre << ','
            << values.areaWeightedRoundedPressureGradientRmsPaPerMetre << ','
            << values.areaWeightedPressureGradientDifferenceRmsPaPerMetre << ','
            << values.areaWeightedFloatSurfaceWindRmsMetresPerSecond << ','
            << values.areaWeightedRoundedSurfaceWindRmsMetresPerSecond << ','
            << values.areaWeightedSurfaceWindDifferenceRmsMetresPerSecond << ','
            << values.areaWeightedFloatUpperHeightGradientRms << ','
            << values.areaWeightedRoundedUpperHeightGradientRms << ','
            << values.areaWeightedUpperHeightGradientDifferenceRms << ','
            << values.areaWeightedFloatUpperWindRmsMetresPerSecond << ','
            << values.areaWeightedRoundedUpperWindRmsMetresPerSecond << ','
            << values.areaWeightedUpperWindDifferenceRmsMetresPerSecond << '\n';
    };

    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
    {
        writerow(season, "base", diagnostics[season].base);
        writerow(season, "final", diagnostics[season].final);
    }
}

void writepressuredecompositioncomparison(
    const filesystem::path& outputdir,
    planet& world)
{
    climatereference::MonthlyGrid reference;
    string error;
    const bool loaded = climatereference::loadMonthlyGrid(
        getappenvironment().referenceClimateDirectory / "era5_slp_anom_monthly.uwclim",
        "slp_anom",
        reference,
        &error);
    ofstream output(outputdir / "pressure_decomposition_comparison.csv");

    if (!output.is_open())
        return;

    if (!loaded || reference.width != world.width() + 1 ||
        reference.height != world.height() + 1)
    {
        output << "status,error\nerror," << csvescape(error) << '\n';
        return;
    }

    constexpr array<const char*, CLIMATESEASONCOUNT> seasonnames = {
        "january", "april", "july", "october"
    };
    constexpr array<int, CLIMATESEASONCOUNT> referencemonths = { 0, 3, 6, 9 };
    const int width = world.width();
    const int height = world.height();
    output << "season,simulated_total_rms_hpa,reference_total_rms_hpa,total_correlation,simulated_zonal_rms_hpa,reference_zonal_rms_hpa,zonal_correlation,simulated_stationary_rms_hpa,reference_stationary_rms_hpa,stationary_correlation\n";
    output << fixed << setprecision(9);

    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
    {
        vector<double> simulatedzonal(height + 1, 0.0);
        vector<double> referencezonal(height + 1, 0.0);

        for (int y = 0; y <= height; y++)
        {
            int referencecount = 0;

            for (int x = 0; x <= width; x++)
            {
                simulatedzonal[y] += world.seasonalpressure(season, x, y);
                const double referencevalue =
                    reference.value(referencemonths[season], x, y);

                if (std::isfinite(referencevalue))
                {
                    referencezonal[y] += referencevalue;
                    referencecount++;
                }
            }

            simulatedzonal[y] /= static_cast<double>(width + 1);

            if (referencecount > 0)
                referencezonal[y] /= static_cast<double>(referencecount);
        }

        weightedcorrelationmoments totalmoments;
        weightedcorrelationmoments zonalmoments;
        weightedcorrelationmoments stationarymoments;
        double zonalweight = 0.0;
        double simulatedzonalsquared = 0.0;
        double referencezonalsquared = 0.0;
        double cellweight = 0.0;
        double simulatedtotalsquared = 0.0;
        double referencetotalsquared = 0.0;
        double simulatedstationarysquared = 0.0;
        double referencestationarysquared = 0.0;

        for (int y = 0; y <= height; y++)
        {
            const double weight = gridcellareaweight(y, height);
            zonalweight += weight;
            simulatedzonalsquared += weight * simulatedzonal[y] * simulatedzonal[y];
            referencezonalsquared += weight * referencezonal[y] * referencezonal[y];
            zonalmoments.add(simulatedzonal[y], referencezonal[y], weight);

            for (int x = 0; x <= width; x++)
            {
                const double referencevalue =
                    reference.value(referencemonths[season], x, y);

                if (!std::isfinite(referencevalue))
                    continue;

                const double simulatedvalue = world.seasonalpressure(season, x, y);
                const double simulatedstationary = simulatedvalue - simulatedzonal[y];
                const double referencestationary = referencevalue - referencezonal[y];
                cellweight += weight;
                simulatedtotalsquared += weight * simulatedvalue * simulatedvalue;
                referencetotalsquared += weight * referencevalue * referencevalue;
                simulatedstationarysquared +=
                    weight * simulatedstationary * simulatedstationary;
                referencestationarysquared +=
                    weight * referencestationary * referencestationary;
                totalmoments.add(simulatedvalue, referencevalue, weight);
                stationarymoments.add(
                    simulatedstationary, referencestationary, weight);
            }
        }

        output
            << seasonnames[season] << ','
            << (cellweight > 0.0 ? sqrt(simulatedtotalsquared / cellweight) : 0.0) << ','
            << (cellweight > 0.0 ? sqrt(referencetotalsquared / cellweight) : 0.0) << ','
            << totalmoments.correlation() << ','
            << (zonalweight > 0.0 ? sqrt(simulatedzonalsquared / zonalweight) : 0.0) << ','
            << (zonalweight > 0.0 ? sqrt(referencezonalsquared / zonalweight) : 0.0) << ','
            << zonalmoments.correlation() << ','
            << (cellweight > 0.0 ? sqrt(simulatedstationarysquared / cellweight) : 0.0) << ','
            << (cellweight > 0.0 ? sqrt(referencestationarysquared / cellweight) : 0.0) << ','
            << stationarymoments.correlation() << '\n';
    }
}

void writetemperaturethresholdcomparison(
    const filesystem::path& outputdir,
    planet& world)
{
    climatereference::MonthlyGrid reference;
    string error;
    const bool loaded = climatereference::loadMonthlyGrid(
        getappenvironment().referenceClimateDirectory / "worldclim_tavg_monthly.uwclim",
        "tavg",
        reference,
        &error);
    ofstream output(outputdir / "temperature_threshold_comparison.csv");

    if (!output.is_open())
        return;

    if (!loaded || reference.width != world.width() + 1 ||
        reference.height != world.height() + 1)
    {
        output << "status,error\nerror," << csvescape(error) << '\n';
        return;
    }

    struct thresholdspecification
    {
        const char* name;
        bool warmest;
        double threshold;
        bool trueabove;
        bool inclusive;
    };

    constexpr array<thresholdspecification, 5> specifications = {
        thresholdspecification{ "polar_warmest_le_10", true, 10.0, false, true },
        { "icecap_warmest_lt_0", true, 0.0, false, false },
        { "hot_summer_warmest_ge_22", true, 22.0, true, true },
        { "continental_coldest_le_neg3", false, -3.0, false, true },
        { "tropical_coldest_ge_18", false, 18.0, true, true }
    };
    output << "threshold,area_weight,simulated_true_fraction,reference_true_fraction,agreement_fraction,false_positive_fraction,false_negative_fraction\n";
    output << fixed << setprecision(9);

    for (const thresholdspecification& specification : specifications)
    {
        double areaweight = 0.0;
        double simulatedtrueweight = 0.0;
        double referencetrueweight = 0.0;
        double agreementweight = 0.0;
        double falsepositiveweight = 0.0;
        double falsenegativeweight = 0.0;

        auto evaluate = [&](double value)
        {
            if (specification.trueabove)
            {
                return specification.inclusive ?
                    value >= specification.threshold : value > specification.threshold;
            }

            return specification.inclusive ?
                value <= specification.threshold : value < specification.threshold;
        };

        for (int y = 0; y <= world.height(); y++)
        {
            const double weight = gridcellareaweight(y, world.height());

            for (int x = 0; x <= world.width(); x++)
            {
                if (!isvalidationland(world, x, y))
                    continue;

                double simulatedminimum = numeric_limits<double>::infinity();
                double simulatedmaximum = -numeric_limits<double>::infinity();
                double referenceminimum = numeric_limits<double>::infinity();
                double referencemaximum = -numeric_limits<double>::infinity();

                for (int season = 0; season < CLIMATESEASONCOUNT; season++)
                {
                    const double value = world.seasonaltemp(season, x, y);
                    simulatedminimum = min(simulatedminimum, value);
                    simulatedmaximum = max(simulatedmaximum, value);
                }

                for (int month = 0; month < reference.monthCount; month++)
                {
                    const double value = reference.value(month, x, y);

                    if (!std::isfinite(value))
                        continue;

                    referenceminimum = min(referenceminimum, value);
                    referencemaximum = max(referencemaximum, value);
                }

                if (!std::isfinite(referenceminimum) ||
                    !std::isfinite(referencemaximum))
                {
                    continue;
                }

                const bool simulated = evaluate(
                    specification.warmest ? simulatedmaximum : simulatedminimum);
                const bool observed = evaluate(
                    specification.warmest ? referencemaximum : referenceminimum);
                areaweight += weight;

                if (simulated)
                    simulatedtrueweight += weight;

                if (observed)
                    referencetrueweight += weight;

                if (simulated == observed)
                    agreementweight += weight;
                else if (simulated)
                    falsepositiveweight += weight;
                else
                    falsenegativeweight += weight;
            }
        }

        output
            << specification.name << ','
            << areaweight << ','
            << (areaweight > 0.0 ? simulatedtrueweight / areaweight : 0.0) << ','
            << (areaweight > 0.0 ? referencetrueweight / areaweight : 0.0) << ','
            << (areaweight > 0.0 ? agreementweight / areaweight : 0.0) << ','
            << (areaweight > 0.0 ? falsepositiveweight / areaweight : 0.0) << ','
            << (areaweight > 0.0 ? falsenegativeweight / areaweight : 0.0) << '\n';
    }
}

void writeprecipitationgrid(const filesystem::path& filepath, planet& world, int mode)
{
    ofstream outfile(filepath);

    if (!outfile.is_open())
        return;

    const int width = world.width();
    const int height = world.height();

    outfile << "y,latitude";

    for (int x = 0; x <= width; x++)
        outfile << ",x" << x;

    outfile << '\n';
    outfile << fixed << setprecision(2);

    for (int y = 0; y <= height; y++)
    {
        outfile << y << ',' << latitudeforrow(y, height);

        for (int x = 0; x <= width; x++)
        {
            double value = 0.0;

            if (mode == 0)
                value = world.averainfloat(x, y);

            if (mode == 1)
                value = world.seasonalrainfloat(seasonjanuary, x, y);

            if (mode == 2)
                value = world.seasonalrainfloat(seasonjuly, x, y);

            if (mode == 3)
                value = isvalidationland(world, x, y) ? 1 : 0;

            outfile << ',' << value;
        }

        outfile << '\n';
    }
}

bool persistclimatebenchmarkdiagnostics(
    long seed,
    int runid,
    const string& timestamp,
    const string& information,
    const vector<string>& codes,
    const vector<long long>& simulationcounts)
{
    constexpr array<const char*, 16> diagnosticfiles = {
        "precipitation_summary.txt",
        "climate_energy_budget.csv",
        "climate_water_budget.csv",
        "climate_water_budget_area_weighted.csv",
        "climate_precipitation_distribution.csv",
        "climate_precipitation_processes.csv",
        "climate_hydrology_spinup.csv",
        "climate_condensation_activity.csv",
        "climate_atmosphere_budget.csv",
        "climate_circulation_precision.csv",
        "pressure_decomposition_comparison.csv",
        "temperature_threshold_comparison.csv",
        "monthly_climate_reference_comparison.csv",
        "monthly_physical_reference_comparison.csv",
        "ipcc_region_climate_comparison.csv",
        "annual_imerg_precipitation_comparison.csv"
    };
    const filesystem::path sourcedirectory = climatevalidationoutputdirectory(seed);
    const filesystem::path targetdirectory = sourcedirectory.parent_path() / "runs" / to_string(runid);
    error_code error;
    filesystem::create_directories(targetdirectory, error);

    if (error)
        return false;

    vector<string> copiedfiles;

    for (const char* filename : diagnosticfiles)
    {
        const filesystem::path source = sourcedirectory / filename;

        if (!filesystem::exists(source))
            continue;

        filesystem::copy_file(
            source,
            targetdirectory / filename,
            filesystem::copy_options::overwrite_existing,
            error);

        if (error)
            return false;

        copiedfiles.push_back(filename);
    }

    const filesystem::path rawpixelspath = targetdirectory / "climate_benchmark_raw_pixels.csv";
    ofstream rawpixels(rawpixelspath);

    if (!rawpixels.is_open() || codes.size() != simulationcounts.size())
        return false;

    rawpixels << "ID,DATETIME,INFORMATION";
    for (const string& code : codes)
        rawpixels << ',' << code;
    rawpixels << '\n' << runid << ',' << timestamp << ',' << csvescape(information);
    for (const long long count : simulationcounts)
        rawpixels << ',' << count;
    rawpixels << '\n';
    rawpixels.close();
    copiedfiles.push_back(rawpixelspath.filename().string());

    ofstream manifest(targetdirectory / "run_manifest.txt");

    if (!manifest.is_open())
        return false;

    manifest << "run_id=" << runid << '\n';
    manifest << "seed=" << seed << '\n';
    manifest << "datetime=" << timestamp << '\n';
    manifest << "information=" << information << '\n';
    manifest << "source_directory=" << sourcedirectory.generic_string() << '\n';
    manifest << "copied_files=";

    for (size_t index = 0; index < copiedfiles.size(); index++)
    {
        if (index > 0)
            manifest << ',';
        manifest << copiedfiles[index];
    }

    manifest << '\n';
    return true;
}
}

bool recordclimatebenchmarkrun(planet& world, const string& information, bool updateworkbook, int* runid)
{
    const int nextrunid = nextclimatebenchmarkrunid();
    const string timestamp = climatebenchmarktimestamp();
    const vector<string> codes = orderedclimatecodes();
    const vector<long long> simulationcounts = collectsimulatedclimatecounts(world);
    const vector<long long> referencecounts = referenceclimatecounts();
    const climatespatialmetrics spatial = compareclimatespatially(world);
    const benchmarkphysicsmetrics physics = collectbenchmarkphysicsmetrics(world);
    long long totalabsoluteerror = 0;
    long long totalreference = 0;

    for (size_t index = 0; index < simulationcounts.size() && index < referencecounts.size(); index++)
    {
        totalabsoluteerror += llabs(simulationcounts[index] - referencecounts[index]);
        totalreference += referencecounts[index];
    }

    const double weightedrelativeerror = totalreference > 0
        ? static_cast<double>(totalabsoluteerror) / static_cast<double>(totalreference)
        : 0.0;

    if (nextrunid < 2)
        return false;

    if (exportclimatebenchmarkimages(world, nextrunid) == false)
        return false;

    if (persistclimatebenchmarkdiagnostics(
        world.seed(), nextrunid, timestamp, information, codes, simulationcounts) == false)
    {
        cerr << "Failed to persist climate benchmark diagnostics.\n";
        return false;
    }

    if (appendclimatebenchmarkrunlog(
        nextrunid, timestamp, information, weightedrelativeerror, spatial, physics) == false)
    {
        cerr << "Failed to update climate benchmark run log.\n";
        return false;
    }

    if (updateworkbook)
    {
        string workbookfailure;

        if (updateclimatebenchmarkworkbook(
            nextrunid, codes, simulationcounts, &workbookfailure) == false)
        {
            const filesystem::path statuspath =
                climatevalidationoutputdirectory(world.seed()).parent_path() /
                "runs" / to_string(nextrunid) / "workbook_update_error.txt";
            ofstream statusfile(statuspath);

            if (statusfile.is_open())
                statusfile << workbookfailure << '\n';

            cerr
                << "Climate benchmark workbook update failed; run " << nextrunid
                << " remains recorded in JSON/CSV with all images and diagnostics.\n";
        }
    }

    if (runid != nullptr)
        *runid = nextrunid;

    return true;
}

void printclimaterelativeerrorreport(planet& world)
{
    const vector<string> codes = orderedclimatecodes();
    const vector<long long> references = referenceclimatecounts();
    const vector<long long> simulated = collectsimulatedclimatecounts(world);

    if (references.size() != codes.size() || simulated.size() != codes.size())
        return;

    double meanrelativeerror = 0.0;
    long long totalsimulated = 0;
    long long totalreference = 0;
    long long totalabsoluteerror = 0;
    double maxrelativeerror = -1.0;
    string maxcode;
    double adjustedmeanrelativeerror = 0.0;
    double adjustedmaxrelativeerror = -1.0;
    string adjustedmaxcode;
    int adjustedcount = 0;

    cout << "Climate relative error (raw pixels):" << '\n';
    cout << fixed << setprecision(4);

    for (size_t index = 0; index < codes.size(); index++)
    {
        const double relativeerror = saferelativeerror(simulated[index], references[index]);
        const long long absoluteerror = llabs(simulated[index] - references[index]);

        cout
            << codes[index]
            << " simulated=" << simulated[index]
            << " reference=" << references[index]
            << " relative_error=" << relativeerror
            << '\n';

        meanrelativeerror += relativeerror;
        totalsimulated += simulated[index];
        totalreference += references[index];
        totalabsoluteerror += absoluteerror;

        if (relativeerror > maxrelativeerror)
        {
            maxrelativeerror = relativeerror;
            maxcode = codes[index];
        }

        if (codes[index] == "Aw" || codes[index] == "As")
            continue;

        adjustedmeanrelativeerror += relativeerror;
        adjustedcount++;

        if (relativeerror > adjustedmaxrelativeerror)
        {
            adjustedmaxrelativeerror = relativeerror;
            adjustedmaxcode = codes[index];
        }
    }

    meanrelativeerror /= static_cast<double>(codes.size());

    const double awasrelativeerror = saferelativeerror(simulated[2] + simulated[3], references[2] + references[3]);
    adjustedmeanrelativeerror += awasrelativeerror;
    adjustedcount++;

    if (awasrelativeerror > adjustedmaxrelativeerror)
    {
        adjustedmaxrelativeerror = awasrelativeerror;
        adjustedmaxcode = "Aw+As";
    }

    adjustedmeanrelativeerror /= static_cast<double>(adjustedcount);

    const double weightedrelativeerror = totalreference > 0
        ? static_cast<double>(totalabsoluteerror) / static_cast<double>(totalreference)
        : 0.0;

    cout << "Aw+As combined"
         << " simulated=" << simulated[2] + simulated[3]
         << " reference=" << references[2] + references[3]
         << " relative_error=" << awasrelativeerror
         << '\n';
    cout << "Climate relative error summary:" << '\n';
    cout << "mean_relative_error=" << meanrelativeerror << '\n';
    cout << "mean_relative_error_adjusted=" << adjustedmeanrelativeerror << '\n';
    cout << "weighted_relative_error=" << weightedrelativeerror << '\n';
    cout << "max_relative_error=" << maxrelativeerror << " (" << maxcode << ")" << '\n';
    cout << "max_relative_error_adjusted=" << adjustedmaxrelativeerror << " (" << adjustedmaxcode << ")" << '\n';
    cout << "simulated_land_total=" << totalsimulated << '\n';
    cout << "reference_land_total=" << totalreference << '\n';

    const climatespatialmetrics spatial = compareclimatespatially(world);

    cout << "Climate spatial agreement:" << '\n';

    if (spatial.referencefound == false)
    {
        cout << "spatial_status=reference_not_found" << '\n';
        return;
    }

    if (spatial.dimensionsmatch == false)
    {
        cout << "spatial_status=dimension_mismatch" << '\n';
        return;
    }

    cout << "spatial_status=ok" << '\n';
    cout << "spatial_compared_cells=" << spatial.comparedcells << '\n';
    cout << "spatial_exact_accuracy=" << spatial.exactaccuracy << '\n';
    cout << "spatial_group_accuracy=" << spatial.groupaccuracy << '\n';
    cout << "spatial_kappa=" << spatial.kappa << '\n';
    cout << "area_weighted_spatial_exact_accuracy=" << spatial.areaweightedexactaccuracy << '\n';
    cout << "area_weighted_spatial_group_accuracy=" << spatial.areaweightedgroupaccuracy << '\n';
    cout << "area_weighted_spatial_kappa=" << spatial.areaweightedkappa << '\n';
    printreferencedfcdriverreport(world);
    printsimulatedclassconfusionreport(world, 6);
    printsimulatedclassconfusionreport(world, 1);
}

void exportclimatevalidationreport(planet& world)
{
    if (isworldgendebugrunactive() == false)
        return;

    const int width = world.width();
    const int height = world.height();
    const filesystem::path outputdir = climatevalidationoutputdirectory();

    ofstream energybudgetfile(outputdir / "climate_energy_budget.csv");

    if (energybudgetfile.is_open())
    {
        const climateenergy::AnnualEnergyBudget& budget = climateenergy::lastAnnualEnergyBudget();
        energybudgetfile << "incoming_solar_wm2,absorbed_solar_wm2,outgoing_longwave_wm2,atmospheric_transport_wm2,storage_tendency_wm2,residual_wm2,calibrated_longwave_intercept_wm2,area_weighted_mean_temperature_c,area_weighted_land_temperature_c,area_weighted_ocean_temperature_c,area_weighted_land_elevation_cooling_c,area_weighted_permanent_land_ice_fraction,local_ice_coupling_iterations,local_ice_coupling_maximum_row_residual,local_ice_coupling_temperature_error_c,local_ice_coupling_converged,northern_50_70_warmest_temperature_c,southern_60_90_warmest_temperature_c,southern_60_90_icecap_thermal_fraction\n";
        energybudgetfile << fixed << setprecision(9)
            << budget.incomingSolarWm2 << ','
            << budget.absorbedSolarWm2 << ','
            << budget.outgoingLongwaveWm2 << ','
            << budget.atmosphericTransportWm2 << ','
            << budget.storageTendencyWm2 << ','
            << budget.residualWm2 << ','
            << budget.calibratedLongwaveInterceptWm2 << ','
            << budget.areaWeightedMeanTemperatureC << ','
            << budget.areaWeightedLandTemperatureC << ','
            << budget.areaWeightedOceanTemperatureC << ','
            << budget.areaWeightedLandElevationCoolingC << ','
            << budget.areaWeightedPermanentLandIceFraction << ','
            << budget.localIceCouplingIterations << ','
            << budget.localIceCouplingMaximumRowResidual << ','
            << budget.localIceCouplingTemperatureErrorC << ','
            << (budget.localIceCouplingConverged ? 1 : 0) << ','
            << budget.northern5070WarmestTemperatureC << ','
            << budget.southern6090WarmestTemperatureC << ','
            << budget.southern6090IcecapThermalFraction << '\n';
    }

    ofstream waterbudgetfile(outputdir / "climate_water_budget.csv");

    if (waterbudgetfile.is_open())
    {
        waterbudgetfile << "quarter,initial_atmospheric_storage,initial_soil_storage,initial_snow_storage,ocean_evaporation,land_evaporation,ocean_precipitation,land_precipitation,runoff,atmospheric_storage,soil_storage,snow_storage,residual,relative_residual\n";
        waterbudgetfile << fixed << setprecision(9);

        const auto& budgets = climatephysics::lastWaterBudgets();

        for (int season = 0; season < CLIMATESEASONCOUNT; season++)
        {
            const climatephysics::WaterBudget& budget = budgets[season];
            waterbudgetfile
                << season << ','
                << budget.initialAtmosphericStorage << ','
                << budget.initialSoilStorage << ','
                << budget.initialSnowStorage << ','
                << budget.oceanEvaporation << ','
                << budget.landEvaporation << ','
                << budget.oceanPrecipitation << ','
                << budget.landPrecipitation << ','
                << budget.runoff << ','
                << budget.atmosphericStorage << ','
                << budget.soilStorage << ','
                << budget.snowStorage << ','
                << budget.residual() << ','
                << budget.relativeResidual() << '\n';
        }
    }

    ofstream areaweightedwaterbudgetfile(outputdir / "climate_water_budget_area_weighted.csv");

    if (areaweightedwaterbudgetfile.is_open())
    {
        areaweightedwaterbudgetfile << "quarter,initial_atmospheric_storage,initial_soil_storage,initial_snow_storage,ocean_evaporation,land_evaporation,ocean_precipitation,land_precipitation,runoff,atmospheric_storage,soil_storage,snow_storage,residual,relative_residual\n";
        areaweightedwaterbudgetfile << fixed << setprecision(9);

        const auto& budgets = climatephysics::lastAreaWeightedWaterBudgets();

        for (int season = 0; season < CLIMATESEASONCOUNT; season++)
        {
            const climatephysics::WaterBudget& budget = budgets[season];
            areaweightedwaterbudgetfile
                << season << ','
                << budget.initialAtmosphericStorage << ','
                << budget.initialSoilStorage << ','
                << budget.initialSnowStorage << ','
                << budget.oceanEvaporation << ','
                << budget.landEvaporation << ','
                << budget.oceanPrecipitation << ','
                << budget.landPrecipitation << ','
                << budget.runoff << ','
                << budget.atmosphericStorage << ','
                << budget.soilStorage << ','
                << budget.snowStorage << ','
                << budget.residual() << ','
                << budget.relativeResidual() << '\n';
        }
    }

    ofstream precipitationdistributionfile(outputdir / "climate_precipitation_distribution.csv");

    if (precipitationdistributionfile.is_open())
    {
        precipitationdistributionfile << "scope,cells,area_weight,raw_zero_fraction,stored_zero_fraction,below_1mm_fraction,area_weighted_raw_zero_fraction,area_weighted_stored_zero_fraction,area_weighted_below_1mm_fraction,mean_monthly_precipitation_mm,area_weighted_mean_monthly_precipitation_mm,maximum_monthly_precipitation_mm,maximum_x,maximum_y,maximum_month,rounded_precipitation_at_maximum_mm,wettest_10_percent_share\n";
        precipitationdistributionfile << fixed << setprecision(9);
        const auto& diagnostics = climatephysics::lastPrecipitationDistributionDiagnostics();

        auto writerow = [&](const char* scope, const climatephysics::PrecipitationDistributionScope& values)
        {
            precipitationdistributionfile
                << scope << ','
                << values.cells << ','
                << values.areaWeight << ','
                << values.rawZeroFraction << ','
                << values.storedZeroFraction << ','
                << values.belowOneMillimetreFraction << ','
                << values.areaWeightedRawZeroFraction << ','
                << values.areaWeightedStoredZeroFraction << ','
                << values.areaWeightedBelowOneMillimetreFraction << ','
                << values.meanMonthlyPrecipitationMm << ','
                << values.areaWeightedMeanMonthlyPrecipitationMm << ','
                << values.maximumMonthlyPrecipitationMm << ','
                << values.maximumMonthlyPrecipitationX << ','
                << values.maximumMonthlyPrecipitationY << ','
                << values.maximumMonthlyPrecipitationMonth << ','
                << values.roundedPrecipitationAtMaximumMm << ','
                << values.wettestTenPercentShare << '\n';
        };

        writerow("land", diagnostics.land);
        writerow("ocean", diagnostics.ocean);
    }

    ofstream spinupfile(outputdir / "climate_hydrology_spinup.csv");

    if (spinupfile.is_open())
    {
        const climatephysics::HydrologySpinupDiagnostics& diagnostics =
            climatephysics::lastHydrologySpinupDiagnostics();
        spinupfile << "cycles_completed,converged,relative_storage_change,relative_atmospheric_storage_change,relative_soil_storage_change,relative_snow_storage_change,relative_snow_cover_change,atmospheric_storage,soil_storage,snow_storage\n";
        spinupfile << fixed << setprecision(9)
            << diagnostics.cyclesCompleted << ','
            << (diagnostics.converged ? 1 : 0) << ','
            << diagnostics.relativeStorageChange << ','
            << diagnostics.relativeAtmosphericStorageChange << ','
            << diagnostics.relativeSoilStorageChange << ','
            << diagnostics.relativeSnowStorageChange << ','
            << diagnostics.relativeSnowCoverChange << ','
            << diagnostics.atmosphericStorage << ','
            << diagnostics.soilStorage << ','
            << diagnostics.snowStorage << '\n';
    }

    ofstream condensationfile(outputdir / "climate_condensation_activity.csv");

    if (condensationfile.is_open())
    {
        constexpr array<const char*, CLIMATESEASONCOUNT> seasonnames = {
            "january-march", "april-june", "july-september", "october-december"
        };
        const auto& activity = climatephysics::lastCondensationActivityDiagnostics();
        const auto& budgets = climatephysics::lastAreaWeightedWaterBudgets();
        constexpr double representativeperioddays = 365.0 / CLIMATESEASONCOUNT;
        condensationfile << "season,cell_step_area_weight,active_cell_step_fraction,active_atmospheric_water_fraction,condensable_excess_water_fraction,precipitation_to_excess_fraction,atmospheric_water_residence_time_days\n";
        condensationfile << fixed << setprecision(9);

        for (int season = 0; season < CLIMATESEASONCOUNT; season++)
        {
            const double precipitation = budgets[season].oceanPrecipitation +
                budgets[season].landPrecipitation;
            condensationfile
                << seasonnames[season] << ','
                << activity[season].cellStepAreaWeight << ','
                << (activity[season].cellStepAreaWeight > 0.0 ?
                    activity[season].activeCellStepAreaWeight /
                        activity[season].cellStepAreaWeight : 0.0) << ','
                << (activity[season].atmosphericWaterAreaWeighted > 0.0 ?
                    activity[season].activeAtmosphericWaterAreaWeighted /
                        activity[season].atmosphericWaterAreaWeighted : 0.0) << ','
                << (activity[season].atmosphericWaterAreaWeighted > 0.0 ?
                    activity[season].excessWaterAreaWeighted /
                        activity[season].atmosphericWaterAreaWeighted : 0.0) << ','
                << (activity[season].excessWaterAreaWeighted > 0.0 ?
                    activity[season].precipitationAreaWeighted /
                        activity[season].excessWaterAreaWeighted : 0.0) << ','
                << (precipitation > 0.0 ?
                    representativeperioddays * budgets[season].atmosphericStorage /
                        precipitation : 0.0) << '\n';
        }
    }

    ofstream precipitationprocessfile(outputdir / "climate_precipitation_processes.csv");

    if (precipitationprocessfile.is_open())
    {
        constexpr array<const char*, CLIMATESEASONCOUNT> quarternames = {
            "january-march", "april-june", "july-september", "october-december"
        };
        const auto& processes = climatephysics::lastPrecipitationProcessDiagnostics();
        const auto& budgets = climatephysics::lastAreaWeightedWaterBudgets();
        precipitationprocessfile
            << "quarter,stratiform_precipitation_mm_area_weighted,orographic_precipitation_mm_area_weighted,convective_precipitation_mm_area_weighted,total_precipitation_mm_area_weighted,water_budget_precipitation_mm_area_weighted,component_budget_difference_mm_area_weighted,reevaporated_precipitation_mm_area_weighted,snowfall_mm_area_weighted,upward_moisture_transfer_mm_area_weighted,downward_moisture_transfer_mm_area_weighted,cloud_fraction_area_weighted,positive_moisture_flux_convergence_mm_area_weighted,negative_moisture_flux_convergence_mm_area_weighted,net_moisture_flux_convergence_mm_area_weighted\n";
        precipitationprocessfile << fixed << setprecision(9);

        for (int season = 0; season < CLIMATESEASONCOUNT; season++)
        {
            const auto& process = processes[season];
            const double processtotal = process.stratiformPrecipitation +
                process.orographicPrecipitation + process.convectivePrecipitation;
            const double budgettotal = budgets[season].oceanPrecipitation +
                budgets[season].landPrecipitation;
            precipitationprocessfile
                << quarternames[season] << ','
                << process.stratiformPrecipitation << ','
                << process.orographicPrecipitation << ','
                << process.convectivePrecipitation << ','
                << processtotal << ','
                << budgettotal << ','
                << processtotal - budgettotal << ','
                << process.reevaporatedPrecipitation << ','
                << process.snowfall << ','
                << process.upwardMoistureTransfer << ','
                << process.downwardMoistureTransfer << ','
                << process.cloudFractionAreaWeighted << ','
                << process.positiveMoistureFluxConvergence << ','
                << process.negativeMoistureFluxConvergence << ','
                << process.positiveMoistureFluxConvergence +
                    process.negativeMoistureFluxConvergence << '\n';
        }
    }

    writemonthlyreferencecomparison(outputdir, world);
    writeimergreferencecomparison(outputdir, world);
    writephysicalreferencecomparison(outputdir, world);
    writeipccregioncomparison(outputdir, world);
    writeatmosphericbudget(outputdir, world);
    writecirculationprecision(outputdir);
    writepressuredecompositioncomparison(outputdir, world);
    writetemperaturethresholdcomparison(outputdir, world);

    vector<zonalstats> rows(height + 1);
    map<string, int> climatecounts;

    double globalannualrain = 0.0;
    double landannualrain = 0.0;
    double oceanannualrain = 0.0;
    double northernannualrain = 0.0;
    double southernannualrain = 0.0;
    int globalcells = 0;
    int landcells = 0;
    int oceancells = 0;
    int northerncells = 0;
    int southerncells = 0;

    for (int y = 0; y <= height; y++)
    {
        zonalstats& row = rows[y];

        for (int x = 0; x <= width; x++)
        {
            const bool sea = world.sea(x, y) == 1;
            const double annualrain = world.averainfloat(x, y);
            const double januaryrain =
                world.seasonalrainfloat(seasonjanuary, x, y);
            const double julyrain =
                world.seasonalrainfloat(seasonjuly, x, y);

            row.cells++;
            row.annualrain = row.annualrain + annualrain;
            row.januaryrain = row.januaryrain + januaryrain;
            row.julyrain = row.julyrain + julyrain;
            row.januarypressure = row.januarypressure + world.seasonalpressure(seasonjanuary, x, y);
            row.julypressure = row.julypressure + world.seasonalpressure(seasonjuly, x, y);
            row.januaryuwind = row.januaryuwind + world.seasonaluwind(seasonjanuary, x, y);
            row.januaryvwind = row.januaryvwind + world.seasonalvwind(seasonjanuary, x, y);
            row.julyuwind = row.julyuwind + world.seasonaluwind(seasonjuly, x, y);
            row.julyvwind = row.julyvwind + world.seasonalvwind(seasonjuly, x, y);
            row.januaryevaporation = row.januaryevaporation + world.seasonalevaporation(seasonjanuary, x, y);
            row.julyevaporation = row.julyevaporation + world.seasonalevaporation(seasonjuly, x, y);
            row.januarymoisture = row.januarymoisture + world.seasonalmoisture(seasonjanuary, x, y);
            row.julymoisture = row.julymoisture + world.seasonalmoisture(seasonjuly, x, y);

            globalcells++;
            globalannualrain = globalannualrain + annualrain;

            if (y < height / 2)
            {
                northernannualrain = northernannualrain + annualrain;
                northerncells++;
            }
            else if (y > height / 2)
            {
                southernannualrain = southernannualrain + annualrain;
                southerncells++;
            }

            if (sea)
            {
                row.oceancells++;
                row.oceanannualrain = row.oceanannualrain + annualrain;
                row.januarysst = row.januarysst + world.seasonalsst(seasonjanuary, x, y);
                row.julysst = row.julysst + world.seasonalsst(seasonjuly, x, y);
                row.januarycurrentu = row.januarycurrentu + world.seasonalcurrentu(seasonjanuary, x, y);
                row.januarycurrentv = row.januarycurrentv + world.seasonalcurrentv(seasonjanuary, x, y);
                row.julycurrentu = row.julycurrentu + world.seasonalcurrentu(seasonjuly, x, y);
                row.julycurrentv = row.julycurrentv + world.seasonalcurrentv(seasonjuly, x, y);

                oceancells++;
                oceanannualrain = oceanannualrain + annualrain;
            }
            else
            {
                row.landcells++;
                row.landannualrain = row.landannualrain + annualrain;

                landcells++;
                landannualrain = landannualrain + annualrain;

                const short climate = static_cast<short>(world.climate(x, y));

                if (climate > 0)
                    climatecounts[getclimatecode(climate)]++;
            }
        }
    }

    ofstream zonalfile(outputdir / "precipitation_zonal.csv");

    if (zonalfile.is_open())
    {
        zonalfile << "y,latitude,cells,land_cells,ocean_cells,mean_annual_rain,mean_jan_rain,mean_jul_rain,land_mean_annual_rain,ocean_mean_annual_rain,mean_jan_pressure,mean_jul_pressure,mean_jan_u_wind,mean_jan_v_wind,mean_jul_u_wind,mean_jul_v_wind,ocean_mean_jan_sst,ocean_mean_jul_sst,ocean_mean_jan_current_u,ocean_mean_jan_current_v,ocean_mean_jul_current_u,ocean_mean_jul_current_v,mean_jan_evaporation,mean_jul_evaporation,mean_jan_moisture,mean_jul_moisture\n";
        zonalfile << fixed << setprecision(4);

        for (int y = 0; y <= height; y++)
        {
            const zonalstats& row = rows[y];

            zonalfile
                << y << ','
                << latitudeforrow(y, height) << ','
                << row.cells << ','
                << row.landcells << ','
                << row.oceancells << ','
                << safeaverage(row.annualrain, row.cells) << ','
                << safeaverage(row.januaryrain, row.cells) << ','
                << safeaverage(row.julyrain, row.cells) << ','
                << safeaverage(row.landannualrain, row.landcells) << ','
                << safeaverage(row.oceanannualrain, row.oceancells) << ','
                << safeaverage(row.januarypressure, row.cells) << ','
                << safeaverage(row.julypressure, row.cells) << ','
                << safeaverage(row.januaryuwind, row.cells) << ','
                << safeaverage(row.januaryvwind, row.cells) << ','
                << safeaverage(row.julyuwind, row.cells) << ','
                << safeaverage(row.julyvwind, row.cells) << ','
                << safeaverage(row.januarysst, row.oceancells) << ','
                << safeaverage(row.julysst, row.oceancells) << ','
                << safeaverage(row.januarycurrentu, row.oceancells) << ','
                << safeaverage(row.januarycurrentv, row.oceancells) << ','
                << safeaverage(row.julycurrentu, row.oceancells) << ','
                << safeaverage(row.julycurrentv, row.oceancells) << ','
                << safeaverage(row.januaryevaporation, row.cells) << ','
                << safeaverage(row.julyevaporation, row.cells) << ','
                << safeaverage(row.januarymoisture, row.cells) << ','
                << safeaverage(row.julymoisture, row.cells) << '\n';
        }
    }

    int annualitczrow = height / 2;
    int januaryitczrow = height / 2;
    int julyitczrow = height / 2;
    int northdryrow = height / 4;
    int southdryrow = (height * 3) / 4;
    double annualitczrain = -1.0;
    double januaryitczrain = -1.0;
    double julyitczrain = -1.0;
    double northdryrain = 1e30;
    double southdryrain = 1e30;

    for (int y = 0; y <= height; y++)
    {
        const float latitude = latitudeforrow(y, height);
        const zonalstats& row = rows[y];
        const double meanannualrain = safeaverage(row.annualrain, row.cells);
        const double meanjanuaryrain = safeaverage(row.januaryrain, row.cells);
        const double meanjulyrain = safeaverage(row.julyrain, row.cells);

        if (fabs(latitude) <= 30.0f)
        {
            if (meanannualrain > annualitczrain)
            {
                annualitczrain = meanannualrain;
                annualitczrow = y;
            }

            if (meanjanuaryrain > januaryitczrain)
            {
                januaryitczrain = meanjanuaryrain;
                januaryitczrow = y;
            }

            if (meanjulyrain > julyitczrain)
            {
                julyitczrain = meanjulyrain;
                julyitczrow = y;
            }
        }

        if (latitude >= 10.0f && latitude <= 45.0f && meanannualrain < northdryrain)
        {
            northdryrain = meanannualrain;
            northdryrow = y;
        }

        if (latitude <= -10.0f && latitude >= -45.0f && meanannualrain < southdryrain)
        {
            southdryrain = meanannualrain;
            southdryrow = y;
        }
    }

    const comparisonmetrics comparison = compareannualprecipitation(outputdir, world, rows);
    ofstream summaryfile(outputdir / "precipitation_summary.txt");

    if (summaryfile.is_open())
    {
        summaryfile << fixed << setprecision(4);
        summaryfile << "seed=" << worldgenerationdebugseed() << '\n';
        summaryfile << "width=" << width << '\n';
        summaryfile << "height=" << height << '\n';
        summaryfile << "global_mean_annual_precip=" << safeaverage(globalannualrain, globalcells) << '\n';
        summaryfile << "land_mean_annual_precip=" << safeaverage(landannualrain, landcells) << '\n';
        summaryfile << "ocean_mean_annual_precip=" << safeaverage(oceanannualrain, oceancells) << '\n';
        summaryfile << "northern_mean_annual_precip=" << safeaverage(northernannualrain, northerncells) << '\n';
        summaryfile << "southern_mean_annual_precip=" << safeaverage(southernannualrain, southerncells) << '\n';
        summaryfile << "annual_itcz_latitude=" << latitudeforrow(annualitczrow, height) << '\n';
        summaryfile << "january_itcz_latitude=" << latitudeforrow(januaryitczrow, height) << '\n';
        summaryfile << "july_itcz_latitude=" << latitudeforrow(julyitczrow, height) << '\n';
        summaryfile << "north_subtropical_dry_latitude=" << latitudeforrow(northdryrow, height) << '\n';
        summaryfile << "south_subtropical_dry_latitude=" << latitudeforrow(southdryrow, height) << '\n';
        summaryfile << "annual_itcz_zonal_precip=" << annualitczrain << '\n';
        summaryfile << "north_subtropical_dry_zonal_precip=" << northdryrain << '\n';
        summaryfile << "south_subtropical_dry_zonal_precip=" << southdryrain << '\n';
        summaryfile << "grid_files=annual_precipitation_grid.csv,january_precipitation_grid.csv,july_precipitation_grid.csv,land_mask_grid.csv\n";
        summaryfile << "physics_diagnostics=climate_energy_budget.csv,climate_water_budget.csv,climate_water_budget_area_weighted.csv,climate_precipitation_distribution.csv,climate_hydrology_spinup.csv,climate_condensation_activity.csv,climate_precipitation_processes.csv,climate_atmosphere_budget.csv,climate_circulation_precision.csv,pressure_decomposition_comparison.csv,temperature_threshold_comparison.csv,monthly_climate_reference_comparison.csv,monthly_physical_reference_comparison.csv,ipcc_region_climate_comparison.csv,annual_imerg_precipitation_comparison.csv\n";
        summaryfile << "reference_grid_path=" << getappenvironment().referencePrecipitationGridPath.string() << '\n';
        summaryfile << "reference_found=" << (comparison.referencefound ? 1 : 0) << '\n';
        summaryfile << "reference_dimensions_match=" << (comparison.dimensionsmatch ? 1 : 0) << '\n';

        if (comparison.dimensionsmatch)
        {
            summaryfile << "reference_compared_cells=" << comparison.comparedcells << '\n';
            summaryfile << "reference_mean_bias=" << comparison.meanbias << '\n';
            summaryfile << "reference_mae=" << comparison.meanabsoluteerror << '\n';
            summaryfile << "reference_rmse=" << comparison.rmse << '\n';
            summaryfile << "reference_correlation=" << comparison.correlation << '\n';
            summaryfile << "reference_tropical_mean_bias=" << comparison.tropicalmeanbias << '\n';
        }
    }

    ofstream climatefile(outputdir / "climate_counts.csv");

    if (climatefile.is_open())
    {
        climatefile << "code,name,cells\n";

        for (const auto& entry : climatecounts)
        {
            const short climate = climatefromcode(entry.first);

            climatefile
                << entry.first << ','
                << csvescape(getclimatename(climate)) << ','
                << entry.second << '\n';
        }
    }

    writeprecipitationgrid(outputdir / "annual_precipitation_grid.csv", world, 0);
    writeprecipitationgrid(outputdir / "january_precipitation_grid.csv", world, 1);
    writeprecipitationgrid(outputdir / "july_precipitation_grid.csv", world, 2);
    writeprecipitationgrid(outputdir / "land_mask_grid.csv", world, 3);
}
