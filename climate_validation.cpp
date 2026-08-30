#include "app_environment.hpp"
#include "climate_atmosphere.hpp"
#include "climate_energy.hpp"
#include "climate_physics.hpp"
#include "climate_reference.hpp"
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
            if (world.sea(x, y) == 1)
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
    surfacewindspeed
};

enum class referencescope
{
    global,
    land,
    ocean
};

double simulatedmonthlyvalue(planet& world, monthlyreferencefield field, int season, int x, int y)
{
    if (field == monthlyreferencefield::temperature)
        return world.seasonaltemp(season, x, y);

    if (field == monthlyreferencefield::precipitation)
        return world.seasonalrain(season, x, y);

    const double u = world.seasonaluwind(season, x, y);
    const double v = world.seasonalvwind(season, x, y);
    return std::sqrt(u * u + v * v);
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
            if (scope == referencescope::land && world.sea(x, y) == 1)
                continue;

            if (scope == referencescope::ocean && world.sea(x, y) == 0)
                continue;

            double simulated = 0.0;
            double observed = 0.0;

            if (season >= 0)
            {
                simulated = simulatedmonthlyvalue(world, field, season, x, y);
                observed = reference.value(referenceMonths[season], x, y);
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
    constexpr array<int, CLIMATESEASONCOUNT> referenceMonths = { 1, 4, 7, 10 };

    auto writecomparisons = [&](const char* variable, const climatereference::MonthlyGrid& grid, monthlyreferencefield field)
    {
        for (int season = -1; season < CLIMATESEASONCOUNT; season++)
        {
            const comparisonmetrics metrics = comparemonthlyfield(world, grid, season, field);
            output
                << variable << ','
                << (season < 0 ? "annual_mean" : periodNames[season]) << ','
                << (season < 0 ? 0 : referenceMonths[season]) << ','
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

filesystem::path climatevalidationoutputdirectory()
{
    filesystem::path outputroot = getappenvironment().profilingWorkbookPath.parent_path();

    if (outputroot.empty())
        outputroot = filesystem::current_path();

    outputroot /= "validation";
    outputroot /= "seed_" + to_string(worldgenerationdebugseed());
    filesystem::create_directories(outputroot);

    return outputroot;
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
            if (world.sea(x, y) == 1)
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
            if (world.sea(x, y) == 1)
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
                rains[season] = static_cast<float>(world.seasonalrain(season, x, y));
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

            if (driestcoldrain < wettestwarmrain / tuning::climate::koppen::continentalWinterDrynessDivisor)
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

bool appendclimatebenchmarkrunlog(
    int runid,
    const string& timestamp,
    const string& information,
    double weightedrelativeerror,
    const climatespatialmetrics& spatial)
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
        << "      \"metrics\": {\n"
        << fixed << setprecision(10)
        << "        \"weighted_relative_error\": " << weightedrelativeerror << ",\n"
        << "        \"spatial_compared_cells\": " << spatial.comparedcells << ",\n"
        << "        \"spatial_exact_accuracy\": " << spatial.exactaccuracy << ",\n"
        << "        \"spatial_group_accuracy\": " << spatial.groupaccuracy << ",\n"
        << "        \"spatial_kappa\": " << spatial.kappa << ",\n"
        << "        \"area_weighted_spatial_exact_accuracy\": " << spatial.areaweightedexactaccuracy << ",\n"
        << "        \"area_weighted_spatial_group_accuracy\": " << spatial.areaweightedgroupaccuracy << ",\n"
        << "        \"area_weighted_spatial_kappa\": " << spatial.areaweightedkappa << "\n"
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

bool updateclimatebenchmarkworkbook(int runid, const vector<string>& codes, const vector<long long>& simulationcounts)
{
    const filesystem::path workbookpath = filesystem::absolute(getappenvironment().climateWorkbookPath).lexically_normal();

    if (filesystem::exists(workbookpath) == false || simulationcounts.size() != codes.size())
        return false;

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
            cerr << "Climate benchmark workbook error: " << errormessage << '\n';
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
            if (world.sea(x, y) == 1)
                continue;

            const double simulated = static_cast<double>(world.averain(x, y));
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

    outfile << "season,area_weighted_mean_pressure_anomaly_hpa,rms_pressure_anomaly_hpa,pressure_limiter_fraction,area_weighted_mean_surface_divergence_per_day,area_weighted_mean_upper_divergence_per_day,area_weighted_mean_vertical_velocity_hpa_per_day,rms_surface_wind_m_s,surface_wind_limiter_fraction,rms_upper_wind_m_s,rms_vertical_velocity_hpa_per_day,vertical_limiter_fraction,rms_column_divergence_per_day\n";
    outfile << fixed << setprecision(9);

    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
    {
        double weighttotal = 0.0;
        double pressuretotal = 0.0;
        double pressuresquared = 0.0;
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

        for (int y = 0; y <= height; y++)
        {
            const int ynorth = y > 0 ? y - 1 : y;
            const int ysouth = y < height ? y + 1 : y;
            const double weight = gridcellareaweight(y, height);
            const auto spacing = climateatmosphere::cellSpacingMetres(
                latitudeforrow(y, height),
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
                        (world.seasonalvwind(season, x, ysouth) - world.seasonalvwind(season, x, ynorth)) /
                        (2.0 * spacing.meridionalMetres)) * tuning::climate::circulation::secondsPerDay;
                const double upperdivergence =
                    ((world.seasonalupperuwind(season, xeast, y) - world.seasonalupperuwind(season, xwest, y)) /
                        (2.0 * spacing.zonalMetres) +
                        (world.seasonaluppervwind(season, x, ysouth) - world.seasonaluppervwind(season, x, ynorth)) /
                        (2.0 * spacing.meridionalMetres)) * tuning::climate::circulation::secondsPerDay;
                const double verticalvelocity =
                    static_cast<double>(world.seasonalverticalvelocity(season, x, y)) / verticalstoragescale;
                const double columndivergence = surfacedivergence + upperdivergence;

                weighttotal += weight;
                pressuretotal += weight * world.seasonalpressure(season, x, y);
                pressuresquared += weight * world.seasonalpressure(season, x, y) *
                    world.seasonalpressure(season, x, y);
                surfacedivergencetotal += weight * surfacedivergence;
                upperdivergencetotal += weight * upperdivergence;
                verticaltotal += weight * verticalvelocity;
                surfacewindsquared += weight * (surfaceu * surfaceu + surfacev * surfacev);
                upperwindsquared += weight * (upperu * upperu + upperv * upperv);
                verticalsquared += weight * verticalvelocity * verticalvelocity;
                columndivergencesquared += weight * columndivergence * columndivergence;

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
            << pressurelimitedweight * inverseweight << ','
            << surfacedivergencetotal * inverseweight << ','
            << upperdivergencetotal * inverseweight << ','
            << verticaltotal * inverseweight << ','
            << sqrt(surfacewindsquared * inverseweight) << ','
            << surfacewindlimitedweight * inverseweight << ','
            << sqrt(upperwindsquared * inverseweight) << ','
            << sqrt(verticalsquared * inverseweight) << ','
            << verticallimitedweight * inverseweight << ','
            << sqrt(columndivergencesquared * inverseweight) << '\n';
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
            int value = 0;

            if (mode == 0)
                value = world.averain(x, y);

            if (mode == 1)
                value = world.janrain(x, y);

            if (mode == 2)
                value = world.julrain(x, y);

            if (mode == 3)
                value = world.sea(x, y) == 0 ? 1 : 0;

            outfile << ',' << value;
        }

        outfile << '\n';
    }
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

    if (updateworkbook && updateclimatebenchmarkworkbook(nextrunid, codes, simulationcounts) == false)
    {
        cerr << "Failed to update climate benchmark workbook.\n";
        return false;
    }

    if (exportclimatebenchmarkimages(world, nextrunid) == false)
        return false;

    if (appendclimatebenchmarkrunlog(nextrunid, timestamp, information, weightedrelativeerror, spatial) == false)
    {
        cerr << "Failed to update climate benchmark run log.\n";
        return false;
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
        energybudgetfile << "incoming_solar_wm2,absorbed_solar_wm2,outgoing_longwave_wm2,atmospheric_transport_wm2,storage_tendency_wm2,residual_wm2,calibrated_longwave_intercept_wm2,area_weighted_mean_temperature_c,area_weighted_land_temperature_c,area_weighted_ocean_temperature_c,area_weighted_land_elevation_cooling_c\n";
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
            << budget.areaWeightedLandElevationCoolingC << '\n';
    }

    ofstream waterbudgetfile(outputdir / "climate_water_budget.csv");

    if (waterbudgetfile.is_open())
    {
        waterbudgetfile << "season,initial_atmospheric_storage,initial_soil_storage,ocean_evaporation,land_evaporation,ocean_precipitation,land_precipitation,runoff,atmospheric_storage,soil_storage,residual,relative_residual\n";
        waterbudgetfile << fixed << setprecision(9);

        const auto& budgets = climatephysics::lastWaterBudgets();

        for (int season = 0; season < CLIMATESEASONCOUNT; season++)
        {
            const climatephysics::WaterBudget& budget = budgets[season];
            waterbudgetfile
                << season << ','
                << budget.initialAtmosphericStorage << ','
                << budget.initialSoilStorage << ','
                << budget.oceanEvaporation << ','
                << budget.landEvaporation << ','
                << budget.oceanPrecipitation << ','
                << budget.landPrecipitation << ','
                << budget.runoff << ','
                << budget.atmosphericStorage << ','
                << budget.soilStorage << ','
                << budget.residual() << ','
                << budget.relativeResidual() << '\n';
        }
    }

    ofstream spinupfile(outputdir / "climate_hydrology_spinup.csv");

    if (spinupfile.is_open())
    {
        const climatephysics::HydrologySpinupDiagnostics& diagnostics =
            climatephysics::lastHydrologySpinupDiagnostics();
        spinupfile << "cycles_completed,converged,relative_storage_change,atmospheric_storage,soil_storage\n";
        spinupfile << fixed << setprecision(9)
            << diagnostics.cyclesCompleted << ','
            << (diagnostics.converged ? 1 : 0) << ','
            << diagnostics.relativeStorageChange << ','
            << diagnostics.atmosphericStorage << ','
            << diagnostics.soilStorage << '\n';
    }

    writemonthlyreferencecomparison(outputdir, world);
    writeimergreferencecomparison(outputdir, world);
    writeatmosphericbudget(outputdir, world);

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
            const int annualrain = world.averain(x, y);
            const int januaryrain = world.janrain(x, y);
            const int julyrain = world.julrain(x, y);

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
        summaryfile << "physics_diagnostics=climate_energy_budget.csv,climate_water_budget.csv,climate_hydrology_spinup.csv,climate_atmosphere_budget.csv,monthly_climate_reference_comparison.csv,annual_imerg_precipitation_comparison.csv\n";
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
