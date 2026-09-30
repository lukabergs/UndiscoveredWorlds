#include <algorithm>
#include <array>
#include <cmath>
#include "planet.hpp"
#include "region.hpp"
#include "climate_koppen.hpp"
#include "climate_classification.hpp"

using namespace std;

namespace
{
struct seasonalclimatesummary
{
    std::array<float, CLIMATESEASONCOUNT> temps{};
    std::array<float, CLIMATESEASONCOUNT> rains{};
    float mintemp = 0.0f;
    float maxtemp = 0.0f;
    float meanannualtemp = 0.0f;
    float annualrain = 0.0f;
    float warmhalffraction = 0.0f;
    float minrain = 0.0f;
    float driestwarmrain = 0.0f;
    float wettestwarmrain = 0.0f;
    float driestcoldrain = 0.0f;
    float wettestcoldrain = 0.0f;
    float secondwarmesttemp = 0.0f;
    bool driestseasoniswarm = false;
};

bool isclassificationwater(planet& world, int x, int y);
seasonalclimatesummary summariseseasonalclimate(planet& world, int x, int y);
short calculateplanetclimate(planet& world, int x, int y);
int calculateplanetbiome(planet& world, int x, int y);
float clampholdridgebiotemperature(float temperature);
float calculateholdridgebiotemperature(const seasonalclimatesummary& summary);
int calculateholdridgebiomeinternal(float janTemp, float aprTemp, float julTemp, float octTemp, float janRain, float aprRain, float julRain, float octRain);
int calculateholdridgebiomefromseasonal(const seasonalclimatesummary& summary);
short calculateclimatefromseasonal(int elev, int sealevel, const seasonalclimatesummary& summary);

bool isclassificationwater(planet& world, int x, int y)
{
    return world.sea(x, y) == 1 || world.truelake(x, y) != 0 || world.riftlakesurface(x, y) != 0;
}

seasonalclimatesummary summariseseasonalclimate(planet& world, int x, int y)
{
    seasonalclimatesummary summary;
    std::array<int, CLIMATESEASONCOUNT> order = { 0, 1, 2, 3 };

    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
    {
        summary.temps[season] = static_cast<float>(world.seasonaltemp(season, x, y));
        summary.rains[season] = world.seasonalrainfloat(season, x, y);
        summary.meanannualtemp += summary.temps[season];
        summary.annualrain += summary.rains[season];
    }

    summary.meanannualtemp = summary.meanannualtemp / static_cast<float>(CLIMATESEASONCOUNT);
    summary.annualrain = summary.annualrain * 3.0f;
    summary.mintemp = summary.temps[0];
    summary.maxtemp = summary.temps[0];
    summary.minrain = summary.rains[0];
    int driestseason = 0;

    for (int season = 1; season < CLIMATESEASONCOUNT; season++)
    {
        if (summary.temps[season] < summary.mintemp)
            summary.mintemp = summary.temps[season];

        if (summary.temps[season] > summary.maxtemp)
            summary.maxtemp = summary.temps[season];

        if (summary.rains[season] < summary.minrain)
        {
            summary.minrain = summary.rains[season];
            driestseason = season;
        }
    }

    for (int left = 0; left < CLIMATESEASONCOUNT - 1; left++)
    {
        for (int right = left + 1; right < CLIMATESEASONCOUNT; right++)
        {
            if (summary.temps[order[left]] > summary.temps[order[right]])
                std::swap(order[left], order[right]);
        }
    }

    const int cold0 = order[0];
    const int cold1 = order[1];
    const int warm0 = order[2];
    const int warm1 = order[3];
    const float coldrain0 = summary.rains[cold0];
    const float coldrain1 = summary.rains[cold1];
    const float warmrain0 = summary.rains[warm0];
    const float warmrain1 = summary.rains[warm1];
    const float totalseasonalrain = summary.rains[0] + summary.rains[1] + summary.rains[2] + summary.rains[3];

    summary.secondwarmesttemp = summary.temps[warm0];
    summary.driestcoldrain = std::min(coldrain0, coldrain1);
    summary.wettestcoldrain = std::max(coldrain0, coldrain1);
    summary.driestwarmrain = std::min(warmrain0, warmrain1);
    summary.wettestwarmrain = std::max(warmrain0, warmrain1);
    summary.driestseasoniswarm = (driestseason == warm0 || driestseason == warm1);

    if (totalseasonalrain > 0.0f)
        summary.warmhalffraction = (warmrain0 + warmrain1) / totalseasonalrain;

    return summary;
}

short calculateplanetclimate(planet& world, int x, int y)
{
    const auto& monthly = world.climatesimulation().monthly;
    if (monthly.columns == world.width() + 1 && monthly.rows == world.height() + 1 &&
        monthly.hasTemperature() && monthly.hasRainfall())
    {
        std::array<float, 12> temperature, rain;
        const int cell = y * monthly.columns + x;
        for (int month = 0; month < 12; ++month)
        {
            const auto calendar = climatehydrology::calendarMonth(month);
            const auto temperatureAdjustment = [&](int season)
            { return world.seasonaltemp(season, x, y) - monthly.referenceTemperatureC[season][cell]; };
            temperature[month] = monthly.temperatureC[month][cell] + climatehydrology::interpolateSeasonal(
                temperatureAdjustment(calendar.firstSeason), temperatureAdjustment(calendar.secondSeason), calendar.interpolation);
            const int season = month / 3;
            const float reference = monthly.referenceRainfallMm[season][cell];
            rain[month] = reference > 0.0f ? monthly.rainfallMm[month][cell] *
                world.seasonalrainfloat(season, x, y) / reference : world.seasonalrainfloat(season, x, y);
        }
        return climatekoppen::classifyMonthly(temperature, rain, y < monthly.rows / 2);
    }
    const int elev = world.map(x, y);
    const int sealevel = world.sealevel();
    const seasonalclimatesummary summary = summariseseasonalclimate(world, x, y);
    return calculateclimatefromseasonal(elev, sealevel, summary);
}

int calculateplanetbiome(planet& world, int x, int y)
{
    const seasonalclimatesummary summary = summariseseasonalclimate(world, x, y);
    return calculateholdridgebiomefromseasonal(summary);
}

float clampholdridgebiotemperature(float temperature)
{
    return clamp(temperature, 0.0f, 30.0f);
}

float calculateholdridgebiotemperature(const seasonalclimatesummary& summary)
{
    float total = 0.0f;

    for (float temp : summary.temps)
        total = total + clampholdridgebiotemperature(temp);

    return total / static_cast<float>(CLIMATESEASONCOUNT);
}

int calculateholdridgebiomeinternal(float janTemp, float aprTemp, float julTemp, float octTemp, float janRain, float aprRain, float julRain, float octRain)
{
    seasonalclimatesummary summary;
    summary.temps[seasonjanuary] = janTemp;
    summary.temps[seasonapril] = aprTemp;
    summary.temps[seasonjuly] = julTemp;
    summary.temps[seasonoctober] = octTemp;
    summary.rains[seasonjanuary] = janRain;
    summary.rains[seasonapril] = aprRain;
    summary.rains[seasonjuly] = julRain;
    summary.rains[seasonoctober] = octRain;
    summary.annualrain = std::max(1.0f, (janRain + aprRain + julRain + octRain) * 3.0f);

    const float biotemperature = calculateholdridgebiotemperature(summary);
    const float annualprecipitation = summary.annualrain;
    const float petratio = biotemperature <= 0.0f ? 1000.0f : (biotemperature * 58.93f) / annualprecipitation;

    if (biotemperature <= 0.0f)
        return biomeice;

    if (biotemperature < 1.5f)
    {
        if (petratio > 16.0f) return biomepolardesert;
        if (petratio > 4.0f) return biomepolardrytundra;
        if (petratio > 2.0f) return biomepolarmoisttundra;
        if (petratio > 0.5f) return biomepolarwettundra;
        return biomepolarraindtundra;
    }

    if (biotemperature < 3.0f)
    {
        if (petratio > 16.0f) return biomesubpolardesert;
        if (petratio > 4.0f) return biomesubpolardrytundra;
        if (petratio > 2.0f) return biomesubpolarmoisttundra;
        if (petratio > 0.5f) return biomesubpolarwettundra;
        return biomesubpolarraindtundra;
    }

    if (biotemperature < 6.0f)
    {
        if (petratio > 16.0f) return biomeborealdesert;
        if (petratio > 4.0f) return biomeborealdrybush;
        if (petratio > 1.0f) return biomeborealmoistforest;
        if (petratio > 0.25f) return biomeborealwetforest;
        return biomeborealrainforest;
    }

    if (biotemperature < 12.0f)
    {
        if (petratio > 16.0f) return biomecooltemperatedesert;
        if (petratio > 8.0f) return biomecooltemperatedesertbush;
        if (petratio > 2.0f) return biomecooltemperatesteppe;
        if (petratio > 1.0f) return biomecooltemperatemoistforest;
        if (petratio > 0.25f) return biomecooltemperatewetforest;
        return biomecooltemperaterainforest;
    }

    if (biotemperature < 18.0f)
    {
        if (petratio > 16.0f) return biomewarmtemperatedesert;
        if (petratio > 8.0f) return biomewarmtemperatedesertbush;
        if (petratio > 4.0f) return biomewarmtemperatethornsteppe;
        if (petratio > 2.0f) return biomewarmtemperatedryforest;
        if (petratio > 1.0f) return biomewarmtemperatemoistforest;
        if (petratio > 0.25f) return biomewarmtemperatewetforest;
        return biomewarmtemperaterainforest;
    }

    if (biotemperature < 24.0f)
    {
        if (petratio > 16.0f) return biomesubtropicaldesert;
        if (petratio > 8.0f) return biomesubtropicaldesertbush;
        if (petratio > 4.0f) return biomesubtropicalthornsteppe;
        if (petratio > 2.0f) return biomesubtropicaldryforest;
        if (petratio > 1.0f) return biomesubtropicalmoistforest;
        if (petratio > 0.25f) return biomesubtropicalwetforest;
        return biomesubtropicalrainforest;
    }

    if (petratio > 16.0f) return biometropicaldesert;
    if (petratio > 8.0f) return biometropicaldesertbush;
    if (petratio > 4.0f) return biometropicalthornsteppe;
    if (petratio > 2.0f) return biometropicalverydryforest;
    if (petratio > 1.0f) return biometropicaldryforest;
    if (petratio > 0.5f) return biometropicalmoistforest;
    if (petratio > 0.25f) return biometropicalwetforest;
    return biometropicalrainforest;
}

int calculateholdridgebiomefromseasonal(const seasonalclimatesummary& summary)
{
    return calculateholdridgebiomeinternal(
        summary.temps[seasonjanuary],
        summary.temps[seasonapril],
        summary.temps[seasonjuly],
        summary.temps[seasonoctober],
        summary.rains[seasonjanuary],
        summary.rains[seasonapril],
        summary.rains[seasonjuly],
        summary.rains[seasonoctober]);
}

short calculateclimatefromseasonal(int elev, int sealevel, const seasonalclimatesummary& summary)
{
    string group, preptype, heattype;

    if (summary.maxtemp <= 10.0f)
    {
        if (summary.maxtemp >= 0.0f)
            group = "ET";

        if (summary.maxtemp < 0.0f)
            group = "EF";
    }

    if (group == "")
    {
        float precthreshold = summary.meanannualtemp * 20.0f;

        if (summary.warmhalffraction >= 0.7f)
            precthreshold = precthreshold + 280.0f;
        else if (summary.warmhalffraction >= 0.3f)
            precthreshold = precthreshold + 140.0f;

        if (summary.annualrain < precthreshold * 0.5f)
            group = "BW";

        if (group == "" && summary.annualrain <= precthreshold)
            group = "BS";
    }

    if (group == "" && summary.mintemp >= 18.0f)
        group = "A";

    if (group == "" && summary.mintemp > -3.0f && summary.mintemp < 18.0f)
        group = "C";

    if (group == "" && summary.mintemp <= -3.0f)
        group = "D";

    if (group == "A")
    {
        if (summary.minrain >= 60.0f)
            preptype = "f";

        if (preptype == "" && summary.minrain >= (100.0f - (summary.annualrain / 25.0f)))
            preptype = "m";

        if (preptype == "")
            preptype = summary.driestseasoniswarm ? "s" : "w";
    }

    if (group == "C" || group == "D")
    {
        if (summary.driestwarmrain < summary.wettestcoldrain / 2.5f && summary.driestwarmrain < 35.0f)
            preptype = "s";

        if (preptype == "" && climatekoppen::isWinterDry(
                summary.driestcoldrain, summary.wettestwarmrain))
            preptype = "w";

        if (preptype == "")
            preptype = "f";
    }

    if (group == "BW" || group == "BS")
    {
        if (summary.meanannualtemp >= 18.0f)
            heattype = "h";

        if (summary.meanannualtemp < 18.0f)
            heattype = "k";
    }

    if (group != "A" && group != "BW" && group != "BS" && group != "ET" && group != "EF")
    {
        if (summary.maxtemp >= 22.0f)
            heattype = "a";
        else if (summary.secondwarmesttemp >= 8.0f)
            heattype = "b";
        else if (summary.mintemp <= -38.0f)
            heattype = "d";
        else
            heattype = "c";
    }

    const string climate = group + preptype + heattype;

    if (climate == "Af") return 1;
    if (climate == "Am") return 2;
    if (climate == "Aw") return 3;
    if (climate == "As") return 4;
    if (climate == "BWh") return 5;
    if (climate == "BWk") return 6;
    if (climate == "BSh") return 7;
    if (climate == "BSk") return 8;
    if (climate == "Csa") return 9;
    if (climate == "Csb") return 10;
    if (climate == "Csc") return 11;
    if (climate == "Cwa") return 12;
    if (climate == "Cwb") return 13;
    if (climate == "Cwc") return 14;
    if (climate == "Cfa") return 15;
    if (climate == "Cfb") return 16;
    if (climate == "Cfc") return 17;
    if (climate == "Dsa") return 18;
    if (climate == "Dsb") return 19;
    if (climate == "Dsc") return 20;
    if (climate == "Dsd") return 21;
    if (climate == "Dwa") return 22;
    if (climate == "Dwb") return 23;
    if (climate == "Dwc") return 24;
    if (climate == "Dwd") return 25;
    if (climate == "Dfa") return 26;
    if (climate == "Dfb") return 27;
    if (climate == "Dfc") return 28;
    if (climate == "Dfd") return 29;
    if (climate == "ET") return 30;
    if (climate == "EF") return 31;

    return 0;
}
}

// This creates the climate map.

void createclimatemap(planet& world)
{
    int width = world.width();
    int height = world.height();

    short climate;

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            climate = getclimate(world, i, j);

            world.setclimate(i, j, climate);
        }
    }
}

void createbiomemap(planet& world)
{
    const int width = world.width();
    const int height = world.height();

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (isclassificationwater(world, i, j))
                world.setbiome(i, j, 0);
            else
                world.setbiome(i, j, calculateplanetbiome(world, i, j));
        }
    }
}

// This returns the climate type of the given point.

short getclimate(planet& world, int x, int y)
{
    if (isclassificationwater(world, x, y))
        return (0);

    return calculateplanetclimate(world, x, y);
}

// The same thing, but for the regional map.

short getclimate(region& region, int x, int y)
{
    if (region.sea(x, y) == 1 || region.truelake(x, y) != 0)
        return (0);

    int elev = region.map(x, y);
    int sealevel = region.sealevel();

    float wrain = (float)region.winterrain(x, y);
    float srain = (float)region.summerrain(x, y);
    float mintemp = (float)region.mintemp(x, y);
    float maxtemp = (float)region.maxtemp(x, y);

    short climate = calculateclimate(elev, sealevel, wrain, srain, mintemp, maxtemp);

    return climate;
}

int calculateholdridgebiome(float janTemp, float aprTemp, float julTemp, float octTemp, float janRain, float aprRain, float julRain, float octRain)
{
    return calculateholdridgebiomeinternal(janTemp, aprTemp, julTemp, octTemp, janRain, aprRain, julRain, octRain);
}

// This does the actual climate calculations for the above two functions.

short calculateclimate(int elev, int sealevel, float wrain, float srain, float mintemp, float maxtemp)
{
    float totalannualrain = 0;

    if (srain > wrain) // If there's more rain in summer, assume a shorter rainy season. The greater the imbalance, the shorter the season.
    {
        float factor = (float)wrain / (float)srain; // 0-1. The higher it is, the more balanced the distribution of rain.

        factor = factor * 12.0f;

        if (factor < 4.0)
            factor = 4.0f;

        int sfactor = (int)factor;

        int wfactor = 12 - sfactor;
        
        totalannualrain = wrain * wfactor + srain * sfactor;
    }
    else
        totalannualrain = (wrain + srain) * 6;

    float meanannualrain = (wrain + srain) / 2;
    float meanannualtemp = (mintemp + maxtemp) / 2;

    float minrain = wrain;
    float maxrain = srain;

    if (wrain > srain)
    {
        minrain = srain;
        maxrain = wrain;
    }

    string group, preptype, heattype; // These three variables define each climate type.

    // First, establish the group (first letter).

    if (maxtemp <= 10)
    {
        if (maxtemp >= 0)
            group = "ET";

        if (maxtemp < 0)
            group = "EF";
    }

    if (group == "")
    {
        float precthreshold = (maxtemp + mintemp) * 10;
        float checkpercent = (srain * 6) / totalannualrain;

        if (checkpercent >= 0.7f)
            precthreshold = precthreshold + 280.0f;

        if (checkpercent >= 0.3 && checkpercent < 0.7f)
            precthreshold = precthreshold + 140.0f;

        if (totalannualrain < precthreshold * 0.5f)
            group = "BW";

        if (totalannualrain >= precthreshold * 0.5 && totalannualrain <= precthreshold)
            group = "BS";
    }

    if (group == "" && mintemp >= 18)
        group = "A";

    if (group == "" && mintemp > -3 && mintemp < 18)
        group = "C";

    if (group == "" && mintemp <= -3)
        group = "D";

    // Now, establish the precipitation type (second letter).

    if (group == "A")
    {
        if (minrain >= 60)
            preptype = "f";

        if (preptype == "" && minrain >= (100 - (totalannualrain / 25))) //  meanannualrain>=(100-minrain))
            preptype = "m";

        if (preptype == "" && srain < 60)
            preptype = "s";

        if (preptype == "" && wrain < 60)
            preptype = "w";
    }

    if (group == "C" || group == "D")
    {
        if (srain < wrain / 3 && srain < 40)
            preptype = "s";

        if (preptype == "" && wrain < srain / 10)
            preptype = "w";

        if (preptype == "")
            preptype = "f";
    }

    // Now, establish the heat type (third letter).

    if (group == "BW" || group == "BS")
    {
        if (meanannualtemp >= 18)
            heattype = "h";

        if (meanannualtemp < 18)
            heattype = "k";
    }

    if (group != "A" && group != "BW" && group != "BS" && group != "ET" && group != "EF")
    {
        if (heattype == "" && maxtemp >= 22)
            heattype = "a";

        if (heattype == "" && maxtemp < 22 && maxtemp >= 14) // Should be: at least four months are >=10, but we can't track that accurately on our model.
            heattype = "b";

        if (heattype == "" && mintemp <= -38)
            heattype = "d";

        if (heattype == "")
            heattype = "c";

    }

    string climate = group + preptype + heattype;

    if (climate == "Af") // Af
        return 1;

    if (climate == "Am") // Am
        return 2;

    if (climate == "Aw") // Aw
        return 3;

    if (climate == "As") // As
        return 4;

    if (climate == "BWh") // Bwh
        return 5;

    if (climate == "BWk") // BWk
        return 6;

    if (climate == "BSh") // BSh
        return 7;

    if (climate == "BSk") // BSk
        return 8;

    if (climate == "Csa") // Csa
        return 9;

    if (climate == "Csb") // Csb
        return 10;

    if (climate == "Csc") // Csc
        return 11;

    if (climate == "Cwa") // Cwa
        return 12;

    if (climate == "Cwb") // Cwb
        return 13;

    if (climate == "Cwc") // Cwc
        return 14;

    if (climate == "Cfa") // Cfa
        return 15;

    if (climate == "Cfb") // Cfb
        return 16;

    if (climate == "Cfc") // Cfc
        return 17;

    if (climate == "Dsa") // Dsa
        return 18;

    if (climate == "Dsb") // Dsb
        return 19;

    if (climate == "Dsc") // Dsc
        return 20;

    if (climate == "Dsd") // Dsd
        return 21;

    if (climate == "Dwa") // Dwa
        return 22;

    if (climate == "Dwb") // Dwb
        return 23;

    if (climate == "Dwc") // Dwc
        return 24;

    if (climate == "Dwd") // Dwd
        return 25;

    if (climate == "Dfa") // Dfa
        return 26;

    if (climate == "Dfb") // Dfb
        return 27;

    if (climate == "Dfc") // Dfc
        return 28;

    if (climate == "Dfd") // Dfd
        return 29;

    if (climate == "ET") // ET
        return 30;

    if (climate == "EF") // EF
        return 31;

    return 0;
}

// This function gives the names of the climate types.

string getclimatename(short climate)
{
    if (climate == 1) // Af
        return "Tropical rainforest";

    if (climate == 2) // Am
        return "Monsoon";

    if (climate == 3) // Aw
        return "Savannah";

    if (climate == 4) // As
        return "Savannah";

    if (climate == 5) // Bwh
        return "Hot desert";

    if (climate == 6) // BWk
        return "Cold desert";

    if (climate == 7) // BSh
        return "Hot semi-arid";

    if (climate == 8) // BSk
        return "Cold steppe";

    if (climate == 9) // Csa
        return "Hot, dry-summer Mediterranean";

    if (climate == 10) // Csb
        return "Warm, dry-summer Mediterranean";

    if (climate == 11) // Csc
        return "Cold, dry-summer Mediterranean";

    if (climate == 12) // Cwa
        return "Dry-winter humid subtropical";

    if (climate == 13) // Cwb
        return "Dry-winter subtropical highland";

    if (climate == 14) // Cwc
        return "Dry-winter subpolar oceanic";

    if (climate == 15) // Cfa
        return "Humid subtropical";

    if (climate == 16) // Cfb
        return "Temperate oceanic";

    if (climate == 17) // Cfc
        return "Subpolar oceanic";

    if (climate == 18) // Dsa
        return "Mediterranean-influenced hot-summer humid continental";

    if (climate == 19) // Dsb
        return "Mediterranean-influenced warm-summer humid continental";

    if (climate == 20) // Dsc
        return "Mediterranean-influenced subarctic";

    if (climate == 21) // Dsd
        return "Mediterranean-influenced extremely cold subarctic";

    if (climate == 22) // Dwa
        return "Monsoon-influenced hot-summer humid continental";

    if (climate == 23) // Dwb
        return "Monsoon-influenced warm-summer humid continental";

    if (climate == 24) // Dwc
        return "Monsoon-influenced subarctic";

    if (climate == 25) // Dwd
        return "Monsoon-influenced extremely cold subarctic";

    if (climate == 26) // Dfa
        return "Hot-summer humid continental";

    if (climate == 27) // Dfb
        return "Warm-summer humid continental";

    if (climate == 28) // Dfc
        return "Subarctic";

    if (climate == 29) // Dfd
        return "Extremely cold subarctic";

    if (climate == 30) // ET
        return "Tundra";

    if (climate == 31) // EF
        return "Frost";

    return "";
}

// This function gives the codes of the climate types.

string getclimatecode(short climate)
{
    if (climate == 1) // Af
        return "Af";

    if (climate == 2) // Am
        return "Am";

    if (climate == 3) // Aw
        return "Aw";

    if (climate == 4) // As
        return "As";

    if (climate == 5) // Bwh
        return "BWh";

    if (climate == 6) // BWk
        return "BWk";

    if (climate == 7) // BSh
        return "BSh";

    if (climate == 8) // BSk
        return "BSk";

    if (climate == 9) // Csa
        return "Csa";

    if (climate == 10) // Csb
        return "Csb";

    if (climate == 11) // Csc
        return "Csc";

    if (climate == 12) // Cwa
        return "Cwa";

    if (climate == 13) // Cwb
        return "Cwb";

    if (climate == 14) // Cwc
        return "Cwc";

    if (climate == 15) // Cfa
        return "Cfa";

    if (climate == 16) // Cfb
        return "Cfb";

    if (climate == 17) // Cfc
        return "Cfc";

    if (climate == 18) // Dsa
        return "Dsa";

    if (climate == 19) // Dsb
        return "Dsb";

    if (climate == 20) // Dsc
        return "Dsc";

    if (climate == 21) // Dsd
        return "Dsd";

    if (climate == 22) // Dwa
        return "Dwa";

    if (climate == 23) // Dwb
        return "Dwb";

    if (climate == 24) // Dwc
        return "Dwc";

    if (climate == 25) // Dwd
        return "Dwd";

    if (climate == 26) // Dfa
        return "Dfa";

    if (climate == 27) // Dfb
        return "Dfb";

    if (climate == 28) // Dfc
        return "Dfc";

    if (climate == 29) // Dfd
        return "Dfd";

    if (climate == 30) // ET
        return "ET";

    if (climate == 31) // EF
        return "EF";

    return "";
}
