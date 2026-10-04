//
//  planet.cpp
//  Undiscovered Worlds
//
//  Created by Jonathan Hill on 22/07/2019.
//
//  Please see functions.hpp for notes.

#include <iostream>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <utility>

#include "planet.hpp"
#include "climate_weather.hpp"
#include "climate_grid.hpp"
#include "climate_hydrology.hpp"
#include "generation_tuning.hpp"
#include "functions.hpp"

namespace
{
BoundaryType decodeboundarytype(std::uint8_t value)
{
    switch (value)
    {
    case static_cast<std::uint8_t>(BoundaryType::convergent):
        return BoundaryType::convergent;
    case static_cast<std::uint8_t>(BoundaryType::divergent):
        return BoundaryType::divergent;
    case static_cast<std::uint8_t>(BoundaryType::transform):
        return BoundaryType::transform;
    case static_cast<std::uint8_t>(BoundaryType::passive_margin):
        return BoundaryType::passive_margin;
    case static_cast<std::uint8_t>(BoundaryType::none):
    default:
        return BoundaryType::none;
    }
}

GeologicRegime decodegeologicregime(std::uint8_t value)
{
    switch (value)
    {
    case static_cast<std::uint8_t>(GeologicRegime::convergent_arc):
        return GeologicRegime::convergent_arc;
    case static_cast<std::uint8_t>(GeologicRegime::continent_collision):
        return GeologicRegime::continent_collision;
    case static_cast<std::uint8_t>(GeologicRegime::divergent_rift):
        return GeologicRegime::divergent_rift;
    case static_cast<std::uint8_t>(GeologicRegime::transform):
        return GeologicRegime::transform;
    case static_cast<std::uint8_t>(GeologicRegime::passive_margin):
        return GeologicRegime::passive_margin;
    case static_cast<std::uint8_t>(GeologicRegime::mid_ocean_ridge):
        return GeologicRegime::mid_ocean_ridge;
    case static_cast<std::uint8_t>(GeologicRegime::trench_adjacent):
        return GeologicRegime::trench_adjacent;
    case static_cast<std::uint8_t>(GeologicRegime::stable):
    default:
        return GeologicRegime::stable;
    }
}

DeformingRegionType decodedeformingregiontype(std::uint8_t value)
{
    switch (value)
    {
    case static_cast<std::uint8_t>(DeformingRegionType::continental_rift):
        return DeformingRegionType::continental_rift;
    case static_cast<std::uint8_t>(DeformingRegionType::diffuse_collision):
        return DeformingRegionType::diffuse_collision;
    case static_cast<std::uint8_t>(DeformingRegionType::none):
    default:
        return DeformingRegionType::none;
    }
}
}

//#define ENABLE_PROFILER
#ifdef ENABLE_PROFILER
#include "profiler.h"
#endif

planet::planet() //constructor
{
    body_.itswidth = 599;
    body_.itsheight = 299;
    geosphere_.itstectonictimeoriginstep = 0;
    geosphere_.itstectonictimemyr = 0.0f;
    geosphere_.itstectonicdeltatimemyr = 0.0f;
    geosphere_.itstectonicplatecount = 0;
    geosphere_.itstectonicsealevelm = 0;
    resizegrids();
}

planet::~planet()
{
}

void planet::resizegrids()
{
    climate_.resize(body_.itswidth, body_.itsheight);
    const int columns = body_.itswidth + 1;
    const int rows = body_.itsheight + 1;
    const bool resized = geosphere_.mapnom.columns() != columns || geosphere_.mapnom.rows() != rows;
    hydrosphere_.rivermapdir.resize(columns, rows);
    hydrosphere_.rivermapjan.resize(columns, rows);
    hydrosphere_.rivermapjul.resize(columns, rows);
    hydrosphere_.lakemap.resize(columns, rows);
    geosphere_.roughnessmap.resize(columns, rows);
    geosphere_.mountainridges.resize(columns, rows);
    geosphere_.mountainheights.resize(columns, rows);
    geosphere_.craterrims.resize(columns, rows);
    geosphere_.cratercentres.resize(columns, rows);
    geosphere_.mapnom.resize(columns, rows);
    hydrosphere_.tidalmap.resize(columns, rows);
    hydrosphere_.riftlakemapsurface.resize(columns, rows);
    hydrosphere_.riftlakemapbed.resize(columns, rows);
    hydrosphere_.lakestartmap.resize(columns, rows);
    hydrosphere_.specials.resize(columns, rows);
    geosphere_.geologicregimemap.resize(columns, rows);
    geosphere_.tectonicconvergencemap.resize(columns, rows);
    geosphere_.tectonicdivergencemap.resize(columns, rows);
    geosphere_.tectonicshearmap.resize(columns, rows);
    geosphere_.tectoniccrustagemyrmap.resize(columns, rows);
    geosphere_.tectoniccrustthicknessmap.resize(columns, rows);
    geosphere_.tectoniccrustclassmap.resize(columns, rows);
    geosphere_.tectonicuplifttendencymap.resize(columns, rows);
    geosphere_.tectonicsubsidencetendencymap.resize(columns, rows);
    geosphere_.tectonicaccumulatedstrainmap.resize(columns, rows);
    geosphere_.tectonicboundarytypemap.resize(columns, rows);
    geosphere_.tectonicboundarydistancemap.resize(columns, rows);
    geosphere_.tectonicboundarysegmentidmap.resize(columns, rows);
    geosphere_.tectonicnearestboundaryidmap.resize(columns, rows);
    geosphere_.tectonicboundaryhistorymap.resize(columns, rows);
    geosphere_.tectonicdeformingregionidmap.resize(columns, rows);
    geosphere_.tectonicdeformingregiontypemap.resize(columns, rows);
    geosphere_.tectonicdeformationratemap.resize(columns, rows);
    geosphere_.tectonicdeformationvelocityxmap.resize(columns, rows);
    geosphere_.tectonicdeformationvelocityymap.resize(columns, rows);
    hydrosphere_.basinclassmap.resize(columns, rows);
    geosphere_.erosionpotentialmap.resize(columns, rows);
    geosphere_.depositionpotentialmap.resize(columns, rows);
    geosphere_.floodplainfertilitymap.resize(columns, rows);
    geosphere_.metalorereservemap.resize(columns, rows);
    geosphere_.placerreservemap.resize(columns, rows);
    geosphere_.evaporitereservemap.resize(columns, rows);
    geosphere_.volcanicreservemap.resize(columns, rows);
    biosphere_.fisheryreservemap.resize(columns, rows);
    anthroposphere_.settlement_suitability.resize(columns, rows);
    anthroposphere_.social_infrastructure.resize(columns, rows);
    anthroposphere_.agricultural_capacity.resize(columns, rows);
    anthroposphere_.route_traffic.resize(columns, rows);
    anthroposphere_.river_access.resize(columns, rows);
    anthroposphere_.harbor_score.resize(columns, rows);
    anthroposphere_.owner_settlement_id.resize(columns, rows);
    anthroposphere_.owner_polity_id.resize(columns, rows);
    geosphere_.extraelevmap.resize(columns, rows);
    hydrosphere_.deltamapdir.resize(columns, rows);
    hydrosphere_.deltamapjan.resize(columns, rows);
    hydrosphere_.deltamapjul.resize(columns, rows);
    geosphere_.islandmap.resize(columns, rows);
    geosphere_.mountainislandmap.resize(columns, rows);
    geosphere_.noisemap.resize(columns, rows);
    geosphere_.oceanridgemap.resize(columns, rows);
    geosphere_.oceanridgeheightmap.resize(columns, rows);
    geosphere_.oceanriftmap.resize(columns, rows);
    geosphere_.oceanridgeoffsetmap.resize(columns, rows);
    geosphere_.oceanridgeanglemap.resize(columns, rows);
    geosphere_.volcanomap.resize(columns, rows);
    geosphere_.stratomap.resize(columns, rows);
    presentation_.noshademap.resize(columns, rows);
    presentation_.testmap.resize(columns, rows);
    geosphere_.shelfmap.resize(columns, rows);
    // Resized grids are zero; give them the values clear() sets, as the former
    // fixed arrays kept them across size changes (imports size the world after clear()).
    if (resized)
    {
        geosphere_.geologicregimemap.fill(static_cast<std::uint8_t>(GeologicRegime::stable));
        geosphere_.tectoniccrustclassmap.fill(static_cast<std::uint8_t>(CrustClass::none));
        geosphere_.tectonicboundarytypemap.fill(static_cast<std::uint8_t>(BoundaryType::none));
        geosphere_.tectonicdeformingregiontypemap.fill(static_cast<std::uint8_t>(DeformingRegionType::none));
        hydrosphere_.basinclassmap.fill(static_cast<std::uint8_t>(BasinClass::none));
        anthroposphere_.owner_settlement_id.fill(-1);
        anthroposphere_.owner_polity_id.fill(-1);
    }
}

void planet::cleartectonicprovenanceinternal()
{
    geosphere_.itstectonictimeoriginstep = 0;
    geosphere_.itstectonictimemyr = 0.0f;
    geosphere_.itstectonicdeltatimemyr = 0.0f;
    geosphere_.itstectonicplatecount = 0;
    geosphere_.itstectonicsealevelm = 0;
    geosphere_.tectonicboundarysegmentlist.clear();
    geosphere_.tectonicdeformingregionlist.clear();
    geosphere_.mineraldepositlist.clear();
    geosphere_.mineralbudgetstate = {};

    parallelforrows(0, body_.itswidth, [&](int startx, int endx)
    {
        for (int i = startx; i <= endx; i++)
        {
            for (int j = 0; j <= body_.itsheight; j++)
            {
                geosphere_.geologicregimemap[i][j] = static_cast<std::uint8_t>(GeologicRegime::stable);
                geosphere_.tectonicconvergencemap[i][j] = 0;
                geosphere_.tectonicdivergencemap[i][j] = 0;
                geosphere_.tectonicshearmap[i][j] = 0;
                geosphere_.tectoniccrustagemyrmap[i][j] = 0.0f;
                geosphere_.tectoniccrustthicknessmap[i][j] = 0.0f;
                geosphere_.tectoniccrustclassmap[i][j] = static_cast<std::uint8_t>(CrustClass::none);
                geosphere_.tectonicuplifttendencymap[i][j] = 0.0f;
                geosphere_.tectonicsubsidencetendencymap[i][j] = 0.0f;
                geosphere_.tectonicaccumulatedstrainmap[i][j] = 0.0f;
                geosphere_.tectonicboundarytypemap[i][j] = static_cast<std::uint8_t>(BoundaryType::none);
                geosphere_.tectonicboundarydistancemap[i][j] = 0;
                geosphere_.tectonicboundarysegmentidmap[i][j] = 0;
                geosphere_.tectonicnearestboundaryidmap[i][j] = 0;
                geosphere_.tectonicboundaryhistorymap[i][j] = 0.0f;
                geosphere_.tectonicdeformingregionidmap[i][j] = 0;
                geosphere_.tectonicdeformingregiontypemap[i][j] = static_cast<std::uint8_t>(DeformingRegionType::none);
                geosphere_.tectonicdeformationratemap[i][j] = 0.0f;
                geosphere_.tectonicdeformationvelocityxmap[i][j] = 0.0f;
                geosphere_.tectonicdeformationvelocityymap[i][j] = 0.0f;
            }
        }
    }, 64);
}

void planet::cleartectonicprovenance()
{
    cleartectonicprovenanceinternal();
}

const std::vector<std::uint8_t>& planet::weatheranomalystate() const
{
    return climate_.weatheranomalystate();
}

void planet::setweatheranomalystate(std::vector<std::uint8_t> state)
{
    climate_.setweatheranomalystate(std::move(state));
}

bool planet::advanceweather(float elapsedSeconds, int horizontalCells)
{
    climate_.bindinputs(climateinputsview());
    return climate_.advanceweather(elapsedSeconds, horizontalCells);
}

void planet::writeshortvectordata(ofstream& outfile, const std::vector<short>& arr)
{
    outfile.write(reinterpret_cast<const char*>(arr.data()), arr.size() * sizeof(short));
}

// The save format keeps the climate fields in their former int16 encodings.
void planet::readlegacyvectordata(ifstream& infile, std::vector<float>& arr, float scale)
{
    std::vector<short> stored(arr.size());
    readshortvectordata(infile, stored);
    for (std::size_t index = 0; index < arr.size(); index++)
        arr[index] = static_cast<float>(stored[index]) / scale;
}

void planet::writelegacydata(ofstream& outfile, const worldgrid<float>& arr)
{
    worldgrid<short> stored(arr.columns(), arr.rows());
    const std::vector<short> encoded = climateworld::legacyfield(arr);
    std::copy(encoded.begin(), encoded.end(), stored.data());
    writedata(outfile, stored);
}

void planet::readlegacydata(ifstream& infile, worldgrid<float>& arr)
{
    worldgrid<short> stored(arr.columns(), arr.rows());
    readdata(infile, stored);
    for (std::size_t index = 0; index < arr.size(); index++)
        arr.data()[index] = static_cast<float>(stored.data()[index]);
}

void planet::readshortvectordata(ifstream& infile, std::vector<short>& arr)
{
    infile.read(reinterpret_cast<char*>(arr.data()), arr.size() * sizeof(short));
}

void planet::syncseasonalclimatefromlegacy()
{
    climate_.bindinputs(climateinputsview());
    climate_.syncseasonalclimatefromlegacy();
}

bool planet::outline(int x, int y) const
{
    if (y<1 || y>body_.itsheight - 1)
        return 0;

    if (sea(x, y) == 0)
    {
        if (seawrap(x - 1, y) == 1)
            return 1;

        if (seawrap(x, y - 1) == 1)
            return 1;

        if (seawrap(x + 1, y) == 1)
            return 1;

        if (seawrap(x, y + 1) == 1)
            return 1;
    }

    return 0;
}

bool planet::coast(int x, int y) const
{
    if (y<1 || y>body_.itsheight - 1 || x<0 || x>body_.itswidth)
        return 0;

    if (sea(x, y) == 1)
    {
        if (seawrap(x - 1, y) == 0)
            return 1;

        if (seawrap(x, y - 1) == 0)
            return 1;

        if (seawrap(x + 1, y) == 0)
            return 1;

        if (seawrap(x, y + 1) == 0)
            return 1;
    }

    return 0;
}

void planet::longitude(int x, int &degrees, int &minutes, int &seconds, bool& negative) const
{
    if ( x<0 || x>body_.itswidth)
        return;

    float fx = (float)x;

    float worldwidth = (float)body_.itswidth + 1.0f;

    float pixelsperlong = worldwidth / 360.0f;

    float longitude = fx / pixelsperlong;

    longitude = longitude - 180.0f;

    degrees = (int)longitude;

    float dec = (longitude - (float)degrees) * 60.0f;

    minutes = (int)dec;

    float dec2 = (dec - (float)minutes) * 60;

    seconds = (int)dec2;

    if (degrees < 0)
        degrees = 0 - degrees;

    if (minutes < 0)
        minutes = 0 - minutes;

    if (seconds < 0)
        seconds = 0 - seconds;

    negative = 0;

    if (fx < worldwidth / 2.f)
        negative = 1;

    return;
}

void planet::latitude(int y, int& degrees, int& minutes, int& seconds, bool &negative) const
{
    if (y<0 || y>body_.itsheight)
        return;

    bool hemisphere = 0; // 1 for southern

    float fy = (float)y;

    float worldheight = (float)body_.itsheight;
    float worldhalfheight = worldheight / 2.0f;

    float pixelsperlat = worldhalfheight / 90.0f;

    float latitude;

    if (fy <= worldhalfheight)
        latitude = (worldhalfheight - fy) / pixelsperlat;
    else
    {
        latitude = (fy - worldhalfheight) / pixelsperlat;
        hemisphere = 1;
    }

    degrees = (int)latitude;

    float dec = (latitude - (float)degrees) * 60.0f;

    minutes = (int)dec;

    float dec2 = (dec - (float)minutes) * 60;

    seconds = (int)dec2;

    negative = 0;

    if (hemisphere == 1)
        negative = 1;

    return;
}

int planet::reverselatitude(int lat) const
{
    float flat = 90.0f - (float)lat;

    float equator = (float)body_.itsheight / 2.0f;

    float pixelsperlat = equator / 90.0f;

    float y = flat * pixelsperlat;

    return (int)y;
}

// slightly more complicated accessor functions

bool planet::seawrap(int x, int y) const
{
    UW_FIELD_READ(base_elevation);
    UW_FIELD_READ(sea_level);
    UW_FIELD_READ(lake_surface);
    x = wrapx(x);
    y = clipy(y);

    if (geosphere_.mapnom[x][y] <= geosphere_.itssealevel && hydrosphere_.lakemap[x][y] == 0)
        return 1;

    else
        return 0;
}

bool planet::outlinewrap(int x, int y) const
{
    x = wrapx(x);

    if (y<1 || y>body_.itsheight - 1)
        return 0;

    if (sea(x, y) == 0)
    {
        if (seawrap(x - 1, y) == 1)
            return 1;

        if (seawrap(x, y - 1) == 1)
            return 1;

        if (seawrap(x + 1, y) == 1)
            return 1;

        if (seawrap(x, y + 1) == 1)
            return 1;
    }

    return 0;
}

int planet::mountainheightwrap(int x, int y) const
{
    UW_FIELD_READ(mountain_height);
    x = wrapx(x);
    y = clipy(y);

    return geosphere_.mountainheights[x][y];
}

// Other public functions.

void planet::clear()
{
    resizegrids();
    climate_.weatheranomalystatebytes.clear();
    geosphere_.itstectonictimeoriginstep = 0;
    geosphere_.itstectonictimemyr = 0.0f;
    geosphere_.itstectonicdeltatimemyr = 0.0f;
    geosphere_.itstectonicplatecount = 0;
    geosphere_.itstectonicsealevelm = 0;
    geosphere_.tectonicboundarysegmentlist.clear();
    geosphere_.tectonicdeformingregionlist.clear();
    geosphere_.mineraldepositlist.clear();
    geosphere_.mineralbudgetstate = {};

    parallelforrows(0, body_.itswidth, [&](int startx, int endx)
    {
        for (int i = startx; i <= endx; i++) // Set all the maps to 0.
        {
            for (int j = 0; j <= body_.itsheight; j++)
            {
                climate_.jantempmap[i][j] = 0;
                climate_.jultempmap[i][j] = 0;
                climate_.climatemap[i][j] = 0;
                climate_.biomemap[i][j] = 0;
                climate_.janrainmap[i][j] = 0;
                climate_.julrainmap[i][j] = 0;
                climate_.janmountainrainmap[i][j] = 0;
                climate_.julmountainrainmap[i][j] = 0;
                climate_.janmountainraindirmap[i][j] = 0;
                climate_.julmountainraindirmap[i][j] = 0;
                climate_.seaicemap[i][j] = 0;
                hydrosphere_.rivermapdir[i][j] = 0;
                hydrosphere_.rivermapjan[i][j] = 0;
                hydrosphere_.rivermapjul[i][j] = 0;
                climate_.windmap[i][j] = 0;
                hydrosphere_.lakemap[i][j] = 0;
                geosphere_.roughnessmap[i][j] = 0;
                geosphere_.mountainridges[i][j] = 0;
                geosphere_.mountainheights[i][j] = 0;
                geosphere_.craterrims[i][j] = 0;
                geosphere_.cratercentres[i][j] = 0;
                geosphere_.mapnom[i][j] = 0;
                hydrosphere_.tidalmap[i][j] = 0;
                hydrosphere_.riftlakemapsurface[i][j] = 0;
                hydrosphere_.riftlakemapbed[i][j] = 0;
                hydrosphere_.specials[i][j] = 0;
                geosphere_.geologicregimemap[i][j] = static_cast<std::uint8_t>(GeologicRegime::stable);
                geosphere_.tectonicconvergencemap[i][j] = 0;
                geosphere_.tectonicdivergencemap[i][j] = 0;
                geosphere_.tectonicshearmap[i][j] = 0;
                geosphere_.tectoniccrustagemyrmap[i][j] = 0.0f;
                geosphere_.tectoniccrustthicknessmap[i][j] = 0.0f;
                geosphere_.tectoniccrustclassmap[i][j] = static_cast<std::uint8_t>(CrustClass::none);
                geosphere_.tectonicuplifttendencymap[i][j] = 0.0f;
                geosphere_.tectonicsubsidencetendencymap[i][j] = 0.0f;
                geosphere_.tectonicaccumulatedstrainmap[i][j] = 0.0f;
                geosphere_.tectonicboundarytypemap[i][j] = static_cast<std::uint8_t>(BoundaryType::none);
                geosphere_.tectonicboundarydistancemap[i][j] = 0;
                geosphere_.tectonicboundarysegmentidmap[i][j] = 0;
                geosphere_.tectonicnearestboundaryidmap[i][j] = 0;
                geosphere_.tectonicboundaryhistorymap[i][j] = 0.0f;
                geosphere_.tectonicdeformingregionidmap[i][j] = 0;
                geosphere_.tectonicdeformingregiontypemap[i][j] = static_cast<std::uint8_t>(DeformingRegionType::none);
                geosphere_.tectonicdeformationratemap[i][j] = 0.0f;
                geosphere_.tectonicdeformationvelocityxmap[i][j] = 0.0f;
                geosphere_.tectonicdeformationvelocityymap[i][j] = 0.0f;
                hydrosphere_.basinclassmap[i][j] = static_cast<std::uint8_t>(BasinClass::none);
                geosphere_.erosionpotentialmap[i][j] = 0;
                geosphere_.depositionpotentialmap[i][j] = 0;
                geosphere_.floodplainfertilitymap[i][j] = 0;
                geosphere_.metalorereservemap[i][j] = 0;
                geosphere_.placerreservemap[i][j] = 0;
                geosphere_.evaporitereservemap[i][j] = 0;
                geosphere_.volcanicreservemap[i][j] = 0;
                biosphere_.fisheryreservemap[i][j] = 0;
                anthroposphere_.settlement_suitability[i][j] = 0;
                anthroposphere_.social_infrastructure[i][j] = 0;
                anthroposphere_.agricultural_capacity[i][j] = 0;
                anthroposphere_.route_traffic[i][j] = 0;
                anthroposphere_.river_access[i][j] = 0;
                anthroposphere_.harbor_score[i][j] = 0;
                anthroposphere_.owner_settlement_id[i][j] = -1;
                anthroposphere_.owner_polity_id[i][j] = -1;
                geosphere_.extraelevmap[i][j] = 0;
                hydrosphere_.deltamapdir[i][j] = 0;
                hydrosphere_.deltamapjan[i][j] = 0;
                hydrosphere_.deltamapjul[i][j] = 0;
                geosphere_.oceanridgemap[i][j] = 0;
                geosphere_.oceanridgeheightmap[i][j] = 0;
                geosphere_.oceanriftmap[i][j] = 0;
                geosphere_.oceanridgeoffsetmap[i][j] = 0;
                geosphere_.islandmap[i][j] = 0;
                presentation_.noshademap[i][j] = 0;
                geosphere_.mountainislandmap[i][j] = 0;
                geosphere_.volcanomap[i][j] = 0;
                geosphere_.stratomap[i][j] = 0;
                geosphere_.noisemap[i][j] = 0;
                presentation_.testmap[i][j] = 0;
            }
        }
    }, 64);
    geosphere_.shelfmap.fill(false);
    // Subtropical-high rows span the full 2048-column capacity in saves.
    climate_.horselats.fill(0);

    for (int i = 0; i < MAXCRATERS; i++)
    {
        geosphere_.cratercentreslist[i].w = 0;
        geosphere_.cratercentreslist[i].x = 0;
        geosphere_.cratercentreslist[i].y = 0;
        geosphere_.cratercentreslist[i].z = 0;
    }

    clearsocialstate();
}

void planet::smoothnom(int amount)
{
    smooth(geosphere_.mapnom, amount, 1, 1);
}

void planet::smoothextraelev(int amount)
{
    smoothoverland(geosphere_.extraelevmap, amount, 0);
}

void planet::shiftterrain(int offset)
{
    // Map columns move left by the offset; deposits move with them.
    const double columns = static_cast<double>(body_.itswidth + 1);
    for (mineralsystems::Deposit& deposit : geosphere_.mineraldepositlist)
    {
        double u = std::fmod(deposit.u - offset / columns, 1.0);
        u = u < 0.0 ? u + 1.0 : u;
        deposit.u = u >= 1.0 ? 0.0 : u;
    }

    shift(geosphere_.mapnom, offset);
    shift(geosphere_.geologicregimemap, offset);
    shift(geosphere_.tectonicconvergencemap, offset);
    shift(geosphere_.tectonicdivergencemap, offset);
    shift(geosphere_.tectonicshearmap, offset);
    shift(geosphere_.tectoniccrustagemyrmap, offset);
    shift(geosphere_.tectoniccrustthicknessmap, offset);
    shift(geosphere_.tectoniccrustclassmap, offset);
    shift(geosphere_.tectonicuplifttendencymap, offset);
    shift(geosphere_.tectonicsubsidencetendencymap, offset);
    shift(geosphere_.tectonicaccumulatedstrainmap, offset);
    shift(geosphere_.tectonicboundarytypemap, offset);
    shift(geosphere_.tectonicboundarydistancemap, offset);
    shift(geosphere_.tectonicboundarysegmentidmap, offset);
    shift(geosphere_.tectonicnearestboundaryidmap, offset);
    shift(geosphere_.tectonicboundaryhistorymap, offset);
    shift(geosphere_.tectonicdeformingregionidmap, offset);
    shift(geosphere_.tectonicdeformingregiontypemap, offset);
    shift(geosphere_.tectonicdeformationratemap, offset);
    shift(geosphere_.tectonicdeformationvelocityxmap, offset);
    shift(geosphere_.tectonicdeformationvelocityymap, offset);
    shift(hydrosphere_.basinclassmap, offset);
    shift(geosphere_.erosionpotentialmap, offset);
    shift(geosphere_.depositionpotentialmap, offset);
    shift(geosphere_.floodplainfertilitymap, offset);
    shift(geosphere_.metalorereservemap, offset);
    shift(geosphere_.placerreservemap, offset);
    shift(geosphere_.evaporitereservemap, offset);
    shift(geosphere_.volcanicreservemap, offset);
    shift(biosphere_.fisheryreservemap, offset);
    shift(anthroposphere_.settlement_suitability, offset);
    shift(anthroposphere_.social_infrastructure, offset);
    shift(anthroposphere_.agricultural_capacity, offset);
    shift(anthroposphere_.route_traffic, offset);
    shift(anthroposphere_.river_access, offset);
    shift(anthroposphere_.harbor_score, offset);
    shift(anthroposphere_.owner_settlement_id, offset);
    shift(anthroposphere_.owner_polity_id, offset);
    shift(geosphere_.mountainheights, offset);
    shift(geosphere_.mountainridges, offset);
    shift(geosphere_.craterrims, offset);
    shift(geosphere_.cratercentres, offset);
    shift(geosphere_.extraelevmap, offset);
    shift(geosphere_.oceanridgemap, offset);
    shift(geosphere_.oceanridgeheightmap, offset);
    shift(geosphere_.oceanriftmap, offset);
    shift(geosphere_.oceanridgeanglemap, offset);
    shift(geosphere_.mountainislandmap, offset);
    shift(presentation_.noshademap, offset);
    shift(geosphere_.volcanomap, offset);
    shift(geosphere_.stratomap, offset);
    shift(presentation_.testmap, offset);

    for (int i = 0; i < geosphere_.itscraterno; i++)
    {
        geosphere_.cratercentreslist[i].x = geosphere_.cratercentreslist[i].x - offset;

        if (geosphere_.cratercentreslist[i].x < 0)
            geosphere_.cratercentreslist[i].x = geosphere_.cratercentreslist[i].x + body_.itswidth;
    }

    for (Settlement& settlement : anthroposphere_.settlementlist)
    {
        settlement.x -= offset;

        while (settlement.x < 0)
            settlement.x += body_.itswidth;

        while (settlement.x > body_.itswidth)
            settlement.x -= body_.itswidth;
    }

    const float worldspan = static_cast<float>(body_.itswidth + 1);
    const float offsetf = static_cast<float>(offset);
    const auto shiftcentroidx = [worldspan, offsetf](float x) -> float
    {
        float shifted = fmodf(x - offsetf, worldspan);

        if (shifted < 0.0f)
            shifted += worldspan;

        return shifted;
    };

    for (TectonicBoundarySegment& segment : geosphere_.tectonicboundarysegmentlist)
        segment.centroidX = shiftcentroidx(segment.centroidX);

    for (TectonicDeformingRegion& region : geosphere_.tectonicdeformingregionlist)
        region.centroidX = shiftcentroidx(region.centroidX);
}

void planet::setmaxriverflow()
{
    int largest = 0;
    int current = 0;

    for (int i = 0; i <= body_.itswidth; i++)
    {
        for (int j = 0; j <= body_.itsheight; j++)
        {
            current = riveraveflow(i, j);

            if (current > largest)
                largest = current;
        }
    }
    hydrosphere_.itsmaxriverflow = largest;
}

void planet::saveworld(string filename)
{
#ifdef ENABLE_PROFILER
    highres_timer_t timer("Save World"); // 26.5s => 10.9s
#endif
    ofstream outfile;
    outfile.open(filename, ios::out);

    writevariable(outfile, itssaveversion);
    writevariable(outfile, body_.itssize);
    writevariable(outfile, body_.itswidth);
    writevariable(outfile, body_.itsheight);
    writevariable(outfile, body_.itsseed);
    writevariable(outfile, body_.itsrotation);
    writevariable(outfile, body_.itstilt);
    writevariable(outfile, body_.itseccentricity);
    writevariable(outfile, body_.itsperihelion);
    writevariable(outfile, body_.itsgravity);
    writevariable(outfile, body_.itslunar);
    writevariable(outfile, atmosphere_.itstempdecrease);
    writevariable(outfile, atmosphere_.itsnorthpolaradjust);
    writevariable(outfile, atmosphere_.itssouthpolaradjust);
    writevariable(outfile, atmosphere_.itsaveragetemp);
    writevariable(outfile, climate_.itsnorthpolartemp);
    writevariable(outfile, climate_.itssouthpolartemp);
    writevariable(outfile, climate_.itseqtemp);
    writevariable(outfile, atmosphere_.itswaterpickup);
    writevariable(outfile, hydrosphere_.itsriverfactor);
    writevariable(outfile, hydrosphere_.itsriverlandreduce);
    writevariable(outfile, hydrosphere_.itsestuarylimit);
    writevariable(outfile, cryosphere_.itsglacialtemp);
    writevariable(outfile, cryosphere_.itsglaciertemp);
    writevariable(outfile, geosphere_.itsmountainreduce);
    writevariable(outfile, atmosphere_.itsclimateno);
    writevariable(outfile, geosphere_.itsmaxheight);
    writevariable(outfile, geosphere_.itssealevel);
    writevariable(outfile, geosphere_.itslandtotal);
    writevariable(outfile, geosphere_.itsseatotal);
    writevariable(outfile, geosphere_.itscraterno);

    writevariable(outfile, presentation_.itslandshading);
    writevariable(outfile, presentation_.itslakeshading);
    writevariable(outfile, presentation_.itsseashading);
    writevariable(outfile, presentation_.itsshadingdir);
    writevariable(outfile, presentation_.itssnowchange);
    writevariable(outfile, presentation_.itsseaiceappearance);
    writevariable(outfile, presentation_.itslandmarbling);
    writevariable(outfile, presentation_.itslakemarbling);
    writevariable(outfile, presentation_.itsseamarbling);
    writevariable(outfile, presentation_.itsminriverflowglobal);
    writevariable(outfile, presentation_.itsminriverflowregional);
    writevariable(outfile, presentation_.itsmangroves);
    writevariable(outfile, presentation_.itscolourcliffs);
    writevariable(outfile, presentation_.itsseaice1);
    writevariable(outfile, presentation_.itsseaice2);
    writevariable(outfile, presentation_.itsseaice3);
    writevariable(outfile, presentation_.itsocean1);
    writevariable(outfile, presentation_.itsocean2);
    writevariable(outfile, presentation_.itsocean3);
    writevariable(outfile, presentation_.itsdeepocean1);
    writevariable(outfile, presentation_.itsdeepocean2);
    writevariable(outfile, presentation_.itsdeepocean3);
    writevariable(outfile, presentation_.itsbase1);
    writevariable(outfile, presentation_.itsbase2);
    writevariable(outfile, presentation_.itsbase3);
    writevariable(outfile, presentation_.itsbasetemp1);
    writevariable(outfile, presentation_.itsbasetemp2);
    writevariable(outfile, presentation_.itsbasetemp3);
    writevariable(outfile, presentation_.itshighbase1);
    writevariable(outfile, presentation_.itshighbase2);
    writevariable(outfile, presentation_.itshighbase3);
    writevariable(outfile, presentation_.itsdesert1);
    writevariable(outfile, presentation_.itsdesert2);
    writevariable(outfile, presentation_.itsdesert3);
    writevariable(outfile, presentation_.itshighdesert1);
    writevariable(outfile, presentation_.itshighdesert2);
    writevariable(outfile, presentation_.itshighdesert3);
    writevariable(outfile, presentation_.itscolddesert1);
    writevariable(outfile, presentation_.itscolddesert2);
    writevariable(outfile, presentation_.itscolddesert3);
    writevariable(outfile, presentation_.itsgrass1);
    writevariable(outfile, presentation_.itsgrass2);
    writevariable(outfile, presentation_.itsgrass3);
    writevariable(outfile, presentation_.itscold1);
    writevariable(outfile, presentation_.itscold2);
    writevariable(outfile, presentation_.itscold3);
    writevariable(outfile, presentation_.itstundra1);
    writevariable(outfile, presentation_.itstundra2);
    writevariable(outfile, presentation_.itstundra3);
    writevariable(outfile, presentation_.itseqtundra1);
    writevariable(outfile, presentation_.itseqtundra2);
    writevariable(outfile, presentation_.itseqtundra3);
    writevariable(outfile, presentation_.itssaltpan1);
    writevariable(outfile, presentation_.itssaltpan2);
    writevariable(outfile, presentation_.itssaltpan3);
    writevariable(outfile, presentation_.itserg1);
    writevariable(outfile, presentation_.itserg2);
    writevariable(outfile, presentation_.itserg3);
    writevariable(outfile, presentation_.itswetlands1);
    writevariable(outfile, presentation_.itswetlands2);
    writevariable(outfile, presentation_.itswetlands3);
    writevariable(outfile, presentation_.itslake1);
    writevariable(outfile, presentation_.itslake2);
    writevariable(outfile, presentation_.itslake3);
    writevariable(outfile, presentation_.itsriver1);
    writevariable(outfile, presentation_.itsriver2);
    writevariable(outfile, presentation_.itsriver3);
    writevariable(outfile, presentation_.itsglacier1);
    writevariable(outfile, presentation_.itsglacier2);
    writevariable(outfile, presentation_.itsglacier3);
    writevariable(outfile, presentation_.itssand1);
    writevariable(outfile, presentation_.itssand2);
    writevariable(outfile, presentation_.itssand3);
    writevariable(outfile, presentation_.itsmud1);
    writevariable(outfile, presentation_.itsmud2);
    writevariable(outfile, presentation_.itsmud3);
    writevariable(outfile, presentation_.itsshingle1);
    writevariable(outfile, presentation_.itsshingle2);
    writevariable(outfile, presentation_.itsshingle3);
    writevariable(outfile, presentation_.itsmangrove1);
    writevariable(outfile, presentation_.itsmangrove2);
    writevariable(outfile, presentation_.itsmangrove3);
    writevariable(outfile, presentation_.itshighlight1);
    writevariable(outfile, presentation_.itshighlight2);
    writevariable(outfile, presentation_.itshighlight3);
    writevariable(outfile, presentation_.itsshowmapoutline);
    writevariable(outfile, presentation_.itsoutline1);
    writevariable(outfile, presentation_.itsoutline2);
    writevariable(outfile, presentation_.itsoutline3);
    writevariable(outfile, presentation_.itselevationlow1);
    writevariable(outfile, presentation_.itselevationlow2);
    writevariable(outfile, presentation_.itselevationlow3);
    writevariable(outfile, presentation_.itselevationhigh1);
    writevariable(outfile, presentation_.itselevationhigh2);
    writevariable(outfile, presentation_.itselevationhigh3);
    writevariable(outfile, presentation_.itstemperaturecold1);
    writevariable(outfile, presentation_.itstemperaturecold2);
    writevariable(outfile, presentation_.itstemperaturecold3);
    writevariable(outfile, presentation_.itstemperaturetemperate1);
    writevariable(outfile, presentation_.itstemperaturetemperate2);
    writevariable(outfile, presentation_.itstemperaturetemperate3);
    writevariable(outfile, presentation_.itstemperaturehot1);
    writevariable(outfile, presentation_.itstemperaturehot2);
    writevariable(outfile, presentation_.itstemperaturehot3);
    writevariable(outfile, presentation_.itsprecipitationdry1);
    writevariable(outfile, presentation_.itsprecipitationdry2);
    writevariable(outfile, presentation_.itsprecipitationdry3);
    writevariable(outfile, presentation_.itsprecipitationwet1);
    writevariable(outfile, presentation_.itsprecipitationwet2);
    writevariable(outfile, presentation_.itsprecipitationwet3);
    for (int i = 0; i < CLIMATEMAPSEACOLOURCOUNT; i++)
    {
        for (int j = 0; j < 3; j++)
            writevariable(outfile, presentation_.itsclimatemapseacolours[i][j]);
    }
    for (int i = 0; i < CLIMATEMAPCOLOURCOUNT; i++)
    {
        for (int j = 0; j < 3; j++)
            writevariable(outfile, presentation_.itsclimatemapcolours[i][j]);
    }
    for (int i = 0; i < RIVERMAPCOLOURCOUNT; i++)
    {
        for (int j = 0; j < 3; j++)
            writevariable(outfile, presentation_.itsrivermapcolours[i][j]);
    }
    for (int i = 0; i < RIVERMAPFEATURECOUNT; i++)
        writevariable(outfile, presentation_.itsshowrivermapfeatures[i]);

    for (int gradient = 0; gradient < MAPGRADIENTTYPECOUNT; gradient++)
    {
        writevariable(outfile, presentation_.itsmapgradientstopcounts[gradient]);
        writevariable(outfile, presentation_.itsmapgradientdiscrete[gradient]);

        for (int stop = 0; stop < MAPGRADIENTMAXSTOPS; stop++)
        {
            writevariable(outfile, presentation_.itsmapgradientpositions[gradient][stop]);

            for (int channel = 0; channel < 3; channel++)
                writevariable(outfile, presentation_.itsmapgradientcolours[gradient][stop][channel]);
        }
    }

    for (int i = 0; i < BIOMEMAPCOLOURCOUNT; i++)
    {
        for (int j = 0; j < 3; j++)
            writevariable(outfile, presentation_.itsbiomemapcolours[i][j]);
    }

    writelegacydata(outfile, climate_.jantempmap);
    writelegacydata(outfile, climate_.jultempmap);
    writedata(outfile, climate_.climatemap);
    writedata(outfile, climate_.biomemap);
    writelegacydata(outfile, climate_.janrainmap);
    writelegacydata(outfile, climate_.julrainmap);
    writelegacydata(outfile, climate_.janmountainrainmap);
    writelegacydata(outfile, climate_.julmountainrainmap);
    writedata(outfile, climate_.janmountainraindirmap);
    writedata(outfile, climate_.julmountainraindirmap);
    writedata(outfile, climate_.seaicemap);
    writedata(outfile, hydrosphere_.rivermapdir);
    writedata(outfile, hydrosphere_.rivermapjan);
    writedata(outfile, hydrosphere_.rivermapjul);
    writedata(outfile, climate_.windmap);
    writedata(outfile, hydrosphere_.lakemap);
    writedata(outfile, geosphere_.roughnessmap);
    writedata(outfile, geosphere_.mountainridges);
    writedata(outfile, geosphere_.mountainheights);
    writedata(outfile, geosphere_.craterrims);
    writedata(outfile, geosphere_.cratercentres);
    writedata(outfile, geosphere_.mapnom);
    writedata(outfile, hydrosphere_.tidalmap);
    writedata(outfile, hydrosphere_.riftlakemapsurface);
    writedata(outfile, hydrosphere_.riftlakemapbed);
    writedata(outfile, hydrosphere_.lakestartmap);
    writedata(outfile, hydrosphere_.specials);
    writedata(outfile, geosphere_.geologicregimemap);
    writedata(outfile, geosphere_.tectonicconvergencemap);
    writedata(outfile, geosphere_.tectonicdivergencemap);
    writedata(outfile, geosphere_.tectonicshearmap);
    writedata(outfile, hydrosphere_.basinclassmap);
    writedata(outfile, geosphere_.erosionpotentialmap);
    writedata(outfile, geosphere_.depositionpotentialmap);
    writedata(outfile, geosphere_.floodplainfertilitymap);
    writedata(outfile, geosphere_.metalorereservemap);
    writedata(outfile, geosphere_.placerreservemap);
    writedata(outfile, geosphere_.evaporitereservemap);
    writedata(outfile, geosphere_.volcanicreservemap);
    writedata(outfile, biosphere_.fisheryreservemap);
    writedata(outfile, geosphere_.extraelevmap);
    writedata(outfile, hydrosphere_.deltamapdir);
    writedata(outfile, hydrosphere_.deltamapjan);
    writedata(outfile, hydrosphere_.deltamapjul);
    writedata(outfile, geosphere_.islandmap);
    writedata(outfile, geosphere_.mountainislandmap);
    writedata(outfile, geosphere_.oceanridgemap);
    writedata(outfile, geosphere_.oceanridgeheightmap);
    writedata(outfile, geosphere_.oceanriftmap);
    writedata(outfile, geosphere_.oceanridgeoffsetmap);
    writedata(outfile, geosphere_.oceanridgeanglemap);
    writedata(outfile, geosphere_.volcanomap);
    writedata(outfile, geosphere_.stratomap);
    writedata(outfile, presentation_.noshademap);
    writedata(outfile, geosphere_.noisemap);
    writedata(outfile, presentation_.testmap);
    writedata(outfile, anthroposphere_.settlement_suitability);
    writedata(outfile, anthroposphere_.social_infrastructure);
    writedata(outfile, anthroposphere_.agricultural_capacity);
    writedata(outfile, anthroposphere_.route_traffic);
    writedata(outfile, anthroposphere_.river_access);
    writedata(outfile, anthroposphere_.harbor_score);
    writedata(outfile, anthroposphere_.owner_settlement_id);
    writedata(outfile, anthroposphere_.owner_polity_id);

    const auto writestring = [&](const std::string& value)
    {
        const int length = static_cast<int>(value.size());
        writevariable(outfile, length);

        if (length > 0)
            outfile.write(value.data(), static_cast<std::streamsize>(length));
    };

    const int settlementcount = static_cast<int>(anthroposphere_.settlementlist.size());
    writevariable(outfile, settlementcount);

    for (const Settlement& settlement : anthroposphere_.settlementlist)
    {
        writevariable(outfile, settlement.id);
        writevariable(outfile, settlement.x);
        writevariable(outfile, settlement.y);
        writestring(settlement.name);
        writevariable(outfile, settlement.urbanPopulation);
        writevariable(outfile, settlement.ruralPopulation);
        writevariable(outfile, settlement.carryingCapacity);
        writevariable(outfile, settlement.infrastructure);
        writevariable(outfile, settlement.marketStrength);
        writevariable(outfile, settlement.harbor);
        writevariable(outfile, settlement.riverAccess);
        writevariable(outfile, settlement.polityId);
    }

    const int politycount = static_cast<int>(anthroposphere_.politylist.size());
    writevariable(outfile, politycount);

    for (const Polity& polity : anthroposphere_.politylist)
    {
        writevariable(outfile, polity.id);
        writestring(polity.name);
        writevariable(outfile, polity.capitalSettlementId);
        writevariable(outfile, polity.cohesion);
        writevariable(outfile, polity.population);
        writevariable(outfile, polity.infrastructure);
        writevariable(outfile, polity.wealth);
        writevariable(outfile, polity.militaryPressure);
    }

    const int routecount = static_cast<int>(anthroposphere_.routeedgelist.size());
    writevariable(outfile, routecount);

    for (const RouteEdge& route : anthroposphere_.routeedgelist)
    {
        writevariable(outfile, route.fromSettlementId);
        writevariable(outfile, route.toSettlementId);
        const std::uint8_t mode = static_cast<std::uint8_t>(route.mode);
        writevariable(outfile, mode);
        writevariable(outfile, route.cost);
        writevariable(outfile, route.capacity);
        writevariable(outfile, route.traffic);
    }

    const int eventcount = static_cast<int>(anthroposphere_.historyeventlist.size());
    writevariable(outfile, eventcount);

    for (const HistoryEvent& event : anthroposphere_.historyeventlist)
    {
        writevariable(outfile, event.year);
        writestring(event.type);
        writevariable(outfile, event.primarySettlementId);
        writevariable(outfile, event.primaryPolityId);
        writevariable(outfile, event.secondarySettlementId);
        writevariable(outfile, event.secondaryPolityId);
        writestring(event.summary);
    }

    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
    {
        writeshortvectordata(outfile, climateworld::legacyfield(climate_.seasonaltempmaps[season]));
        writeshortvectordata(outfile, climateworld::legacyfield(climate_.seasonalrainmaps[season]));
        writeshortvectordata(outfile, climateworld::legacyfield(climate_.seasonalpressuremaps[season]));
        writeshortvectordata(outfile, climateworld::legacyfield(climate_.seasonaluwindmaps[season]));
        writeshortvectordata(outfile, climateworld::legacyfield(climate_.seasonalvwindmaps[season]));
        writeshortvectordata(outfile, climateworld::legacyfield(climate_.seasonalupperheightmaps[season]));
        writeshortvectordata(outfile, climateworld::legacyfield(climate_.seasonalupperuwindmaps[season]));
        writeshortvectordata(outfile, climateworld::legacyfield(climate_.seasonaluppervwindmaps[season]));
        writeshortvectordata(outfile, climateworld::legacyfield(climate_.seasonalverticalvelocitymaps[season], climatelegacyscale::verticalvelocity));
        writeshortvectordata(outfile, climateworld::legacyfield(climate_.seasonalcurrentumaps[season], climatelegacyscale::current));
        writeshortvectordata(outfile, climateworld::legacyfield(climate_.seasonalcurrentvmaps[season], climatelegacyscale::current));
        writeshortvectordata(outfile, climateworld::legacyfield(climate_.seasonalsstmaps[season]));
        writeshortvectordata(outfile, climateworld::legacyfield(climate_.seasonalevaporationmaps[season]));
        writeshortvectordata(outfile, climateworld::legacyfield(climate_.seasonalmaritimeinfluencemaps[season], climatelegacyscale::maritimeinfluence));
        writeshortvectordata(outfile, climateworld::legacyfield(climate_.seasonalmaritimethermalanomalymaps[season], climatelegacyscale::maritimethermalanomaly));
        writeshortvectordata(outfile, climate_.seasonalmaritimefetchmaps[season]);
        writeshortvectordata(outfile, climateworld::legacyfield(climate_.seasonalmoisturemaps[season]));
        writeshortvectordata(outfile, climateworld::legacyfield(climate_.seasonalconvergencemaps[season], climatelegacyscale::convergence));
        writeshortvectordata(outfile, climateworld::legacyfield(climate_.seasonalupliftmaps[season], climatelegacyscale::orographicmotion));
        writeshortvectordata(outfile, climateworld::legacyfield(climate_.seasonalsubsidencemaps[season], climatelegacyscale::orographicmotion));
    }

    for (int i = 0; i < ARRAYWIDTH; i++)
    {
        for (int j = 0; j < 6; j++)
            writevariable(outfile, climate_.horselats[i][j]);
    }

    for (int i = 0; i < geosphere_.itscraterno; i++)
    {
        writevariable(outfile, geosphere_.cratercentreslist[i].w);
        writevariable(outfile, geosphere_.cratercentreslist[i].x);
        writevariable(outfile, geosphere_.cratercentreslist[i].y);
        writevariable(outfile, geosphere_.cratercentreslist[i].z);
    }

    writevariable(outfile, geosphere_.itstectonictimeoriginstep);
    writevariable(outfile, geosphere_.itstectonictimemyr);
    writevariable(outfile, geosphere_.itstectonicdeltatimemyr);
    writedata(outfile, geosphere_.tectoniccrustagemyrmap);
    writedata(outfile, geosphere_.tectoniccrustthicknessmap);
    writedata(outfile, geosphere_.tectoniccrustclassmap);
    writedata(outfile, geosphere_.tectonicuplifttendencymap);
    writedata(outfile, geosphere_.tectonicsubsidencetendencymap);
    writedata(outfile, geosphere_.tectonicaccumulatedstrainmap);
    writedata(outfile, geosphere_.tectonicboundarytypemap);
    writedata(outfile, geosphere_.tectonicboundarydistancemap);
    writedata(outfile, geosphere_.tectonicboundarysegmentidmap);
    writedata(outfile, geosphere_.tectonicnearestboundaryidmap);
    writedata(outfile, geosphere_.tectonicboundaryhistorymap);
    writedata(outfile, geosphere_.tectonicdeformingregionidmap);
    writedata(outfile, geosphere_.tectonicdeformingregiontypemap);
    writedata(outfile, geosphere_.tectonicdeformationratemap);
    writedata(outfile, geosphere_.tectonicdeformationvelocityxmap);
    writedata(outfile, geosphere_.tectonicdeformationvelocityymap);
    writevariable(outfile, geosphere_.itstectonicplatecount);
    writevariable(outfile, geosphere_.itstectonicsealevelm);

    int boundarysegmentcount = static_cast<int>(geosphere_.tectonicboundarysegmentlist.size());
    writevariable(outfile, boundarysegmentcount);
    for (const TectonicBoundarySegment& segment : geosphere_.tectonicboundarysegmentlist)
    {
        writevariable(outfile, segment.id);
        writevariable(outfile, segment.leftPlateId);
        writevariable(outfile, segment.rightPlateId);
        writevariable(outfile, segment.cellCount);
        writevariable(outfile, segment.persistenceSteps);
        writevariable(outfile, segment.centroidX);
        writevariable(outfile, segment.centroidY);
        writevariable(outfile, segment.lengthCells);
        writevariable(outfile, segment.averageNormalMotion);
        writevariable(outfile, segment.averageShearMotion);
        writevariable(outfile, segment.averageConvergenceScore);
        writevariable(outfile, segment.averageDivergenceScore);
        writevariable(outfile, segment.averageShearScore);
        std::uint8_t boundarytype = static_cast<std::uint8_t>(segment.boundaryType);
        std::uint8_t geologicregime = static_cast<std::uint8_t>(segment.geologicRegime);
        writevariable(outfile, boundarytype);
        writevariable(outfile, geologicregime);
        writevariable(outfile, segment.ageMyr);
    }

    int deformingregioncount = static_cast<int>(geosphere_.tectonicdeformingregionlist.size());
    writevariable(outfile, deformingregioncount);
    for (const TectonicDeformingRegion& region : geosphere_.tectonicdeformingregionlist)
    {
        writevariable(outfile, region.id);
        writevariable(outfile, region.boundarySegmentId);
        writevariable(outfile, region.primaryPlateId);
        writevariable(outfile, region.secondaryPlateId);
        writevariable(outfile, region.cellCount);
        writevariable(outfile, region.persistenceSteps);
        writevariable(outfile, region.centroidX);
        writevariable(outfile, region.centroidY);
        writevariable(outfile, region.averageDeformationRate);
        writevariable(outfile, region.averageInterpolatedVelocityX);
        writevariable(outfile, region.averageInterpolatedVelocityY);
        writevariable(outfile, region.averageNormalMotion);
        writevariable(outfile, region.averageShearMotion);
        std::uint8_t regiontype = static_cast<std::uint8_t>(region.type);
        writevariable(outfile, regiontype);
        writevariable(outfile, region.ageMyr);
    }

    writevariable(outfile, static_cast<int>(climate_.weatheranomalystatebytes.size()));
    for (std::uint8_t byte : climate_.weatheranomalystatebytes)
        writevariable(outfile, static_cast<int>(byte));

    writemineraldeposits(outfile);

    if (!outfile.good())
    {
        cerr << "Error writing world '" << filename << "'" << endl;
    }
}

// Optional trailing section: loaders that predate it stop before the tag, and
// saves without it load with no deposits.
static const char* const kMineralDepositSection = "mineral-deposits-v1";

void planet::writemineraldeposits(ofstream& outfile)
{
    using namespace mineralsystems;
    const std::streamsize precision = outfile.precision(12);
    outfile << kMineralDepositSection << '\n';
    writevariable(outfile, static_cast<int>(geosphere_.mineraldepositlist.size()));
    for (const Deposit& deposit : geosphere_.mineraldepositlist)
    {
        const FormationRecord& formation = deposit.formation;
        writevariable(outfile, static_cast<long>(deposit.id));
        writevariable(outfile, static_cast<int>(deposit.system));
        writevariable(outfile, static_cast<int>(deposit.origin));
        writevariable(outfile, static_cast<int>(deposit.fate));
        writevariable(outfile, deposit.fateTimeMyr);
        writevariable(outfile, deposit.u);
        writevariable(outfile, deposit.v);
        writevariable(outfile, formation.timeMyr);
        writevariable(outfile, static_cast<int>(formation.geologicRegime));
        writevariable(outfile, static_cast<int>(formation.boundaryType));
        writevariable(outfile, static_cast<int>(formation.hostCrustClass));
        writevariable(outfile, static_cast<long>(formation.boundarySegmentId));
        writevariable(outfile, static_cast<long>(formation.provinceId));
        writevariable(outfile, formation.boundaryDistanceKm);
        writevariable(outfile, formation.boundaryAgeMyr);
        writevariable(outfile, formation.settingIntensity);
        writevariable(outfile, deposit.oreMt);
        writevariable(outfile, deposit.goldGradeGt);
        writevariable(outfile, deposit.copperGradePct);
        writevariable(outfile, deposit.formationDepthM);
        writevariable(outfile, deposit.verticalExtentM);
        writevariable(outfile, deposit.exhumationM);
        writevariable(outfile, deposit.erodedFraction);
        writevariable(outfile, deposit.recentLiberatedGoldT);
        writevariable(outfile, static_cast<int>(deposit.sources.size()));
        for (const SourceShare& share : deposit.sources)
        {
            writevariable(outfile, static_cast<long>(share.depositId));
            writevariable(outfile, share.goldTonnes);
        }
        const SurfaceState& surface = deposit.surface;
        writevariable(outfile, static_cast<int>(surface.evaluated));
        writevariable(outfile, surface.worldX);
        writevariable(outfile, surface.worldY);
        writevariable(outfile, static_cast<int>(surface.submerged));
        writevariable(outfile, static_cast<int>(surface.iceCovered));
        writevariable(outfile, surface.elevationM);
        writevariable(outfile, surface.localReliefM);
        writevariable(outfile, surface.annualRainfallMm);
    }

    const MineralBudget& budget = geosphere_.mineralbudgetstate;
    for (int system = 0; system < kDepositSystemCount; ++system)
    {
        writevariable(outfile, budget.formedGoldT[system]);
        writevariable(outfile, budget.formedCount[system]);
    }
    for (double value : { budget.inheritedGoldT, budget.erodedGoldT, budget.subductedGoldT, budget.placerInjectedGoldT,
             budget.placerDepositGoldT, budget.placerDispersedGoldT, budget.placerLakeGoldT, budget.placerExportedGoldT })
        writevariable(outfile, value);
    for (int value : { budget.inheritedCount, budget.erodedCount, budget.subductedCount })
        writevariable(outfile, value);
    outfile.precision(precision);
}

bool planet::readmineraldeposits(ifstream& infile)
{
    using namespace mineralsystems;
    geosphere_.mineraldepositlist.clear();
    geosphere_.mineralbudgetstate = {};

    string tag;
    if (!getline(infile, tag) || tag != kMineralDepositSection)
    {
        infile.clear(); // Worlds saved before the section end here.
        return true;
    }

    const auto readint = [&](int minimum, int maximum, int& value)
    {
        readvariable(infile, value);
        return value >= minimum && value <= maximum;
    };
    const auto readid = [&](std::uint32_t& value)
    {
        long read = 0;
        readvariable(infile, read);
        value = static_cast<std::uint32_t>(std::max(0L, read));
    };

    int count = 0;
    if (!readint(0, 50000000, count))
        return false;
    geosphere_.mineraldepositlist.resize(static_cast<size_t>(count));
    for (Deposit& deposit : geosphere_.mineraldepositlist)
    {
        FormationRecord& formation = deposit.formation;
        int system = 0, origin = 0, fate = 0, regime = 0, boundary = 0, crust = 0;
        readid(deposit.id);
        if (!readint(0, kDepositSystemCount - 1, system) || !readint(0, 2, origin) || !readint(0, 2, fate))
            return false;
        deposit.system = static_cast<DepositSystem>(system);
        deposit.origin = static_cast<DepositOrigin>(origin);
        deposit.fate = static_cast<DepositFate>(fate);
        readvariable(infile, deposit.fateTimeMyr);
        readvariable(infile, deposit.u);
        readvariable(infile, deposit.v);
        readvariable(infile, formation.timeMyr);
        if (!readint(0, 255, regime) || !readint(0, 255, boundary) || !readint(0, 255, crust))
            return false;
        formation.geologicRegime = static_cast<std::uint8_t>(regime);
        formation.boundaryType = static_cast<std::uint8_t>(boundary);
        formation.hostCrustClass = static_cast<std::uint8_t>(crust);
        readid(formation.boundarySegmentId);
        readid(formation.provinceId);
        readvariable(infile, formation.boundaryDistanceKm);
        readvariable(infile, formation.boundaryAgeMyr);
        readvariable(infile, formation.settingIntensity);
        readvariable(infile, deposit.oreMt);
        readvariable(infile, deposit.goldGradeGt);
        readvariable(infile, deposit.copperGradePct);
        readvariable(infile, deposit.formationDepthM);
        readvariable(infile, deposit.verticalExtentM);
        readvariable(infile, deposit.exhumationM);
        readvariable(infile, deposit.erodedFraction);
        readvariable(infile, deposit.recentLiberatedGoldT);
        int sources = 0;
        if (!readint(0, 1000, sources))
            return false;
        deposit.sources.resize(static_cast<size_t>(sources));
        for (SourceShare& share : deposit.sources)
        {
            readid(share.depositId);
            readvariable(infile, share.goldTonnes);
        }
        SurfaceState& surface = deposit.surface;
        int evaluated = 0, submerged = 0, ice = 0;
        readvariable(infile, evaluated);
        readvariable(infile, surface.worldX);
        readvariable(infile, surface.worldY);
        readvariable(infile, submerged);
        readvariable(infile, ice);
        surface.evaluated = evaluated != 0;
        surface.submerged = submerged != 0;
        surface.iceCovered = ice != 0;
        readvariable(infile, surface.elevationM);
        readvariable(infile, surface.localReliefM);
        readvariable(infile, surface.annualRainfallMm);
    }

    MineralBudget& budget = geosphere_.mineralbudgetstate;
    for (int system = 0; system < kDepositSystemCount; ++system)
    {
        readvariable(infile, budget.formedGoldT[system]);
        readvariable(infile, budget.formedCount[system]);
    }
    for (double* value : { &budget.inheritedGoldT, &budget.erodedGoldT, &budget.subductedGoldT, &budget.placerInjectedGoldT,
             &budget.placerDepositGoldT, &budget.placerDispersedGoldT, &budget.placerLakeGoldT, &budget.placerExportedGoldT })
        readvariable(infile, *value);
    for (int* value : { &budget.inheritedCount, &budget.erodedCount, &budget.subductedCount })
        readvariable(infile, *value);
    return true;
}

bool planet::loadworld(string filename)
{
#ifdef ENABLE_PROFILER
    highres_timer_t timer("Load World"); // 9.1s => 8.8s
#endif
    ifstream infile;
    infile.open(filename, ios::in);

    int val;
    readvariable(infile, val);

    const int fileversion = val;

    if (fileversion != itssaveversion && fileversion != 17) // Incompatible file format!
        return 0;

    readvariable(infile, body_.itssize);
    readvariable(infile, body_.itswidth);
    readvariable(infile, body_.itsheight);
    const int columns = body_.itswidth + 1;
    const int rows = body_.itsheight + 1;
    if (columns < 2 || columns > ARRAYWIDTH || columns % 2 != 0 ||
        rows != columns / 2)
    {
        return false;
    }
    readvariable(infile, body_.itsseed);
    resizegrids();
    readvariable(infile, body_.itsrotation);
    readvariable(infile, body_.itstilt);
    readvariable(infile, body_.itseccentricity);
    readvariable(infile, body_.itsperihelion);
    readvariable(infile, body_.itsgravity);
    readvariable(infile, body_.itslunar);
    readvariable(infile, atmosphere_.itstempdecrease);
    readvariable(infile, atmosphere_.itsnorthpolaradjust);
    readvariable(infile, atmosphere_.itssouthpolaradjust);
    readvariable(infile, atmosphere_.itsaveragetemp);
    readvariable(infile, climate_.itsnorthpolartemp);
    readvariable(infile, climate_.itssouthpolartemp);
    readvariable(infile, climate_.itseqtemp);
    readvariable(infile, atmosphere_.itswaterpickup);
    readvariable(infile, hydrosphere_.itsriverfactor);
    readvariable(infile, hydrosphere_.itsriverlandreduce);
    readvariable(infile, hydrosphere_.itsestuarylimit);
    readvariable(infile, cryosphere_.itsglacialtemp);
    readvariable(infile, cryosphere_.itsglaciertemp);
    readvariable(infile, geosphere_.itsmountainreduce);
    readvariable(infile, atmosphere_.itsclimateno);
    readvariable(infile, geosphere_.itsmaxheight);
    readvariable(infile, geosphere_.itssealevel);
    readvariable(infile, geosphere_.itslandtotal);
    readvariable(infile, geosphere_.itsseatotal);
    readvariable(infile, geosphere_.itscraterno);

    readvariable(infile, presentation_.itslandshading);
    readvariable(infile, presentation_.itslakeshading);
    readvariable(infile, presentation_.itsseashading);
    readvariable(infile, presentation_.itsshadingdir);
    readvariable(infile, presentation_.itssnowchange);
    readvariable(infile, presentation_.itsseaiceappearance);
    readvariable(infile, presentation_.itslandmarbling);
    readvariable(infile, presentation_.itslakemarbling);
    readvariable(infile, presentation_.itsseamarbling);
    readvariable(infile, presentation_.itsminriverflowglobal);
    readvariable(infile, presentation_.itsminriverflowregional);
    readvariable(infile, presentation_.itsmangroves);
    readvariable(infile, presentation_.itscolourcliffs);
    readvariable(infile, presentation_.itsseaice1);
    readvariable(infile, presentation_.itsseaice2);
    readvariable(infile, presentation_.itsseaice3);
    readvariable(infile, presentation_.itsocean1);
    readvariable(infile, presentation_.itsocean2);
    readvariable(infile, presentation_.itsocean3);
    readvariable(infile, presentation_.itsdeepocean1);
    readvariable(infile, presentation_.itsdeepocean2);
    readvariable(infile, presentation_.itsdeepocean3);
    readvariable(infile, presentation_.itsbase1);
    readvariable(infile, presentation_.itsbase2);
    readvariable(infile, presentation_.itsbase3);
    readvariable(infile, presentation_.itsbasetemp1);
    readvariable(infile, presentation_.itsbasetemp2);
    readvariable(infile, presentation_.itsbasetemp3);
    readvariable(infile, presentation_.itshighbase1);
    readvariable(infile, presentation_.itshighbase2);
    readvariable(infile, presentation_.itshighbase3);
    readvariable(infile, presentation_.itsdesert1);
    readvariable(infile, presentation_.itsdesert2);
    readvariable(infile, presentation_.itsdesert3);
    readvariable(infile, presentation_.itshighdesert1);
    readvariable(infile, presentation_.itshighdesert2);
    readvariable(infile, presentation_.itshighdesert3);
    readvariable(infile, presentation_.itscolddesert1);
    readvariable(infile, presentation_.itscolddesert2);
    readvariable(infile, presentation_.itscolddesert3);
    readvariable(infile, presentation_.itsgrass1);
    readvariable(infile, presentation_.itsgrass2);
    readvariable(infile, presentation_.itsgrass3);
    readvariable(infile, presentation_.itscold1);
    readvariable(infile, presentation_.itscold2);
    readvariable(infile, presentation_.itscold3);
    readvariable(infile, presentation_.itstundra1);
    readvariable(infile, presentation_.itstundra2);
    readvariable(infile, presentation_.itstundra3);
    readvariable(infile, presentation_.itseqtundra1);
    readvariable(infile, presentation_.itseqtundra2);
    readvariable(infile, presentation_.itseqtundra3);
    readvariable(infile, presentation_.itssaltpan1);
    readvariable(infile, presentation_.itssaltpan2);
    readvariable(infile, presentation_.itssaltpan3);
    readvariable(infile, presentation_.itserg1);
    readvariable(infile, presentation_.itserg2);
    readvariable(infile, presentation_.itserg3);
    readvariable(infile, presentation_.itswetlands1);
    readvariable(infile, presentation_.itswetlands2);
    readvariable(infile, presentation_.itswetlands3);
    readvariable(infile, presentation_.itslake1);
    readvariable(infile, presentation_.itslake2);
    readvariable(infile, presentation_.itslake3);
    readvariable(infile, presentation_.itsriver1);
    readvariable(infile, presentation_.itsriver2);
    readvariable(infile, presentation_.itsriver3);
    readvariable(infile, presentation_.itsglacier1);
    readvariable(infile, presentation_.itsglacier2);
    readvariable(infile, presentation_.itsglacier3);
    readvariable(infile, presentation_.itssand1);
    readvariable(infile, presentation_.itssand2);
    readvariable(infile, presentation_.itssand3);
    readvariable(infile, presentation_.itsmud1);
    readvariable(infile, presentation_.itsmud2);
    readvariable(infile, presentation_.itsmud3);
    readvariable(infile, presentation_.itsshingle1);
    readvariable(infile, presentation_.itsshingle2);
    readvariable(infile, presentation_.itsshingle3);
    readvariable(infile, presentation_.itsmangrove1);
    readvariable(infile, presentation_.itsmangrove2);
    readvariable(infile, presentation_.itsmangrove3);
    readvariable(infile, presentation_.itshighlight1);
    readvariable(infile, presentation_.itshighlight2);
    readvariable(infile, presentation_.itshighlight3);

    if (fileversion >= 2)
    {
        readvariable(infile, presentation_.itsshowmapoutline);
        readvariable(infile, presentation_.itsoutline1);
        readvariable(infile, presentation_.itsoutline2);
        readvariable(infile, presentation_.itsoutline3);
        readvariable(infile, presentation_.itselevationlow1);
        readvariable(infile, presentation_.itselevationlow2);
        readvariable(infile, presentation_.itselevationlow3);
        readvariable(infile, presentation_.itselevationhigh1);
        readvariable(infile, presentation_.itselevationhigh2);
        readvariable(infile, presentation_.itselevationhigh3);
        readvariable(infile, presentation_.itstemperaturecold1);
        readvariable(infile, presentation_.itstemperaturecold2);
        readvariable(infile, presentation_.itstemperaturecold3);
        readvariable(infile, presentation_.itstemperaturetemperate1);
        readvariable(infile, presentation_.itstemperaturetemperate2);
        readvariable(infile, presentation_.itstemperaturetemperate3);
        readvariable(infile, presentation_.itstemperaturehot1);
        readvariable(infile, presentation_.itstemperaturehot2);
        readvariable(infile, presentation_.itstemperaturehot3);
        readvariable(infile, presentation_.itsprecipitationdry1);
        readvariable(infile, presentation_.itsprecipitationdry2);
        readvariable(infile, presentation_.itsprecipitationdry3);
        readvariable(infile, presentation_.itsprecipitationwet1);
        readvariable(infile, presentation_.itsprecipitationwet2);
        readvariable(infile, presentation_.itsprecipitationwet3);

        if (fileversion >= 3)
        {
            for (int i = 0; i < CLIMATEMAPSEACOLOURCOUNT; i++)
            {
                for (int j = 0; j < 3; j++)
                    readvariable(infile, presentation_.itsclimatemapseacolours[i][j]);
            }

            for (int i = 0; i < CLIMATEMAPCOLOURCOUNT; i++)
            {
                for (int j = 0; j < 3; j++)
                    readvariable(infile, presentation_.itsclimatemapcolours[i][j]);
            }

            for (int i = 0; i < RIVERMAPCOLOURCOUNT; i++)
            {
                for (int j = 0; j < 3; j++)
                    readvariable(infile, presentation_.itsrivermapcolours[i][j]);
            }

            for (int i = 0; i < RIVERMAPFEATURECOUNT; i++)
                readvariable(infile, presentation_.itsshowrivermapfeatures[i]);

            if (fileversion >= 4)
            {
                for (int gradient = 0; gradient < MAPGRADIENTTYPECOUNT; gradient++)
                {
                    readvariable(infile, presentation_.itsmapgradientstopcounts[gradient]);
                    readvariable(infile, presentation_.itsmapgradientdiscrete[gradient]);

                    for (int stop = 0; stop < MAPGRADIENTMAXSTOPS; stop++)
                    {
                        readvariable(infile, presentation_.itsmapgradientpositions[gradient][stop]);

                        for (int channel = 0; channel < 3; channel++)
                            readvariable(infile, presentation_.itsmapgradientcolours[gradient][stop][channel]);
                    }
                }
            }
            else
                initialisegradientmapappearance(*this);

            if (fileversion >= 5)
            {
                for (int i = 0; i < BIOMEMAPCOLOURCOUNT; i++)
                {
                    for (int j = 0; j < 3; j++)
                        readvariable(infile, presentation_.itsbiomemapcolours[i][j]);
                }
            }
            else
            {
                for (int i = 0; i < BIOMEMAPCOLOURCOUNT; i++)
                {
                    for (int j = 0; j < 3; j++)
                        presentation_.itsbiomemapcolours[i][j] = defaultbiomemapcolours[i][j];
                }
            }
        }
        else
        {
            static const int defaultclimatecolours[CLIMATEMAPCOLOURCOUNT][3] =
            {
                { 0, 0, 0 }, { 0, 0, 254 }, { 1, 119, 255 }, { 70, 169, 250 }, { 70, 169, 250 }, { 249, 15, 0 }, { 251, 150, 149 }, { 245, 163, 1 },
                { 254, 219, 99 }, { 255, 255, 0 }, { 198, 199, 1 }, { 184, 184, 114 }, { 138, 255, 162 }, { 86, 199, 112 }, { 30, 150, 66 }, { 192, 254, 109 },
                { 76, 255, 93 }, { 19, 203, 74 }, { 255, 8, 245 }, { 204, 3, 192 }, { 154, 51, 144 }, { 153, 100, 146 }, { 172, 178, 249 }, { 91, 121, 213 },
                { 78, 83, 175 }, { 54, 3, 130 }, { 0, 255, 245 }, { 32, 200, 250 }, { 0, 126, 125 }, { 0, 69, 92 }, { 178, 178, 178 }, { 104, 104, 104 }
            };
            static const int defaultclimateseacolours[CLIMATEMAPSEACOLOURCOUNT][3] =
            {
                { 13, 49, 109 }, { 228, 228, 255 }, { 243, 243, 255 }
            };

            for (int i = 0; i < CLIMATEMAPSEACOLOURCOUNT; i++)
            {
                for (int j = 0; j < 3; j++)
                    presentation_.itsclimatemapseacolours[i][j] = defaultclimateseacolours[i][j];
            }

            for (int i = 0; i < CLIMATEMAPCOLOURCOUNT; i++)
            {
                for (int j = 0; j < 3; j++)
                    presentation_.itsclimatemapcolours[i][j] = defaultclimatecolours[i][j];
            }

            presentation_.itsrivermapcolours[rivermapbackground] = { 255, 255, 255 };
            presentation_.itsrivermapcolours[rivermaplowflow] = { 255, 255, 255 };
            presentation_.itsrivermapcolours[rivermaphighflow] = { presentation_.itsriver1, presentation_.itsriver2, presentation_.itsriver3 };
            presentation_.itsrivermapcolours[rivermaplake] = { presentation_.itslake1, presentation_.itslake2, presentation_.itslake3 };
            presentation_.itsrivermapcolours[rivermapsaltpan] = { presentation_.itssaltpan1, presentation_.itssaltpan2, presentation_.itssaltpan3 };
            presentation_.itsrivermapcolours[rivermapwetlands] = { presentation_.itswetlands1, presentation_.itswetlands2, presentation_.itswetlands3 };
            presentation_.itsrivermapcolours[rivermapmud] = { presentation_.itsmud1, presentation_.itsmud2, presentation_.itsmud3 };
            presentation_.itsrivermapcolours[rivermapsand] = { presentation_.itssand1, presentation_.itssand2, presentation_.itssand3 };
            presentation_.itsrivermapcolours[rivermapshingle] = { presentation_.itsshingle1, presentation_.itsshingle2, presentation_.itsshingle3 };
            presentation_.itsrivermapcolours[rivermapvolcano] = { 240, 0, 0 };

            for (int i = 0; i < RIVERMAPFEATURECOUNT; i++)
                presentation_.itsshowrivermapfeatures[i] = true;

            initialisegradientmapappearance(*this);

            for (int i = 0; i < BIOMEMAPCOLOURCOUNT; i++)
            {
                for (int j = 0; j < 3; j++)
                    presentation_.itsbiomemapcolours[i][j] = defaultbiomemapcolours[i][j];
            }
        }
    }
    else
        setdefaultnonreliefmapappearance(*this);

    readlegacydata(infile, climate_.jantempmap);
    readlegacydata(infile, climate_.jultempmap);
    readdata(infile, climate_.climatemap);
    readdata(infile, climate_.biomemap);
    readlegacydata(infile, climate_.janrainmap);
    readlegacydata(infile, climate_.julrainmap);
    readlegacydata(infile, climate_.janmountainrainmap);
    readlegacydata(infile, climate_.julmountainrainmap);
    readdata(infile, climate_.janmountainraindirmap);
    readdata(infile, climate_.julmountainraindirmap);
    readdata(infile, climate_.seaicemap);
    readdata(infile, hydrosphere_.rivermapdir);
    readdata(infile, hydrosphere_.rivermapjan);
    readdata(infile, hydrosphere_.rivermapjul);
    readdata(infile, climate_.windmap);
    readdata(infile, hydrosphere_.lakemap);
    readdata(infile, geosphere_.roughnessmap);
    readdata(infile, geosphere_.mountainridges);
    readdata(infile, geosphere_.mountainheights);
    readdata(infile, geosphere_.craterrims);
    readdata(infile, geosphere_.cratercentres);
    readdata(infile, geosphere_.mapnom);
    readdata(infile, hydrosphere_.tidalmap);
    readdata(infile, hydrosphere_.riftlakemapsurface);
    readdata(infile, hydrosphere_.riftlakemapbed);
    readdata(infile, hydrosphere_.lakestartmap);
    readdata(infile, hydrosphere_.specials);
    if (fileversion >= 12)
    {
        readdata(infile, geosphere_.geologicregimemap);
        readdata(infile, geosphere_.tectonicconvergencemap);
        readdata(infile, geosphere_.tectonicdivergencemap);
        readdata(infile, geosphere_.tectonicshearmap);
        readdata(infile, hydrosphere_.basinclassmap);
        readdata(infile, geosphere_.erosionpotentialmap);
        readdata(infile, geosphere_.depositionpotentialmap);
        readdata(infile, geosphere_.floodplainfertilitymap);
        readdata(infile, geosphere_.metalorereservemap);
        readdata(infile, geosphere_.placerreservemap);
        readdata(infile, geosphere_.evaporitereservemap);
        readdata(infile, geosphere_.volcanicreservemap);
        readdata(infile, biosphere_.fisheryreservemap);
    }
    else
    {
        cleartectonicprovenanceinternal();

        parallelforrows(0, body_.itswidth, [&](int startx, int endx)
        {
            for (int i = startx; i <= endx; i++)
            {
                for (int j = 0; j <= body_.itsheight; j++)
                {
                    hydrosphere_.basinclassmap[i][j] = static_cast<std::uint8_t>(BasinClass::none);
                    geosphere_.erosionpotentialmap[i][j] = 0;
                    geosphere_.depositionpotentialmap[i][j] = 0;
                    geosphere_.floodplainfertilitymap[i][j] = 0;
                    geosphere_.metalorereservemap[i][j] = 0;
                    geosphere_.placerreservemap[i][j] = 0;
                    geosphere_.evaporitereservemap[i][j] = 0;
                    geosphere_.volcanicreservemap[i][j] = 0;
                    biosphere_.fisheryreservemap[i][j] = 0;
                }
            }
        }, 64);
    }
    readdata(infile, geosphere_.extraelevmap);
    readdata(infile, hydrosphere_.deltamapdir);
    readdata(infile, hydrosphere_.deltamapjan);
    readdata(infile, hydrosphere_.deltamapjul);
    readdata(infile, geosphere_.islandmap);
    readdata(infile, geosphere_.mountainislandmap);
    readdata(infile, geosphere_.oceanridgemap);
    readdata(infile, geosphere_.oceanridgeheightmap);
    readdata(infile, geosphere_.oceanriftmap);
    readdata(infile, geosphere_.oceanridgeoffsetmap);
    readdata(infile, geosphere_.oceanridgeanglemap);
    readdata(infile, geosphere_.volcanomap);
    readdata(infile, geosphere_.stratomap);
    readdata(infile, presentation_.noshademap);
    readdata(infile, geosphere_.noisemap);
    readdata(infile, presentation_.testmap);

    auto readstring = [&](std::string& value)
    {
        int length = 0;
        readvariable(infile, length);

        if (length <= 0)
        {
            value.clear();
            return;
        }

        value.assign(static_cast<size_t>(length), '\0');
        infile.read(value.data(), static_cast<std::streamsize>(length));
    };

    if (fileversion >= 14)
    {
        readdata(infile, anthroposphere_.settlement_suitability);
        readdata(infile, anthroposphere_.social_infrastructure);
        readdata(infile, anthroposphere_.agricultural_capacity);
        readdata(infile, anthroposphere_.route_traffic);
        readdata(infile, anthroposphere_.river_access);
        readdata(infile, anthroposphere_.harbor_score);
        readdata(infile, anthroposphere_.owner_settlement_id);
        readdata(infile, anthroposphere_.owner_polity_id);

        clearsocialstate();

        int settlementcount = 0;
        readvariable(infile, settlementcount);
        settlementcount = std::max(0, settlementcount);
        anthroposphere_.settlementlist.resize(static_cast<size_t>(settlementcount));

        for (Settlement& settlement : anthroposphere_.settlementlist)
        {
            readvariable(infile, settlement.id);
            readvariable(infile, settlement.x);
            readvariable(infile, settlement.y);
            readstring(settlement.name);
            readvariable(infile, settlement.urbanPopulation);
            readvariable(infile, settlement.ruralPopulation);
            readvariable(infile, settlement.carryingCapacity);
            readvariable(infile, settlement.infrastructure);
            readvariable(infile, settlement.marketStrength);
            readvariable(infile, settlement.harbor);
            readvariable(infile, settlement.riverAccess);
            readvariable(infile, settlement.polityId);
        }

        int politycount = 0;
        readvariable(infile, politycount);
        politycount = std::max(0, politycount);
        anthroposphere_.politylist.resize(static_cast<size_t>(politycount));

        for (Polity& polity : anthroposphere_.politylist)
        {
            readvariable(infile, polity.id);
            readstring(polity.name);
            readvariable(infile, polity.capitalSettlementId);
            readvariable(infile, polity.cohesion);
            readvariable(infile, polity.population);
            readvariable(infile, polity.infrastructure);
            readvariable(infile, polity.wealth);
            readvariable(infile, polity.militaryPressure);
        }

        int routecount = 0;
        readvariable(infile, routecount);
        routecount = std::max(0, routecount);
        anthroposphere_.routeedgelist.resize(static_cast<size_t>(routecount));

        for (RouteEdge& route : anthroposphere_.routeedgelist)
        {
            readvariable(infile, route.fromSettlementId);
            readvariable(infile, route.toSettlementId);
            std::uint8_t mode = 0;
            readvariable(infile, mode);
            if (mode > static_cast<std::uint8_t>(RouteMode::sea))
                mode = static_cast<std::uint8_t>(RouteMode::land);
            route.mode = static_cast<RouteMode>(mode);
            readvariable(infile, route.cost);
            readvariable(infile, route.capacity);
            readvariable(infile, route.traffic);
        }

        int eventcount = 0;
        readvariable(infile, eventcount);
        eventcount = std::max(0, eventcount);
        anthroposphere_.historyeventlist.resize(static_cast<size_t>(eventcount));

        for (HistoryEvent& event : anthroposphere_.historyeventlist)
        {
            readvariable(infile, event.year);
            readstring(event.type);
            readvariable(infile, event.primarySettlementId);
            readvariable(infile, event.primaryPolityId);
            readvariable(infile, event.secondarySettlementId);
            readvariable(infile, event.secondaryPolityId);
            readstring(event.summary);
        }
    }
    else
    {
        parallelforrows(0, body_.itswidth, [&](int startx, int endx)
        {
            for (int i = startx; i <= endx; i++)
            {
                for (int j = 0; j <= body_.itsheight; j++)
                {
                    anthroposphere_.settlement_suitability[i][j] = 0;
                    anthroposphere_.social_infrastructure[i][j] = 0;
                    anthroposphere_.agricultural_capacity[i][j] = 0;
                    anthroposphere_.route_traffic[i][j] = 0;
                    anthroposphere_.river_access[i][j] = 0;
                    anthroposphere_.harbor_score[i][j] = 0;
                    anthroposphere_.owner_settlement_id[i][j] = -1;
                    anthroposphere_.owner_polity_id[i][j] = -1;
                }
            }
        }, 64);

        clearsocialstate();
    }

    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
    {
        readlegacyvectordata(infile, climate_.seasonaltempmaps[season]);
        readlegacyvectordata(infile, climate_.seasonalrainmaps[season]);
        climate_.seasonalrainfloatmaps[season] = climate_.seasonalrainmaps[season];
        readlegacyvectordata(infile, climate_.seasonalpressuremaps[season]);
        readlegacyvectordata(infile, climate_.seasonaluwindmaps[season]);
        readlegacyvectordata(infile, climate_.seasonalvwindmaps[season]);
        if (fileversion >= 13)
        {
            readlegacyvectordata(infile, climate_.seasonalupperheightmaps[season]);
            readlegacyvectordata(infile, climate_.seasonalupperuwindmaps[season]);
            readlegacyvectordata(infile, climate_.seasonaluppervwindmaps[season]);
            readlegacyvectordata(infile, climate_.seasonalverticalvelocitymaps[season], climatelegacyscale::verticalvelocity);
        }
        else
        {
            std::fill(climate_.seasonalupperheightmaps[season].begin(), climate_.seasonalupperheightmaps[season].end(), 0.0f);
            std::fill(climate_.seasonalupperuwindmaps[season].begin(), climate_.seasonalupperuwindmaps[season].end(), 0.0f);
            std::fill(climate_.seasonaluppervwindmaps[season].begin(), climate_.seasonaluppervwindmaps[season].end(), 0.0f);
            std::fill(climate_.seasonalverticalvelocitymaps[season].begin(), climate_.seasonalverticalvelocitymaps[season].end(), 0.0f);
        }
        readlegacyvectordata(infile, climate_.seasonalcurrentumaps[season], climatelegacyscale::current);
        readlegacyvectordata(infile, climate_.seasonalcurrentvmaps[season], climatelegacyscale::current);
        readlegacyvectordata(infile, climate_.seasonalsstmaps[season]);
        readlegacyvectordata(infile, climate_.seasonalevaporationmaps[season]);
        readlegacyvectordata(infile, climate_.seasonalmaritimeinfluencemaps[season], climatelegacyscale::maritimeinfluence);
        readlegacyvectordata(infile, climate_.seasonalmaritimethermalanomalymaps[season], climatelegacyscale::maritimethermalanomaly);
        readshortvectordata(infile, climate_.seasonalmaritimefetchmaps[season]);
        readlegacyvectordata(infile, climate_.seasonalmoisturemaps[season]);
        readlegacyvectordata(infile, climate_.seasonalconvergencemaps[season], climatelegacyscale::convergence);
        readlegacyvectordata(infile, climate_.seasonalupliftmaps[season], climatelegacyscale::orographicmotion);
        readlegacyvectordata(infile, climate_.seasonalsubsidencemaps[season], climatelegacyscale::orographicmotion);
    }

    for (int i = 0; i < ARRAYWIDTH; i++)
    {
        for (int j = 0; j < 6; j++)
            readvariable(infile, climate_.horselats[i][j]);
    }

    for (int i = 0; i < geosphere_.itscraterno; i++)
    {
        readvariable(infile, geosphere_.cratercentreslist[i].w);
        readvariable(infile, geosphere_.cratercentreslist[i].x);
        readvariable(infile, geosphere_.cratercentreslist[i].y);
        readvariable(infile, geosphere_.cratercentreslist[i].z);
    }

    if (fileversion >= 15)
    {
        readvariable(infile, geosphere_.itstectonictimeoriginstep);
        readvariable(infile, geosphere_.itstectonictimemyr);
        readvariable(infile, geosphere_.itstectonicdeltatimemyr);
        readdata(infile, geosphere_.tectoniccrustagemyrmap);
        readdata(infile, geosphere_.tectoniccrustthicknessmap);
        readdata(infile, geosphere_.tectoniccrustclassmap);
        readdata(infile, geosphere_.tectonicuplifttendencymap);
        readdata(infile, geosphere_.tectonicsubsidencetendencymap);
        readdata(infile, geosphere_.tectonicaccumulatedstrainmap);
        readdata(infile, geosphere_.tectonicboundarytypemap);
        readdata(infile, geosphere_.tectonicboundarydistancemap);
        readdata(infile, geosphere_.tectonicboundarysegmentidmap);
        readdata(infile, geosphere_.tectonicnearestboundaryidmap);
        readdata(infile, geosphere_.tectonicboundaryhistorymap);
        readdata(infile, geosphere_.tectonicdeformingregionidmap);
        readdata(infile, geosphere_.tectonicdeformingregiontypemap);
        readdata(infile, geosphere_.tectonicdeformationratemap);
        readdata(infile, geosphere_.tectonicdeformationvelocityxmap);
        readdata(infile, geosphere_.tectonicdeformationvelocityymap);

        if (fileversion >= 16)
        {
            if (fileversion == 17)
            {
                int retiredcounter = 0; // Always zero; dropped from version 18.
                readvariable(infile, retiredcounter);
            }
            readvariable(infile, geosphere_.itstectonicplatecount);
            readvariable(infile, geosphere_.itstectonicsealevelm);

            int boundarysegmentcount = 0;
            readvariable(infile, boundarysegmentcount);
            boundarysegmentcount = std::max(0, boundarysegmentcount);
            geosphere_.tectonicboundarysegmentlist.resize(static_cast<size_t>(boundarysegmentcount));
            for (TectonicBoundarySegment& segment : geosphere_.tectonicboundarysegmentlist)
            {
                readvariable(infile, segment.id);
                readvariable(infile, segment.leftPlateId);
                readvariable(infile, segment.rightPlateId);
                readvariable(infile, segment.cellCount);
                readvariable(infile, segment.persistenceSteps);
                readvariable(infile, segment.centroidX);
                readvariable(infile, segment.centroidY);
                readvariable(infile, segment.lengthCells);
                readvariable(infile, segment.averageNormalMotion);
                readvariable(infile, segment.averageShearMotion);
                readvariable(infile, segment.averageConvergenceScore);
                readvariable(infile, segment.averageDivergenceScore);
                readvariable(infile, segment.averageShearScore);
                std::uint8_t boundarytype = 0;
                std::uint8_t geologicregime = 0;
                readvariable(infile, boundarytype);
                readvariable(infile, geologicregime);
                segment.boundaryType = decodeboundarytype(boundarytype);
                segment.geologicRegime = decodegeologicregime(geologicregime);
                readvariable(infile, segment.ageMyr);
                segment.id = std::max(0, segment.id);
                segment.leftPlateId = std::max(0, segment.leftPlateId);
                segment.rightPlateId = std::max(0, segment.rightPlateId);
                segment.cellCount = std::max(0, segment.cellCount);
                segment.persistenceSteps = std::max(0, segment.persistenceSteps);
                segment.averageConvergenceScore = std::clamp(segment.averageConvergenceScore, 0, 255);
                segment.averageDivergenceScore = std::clamp(segment.averageDivergenceScore, 0, 255);
                segment.averageShearScore = std::clamp(segment.averageShearScore, 0, 255);
                segment.ageMyr = std::max(0.0, segment.ageMyr);
            }

            int deformingregioncount = 0;
            readvariable(infile, deformingregioncount);
            deformingregioncount = std::max(0, deformingregioncount);
            geosphere_.tectonicdeformingregionlist.resize(static_cast<size_t>(deformingregioncount));
            for (TectonicDeformingRegion& region : geosphere_.tectonicdeformingregionlist)
            {
                readvariable(infile, region.id);
                readvariable(infile, region.boundarySegmentId);
                readvariable(infile, region.primaryPlateId);
                readvariable(infile, region.secondaryPlateId);
                readvariable(infile, region.cellCount);
                readvariable(infile, region.persistenceSteps);
                readvariable(infile, region.centroidX);
                readvariable(infile, region.centroidY);
                readvariable(infile, region.averageDeformationRate);
                readvariable(infile, region.averageInterpolatedVelocityX);
                readvariable(infile, region.averageInterpolatedVelocityY);
                readvariable(infile, region.averageNormalMotion);
                readvariable(infile, region.averageShearMotion);
                std::uint8_t regiontype = 0;
                readvariable(infile, regiontype);
                region.type = decodedeformingregiontype(regiontype);
                readvariable(infile, region.ageMyr);
                region.id = std::max(0, region.id);
                region.boundarySegmentId = std::max(0, region.boundarySegmentId);
                region.primaryPlateId = std::max(0, region.primaryPlateId);
                region.secondaryPlateId = std::max(0, region.secondaryPlateId);
                region.cellCount = std::max(0, region.cellCount);
                region.persistenceSteps = std::max(0, region.persistenceSteps);
                region.ageMyr = std::max(0.0, region.ageMyr);
            }
        }
        else
        {
            geosphere_.itstectonicplatecount = 0;
            geosphere_.itstectonicsealevelm = 0;
            geosphere_.tectonicboundarysegmentlist.clear();
            geosphere_.tectonicdeformingregionlist.clear();
        }
    }
    else
    {
        geosphere_.itstectonictimeoriginstep = 0;
        geosphere_.itstectonictimemyr = 0.0f;
        geosphere_.itstectonicdeltatimemyr = 0.0f;
        geosphere_.itstectonicplatecount = 0;
        geosphere_.itstectonicsealevelm = 0;
        geosphere_.tectonicboundarysegmentlist.clear();
        geosphere_.tectonicdeformingregionlist.clear();

        parallelforrows(0, body_.itswidth, [&](int startx, int endx)
        {
            for (int i = startx; i <= endx; i++)
            {
                for (int j = 0; j <= body_.itsheight; j++)
                {
                    geosphere_.tectoniccrustagemyrmap[i][j] = 0.0f;
                    geosphere_.tectoniccrustthicknessmap[i][j] = 0.0f;
                    geosphere_.tectoniccrustclassmap[i][j] = static_cast<std::uint8_t>(CrustClass::none);
                    geosphere_.tectonicuplifttendencymap[i][j] = 0.0f;
                    geosphere_.tectonicsubsidencetendencymap[i][j] = 0.0f;
                    geosphere_.tectonicaccumulatedstrainmap[i][j] = 0.0f;
                    geosphere_.tectonicboundarytypemap[i][j] = static_cast<std::uint8_t>(BoundaryType::none);
                    geosphere_.tectonicboundarydistancemap[i][j] = 0;
                    geosphere_.tectonicboundarysegmentidmap[i][j] = 0;
                    geosphere_.tectonicnearestboundaryidmap[i][j] = 0;
                    geosphere_.tectonicboundaryhistorymap[i][j] = 0.0f;
                    geosphere_.tectonicdeformingregionidmap[i][j] = 0;
                    geosphere_.tectonicdeformingregiontypemap[i][j] = static_cast<std::uint8_t>(DeformingRegionType::none);
                    geosphere_.tectonicdeformationratemap[i][j] = 0.0f;
                    geosphere_.tectonicdeformationvelocityxmap[i][j] = 0.0f;
                    geosphere_.tectonicdeformationvelocityymap[i][j] = 0.0f;
                }
            }
        }, 64);
    }

    int weatherstatesize = 0;
    readvariable(infile, weatherstatesize);
    if (weatherstatesize < 0 || weatherstatesize > 16 * 1024 * 1024)
        return false;
    climate_.weatheranomalystatebytes.resize(static_cast<std::size_t>(weatherstatesize));
    for (std::uint8_t& byte : climate_.weatheranomalystatebytes)
    {
        int value = 0;
        readvariable(infile, value);
        if (value < 0 || value > 255)
            return false;
        byte = static_cast<std::uint8_t>(value);
    }
    climateweather::ShallowWaterState loadedweather;
    if (!climate_.weatheranomalystatebytes.empty() && !climateweather::deserializeState(climate_.weatheranomalystatebytes, loadedweather))
        return false;

    if (!readmineraldeposits(infile))
        return false;

    setmaxriverflow();

    if (!infile.good())
    {
        cerr << "Error reading world '" << filename << "'" << endl;
    }

    return 1;
}

// Private member functions.

int planet::wrapx(int x) const // This wraps X coordinates so they point to proper locations on the map.
{
    while (x > body_.itswidth) // If it's too large, wrap it.
    {
        x = x - body_.itswidth;
    }

    while (x < 0) // If it's too small, wrap it.
    {
        x = x + body_.itswidth;
    }

    return(x);
}

int planet::clipy(int y) const // This clips Y coordinates so they can't be off the map.
{
    if (y < 0)
        y = 0;

    if (y > body_.itsheight)
        y = body_.itsheight;

    return(y);
}

namespace
{
// Rows beyond the poles count as zero in the averages, as they did when grids
// were zero-filled 2048 x 1024 arrays; the averages keep their former results.
template<typename T>
T smoothingcell(const worldgrid<T>& grid, int x, int y)
{
    return y < 0 || y >= grid.rows() ? T{} : grid[x][y];
}
}

void planet::smooth(worldgrid<int>& arr, int amount, bool vary, bool avoidmountains) // This smoothes the given array by the given amount.
{
    for (int i = 0; i <= body_.itswidth; i++)
    {
        for (int j = 1; j < body_.itsheight; j++)
        {
            if (avoidmountains == 0 || geosphere_.mountainheights[i][j] == 0)
            {
                int crount = 0;
                int ave = 0;

                for (int k = i - amount; k <= i + amount; k++)
                {
                    int kk = k;

                    if (kk < 0)
                        kk = body_.itswidth;

                    if (kk > body_.itswidth)
                        kk = 0;

                    for (int l = j - amount; l <= j + amount; l++)
                    {
                        ave = ave + smoothingcell(geosphere_.mapnom, kk, l);
                        crount++;
                    }
                }
                ave = ave / crount;

                if (ave > 0 && ave < geosphere_.itsmaxheight)
                    geosphere_.mapnom[i][j] = ave;
            }
        }
    }
}

void planet::smooth(worldgrid<short>& arr, int amount, bool vary, bool avoidmountains) // This smoothes the given array by the given amount.
{
    for (int i = 0; i <= body_.itswidth; i++)
    {
        for (int j = 1; j < body_.itsheight; j++)
        {
            if (avoidmountains == 0 || geosphere_.mountainheights[i][j] == 0)
            {
                int crount = 0;
                int ave = 0;

                for (int k = i - amount; k <= i + amount; k++)
                {
                    int kk = k;

                    if (kk < 0)
                        kk = body_.itswidth;

                    if (kk > body_.itswidth)
                        kk = 0;

                    for (int l = j - amount; l <= j + amount; l++)
                    {
                        ave = ave + smoothingcell(geosphere_.mapnom, kk, l);
                        crount++;
                    }
                }
                ave = ave / crount;

                if (ave > 0 && ave < geosphere_.itsmaxheight)
                    geosphere_.mapnom[i][j] = ave;
            }
        }
    }
}

// This does the same, but only over land.

void planet::smoothoverland(worldgrid<int>& arr, int amount, bool uponly)
{
    for (int i = 0; i <= body_.itswidth; i++)
    {
        for (int j = 1; j < body_.itsheight; j++)
        {
            if (sea(i, j) == 0)
            {
                int crount = 0;
                int ave = 0;

                for (int k = i - amount; k <= i + amount; k++)
                {
                    int kk = k;

                    if (kk < 0)
                        kk = body_.itswidth;

                    if (kk > body_.itswidth)
                        kk = 0;

                    for (int l = j - amount; l <= j + amount; l++)
                    {
                        ave = ave + smoothingcell(arr, kk, l);
                        crount++;
                    }
                }

                if (crount > 0)
                {
                    ave = ave / crount;

                    if (ave > 0 && ave < geosphere_.itsmaxheight)
                    {
                        if (uponly == 0)
                            arr[i][j] = ave;
                        else
                        {
                            if (ave > arr[i][j])
                                arr[i][j] = ave;
                        }
                    }
                }
            }
        }
    }
}

void planet::smoothoverland(worldgrid<short>& arr, int amount, bool uponly)
{
    for (int i = 0; i <= body_.itswidth; i++)
    {
        for (int j = 1; j < body_.itsheight; j++)
        {
            if (sea(i, j) == 0)
            {
                int crount = 0;
                int ave = 0;

                for (int k = i - amount; k <= i + amount; k++)
                {
                    int kk = k;

                    if (kk < 0)
                        kk = body_.itswidth;

                    if (kk > body_.itswidth)
                        kk = 0;

                    for (int l = j - amount; l <= j + amount; l++)
                    {
                        ave = ave + (int)smoothingcell(arr, kk, l);
                        crount++;
                    }
                }

                if (crount > 0)
                {
                    ave = ave / crount;

                    if (ave > 0 && ave < geosphere_.itsmaxheight)
                    {
                        if (uponly == 0)
                            arr[i][j] = (short)ave;
                        else
                        {
                            if (ave > arr[i][j])
                                arr[i][j] = (short)ave;
                        }
                    }
                }
            }
        }
    }
}

// This function shifts everything in an array to the left by a given number of pixels.

template<typename T> void planet::shift(worldgrid<T>& arr, int offset)
{
    vector<vector<T>> dummy(ARRAYWIDTH, vector<T>(ARRAYHEIGHT, 0));

    //T dummy[ARRAYWIDTH][ARRAYHEIGHT];

    for (int i = 0; i <= body_.itswidth; i++)
    {
        for (int j = 0; j <= body_.itsheight; j++)
            dummy[i][j] = arr[i][j];
    }

    for (int i = 0; i <= body_.itswidth; i++)
    {
        int ii = i + offset;

        if (ii<0 || ii>body_.itswidth)
            ii = wrap(ii, body_.itswidth);

        for (int j = 0; j <= body_.itsheight; j++)
            arr[i][j] = dummy[ii][j];
    }
}

// Functions for saving member variables.

template<typename T> void write_val(T const val, ostream& out) { // default
    out << val;
}
void write_int_val(int val, ostream& out) {
    if (val < 0) { out.put('-'); val = -val; } // negative
    if (val < 10) { out.put('0' + char(val)); } // 1 digit
    else if (val < 100) { out.put('0' + char(val / 10)); out.put('0' + char(val % 10)); } // 2 digits
    else if (val < 1000) { out.put('0' + char(val / 100)); out.put('0' + char((val / 10) % 10)); out.put('0' + char(val % 10)); } // 3 digits
    else { out << val; } // 4+ digits
}
void write_val(int   const val, ostream& out) { write_int_val(val, out); }
void write_val(short const val, ostream& out) { write_int_val(val, out); }

void write_val(bool const val, ostream& out) {
    out.put(val ? '1' : '0');
}

template<typename T> void planet::writevariable(ofstream& outfile, T val)
{
    write_val(val, outfile);
    outfile.put('\n');
}

// Functions for saving member arrays.

// Grids are written column by column over the world size, the layout of the former arrays.
template<typename T> void planet::writedata(ofstream& outfile, const worldgrid<T>& arr)
{
    for (int i = 0; i <= body_.itswidth; i++)
    {
        for (int j = 0; j <= body_.itsheight; j++)
        {
            write_val(arr[i][j], outfile);
            outfile.put('\n');
        }
    }
}

// Functions for loading member variables.

void read_val(string const& line, int& val)
{
    val = stoi(line);
}
void read_val(string const& line, bool& val)
{
    val = stob(line);
}
void read_val(string const& line, short& val)
{
    val = stos(line);
}
void read_val(string const& line, unsigned short& val)
{
    val = stous(line);
}
void read_val(string const& line, float& val)
{
    val = stof(line);
}
void read_val(string const& line, double& val)
{
    val = stod(line);
}
void read_val(string const& line, long& val)
{
    val = stol(line);
}
void read_val(string const& line, char& val)
{
    val = stoc(line);
}
void read_val(string const& line, unsigned char& val)
{
    val = stouc(line);
}

template<typename T> void planet::readvariable(ifstream& infile, T& val)
{
    getline(infile, line_for_file_read);
    read_val(line_for_file_read, val);
}

// Functions for loading member arrays.

template<typename T> void planet::readdata(ifstream& infile, worldgrid<T>& arr)
{
    for (int i = 0; i <= body_.itswidth; i++)
    {
        for (int j = 0; j <= body_.itsheight; j++)
        {
            getline(infile, line_for_file_read);
            read_val(line_for_file_read, arr[i][j]);
        }
    }
}

climateinputs planet::climateinputsview() const
{
    climateinputs inputs;
    inputs.width = body_.itswidth;
    inputs.height = body_.itsheight;
    inputs.seed = body_.itsseed;
    inputs.rotation = body_.itsrotation;
    inputs.tilt = body_.itstilt;
    inputs.eccentricity = body_.itseccentricity;
    inputs.perihelion = body_.itsperihelion;
    inputs.gravity = body_.itsgravity;
    inputs.lunar = body_.itslunar;
    inputs.tempdecrease = atmosphere_.itstempdecrease;
    inputs.northpolaradjust = atmosphere_.itsnorthpolaradjust;
    inputs.southpolaradjust = atmosphere_.itssouthpolaradjust;
    inputs.averagetemp = atmosphere_.itsaveragetemp;
    inputs.waterpickup = atmosphere_.itswaterpickup;
    inputs.maxelevation = geosphere_.itsmaxheight;
    inputs.sealevel = geosphere_.itssealevel;
    inputs.landtotal = geosphere_.itslandtotal;
    inputs.seatotal = geosphere_.itsseatotal;
    inputs.rows = body_.itsheight + 1;
    inputs.nom = &geosphere_.mapnom[0][0];
    inputs.mountainheights = &geosphere_.mountainheights[0][0];
    inputs.extraelev = &geosphere_.extraelevmap[0][0];
    inputs.craterrims = &geosphere_.craterrims[0][0];
    inputs.cratercentres = &geosphere_.cratercentres[0][0];
    inputs.volcano = &geosphere_.volcanomap[0][0];
    inputs.oceanridgeheight = &geosphere_.oceanridgeheightmap[0][0];
    inputs.lakesurface = &hydrosphere_.lakemap[0][0];
    inputs.riftlakesurface = &hydrosphere_.riftlakemapsurface[0][0];
    inputs.special = &hydrosphere_.specials[0][0];
    return inputs;
}

// Binding on every call keeps the inputs pointing at this world after copies.
climateworld& planet::climate()
{
    climate_.bindinputs(climateinputsview());
    return climate_;
}

const climateworld& planet::climate() const
{
    climate_.bindinputs(climateinputsview());
    return climate_;
}
