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

#include "planet.hpp"
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
    itswidth = 599;
    itsheight = 399;
    itstectonictimeoriginstep = 0;
    itstectonictimemyr = 0.0f;
    itstectonicdeltatimemyr = 0.0f;
    itstectoniccyclecount = 0;
    itstectonicplatecount = 0;
    itstectonicsealevelm = 0;
    resizeseasonalclimatefields();
}

planet::~planet()
{
}

void planet::resizeseasonalclimatefields()
{
    const int cellcount = (itswidth + 1) * (itsheight + 1);

    auto resizefield = [cellcount](std::array<std::vector<short>, CLIMATESEASONCOUNT>& field)
    {
        for (std::vector<short>& season : field)
            season.assign(cellcount, 0);
    };
    auto resizefloatfield = [cellcount](
        std::array<std::vector<float>, CLIMATESEASONCOUNT>& field)
    {
        for (std::vector<float>& season : field)
            season.assign(cellcount, 0.0f);
    };

    resizefield(seasonaltempmaps);
    resizefield(seasonalrainmaps);
    resizefloatfield(seasonalrainfloatmaps);
    resizefield(seasonalpressuremaps);
    resizefield(seasonaluwindmaps);
    resizefield(seasonalvwindmaps);
    resizefield(seasonalupperheightmaps);
    resizefield(seasonalupperuwindmaps);
    resizefield(seasonaluppervwindmaps);
    resizefield(seasonalverticalvelocitymaps);
    resizefield(seasonalcurrentumaps);
    resizefield(seasonalcurrentvmaps);
    resizefield(seasonalsstmaps);
    resizefield(seasonalevaporationmaps);
    resizefield(seasonalmaritimeinfluencemaps);
    resizefield(seasonalmaritimethermalanomalymaps);
    resizefield(seasonalmaritimefetchmaps);
    resizefield(seasonalmoisturemaps);
    resizefield(seasonalconvergencemaps);
    resizefield(seasonalupliftmaps);
    resizefield(seasonalsubsidencemaps);
}

void planet::cleartectonicprovenanceinternal()
{
    itstectonictimeoriginstep = 0;
    itstectonictimemyr = 0.0f;
    itstectonicdeltatimemyr = 0.0f;
    itstectoniccyclecount = 0;
    itstectonicplatecount = 0;
    itstectonicsealevelm = 0;
    tectonicboundarysegmentlist.clear();
    tectonicdeformingregionlist.clear();

    parallelforrows(0, ARRAYWIDTH - 1, [&](int startx, int endx)
    {
        for (int i = startx; i <= endx; i++)
        {
            for (int j = 0; j < ARRAYHEIGHT; j++)
            {
                geologicregimemap[i][j] = static_cast<std::uint8_t>(GeologicRegime::stable);
                tectonicconvergencemap[i][j] = 0;
                tectonicdivergencemap[i][j] = 0;
                tectonicshearmap[i][j] = 0;
                tectoniccrustagemyrmap[i][j] = 0.0f;
                tectoniccrustthicknessmap[i][j] = 0.0f;
                tectoniccrustclassmap[i][j] = static_cast<std::uint8_t>(CrustClass::none);
                tectonicuplifttendencymap[i][j] = 0.0f;
                tectonicsubsidencetendencymap[i][j] = 0.0f;
                tectonicaccumulatedstrainmap[i][j] = 0.0f;
                tectonicboundarytypemap[i][j] = static_cast<std::uint8_t>(BoundaryType::none);
                tectonicboundarydistancemap[i][j] = 0;
                tectonicboundarysegmentidmap[i][j] = 0;
                tectonicnearestboundaryidmap[i][j] = 0;
                tectonicboundaryhistorymap[i][j] = 0.0f;
                tectonicdeformingregionidmap[i][j] = 0;
                tectonicdeformingregiontypemap[i][j] = static_cast<std::uint8_t>(DeformingRegionType::none);
                tectonicdeformationratemap[i][j] = 0.0f;
                tectonicdeformationvelocityxmap[i][j] = 0.0f;
                tectonicdeformationvelocityymap[i][j] = 0.0f;
            }
        }
    }, 64);
}

void planet::cleartectonicprovenance()
{
    cleartectonicprovenanceinternal();
}

void planet::writeshortvectordata(ofstream& outfile, const std::vector<short>& arr)
{
    outfile.write(reinterpret_cast<const char*>(arr.data()), arr.size() * sizeof(short));
}

void planet::readshortvectordata(ifstream& infile, std::vector<short>& arr)
{
    infile.read(reinterpret_cast<char*>(arr.data()), arr.size() * sizeof(short));
}

void planet::syncseasonalclimatefromlegacy()
{
    bool hasexplicittransitiontemp = false;
    bool hasexplicittransitionrain = false;

    for (const short value : seasonaltempmaps[seasonapril])
    {
        if (value != 0)
        {
            hasexplicittransitiontemp = true;
            break;
        }
    }

    if (hasexplicittransitiontemp == false)
    {
        for (const short value : seasonaltempmaps[seasonoctober])
        {
            if (value != 0)
            {
                hasexplicittransitiontemp = true;
                break;
            }
        }
    }

    for (const short value : seasonalrainmaps[seasonapril])
    {
        if (value != 0)
        {
            hasexplicittransitionrain = true;
            break;
        }
    }

    if (hasexplicittransitionrain == false)
    {
        for (const short value : seasonalrainmaps[seasonoctober])
        {
            if (value != 0)
            {
                hasexplicittransitionrain = true;
                break;
            }
        }
    }

    for (int y = 0; y <= itsheight; y++)
    {
        for (int x = 0; x <= itswidth; x++)
        {
            const int index = seasonalclimateindex(x, y);

            const float januarytemp = static_cast<float>(jantempmap[x][y]);
            const float julytemp = static_cast<float>(jultempmap[x][y]);
            const float januaryrain = static_cast<float>(janrainmap[x][y]);
            const float julyrain = static_cast<float>(julrainmap[x][y]);

            seasonaltempmaps[seasonjanuary][index] = static_cast<short>(jantempmap[x][y]);
            seasonaltempmaps[seasonjuly][index] = static_cast<short>(jultempmap[x][y]);
            seasonalrainmaps[seasonjanuary][index] = static_cast<short>(janrainmap[x][y]);
            seasonalrainmaps[seasonjuly][index] = static_cast<short>(julrainmap[x][y]);

            float summertemp = julytemp;
            float wintertemp = januarytemp;

            if (itsperihelion == 1)
            {
                summertemp = januarytemp;
                wintertemp = julytemp;
            }

            const float winterstrength = 0.5f + itseccentricity * 0.5f;
            const float summerstrength = 1.0f - winterstrength;
            float transitiontemp = summertemp * summerstrength + wintertemp * winterstrength;

            float fourseason = itstilt * 0.294592f - 2.45428f;
            float lat = static_cast<float>(y);

            if (y > itsheight / 2.0f)
                lat = static_cast<float>(itsheight - y);

            const float fourseasonstrength = lat / (static_cast<float>(itsheight) / 2.0f);
            const float transitiontempdiff = (fourseason * fourseasonstrength) / 2.0f;

            if (hasexplicittransitiontemp == false)
            {
                seasonaltempmaps[seasonapril][index] = static_cast<short>(transitiontemp + transitiontempdiff);
                seasonaltempmaps[seasonoctober][index] = static_cast<short>(transitiontemp + transitiontempdiff);
            }

            if (hasexplicittransitionrain == false)
            {
                float apriljanfactor = 0.5f;
                float apriljulfactor = 0.5f;

                if (jultempmap[x][y] > jantempmap[x][y] && julrainmap[x][y] > janrainmap[x][y] && julyrain > 0.0f)
                {
                    const float monsoonfactor = 1.0f - januaryrain / julyrain;
                    apriljanfactor = monsoonfactor * 0.9f;
                    apriljulfactor = 1.0f - apriljanfactor;
                }

                if (jultempmap[x][y] < jantempmap[x][y] && julrainmap[x][y] < janrainmap[x][y] && januaryrain > 0.0f)
                {
                    const float monsoonfactor = 1.0f - julyrain / januaryrain;
                    apriljanfactor = monsoonfactor * 0.7f;
                    apriljulfactor = 1.0f - apriljanfactor;
                }

                seasonalrainmaps[seasonapril][index] = static_cast<short>(januaryrain * apriljanfactor + julyrain * apriljulfactor);

                float octoberjanfactor = 0.5f;
                float octoberjulfactor = 0.5f;

                if (jultempmap[x][y] > jantempmap[x][y] && julrainmap[x][y] > janrainmap[x][y] && julyrain > 0.0f)
                {
                    const float monsoonfactor = 1.0f - januaryrain / julyrain;
                    octoberjulfactor = monsoonfactor * 0.7f;
                    octoberjanfactor = 1.0f - octoberjulfactor;
                }

                if (jultempmap[x][y] < jantempmap[x][y] && julrainmap[x][y] < janrainmap[x][y] && januaryrain > 0.0f)
                {
                    const float monsoonfactor = 1.0f - julyrain / januaryrain;
                    octoberjulfactor = monsoonfactor * 0.9f;
                    octoberjanfactor = 1.0f - octoberjulfactor;
                }

                seasonalrainmaps[seasonoctober][index] = static_cast<short>(januaryrain * octoberjanfactor + julyrain * octoberjulfactor);
            }
        }
    }
}

bool planet::outline(int x, int y) const
{
    if (y<1 || y>itsheight - 1)
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
    if (y<1 || y>itsheight - 1 || x<0 || x>itswidth)
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
    if ( x<0 || x>itswidth)
        return;

    float fx = (float)x;

    float worldwidth = (float)itswidth + 1.0f;

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
    if (y<0 || y>itsheight)
        return;

    bool hemisphere = 0; // 1 for southern

    float fy = (float)y;

    float worldheight = (float)itsheight;
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

    float equator = (float)itsheight / 2.0f;

    float pixelsperlat = equator / 90.0f;

    float y = flat * pixelsperlat;

    return (int)y;
}

// slightly more complicated accessor functions

bool planet::seawrap(int x, int y) const
{
    x = wrapx(x);
    y = clipy(y);

    if (mapnom[x][y] <= itssealevel && lakemap[x][y] == 0)
        return 1;

    else
        return 0;
}

bool planet::outlinewrap(int x, int y) const
{
    x = wrapx(x);

    if (y<1 || y>itsheight - 1)
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
    x = wrapx(x);
    y = clipy(y);

    return mountainheights[x][y];
}

// Other public functions.

void planet::clear()
{
    resizeseasonalclimatefields();
    itstectonictimeoriginstep = 0;
    itstectonictimemyr = 0.0f;
    itstectonicdeltatimemyr = 0.0f;
    itstectoniccyclecount = 0;
    itstectonicplatecount = 0;
    itstectonicsealevelm = 0;
    tectonicboundarysegmentlist.clear();
    tectonicdeformingregionlist.clear();

    parallelforrows(0, ARRAYWIDTH - 1, [&](int startx, int endx)
    {
        for (int i = startx; i <= endx; i++) // Set all the maps to 0.
        {
            for (int j = 0; j < ARRAYHEIGHT; j++)
            {
                jantempmap[i][j] = 0;
                jultempmap[i][j] = 0;
                climatemap[i][j] = 0;
                biomemap[i][j] = 0;
                janrainmap[i][j] = 0;
                julrainmap[i][j] = 0;
                janmountainrainmap[i][j] = 0;
                julmountainrainmap[i][j] = 0;
                janmountainraindirmap[i][j] = 0;
                julmountainraindirmap[i][j] = 0;
                seaicemap[i][j] = 0;
                rivermapdir[i][j] = 0;
                rivermapjan[i][j] = 0;
                rivermapjul[i][j] = 0;
                windmap[i][j] = 0;
                lakemap[i][j] = 0;
                roughnessmap[i][j] = 0;
                mountainridges[i][j] = 0;
                mountainheights[i][j] = 0;
                craterrims[i][j] = 0;
                cratercentres[i][j] = 0;
                mapnom[i][j] = 0;
                tidalmap[i][j] = 0;
                riftlakemapsurface[i][j] = 0;
                riftlakemapbed[i][j] = 0;
                specials[i][j] = 0;
                geologicregimemap[i][j] = static_cast<std::uint8_t>(GeologicRegime::stable);
                tectonicconvergencemap[i][j] = 0;
                tectonicdivergencemap[i][j] = 0;
                tectonicshearmap[i][j] = 0;
                tectoniccrustagemyrmap[i][j] = 0.0f;
                tectoniccrustthicknessmap[i][j] = 0.0f;
                tectoniccrustclassmap[i][j] = static_cast<std::uint8_t>(CrustClass::none);
                tectonicuplifttendencymap[i][j] = 0.0f;
                tectonicsubsidencetendencymap[i][j] = 0.0f;
                tectonicaccumulatedstrainmap[i][j] = 0.0f;
                tectonicboundarytypemap[i][j] = static_cast<std::uint8_t>(BoundaryType::none);
                tectonicboundarydistancemap[i][j] = 0;
                tectonicboundarysegmentidmap[i][j] = 0;
                tectonicnearestboundaryidmap[i][j] = 0;
                tectonicboundaryhistorymap[i][j] = 0.0f;
                tectonicdeformingregionidmap[i][j] = 0;
                tectonicdeformingregiontypemap[i][j] = static_cast<std::uint8_t>(DeformingRegionType::none);
                tectonicdeformationratemap[i][j] = 0.0f;
                tectonicdeformationvelocityxmap[i][j] = 0.0f;
                tectonicdeformationvelocityymap[i][j] = 0.0f;
                basinclassmap[i][j] = static_cast<std::uint8_t>(BasinClass::none);
                erosionpotentialmap[i][j] = 0;
                depositionpotentialmap[i][j] = 0;
                floodplainfertilitymap[i][j] = 0;
                metalorereservemap[i][j] = 0;
                placerreservemap[i][j] = 0;
                evaporitereservemap[i][j] = 0;
                volcanicreservemap[i][j] = 0;
                fisheryreservemap[i][j] = 0;
                settlement_suitability[i][j] = 0;
                social_infrastructure[i][j] = 0;
                agricultural_capacity[i][j] = 0;
                route_traffic[i][j] = 0;
                river_access[i][j] = 0;
                harbor_score[i][j] = 0;
                owner_settlement_id[i][j] = -1;
                owner_polity_id[i][j] = -1;
                extraelevmap[i][j] = 0;
                deltamapdir[i][j] = 0;
                deltamapjan[i][j] = 0;
                deltamapjul[i][j] = 0;
                oceanridgemap[i][j] = 0;
                oceanridgeheightmap[i][j] = 0;
                oceanriftmap[i][j] = 0;
                oceanridgeoffsetmap[i][j] = 0;
                islandmap[i][j] = 0;
                noshademap[i][j] = 0;
                mountainislandmap[i][j] = 0;
                volcanomap[i][j] = 0;
                stratomap[i][j] = 0;
                noisemap[i][j] = 0;
                testmap[i][j] = 0;
            }

            for (int j = 0; j < 6; j++)
            {
                horselats[i][j] = 0;
            }
        }
    }, 64);

    for (int i = 0; i < MAXCRATERS; i++)
    {
        cratercentreslist[i].w = 0;
        cratercentreslist[i].x = 0;
        cratercentreslist[i].y = 0;
        cratercentreslist[i].z = 0;
    }

    clearsocialstate();
}

void planet::smoothnom(int amount)
{
    smooth(mapnom, amount, 1, 1);
}

void planet::smoothextraelev(int amount)
{
    smoothoverland(extraelevmap, amount, 0);
}

void planet::shiftterrain(int offset)
{
    shift(mapnom, offset);
    shift(geologicregimemap, offset);
    shift(tectonicconvergencemap, offset);
    shift(tectonicdivergencemap, offset);
    shift(tectonicshearmap, offset);
    shift(tectoniccrustagemyrmap, offset);
    shift(tectoniccrustthicknessmap, offset);
    shift(tectoniccrustclassmap, offset);
    shift(tectonicuplifttendencymap, offset);
    shift(tectonicsubsidencetendencymap, offset);
    shift(tectonicaccumulatedstrainmap, offset);
    shift(tectonicboundarytypemap, offset);
    shift(tectonicboundarydistancemap, offset);
    shift(tectonicboundarysegmentidmap, offset);
    shift(tectonicnearestboundaryidmap, offset);
    shift(tectonicboundaryhistorymap, offset);
    shift(tectonicdeformingregionidmap, offset);
    shift(tectonicdeformingregiontypemap, offset);
    shift(tectonicdeformationratemap, offset);
    shift(tectonicdeformationvelocityxmap, offset);
    shift(tectonicdeformationvelocityymap, offset);
    shift(basinclassmap, offset);
    shift(erosionpotentialmap, offset);
    shift(depositionpotentialmap, offset);
    shift(floodplainfertilitymap, offset);
    shift(metalorereservemap, offset);
    shift(placerreservemap, offset);
    shift(evaporitereservemap, offset);
    shift(volcanicreservemap, offset);
    shift(fisheryreservemap, offset);
    shift(settlement_suitability, offset);
    shift(social_infrastructure, offset);
    shift(agricultural_capacity, offset);
    shift(route_traffic, offset);
    shift(river_access, offset);
    shift(harbor_score, offset);
    shift(owner_settlement_id, offset);
    shift(owner_polity_id, offset);
    shift(mountainheights, offset);
    shift(mountainridges, offset);
    shift(craterrims, offset);
    shift(cratercentres, offset);
    shift(extraelevmap, offset);
    shift(oceanridgemap, offset);
    shift(oceanridgeheightmap, offset);
    shift(oceanriftmap, offset);
    shift(oceanridgeanglemap, offset);
    shift(mountainislandmap, offset);
    shift(noshademap, offset);
    shift(volcanomap, offset);
    shift(stratomap, offset);
    shift(testmap, offset);

    for (int i = 0; i < itscraterno; i++)
    {
        cratercentreslist[i].x = cratercentreslist[i].x - offset;

        if (cratercentreslist[i].x < 0)
            cratercentreslist[i].x = cratercentreslist[i].x + itswidth;
    }

    for (Settlement& settlement : settlementlist)
    {
        settlement.x -= offset;

        while (settlement.x < 0)
            settlement.x += itswidth;

        while (settlement.x > itswidth)
            settlement.x -= itswidth;
    }

    const float worldspan = static_cast<float>(itswidth + 1);
    const float offsetf = static_cast<float>(offset);
    const auto shiftcentroidx = [worldspan, offsetf](float x) -> float
    {
        float shifted = fmodf(x - offsetf, worldspan);

        if (shifted < 0.0f)
            shifted += worldspan;

        return shifted;
    };

    for (TectonicBoundarySegment& segment : tectonicboundarysegmentlist)
        segment.centroidX = shiftcentroidx(segment.centroidX);

    for (TectonicDeformingRegion& region : tectonicdeformingregionlist)
        region.centroidX = shiftcentroidx(region.centroidX);
}

void planet::smoothrainmaps(int amount)
{
    smoothoverland(janrainmap, amount, 0);
    smoothoverland(julrainmap, amount, 0);
}

void planet::setmaxriverflow()
{
    int largest = 0;
    int current = 0;

    for (int i = 0; i <= itswidth; i++)
    {
        for (int j = 0; j <= itsheight; j++)
        {
            current = riveraveflow(i, j);

            if (current > largest)
                largest = current;
        }
    }
    itsmaxriverflow = largest;
}

void planet::saveworld(string filename)
{
#ifdef ENABLE_PROFILER
    highres_timer_t timer("Save World"); // 26.5s => 10.9s
#endif
    ofstream outfile;
    outfile.open(filename, ios::out);

    writevariable(outfile, itssaveversion);
    writevariable(outfile, itssize);
    writevariable(outfile, itswidth);
    writevariable(outfile, itsheight);
    writevariable(outfile, itsseed);
    writevariable(outfile, itsrotation);
    writevariable(outfile, itstilt);
    writevariable(outfile, itseccentricity);
    writevariable(outfile, itsperihelion);
    writevariable(outfile, itsgravity);
    writevariable(outfile, itslunar);
    writevariable(outfile, itstempdecrease);
    writevariable(outfile, itsnorthpolaradjust);
    writevariable(outfile, itssouthpolaradjust);
    writevariable(outfile, itsaveragetemp);
    writevariable(outfile, itsnorthpolartemp);
    writevariable(outfile, itssouthpolartemp);
    writevariable(outfile, itseqtemp);
    writevariable(outfile, itswaterpickup);
    writevariable(outfile, itsriverfactor);
    writevariable(outfile, itsriverlandreduce);
    writevariable(outfile, itsestuarylimit);
    writevariable(outfile, itsglacialtemp);
    writevariable(outfile, itsglaciertemp);
    writevariable(outfile, itsmountainreduce);
    writevariable(outfile, itsclimateno);
    writevariable(outfile, itsmaxheight);
    writevariable(outfile, itssealevel);
    writevariable(outfile, itslandtotal);
    writevariable(outfile, itsseatotal);
    writevariable(outfile, itscraterno);

    writevariable(outfile, itslandshading);
    writevariable(outfile, itslakeshading);
    writevariable(outfile, itsseashading);
    writevariable(outfile, itsshadingdir);
    writevariable(outfile, itssnowchange);
    writevariable(outfile, itsseaiceappearance);
    writevariable(outfile, itslandmarbling);
    writevariable(outfile, itslakemarbling);
    writevariable(outfile, itsseamarbling);
    writevariable(outfile, itsminriverflowglobal);
    writevariable(outfile, itsminriverflowregional);
    writevariable(outfile, itsmangroves);
    writevariable(outfile, itscolourcliffs);
    writevariable(outfile, itsseaice1);
    writevariable(outfile, itsseaice2);
    writevariable(outfile, itsseaice3);
    writevariable(outfile, itsocean1);
    writevariable(outfile, itsocean2);
    writevariable(outfile, itsocean3);
    writevariable(outfile, itsdeepocean1);
    writevariable(outfile, itsdeepocean2);
    writevariable(outfile, itsdeepocean3);
    writevariable(outfile, itsbase1);
    writevariable(outfile, itsbase2);
    writevariable(outfile, itsbase3);
    writevariable(outfile, itsbasetemp1);
    writevariable(outfile, itsbasetemp2);
    writevariable(outfile, itsbasetemp3);
    writevariable(outfile, itshighbase1);
    writevariable(outfile, itshighbase2);
    writevariable(outfile, itshighbase3);
    writevariable(outfile, itsdesert1);
    writevariable(outfile, itsdesert2);
    writevariable(outfile, itsdesert3);
    writevariable(outfile, itshighdesert1);
    writevariable(outfile, itshighdesert2);
    writevariable(outfile, itshighdesert3);
    writevariable(outfile, itscolddesert1);
    writevariable(outfile, itscolddesert2);
    writevariable(outfile, itscolddesert3);
    writevariable(outfile, itsgrass1);
    writevariable(outfile, itsgrass2);
    writevariable(outfile, itsgrass3);
    writevariable(outfile, itscold1);
    writevariable(outfile, itscold2);
    writevariable(outfile, itscold3);
    writevariable(outfile, itstundra1);
    writevariable(outfile, itstundra2);
    writevariable(outfile, itstundra3);
    writevariable(outfile, itseqtundra1);
    writevariable(outfile, itseqtundra2);
    writevariable(outfile, itseqtundra3);
    writevariable(outfile, itssaltpan1);
    writevariable(outfile, itssaltpan2);
    writevariable(outfile, itssaltpan3);
    writevariable(outfile, itserg1);
    writevariable(outfile, itserg2);
    writevariable(outfile, itserg3);
    writevariable(outfile, itswetlands1);
    writevariable(outfile, itswetlands2);
    writevariable(outfile, itswetlands3);
    writevariable(outfile, itslake1);
    writevariable(outfile, itslake2);
    writevariable(outfile, itslake3);
    writevariable(outfile, itsriver1);
    writevariable(outfile, itsriver2);
    writevariable(outfile, itsriver3);
    writevariable(outfile, itsglacier1);
    writevariable(outfile, itsglacier2);
    writevariable(outfile, itsglacier3);
    writevariable(outfile, itssand1);
    writevariable(outfile, itssand2);
    writevariable(outfile, itssand3);
    writevariable(outfile, itsmud1);
    writevariable(outfile, itsmud2);
    writevariable(outfile, itsmud3);
    writevariable(outfile, itsshingle1);
    writevariable(outfile, itsshingle2);
    writevariable(outfile, itsshingle3);
    writevariable(outfile, itsmangrove1);
    writevariable(outfile, itsmangrove2);
    writevariable(outfile, itsmangrove3);
    writevariable(outfile, itshighlight1);
    writevariable(outfile, itshighlight2);
    writevariable(outfile, itshighlight3);
    writevariable(outfile, itsshowmapoutline);
    writevariable(outfile, itsoutline1);
    writevariable(outfile, itsoutline2);
    writevariable(outfile, itsoutline3);
    writevariable(outfile, itselevationlow1);
    writevariable(outfile, itselevationlow2);
    writevariable(outfile, itselevationlow3);
    writevariable(outfile, itselevationhigh1);
    writevariable(outfile, itselevationhigh2);
    writevariable(outfile, itselevationhigh3);
    writevariable(outfile, itstemperaturecold1);
    writevariable(outfile, itstemperaturecold2);
    writevariable(outfile, itstemperaturecold3);
    writevariable(outfile, itstemperaturetemperate1);
    writevariable(outfile, itstemperaturetemperate2);
    writevariable(outfile, itstemperaturetemperate3);
    writevariable(outfile, itstemperaturehot1);
    writevariable(outfile, itstemperaturehot2);
    writevariable(outfile, itstemperaturehot3);
    writevariable(outfile, itsprecipitationdry1);
    writevariable(outfile, itsprecipitationdry2);
    writevariable(outfile, itsprecipitationdry3);
    writevariable(outfile, itsprecipitationwet1);
    writevariable(outfile, itsprecipitationwet2);
    writevariable(outfile, itsprecipitationwet3);
    for (int i = 0; i < CLIMATEMAPSEACOLOURCOUNT; i++)
    {
        for (int j = 0; j < 3; j++)
            writevariable(outfile, itsclimatemapseacolours[i][j]);
    }
    for (int i = 0; i < CLIMATEMAPCOLOURCOUNT; i++)
    {
        for (int j = 0; j < 3; j++)
            writevariable(outfile, itsclimatemapcolours[i][j]);
    }
    for (int i = 0; i < RIVERMAPCOLOURCOUNT; i++)
    {
        for (int j = 0; j < 3; j++)
            writevariable(outfile, itsrivermapcolours[i][j]);
    }
    for (int i = 0; i < RIVERMAPFEATURECOUNT; i++)
        writevariable(outfile, itsshowrivermapfeatures[i]);

    for (int gradient = 0; gradient < MAPGRADIENTTYPECOUNT; gradient++)
    {
        writevariable(outfile, itsmapgradientstopcounts[gradient]);
        writevariable(outfile, itsmapgradientdiscrete[gradient]);

        for (int stop = 0; stop < MAPGRADIENTMAXSTOPS; stop++)
        {
            writevariable(outfile, itsmapgradientpositions[gradient][stop]);

            for (int channel = 0; channel < 3; channel++)
                writevariable(outfile, itsmapgradientcolours[gradient][stop][channel]);
        }
    }

    for (int i = 0; i < BIOMEMAPCOLOURCOUNT; i++)
    {
        for (int j = 0; j < 3; j++)
            writevariable(outfile, itsbiomemapcolours[i][j]);
    }

    writedata(outfile, jantempmap);
    writedata(outfile, jultempmap);
    writedata(outfile, climatemap);
    writedata(outfile, biomemap);
    writedata(outfile, janrainmap);
    writedata(outfile, julrainmap);
    writedata(outfile, janmountainrainmap);
    writedata(outfile, julmountainrainmap);
    writedata(outfile, janmountainraindirmap);
    writedata(outfile, julmountainraindirmap);
    writedata(outfile, seaicemap);
    writedata(outfile, rivermapdir);
    writedata(outfile, rivermapjan);
    writedata(outfile, rivermapjul);
    writedata(outfile, windmap);
    writedata(outfile, lakemap);
    writedata(outfile, roughnessmap);
    writedata(outfile, mountainridges);
    writedata(outfile, mountainheights);
    writedata(outfile, craterrims);
    writedata(outfile, cratercentres);
    writedata(outfile, mapnom);
    writedata(outfile, tidalmap);
    writedata(outfile, riftlakemapsurface);
    writedata(outfile, riftlakemapbed);
    writedata(outfile, lakestartmap);
    writedata(outfile, specials);
    writedata(outfile, geologicregimemap);
    writedata(outfile, tectonicconvergencemap);
    writedata(outfile, tectonicdivergencemap);
    writedata(outfile, tectonicshearmap);
    writedata(outfile, basinclassmap);
    writedata(outfile, erosionpotentialmap);
    writedata(outfile, depositionpotentialmap);
    writedata(outfile, floodplainfertilitymap);
    writedata(outfile, metalorereservemap);
    writedata(outfile, placerreservemap);
    writedata(outfile, evaporitereservemap);
    writedata(outfile, volcanicreservemap);
    writedata(outfile, fisheryreservemap);
    writedata(outfile, extraelevmap);
    writedata(outfile, deltamapdir);
    writedata(outfile, deltamapjan);
    writedata(outfile, deltamapjul);
    writedata(outfile, islandmap);
    writedata(outfile, mountainislandmap);
    writedata(outfile, oceanridgemap);
    writedata(outfile, oceanridgeheightmap);
    writedata(outfile, oceanriftmap);
    writedata(outfile, oceanridgeoffsetmap);
    writedata(outfile, oceanridgeanglemap);
    writedata(outfile, volcanomap);
    writedata(outfile, stratomap);
    writedata(outfile, noshademap);
    writedata(outfile, noisemap);
    writedata(outfile, testmap);
    writedata(outfile, settlement_suitability);
    writedata(outfile, social_infrastructure);
    writedata(outfile, agricultural_capacity);
    writedata(outfile, route_traffic);
    writedata(outfile, river_access);
    writedata(outfile, harbor_score);
    writedata(outfile, owner_settlement_id);
    writedata(outfile, owner_polity_id);

    const auto writestring = [&](const std::string& value)
    {
        const int length = static_cast<int>(value.size());
        writevariable(outfile, length);

        if (length > 0)
            outfile.write(value.data(), static_cast<std::streamsize>(length));
    };

    const int settlementcount = static_cast<int>(settlementlist.size());
    writevariable(outfile, settlementcount);

    for (const Settlement& settlement : settlementlist)
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

    const int politycount = static_cast<int>(politylist.size());
    writevariable(outfile, politycount);

    for (const Polity& polity : politylist)
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

    const int routecount = static_cast<int>(routeedgelist.size());
    writevariable(outfile, routecount);

    for (const RouteEdge& route : routeedgelist)
    {
        writevariable(outfile, route.fromSettlementId);
        writevariable(outfile, route.toSettlementId);
        const std::uint8_t mode = static_cast<std::uint8_t>(route.mode);
        writevariable(outfile, mode);
        writevariable(outfile, route.cost);
        writevariable(outfile, route.capacity);
        writevariable(outfile, route.traffic);
    }

    const int eventcount = static_cast<int>(historyeventlist.size());
    writevariable(outfile, eventcount);

    for (const HistoryEvent& event : historyeventlist)
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
        writeshortvectordata(outfile, seasonaltempmaps[season]);
        writeshortvectordata(outfile, seasonalrainmaps[season]);
        writeshortvectordata(outfile, seasonalpressuremaps[season]);
        writeshortvectordata(outfile, seasonaluwindmaps[season]);
        writeshortvectordata(outfile, seasonalvwindmaps[season]);
        writeshortvectordata(outfile, seasonalupperheightmaps[season]);
        writeshortvectordata(outfile, seasonalupperuwindmaps[season]);
        writeshortvectordata(outfile, seasonaluppervwindmaps[season]);
        writeshortvectordata(outfile, seasonalverticalvelocitymaps[season]);
        writeshortvectordata(outfile, seasonalcurrentumaps[season]);
        writeshortvectordata(outfile, seasonalcurrentvmaps[season]);
        writeshortvectordata(outfile, seasonalsstmaps[season]);
        writeshortvectordata(outfile, seasonalevaporationmaps[season]);
        writeshortvectordata(outfile, seasonalmaritimeinfluencemaps[season]);
        writeshortvectordata(outfile, seasonalmaritimethermalanomalymaps[season]);
        writeshortvectordata(outfile, seasonalmaritimefetchmaps[season]);
        writeshortvectordata(outfile, seasonalmoisturemaps[season]);
        writeshortvectordata(outfile, seasonalconvergencemaps[season]);
        writeshortvectordata(outfile, seasonalupliftmaps[season]);
        writeshortvectordata(outfile, seasonalsubsidencemaps[season]);
    }

    for (int i = 0; i < ARRAYWIDTH; i++)
    {
        for (int j = 0; j < 6; j++)
            writevariable(outfile, horselats[i][j]);
    }

    for (int i = 0; i < itscraterno; i++)
    {
        writevariable(outfile, cratercentreslist[i].w);
        writevariable(outfile, cratercentreslist[i].x);
        writevariable(outfile, cratercentreslist[i].y);
        writevariable(outfile, cratercentreslist[i].z);
    }

    writevariable(outfile, itstectonictimeoriginstep);
    writevariable(outfile, itstectonictimemyr);
    writevariable(outfile, itstectonicdeltatimemyr);
    writedata(outfile, tectoniccrustagemyrmap);
    writedata(outfile, tectoniccrustthicknessmap);
    writedata(outfile, tectoniccrustclassmap);
    writedata(outfile, tectonicuplifttendencymap);
    writedata(outfile, tectonicsubsidencetendencymap);
    writedata(outfile, tectonicaccumulatedstrainmap);
    writedata(outfile, tectonicboundarytypemap);
    writedata(outfile, tectonicboundarydistancemap);
    writedata(outfile, tectonicboundarysegmentidmap);
    writedata(outfile, tectonicnearestboundaryidmap);
    writedata(outfile, tectonicboundaryhistorymap);
    writedata(outfile, tectonicdeformingregionidmap);
    writedata(outfile, tectonicdeformingregiontypemap);
    writedata(outfile, tectonicdeformationratemap);
    writedata(outfile, tectonicdeformationvelocityxmap);
    writedata(outfile, tectonicdeformationvelocityymap);
    writevariable(outfile, itstectoniccyclecount);
    writevariable(outfile, itstectonicplatecount);
    writevariable(outfile, itstectonicsealevelm);

    int boundarysegmentcount = static_cast<int>(tectonicboundarysegmentlist.size());
    writevariable(outfile, boundarysegmentcount);
    for (const TectonicBoundarySegment& segment : tectonicboundarysegmentlist)
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

    int deformingregioncount = static_cast<int>(tectonicdeformingregionlist.size());
    writevariable(outfile, deformingregioncount);
    for (const TectonicDeformingRegion& region : tectonicdeformingregionlist)
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

    if (!outfile.good())
    {
        cerr << "Error writing world '" << filename << "'" << endl;
    }
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

    if (fileversion != itssaveversion) // Incompatible file format!
        return 0;

    readvariable(infile, itssize);
    readvariable(infile, itswidth);
    readvariable(infile, itsheight);
    readvariable(infile, itsseed);
    resizeseasonalclimatefields();
    readvariable(infile, itsrotation);
    readvariable(infile, itstilt);
    readvariable(infile, itseccentricity);
    readvariable(infile, itsperihelion);
    readvariable(infile, itsgravity);
    readvariable(infile, itslunar);
    readvariable(infile, itstempdecrease);
    readvariable(infile, itsnorthpolaradjust);
    readvariable(infile, itssouthpolaradjust);
    readvariable(infile, itsaveragetemp);
    readvariable(infile, itsnorthpolartemp);
    readvariable(infile, itssouthpolartemp);
    readvariable(infile, itseqtemp);
    readvariable(infile, itswaterpickup);
    readvariable(infile, itsriverfactor);
    readvariable(infile, itsriverlandreduce);
    readvariable(infile, itsestuarylimit);
    readvariable(infile, itsglacialtemp);
    readvariable(infile, itsglaciertemp);
    readvariable(infile, itsmountainreduce);
    readvariable(infile, itsclimateno);
    readvariable(infile, itsmaxheight);
    readvariable(infile, itssealevel);
    readvariable(infile, itslandtotal);
    readvariable(infile, itsseatotal);
    readvariable(infile, itscraterno);

    readvariable(infile, itslandshading);
    readvariable(infile, itslakeshading);
    readvariable(infile, itsseashading);
    readvariable(infile, itsshadingdir);
    readvariable(infile, itssnowchange);
    readvariable(infile, itsseaiceappearance);
    readvariable(infile, itslandmarbling);
    readvariable(infile, itslakemarbling);
    readvariable(infile, itsseamarbling);
    readvariable(infile, itsminriverflowglobal);
    readvariable(infile, itsminriverflowregional);
    readvariable(infile, itsmangroves);
    readvariable(infile, itscolourcliffs);
    readvariable(infile, itsseaice1);
    readvariable(infile, itsseaice2);
    readvariable(infile, itsseaice3);
    readvariable(infile, itsocean1);
    readvariable(infile, itsocean2);
    readvariable(infile, itsocean3);
    readvariable(infile, itsdeepocean1);
    readvariable(infile, itsdeepocean2);
    readvariable(infile, itsdeepocean3);
    readvariable(infile, itsbase1);
    readvariable(infile, itsbase2);
    readvariable(infile, itsbase3);
    readvariable(infile, itsbasetemp1);
    readvariable(infile, itsbasetemp2);
    readvariable(infile, itsbasetemp3);
    readvariable(infile, itshighbase1);
    readvariable(infile, itshighbase2);
    readvariable(infile, itshighbase3);
    readvariable(infile, itsdesert1);
    readvariable(infile, itsdesert2);
    readvariable(infile, itsdesert3);
    readvariable(infile, itshighdesert1);
    readvariable(infile, itshighdesert2);
    readvariable(infile, itshighdesert3);
    readvariable(infile, itscolddesert1);
    readvariable(infile, itscolddesert2);
    readvariable(infile, itscolddesert3);
    readvariable(infile, itsgrass1);
    readvariable(infile, itsgrass2);
    readvariable(infile, itsgrass3);
    readvariable(infile, itscold1);
    readvariable(infile, itscold2);
    readvariable(infile, itscold3);
    readvariable(infile, itstundra1);
    readvariable(infile, itstundra2);
    readvariable(infile, itstundra3);
    readvariable(infile, itseqtundra1);
    readvariable(infile, itseqtundra2);
    readvariable(infile, itseqtundra3);
    readvariable(infile, itssaltpan1);
    readvariable(infile, itssaltpan2);
    readvariable(infile, itssaltpan3);
    readvariable(infile, itserg1);
    readvariable(infile, itserg2);
    readvariable(infile, itserg3);
    readvariable(infile, itswetlands1);
    readvariable(infile, itswetlands2);
    readvariable(infile, itswetlands3);
    readvariable(infile, itslake1);
    readvariable(infile, itslake2);
    readvariable(infile, itslake3);
    readvariable(infile, itsriver1);
    readvariable(infile, itsriver2);
    readvariable(infile, itsriver3);
    readvariable(infile, itsglacier1);
    readvariable(infile, itsglacier2);
    readvariable(infile, itsglacier3);
    readvariable(infile, itssand1);
    readvariable(infile, itssand2);
    readvariable(infile, itssand3);
    readvariable(infile, itsmud1);
    readvariable(infile, itsmud2);
    readvariable(infile, itsmud3);
    readvariable(infile, itsshingle1);
    readvariable(infile, itsshingle2);
    readvariable(infile, itsshingle3);
    readvariable(infile, itsmangrove1);
    readvariable(infile, itsmangrove2);
    readvariable(infile, itsmangrove3);
    readvariable(infile, itshighlight1);
    readvariable(infile, itshighlight2);
    readvariable(infile, itshighlight3);

    if (fileversion >= 2)
    {
        readvariable(infile, itsshowmapoutline);
        readvariable(infile, itsoutline1);
        readvariable(infile, itsoutline2);
        readvariable(infile, itsoutline3);
        readvariable(infile, itselevationlow1);
        readvariable(infile, itselevationlow2);
        readvariable(infile, itselevationlow3);
        readvariable(infile, itselevationhigh1);
        readvariable(infile, itselevationhigh2);
        readvariable(infile, itselevationhigh3);
        readvariable(infile, itstemperaturecold1);
        readvariable(infile, itstemperaturecold2);
        readvariable(infile, itstemperaturecold3);
        readvariable(infile, itstemperaturetemperate1);
        readvariable(infile, itstemperaturetemperate2);
        readvariable(infile, itstemperaturetemperate3);
        readvariable(infile, itstemperaturehot1);
        readvariable(infile, itstemperaturehot2);
        readvariable(infile, itstemperaturehot3);
        readvariable(infile, itsprecipitationdry1);
        readvariable(infile, itsprecipitationdry2);
        readvariable(infile, itsprecipitationdry3);
        readvariable(infile, itsprecipitationwet1);
        readvariable(infile, itsprecipitationwet2);
        readvariable(infile, itsprecipitationwet3);

        if (fileversion >= 3)
        {
            for (int i = 0; i < CLIMATEMAPSEACOLOURCOUNT; i++)
            {
                for (int j = 0; j < 3; j++)
                    readvariable(infile, itsclimatemapseacolours[i][j]);
            }

            for (int i = 0; i < CLIMATEMAPCOLOURCOUNT; i++)
            {
                for (int j = 0; j < 3; j++)
                    readvariable(infile, itsclimatemapcolours[i][j]);
            }

            for (int i = 0; i < RIVERMAPCOLOURCOUNT; i++)
            {
                for (int j = 0; j < 3; j++)
                    readvariable(infile, itsrivermapcolours[i][j]);
            }

            for (int i = 0; i < RIVERMAPFEATURECOUNT; i++)
                readvariable(infile, itsshowrivermapfeatures[i]);

            if (fileversion >= 4)
            {
                for (int gradient = 0; gradient < MAPGRADIENTTYPECOUNT; gradient++)
                {
                    readvariable(infile, itsmapgradientstopcounts[gradient]);
                    readvariable(infile, itsmapgradientdiscrete[gradient]);

                    for (int stop = 0; stop < MAPGRADIENTMAXSTOPS; stop++)
                    {
                        readvariable(infile, itsmapgradientpositions[gradient][stop]);

                        for (int channel = 0; channel < 3; channel++)
                            readvariable(infile, itsmapgradientcolours[gradient][stop][channel]);
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
                        readvariable(infile, itsbiomemapcolours[i][j]);
                }
            }
            else
            {
                for (int i = 0; i < BIOMEMAPCOLOURCOUNT; i++)
                {
                    for (int j = 0; j < 3; j++)
                        itsbiomemapcolours[i][j] = defaultbiomemapcolours[i][j];
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
                    itsclimatemapseacolours[i][j] = defaultclimateseacolours[i][j];
            }

            for (int i = 0; i < CLIMATEMAPCOLOURCOUNT; i++)
            {
                for (int j = 0; j < 3; j++)
                    itsclimatemapcolours[i][j] = defaultclimatecolours[i][j];
            }

            itsrivermapcolours[rivermapbackground] = { 255, 255, 255 };
            itsrivermapcolours[rivermaplowflow] = { 255, 255, 255 };
            itsrivermapcolours[rivermaphighflow] = { itsriver1, itsriver2, itsriver3 };
            itsrivermapcolours[rivermaplake] = { itslake1, itslake2, itslake3 };
            itsrivermapcolours[rivermapsaltpan] = { itssaltpan1, itssaltpan2, itssaltpan3 };
            itsrivermapcolours[rivermapwetlands] = { itswetlands1, itswetlands2, itswetlands3 };
            itsrivermapcolours[rivermapmud] = { itsmud1, itsmud2, itsmud3 };
            itsrivermapcolours[rivermapsand] = { itssand1, itssand2, itssand3 };
            itsrivermapcolours[rivermapshingle] = { itsshingle1, itsshingle2, itsshingle3 };
            itsrivermapcolours[rivermapvolcano] = { 240, 0, 0 };

            for (int i = 0; i < RIVERMAPFEATURECOUNT; i++)
                itsshowrivermapfeatures[i] = true;

            initialisegradientmapappearance(*this);

            for (int i = 0; i < BIOMEMAPCOLOURCOUNT; i++)
            {
                for (int j = 0; j < 3; j++)
                    itsbiomemapcolours[i][j] = defaultbiomemapcolours[i][j];
            }
        }
    }
    else
        setdefaultnonreliefmapappearance(*this);

    readdata(infile, jantempmap);
    readdata(infile, jultempmap);
    readdata(infile, climatemap);
    readdata(infile, biomemap);
    readdata(infile, janrainmap);
    readdata(infile, julrainmap);
    readdata(infile, janmountainrainmap);
    readdata(infile, julmountainrainmap);
    readdata(infile, janmountainraindirmap);
    readdata(infile, julmountainraindirmap);
    readdata(infile, seaicemap);
    readdata(infile, rivermapdir);
    readdata(infile, rivermapjan);
    readdata(infile, rivermapjul);
    readdata(infile, windmap);
    readdata(infile, lakemap);
    readdata(infile, roughnessmap);
    readdata(infile, mountainridges);
    readdata(infile, mountainheights);
    readdata(infile, craterrims);
    readdata(infile, cratercentres);
    readdata(infile, mapnom);
    readdata(infile, tidalmap);
    readdata(infile, riftlakemapsurface);
    readdata(infile, riftlakemapbed);
    readdata(infile, lakestartmap);
    readdata(infile, specials);
    if (fileversion >= 12)
    {
        readdata(infile, geologicregimemap);
        readdata(infile, tectonicconvergencemap);
        readdata(infile, tectonicdivergencemap);
        readdata(infile, tectonicshearmap);
        readdata(infile, basinclassmap);
        readdata(infile, erosionpotentialmap);
        readdata(infile, depositionpotentialmap);
        readdata(infile, floodplainfertilitymap);
        readdata(infile, metalorereservemap);
        readdata(infile, placerreservemap);
        readdata(infile, evaporitereservemap);
        readdata(infile, volcanicreservemap);
        readdata(infile, fisheryreservemap);
    }
    else
    {
        cleartectonicprovenanceinternal();

        parallelforrows(0, ARRAYWIDTH - 1, [&](int startx, int endx)
        {
            for (int i = startx; i <= endx; i++)
            {
                for (int j = 0; j < ARRAYHEIGHT; j++)
                {
                    basinclassmap[i][j] = static_cast<std::uint8_t>(BasinClass::none);
                    erosionpotentialmap[i][j] = 0;
                    depositionpotentialmap[i][j] = 0;
                    floodplainfertilitymap[i][j] = 0;
                    metalorereservemap[i][j] = 0;
                    placerreservemap[i][j] = 0;
                    evaporitereservemap[i][j] = 0;
                    volcanicreservemap[i][j] = 0;
                    fisheryreservemap[i][j] = 0;
                }
            }
        }, 64);
    }
    readdata(infile, extraelevmap);
    readdata(infile, deltamapdir);
    readdata(infile, deltamapjan);
    readdata(infile, deltamapjul);
    readdata(infile, islandmap);
    readdata(infile, mountainislandmap);
    readdata(infile, oceanridgemap);
    readdata(infile, oceanridgeheightmap);
    readdata(infile, oceanriftmap);
    readdata(infile, oceanridgeoffsetmap);
    readdata(infile, oceanridgeanglemap);
    readdata(infile, volcanomap);
    readdata(infile, stratomap);
    readdata(infile, noshademap);
    readdata(infile, noisemap);
    readdata(infile, testmap);

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
        readdata(infile, settlement_suitability);
        readdata(infile, social_infrastructure);
        readdata(infile, agricultural_capacity);
        readdata(infile, route_traffic);
        readdata(infile, river_access);
        readdata(infile, harbor_score);
        readdata(infile, owner_settlement_id);
        readdata(infile, owner_polity_id);

        clearsocialstate();

        int settlementcount = 0;
        readvariable(infile, settlementcount);
        settlementcount = std::max(0, settlementcount);
        settlementlist.resize(static_cast<size_t>(settlementcount));

        for (Settlement& settlement : settlementlist)
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
        politylist.resize(static_cast<size_t>(politycount));

        for (Polity& polity : politylist)
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
        routeedgelist.resize(static_cast<size_t>(routecount));

        for (RouteEdge& route : routeedgelist)
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
        historyeventlist.resize(static_cast<size_t>(eventcount));

        for (HistoryEvent& event : historyeventlist)
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
        parallelforrows(0, ARRAYWIDTH - 1, [&](int startx, int endx)
        {
            for (int i = startx; i <= endx; i++)
            {
                for (int j = 0; j < ARRAYHEIGHT; j++)
                {
                    settlement_suitability[i][j] = 0;
                    social_infrastructure[i][j] = 0;
                    agricultural_capacity[i][j] = 0;
                    route_traffic[i][j] = 0;
                    river_access[i][j] = 0;
                    harbor_score[i][j] = 0;
                    owner_settlement_id[i][j] = -1;
                    owner_polity_id[i][j] = -1;
                }
            }
        }, 64);

        clearsocialstate();
    }

    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
    {
        readshortvectordata(infile, seasonaltempmaps[season]);
        readshortvectordata(infile, seasonalrainmaps[season]);
        std::transform(
            seasonalrainmaps[season].begin(),
            seasonalrainmaps[season].end(),
            seasonalrainfloatmaps[season].begin(),
            [](short value) { return static_cast<float>(value); });
        readshortvectordata(infile, seasonalpressuremaps[season]);
        readshortvectordata(infile, seasonaluwindmaps[season]);
        readshortvectordata(infile, seasonalvwindmaps[season]);
        if (fileversion >= 13)
        {
            readshortvectordata(infile, seasonalupperheightmaps[season]);
            readshortvectordata(infile, seasonalupperuwindmaps[season]);
            readshortvectordata(infile, seasonaluppervwindmaps[season]);
            readshortvectordata(infile, seasonalverticalvelocitymaps[season]);
        }
        else
        {
            std::fill(seasonalupperheightmaps[season].begin(), seasonalupperheightmaps[season].end(), 0);
            std::fill(seasonalupperuwindmaps[season].begin(), seasonalupperuwindmaps[season].end(), 0);
            std::fill(seasonaluppervwindmaps[season].begin(), seasonaluppervwindmaps[season].end(), 0);
            std::fill(seasonalverticalvelocitymaps[season].begin(), seasonalverticalvelocitymaps[season].end(), 0);
        }
        readshortvectordata(infile, seasonalcurrentumaps[season]);
        readshortvectordata(infile, seasonalcurrentvmaps[season]);
        readshortvectordata(infile, seasonalsstmaps[season]);
        readshortvectordata(infile, seasonalevaporationmaps[season]);
        readshortvectordata(infile, seasonalmaritimeinfluencemaps[season]);
        readshortvectordata(infile, seasonalmaritimethermalanomalymaps[season]);
        readshortvectordata(infile, seasonalmaritimefetchmaps[season]);
        readshortvectordata(infile, seasonalmoisturemaps[season]);
        readshortvectordata(infile, seasonalconvergencemaps[season]);
        readshortvectordata(infile, seasonalupliftmaps[season]);
        readshortvectordata(infile, seasonalsubsidencemaps[season]);
    }

    for (int i = 0; i < ARRAYWIDTH; i++)
    {
        for (int j = 0; j < 6; j++)
            readvariable(infile, horselats[i][j]);
    }

    for (int i = 0; i < itscraterno; i++)
    {
        readvariable(infile, cratercentreslist[i].w);
        readvariable(infile, cratercentreslist[i].x);
        readvariable(infile, cratercentreslist[i].y);
        readvariable(infile, cratercentreslist[i].z);
    }

    if (fileversion >= 15)
    {
        readvariable(infile, itstectonictimeoriginstep);
        readvariable(infile, itstectonictimemyr);
        readvariable(infile, itstectonicdeltatimemyr);
        readdata(infile, tectoniccrustagemyrmap);
        readdata(infile, tectoniccrustthicknessmap);
        readdata(infile, tectoniccrustclassmap);
        readdata(infile, tectonicuplifttendencymap);
        readdata(infile, tectonicsubsidencetendencymap);
        readdata(infile, tectonicaccumulatedstrainmap);
        readdata(infile, tectonicboundarytypemap);
        readdata(infile, tectonicboundarydistancemap);
        readdata(infile, tectonicboundarysegmentidmap);
        readdata(infile, tectonicnearestboundaryidmap);
        readdata(infile, tectonicboundaryhistorymap);
        readdata(infile, tectonicdeformingregionidmap);
        readdata(infile, tectonicdeformingregiontypemap);
        readdata(infile, tectonicdeformationratemap);
        readdata(infile, tectonicdeformationvelocityxmap);
        readdata(infile, tectonicdeformationvelocityymap);

        if (fileversion >= 16)
        {
            readvariable(infile, itstectoniccyclecount);
            readvariable(infile, itstectonicplatecount);
            readvariable(infile, itstectonicsealevelm);

            int boundarysegmentcount = 0;
            readvariable(infile, boundarysegmentcount);
            boundarysegmentcount = std::max(0, boundarysegmentcount);
            tectonicboundarysegmentlist.resize(static_cast<size_t>(boundarysegmentcount));
            for (TectonicBoundarySegment& segment : tectonicboundarysegmentlist)
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
            tectonicdeformingregionlist.resize(static_cast<size_t>(deformingregioncount));
            for (TectonicDeformingRegion& region : tectonicdeformingregionlist)
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
            itstectoniccyclecount = 0;
            itstectonicplatecount = 0;
            itstectonicsealevelm = 0;
            tectonicboundarysegmentlist.clear();
            tectonicdeformingregionlist.clear();
        }
    }
    else
    {
        itstectonictimeoriginstep = 0;
        itstectonictimemyr = 0.0f;
        itstectonicdeltatimemyr = 0.0f;
        itstectoniccyclecount = 0;
        itstectonicplatecount = 0;
        itstectonicsealevelm = 0;
        tectonicboundarysegmentlist.clear();
        tectonicdeformingregionlist.clear();

        parallelforrows(0, ARRAYWIDTH - 1, [&](int startx, int endx)
        {
            for (int i = startx; i <= endx; i++)
            {
                for (int j = 0; j < ARRAYHEIGHT; j++)
                {
                    tectoniccrustagemyrmap[i][j] = 0.0f;
                    tectoniccrustthicknessmap[i][j] = 0.0f;
                    tectoniccrustclassmap[i][j] = static_cast<std::uint8_t>(CrustClass::none);
                    tectonicuplifttendencymap[i][j] = 0.0f;
                    tectonicsubsidencetendencymap[i][j] = 0.0f;
                    tectonicaccumulatedstrainmap[i][j] = 0.0f;
                    tectonicboundarytypemap[i][j] = static_cast<std::uint8_t>(BoundaryType::none);
                    tectonicboundarydistancemap[i][j] = 0;
                    tectonicboundarysegmentidmap[i][j] = 0;
                    tectonicnearestboundaryidmap[i][j] = 0;
                    tectonicboundaryhistorymap[i][j] = 0.0f;
                    tectonicdeformingregionidmap[i][j] = 0;
                    tectonicdeformingregiontypemap[i][j] = static_cast<std::uint8_t>(DeformingRegionType::none);
                    tectonicdeformationratemap[i][j] = 0.0f;
                    tectonicdeformationvelocityxmap[i][j] = 0.0f;
                    tectonicdeformationvelocityymap[i][j] = 0.0f;
                }
            }
        }, 64);
    }

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
    while (x > itswidth) // If it's too large, wrap it.
    {
        x = x - itswidth;
    }

    while (x < 0) // If it's too small, wrap it.
    {
        x = x + itswidth;
    }

    return(x);
}

int planet::clipy(int y) const // This clips Y coordinates so they can't be off the map.
{
    if (y < 0)
        y = 0;

    if (y > itsheight)
        y = itsheight;

    return(y);
}

void planet::smooth(int arr[][ARRAYHEIGHT], int amount, bool vary, bool avoidmountains) // This smoothes the given array by the given amount.
{
    for (int i = 0; i <= itswidth; i++)
    {
        for (int j = 1; j < itsheight; j++)
        {
            if (avoidmountains == 0 || mountainheights[i][j] == 0)
            {
                int crount = 0;
                int ave = 0;

                for (int k = i - amount; k <= i + amount; k++)
                {
                    int kk = k;

                    if (kk < 0)
                        kk = itswidth;

                    if (kk > itswidth)
                        kk = 0;

                    for (int l = j - amount; l <= j + amount; l++)
                    {
                        ave = ave + mapnom[kk][l];
                        crount++;
                    }
                }
                ave = ave / crount;

                if (ave > 0 && ave < itsmaxheight)
                    mapnom[i][j] = ave;
            }
        }
    }
}

void planet::smooth(short arr[][ARRAYHEIGHT], int amount, bool vary, bool avoidmountains) // This smoothes the given array by the given amount.
{
    for (int i = 0; i <= itswidth; i++)
    {
        for (int j = 1; j < itsheight; j++)
        {
            if (avoidmountains == 0 || mountainheights[i][j] == 0)
            {
                int crount = 0;
                int ave = 0;

                for (int k = i - amount; k <= i + amount; k++)
                {
                    int kk = k;

                    if (kk < 0)
                        kk = itswidth;

                    if (kk > itswidth)
                        kk = 0;

                    for (int l = j - amount; l <= j + amount; l++)
                    {
                        ave = ave + mapnom[kk][l];
                        crount++;
                    }
                }
                ave = ave / crount;

                if (ave > 0 && ave < itsmaxheight)
                    mapnom[i][j] = ave;
            }
        }
    }
}

// This does the same, but only over land.

void planet::smoothoverland(int arr[][ARRAYHEIGHT], int amount, bool uponly)
{
    for (int i = 0; i <= itswidth; i++)
    {
        for (int j = 1; j < itsheight; j++)
        {
            if (sea(i, j) == 0)
            {
                int crount = 0;
                int ave = 0;

                for (int k = i - amount; k <= i + amount; k++)
                {
                    int kk = k;

                    if (kk < 0)
                        kk = itswidth;

                    if (kk > itswidth)
                        kk = 0;

                    for (int l = j - amount; l <= j + amount; l++)
                    {
                        ave = ave + arr[kk][l];
                        crount++;
                    }
                }

                if (crount > 0)
                {
                    ave = ave / crount;

                    if (ave > 0 && ave < itsmaxheight)
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

void planet::smoothoverland(short arr[][ARRAYHEIGHT], int amount, bool uponly)
{
    for (int i = 0; i <= itswidth; i++)
    {
        for (int j = 1; j < itsheight; j++)
        {
            if (sea(i, j) == 0)
            {
                int crount = 0;
                int ave = 0;

                for (int k = i - amount; k <= i + amount; k++)
                {
                    int kk = k;

                    if (kk < 0)
                        kk = itswidth;

                    if (kk > itswidth)
                        kk = 0;

                    for (int l = j - amount; l <= j + amount; l++)
                    {
                        ave = ave + (int)arr[kk][l];
                        crount++;
                    }
                }

                if (crount > 0)
                {
                    ave = ave / crount;

                    if (ave > 0 && ave < itsmaxheight)
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

template<typename T> void planet::shift(T arr[][ARRAYHEIGHT], int offset)
{
    vector<vector<T>> dummy(ARRAYWIDTH, vector<T>(ARRAYHEIGHT, 0));

    //T dummy[ARRAYWIDTH][ARRAYHEIGHT];

    for (int i = 0; i <= itswidth; i++)
    {
        for (int j = 0; j <= itsheight; j++)
            dummy[i][j] = arr[i][j];
    }

    for (int i = 0; i <= itswidth; i++)
    {
        int ii = i + offset;

        if (ii<0 || ii>itswidth)
            ii = wrap(ii, itswidth);

        for (int j = 0; j <= itsheight; j++)
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

template<typename T> void planet::writedata(ofstream& outfile, T const arr[ARRAYWIDTH][ARRAYHEIGHT])
{
    for (int i = 0; i <= itswidth; i++)
    {
        for (int j = 0; j <= itsheight; j++)
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

template<typename T> void planet::readdata(ifstream& infile, T arr[ARRAYWIDTH][ARRAYHEIGHT])
{
    for (int i = 0; i <= itswidth; i++)
    {
        for (int j = 0; j <= itsheight; j++)
        {
            getline(infile, line_for_file_read);
            read_val(line_for_file_read, arr[i][j]);
        }
    }
}
