#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

#include "climate_atmosphere.hpp"
#include "climate_physics.hpp"
#include "generation_tuning.hpp"
#include "planet.hpp"
#include "functions.hpp"

using namespace std;

namespace
{
using floatgrid = vector<vector<float>>;

constexpr std::array<float, CLIMATESEASONCOUNT> seasonlatitudephase = { -1.0f, 0.0f, 1.0f, 0.0f };

int wrapx(int x, int width)
{
    if (x < 0 || x > width)
        return wrap(x, width);

    return x;
}

float latitudeforrow(int y, int height)
{
    if (height <= 0)
        return 0.0f;

    return 90.0f - (180.0f * static_cast<float>(y) / static_cast<float>(height));
}

float coastalweight(int distance, int maxdistance)
{
    if (distance <= 0 || distance > maxdistance)
        return 0.0f;

    return 1.0f - (static_cast<float>(distance - 1) / static_cast<float>(maxdistance));
}

bool landindir(planet& world, int x, int y, int dx, int dy, int maxdistance, int& nearestdistance)
{
    const int width = world.width();
    const int height = world.height();
    const int sealevel = world.sealevel();

    nearestdistance = maxdistance + 1;

    for (int distance = 1; distance <= maxdistance; distance++)
    {
        const int xx = wrapx(x + dx * distance, width);
        const int yy = y + dy * distance;

        if (yy < 0 || yy > height)
            break;

        if (world.nom(xx, yy) > sealevel)
        {
            nearestdistance = distance;
            return true;
        }
    }

    return false;
}

float sampleoceanfield(const floatgrid& field, planet& world, int x, int y)
{
    const int width = world.width();
    const int height = world.height();
    const int sealevel = world.sealevel();

    if (y >= 0 && y <= height && world.nom(x, y) <= sealevel)
        return field[x][y];

    for (int radius = 1; radius <= 2; radius++)
    {
        for (int dy = -radius; dy <= radius; dy++)
        {
            const int yy = y + dy;

            if (yy < 0 || yy > height)
                continue;

            for (int dx = -radius; dx <= radius; dx++)
            {
                const int xx = wrapx(x + dx, width);

                if (world.nom(xx, yy) <= sealevel)
                    return field[xx][yy];
            }
        }
    }

    return field[wrapx(x, width)][std::clamp(y, 0, height)];
}

void smoothallfield(planet& world, floatgrid& field, int iterations);

float samplewrappedfield(const floatgrid& field, planet& world, float x, float y)
{
    const int width = world.width();
    const int height = world.height();
    const float span = static_cast<float>(width + 1);

    while (x < 0.0f)
        x += span;

    while (x >= span)
        x -= span;

    y = std::clamp(y, 0.0f, static_cast<float>(height));

    const int x0 = static_cast<int>(std::floor(x));
    const int y0 = static_cast<int>(std::floor(y));
    const int x1 = wrapx(x0 + 1, width);
    const int y1 = std::min(y0 + 1, height);
    const float fracx = x - static_cast<float>(x0);
    const float fracy = y - static_cast<float>(y0);
    const float v00 = field[x0][y0];
    const float v10 = field[x1][y0];
    const float v01 = field[x0][y1];
    const float v11 = field[x1][y1];

    return
        v00 * (1.0f - fracx) * (1.0f - fracy) +
        v10 * fracx * (1.0f - fracy) +
        v01 * (1.0f - fracx) * fracy +
        v11 * fracx * fracy;
}

void smoothseasonalfield(planet& world, floatgrid& field, int iterations)
{
    const int width = world.width();
    const int height = world.height();
    const int sealevel = world.sealevel();

    floatgrid scratch = field;

    for (int iteration = 0; iteration < iterations; iteration++)
    {
        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int x = 0; x <= width; x++)
                {
                    if (world.nom(x, y) > sealevel)
                    {
                        scratch[x][y] = 0.0f;
                        continue;
                    }

                    float total = 0.0f;
                    float weighttotal = 0.0f;

                    for (int dy = -1; dy <= 1; dy++)
                    {
                        const int yy = y + dy;

                        if (yy < 0 || yy > height)
                            continue;

                        for (int dx = -1; dx <= 1; dx++)
                        {
                            const int xx = wrapx(x + dx, width);

                            if (world.nom(xx, yy) > sealevel)
                                continue;

                            const float weight = (dx == 0 && dy == 0) ? 2.0f : 1.0f;
                            total += field[xx][yy] * weight;
                            weighttotal += weight;
                        }
                    }

                    scratch[x][y] = (weighttotal > 0.0f) ? total / weighttotal : field[x][y];
                }
            }
        });

        field.swap(scratch);
    }
}

void applytopographicwindeffects(planet& world, const floatgrid& macroterrain, floatgrid& windu, floatgrid& windv)
{
    const int width = world.width();
    const int height = world.height();
    const float maxvectorwind = tuning::climate::atmosphere::maxVectorWind;

    auto computegradient = [&](int x, int y)
    {
        const int xwest = wrapx(x - 1, width);
        const int xeast = wrapx(x + 1, width);
        const int ynorth = (y > 0) ? y - 1 : y;
        const int ysouth = (y < height) ? y + 1 : y;

        const float gradx = (macroterrain[xeast][y] - macroterrain[xwest][y]) / 2.0f;
        const float grady = (macroterrain[x][ysouth] - macroterrain[x][ynorth]) / 2.0f;

        return std::pair<float, float>(gradx, grady);
    };

    for (int iteration = 0; iteration < tuning::climate::atmosphere::topographyIterations; iteration++)
    {
        floatgrid nextu = windu;
        floatgrid nextv = windv;

        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int x = 0; x <= width; x++)
                {
                    const float u = windu[x][y];
                    const float v = windv[x][y];
                    const float speed = std::sqrt(u * u + v * v);

                    if (speed < tuning::climate::atmosphere::topographyMinimumWindSpeed)
                        continue;

                    const auto [gradx, grady] = computegradient(x, y);
                    const float gradmag = std::sqrt(gradx * gradx + grady * grady);
                    const float terrainhere = macroterrain[x][y];

                    if (gradmag < tuning::climate::atmosphere::topographyMinimumRelief && terrainhere < tuning::climate::atmosphere::topographyMinimumRelief)
                        continue;

                    const float nx = gradx / std::max(gradmag, 0.0001f);
                    const float ny = grady / std::max(gradmag, 0.0001f);
                    const float tx = -ny;
                    const float ty = nx;
                    const float dirx = u / speed;
                    const float diry = v / speed;

                    float lookaheadrise = 0.0f;

                    for (int step = 1; step <= tuning::climate::atmosphere::topographyLookaheadDistance; step++)
                    {
                        const float samplex = static_cast<float>(x) + dirx * static_cast<float>(step);
                        const float sampley = static_cast<float>(y) + diry * static_cast<float>(step);
                        const float rise = samplewrappedfield(macroterrain, world, samplex, sampley) - terrainhere;

                        if (rise > lookaheadrise)
                            lookaheadrise = rise;
                    }

                    const float positivecontour = samplewrappedfield(macroterrain, world,
                        static_cast<float>(x) + tx * tuning::climate::atmosphere::topographySideSampleDistance,
                        static_cast<float>(y) + ty * tuning::climate::atmosphere::topographySideSampleDistance);
                    const float negativecontour = samplewrappedfield(macroterrain, world,
                        static_cast<float>(x) - tx * tuning::climate::atmosphere::topographySideSampleDistance,
                        static_cast<float>(y) - ty * tuning::climate::atmosphere::topographySideSampleDistance);

                    const float steeringdirection = (positivecontour <= negativecontour) ? 1.0f : -1.0f;
                    const float gradientbarrier = std::clamp(gradmag / tuning::climate::atmosphere::topographyGradientScale, 0.0f, 1.0f);
                    const float lookaheadbarrier = std::clamp(lookaheadrise / tuning::climate::atmosphere::topographyLookaheadRiseScale, 0.0f, 1.0f);
                    const float barrier = std::max(gradientbarrier, lookaheadbarrier);

                    if (barrier <= 0.0f)
                        continue;

                    float crossridge = u * nx + v * ny;
                    float alongridge = u * tx + v * ty;

                    if (crossridge > 0.0f)
                    {
                        const float blockedfraction = barrier * (1.0f - tuning::climate::atmosphere::blockedComponentFactor);
                        const float blockedcross = crossridge * blockedfraction;

                        crossridge = crossridge - blockedcross;
                        alongridge = alongridge + blockedcross * tuning::climate::atmosphere::topographyDeflectionFactor * steeringdirection;
                        crossridge = crossridge * (1.0f - barrier * tuning::climate::atmosphere::topographyChannelFactor);
                    }
                    else if (crossridge < 0.0f)
                    {
                        crossridge = crossridge * (1.0f + barrier * tuning::climate::atmosphere::topographyDownslopeAcceleration);
                    }

                    float newu = tx * alongridge + nx * crossridge;
                    float newv = ty * alongridge + ny * crossridge;
                    const float newspeed = std::sqrt(newu * newu + newv * newv);

                    if (newspeed > 0.0f)
                    {
                        const float roughnessdrag = 1.0f - barrier * tuning::climate::atmosphere::topographySpeedReduction;
                        const float dragfactor = std::max(0.0f, roughnessdrag);
                        const float cappedspeed = std::min(maxvectorwind, newspeed * dragfactor);
                        const float speedscale = cappedspeed / newspeed;

                        newu = newu * speedscale;
                        newv = newv * speedscale;
                    }

                    nextu[x][y] = std::clamp(newu, -maxvectorwind, maxvectorwind);
                    nextv[x][y] = std::clamp(newv, -maxvectorwind, maxvectorwind);
                }
            }
        });

        windu.swap(nextu);
        windv.swap(nextv);
        smoothallfield(world, windu, 1);
        smoothallfield(world, windv, 1);
    }
}

void storeterrainverticalmotion(planet& world, int season, const floatgrid& macroterrain, const floatgrid& windu, const floatgrid& windv)
{
    const int width = world.width();
    const int height = world.height();

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int y = startrow; y <= endrow; y++)
        {
            const int ynorth = (y > 0) ? y - 1 : y;
            const int ysouth = (y < height) ? y + 1 : y;

            for (int x = 0; x <= width; x++)
            {
                const float u = windu[x][y];
                const float v = windv[x][y];
                const float speed = std::sqrt(u * u + v * v);

                if (speed < tuning::climate::atmosphere::topographyMinimumWindSpeed)
                {
                    world.setseasonaluplift(season, x, y, 0);
                    world.setseasonalsubsidence(season, x, y, 0);
                    continue;
                }

                const int xwest = wrapx(x - 1, width);
                const int xeast = wrapx(x + 1, width);
                const float gradx = (macroterrain[xeast][y] - macroterrain[xwest][y]) / 2.0f;
                const float grady = (macroterrain[x][ysouth] - macroterrain[x][ynorth]) / 2.0f;
                const float gradmag = std::sqrt(gradx * gradx + grady * grady);
                const float terrainhere = macroterrain[x][y];

                if (gradmag < tuning::climate::atmosphere::topographyMinimumRelief && terrainhere < tuning::climate::atmosphere::topographyMinimumRelief)
                {
                    world.setseasonaluplift(season, x, y, 0);
                    world.setseasonalsubsidence(season, x, y, 0);
                    continue;
                }

                const float dirx = u / speed;
                const float diry = v / speed;
                float lookaheadrise = 0.0f;
                float lookaheaddrop = 0.0f;

                for (int step = 1; step <= tuning::climate::atmosphere::topographyLookaheadDistance; step++)
                {
                    const float samplex = static_cast<float>(x) + dirx * static_cast<float>(step);
                    const float sampley = static_cast<float>(y) + diry * static_cast<float>(step);
                    const float change = samplewrappedfield(macroterrain, world, samplex, sampley) - terrainhere;

                    if (change > lookaheadrise)
                        lookaheadrise = change;

                    if (-change > lookaheaddrop)
                        lookaheaddrop = -change;
                }

                const float nx = gradx / std::max(gradmag, 0.0001f);
                const float ny = grady / std::max(gradmag, 0.0001f);
                const float crossridge = u * nx + v * ny;
                const float gradientbarrier = std::clamp(gradmag / tuning::climate::atmosphere::topographyGradientScale, 0.0f, 1.0f);
                const float lookaheadbarrier = std::clamp(lookaheadrise / tuning::climate::atmosphere::topographyLookaheadRiseScale, 0.0f, 1.0f);
                const float leefactor = std::clamp(lookaheaddrop / tuning::climate::atmosphere::topographyLookaheadRiseScale, 0.0f, 1.0f);
                const float barrier = std::max(gradientbarrier, lookaheadbarrier);
                const float uplift = std::max(0.0f, crossridge) * barrier / tuning::climate::atmosphere::topographyVerticalMotionWindScale;
                const float subsidence = std::max(0.0f, -crossridge) * std::max(gradientbarrier, leefactor) / tuning::climate::atmosphere::topographyVerticalMotionWindScale;

                world.setseasonaluplift(season, x, y, static_cast<int>(std::round(uplift * tuning::climate::atmosphere::topographyVerticalMotionStorageScale)));
                world.setseasonalsubsidence(season, x, y, static_cast<int>(std::round(subsidence * tuning::climate::atmosphere::topographyVerticalMotionStorageScale)));
            }
        }
    });
}
}

void createoceancurrentmap(planet& world)
{
    const int width = world.width();
    const int height = world.height();
    const int sealevel = world.sealevel();
    const int coastalsearchdistance = tuning::climate::oceancurrents::coastalSearchDistance;
    const float equatorialband = tuning::climate::oceancurrents::equatorialBand;
    const float midlatitudeband = tuning::climate::oceancurrents::midLatitudeBand;
    const float polarband = tuning::climate::oceancurrents::polarBand;
    const float countercurrentband = tuning::climate::oceancurrents::counterCurrentBand;
    const float retainedbasestrength = tuning::climate::oceancurrents::retainedBaseStrength;
    const float smoothingblend = tuning::climate::oceancurrents::smoothingBlend;
    const float blockedcomponentfactor = tuning::climate::oceancurrents::blockedComponentFactor;
    const float maxcurrentspeed = tuning::climate::oceancurrents::equatorialSpeed;

    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
    {
        const float latitudeshift = seasonlatitudephase[season] * world.tilt() * tuning::climate::oceancurrents::seasonalShiftFactor;

        floatgrid baseu(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid basev(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid currentu(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid currentv(width + 1, vector<float>(height + 1, 0.0f));

        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int x = 0; x <= width; x++)
                {
                    if (world.nom(x, y) > sealevel)
                        continue;

                    const float latitude = latitudeforrow(y, height);
                    const float effectivelatitude = latitude - latitudeshift;
                    const float absolutelatitude = std::fabs(effectivelatitude);

                    float u = 0.0f;
                    float v = 0.0f;

                    if (absolutelatitude <= countercurrentband)
                        u = tuning::climate::oceancurrents::counterCurrentSpeed;
                    else if (absolutelatitude < equatorialband)
                        u = -tuning::climate::oceancurrents::equatorialSpeed;
                    else if (absolutelatitude < midlatitudeband)
                        u = tuning::climate::oceancurrents::midLatitudeSpeed;
                    else
                        u = -tuning::climate::oceancurrents::polarSpeed;

                    const float polewardsign = (latitude >= 0.0f) ? -1.0f : 1.0f;

                    int westdistance = 0;
                    int eastdistance = 0;

                    const bool westland = landindir(world, x, y, -1, 0, coastalsearchdistance, westdistance);
                    const bool eastland = landindir(world, x, y, 1, 0, coastalsearchdistance, eastdistance);

                    if (absolutelatitude >= countercurrentband && absolutelatitude < midlatitudeband)
                    {
                        if (westland)
                            v += polewardsign * tuning::climate::oceancurrents::westernBoundarySpeed * coastalweight(westdistance, coastalsearchdistance);

                        if (eastland)
                            v -= polewardsign * tuning::climate::oceancurrents::easternBoundarySpeed * coastalweight(eastdistance, coastalsearchdistance);
                    }
                    else if (absolutelatitude >= midlatitudeband && absolutelatitude < polarband)
                    {
                        if (westland)
                            v -= polewardsign * tuning::climate::oceancurrents::subpolarBoundarySpeed * coastalweight(westdistance, coastalsearchdistance);

                        if (eastland)
                            v += polewardsign * tuning::climate::oceancurrents::subpolarBoundarySpeed * coastalweight(eastdistance, coastalsearchdistance);
                    }

                    baseu[x][y] = u;
                    basev[x][y] = v;
                    currentu[x][y] = u;
                    currentv[x][y] = v;
                }
            }
        });

        floatgrid nextu = currentu;
        floatgrid nextv = currentv;

        for (int iteration = 0; iteration < tuning::climate::oceancurrents::smoothingIterations; iteration++)
        {
            parallelforrows(0, height, [&](int startrow, int endrow)
            {
                for (int y = startrow; y <= endrow; y++)
                {
                    for (int x = 0; x <= width; x++)
                    {
                        if (world.nom(x, y) > sealevel)
                        {
                            nextu[x][y] = 0.0f;
                            nextv[x][y] = 0.0f;
                            continue;
                        }

                        float avgu = 0.0f;
                        float avgv = 0.0f;
                        float weighttotal = 0.0f;

                        for (int dy = -1; dy <= 1; dy++)
                        {
                            const int yy = y + dy;

                            if (yy < 0 || yy > height)
                                continue;

                            for (int dx = -1; dx <= 1; dx++)
                            {
                                const int xx = wrapx(x + dx, width);

                                if (world.nom(xx, yy) > sealevel)
                                    continue;

                                const float weight = (dx == 0 && dy == 0) ? 2.0f : 1.0f;
                                avgu += currentu[xx][yy] * weight;
                                avgv += currentv[xx][yy] * weight;
                                weighttotal += weight;
                            }
                        }

                        if (weighttotal > 0.0f)
                        {
                            avgu = avgu / weighttotal;
                            avgv = avgv / weighttotal;
                        }

                        float blendedu = currentu[x][y] * (1.0f - smoothingblend) + avgu * smoothingblend;
                        float blendedv = currentv[x][y] * (1.0f - smoothingblend) + avgv * smoothingblend;

                        blendedu = blendedu * retainedbasestrength + baseu[x][y] * (1.0f - retainedbasestrength);
                        blendedv = blendedv * retainedbasestrength + basev[x][y] * (1.0f - retainedbasestrength);

                        const bool eastblocked = world.nom(wrapx(x + 1, width), y) > sealevel;
                        const bool westblocked = world.nom(wrapx(x - 1, width), y) > sealevel;
                        const bool northblocked = (y == 0) || world.nom(x, y - 1) > sealevel;
                        const bool southblocked = (y == height) || world.nom(x, y + 1) > sealevel;

                        if (eastblocked && blendedu > 0.0f)
                            blendedu = blendedu * blockedcomponentfactor;

                        if (westblocked && blendedu < 0.0f)
                            blendedu = blendedu * blockedcomponentfactor;

                        if (northblocked && blendedv < 0.0f)
                            blendedv = blendedv * blockedcomponentfactor;

                        if (southblocked && blendedv > 0.0f)
                            blendedv = blendedv * blockedcomponentfactor;

                        nextu[x][y] = std::clamp(blendedu, -maxcurrentspeed, maxcurrentspeed);
                        nextv[x][y] = std::clamp(blendedv, -maxcurrentspeed, maxcurrentspeed);
                    }
                }
            });

            currentu.swap(nextu);
            currentv.swap(nextv);
        }

        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int x = 0; x <= width; x++)
                {
                    if (world.nom(x, y) > sealevel)
                    {
                        world.setseasonalcurrentu(season, x, y, 0);
                        world.setseasonalcurrentv(season, x, y, 0);
                        continue;
                    }

                    world.setseasonalcurrentu(season, x, y, static_cast<int>(std::round(currentu[x][y])));
                    world.setseasonalcurrentv(season, x, y, static_cast<int>(std::round(currentv[x][y])));
                }
            }
        });
    }
}

void createsurfacetemperaturemap(planet& world)
{
    const int width = world.width();
    const int height = world.height();
    const int sealevel = world.sealevel();
    const int coastalsearchdistance = tuning::climate::oceancurrents::coastalSearchDistance;
    const float maxcurrentspeed = tuning::climate::oceancurrents::equatorialSpeed;

    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
    {
        floatgrid basetemperatures(width + 1, vector<float>(height + 1, 0.0f));

        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int x = 0; x <= width; x++)
                {
                    if (world.nom(x, y) <= sealevel)
                        basetemperatures[x][y] = static_cast<float>(world.seasonaltemp(season, x, y));
                }
            }
        });

        smoothseasonalfield(world, basetemperatures, tuning::climate::sst::smoothingIterations);

        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int x = 0; x <= width; x++)
                {
                    if (world.nom(x, y) > sealevel)
                    {
                        world.setseasonalsst(season, x, y, 0);
                        world.setseasonalevaporation(season, x, y, 0);
                        world.setseasonalmoisture(season, x, y, 0);
                        continue;
                    }

                    const float currentu = static_cast<float>(world.seasonalcurrentu(season, x, y));
                    const float currentv = static_cast<float>(world.seasonalcurrentv(season, x, y));
                    const float magnitude = std::sqrt(currentu * currentu + currentv * currentv);

                    const int sourcex = wrapx(x - static_cast<int>(std::round((currentu / maxcurrentspeed) * tuning::climate::sst::advectionSampleDistance)), width);
                    const int sourcey = std::clamp(y - static_cast<int>(std::round((currentv / maxcurrentspeed) * tuning::climate::sst::advectionSampleDistance)), 0, height);

                    float sst = basetemperatures[x][y];
                    const float sourcetemperature = sampleoceanfield(basetemperatures, world, sourcex, sourcey);
                    sst = sst + (sourcetemperature - sst) * tuning::climate::sst::advectionBlend;

                    const float latitude = latitudeforrow(y, height);
                    const float polewardflow = ((latitude >= 0.0f) ? -currentv : currentv) / maxcurrentspeed;
                    const float equatorwardflow = ((latitude >= 0.0f) ? currentv : -currentv) / maxcurrentspeed;

                    int westdistance = 0;
                    int eastdistance = 0;

                    const bool westland = landindir(world, x, y, -1, 0, coastalsearchdistance, westdistance);
                    const bool eastland = landindir(world, x, y, 1, 0, coastalsearchdistance, eastdistance);

                    if (westland)
                        sst += std::max(0.0f, polewardflow) * tuning::climate::sst::westernBoundaryWarming * coastalweight(westdistance, coastalsearchdistance);

                    if (eastland)
                        sst -= std::max(0.0f, equatorwardflow) * tuning::climate::sst::easternBoundaryCooling * coastalweight(eastdistance, coastalsearchdistance);

                    sst = std::clamp(sst, tuning::climate::sst::minimumSst, tuning::climate::sst::maximumSst);

                    const float evaporation = std::max(0.0f, (sst + 10.0f) * tuning::climate::sst::evaporationScale + magnitude * tuning::climate::sst::evaporationCurrentBoost);

                    world.setseasonalsst(season, x, y, static_cast<int>(std::round(sst)));
                    world.setseasonalevaporation(season, x, y, static_cast<int>(std::round(evaporation)));
                    world.setseasonalmoisture(season, x, y, static_cast<int>(std::round(evaporation)));
                }
            }
        });
    }
}

namespace
{
float pressuresurfacetemperature(planet& world, int season, int x, int y)
{
    if (world.sea(x, y) == 1 && world.seasonalsst(season, x, y) != 0)
        return static_cast<float>(world.seasonalsst(season, x, y));

    return static_cast<float>(world.seasonaltemp(season, x, y));
}

void smoothallfield(planet& world, floatgrid& field, int iterations)
{
    const int width = world.width();
    const int height = world.height();
    floatgrid scratch = field;

    for (int iteration = 0; iteration < iterations; iteration++)
    {
        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int x = 0; x <= width; x++)
                {
                    float total = 0.0f;
                    float weighttotal = 0.0f;

                    for (int dy = -1; dy <= 1; dy++)
                    {
                        const int yy = y + dy;

                        if (yy < 0 || yy > height)
                            continue;

                        for (int dx = -1; dx <= 1; dx++)
                        {
                            const int xx = wrapx(x + dx, width);
                            const float weight = (dx == 0 && dy == 0) ? 2.0f : 1.0f;
                            total += field[xx][yy] * weight;
                            weighttotal += weight;
                        }
                    }

                    scratch[x][y] = (weighttotal > 0.0f) ? total / weighttotal : field[x][y];
                }
            }
        });

        field.swap(scratch);
    }
}

floatgrid buildcontinentalityfield(planet& world, int smoothingiterations, float exponent)
{
    const int width = world.width();
    const int height = world.height();
    floatgrid continentality(width + 1, vector<float>(height + 1, 0.0f));

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int y = startrow; y <= endrow; y++)
        {
            for (int x = 0; x <= width; x++)
                continentality[x][y] = world.sea(x, y) == 1 ? 0.0f : 1.0f;
        }
    });

    smoothallfield(world, continentality, smoothingiterations);

    if (std::fabs(exponent - 1.0f) > 0.001f)
    {
        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int x = 0; x <= width; x++)
                    continentality[x][y] = std::pow(std::clamp(continentality[x][y], 0.0f, 1.0f), exponent);
            }
        });
    }

    return continentality;
}
}

void createpressuremap(planet& world)
{
    const int width = world.width();
    const int height = world.height();
    const float thicknessresponse = climateatmosphere::hypsometricHeightResponseMetresPerKelvin(
        tuning::climate::circulation::surfaceReferencePressurePa,
        tuning::climate::circulation::upperReferencePressurePa);

    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
    {
        floatgrid surface(width + 1, vector<float>(height + 1, 0.0f));

        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int x = 0; x <= width; x++)
                    surface[x][y] = pressuresurfacetemperature(world, season, x, y);
            }
        });

        double temperaturesum = 0.0;
        double weighttotal = 0.0;

        for (int y = 0; y <= height; y++)
        {
            constexpr double pi = 3.14159265358979323846;
            const double weight = std::max(
                0.0,
                std::cos(static_cast<double>(latitudeforrow(y, height)) * pi / 180.0));

            for (int x = 0; x <= width; x++)
            {
                temperaturesum += weight * static_cast<double>(surface[x][y]);
                weighttotal += weight;
            }
        }

        const float globalmean = weighttotal > 0.0 ?
            static_cast<float>(temperaturesum / weighttotal) : 0.0f;
        floatgrid pressure(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid upperheight(width + 1, vector<float>(height + 1, 0.0f));

        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int x = 0; x <= width; x++)
                {
                    pressure[x][y] = 0.0f;
                    upperheight[x][y] = (surface[x][y] - globalmean) * thicknessresponse;
                }
            }
        });

        smoothallfield(world, upperheight, tuning::climate::circulation::thermalHeightSmoothingIterations);

        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int x = 0; x <= width; x++)
                {
                    world.setseasonalpressure(season, x, y, static_cast<int>(std::round(pressure[x][y])));
                    world.setseasonalupperheight(season, x, y, static_cast<int>(std::round(upperheight[x][y])));
                }
            }
        });
    }
}

void updatehorsebeltsfrompressure(planet& world)
{
    const int width = world.width();
    const int height = world.height();

    auto rowfromlatitude = [height](float latitude)
    {
        return std::clamp(static_cast<int>(std::round((90.0f - latitude) * static_cast<float>(height) / 180.0f)), 0, height);
    };

    const int defaultnorthouter = rowfromlatitude(32.0f);
    const int defaultnorthinner = rowfromlatitude(26.0f);
    const int defaultsouthinner = rowfromlatitude(-26.0f);
    const int defaultsouthouter = rowfromlatitude(-32.0f);
    const int minband = std::max(2, height / 60);

    vector<int> northouter(width + 1, defaultnorthouter);
    vector<int> northinner(width + 1, defaultnorthinner);
    vector<int> southinner(width + 1, defaultsouthinner);
    vector<int> southouter(width + 1, defaultsouthouter);

    for (int x = 0; x <= width; x++)
    {
        int northpeakpoleward = height;
        int northpeakequatorward = 0;
        int southpeakequatorward = height;
        int southpeakpoleward = 0;
        bool foundnorth = false;
        bool foundsouth = false;

        for (int season = 0; season < CLIMATESEASONCOUNT; season++)
        {
            float bestnorth = -1000000.0f;
            float bestsouth = -1000000.0f;
            int bestnorthy = defaultnorthinner;
            int bestsouthy = defaultsouthinner;

            for (int y = 0; y <= height; y++)
            {
                const float latitude = latitudeforrow(y, height);
                float pressure = 0.0f;

                for (int dx = -1; dx <= 1; dx++)
                {
                    const int xx = wrapx(x + dx, width);
                    pressure += static_cast<float>(world.seasonalpressure(season, xx, y));
                }

                pressure = pressure / 3.0f;

                if (latitude >= 15.0f && latitude <= 45.0f)
                {
                    if (pressure > bestnorth)
                    {
                        bestnorth = pressure;
                        bestnorthy = y;
                    }
                }

                if (latitude <= -15.0f && latitude >= -45.0f)
                {
                    if (pressure > bestsouth)
                    {
                        bestsouth = pressure;
                        bestsouthy = y;
                    }
                }
            }

            northpeakpoleward = std::min(northpeakpoleward, bestnorthy);
            northpeakequatorward = std::max(northpeakequatorward, bestnorthy);
            southpeakequatorward = std::min(southpeakequatorward, bestsouthy);
            southpeakpoleward = std::max(southpeakpoleward, bestsouthy);
            foundnorth = true;
            foundsouth = true;
        }

        if (foundnorth)
        {
            if (northpeakequatorward - northpeakpoleward < minband)
            {
                const int centre = (northpeakequatorward + northpeakpoleward) / 2;
                northpeakpoleward = std::max(rowfromlatitude(45.0f), centre - minband / 2);
                northpeakequatorward = std::min(rowfromlatitude(15.0f), northpeakpoleward + minband);
            }

            northouter[x] = northpeakpoleward;
            northinner[x] = northpeakequatorward;
        }

        if (foundsouth)
        {
            if (southpeakpoleward - southpeakequatorward < minband)
            {
                const int centre = (southpeakpoleward + southpeakequatorward) / 2;
                southpeakequatorward = std::max(rowfromlatitude(-15.0f), centre - minband / 2);
                southpeakpoleward = std::min(rowfromlatitude(-45.0f), southpeakequatorward + minband);
            }

            southinner[x] = southpeakequatorward;
            southouter[x] = southpeakpoleward;
        }
    }

    auto smoothband = [width](vector<int>& band)
    {
        vector<int> smoothed = band;

        for (int x = 0; x <= width; x++)
        {
            int total = 0;
            int count = 0;

            for (int dx = -2; dx <= 2; dx++)
            {
                total += band[wrapx(x + dx, width)];
                count++;
            }

            smoothed[x] = total / count;
        }

        band.swap(smoothed);
    };

    smoothband(northouter);
    smoothband(northinner);
    smoothband(southinner);
    smoothband(southouter);

    for (int x = 0; x <= width; x++)
    {
        if (northinner[x] <= northouter[x])
            northinner[x] = std::min(height, northouter[x] + minband);

        if (southouter[x] <= southinner[x])
            southouter[x] = std::min(height, southinner[x] + minband);

        world.sethorse(x, 1, northouter[x]);
        world.sethorse(x, 2, northinner[x]);
        world.sethorse(x, 3, southinner[x]);
        world.sethorse(x, 4, southouter[x]);
    }
}

void createvectorwindmap(planet& world)
{
    const int width = world.width();
    const int height = world.height();
    const int sealevel = world.sealevel();
    const float maxvectorwind = tuning::climate::atmosphere::maxVectorWind;
    const floatgrid continentality = buildcontinentalityfield(world, tuning::climate::atmosphere::landmaskSmoothingIterations, 1.0f);
    floatgrid macroterrain(width + 1, vector<float>(height + 1, 0.0f));

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int y = startrow; y <= endrow; y++)
        {
            for (int x = 0; x <= width; x++)
                macroterrain[x][y] = static_cast<float>(std::max(0, world.nom(x, y) - sealevel));
        }
    });

    smoothallfield(world, macroterrain, tuning::climate::atmosphere::topographySmoothingIterations);

    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
    {
        floatgrid pressure(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid basepressure(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid nextpressure(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid baseupperheight(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid upperheight(width + 1, vector<float>(height + 1, 0.0f));

        for (int y = 0; y <= height; y++)
        {
            for (int x = 0; x <= width; x++)
            {
                pressure[x][y] = static_cast<float>(world.seasonalpressure(season, x, y));
                basepressure[x][y] = pressure[x][y];
                nextpressure[x][y] = pressure[x][y];
                baseupperheight[x][y] = static_cast<float>(world.seasonalupperheight(season, x, y));
                upperheight[x][y] = baseupperheight[x][y];
            }
        }

        floatgrid windu(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid windv(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid upperu(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid upperv(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid vertical(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid nextupperheight = upperheight;
        floatgrid nextvertical = vertical;
        floatgrid surfacedivergence(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid upperdivergence(width + 1, vector<float>(height + 1, 0.0f));

        auto updatesurfaceflow = [&]()
        {
            parallelforrows(0, height, [&](int startrow, int endrow)
            {
                for (int y = startrow; y <= endrow; y++)
                {
                    const float latitude = latitudeforrow(y, height);
                    const auto spacing = climateatmosphere::cellSpacingMetres(
                        latitude,
                        width + 1,
                        height + 1,
                        tuning::climate::atmosphere::referencePlanetRadiusMetres);

                    for (int x = 0; x <= width; x++)
                    {
                        const int ynorth = (y > 0) ? y - 1 : y;
                        const int ysouth = (y < height) ? y + 1 : y;
                        const int xwest = wrapx(x - 1, width);
                        const int xeast = wrapx(x + 1, width);
                        const float pressuregradienteast =
                            (pressure[xeast][y] - pressure[xwest][y]) *
                            tuning::climate::atmosphere::pressurePascalsPerHectopascal /
                            (2.0f * spacing.zonalMetres);
                        const float pressuregradientnorth =
                            (pressure[x][ynorth] - pressure[x][ysouth]) *
                            tuning::climate::atmosphere::pressurePascalsPerHectopascal /
                            (2.0f * spacing.meridionalMetres);
                        const float continental = std::clamp(continentality[x][y], 0.0f, 1.0f);
                        const float relief = std::clamp(macroterrain[x][y] / 4500.0f, 0.0f, 1.0f);
                        const float baseDragTime =
                            tuning::climate::atmosphere::oceanBoundaryLayerDragTimeSeconds +
                            continental * (tuning::climate::atmosphere::landBoundaryLayerDragTimeSeconds -
                                tuning::climate::atmosphere::oceanBoundaryLayerDragTimeSeconds);
                        const float dragTime = baseDragTime + relief *
                            (tuning::climate::atmosphere::highReliefDragTimeSeconds - baseDragTime);
                        const auto wind = climateatmosphere::steadyRayleighCoriolisWind(
                            -pressuregradienteast / tuning::climate::atmosphere::surfaceAirDensityKgM3,
                            -pressuregradientnorth / tuning::climate::atmosphere::surfaceAirDensityKgM3,
                            latitude,
                            dragTime,
                            tuning::climate::atmosphere::rotationRatePerSecond,
                            world.rotation() ? 1.0f : -1.0f);

                        windu[x][y] = std::clamp(wind.eastMetresPerSecond, -maxvectorwind, maxvectorwind);
                        windv[x][y] = std::clamp(wind.southMetresPerSecond, -maxvectorwind, maxvectorwind);
                    }
                }
            });

            smoothallfield(world, windu, tuning::climate::atmosphere::smoothingIterations);
            smoothallfield(world, windv, tuning::climate::atmosphere::smoothingIterations);

            parallelforrows(0, height, [&](int startrow, int endrow)
            {
                for (int y = startrow; y <= endrow; y++)
                {
                    const int ynorth = (y > 0) ? y - 1 : y;
                    const int ysouth = (y < height) ? y + 1 : y;
                    const auto spacing = climateatmosphere::cellSpacingMetres(
                        latitudeforrow(y, height),
                        width + 1,
                        height + 1,
                        tuning::climate::atmosphere::referencePlanetRadiusMetres);

                    for (int x = 0; x <= width; x++)
                    {
                        const int xwest = wrapx(x - 1, width);
                        const int xeast = wrapx(x + 1, width);
                        surfacedivergence[x][y] =
                            ((windu[xeast][y] - windu[xwest][y]) / (2.0f * spacing.zonalMetres) +
                                (windv[x][ysouth] - windv[x][ynorth]) / (2.0f * spacing.meridionalMetres)) *
                            tuning::climate::circulation::secondsPerDay;
                    }
                }
            });

            smoothallfield(
                world,
                surfacedivergence,
                tuning::climate::circulation::divergenceSmoothingIterations);
        };

        updatesurfaceflow();

        for (int iteration = 0; iteration < tuning::climate::circulation::iterations; iteration++)
        {
            parallelforrows(0, height, [&](int startrow, int endrow)
            {
                for (int y = startrow; y <= endrow; y++)
                {
                    const int ynorth = (y > 0) ? y - 1 : y;
                    const int ysouth = (y < height) ? y + 1 : y;
                    const float latitude = latitudeforrow(y, height);
                    const auto spacing = climateatmosphere::cellSpacingMetres(
                        latitude,
                        width + 1,
                        height + 1,
                        tuning::climate::atmosphere::referencePlanetRadiusMetres);

                    for (int x = 0; x <= width; x++)
                    {
                        const int xwest = wrapx(x - 1, width);
                        const int xeast = wrapx(x + 1, width);
                        const float heightgradienteast =
                            (upperheight[xeast][y] - upperheight[xwest][y]) /
                            (2.0f * spacing.zonalMetres);
                        const float heightgradientnorth =
                            (upperheight[x][ynorth] - upperheight[x][ysouth]) /
                            (2.0f * spacing.meridionalMetres);
                        const auto wind = climateatmosphere::steadyRayleighCoriolisWind(
                            -tuning::climate::atmosphere::gravityMetresPerSecondSquared * heightgradienteast,
                            -tuning::climate::atmosphere::gravityMetresPerSecondSquared * heightgradientnorth,
                            latitude,
                            tuning::climate::circulation::upperLayerDragTimeSeconds,
                            tuning::climate::atmosphere::rotationRatePerSecond,
                            world.rotation() ? 1.0f : -1.0f);

                        upperu[x][y] = std::clamp(wind.eastMetresPerSecond, -maxvectorwind, maxvectorwind);
                        upperv[x][y] = std::clamp(wind.southMetresPerSecond, -maxvectorwind, maxvectorwind);
                    }
                }
            });

            smoothallfield(world, upperu, tuning::climate::circulation::windSmoothingIterations);
            smoothallfield(world, upperv, tuning::climate::circulation::windSmoothingIterations);

            parallelforrows(0, height, [&](int startrow, int endrow)
            {
                for (int y = startrow; y <= endrow; y++)
                {
                    const int ynorth = (y > 0) ? y - 1 : y;
                    const int ysouth = (y < height) ? y + 1 : y;
                    const auto spacing = climateatmosphere::cellSpacingMetres(
                        latitudeforrow(y, height),
                        width + 1,
                        height + 1,
                        tuning::climate::atmosphere::referencePlanetRadiusMetres);

                    for (int x = 0; x <= width; x++)
                    {
                        const int xwest = wrapx(x - 1, width);
                        const int xeast = wrapx(x + 1, width);
                        upperdivergence[x][y] =
                            ((upperu[xeast][y] - upperu[xwest][y]) / (2.0f * spacing.zonalMetres) +
                                (upperv[x][ysouth] - upperv[x][ynorth]) / (2.0f * spacing.meridionalMetres)) *
                            tuning::climate::circulation::secondsPerDay;
                    }
                }
            });

            smoothallfield(
                world,
                upperdivergence,
                tuning::climate::circulation::divergenceSmoothingIterations);

            parallelforrows(0, height, [&](int startrow, int endrow)
            {
                for (int y = startrow; y <= endrow; y++)
                {
                    const int ynorth = (y > 0) ? y - 1 : y;
                    const int ysouth = (y < height) ? y + 1 : y;

                    for (int x = 0; x <= width; x++)
                    {
                        const int xwest = wrapx(x - 1, width);
                        const int xeast = wrapx(x + 1, width);
                        const float targetvertical =
                            (upperdivergence[x][y] - surfacedivergence[x][y]) *
                            tuning::climate::circulation::layerPressureDepthHpa * 0.5f;
                        const float clampedvertical = std::clamp(
                            targetvertical,
                            -tuning::climate::circulation::maximumVerticalVelocity,
                            tuning::climate::circulation::maximumVerticalVelocity);
                        const float updatedvertical = vertical[x][y] +
                            (clampedvertical - vertical[x][y]) * tuning::climate::circulation::verticalRelaxation;
                        const float upperlaplacian =
                            (upperheight[xwest][y] + upperheight[xeast][y] + upperheight[x][ynorth] + upperheight[x][ysouth]) * 0.25f - upperheight[x][y];

                        nextvertical[x][y] = updatedvertical;
                        nextupperheight[x][y] = upperheight[x][y] + tuning::climate::circulation::upperHeightRelaxation *
                            ((baseupperheight[x][y] - upperheight[x][y]) +
                                upperlaplacian * tuning::climate::circulation::upperHeightDiffusion);
                    }
                }
            });

            upperheight.swap(nextupperheight);
            vertical.swap(nextvertical);

            parallelforrows(0, height, [&](int startrow, int endrow)
            {
                for (int y = startrow; y <= endrow; y++)
                {
                    const int ynorth = (y > 0) ? y - 1 : y;
                    const int ysouth = (y < height) ? y + 1 : y;

                    for (int x = 0; x <= width; x++)
                    {
                        const int xwest = wrapx(x - 1, width);
                        const int xeast = wrapx(x + 1, width);
                        const float pressurelaplacian =
                            (pressure[xwest][y] + pressure[xeast][y] +
                                pressure[x][ynorth] + pressure[x][ysouth]) * 0.25f -
                            pressure[x][y];
                        const float columndivergence =
                            (surfacedivergence[x][y] + upperdivergence[x][y]) * 0.5f;
                        const float pressuretendency = std::clamp(
                            -tuning::climate::circulation::surfacePressureReferenceHpa * columndivergence,
                            -tuning::climate::circulation::maximumSurfacePressureTendencyHpaPerDay,
                            tuning::climate::circulation::maximumSurfacePressureTendencyHpaPerDay);
                        const float restoringtendency =
                            (basepressure[x][y] - pressure[x][y]) /
                            tuning::climate::circulation::surfacePressureRestoringTimeDays;
                        nextpressure[x][y] = std::clamp(
                            pressure[x][y] + tuning::climate::circulation::timeStepDays *
                                (pressuretendency + restoringtendency) +
                                pressurelaplacian * tuning::climate::circulation::surfacePressureDiffusion,
                            -tuning::climate::circulation::maximumSurfacePressureAnomalyHpa,
                            tuning::climate::circulation::maximumSurfacePressureAnomalyHpa);
                    }
                }
            });

            double pressuretotal = 0.0;
            double pressureweight = 0.0;

            for (int y = 0; y <= height; y++)
            {
                constexpr double pi = 3.14159265358979323846;
                const double areaweight = std::max(
                    0.0,
                    std::cos(static_cast<double>(latitudeforrow(y, height)) * pi / 180.0));

                for (int x = 0; x <= width; x++)
                {
                    pressuretotal += areaweight * static_cast<double>(nextpressure[x][y]);
                    pressureweight += areaweight;
                }
            }

            const float meanpressureanomaly = pressureweight > 0.0 ?
                static_cast<float>(pressuretotal / pressureweight) : 0.0f;

            parallelforrows(0, height, [&](int startrow, int endrow)
            {
                for (int y = startrow; y <= endrow; y++)
                {
                    for (int x = 0; x <= width; x++)
                        pressure[x][y] = nextpressure[x][y] - meanpressureanomaly;
                }
            });

            updatesurfaceflow();
        }

        applytopographicwindeffects(world, macroterrain, windu, windv);
        storeterrainverticalmotion(world, season, macroterrain, windu, windv);
        smoothallfield(world, upperheight, tuning::climate::circulation::windSmoothingIterations);

        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int x = 0; x <= width; x++)
                {
                    const float dynamicvertical = std::clamp(
                        vertical[x][y],
                        -tuning::climate::circulation::maximumVerticalVelocity,
                        tuning::climate::circulation::maximumVerticalVelocity);

                    world.setseasonalpressure(season, x, y, static_cast<int>(std::round(pressure[x][y])));
                    world.setseasonalupperheight(season, x, y, static_cast<int>(std::round(upperheight[x][y])));
                    world.setseasonaluwind(season, x, y, static_cast<int>(std::round(std::clamp(windu[x][y], -maxvectorwind, maxvectorwind))));
                    world.setseasonalvwind(season, x, y, static_cast<int>(std::round(std::clamp(windv[x][y], -maxvectorwind, maxvectorwind))));
                    world.setseasonalupperuwind(season, x, y, static_cast<int>(std::round(std::clamp(upperu[x][y], -maxvectorwind, maxvectorwind))));
                    world.setseasonaluppervwind(season, x, y, static_cast<int>(std::round(std::clamp(upperv[x][y], -maxvectorwind, maxvectorwind))));
                    world.setseasonalverticalvelocity(season, x, y, static_cast<int>(std::round(
                        dynamicvertical * tuning::climate::circulation::verticalVelocityStorageScale)));
                }
            }
        });
    }

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int y = startrow; y <= endrow; y++)
        {
            for (int x = 0; x <= width; x++)
            {
                float averageu = 0.0f;
                float averageupperu = 0.0f;

                for (int season = 0; season < CLIMATESEASONCOUNT; season++)
                {
                    averageu += static_cast<float>(world.seasonaluwind(season, x, y));
                    averageupperu += static_cast<float>(world.seasonalupperuwind(season, x, y));
                }

                averageu = averageu / static_cast<float>(CLIMATESEASONCOUNT);
                averageupperu = averageupperu / static_cast<float>(CLIMATESEASONCOUNT);

                int scalarwind = static_cast<int>(std::round(averageu / tuning::climate::atmosphere::scalarWindDivisor));
                scalarwind = std::clamp(scalarwind, -10, 10);

                if (std::fabs(averageu) < tuning::climate::atmosphere::minimumScalarZonalWind || scalarwind == 0)
                {
                    const float transportu = averageu + (averageupperu - averageu) *
                        tuning::climate::moistureadvection::upperWindTransportFraction;

                    if (transportu > 0.25f)
                        world.setwind(x, y, 101);
                    else if (transportu < -0.25f)
                        world.setwind(x, y, 99);
                    else
                        world.setwind(x, y, 0);
                }
                else
                {
                    if (scalarwind > 0)
                        scalarwind = std::max(1, scalarwind);
                    else
                        scalarwind = std::min(-1, scalarwind);

                    world.setwind(x, y, scalarwind);
                }
            }
        }
    });

    for (int y = 0; y <= height; y++)
        world.setwind(0, y, world.wind(width, y));
}

void createadvectedrainfall(planet& world, vector<vector<int>>& inland, vector<vector<int>>& fractal)
{
    (void)inland;
    (void)fractal;
    const int width = world.width();
    const int height = world.height();
    const int sealevel = world.sealevel();
    const int baselineiterations = tuning::climate::moistureadvection::iterations;
    const int iterations = tuning::climateresolution::scaleDistance(
        baselineiterations, width, height);
    const float periterationfactor = static_cast<float>(baselineiterations) / static_cast<float>(iterations);
    const float timestepseconds = tuning::climate::moistureadvection::advectionTimeStepSeconds * periterationfactor;
    const float condensationefficiency = 1.0f - std::pow(
        1.0f - tuning::climate::moistureadvection::condensationEfficiency, periterationfactor);
    floatgrid soilmoisture(width + 1, vector<float>(height + 1, 0.0f));
    floatgrid moisture(width + 1, vector<float>(height + 1, 0.0f));

    auto runsolver = [&](int season)
    {
        floatgrid totalrain(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid totalseaevaporation(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid totallandevaporation(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid runoff(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid uplift(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid descent(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid dynamicvertical(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid saturationcapacity(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid surfacewindspeed(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid transportwindu(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid transportwindv(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid available(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid zonalmoisture(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid nextmoisture(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid fluxconvergence(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid totalconvergence(width + 1, vector<float>(height + 1, 0.0f));
        climatephysics::WaterBudget budget;

        for (int y = 0; y <= height; y++)
        {
            for (int x = 0; x <= width; x++)
            {
                budget.initialAtmosphericStorage += moisture[x][y];

                if (world.sea(x, y) == 0)
                    budget.initialSoilStorage += soilmoisture[x][y];
            }
        }

        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int x = 0; x <= width; x++)
                {
                    const float u = static_cast<float>(world.seasonaluwind(season, x, y));
                    const float v = static_cast<float>(world.seasonalvwind(season, x, y));
                    const float upperu = static_cast<float>(world.seasonalupperuwind(season, x, y));
                    const float upperv = static_cast<float>(world.seasonaluppervwind(season, x, y));
                    const float temperature = (world.sea(x, y) == 1) ? static_cast<float>(world.seasonalsst(season, x, y)) : static_cast<float>(world.seasonaltemp(season, x, y));
                    const float elevation = static_cast<float>(std::max(0, world.map(x, y) - sealevel));

                    saturationcapacity[x][y] = climatephysics::saturationColumnWater(temperature, elevation);
                    surfacewindspeed[x][y] = std::max(
                        tuning::climate::moistureadvection::minimumSurfaceWind,
                        std::sqrt(u * u + v * v));
                    transportwindu[x][y] = u + (upperu - u) *
                        tuning::climate::moistureadvection::upperWindTransportFraction;
                    transportwindv[x][y] = v + (upperv - v) *
                        tuning::climate::moistureadvection::upperWindTransportFraction;
                    uplift[x][y] = static_cast<float>(world.seasonaluplift(season, x, y)) / tuning::climate::atmosphere::topographyVerticalMotionStorageScale;
                    descent[x][y] = static_cast<float>(world.seasonalsubsidence(season, x, y)) / tuning::climate::atmosphere::topographyVerticalMotionStorageScale;
                    dynamicvertical[x][y] = static_cast<float>(world.seasonalverticalvelocity(season, x, y)) /
                        tuning::climate::circulation::verticalVelocityStorageScale;
                }
            }
        });

        for (int iteration = 0; iteration < iterations; iteration++)
        {
            parallelforrows(0, height, [&](int startrow, int endrow)
            {
                for (int y = startrow; y <= endrow; y++)
                {
                    for (int x = 0; x <= width; x++)
                    {
                        const float retained = std::max(0.0f, moisture[x][y]);
                        const bool sea = world.sea(x, y) == 1;
                        const float temperature = sea ?
                            static_cast<float>(world.seasonalsst(season, x, y)) :
                            static_cast<float>(world.seasonaltemp(season, x, y));
                        const float elevation = static_cast<float>(std::max(0, world.map(x, y) - sealevel));
                        const float relativehumidity = std::clamp(
                            retained / saturationcapacity[x][y], 0.0f, 1.0f);
                        const float potentialevaporation = climatephysics::bulkAerodynamicEvaporationMm(
                            temperature,
                            elevation,
                            surfacewindspeed[x][y],
                            relativehumidity,
                            timestepseconds,
                            tuning::climate::moistureadvection::surfaceExchangeCoefficient);
                        float seaevaporation = 0.0f;
                        float landevaporation = 0.0f;

                        if (sea)
                        {
                            seaevaporation = potentialevaporation;

                            if (world.seaice(x, y) == 1)
                                seaevaporation *= tuning::climate::moistureadvection::seaIceFactor;

                            if (world.seaice(x, y) == 2)
                                seaevaporation = 0.0f;

                            totalseaevaporation[x][y] += seaevaporation;
                        }
                        else if (soilmoisture[x][y] > 0.0f)
                        {
                            const float wateravailability = std::clamp(
                                soilmoisture[x][y] / tuning::climate::moistureadvection::landSoilMoistureCapacity,
                                0.0f,
                                1.0f);
                            landevaporation = std::min(
                                soilmoisture[x][y],
                                potentialevaporation * wateravailability *
                                    tuning::climate::moistureadvection::landSurfaceResistance);
                            soilmoisture[x][y] -= landevaporation;
                            totallandevaporation[x][y] += landevaporation;
                        }

                        const float localavailable = retained + seaevaporation + landevaporation;
                        available[x][y] = localavailable;
                    }
                }
            });

            parallelforrows(0, height, [&](int startrow, int endrow)
            {
                for (int y = startrow; y <= endrow; y++)
                {
                    const float latitudeRadians = latitudeforrow(y, height) * 3.14159265358979323846f / 180.0f;
                    const float zonalcircumference = 2.0f * 3.14159265358979323846f *
                        tuning::climate::moistureadvection::referencePlanetRadiusMetres *
                        std::max(0.05f, std::cos(latitudeRadians));
                    const float cellwidthmetres = zonalcircumference / static_cast<float>(width + 1);

                    for (int x = 0; x <= width; x++)
                        zonalmoisture[x][y] = 0.0f;

                    for (int x = 0; x <= width; x++)
                    {
                        const float u = transportwindu[x][y];
                        const float displacement = std::clamp(
                            u * timestepseconds / cellwidthmetres,
                            -tuning::climate::moistureadvection::maximumAdvectionCellsPerStep,
                            tuning::climate::moistureadvection::maximumAdvectionCellsPerStep);
                        const float transportfraction = std::min(
                            tuning::climate::moistureadvection::transportMaxFraction,
                            std::fabs(displacement));
                        const float targetposition = static_cast<float>(x) +
                            (std::fabs(displacement) < 1.0f ?
                                static_cast<float>((displacement > 0.0f) - (displacement < 0.0f)) :
                                displacement);
                        const int targetbase = static_cast<int>(std::floor(targetposition));
                        const float targetfraction = targetposition - std::floor(targetposition);
                        const int targetx0 = wrapx(targetbase, width);
                        const int targetx1 = wrapx(targetbase + 1, width);
                        const float transported = available[x][y] * transportfraction;

                        zonalmoisture[x][y] += available[x][y] - transported;
                        zonalmoisture[targetx0][y] += transported * (1.0f - targetfraction);
                        zonalmoisture[targetx1][y] += transported * targetfraction;
                    }
                }
            });

            parallelforrows(0, width, [&](int startcolumn, int endcolumn)
            {
                const float meridionalcellheight = 3.14159265358979323846f *
                    tuning::climate::moistureadvection::referencePlanetRadiusMetres /
                    static_cast<float>(height + 1);

                for (int x = startcolumn; x <= endcolumn; x++)
                {
                    for (int y = 0; y <= height; y++)
                        nextmoisture[x][y] = 0.0f;

                    for (int y = 0; y <= height; y++)
                    {
                        const float v = transportwindv[x][y];
                        const float displacement = std::clamp(
                            v * timestepseconds / meridionalcellheight,
                            -tuning::climate::moistureadvection::maximumAdvectionCellsPerStep,
                            tuning::climate::moistureadvection::maximumAdvectionCellsPerStep);
                        const float transportfraction = std::min(
                            tuning::climate::moistureadvection::transportMaxFraction,
                            std::fabs(displacement));
                        const float targetposition = std::clamp(
                            static_cast<float>(y) +
                                (std::fabs(displacement) < 1.0f ?
                                    static_cast<float>((displacement > 0.0f) - (displacement < 0.0f)) :
                                    displacement),
                            0.0f,
                            static_cast<float>(height));
                        const int targety0 = static_cast<int>(std::floor(targetposition));
                        const int targety1 = std::min(targety0 + 1, height);
                        const float targetfraction = targetposition - static_cast<float>(targety0);
                        const float transported = zonalmoisture[x][y] * transportfraction;

                        nextmoisture[x][y] += zonalmoisture[x][y] - transported;
                        nextmoisture[x][targety0] += transported * (1.0f - targetfraction);
                        nextmoisture[x][targety1] += transported * targetfraction;
                    }
                }
            });

            parallelforrows(0, height, [&](int startrow, int endrow)
            {
                for (int y = startrow; y <= endrow; y++)
                {
                    for (int x = 0; x <= width; x++)
                    {

                        const float temperature = static_cast<float>(
                            (world.sea(x, y) == 1) ?
                                world.seasonalsst(season, x, y) :
                                world.seasonaltemp(season, x, y));
                        const float elevation = static_cast<float>(
                            std::max(0, world.map(x, y) - sealevel));
                        const float moisturecapacity = climatephysics::saturationColumnWater(
                            temperature, elevation);
                        const float convergenceamount = std::max(0.0f, nextmoisture[x][y] - available[x][y]);

                        fluxconvergence[x][y] = convergenceamount / moisturecapacity;
                        totalconvergence[x][y] += fluxconvergence[x][y];
                    }
                }
            });

            parallelforrows(0, height, [&](int startrow, int endrow)
            {
                for (int y = startrow; y <= endrow; y++)
                {
                    for (int x = 0; x <= width; x++)
                    {
                        const bool sea = world.sea(x, y) == 1;
                        const float temperature = sea ? static_cast<float>(world.seasonalsst(season, x, y)) : static_cast<float>(world.seasonaltemp(season, x, y));
                        const float elevation = static_cast<float>(std::max(0, world.map(x, y) - sealevel));
                        const float liftingcooling = std::clamp(
                            std::max(0.0f, dynamicvertical[x][y]) * tuning::climate::moistureadvection::dynamicVerticalCooling +
                                uplift[x][y] * tuning::climate::moistureadvection::topographicUpliftCooling,
                            0.0f,
                            tuning::climate::moistureadvection::maximumParcelTemperatureAdjustment);
                        const float subsidencewarming = std::clamp(
                            std::max(0.0f, -dynamicvertical[x][y]) * tuning::climate::moistureadvection::dynamicSubsidenceWarming +
                                descent[x][y] * tuning::climate::moistureadvection::topographicDescentWarming,
                            0.0f,
                            tuning::climate::moistureadvection::maximumParcelTemperatureAdjustment);
                        const float parceltemperature = temperature - liftingcooling + subsidencewarming;
                        const float moisturecapacity = climatephysics::saturationColumnWater(
                            parceltemperature, elevation);
                        const float availablemoisture = std::max(0.0f, nextmoisture[x][y]);
                        const float condensablewater = std::max(0.0f, availablemoisture - moisturecapacity);
                        float precipitation = condensablewater * condensationefficiency;

                        precipitation = std::clamp(precipitation, 0.0f, availablemoisture);
                        totalrain[x][y] += precipitation;
                        moisture[x][y] = std::max(0.0f, availablemoisture - precipitation);

                        if (sea == false)
                        {
                            const float storageavailable = std::max(
                                0.0f,
                                tuning::climate::moistureadvection::landSoilMoistureCapacity - soilmoisture[x][y]);
                            const float infiltration = std::min(
                                storageavailable,
                                precipitation * tuning::climate::moistureadvection::landInfiltrationFraction);

                            soilmoisture[x][y] += infiltration;
                            runoff[x][y] += precipitation - infiltration;
                        }
                    }
                }
            });
        }

        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int x = 0; x <= width; x++)
                {
                    const int rain = std::clamp(
                        static_cast<int>(std::round(
                            totalrain[x][y] * tuning::climate::moistureadvection::rainfallScale)),
                        0,
                        static_cast<int>(std::numeric_limits<short>::max()));
                    world.setseasonalrain(season, x, y, rain);

                    if (season == seasonjanuary)
                        world.setjanrain(x, y, rain);

                    if (season == seasonjuly)
                        world.setjulrain(x, y, rain);

                    world.setseasonalmoisture(season, x, y, static_cast<int>(std::round(moisture[x][y])));
                    world.setseasonalconvergence(season, x, y, static_cast<int>(std::round(
                        (totalconvergence[x][y] / static_cast<float>(iterations)) * tuning::climate::moistureadvection::convergenceStorageScale)));
                }
            }
        });

        for (int y = 0; y <= height; y++)
        {
            for (int x = 0; x <= width; x++)
            {
                if (world.sea(x, y) == 1)
                {
                    budget.oceanEvaporation += totalseaevaporation[x][y];
                    budget.oceanPrecipitation += totalrain[x][y];
                }
                else
                {
                    budget.landEvaporation += totallandevaporation[x][y];
                    budget.landPrecipitation += totalrain[x][y];
                    budget.runoff += runoff[x][y];
                    budget.soilStorage += soilmoisture[x][y];
                }

                budget.atmosphericStorage += moisture[x][y];
            }
        }

        return budget;
    };

    std::array<climatephysics::WaterBudget, CLIMATESEASONCOUNT> finalbudgets{};
    climatephysics::HydrologySpinupDiagnostics diagnostics;

    for (int cycle = 0; cycle < tuning::climate::moistureadvection::maximumSpinupCycles; cycle++)
    {
        const floatgrid initialmoisture = moisture;
        const floatgrid initialsoilmoisture = soilmoisture;

        for (int season = 0; season < CLIMATESEASONCOUNT; season++)
            finalbudgets[season] = runsolver(season);

        double absolutechange = 0.0;
        double referencestorage = 0.0;
        double atmosphericstorage = 0.0;
        double soilstorage = 0.0;

        for (int y = 0; y <= height; y++)
        {
            for (int x = 0; x <= width; x++)
            {
                absolutechange += std::abs(static_cast<double>(moisture[x][y] - initialmoisture[x][y]));
                referencestorage += 0.5 * static_cast<double>(moisture[x][y] + initialmoisture[x][y]);
                atmosphericstorage += moisture[x][y];

                if (world.sea(x, y) == 0)
                {
                    absolutechange += std::abs(static_cast<double>(soilmoisture[x][y] - initialsoilmoisture[x][y]));
                    referencestorage += 0.5 * static_cast<double>(soilmoisture[x][y] + initialsoilmoisture[x][y]);
                    soilstorage += soilmoisture[x][y];
                }
            }
        }

        diagnostics.cyclesCompleted = cycle + 1;
        diagnostics.relativeStorageChange = absolutechange / std::max(1.0, referencestorage);
        diagnostics.atmosphericStorage = atmosphericstorage;
        diagnostics.soilStorage = soilstorage;
        diagnostics.converged = diagnostics.cyclesCompleted >=
                tuning::climate::moistureadvection::minimumSpinupCycles &&
            diagnostics.relativeStorageChange <=
                tuning::climate::moistureadvection::spinupRelativeStorageTolerance;

        std::cout
            << "Climate hydrology spinup cycle=" << diagnostics.cyclesCompleted
            << " relative_storage_change=" << diagnostics.relativeStorageChange
            << " atmospheric_storage=" << diagnostics.atmosphericStorage
            << " soil_storage=" << diagnostics.soilStorage
            << " converged=" << (diagnostics.converged ? 1 : 0)
            << '\n';

        if (diagnostics.converged)
            break;
    }

    climatephysics::setLastHydrologySpinupDiagnostics(diagnostics);

    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
    {
        const climatephysics::WaterBudget& budget = finalbudgets[season];
        climatephysics::setLastWaterBudget(season, budget);
        std::cout
            << "Climate water budget season=" << season
            << " initial_atmospheric_storage=" << budget.initialAtmosphericStorage
            << " initial_soil_storage=" << budget.initialSoilStorage
            << " ocean_evaporation=" << budget.oceanEvaporation
            << " land_evaporation=" << budget.landEvaporation
            << " ocean_precipitation=" << budget.oceanPrecipitation
            << " land_precipitation=" << budget.landPrecipitation
            << " runoff=" << budget.runoff
            << " atmospheric_storage=" << budget.atmosphericStorage
            << " soil_storage=" << budget.soilStorage
            << " residual=" << budget.residual()
            << " relative_residual=" << budget.relativeResidual()
            << '\n';
    }
}
