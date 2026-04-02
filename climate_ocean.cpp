#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

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
float gaussianpressurebell(float x, float centre, float width)
{
    const float diff = (x - centre) / width;
    return std::exp(-(diff * diff));
}

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
}

void createpressuremap(planet& world)
{
    const int width = world.width();
    const int height = world.height();
    const int sealevel = world.sealevel();

    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
    {
        const float latitudeshift = seasonlatitudephase[season] * world.tilt() * tuning::climate::pressure::seasonalShiftFactor;
        floatgrid surface(width + 1, vector<float>(height + 1, 0.0f));
        vector<float> zonalmean(height + 1, 0.0f);

        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int x = 0; x <= width; x++)
                    surface[x][y] = pressuresurfacetemperature(world, season, x, y);
            }
        });

        for (int y = 0; y <= height; y++)
        {
            float total = 0.0f;

            for (int x = 0; x <= width; x++)
                total += surface[x][y];

            zonalmean[y] = total / static_cast<float>(width + 1);
        }

        floatgrid pressure(width + 1, vector<float>(height + 1, 0.0f));

        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int x = 0; x <= width; x++)
                {
                    const float latitude = latitudeforrow(y, height) - latitudeshift;
                    const float absolutelatitude = std::fabs(latitude);
                    const float basepressure =
                        -tuning::climate::pressure::equatorialLow * gaussianpressurebell(absolutelatitude, 0.0f, tuning::climate::pressure::equatorialWidth) +
                        tuning::climate::pressure::subtropicalHigh * gaussianpressurebell(absolutelatitude, tuning::climate::pressure::subtropicalLatitude, tuning::climate::pressure::subtropicalWidth) -
                        tuning::climate::pressure::subpolarLow * gaussianpressurebell(absolutelatitude, tuning::climate::pressure::subpolarLatitude, tuning::climate::pressure::subpolarWidth) +
                        tuning::climate::pressure::polarHigh * gaussianpressurebell(absolutelatitude, tuning::climate::pressure::polarLatitude, tuning::climate::pressure::polarWidth);

                    const bool land = world.nom(x, y) > sealevel;
                    const float thermalresponse = land ? tuning::climate::pressure::landThermalResponse : tuning::climate::pressure::oceanThermalResponse;
                    const float temperatureanomaly = surface[x][y] - zonalmean[y];
                    const float elevation = static_cast<float>(std::max(0, world.nom(x, y) - sealevel));

                    pressure[x][y] = basepressure - temperatureanomaly * thermalresponse - elevation * tuning::climate::pressure::elevationResponse;
                }
            }
        });

        smoothallfield(world, pressure, tuning::climate::pressure::smoothingIterations);

        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int x = 0; x <= width; x++)
                    world.setseasonalpressure(season, x, y, static_cast<int>(std::round(pressure[x][y])));
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

        for (int y = 0; y <= height; y++)
        {
            for (int x = 0; x <= width; x++)
                pressure[x][y] = static_cast<float>(world.seasonalpressure(season, x, y));
        }

        floatgrid windu(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid windv(width + 1, vector<float>(height + 1, 0.0f));

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
                    const float dpdx = (pressure[xeast][y] - pressure[xwest][y]) / 2.0f;
                    const float dpdy = (pressure[x][ysouth] - pressure[x][ynorth]) / 2.0f;
                    const float latitude = latitudeforrow(y, height);
                    const float coriolis = std::clamp(std::fabs(latitude) / tuning::climate::atmosphere::coriolisLatitude, 0.0f, 1.0f);

                    const float directu = -dpdx * tuning::climate::atmosphere::directFlowFactor;
                    const float directv = -dpdy * tuning::climate::atmosphere::directFlowFactor;

                    float geou = 0.0f;
                    float geov = 0.0f;

                    if (latitude >= 0.0f)
                    {
                        geou = dpdy * tuning::climate::atmosphere::geostrophicFactor;
                        geov = -dpdx * tuning::climate::atmosphere::geostrophicFactor;
                    }
                    else
                    {
                        geou = -dpdy * tuning::climate::atmosphere::geostrophicFactor;
                        geov = dpdx * tuning::climate::atmosphere::geostrophicFactor;
                    }

                    float u = directu * (1.0f - coriolis) + geou * coriolis;
                    float v = directv * (1.0f - coriolis * 0.5f) + geov * coriolis;

                    windu[x][y] = std::clamp(u, -maxvectorwind, maxvectorwind);
                    windv[x][y] = std::clamp(v, -maxvectorwind, maxvectorwind);
                }
            }
        });

        smoothallfield(world, windu, tuning::climate::atmosphere::smoothingIterations);
        smoothallfield(world, windv, tuning::climate::atmosphere::smoothingIterations);
        applytopographicwindeffects(world, macroterrain, windu, windv);
        storeterrainverticalmotion(world, season, macroterrain, windu, windv);

        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int x = 0; x <= width; x++)
                {
                    world.setseasonaluwind(season, x, y, static_cast<int>(std::round(std::clamp(windu[x][y], -maxvectorwind, maxvectorwind))));
                    world.setseasonalvwind(season, x, y, static_cast<int>(std::round(std::clamp(windv[x][y], -maxvectorwind, maxvectorwind))));
                }
            }
        });
    }

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int y = startrow; y <= endrow; y++)
        {
            const float latitude = latitudeforrow(y, height);
            bool eastfallback = std::fabs(latitude) >= 22.0f && std::fabs(latitude) <= 60.0f;

            for (int x = 0; x <= width; x++)
            {
                float averageu = 0.0f;

                for (int season = 0; season < CLIMATESEASONCOUNT; season++)
                    averageu += static_cast<float>(world.seasonaluwind(season, x, y));

                averageu = averageu / static_cast<float>(CLIMATESEASONCOUNT);

                int scalarwind = static_cast<int>(std::round(averageu / tuning::climate::atmosphere::scalarWindDivisor));
                scalarwind = std::clamp(scalarwind, -10, 10);

                if (std::fabs(averageu) < tuning::climate::atmosphere::minimumScalarZonalWind || scalarwind == 0)
                {
                    if (averageu > 0.25f)
                        eastfallback = true;
                    else if (averageu < -0.25f)
                        eastfallback = false;

                    world.setwind(x, y, eastfallback ? 101 : 99);
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
    const int maxadvectiondistance = tuning::climate::moistureadvection::maxAdvectionDistance;
    const int iterations = tuning::climate::moistureadvection::iterations;

    auto runsolver = [&](int season)
    {
        floatgrid source(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid moisture(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid totalrain(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid uplift(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid descent(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid pressurelift(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid pressuresubsidence(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid available(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid nextmoisture(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid eastshare(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid westshare(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid northshare(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid southshare(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid convergence(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid totalconvergence(width + 1, vector<float>(height + 1, 0.0f));

        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                const int ynorth = (y > 0) ? y - 1 : y;
                const int ysouth = (y < height) ? y + 1 : y;

                for (int x = 0; x <= width; x++)
                {
                    const float u = static_cast<float>(world.seasonaluwind(season, x, y));
                    const float v = static_cast<float>(world.seasonalvwind(season, x, y));
                    const float magnitude = std::sqrt(u * u + v * v);
                    const float temperature = (world.sea(x, y) == 1) ? static_cast<float>(world.seasonalsst(season, x, y)) : static_cast<float>(world.seasonaltemp(season, x, y));

                    float evaporation = 0.0f;

                    if (world.sea(x, y) == 1)
                    {
                        evaporation = static_cast<float>(world.seasonalevaporation(season, x, y)) * tuning::climate::moistureadvection::sourceScale;

                        if (world.seaice(x, y) == 1)
                            evaporation = evaporation * tuning::climate::moistureadvection::seaIceFactor;

                        if (world.seaice(x, y) == 2)
                            evaporation = 0.0f;
                    }
                    else if (temperature > 0.0f)
                    {
                        evaporation = temperature * tuning::climate::moistureadvection::landEvaporationFactor;
                    }

                    source[x][y] = evaporation;

                    const int xwest = wrapx(x - 1, width);
                    const int xeast = wrapx(x + 1, width);
                    const float uwest = static_cast<float>(world.seasonaluwind(season, xwest, y));
                    const float ueast = static_cast<float>(world.seasonaluwind(season, xeast, y));
                    const float vnorth = static_cast<float>(world.seasonalvwind(season, x, ynorth));
                    const float vsouth = static_cast<float>(world.seasonalvwind(season, x, ysouth));
                    const float divergence = ((ueast - uwest) + (vsouth - vnorth)) / 8.0f;

                    convergence[x][y] = std::max(0.0f, -divergence);
                    pressurelift[x][y] = std::max(0.0f, -static_cast<float>(world.seasonalpressure(season, x, y)));
                    pressuresubsidence[x][y] = std::max(0.0f, static_cast<float>(world.seasonalpressure(season, x, y)));
                    uplift[x][y] = static_cast<float>(world.seasonaluplift(season, x, y)) / tuning::climate::atmosphere::topographyVerticalMotionStorageScale;
                    descent[x][y] = static_cast<float>(world.seasonalsubsidence(season, x, y)) / tuning::climate::atmosphere::topographyVerticalMotionStorageScale;
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
                        const float u = static_cast<float>(world.seasonaluwind(season, x, y));
                        const float v = static_cast<float>(world.seasonalvwind(season, x, y));
                        const float retained = std::max(0.0f, moisture[x][y] * tuning::climate::moistureadvection::carryRetention);
                        const float localavailable = retained + source[x][y];
                        const float zonalshare = std::clamp(std::fabs(u) * tuning::climate::moistureadvection::advectionDistanceScale / static_cast<float>(maxadvectiondistance), 0.0f, 1.0f);
                        const float meridionalshare = std::clamp(std::fabs(v) * tuning::climate::moistureadvection::advectionDistanceScale / static_cast<float>(maxadvectiondistance), 0.0f, 1.0f);
                        float east = (u > 0.0f) ? zonalshare : 0.0f;
                        float west = (u < 0.0f) ? zonalshare : 0.0f;
                        float south = (v > 0.0f) ? meridionalshare : 0.0f;
                        float north = (v < 0.0f) ? meridionalshare : 0.0f;
                        float totalshare = east + west + south + north;

                        if (totalshare > tuning::climate::moistureadvection::transportMaxFraction && totalshare > 0.0f)
                        {
                            const float scale = tuning::climate::moistureadvection::transportMaxFraction / totalshare;
                            east = east * scale;
                            west = west * scale;
                            south = south * scale;
                            north = north * scale;
                        }

                        available[x][y] = localavailable;
                        eastshare[x][y] = east;
                        westshare[x][y] = west;
                        northshare[x][y] = north;
                        southshare[x][y] = south;
                    }
                }
            });

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
                        const float localshare = 1.0f - eastshare[x][y] - westshare[x][y] - northshare[x][y] - southshare[x][y];

                        nextmoisture[x][y] =
                            available[x][y] * localshare +
                            available[xwest][y] * eastshare[xwest][y] +
                            available[xeast][y] * westshare[xeast][y] +
                            available[x][ynorth] * southshare[x][ynorth] +
                            available[x][ysouth] * northshare[x][ysouth];

                        const float moisturecapacity = tuning::climate::moistureadvection::moistureCapacityBase +
                            std::max(0.0f, static_cast<float>((world.sea(x, y) == 1) ? world.seasonalsst(season, x, y) : world.seasonaltemp(season, x, y)) + 8.0f) *
                            tuning::climate::moistureadvection::moistureCapacityTemperatureFactor;
                        const float safecapacity = std::max(5.0f, moisturecapacity);
                        const float convergenceamount = std::max(0.0f, nextmoisture[x][y] - available[x][y]);

                        convergence[x][y] = convergenceamount / safecapacity;
                        totalconvergence[x][y] += convergence[x][y];
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
                        const float u = static_cast<float>(world.seasonaluwind(season, x, y));
                        const float v = static_cast<float>(world.seasonalvwind(season, x, y));
                        const float magnitude = std::sqrt(u * u + v * v);
                        const float temperature = sea ? static_cast<float>(world.seasonalsst(season, x, y)) : static_cast<float>(world.seasonaltemp(season, x, y));
                        const float moisturecapacity = tuning::climate::moistureadvection::moistureCapacityBase +
                            std::max(0.0f, temperature + 8.0f) * tuning::climate::moistureadvection::moistureCapacityTemperatureFactor;
                        const float safecapacity = std::max(5.0f, moisturecapacity);
                        const float availablemoisture = std::max(0.0f, nextmoisture[x][y]);
                        const float humidityratio = (safecapacity > 0.0f) ? availablemoisture / safecapacity : 0.0f;
                        const float saturationfraction = std::clamp(
                            (humidityratio - tuning::climate::moistureadvection::saturationThreshold) /
                            std::max(0.01f, 1.0f - tuning::climate::moistureadvection::saturationThreshold),
                            0.0f, 1.0f);
                        float condensationrate = sea ? tuning::climate::moistureadvection::oceanCondensation : tuning::climate::moistureadvection::landCondensation;

                        if (sea && magnitude < 4.0f)
                            condensationrate += tuning::climate::moistureadvection::calmOceanBoost;

                        condensationrate += saturationfraction * tuning::climate::moistureadvection::humidityCondensationFactor;
                        condensationrate += convergence[x][y] * tuning::climate::moistureadvection::convergenceFactor;
                        condensationrate += pressurelift[x][y] * tuning::climate::moistureadvection::lowPressureFactor;
                        condensationrate += uplift[x][y] * tuning::climate::moistureadvection::upliftFactor;
                        condensationrate -= pressuresubsidence[x][y] * tuning::climate::moistureadvection::highPressureFactor;
                        condensationrate -= descent[x][y] * tuning::climate::moistureadvection::descentFactor;
                        condensationrate = condensationrate * std::max(tuning::climate::moistureadvection::minimumForcedHumidity, std::min(1.0f, humidityratio));
                        condensationrate = std::clamp(condensationrate, 0.0f, 0.92f);

                        float precipitation = availablemoisture * condensationrate;

                        if (availablemoisture > safecapacity)
                            precipitation += (availablemoisture - safecapacity) * tuning::climate::moistureadvection::overflowFactor;

                        precipitation = std::clamp(precipitation, 0.0f, availablemoisture);
                        totalrain[x][y] += precipitation;
                        moisture[x][y] = std::max(0.0f, availablemoisture - precipitation);
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
                    const int rain = static_cast<int>(std::round(totalrain[x][y] * tuning::climate::moistureadvection::rainfallScale));
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
    };

    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
        runsolver(season);
}
