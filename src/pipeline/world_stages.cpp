#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <queue>
#include <stdexcept>
#include <string>

#include "world_stages.hpp"
#include "climate_stages.hpp"

#include "classes.hpp"
#include "fastlem_mountains.hpp"
#include "functions.hpp"
#include "generation_tuning.hpp"
#include "physical_layers.hpp"
#include "world_generation_debug.hpp"
#include "world_processes.hpp"

using namespace std;

namespace
{
vector<vector<int>> makeintgrid(int value = 0)
{
    return vector<vector<int>>(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, value));
}

vector<vector<bool>> makeboolgrid(bool value = false)
{
    return vector<vector<bool>>(ARRAYWIDTH, vector<bool>(ARRAYHEIGHT, value));
}

vector<vector<std::uint8_t>> makefeaturegrid(std::uint8_t value = 0)
{
    return vector<vector<std::uint8_t>>(ARRAYWIDTH, vector<std::uint8_t>(ARRAYHEIGHT, value));
}

vector<vector<vector<int>>> makethreeintgrid()
{
    return vector<vector<vector<int>>>(ARRAYWIDTH, vector<vector<int>>(ARRAYHEIGHT, vector<int>(2, 0)));
}

void reseedterrainpass(planet& world, int salt)
{
    fast_srand(deterministicfastseed(deterministiccontextseed(world.seed(), salt)));
}

float clamp01(float value)
{
    return clamp(value, 0.0f, 1.0f);
}

int wrapx(const planet& world, int x)
{
    return wrap(x, world.width());
}

float normalizeterrainvalue(int current, int minimum, int maximum)
{
    return clamp01(static_cast<float>(current - minimum) / static_cast<float>(max(1, maximum - minimum)));
}

float samplecurveaxis(const vector<WorldCurvePoint>& points, float input, bool samplexfromy)
{
    if (points.empty())
        return clamp01(input);

    const float clampedinput = clamp01(input);

    for (size_t index = 1; index < points.size(); index++)
    {
        const float lowaxis = samplexfromy ? points[index - 1].y : points[index - 1].x;
        const float highaxis = samplexfromy ? points[index].y : points[index].x;

        if (clampedinput > highaxis)
            continue;

        const float span = max(0.0001f, highaxis - lowaxis);
        const float factor = clamp01((clampedinput - lowaxis) / span);
        const float lowvalue = samplexfromy ? points[index - 1].x : points[index - 1].y;
        const float highvalue = samplexfromy ? points[index].x : points[index].y;
        return clamp01(lowvalue + (highvalue - lowvalue) * factor);
    }

    return samplexfromy ? points.back().x : points.back().y;
}

float sampleheightcdf(const vector<float>& values, float normalizedheight)
{
    if (values.empty())
        return clamp01(normalizedheight);

    const float position = clamp01(normalizedheight) * static_cast<float>(values.size() - 1);
    const size_t lowindex = static_cast<size_t>(floorf(position));
    const size_t highindex = min(lowindex + 1, values.size() - 1);
    const float factor = position - static_cast<float>(lowindex);
    return clamp01(values[lowindex] + (values[highindex] - values[lowindex]) * factor);
}

void clearfeaturemask(GenerationScratch& scratch)
{
    scratch.stageFeatureMask = makefeaturegrid(0);
    scratch.stageFeatureCellCount = 0;
}

void markfeaturecell(GenerationScratch& scratch, int x, int y)
{
    if (x < 0 || x >= ARRAYWIDTH || y < 0 || y >= ARRAYHEIGHT)
        return;

    if (scratch.stageFeatureMask.empty())
        clearfeaturemask(scratch);

    if (scratch.stageFeatureMask[x][y] != 0)
        return;

    scratch.stageFeatureMask[x][y] = 1;
    scratch.stageFeatureCellCount++;
}

WorldTectonicsOptions tectonicoptions(const WorldStageParameters& parameters)
{
    WorldTectonicsOptions options;
    options.durationMyr = parameters.tectonicDurationMyr;
    options.plateCount = parameters.tectonicPlateCount;
    options.deltaTimeMyr = parameters.tectonicDeltaTimeMyr;
    options.scenario = parameters.tectonicScenario;
    options.worldGrid = parameters.tectonicWorldGrid;
    options.mineralHistory = parameters.tectonicMineralHistory;
    return options;
}

void rebuildshelves(planet& world, GenerationScratch& scratch)
{
    const int width = world.width();
    const int height = world.height();

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int y = startrow; y <= endrow; y++)
        {
            for (int x = 0; x <= width; x++)
                world.continentalshelves()[x][y] = false;
        }
    });

    makecontinentalshelves(world, world.continentalshelves(), 4);
}

void applyhypsometricremap(planet& world, GenerationScratch& scratch, const WorldStageParameters& parameters)
{
    const int width = world.width();
    const int height = world.height();
    const int minimum = 1;
    const int maximum = world.maxelevation() - 1;
    const vector<float> sourcecdf = scratch.heightCdf;
    const vector<WorldCurvePoint>& targetcurve = parameters.hypsometricCurve;

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int y = startrow; y <= endrow; y++)
        {
            for (int x = 0; x <= width; x++)
            {
                const float normalized = normalizeterrainvalue(world.nom(x, y), minimum, maximum);
                const float quantile = sampleheightcdf(sourcecdf, normalized);
                const float remapped = samplecurveaxis(targetcurve, quantile, false);
                const int newvalue = minimum + static_cast<int>(roundf(remapped * static_cast<float>(maximum - minimum)));
                world.setnom(x, y, clamp(newvalue, minimum, maximum));
            }
        }
    });

    rebuildheightdistribution(world, scratch);
}

void applyterraforming(planet& world, GenerationScratch& scratch, const WorldStageParameters& parameters)
{
    world.setsealevel(clamp(parameters.seaLevel, 1, world.maxelevation() - 2));

    if (parameters.useHypsometricRemap)
        applyhypsometricremap(world, scratch, parameters);

    rebuildshelves(world, scratch);
    rebuildheightdistribution(world, scratch);
}

vector<vector<int>> buildterrainfractal(planet& world, int salt, int maximum)
{
    vector<vector<int>> fractal = makeintgrid();
    reseedterrainpass(world, salt);
    createfractal(fractal, world.width(), world.height(), 8, 0.2f, static_cast<float>(random(3, 6)), 1, maximum, 0, 0);
    warp(fractal, world.width(), world.height(), world.maxelevation(), random(20, 80), 1);
    return fractal;
}

void metadataoceanrefinement(planet& world, GenerationScratch& scratch)
{
    const int width = world.width();
    const int height = world.height();
    const int maxelev = world.maxelevation();
    const int sealevel = world.sealevel();
    vector<vector<int>> seafractal = makeintgrid();

    reseedterrainpass(world, 0x4110);
    createfractal(seafractal, width, height, 8, 0.2f, static_cast<float>(random(3, 6)), 1, maxelev, 0, 0);
    warp(seafractal, width, height, maxelev, random(20, 80), 1);

    const float coastalvarreduce = static_cast<float>(maxelev) / 500.0f;
    const float oceanvarreduce = static_cast<float>(maxelev) / 1000.0f;

    for (int x = 0; x <= width; x++)
    {
        for (int y = 0; y <= height; y++)
        {
            if (world.sea(x, y) == 0)
                continue;

            if (world.continentalshelves()[x][y])
            {
                const float var = static_cast<float>(seafractal[x][y] - maxelev / 2) / coastalvarreduce;
                int newvalue = sealevel - 200 + static_cast<int>(var);
                newvalue = min(newvalue, sealevel - 10);
                world.setnom(x, y, max(1, newvalue));
            }
            else
            {
                int xx = x + width / 2;
                if (xx > width)
                    xx -= width;

                const float var = static_cast<float>(seafractal[xx][y] - maxelev / 2) / oceanvarreduce;
                int newvalue = sealevel - 5000 + static_cast<int>(var);
                newvalue = min(newvalue, sealevel - 3000);
                world.setnom(x, y, max(1, newvalue));
            }
        }
    }
}

int wrapdistx(const planet& world, int left, int right)
{
    const int width = world.width() + 1;
    int delta = abs(left - right);
    delta = min(delta, width - delta);
    return delta;
}

void metadataactivetrenches(planet& world, GenerationScratch& scratch)
{
    const int width = world.width();
    const int height = world.height();
    const int sealevel = world.sealevel();
    clearfeaturemask(scratch);
    reseedterrainpass(world, 0x4114);

    int maxsegmentid = 0;
    for (const TectonicBoundarySegment& segment : world.tectonicboundarysegments())
        maxsegmentid = max(maxsegmentid, segment.id);

    vector<int> segmentlookup(static_cast<size_t>(maxsegmentid + 1), -1);
    vector<int> segmentids;
    vector<vector<pair<int, int>>> segmentcells;

    for (int y = 0; y <= height; y++)
    {
        for (int x = 0; x <= width; x++)
        {
            if (world.sea(x, y) == 0 || world.tectonicboundarytype(x, y) != BoundaryType::convergent)
                continue;

            const CrustClass crustclass = world.tectoniccrustclass(x, y);
            if (crustclass != CrustClass::oceanic && crustclass != CrustClass::transitional)
                continue;

            const int segmentid = world.tectonicboundarysegmentid(x, y);
            if (segmentid <= 0)
                continue;

            if (segmentlookup[segmentid] < 0)
            {
                segmentlookup[segmentid] = static_cast<int>(segmentcells.size());
                segmentids.push_back(segmentid);
                segmentcells.emplace_back();
            }

            segmentcells[segmentlookup[segmentid]].emplace_back(x, y);
        }
    }

    vector<int> order(segmentids.size(), 0);
    for (size_t index = 0; index < order.size(); index++)
        order[index] = static_cast<int>(index);

    for (int index = static_cast<int>(order.size()) - 1; index > 0; index--)
        swap(order[index], order[random(0, index)]);

    const int trenchcount = min(static_cast<int>(order.size()), clamp(static_cast<int>(order.size()) / 3 + 1, 2, 6));

    for (int trenchindex = 0; trenchindex < trenchcount; trenchindex++)
    {
        const vector<pair<int, int>>& cells = segmentcells[order[trenchindex]];
        if (cells.empty())
            continue;

        const int branchcount = clamp(static_cast<int>(cells.size()) / 48, 1, 3);

        for (int branch = 0; branch < branchcount; branch++)
        {
            const pair<int, int>& anchor = cells[random(0, static_cast<int>(cells.size()) - 1)];
            const float maxsegmentdistance = static_cast<float>(random(18, 38));
            const int radius = random(2, 4);

            for (const pair<int, int>& cell : cells)
            {
                const int dx = wrapdistx(world, cell.first, anchor.first);
                const int dy = abs(cell.second - anchor.second);
                const float segmentdistance = sqrtf(static_cast<float>(dx * dx + dy * dy));

                if (segmentdistance > maxsegmentdistance)
                    continue;

                const float branchstrength = 1.0f - clamp01(segmentdistance / maxsegmentdistance);

                for (int ox = -radius; ox <= radius; ox++)
                {
                    for (int oy = -radius; oy <= radius; oy++)
                    {
                        const int nx = wrapx(world, cell.first + ox);
                        const int ny = cell.second + oy;

                        if (ny < 0 || ny > height || world.sea(nx, ny) == 0)
                            continue;

                        const float kerneldistance = sqrtf(static_cast<float>(ox * ox + oy * oy));
                        if (kerneldistance > static_cast<float>(radius) + 0.5f)
                            continue;

                        const float kernelfalloff = 1.0f - clamp01(kerneldistance / (static_cast<float>(radius) + 0.5f));
                        const float proximity = 1.0f - clamp01(static_cast<float>(world.tectonicboundarydistance(nx, ny)) / 16.0f);
                        const float subsidence = clamp01(world.tectonicsubsidencetendency(nx, ny));
                        const float convergence = clamp01(static_cast<float>(world.tectonicconvergence(nx, ny)) / 255.0f);
                        const float randomness = 0.82f + static_cast<float>(random(0, 35)) / 100.0f;
                        const int trenchdepth = static_cast<int>(roundf((1200.0f + 2200.0f * branchstrength) * kernelfalloff * randomness * (0.45f + subsidence * 0.30f + convergence * 0.25f)));

                        if (trenchdepth <= 0)
                            continue;

                        const int target = max(1, min(world.nom(nx, ny), sealevel - 450 - trenchdepth));
                        if (target >= world.nom(nx, ny))
                            continue;

                        world.setnom(nx, ny, target);
                        world.setgeologicregime(nx, ny, GeologicRegime::trench_adjacent);
                        markfeaturecell(scratch, nx, ny);
                    }
                }
            }
        }
    }
}

void applyvolcanicupliftkernel(planet& world, GenerationScratch& scratch, int centrex, int centrey, int radius, int maximumuplift, GeologicRegime regime)
{
    for (int dx = -radius; dx <= radius; dx++)
    {
        for (int dy = -radius; dy <= radius; dy++)
        {
            const int nx = wrapx(world, centrex + dx);
            const int ny = centrey + dy;

            if (ny < 0 || ny > world.height())
                continue;

            const float distance = sqrtf(static_cast<float>(dx * dx + dy * dy));
            if (distance > static_cast<float>(radius) + 0.5f)
                continue;

            const float falloff = 1.0f - clamp01(distance / (static_cast<float>(radius) + 0.5f));
            const int uplift = static_cast<int>(roundf(static_cast<float>(maximumuplift) * falloff * falloff));
            if (uplift <= 0)
                continue;

            const int target = min(world.maxelevation() - 1, world.nom(nx, ny) + uplift);
            if (target <= world.nom(nx, ny))
                continue;

            world.setnom(nx, ny, target);
            if (regime != GeologicRegime::stable)
                world.setgeologicregime(nx, ny, regime);
            markfeaturecell(scratch, nx, ny);
        }
    }
}

void metadatavolcanoes(planet& world, GenerationScratch& scratch)
{
    const int width = world.width();
    const int height = world.height();
    clearfeaturemask(scratch);
    rebuildshelves(world, scratch);
    reseedterrainpass(world, 0x4115);

    vector<vector<int>> volcanodirection = buildterrainfractal(world, 0x4116, world.maxelevation());
    vector<pair<int, int>> selectedseeds;

    auto farenough = [&](int x, int y, int minimumdistance)
    {
        for (const pair<int, int>& seed : selectedseeds)
        {
            const int dx = wrapdistx(world, x, seed.first);
            const int dy = abs(y - seed.second);
            if (dx * dx + dy * dy < minimumdistance * minimumdistance)
                return false;
        }

        return true;
    };

    struct VolcanoCandidate
    {
        int x = 0;
        int y = 0;
        float score = 0.0f;
    };

    vector<VolcanoCandidate> conecandidates;

    for (int y = 0; y <= height; y++)
    {
        for (int x = 0; x <= width; x++)
        {
            if (world.tectonicboundarytype(x, y) != BoundaryType::convergent)
                continue;

            const float proximity = 1.0f - clamp01(static_cast<float>(world.tectonicboundarydistance(x, y)) / 12.0f);
            const float uplift = clamp01(world.tectonicuplifttendency(x, y));
            const float convergence = clamp01(static_cast<float>(world.tectonicconvergence(x, y)) / 255.0f);
            const float history = clamp01(world.tectonicboundaryhistory(x, y));
            const float score = proximity * (0.35f + uplift * 0.35f + convergence * 0.20f + history * 0.10f);

            if (score < 0.34f)
                continue;

            conecandidates.push_back({ x, y, score });
        }
    }

    for (int index = static_cast<int>(conecandidates.size()) - 1; index > 0; index--)
        swap(conecandidates[index], conecandidates[random(0, index)]);

    const int area = (width + 1) * (height + 1);
    const int desiredcones = clamp(area / 200000, 2, 10);

    for (const VolcanoCandidate& candidate : conecandidates)
    {
        if (static_cast<int>(selectedseeds.size()) >= desiredcones)
            break;

        if (!farenough(candidate.x, candidate.y, 24))
            continue;

        const int peakheight = static_cast<int>(roundf(1800.0f + candidate.score * 3800.0f));
        applyvolcanicupliftkernel(world, scratch, candidate.x, candidate.y, world.sea(candidate.x, candidate.y) ? 2 : 1, static_cast<int>(roundf(240.0f + candidate.score * 520.0f)), GeologicRegime::convergent_arc);
        createisolatedvolcano(world, candidate.x, candidate.y, world.continentalshelves(), volcanodirection, peakheight, true);
        selectedseeds.emplace_back(candidate.x, candidate.y);
    }

    auto placenexthotspot = [&](bool oceanpreferred, int targetcount)
    {
        for (int placed = 0, attempts = 0; placed < targetcount && attempts < 512; attempts++)
        {
            const int x = random(0, width);
            const int y = random(0, height);

            if (!farenough(x, y, 34))
                continue;

            if (world.tectonicboundarydistance(x, y) < 18)
                continue;

            if (oceanpreferred != (world.sea(x, y) == 1))
                continue;

            if (world.sea(x, y) == 1 && world.continentalshelves()[x][y])
                continue;

            const int radius = oceanpreferred ? random(4, 6) : random(3, 5);
            const int uplift = oceanpreferred ? random(260, 700) : random(180, 460);
            const int peakheight = oceanpreferred ? random(900, 2400) : random(700, 1600);
            applyvolcanicupliftkernel(world, scratch, x, y, radius, uplift, GeologicRegime::stable);
            createisolatedvolcano(world, x, y, world.continentalshelves(), volcanodirection, peakheight, false);
            selectedseeds.emplace_back(x, y);
            placed++;
        }
    };

    placenexthotspot(true, clamp(area / 280000 + 1, 1, 6));
    placenexthotspot(false, clamp(area / 520000, 0, 3));

    for (int y = 0; y <= height; y++)
    {
        for (int x = 0; x <= width; x++)
        {
            if (world.volcano(x, y) == 0)
                continue;

            markfeaturecell(scratch, x, y);
            if (world.strato(x, y))
                world.setgeologicregime(x, y, GeologicRegime::convergent_arc);
        }
    }
}

vector<pair<int, int>> findcomponentdrainpath(planet& world, const vector<vector<int>>& componentids, int componentid)
{
    const int width = world.width();
    const int height = world.height();
    vector<vector<bool>> visited = makeboolgrid(false);
    vector<vector<int>> parentx = makeintgrid(-1);
    vector<vector<int>> parenty = makeintgrid(-1);
    queue<pair<int, int>> frontier;
    pair<int, int> endpoint = { -1, -1 };

    for (int y = 0; y <= height; y++)
    {
        for (int x = 0; x <= width; x++)
        {
            if (componentids[x][y] == componentid)
            {
                frontier.emplace(x, y);
                visited[x][y] = true;
                parentx[x][y] = x;
                parenty[x][y] = y;
            }
        }
    }

    const int directions[8][2] =
    {
        { 0, -1 }, { 1, -1 }, { 1, 0 }, { 1, 1 },
        { 0, 1 }, { -1, 1 }, { -1, 0 }, { -1, -1 }
    };

    while (frontier.empty() == false && endpoint.first < 0)
    {
        const pair<int, int> current = frontier.front();
        frontier.pop();

        for (const auto& direction : directions)
        {
            const int nx = wrapx(world, current.first + direction[0]);
            const int ny = current.second + direction[1];

            if (ny < 0 || ny > height || visited[nx][ny])
                continue;

            visited[nx][ny] = true;
            parentx[nx][ny] = current.first;
            parenty[nx][ny] = current.second;

            if (world.sea(nx, ny) == 1 && componentids[nx][ny] == 0)
            {
                endpoint = { nx, ny };
                break;
            }

            frontier.emplace(nx, ny);
        }
    }

    vector<pair<int, int>> path;

    if (endpoint.first < 0)
        return path;

    int x = endpoint.first;
    int y = endpoint.second;

    while (parentx[x][y] != x || parenty[x][y] != y)
    {
        path.emplace_back(x, y);
        const int px = parentx[x][y];
        const int py = parenty[x][y];
        x = px;
        y = py;
    }

    path.emplace_back(x, y);
    reverse(path.begin(), path.end());
    return path;
}

void applyinlandseapolicies(planet& world, GenerationScratch& scratch)
{
    if (scratch.inlandSeaComponents.empty())
        detectinlandseas(world, scratch);

    for (const GenerationWaterComponent& component : scratch.inlandSeaComponents)
    {
        if (component.policy == GenerationComponentPolicy::keep)
            continue;

        if (component.policy == GenerationComponentPolicy::fill)
        {
            for (int x = 0; x <= world.width(); x++)
            {
                for (int y = 0; y <= world.height(); y++)
                {
                    if (scratch.inlandSeaComponentIds[x][y] == component.id)
                        world.setnom(x, y, max(world.sealevel() + 1, world.nom(x, y)));
                }
            }
        }
        else
        {
            const vector<pair<int, int>> path = findcomponentdrainpath(world, scratch.inlandSeaComponentIds, component.id);

            for (const pair<int, int>& cell : path)
                world.setnom(cell.first, cell.second, max(1, world.sealevel() - 10));
        }
    }

    rebuildshelves(world, scratch);
    detectinlandseas(world, scratch);
}

void applybasinpolicies(planet& world, GenerationScratch& scratch)
{
    if (scratch.basinComponents.empty())
        detectbasincomponents(world, scratch);

    for (const GenerationBasinComponent& component : scratch.basinComponents)
    {
        if (component.policy == GenerationComponentPolicy::keep)
            continue;

        if (component.policy == GenerationComponentPolicy::fill)
        {
            int targetelevation = 0;

            for (int x = 0; x <= world.width(); x++)
            {
                for (int y = 0; y <= world.height(); y++)
                {
                    if (scratch.basinComponentIds[x][y] == component.id)
                        targetelevation = max(targetelevation, world.nom(x, y));
                }
            }

            for (int x = 0; x <= world.width(); x++)
            {
                for (int y = 0; y <= world.height(); y++)
                {
                    if (scratch.basinComponentIds[x][y] == component.id)
                        world.setnom(x, y, targetelevation);
                }
            }
        }
        else
        {
            const vector<pair<int, int>> path = findcomponentdrainpath(world, scratch.basinComponentIds, component.id);

            for (const pair<int, int>& cell : path)
                world.setnom(cell.first, cell.second, max(1, world.sealevel() - 10));
        }
    }

    detectbasincomponents(world, scratch);
}

void applyterraintexturing(planet& world, GenerationScratch& scratch)
{
    clearfeaturemask(scratch);
    vector<vector<int>> fractal = buildterrainfractal(world, 0x4111, world.maxelevation());
    flip(fractal, world.width(), world.height(), 1, 1);
    shift(fractal, world.width(), world.height(), random(1, world.width()));

    const int sealevel = world.sealevel();
    const int midpoint = world.maxelevation() / 2;

    for (int x = 0; x <= world.width(); x++)
    {
        for (int y = 0; y <= world.height(); y++)
        {
            const int oldvalue = world.nom(x, y);
            const int noise = fractal[x][y] - midpoint;
            const float strength = world.sea(x, y) ? 0.14f : 0.10f;
            int newvalue = oldvalue + static_cast<int>(roundf(static_cast<float>(noise) * strength));

            if (world.sea(x, y))
                newvalue = min(newvalue, sealevel - 10);
            else
                newvalue = max(newvalue, sealevel + 5);

            newvalue = clamp(newvalue, 1, world.maxelevation() - 1);

            if (newvalue == oldvalue)
                continue;

            world.setnom(x, y, newvalue);
            markfeaturecell(scratch, x, y);
        }
    }
}

bool applyfastlemmountains(planet& world)
{
    vector<vector<int>> fractal = buildterrainfractal(world, 0x4112, 12750);
    return generatefastlemmountains(world, fractal);
}

void rebuildterrainroughness(planet& world)
{
    vector<vector<int>> roughness = buildterrainfractal(world, 0x4113, world.maxelevation());

    for (int x = 0; x <= world.width(); x++)
    {
        for (int y = 0; y <= world.height(); y++)
            world.setroughness(x, y, static_cast<float>(roughness[x][y]));
    }
}

void applyriverandbasinstage(planet& world, GenerationScratch& scratch, const GenerationExecutionContext& context)
{
    if (context.dorivers == false)
        return;

    const int width = world.width();
    const int height = world.height();
    int saltlakesplaced = 0;
    scratch.saltLakeMap = makethreeintgrid();
    scratch.noLake = makeintgrid();
    scratch.basinSeeds = makeintgrid();

    createrivermap(world, scratch.mountaindrainage);

    if (context.dolakes || world.seatotal() == 0)
    {
        createsaltlakes(world, saltlakesplaced, scratch.saltLakeMap, scratch.noLake, scratch.basinSeeds, context.smalllake);
        addlandnoise(world);
        depressionfill(world);

        for (int x = 0; x <= width; x++)
        {
            for (int y = 0; y <= height; y++)
            {
                world.setriverdir(x, y, 0);
                world.setriverjan(x, y, 0);
                world.setriverjul(x, y, 0);
            }
        }

        createrivermap(world, scratch.mountaindrainage);
    }

    world.setmaxriverflow();
    detectbasincomponents(world, scratch);
}

void applylakestage(planet& world, GenerationScratch& scratch, const GenerationExecutionContext& context)
{
    if (context.dolakes == false)
    {
        if (world.seatotal() == 0)
            convertsaltlakes(world, scratch.saltLakeMap);

        return;
    }

    convertsaltlakes(world, scratch.saltLakeMap);
    createlakemap(world, scratch.noLake, context.smalllake, context.largelake);
    createriftlakemap(world, scratch.noLake);
    world.setmaxriverflow();
}

void applyaridfeatures(planet& world, const GenerationExecutionContext& context)
{
    createergs(world, context.smalllake, context.largelake, context.landshape);
    createsaltpans(world, context.smalllake, context.largelake);
}

void applydeltawetlandstage(planet& world, const GenerationExecutionContext& context)
{
    if (context.dodeltas)
    {
        createriverdeltas(world);
        checkrivers(world);
    }

    createwetlands(world, context.smalllake);
    removeexcesswetlands(world);
    refineroughnessmap(world);
    removesealakes(world);
    connectlakes(world);
    checkpoleclimates(world.climate());
    world.syncseasonalclimatefromlegacy();
    createclimatemap(world.climate());
    createbiomemap(world.climate());
    exportclimatevalidationreport(world.climate());
}

}

void detectinlandseas(planet& world, GenerationScratch& scratch)
{
    const int width = world.width();
    const int height = world.height();
    scratch.inlandSeaComponentIds = makeintgrid(0);
    scratch.inlandSeaComponents.clear();

    vector<vector<bool>> visited = makeboolgrid(false);
    vector<vector<pair<int, int>>> components(1);
    int largestcomponentid = 0;
    int largestcomponentarea = 0;

    for (int y = 0; y <= height; y++)
    {
        for (int x = 0; x <= width; x++)
        {
            if (visited[x][y] || world.sea(x, y) == 0)
                continue;

            queue<pair<int, int>> frontier;
            vector<pair<int, int>> cells;
            frontier.emplace(x, y);
            visited[x][y] = true;

            while (frontier.empty() == false)
            {
                const pair<int, int> current = frontier.front();
                frontier.pop();
                cells.push_back(current);

                const int directions[8][2] =
                {
                    { 0, -1 }, { 1, -1 }, { 1, 0 }, { 1, 1 },
                    { 0, 1 }, { -1, 1 }, { -1, 0 }, { -1, -1 }
                };

                for (const auto& direction : directions)
                {
                    const int nx = wrapx(world, current.first + direction[0]);
                    const int ny = current.second + direction[1];

                    if (ny < 0 || ny > height || visited[nx][ny] || world.sea(nx, ny) == 0)
                        continue;

                    visited[nx][ny] = true;
                    frontier.emplace(nx, ny);
                }
            }

            components.push_back(std::move(cells));

            if (static_cast<int>(components.back().size()) > largestcomponentarea)
            {
                largestcomponentarea = static_cast<int>(components.back().size());
                largestcomponentid = static_cast<int>(components.size() - 1);
            }
        }
    }

    for (size_t componentindex = 1; componentindex < components.size(); componentindex++)
    {
        if (static_cast<int>(componentindex) == largestcomponentid)
            continue;

        GenerationWaterComponent component;
        component.id = static_cast<int>(scratch.inlandSeaComponents.size() + 1);
        component.seedX = components[componentindex].front().first;
        component.seedY = components[componentindex].front().second;
        component.cellCount = static_cast<int>(components[componentindex].size());
        scratch.inlandSeaComponents.push_back(component);

        for (const pair<int, int>& cell : components[componentindex])
            scratch.inlandSeaComponentIds[cell.first][cell.second] = component.id;
    }
}

void detectbasincomponents(planet& world, GenerationScratch& scratch)
{
    const int width = world.width();
    const int height = world.height();
    scratch.basinComponentIds = makeintgrid(0);
    scratch.basinComponents.clear();
    vector<vector<bool>> visited = makeboolgrid(false);
    int nextid = 1;

    for (int y = 0; y <= height; y++)
    {
        for (int x = 0; x <= width; x++)
        {
            if (visited[x][y] || world.basinclass(x, y) != BasinClass::endorheic)
                continue;

            queue<pair<int, int>> frontier;
            frontier.emplace(x, y);
            visited[x][y] = true;
            GenerationBasinComponent component;
            component.id = nextid;
            component.seedX = x;
            component.seedY = y;
            component.basinClass = BasinClass::endorheic;

            while (frontier.empty() == false)
            {
                const pair<int, int> current = frontier.front();
                frontier.pop();
                component.cellCount++;
                scratch.basinComponentIds[current.first][current.second] = nextid;

                const int directions[8][2] =
                {
                    { 0, -1 }, { 1, -1 }, { 1, 0 }, { 1, 1 },
                    { 0, 1 }, { -1, 1 }, { -1, 0 }, { -1, -1 }
                };

                for (const auto& direction : directions)
                {
                    const int nx = wrapx(world, current.first + direction[0]);
                    const int ny = current.second + direction[1];

                    if (ny < 0 || ny > height || visited[nx][ny] || world.basinclass(nx, ny) != BasinClass::endorheic)
                        continue;

                    visited[nx][ny] = true;
                    frontier.emplace(nx, ny);
                }
            }

            scratch.basinComponents.push_back(component);
            nextid++;
        }
    }
}

bool executegenerationstage(GenerationStageId stageId, planet& world, GenerationScratch& scratch, WorldStageParameters& parameters, const GenerationExecutionContext& context, string* errormessage)
{
    const GenerationStageDefinition& stagedefinition = getgenerationstage(stageId);

    if (beginworldgenstep(stagedefinition.label) == false)
        return true;

    clearfeaturemask(scratch);
    seedprocessstream(world.seed(), stageId);
    ClimateStageContext climatecontext;
    climatecontext.smalllake = context.smalllake;
    climatecontext.landshape = context.landshape;

    switch (stageId)
    {
    case GenerationStageId::geosphere_initialization:
    {
        world.clear();
        scratch.geologySolver.set(initializegeosphere(world, world.continentalshelves(), tectonicoptions(parameters)));
        parameters.seaLevel = world.sealevel();
        rebuildheightdistribution(world, scratch);
        return true;
    }

    case GenerationStageId::plate_tectonics:
    {
        // Continues the initializer's engine; without one (a recipe or workbench
        // copy without that step), the world starts over and the engine is created here.
        std::shared_ptr<GeologySolverState> solver = scratch.geologySolver.take();
        if (solver == nullptr)
            world.clear();
        const WorldTectonicsOptions options = tectonicoptions(parameters);
        // Checkpoints are not part of what the solver state was made for.
        WorldTectonicsOptions evolveroptions = options;
        if (context.evolverCheckpoint)
        {
            evolveroptions.checkpointEveryMyr = parameters.tectonicCheckpointEveryMyr;
            evolveroptions.checkpointsMyr = parameters.tectonicCheckpointsMyr;
            evolveroptions.checkpoint = [&](double timemyr)
            {
                rebuildshelves(world, scratch);
                context.evolverCheckpoint(stageId, timemyr);
            };
        }
        evolvegeosphere(world, world.continentalshelves(), evolveroptions, std::move(solver));
        rebuildshelves(world, scratch);
        parameters.seaLevel = world.sealevel();
        rebuildheightdistribution(world, scratch);
        return true;
    }

    case GenerationStageId::terraforming:
        applyterraforming(world, scratch, parameters);
        return true;

    case GenerationStageId::inland_seas:
        applyinlandseapolicies(world, scratch);
        return true;

    case GenerationStageId::archipelagos:
    {
        vector<vector<bool>> removedland = makeboolgrid(false);
        makearchipelagos(world, removedland, context.landshape);
        return true;
    }

    case GenerationStageId::ocean_cleanup:
        widenchannels(world);
        loweroceans(world);
        rebuildshelves(world, scratch);
        return true;

    case GenerationStageId::coastline_refinement:
        removestraights(world);
        return true;

    case GenerationStageId::ocean_depth_refinement:
        metadataoceanrefinement(world, scratch);
        return true;

    case GenerationStageId::tectonic_trenches:
        metadataactivetrenches(world, scratch);
        return true;

    case GenerationStageId::tectonic_volcanoes:
        metadatavolcanoes(world, scratch);
        return true;

    case GenerationStageId::terrain_texturing:
        applyterraintexturing(world, scratch);
        return true;

    case GenerationStageId::fastlem_mountains:
        if (applyfastlemmountains(world) == false && errormessage != nullptr)
            *errormessage = "FastLEM mountains did not produce enough ridge candidates for this seed.";
        return true;

    case GenerationStageId::mountain_bases:
    {
        vector<vector<bool>> defaultmountains = makeboolgrid(false);
        vector<vector<bool>>& mountainsok = context.okmountains != nullptr ? *context.okmountains : defaultmountains;
        raisemountainbases(world, scratch.mountaindrainage, mountainsok);
        return true;
    }

    case GenerationStageId::terrain_smoothing:
        smoothland(world, 2);
        return true;

    case GenerationStageId::canyon_uplift:
        createextraelev(world);
        return true;

    case GenerationStageId::depression_fill:
        depressionfill(world);
        addlandnoise(world);
        depressionfill(world);
        return true;

    case GenerationStageId::coastline_adjustment:
        normalisecoasts(world, 13, 11, 4);
        normalisecoasts(world, 13, 11, 4);
        return true;

    case GenerationStageId::island_check:
        checkislands(world);
        extendnoshade(world);
        return true;

    case GenerationStageId::terrain_roughness:
        getlandandseatotals(world);
        rebuildterrainroughness(world);
        return true;

    case GenerationStageId::global_temperature:
    case GenerationStageId::ocean_currents:
    case GenerationStageId::sea_surface_temperatures:
    case GenerationStageId::pressure:
    case GenerationStageId::winds:
        executeclimatestage(stageId, world.climate(), climatecontext);
        return true;

    case GenerationStageId::sea_ice_and_tides:
        executeclimatestage(stageId, world.climate(), climatecontext);
        if (world.seatotal() > 0)
            createtidalmap(world);
        return true;

    case GenerationStageId::rainfall:
        // The planet keeps its own land/sea counts for the stages after climate.
        getlandandseatotals(world);
        executeclimatestage(stageId, world.climate(), climatecontext);
        return true;

    case GenerationStageId::fjords:
        if (world.seatotal() > 0)
            addfjordmountains(world);
        return true;

    case GenerationStageId::rivers_and_basins:
        applyriverandbasinstage(world, scratch, context);
        return true;

    case GenerationStageId::basin_editor:
        applybasinpolicies(world, scratch);
        if (context.dorivers)
            createrivermap(world, scratch.mountaindrainage);
        return true;

    case GenerationStageId::lakes:
        applylakestage(world, scratch, context);
        return true;

    case GenerationStageId::post_river_fastlem:
        broadenfastlemterrainfromrivers(world);
        return true;

    case GenerationStageId::mountain_temperature_lapse:
    case GenerationStageId::climates_and_biomes:
        executeclimatestage(stageId, world.climate(), climatecontext);
        return true;

    case GenerationStageId::arid_features:
        applyaridfeatures(world, context);
        return true;

    case GenerationStageId::deltas_wetlands_and_roughness:
        applydeltawetlandstage(world, context);
        return true;

    case GenerationStageId::finalize_layers:
        finalizegeneratedworld(world, scratch, context);
        return true;
    }

    return false;
}


bool stageusesmanualpreview(GenerationStageId stageid)
{
    return stageid == GenerationStageId::inland_seas || stageid == GenerationStageId::basin_editor;
}

void rebuildheightdistribution(const planet& world, GenerationScratch& scratch)
{
    const int width = world.width();
    const int height = world.height();
    const int minimum = 1;
    const int maximum = world.maxelevation() - 1;

    scratch.heightHistogram.assign(kHeightHistogramBins, 0);
    scratch.heightCdf.assign(kHeightHistogramBins, 0.0f);

    int total = 0;

    for (int x = 0; x <= width; x++)
    {
        for (int y = 0; y <= height; y++)
        {
            const float normalized = normalizeterrainvalue(world.nom(x, y), minimum, maximum);
            int bin = static_cast<int>(floorf(normalized * static_cast<float>(kHeightHistogramBins)));
            bin = clamp(bin, 0, kHeightHistogramBins - 1);
            scratch.heightHistogram[bin]++;
            total++;
        }
    }

    if (total <= 0)
        return;

    int running = 0;

    for (int index = 0; index < kHeightHistogramBins; index++)
    {
        running += scratch.heightHistogram[index];
        scratch.heightCdf[index] = static_cast<float>(running) / static_cast<float>(total);
    }
}

void initializegenerationscratch(GenerationScratch& scratch)
{
    scratch.mountaindrainage = makeintgrid();
    scratch.squareroot.assign((MAXCRATERRADIUS * MAXCRATERRADIUS + MAXCRATERRADIUS + 1) * 24, 0);

    for (size_t index = 1; index < scratch.squareroot.size(); index++)
        scratch.squareroot[index] = static_cast<int>(sqrt(static_cast<float>(index)));

    scratch.saltLakeMap = makethreeintgrid();
    scratch.noLake = makeintgrid();
    scratch.basinSeeds = makeintgrid();
    scratch.inlandSeaComponentIds = makeintgrid(0);
    scratch.basinComponentIds = makeintgrid(0);
    scratch.heightHistogram.assign(kHeightHistogramBins, 0);
    scratch.heightCdf.assign(kHeightHistogramBins, 0.0f);
    scratch.inlandSeaComponents.clear();
    scratch.basinComponents.clear();
    clearfeaturemask(scratch);
}

void finalizegeneratedworld(planet& world, GenerationScratch& scratch, const GenerationExecutionContext& context)
{
    generatephysicalworldlayers(world, world.continentalshelves());
    generatesocialworld(world, context.socialoptions);

    if (context.appendclimateworkbook && recordclimatebenchmarkrun(
            world.climate(),
            "",
            climatebenchmarkmapselection::defaultselection()) == false)
        updatereport("Climate benchmark artifact update failed");
}
