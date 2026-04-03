#include <array>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

#include "third_party/plate_tectonics/src/platecapi.hpp"
#ifdef max
#undef max
#endif
#ifdef min
#undef min
#endif
#include "functions.hpp"
#include "generation_tuning.hpp"
#include "planet.hpp"

namespace
{
struct PlateTectonicsHandle
{
    void* pointer = nullptr;

    ~PlateTectonicsHandle()
    {
        if (pointer != nullptr)
            platec_api_destroy(pointer);
    }
};

struct PlateTectonicsSnapshot
{
    int width = 0;
    int height = 0;
    std::vector<uint32_t> platesmap;
    std::vector<std::array<float, 2>> velocities;
    bool valid = false;
};

int snapshotindex(int x, int y, int width)
{
    return y * width + x;
}

float clamp01(float value)
{
    return std::clamp(value, 0.0f, 1.0f);
}

float currentsearatio(planet& world)
{
    const int width = world.width();
    const int height = world.height();
    int seacells = 0;
    const int totalcells = (width + 1) * (height + 1);

    for (int y = 0; y <= height; y++)
    {
        for (int x = 0; x <= width; x++)
            seacells += world.sea(x, y) ? 1 : 0;
    }

    return static_cast<float>(seacells) / static_cast<float>(std::max(1, totalcells));
}

void normalizeheightmap(std::vector<float>& heightmap)
{
    if (heightmap.empty())
        return;

    const auto bounds = std::minmax_element(heightmap.begin(), heightmap.end());
    const float lowest = *bounds.first;
    const float highest = *bounds.second;

    if (highest <= lowest)
    {
        std::fill(heightmap.begin(), heightmap.end(), 0.0f);
        return;
    }

    const float inverserange = 1.0f / (highest - lowest);

    for (float& value : heightmap)
        value = (value - lowest) * inverserange;
}

float findseathreshold(const std::vector<float>& heightmap, float searatio)
{
    float threshold = 0.5f;
    float step = 0.5f;
    const std::size_t cellcount = heightmap.size();

    while (step > 0.0005f)
    {
        std::size_t seacells = 0;

        for (float value : heightmap)
            seacells += value < threshold ? 1U : 0U;

        step *= 0.5f;

        if (static_cast<float>(seacells) / static_cast<float>(cellcount) < searatio)
            threshold += step;
        else
            threshold -= step;
    }

    return clamp01(threshold);
}

void captureplatetectonicssnapshot(void* simulation, int width, int height, PlateTectonicsSnapshot& snapshot)
{
    const uint32_t platecount = platec_api_get_plate_count(simulation);
    const uint32_t* platesmap = platec_api_get_platesmap(simulation);

    if (platecount == 0 || platesmap == nullptr)
    {
        snapshot.valid = false;
        return;
    }

    const std::size_t cellcount = static_cast<std::size_t>(width) * static_cast<std::size_t>(height);
    snapshot.width = width;
    snapshot.height = height;
    snapshot.platesmap.assign(platesmap, platesmap + cellcount);
    snapshot.velocities.assign(platecount, { 0.0f, 0.0f });

    for (uint32_t plateindex = 0; plateindex < platecount; plateindex++)
    {
        snapshot.velocities[plateindex][0] = platec_api_velocity_unity_vector_x(simulation, plateindex);
        snapshot.velocities[plateindex][1] = platec_api_velocity_unity_vector_y(simulation, plateindex);
    }

    snapshot.valid = true;
}

int normalizetectonicsignal(float signal)
{
    const float threshold = tuning::terrain::platetectonics::convergentBoundaryThreshold;
    const float normalized = clamp01((signal - threshold) / (2.0f - threshold));
    return static_cast<int>(std::round(normalized * 100.0f));
}

bool hasadjacentland(const planet& world, int x, int y)
{
    const int width = world.width();
    const int height = world.height();

    for (int dy = -1; dy <= 1; dy++)
    {
        for (int dx = -1; dx <= 1; dx++)
        {
            if (dx == 0 && dy == 0)
                continue;

            const int ny = y + dy;

            if (ny < 0 || ny > height)
                continue;

            const int nx = wrap(x + dx, width);

            if (world.sea(nx, ny) == 0)
                return true;
        }
    }

    return false;
}

bool hasadjacentsea(const planet& world, int x, int y)
{
    const int width = world.width();
    const int height = world.height();

    for (int dy = -1; dy <= 1; dy++)
    {
        for (int dx = -1; dx <= 1; dx++)
        {
            if (dx == 0 && dy == 0)
                continue;

            const int ny = y + dy;

            if (ny < 0 || ny > height)
                continue;

            const int nx = wrap(x + dx, width);

            if (world.sea(nx, ny) != 0)
                return true;
        }
    }

    return false;
}

GeologicRegime classifygeologicregime(const planet& world, int x, int y, int convergence, int divergence, int shear)
{
    const bool sea = world.sea(x, y) != 0;
    const bool adjacentland = hasadjacentland(world, x, y);
    const bool adjacentsea = hasadjacentsea(world, x, y);
    const int maxsignal = std::max({ convergence, divergence, shear });

    if (convergence >= 65)
    {
        if (sea && adjacentland)
            return GeologicRegime::trench_adjacent;

        if (!sea && adjacentsea)
            return GeologicRegime::convergent_arc;

        return GeologicRegime::continent_collision;
    }

    if (divergence >= 55)
        return sea ? GeologicRegime::mid_ocean_ridge : GeologicRegime::divergent_rift;

    if (shear >= 50 && shear >= convergence + 10 && shear >= divergence + 10)
        return GeologicRegime::transform;

    if ((world.outline(x, y) || world.coast(x, y)) && maxsignal <= 20)
        return GeologicRegime::passive_margin;

    return GeologicRegime::stable;
}

void recordboundarymotion(const PlateTectonicsSnapshot& snapshot, int x, int y, int nx, int ny, std::vector<float>& convergence, std::vector<float>& divergence, std::vector<float>& shear)
{
    if (nx < 0 || nx >= snapshot.width || ny < 0 || ny >= snapshot.height)
        return;

    const int index = snapshotindex(x, y, snapshot.width);
    const int neighbourindex = snapshotindex(nx, ny, snapshot.width);
    const uint32_t plate = snapshot.platesmap[index];
    const uint32_t neighbourplate = snapshot.platesmap[neighbourindex];

    if (plate == neighbourplate || plate >= snapshot.velocities.size() || neighbourplate >= snapshot.velocities.size())
        return;

    float normalx = static_cast<float>(nx - x);
    float normaly = static_cast<float>(ny - y);
    const float length = std::sqrt(normalx * normalx + normaly * normaly);

    if (length <= 0.0f)
        return;

    normalx /= length;
    normaly /= length;

    const float relativex = snapshot.velocities[plate][0] - snapshot.velocities[neighbourplate][0];
    const float relativey = snapshot.velocities[plate][1] - snapshot.velocities[neighbourplate][1];
    const float boundarymotion = relativex * normalx + relativey * normaly;
    const float tangentialmotion = std::abs(relativex * (-normaly) + relativey * normalx);

    if (boundarymotion > 0.0f)
    {
        convergence[index] = std::max(convergence[index], boundarymotion);
        convergence[neighbourindex] = std::max(convergence[neighbourindex], boundarymotion);
    }
    else if (boundarymotion < 0.0f)
    {
        const float divergencesignal = -boundarymotion;
        divergence[index] = std::max(divergence[index], divergencesignal);
        divergence[neighbourindex] = std::max(divergence[neighbourindex], divergencesignal);
    }

    shear[index] = std::max(shear[index], tangentialmotion);
    shear[neighbourindex] = std::max(shear[neighbourindex], tangentialmotion);
}
}

void applyplatetectonicssimulation(planet& world, std::vector<std::vector<bool>>& shelves)
{
    const int width = world.width();
    const int height = world.height();
    const int simwidth = width + 1;
    const int simheight = height + 1;
    const std::size_t cellcount = static_cast<std::size_t>(simwidth) * static_cast<std::size_t>(simheight);
    const int sealevel = world.sealevel();
    const int maxelev = world.maxelevation();
    const float searatio = std::clamp(currentsearatio(world) + tuning::terrain::platetectonics::seaLevelBias, 0.05f, 0.95f);

    world.cleartectonicprovenance();

    std::vector<float> inputheightmap(cellcount, 0.0f);
    std::vector<bool> originalsea(cellcount, false);
    int lowest = world.nom(0, 0);
    int highest = world.nom(0, 0);

    for (int y = 0; y <= height; y++)
    {
        for (int x = 0; x <= width; x++)
        {
            const int value = world.nom(x, y);
            lowest = std::min(lowest, value);
            highest = std::max(highest, value);
        }
    }

    const float inverserange = highest > lowest ? 1.0f / static_cast<float>(highest - lowest) : 0.0f;

    for (int y = 0; y <= height; y++)
    {
        for (int x = 0; x <= width; x++)
        {
            const std::size_t index = static_cast<std::size_t>(y) * static_cast<std::size_t>(simwidth) + static_cast<std::size_t>(x);
            const int current = world.nom(x, y);
            inputheightmap[index] = inverserange > 0.0f ? static_cast<float>(current - lowest) * inverserange : 0.0f;
            originalsea[index] = world.sea(x, y) != 0;
        }
    }

    PlateTectonicsHandle simulation;
    simulation.pointer = platec_api_create_from_heightmap(world.seed(), static_cast<uint32_t>(simwidth), static_cast<uint32_t>(simheight),
        inputheightmap.data(), searatio,
        tuning::terrain::platetectonics::erosionPeriod,
        tuning::terrain::platetectonics::foldingRatio,
        tuning::terrain::platetectonics::aggregationOverlapAbsolute,
        tuning::terrain::platetectonics::aggregationOverlapRelative,
        static_cast<uint32_t>(platetectonicscyclecount()),
        tuning::terrain::platetectonics::plateCount);

    if (simulation.pointer == nullptr)
        return;

    PlateTectonicsSnapshot snapshot;

    for (uint32_t step = 0; step < tuning::terrain::platetectonics::maximumSimulationSteps; step++)
    {
        captureplatetectonicssnapshot(simulation.pointer, simwidth, simheight, snapshot);

        if (platec_api_is_finished(simulation.pointer) != 0)
            break;

        platec_api_step(simulation.pointer);
    }

    float* outputheightmap = platec_api_get_heightmap(simulation.pointer);

    if (outputheightmap == nullptr)
        return;

    std::vector<float> transformed(outputheightmap, outputheightmap + cellcount);
    normalizeheightmap(transformed);

    if (tuning::terrain::platetectonics::outputBlend < 1.0f)
    {
        const float blend = clamp01(tuning::terrain::platetectonics::outputBlend);

        for (std::size_t index = 0; index < cellcount; index++)
            transformed[index] = inputheightmap[index] * (1.0f - blend) + transformed[index] * blend;
    }

    const float landbiasedsearatio = clamp01(searatio - tuning::terrain::platetectonics::landRetentionSeaBias);
    const float seathreshold = findseathreshold(transformed, landbiasedsearatio);
    const float searange = std::max(0.0001f, seathreshold);
    const float landrange = std::max(0.0001f, 1.0f - seathreshold);
    const int oceanceiling = std::max(tuning::terrain::platetectonics::minimumOceanDepth, sealevel - tuning::terrain::platetectonics::coastalOceanOffset);
    const int landfloor = std::min(maxelev - 1, sealevel + tuning::terrain::platetectonics::landStartOffset);

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int y = startrow; y <= endrow; y++)
        {
            for (int x = 0; x <= width; x++)
            {
                const std::size_t index = static_cast<std::size_t>(y) * static_cast<std::size_t>(simwidth) + static_cast<std::size_t>(x);
                const float value = transformed[index];

                if (originalsea[index])
                {
                    const float oceanrelative = std::pow(clamp01(value / searange), tuning::terrain::platetectonics::oceanExponent);
                    const int newvalue = tuning::terrain::platetectonics::minimumOceanDepth +
                        static_cast<int>(std::round(static_cast<float>(oceanceiling - tuning::terrain::platetectonics::minimumOceanDepth) * oceanrelative));
                    world.setnom(x, y, std::max(tuning::terrain::platetectonics::minimumOceanDepth, std::min(oceanceiling, newvalue)));
                }
                else if (value <= seathreshold)
                {
                    const float oceanrelative = std::pow(clamp01(value / searange), tuning::terrain::platetectonics::oceanExponent);
                    const int newvalue = tuning::terrain::platetectonics::minimumOceanDepth +
                        static_cast<int>(std::round(static_cast<float>(oceanceiling - tuning::terrain::platetectonics::minimumOceanDepth) * oceanrelative));
                    world.setnom(x, y, std::max(tuning::terrain::platetectonics::minimumOceanDepth, std::min(oceanceiling, newvalue)));
                }
                else
                {
                    const float landrelative = std::pow(clamp01((value - seathreshold) / landrange), tuning::terrain::platetectonics::landExponent);
                    const int newvalue = landfloor +
                        static_cast<int>(std::round(static_cast<float>((maxelev - 1) - landfloor) * landrelative));
                    world.setnom(x, y, std::max(landfloor, std::min(maxelev - 1, newvalue)));
                }

                shelves[x][y] = false;
            }
        }
    });

    if (snapshot.valid)
    {
        std::vector<float> convergence(cellcount, 0.0f);
        std::vector<float> divergence(cellcount, 0.0f);
        std::vector<float> shear(cellcount, 0.0f);

        for (int y = 0; y < simheight; y++)
        {
            for (int x = 0; x < simwidth; x++)
            {
                if (x + 1 < simwidth)
                    recordboundarymotion(snapshot, x, y, x + 1, y, convergence, divergence, shear);

                if (y + 1 < simheight)
                    recordboundarymotion(snapshot, x, y, x, y + 1, convergence, divergence, shear);
            }
        }

        std::vector<std::vector<int>> rawmountains(ARRAYWIDTH, std::vector<int>(ARRAYHEIGHT, 0));
        bool havemountains = false;

        for (int y = 0; y <= height; y++)
        {
            for (int x = 0; x <= width; x++)
            {
                const std::size_t index = static_cast<std::size_t>(y) * static_cast<std::size_t>(simwidth) + static_cast<std::size_t>(x);
                const int convergencescore = normalizetectonicsignal(convergence[index]);
                const int divergencescore = normalizetectonicsignal(divergence[index]);
                const int shearscore = normalizetectonicsignal(shear[index]);

                world.settectonicconvergence(x, y, convergencescore);
                world.settectonicdivergence(x, y, divergencescore);
                world.settectonicshear(x, y, shearscore);
                world.setgeologicregime(x, y, classifygeologicregime(world, x, y, convergencescore, divergencescore, shearscore));

                if (originalsea[index] || world.nom(x, y) <= sealevel + 25)
                    continue;

                const float convergencesignal = clamp01((convergence[index] - tuning::terrain::platetectonics::convergentBoundaryThreshold) /
                    (2.0f - tuning::terrain::platetectonics::convergentBoundaryThreshold));

                if (convergencesignal <= 0.0f)
                    continue;

                const float mountainsignal = std::pow(convergencesignal, 0.8f);
                const int uplift = static_cast<int>(std::round(static_cast<float>(tuning::terrain::platetectonics::collisionUplift) * std::sqrt(convergencesignal)));
                const int peakheight = tuning::terrain::platetectonics::collisionMinimumPeak +
                    static_cast<int>(std::round(static_cast<float>(tuning::terrain::platetectonics::collisionMaximumPeak - tuning::terrain::platetectonics::collisionMinimumPeak) * mountainsignal));

                if (uplift > 0)
                    world.setnom(x, y, std::min(maxelev - 1, world.nom(x, y) + uplift));

                rawmountains[x][y] = peakheight;
                havemountains = true;
            }
        }

        if (havemountains)
        {
            std::vector<std::vector<bool>> dummyok(ARRAYWIDTH, std::vector<bool>(ARRAYHEIGHT, false));
            createmountainsfromraw(world, rawmountains, dummyok);
        }
    }
}
