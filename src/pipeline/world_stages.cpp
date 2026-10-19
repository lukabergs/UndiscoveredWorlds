#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <queue>
#include <stdexcept>
#include <string>

#include "pipeline/world_stages.hpp"
#include "pipeline/climate_stages.hpp"

#include "core/classes.hpp"
#include "core/grid_shapes.hpp"
#include "core/mesh/cell_grid.hpp"
#include "world/functions.hpp"
#include "pipeline/physical_layers.hpp"
#include "pipeline/world_generation_debug.hpp"
#include "pipeline/world_processes.hpp"
#include "simulations/hydrology/flow_routing.hpp"
#include "simulations/hydrology/lake_levels.hpp"
#include "simulations/hydrology/rivers.hpp"
#include "core/deterministic_random.hpp"
#include "simulations/terrain/fluvial_erosion.hpp"
#include "simulations/terrain/terrain_elevation.hpp"
#include "simulations/terrain/terrain_shapes.hpp"
#include "simulations/terrain/terrain_volcanoes.hpp"
#include "simulations/terrain/subgrid_statistics.hpp"
#include "simulations/terrain/terrain_parameters.hpp"
#include "world/world_parameters.hpp"

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

void reseedterrainpass(planet& world, int salt)
{
    fast_srand(deterministicfastseed(deterministiccontextseed(world.seed(), salt)));
}

float clamp01(float value)
{
    return clamp(value, 0.0f, 1.0f);
}

float normalizeterrainvalue(float current, int minimum, int maximum)
{
    return clamp01((current - static_cast<float>(minimum)) / static_cast<float>(max(1, maximum - minimum)));
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

// The scratch's [x][y] arrays hold a cell at its legacy slot (core/grid_shapes.hpp).
void markfeaturecell(GenerationScratch& scratch, const planet& world, mesh::CellId c)
{
    int x = 0;
    int y = 0;
    legacyslot(world.grid(), c, x, y);
    markfeaturecell(scratch, x, y);
}

template <typename T>
T& scratchcell(std::vector<std::vector<T>>& grid, const planet& world, mesh::CellId c)
{
    int x = 0;
    int y = 0;
    legacyslot(world.grid(), c, x, y);
    return grid[x][y];
}

template <typename T>
T scratchcell(const std::vector<std::vector<T>>& grid, const planet& world, mesh::CellId c)
{
    int x = 0;
    int y = 0;
    legacyslot(world.grid(), c, x, y);
    return grid[x][y];
}

mesh::CellId cellcount(const planet& world)
{
    return static_cast<mesh::CellId>(world.cellcount());
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
    world.continentalshelves().fill(false);

    makecontinentalshelves(world, world.continentalshelves(), 4);
}

void applyhypsometricremap(planet& world, GenerationScratch& scratch, const WorldStageParameters& parameters)
{
    const int minimum = 1;
    const int maximum = world.maxelevation() - 1;
    const vector<float> sourcecdf = scratch.heightCdf;
    const vector<WorldCurvePoint>& targetcurve = parameters.hypsometricCurve;

    parallelforrows(0, static_cast<int>(cellcount(world)) - 1, [&](int first, int last)
    {
        for (mesh::CellId c = static_cast<mesh::CellId>(first); c <= static_cast<mesh::CellId>(last); c++)
        {
            const float normalized = normalizeterrainvalue(world.nomfloat(c), minimum, maximum);
            const float quantile = sampleheightcdf(sourcecdf, normalized);
            const float remapped = samplecurveaxis(targetcurve, quantile, false);
            const float newvalue = static_cast<float>(minimum) + remapped * static_cast<float>(maximum - minimum);
            world.setnomfloat(c, clamp(newvalue, static_cast<float>(minimum), static_cast<float>(maximum)));
        }
    }, 1024);

    rebuildheightdistribution(world, scratch);
}

// The whole-metre sea level that leaves the area share `landfraction` of the surface above it
// (base elevation; a cell weighs its area, mesh::CellGrid::area_m2). No connectivity rule: a
// closed basin below the level is below sea like any other cell.
int sealevelforlandfraction(const planet& world, float landfraction)
{
    const mesh::CellGrid& grid = world.grid();

    vector<pair<int, double>> cells;   // whole metres, as basebelowsea compares them; area weight
    cells.reserve(world.cellcount());
    double totalweight = 0.0;

    for (mesh::CellId c = 0; c < cellcount(world); c++)
    {
        const double weight = grid.area_m2(c) / grid.totalarea_m2();
        cells.push_back({ wholemetres(world.nomfloat(c)), weight });
        totalweight += weight;
    }

    if (cells.empty() || totalweight <= 0.0)
        return world.sealevel();

    // Highest first; a level just below a group of equal heights makes that group land. The level
    // whose land share is nearest the target wins (on a tie the higher one).
    stable_sort(cells.begin(), cells.end(), [](const pair<int, double>& a, const pair<int, double>& b) { return a.first > b.first; });

    const double target = static_cast<double>(landfraction) * totalweight;
    double cumulative = 0.0;
    int best = cells.front().first;
    double besterror = target;
    size_t first = 0;

    while (first < cells.size())
    {
        size_t last = first;
        while (last < cells.size() && cells[last].first == cells[first].first)
        {
            cumulative += cells[last].second;
            last++;
        }

        const double error = abs(cumulative - target);
        if (error < besterror)
        {
            best = cells[first].first - 1;
            besterror = error;
        }
        first = last;
    }

    return best;
}

void applyterraforming(planet& world, GenerationScratch& scratch, const WorldStageParameters& parameters)
{
    // world.ocean.land_fraction in (0, 1) places the sea level at that land share of the surface,
    // unless the stage was asked for another level than the one the tectonics left (the
    // workbench's sea-level slider); 0 keeps the requested level (before 2026-10-16 always the
    // tectonic engine's). With world.ocean.land_fraction_after_coast the fraction is placed by
    // terrain.sea_level instead, and this stage keeps the requested level.
    const float landfraction = tuning::world::ocean::landFraction;
    int sealevel = parameters.seaLevel;
    if (landfraction > 0.0f && landfraction < 1.0f && !tuning::world::ocean::landFractionAfterCoast &&
        parameters.seaLevel == world.sealevel())
        sealevel = sealevelforlandfraction(world, landfraction);

    world.setsealevel(clamp(sealevel, 1, world.maxelevation() - 2));

    if (parameters.useHypsometricRemap)
        applyhypsometricremap(world, scratch, parameters);

    rebuildshelves(world, scratch);
    rebuildheightdistribution(world, scratch);
}

// terrain.sea_level: world.ocean.land_fraction in (0, 1) on the finished coast (after the coast
// passes, before the island check, which reads the level), when
// world.ocean.land_fraction_after_coast; otherwise nothing.
void applysealevel(planet& world, GenerationScratch& scratch)
{
    const float landfraction = tuning::world::ocean::landFraction;
    if (!tuning::world::ocean::landFractionAfterCoast || !(landfraction > 0.0f && landfraction < 1.0f))
        return;

    world.setsealevel(clamp(sealevelforlandfraction(world, landfraction), 1, world.maxelevation() - 2));
    rebuildshelves(world, scratch);
    rebuildheightdistribution(world, scratch);
}

// The terrain stages' fractal on the sphere (terrain.sphere_noise.terrain_fractal_*, warped as
// the legacy grain-8 fractal was): one value per cell in [1, maximum]. The pass's stream is
// reseeded with the salt first, as before, and the warp's reach is drawn from it.
vector<float> buildterrainfractal(planet& world, int salt, int maximum)
{
    reseedterrainpass(world, salt);
    const int warpfactor = random(20, 80);
    return terrain::spherefractalwarped(world, deterministiccontextseed(world.seed(), salt), terrain::terrainfractal(),
                                        static_cast<float>(maximum), warpfactor);
}

// The ocean's depth curve (tuning::terrain::ocean_floor): depth in metres at an area share
// 0..1 of the ocean, linear within each band.
float oceandepthatshare(float share)
{
    namespace of = tuning::terrain::ocean_floor;
    const array<float, 9> edges = { 0.0f, 200.0f, 1000.0f, 2000.0f, 3000.0f, 4000.0f, 5000.0f, 6000.0f, of::deepestM };
    const array<float, 8> shares = { of::share0200M, of::share2001000M, of::share10002000M, of::share20003000M,
                                     of::share30004000M, of::share40005000M, of::share50006000M, of::shareOver6000M };
    float total = 0.0f;
    for (const float value : shares)
        total += max(0.0f, value);

    if (total <= 0.0f)
        return edges.back();

    float lower = 0.0f;
    for (size_t band = 0; band < shares.size(); band++)
    {
        const float bandshare = max(0.0f, shares[band]) / total;
        if (bandshare > 0.0f && (share < lower + bandshare || band == shares.size() - 1))
            return edges[band] + clamp01((share - lower) / bandshare) * (edges[band + 1] - edges[band]);

        lower += bandshare;
    }

    return edges.back();
}

// Every sea cell takes the depth of the ocean curve at its area-weighted rank of the engine's
// depth: the engine's structure (shallow margins and ridges, deep old basins) stays, and the
// depth distribution no longer depends on the grid width or on the shelf mask.
void metadataoceanrefinement(planet& world, GenerationScratch&)
{
    const mesh::CellGrid& grid = world.grid();
    const float sealevel = static_cast<float>(world.sealevel());

    struct SeaCell
    {
        float elevation;
        double weight;
        mesh::CellId cell;
    };

    vector<SeaCell> cells;
    double totalweight = 0.0;

    for (mesh::CellId c = 0; c < cellcount(world); c++)
    {
        if (world.sea(c) == 0)
            continue;

        // A cell weighs its area.
        const double weight = grid.area_m2(c) / grid.totalarea_m2();
        cells.push_back({ world.nomfloat(c), weight, c });
        totalweight += weight;
    }

    if (cells.empty() || totalweight <= 0.0)
        return;

    // Shallowest first.
    stable_sort(cells.begin(), cells.end(), [](const SeaCell& a, const SeaCell& b) { return a.elevation > b.elevation; });

    const float minimumdepth = tuning::terrain::ocean_floor::minimumDepthM;
    double cumulative = 0.0;
    size_t first = 0;

    while (first < cells.size())
    {
        // Cells of equal engine height form one group and all take the share at the group's
        // middle, so a tie is not broken by the order of the cells.
        size_t last = first;
        double groupweight = 0.0;
        while (last < cells.size() && cells[last].elevation == cells[first].elevation)
        {
            groupweight += cells[last].weight;
            last++;
        }

        const float share = static_cast<float>((cumulative + 0.5 * groupweight) / totalweight);
        const float depth = max(oceandepthatshare(share), minimumdepth);

        for (size_t index = first; index < last; index++)
            world.setnomfloat(cells[index].cell, max(1.0f, sealevel - depth));

        cumulative += groupweight;
        first = last;
    }
}

// Trenches along convergent oceanic boundary segments: a branch reaches up to 18 to 38 legacy
// pixels (terrain.sphere_noise.legacy_columns, in km) along its segment from a random anchor,
// each segment cell cutting the sea within 2 to 4 pixels (plus half a pixel) of it.
void metadataactivetrenches(planet& world, GenerationScratch& scratch)
{
    const int sealevel = world.sealevel();
    const double pixel_km = terrain::legacycellkm(world);
    clearfeaturemask(scratch);
    reseedterrainpass(world, 0x4114);

    int maxsegmentid = 0;
    for (const TectonicBoundarySegment& segment : world.tectonicboundarysegments())
        maxsegmentid = max(maxsegmentid, segment.id);

    vector<int> segmentlookup(static_cast<size_t>(maxsegmentid + 1), -1);
    vector<int> segmentids;
    vector<vector<mesh::CellId>> segmentcells;

    // A trench is cut below the floor as it was before this stage (at least 450 m deep), so it
    // shows on the ocean curve's deep floor and overlapping kernels do not add up.
    vector<float> trenchfloor(world.cellcount(), 0.0f);
    for (mesh::CellId c = 0; c < cellcount(world); c++)
        trenchfloor[c] = min(world.nomfloat(c), static_cast<float>(sealevel - 450));

    for (mesh::CellId c = 0; c < cellcount(world); c++)
    {
        if (world.sea(c) == 0 || world.tectonicboundarytype(c) != BoundaryType::convergent)
            continue;

        const CrustClass crustclass = world.tectoniccrustclass(c);
        if (crustclass != CrustClass::oceanic && crustclass != CrustClass::transitional)
            continue;

        const int segmentid = world.tectonicboundarysegmentid(c);
        if (segmentid <= 0)
            continue;

        if (segmentlookup[segmentid] < 0)
        {
            segmentlookup[segmentid] = static_cast<int>(segmentcells.size());
            segmentids.push_back(segmentid);
            segmentcells.emplace_back();
        }

        segmentcells[segmentlookup[segmentid]].push_back(c);
    }

    vector<int> order(segmentids.size(), 0);
    for (size_t index = 0; index < order.size(); index++)
        order[index] = static_cast<int>(index);

    for (int index = static_cast<int>(order.size()) - 1; index > 0; index--)
        swap(order[index], order[random(0, index)]);

    const int trenchcount = min(static_cast<int>(order.size()), clamp(static_cast<int>(order.size()) / 3 + 1, 2, 6));

    for (int trenchindex = 0; trenchindex < trenchcount; trenchindex++)
    {
        const vector<mesh::CellId>& cells = segmentcells[order[trenchindex]];
        if (cells.empty())
            continue;

        // One branch per 48 legacy pixels of segment: its cells' widths (the square roots of
        // their areas) summed, in legacy pixels.
        double lengthpixels = 0.0;
        for (const mesh::CellId c : cells)
            lengthpixels += sqrt(world.grid().area_m2(c)) / 1000.0 / pixel_km;
        const int branchcount = clamp(static_cast<int>(lengthpixels / 48.0), 1, 3);

        for (int branch = 0; branch < branchcount; branch++)
        {
            const mesh::CellId anchor = cells[random(0, static_cast<int>(cells.size()) - 1)];
            const float maxsegmentdistance = static_cast<float>(random(18, 38));
            const int radius = random(2, 4);
            const float reach = static_cast<float>(radius) + 0.5f;

            for (const mesh::CellId c : cells)
            {
                // Distances in legacy pixels.
                const float segmentdistance = static_cast<float>(terrain::distancekm(world, c, anchor) / pixel_km);

                if (segmentdistance > maxsegmentdistance)
                    continue;

                const float branchstrength = 1.0f - clamp01(segmentdistance / maxsegmentdistance);

                for (const mesh::CellId n : terrain::cellswithinkm(world, c, reach * pixel_km))
                {
                    if (world.sea(n) == 0)
                        continue;

                    const float kerneldistance = static_cast<float>(terrain::distancekm(world, c, n) / pixel_km);
                    if (kerneldistance > reach)
                        continue;

                    const float kernelfalloff = 1.0f - clamp01(kerneldistance / reach);
                    const float proximity = 1.0f - clamp01(static_cast<float>(world.tectonicboundarydistance(n)) / 16.0f);
                    (void)proximity;  // read as before (the stage's declared access), unused
                    const float subsidence = clamp01(world.tectonicsubsidencetendency(n));
                    const float convergence = clamp01(static_cast<float>(world.tectonicconvergence(n)) / 255.0f);
                    const float randomness = 0.82f + static_cast<float>(random(0, 35)) / 100.0f;
                    const int trenchdepth = static_cast<int>(roundf((1200.0f + 2200.0f * branchstrength) * kernelfalloff * randomness * (0.45f + subsidence * 0.30f + convergence * 0.25f)));

                    if (trenchdepth <= 0)
                        continue;

                    const float target = max(1.0f, min(world.nomfloat(n), trenchfloor[n] - static_cast<float>(trenchdepth)));
                    if (target >= world.nomfloat(n))
                        continue;

                    world.setnomfloat(n, target);
                    world.setgeologicregime(n, GeologicRegime::trench_adjacent);
                    markfeaturecell(scratch, world, n);
                }
            }
        }
    }
}

// Raises the cells within radius legacy pixels (plus half a pixel) of a centre cell.
void applyvolcanicupliftkernel(planet& world, GenerationScratch& scratch, mesh::CellId centre, int radius, int maximumuplift, GeologicRegime regime)
{
    const double pixel_km = terrain::legacycellkm(world);
    const float reach = static_cast<float>(radius) + 0.5f;

    for (const mesh::CellId n : terrain::cellswithinkm(world, centre, reach * pixel_km))
    {
        const float distance = static_cast<float>(terrain::distancekm(world, centre, n) / pixel_km);
        if (distance > reach)
            continue;

        const float falloff = 1.0f - clamp01(distance / reach);
        const int uplift = static_cast<int>(roundf(static_cast<float>(maximumuplift) * falloff * falloff));
        if (uplift <= 0)
            continue;

        const float target = min(static_cast<float>(world.maxelevation() - 1), world.nomfloat(n) + static_cast<float>(uplift));
        if (target <= world.nomfloat(n))
            continue;

        world.setnomfloat(n, target);
        if (regime != GeologicRegime::stable)
            world.setgeologicregime(n, regime);
        markfeaturecell(scratch, world, n);
    }
}

// Arc cones on convergent boundaries and hotspots: seeds at least 24 (cones) or 34 (hotspots)
// legacy pixels apart in km, hotspots drawn by area; the counts are those of the 2048-column
// world (per legacy cell, terrain.sphere_noise.legacy_columns), so every resolution has them.
void metadatavolcanoes(planet& world, GenerationScratch& scratch)
{
    clearfeaturemask(scratch);
    rebuildshelves(world, scratch);
    reseedterrainpass(world, 0x4115);

    // The two direction fields (east and north shifts): the legacy fractal and its value half a
    // turn of longitude away, which acted as a second, independent one.
    const vector<float> eastdirection = buildterrainfractal(world, 0x4116, world.maxelevation());
    const vector<float> northdirection = terrain::spherefractalwarped(
        world, deterministiccombine(deterministiccontextseed(world.seed(), 0x4116), 7), terrain::terrainfractal(),
        static_cast<float>(world.maxelevation()), random(20, 80));
    const double pixel_km = terrain::legacycellkm(world);
    vector<mesh::CellId> selectedseeds;

    auto farenough = [&](mesh::CellId c, int minimumdistance)
    {
        for (const mesh::CellId seed : selectedseeds)
        {
            if (terrain::distancekm(world, c, seed) < static_cast<double>(minimumdistance) * pixel_km)
                return false;
        }

        return true;
    };

    struct VolcanoCandidate
    {
        mesh::CellId cell = 0;
        float score = 0.0f;
    };

    vector<VolcanoCandidate> conecandidates;

    for (mesh::CellId c = 0; c < cellcount(world); c++)
    {
        if (world.tectonicboundarytype(c) != BoundaryType::convergent)
            continue;

        const float proximity = 1.0f - clamp01(static_cast<float>(world.tectonicboundarydistance(c)) / 12.0f);
        const float uplift = clamp01(world.tectonicuplifttendency(c));
        const float convergence = clamp01(static_cast<float>(world.tectonicconvergence(c)) / 255.0f);
        const float history = clamp01(world.tectonicboundaryhistory(c));
        const float score = proximity * (0.35f + uplift * 0.35f + convergence * 0.20f + history * 0.10f);

        if (score < 0.34f)
            continue;

        conecandidates.push_back({ c, score });
    }

    for (int index = static_cast<int>(conecandidates.size()) - 1; index > 0; index--)
        swap(conecandidates[index], conecandidates[random(0, index)]);

    const double area = terrain::legacycellcount();
    const int desiredcones = clamp(static_cast<int>(area / 200000.0), 2, 10);

    for (const VolcanoCandidate& candidate : conecandidates)
    {
        if (static_cast<int>(selectedseeds.size()) >= desiredcones)
            break;

        if (!farenough(candidate.cell, 24))
            continue;

        const int peakheight = static_cast<int>(roundf(1800.0f + candidate.score * 3800.0f));
        applyvolcanicupliftkernel(world, scratch, candidate.cell, world.sea(candidate.cell) ? 2 : 1, static_cast<int>(roundf(240.0f + candidate.score * 520.0f)), GeologicRegime::convergent_arc);
        terrain::createisolatedvolcano(world, candidate.cell, world.continentalshelves(), eastdirection, northdirection, peakheight, true);
        selectedseeds.push_back(candidate.cell);
    }

    auto placenexthotspot = [&](bool oceanpreferred, int targetcount)
    {
        for (int placed = 0, attempts = 0; placed < targetcount && attempts < 512; attempts++)
        {
            const mesh::CellId c = terrain::randomcellbyarea(world);

            if (!farenough(c, 34))
                continue;

            if (world.tectonicboundarydistance(c) < 18)
                continue;

            if (oceanpreferred != (world.sea(c) == 1))
                continue;

            if (world.sea(c) == 1 && world.continentalshelves()[c])
                continue;

            const int radius = oceanpreferred ? random(4, 6) : random(3, 5);
            const int uplift = oceanpreferred ? random(260, 700) : random(180, 460);
            const int peakheight = oceanpreferred ? random(900, 2400) : random(700, 1600);
            applyvolcanicupliftkernel(world, scratch, c, radius, uplift, GeologicRegime::stable);
            terrain::createisolatedvolcano(world, c, world.continentalshelves(), eastdirection, northdirection, peakheight, false);
            selectedseeds.push_back(c);
            placed++;
        }
    };

    placenexthotspot(true, clamp(static_cast<int>(area / 280000.0) + 1, 1, 6));
    placenexthotspot(false, clamp(static_cast<int>(area / 520000.0), 0, 3));

    for (mesh::CellId c = 0; c < cellcount(world); c++)
    {
        if (world.volcanofloat(c) == 0)
            continue;

        markfeaturecell(scratch, world, c);
        if (world.strato(c))
            world.setgeologicregime(c, GeologicRegime::convergent_arc);
    }
}

vector<mesh::CellId> findcomponentdrainpath(planet& world, const vector<vector<int>>& componentids, int componentid)
{
    const mesh::CellGrid& grid = world.grid();
    const mesh::CellId count = cellcount(world);
    vector<std::uint8_t> visited(count, 0);
    vector<mesh::CellId> parent(count, mesh::nocell);
    queue<mesh::CellId> frontier;
    mesh::CellId endpoint = mesh::nocell;

    for (mesh::CellId c = 0; c < count; c++)
    {
        if (scratchcell(componentids, world, c) == componentid)
        {
            frontier.push(c);
            visited[c] = 1;
            parent[c] = c;
        }
    }

    while (frontier.empty() == false && endpoint == mesh::nocell)
    {
        const mesh::CellId current = frontier.front();
        frontier.pop();

        for (int slot = 1; slot <= mesh::slotcount; slot++)
        {
            const std::int32_t neighbour = grid.neighbour(current, slot);
            if (neighbour < 0 || visited[neighbour])
                continue;

            const mesh::CellId next = static_cast<mesh::CellId>(neighbour);
            visited[next] = 1;
            parent[next] = current;

            if (world.sea(next) == 1 && scratchcell(componentids, world, next) == 0)
            {
                endpoint = next;
                break;
            }

            frontier.push(next);
        }
    }

    vector<mesh::CellId> path;

    if (endpoint == mesh::nocell)
        return path;

    mesh::CellId c = endpoint;

    while (parent[c] != c)
    {
        path.push_back(c);
        c = parent[c];
    }

    path.push_back(c);
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
            for (mesh::CellId c = 0; c < cellcount(world); c++)
            {
                if (scratchcell(scratch.inlandSeaComponentIds, world, c) == component.id)
                    world.setnomfloat(c, max(static_cast<float>(world.sealevel() + 1), world.nomfloat(c)));
            }
        }
        else
        {
            const vector<mesh::CellId> path = findcomponentdrainpath(world, scratch.inlandSeaComponentIds, component.id);

            for (const mesh::CellId c : path)
                world.setnom(c, max(1, world.sealevel() - 10));
        }
    }

    rebuildshelves(world, scratch);
    detectinlandseas(world, scratch);
}

void applyterraintexturing(planet& world, GenerationScratch& scratch)
{
    clearfeaturemask(scratch);
    // The legacy flip and shift only made the field another, independent sample: a sphere field
    // of its own salt is one already.
    const vector<float> fractal = buildterrainfractal(world, 0x4111, world.maxelevation());
    const mesh::CellGrid& grid = world.grid();

    const int sealevel = world.sealevel();

    // The sea's texture is the fractal at zero mean and unit spread, both by area: the fractal's
    // values cover only part of its range.
    double sum = 0.0;
    double sumsquares = 0.0;
    double count = 0.0;

    for (mesh::CellId c = 0; c < cellcount(world); c++)
    {
        const double weight = grid.area_m2(c) / grid.totalarea_m2();
        const double value = static_cast<double>(fractal[c]);
        sum += weight * value;
        sumsquares += weight * value * value;
        count += weight;
    }

    const double mean = count > 0.0 ? sum / count : 0.0;
    const double spread = count > 0.0 ? sqrt(max(0.0, sumsquares / count - mean * mean)) : 0.0;
    const float seascale = spread > 0.0 ? tuning::terrain::ocean_floor::textureSdM / static_cast<float>(spread) : 0.0f;

    // The land's texture is centred at the same measured mean: centred at half the maximum
    // elevation it lowered almost all land, and the floor below held most of it at the coast.
    const float landscale = tuning::terrain::texturing::landScale;

    for (mesh::CellId c = 0; c < cellcount(world); c++)
    {
        const float oldvalue = world.nomfloat(c);
        const float noise = static_cast<float>(static_cast<double>(fractal[c]) - mean);
        float newvalue = world.sea(c)
            ? oldvalue + noise * seascale
            : oldvalue + noise * landscale;

        if (world.sea(c))
            newvalue = min(newvalue, static_cast<float>(sealevel - 10));
        else
            newvalue = max(newvalue, static_cast<float>(sealevel + 5));

        newvalue = clamp(newvalue, 1.0f, static_cast<float>(world.maxelevation() - 1));

        if (newvalue == oldvalue)
            continue;

        world.setnomfloat(c, newvalue);
        markfeaturecell(scratch, world, c);
    }
}

void rebuildterrainroughness(planet& world)
{
    const vector<float> roughness = buildterrainfractal(world, 0x4113, world.maxelevation());

    for (mesh::CellId c = 0; c < cellcount(world); c++)
        world.setroughness(c, roughness[c]);
}

void applyaridfeatures(planet& world, const GenerationExecutionContext& context)
{
    createergs(world, context.smalllake, context.largelake, context.landshape);
    createsaltpans(world, context.smalllake, context.largelake);
}

void applydeltawetlandstage(planet& world, GenerationScratch& scratch, const GenerationExecutionContext& context)
{
    // The deltas come from the erosion (erosion plan C5b): marine deposition builds land at
    // river mouths and hydrology.rivers writes the distributaries over it. With deltas turned off, the
    // distributaries hydrology.rivers wrote are cleared here, before the wetlands read them
    // (hydrology.rivers does not see the option); a world without erosion has none.
    const planet& view = world;
    const bool eroded = view.erosionduration() > 0.0f;
    if (!context.dodeltas && eroded)
    {
        for (mesh::CellId c = 0; c < cellcount(world); c++)
        {
            world.setdeltadir(c, 0);
            world.setdeltajan(c, 0);
            world.setdeltajul(c, 0);
        }
    }

    createwetlands(world, context.smalllake);
    removeexcesswetlands(world);
    refineroughnessmap(world);
    // The lakes are hydrology.lake_levels' (erosion plan H4); nothing here spreads or joins them.
    // The climate classes and biomes through the climate adapter: the world itself on lat-lon,
    // its shadow on a cube.
    runclimateprocess(GenerationStageId::deltas_wetlands_and_roughness, world, scratch.climateShadow, [](planet& target)
    {
        checkpoleclimates(target.climate());
        target.syncseasonalclimatefromlegacy();
        createclimatemap(target.climate());
        createbiomemap(target.climate());
        writeclimatediagnostics(target.climate());
    });
}

}

void detectinlandseas(planet& world, GenerationScratch& scratch)
{
    const mesh::CellGrid& grid = world.grid();
    const mesh::CellId count = cellcount(world);
    scratch.inlandSeaComponentIds = makeintgrid(0);
    scratch.inlandSeaComponents.clear();

    vector<std::uint8_t> visited(count, 0);
    vector<vector<mesh::CellId>> components(1);
    int largestcomponentid = 0;
    int largestcomponentarea = 0;

    for (mesh::CellId c = 0; c < count; c++)
    {
        if (visited[c] || world.sea(c) == 0)
            continue;

        queue<mesh::CellId> frontier;
        vector<mesh::CellId> cells;
        frontier.push(c);
        visited[c] = 1;

        while (frontier.empty() == false)
        {
            const mesh::CellId current = frontier.front();
            frontier.pop();
            cells.push_back(current);

            for (int slot = 1; slot <= mesh::slotcount; slot++)
            {
                const std::int32_t neighbour = grid.neighbour(current, slot);
                if (neighbour < 0 || visited[neighbour] || world.sea(static_cast<mesh::CellId>(neighbour)) == 0)
                    continue;

                visited[neighbour] = 1;
                frontier.push(static_cast<mesh::CellId>(neighbour));
            }
        }

        components.push_back(std::move(cells));

        if (static_cast<int>(components.back().size()) > largestcomponentarea)
        {
            largestcomponentarea = static_cast<int>(components.back().size());
            largestcomponentid = static_cast<int>(components.size() - 1);
        }
    }

    for (size_t componentindex = 1; componentindex < components.size(); componentindex++)
    {
        if (static_cast<int>(componentindex) == largestcomponentid)
            continue;

        GenerationWaterComponent component;
        component.id = static_cast<int>(scratch.inlandSeaComponents.size() + 1);
        legacyslot(grid, components[componentindex].front(), component.seedX, component.seedY);
        component.cellCount = static_cast<int>(components[componentindex].size());
        scratch.inlandSeaComponents.push_back(component);

        for (const mesh::CellId cell : components[componentindex])
            scratchcell(scratch.inlandSeaComponentIds, world, cell) = component.id;
    }
}

bool executegenerationstage(GenerationStageId stageId, planet& world, GenerationScratch& scratch, WorldStageParameters& parameters, const GenerationExecutionContext& context, string* errormessage)
{
    const GenerationStageDefinition& stagedefinition = getgenerationstage(stageId);

    // A cube world runs only the processes ported to the cell grid (sphere-mesh.md section 4).
    if (world.grid().kind() == mesh::GridKind::cube)
    {
        const WorldProcessDefinition* process = nullptr;
        for (const WorldProcessDefinition& candidate : getworldprocesses())
            if (candidate.stage == stageId)
                process = &candidate;
        if (process == nullptr || !process->cubeready)
        {
            if (errormessage != nullptr)
                *errormessage = "Process " + string(process != nullptr ? process->id : stagedefinition.label) +
                    " is not ready for a cube world (--grid cube); it runs on --grid latlon.";
            return false;
        }
    }

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
        scratch.climateShadow.reset();   // a new world: climate's shadow starts again
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
        {
            world.clear();
            scratch.climateShadow.reset();
        }
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
        parameters.seaLevel = world.sealevel();
        return true;

    case GenerationStageId::inland_seas:
        applyinlandseapolicies(world, scratch);
        return true;

    case GenerationStageId::archipelagos:
    {
        const vector<std::uint8_t> removedland(world.cellcount(), 0);
        terrain::makearchipelagos(world, removedland, context.landshape);
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

    case GenerationStageId::sea_level:
        applysealevel(world, scratch);
        parameters.seaLevel = world.sealevel();
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
        // Through the climate adapter (climate_adapter.hpp): the world itself on lat-lon, its shadow on a cube.
        runclimateprocess(stageId, world, scratch.climateShadow,
            [&](planet& target) { executeclimatestage(stageId, target.climate(), climatecontext); });
        return true;

    case GenerationStageId::sea_ice_and_tides:
        runclimateprocess(stageId, world, scratch.climateShadow, [&](planet& target)
        {
            executeclimatestage(stageId, target.climate(), climatecontext);
            if (target.seatotal() > 0)
                createtidalmap(target);
        });
        return true;

    case GenerationStageId::rainfall:
        // The planet keeps its own land/sea counts for the stages after climate.
        getlandandseatotals(world);
        runclimateprocess(stageId, world, scratch.climateShadow,
            [&](planet& target) { executeclimatestage(stageId, target.climate(), climatecontext); });
        return true;

    case GenerationStageId::flow_routing:
        applyflowroutingstage(world);
        return true;

    case GenerationStageId::fluvial_erosion:
    {
        // The land and sea totals are not recounted here: they are undeclared writes of
        // this process, and the stages after it read the counts of the rainfall stage.
        FluvialErosionCheckpoints checkpoints;
        if (context.evolverCheckpoint)
        {
            checkpoints.everyMyr = parameters.erosionCheckpointEveryMyr;
            checkpoints.atMyr = parameters.erosionCheckpointsMyr;
            checkpoints.checkpoint = [&](double timemyr) { context.evolverCheckpoint(stageId, timemyr); };
        }
        applyfluvialerosionstage(world, checkpoints);
        // The erosion moves coasts: refresh the planet's land and sea counts here, so the
        // climate after it (procedural-world's @full climate) does not have to.
        getlandandseatotals(world);
        return true;
    }

    case GenerationStageId::rivers:
        applyriverstage(world);
        return true;

    case GenerationStageId::lake_levels:
        applylakelevelstage(world);
        return true;

    case GenerationStageId::mountain_temperature_lapse:
    case GenerationStageId::climates_and_biomes:
        runclimateprocess(stageId, world, scratch.climateShadow,
            [&](planet& target) { executeclimatestage(stageId, target.climate(), climatecontext); });
        return true;

    case GenerationStageId::arid_features:
        applyaridfeatures(world, context);
        return true;

    case GenerationStageId::deltas_wetlands_and_roughness:
        applydeltawetlandstage(world, scratch, context);
        return true;

    case GenerationStageId::subgrid_statistics:
        applysubgridstatisticsstage(world);
        return true;

    case GenerationStageId::finalize_layers:
        finalizegeneratedworld(world, scratch, context);
        return true;
    }

    return false;
}


bool stageusesmanualpreview(GenerationStageId stageid)
{
    return stageid == GenerationStageId::inland_seas;
}

void rebuildheightdistribution(const planet& world, GenerationScratch& scratch)
{
    const int minimum = 1;
    const int maximum = world.maxelevation() - 1;

    scratch.heightHistogram.assign(kHeightHistogramBins, 0);
    scratch.heightCdf.assign(kHeightHistogramBins, 0.0f);

    int total = 0;

    for (mesh::CellId c = 0; c < cellcount(world); c++)
    {
        const float normalized = normalizeterrainvalue(world.nomfloat(c), minimum, maximum);
        int bin = static_cast<int>(floorf(normalized * static_cast<float>(kHeightHistogramBins)));
        bin = clamp(bin, 0, kHeightHistogramBins - 1);
        scratch.heightHistogram[bin]++;
        total++;
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

    scratch.inlandSeaComponentIds = makeintgrid(0);
    scratch.heightHistogram.assign(kHeightHistogramBins, 0);
    scratch.heightCdf.assign(kHeightHistogramBins, 0.0f);
    scratch.inlandSeaComponents.clear();
    clearfeaturemask(scratch);
    scratch.climateShadow.reset();
}

void finalizegeneratedworld(planet& world, GenerationScratch& scratch, const GenerationExecutionContext& context)
{
    generatephysicalworldlayers(world, world.continentalshelves());
    generatesocialworld(world, context.socialoptions);
}
