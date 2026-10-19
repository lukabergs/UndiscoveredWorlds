// Tests of the regional pass's cube forms (T8 W2-W2, sphere-mesh.md section 8), suite Landscape:
// Landscape.RegionalCube*. Synthetic cube worlds (level 0 N 32) and cube windows built from rings of
// level-0 cells; two windows that share cells far enough from both domains' edges must give them
// bit-identical values (the overlap property of regional_pass.hpp on the cube).

#include "landscape_fixture.hpp"
#include "regional_test_world.hpp"

#include "core/mesh/cell_grid.hpp"
#include "core/mesh/cube_geometry.hpp"
#include "core/mesh/window_grid.hpp"
#include "simulations/landscape/channel_network.hpp"
#include "simulations/landscape/level_blocks.hpp"
#include "simulations/landscape/regional_coast.hpp"
#include "simulations/landscape/regional_edges.hpp"
#include "simulations/landscape/regional_surface.hpp"
#include "simulations/landscape/regional_cirques.hpp"
#include "simulations/landscape/regional_glacial.hpp"
#include "simulations/landscape/regional_ice.hpp"
#include "simulations/landscape/regional_troughs.hpp"
#include "simulations/landscape/regional_valleys.hpp"
#include "simulations/landscape/regional_water.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <memory>
#include <vector>

namespace
{
constexpr int cuben = 32;

// The level-0 keys within `rings` rings of a cell, ascending.
std::vector<mesh::LevelKey> disc(mesh::LevelKey centre, int rings) { return mesh::dilatekeys(cuben, { centre }, rings); }

// The cube window of the children at k of a set of level-0 keys.
std::shared_ptr<const mesh::CellGrid> windowof(const std::vector<mesh::LevelKey>& parents, int k)
{
    std::vector<mesh::LevelKey> keys;
    for (const mesh::LevelKey parent : parents)
        for (const mesh::LevelKey child : mesh::cubechildkeys(cuben, k, parent))
            keys.push_back(child);
    return mesh::cubewindowgrid(cuben * k, radius_m, keys);
}

// A cube world: smooth relief of a few thousand metres with coasts, seas and one-cell islands.
RegionalGlobalInputs cubeworld()
{
    RegionalGlobalInputs inputs;
    inputs.grid = mesh::makecubegrid(cuben, radius_m);
    inputs.radius_m = radius_m;
    inputs.sealevel = 0.0f;
    inputs.seed = 1;
    inputs.lapserate_c_per_km = 6.5f;
    const std::size_t size = inputs.grid->cellcount();
    inputs.elevation_m.resize(size);
    inputs.network.sea.assign(size, 0);
    inputs.lakesurface_m.assign(size, 0.0f);
    inputs.island.assign(size, 0);
    for (mesh::CellId c = 0; c < size; c++)
    {
        const mesh::Vec3 p = inputs.grid->centre(c);
        const double z = 1500.0 * std::sin(5.0 * p.x + 1.0) * std::cos(4.0 * p.y - 0.5) + 900.0 * std::sin(7.0 * p.z) - 200.0;
        inputs.elevation_m[c] = static_cast<float>(z);
        inputs.network.sea[c] = z < 0.0 ? 1 : 0;
    }
    return inputs;
}

// The coast and the land-share rule at the defaults of configs/parameters/landscape.json, with fewer
// share and sign-keeping passes (the box's margin shrinks with them).
RegionalParameters coastparameters()
{
    RegionalParameters parameters;
    parameters.meanpasses = 2;
    parameters.warp.amplitude_cells = 0.5;
    parameters.warp.wavelength_cells = 2.0;
    parameters.warp.fade_m = 300.0;
    parameters.coast.noiseamplitude_m = 50.0;
    parameters.coast.noisepersistence = 0.6;
    parameters.coast.noiserelieffraction = 0.25;
    parameters.coast.seadrop_m = 300.0;
    parameters.coast.landshare = 0.5;
    parameters.coast.shorefloor_m = 20.0;
    parameters.coast.sharepasses = 2;
    parameters.coast.signedpasses = 2;
    parameters.coast.onesided = 1;
    parameters.coast.islandmargin = 0.125;
    parameters.coast.warpcells = 0.15;
    parameters.coast.warppersistence = 0.6;
    parameters.coast.pocketreach = 2;
    parameters.edges.smoothgate = 3;
    parameters.edges.straightmax = 4;
    parameters.edges.chunkmin = 2;
    parameters.edges.chunkmax = 7;
    return parameters;
}

std::uint32_t floatbits(float value)
{
    std::uint32_t bits = 0;
    std::memcpy(&bits, &value, sizeof bits);
    return bits;
}
}

TEST_F(Landscape, RegionalCubeStartSurfaceOverlaps)
{
    // Two windows at k 4, discs of 11 rings around a cell at a cube corner and one 5 cells along its
    // face edge: on the cells of the level-0 cells within 3 rings of both centres (8 rings inside both
    // domains) the land-share start surface, its base and the sea are bit-identical, the land share
    // holds on the domain, and no cell of the start surface is far from its block's mean.
    constexpr int k = 4;
    const RegionalGlobalInputs inputs = cubeworld();
    const RegionalParameters parameters = coastparameters();
    const mesh::LevelKey corner = mesh::cubekey(cuben, 1, cuben - 1, 0), along = mesh::cubekey(cuben, 1, cuben - 1, 5);
    const auto a = windowof(disc(corner, 11), k);
    const auto b = windowof(disc(along, 11), k);
    const landscape::FineWarp warpa = landscape::regionalwarp(*a, cuben, inputs.seed, parameters);
    const landscape::FineWarp warpb = landscape::regionalwarp(*b, cuben, inputs.seed, parameters);
    landscape::StartSurface sa = landscape::regionalstartsurface(inputs, *a, warpa, parameters, false);
    landscape::StartSurface sb = landscape::regionalstartsurface(inputs, *b, warpb, parameters, false);
    ASSERT_EQ(sa.start_m.size(), a->cellcount());
    std::vector<float> surfacea = sa.start_m, surfaceb = sb.start_m;
    const std::vector<uint8_t> seaa = landscape::regionalsea(inputs, *a, warpa, parameters, surfacea, nullptr, false);
    const std::vector<uint8_t> seab = landscape::regionalsea(inputs, *b, warpb, parameters, surfaceb, nullptr, false);
    const std::vector<mesh::LevelKey> innera = disc(corner, 3), innerb = disc(along, 3);
    int shared = 0, starts = 0, bases = 0, seas = 0;
    for (mesh::CellId c = 0; c < a->cellcount(); c++)
    {
        const mesh::LevelKey key = a->globalkey(c);
        const mesh::LevelKey parent = mesh::cubeparentkey(cuben * k, k, key);
        if (!std::binary_search(innera.begin(), innera.end(), parent) || !std::binary_search(innerb.begin(), innerb.end(), parent))
            continue;
        const mesh::CellId other = b->localofkey(key);
        ASSERT_NE(other, mesh::nocell);
        shared++;
        starts += floatbits(sa.start_m[c]) != floatbits(sb.start_m[other]) ? 1 : 0;
        bases += floatbits(sa.base_m[c]) != floatbits(sb.base_m[other]) ? 1 : 0;
        seas += seaa[c] != seab[other] ? 1 : 0;
    }
    int sea = 0;
    for (const uint8_t value : seaa)
        sea += value;
    std::printf("[cube coast] %d shared cells: %d starts, %d bases, %d sea flags differ; sea %d of %zu cells; share passes %d, "
                "misses %d; block-mean error %.3g m; pockets %d filled, %d left\n",
                shared, starts, bases, seas, sea, seaa.size(), sa.sharepasses, sa.sharemisses, sa.blockmeanerror_m,
                sa.pocketsfilled, sa.pockets);
    EXPECT_GT(shared, 0);
    EXPECT_EQ(starts, 0);
    EXPECT_EQ(bases, 0);
    EXPECT_EQ(seas, 0);
    EXPECT_GT(sea, 0);
    EXPECT_LT(sea, static_cast<int>(seaa.size()));
    EXPECT_LT(sa.blockmeanerror_m, 50.0);
}

namespace
{
// The cube world with lakes (a lake surface 40 m above the cells between 100 and 400 m), ice above
// within 0.1 rad of a cell on face 0 (where the overlap test looks), full glaciation within 0.14 rad
// and an ELA of 600 m.
RegionalGlobalInputs lakeworld()
{
    RegionalGlobalInputs inputs = cubeworld();
    const std::size_t size = inputs.grid->cellcount();
    inputs.icethickness_m.assign(size, 0.0f);
    inputs.glaciation.assign(size, 0.0f);
    inputs.equilibriumline_m.assign(size, 600.0f);
    inputs.network.filledelevation_m = inputs.elevation_m;
    inputs.network.direction.assign(size, 0);
    const mesh::Vec3 focus = inputs.grid->centre(static_cast<mesh::CellId>(mesh::cubekey(cuben, 0, cuben - 2, 21)));
    for (mesh::CellId c = 0; c < size; c++)
    {
        const mesh::Vec3 p = inputs.grid->centre(c);
        const float z = inputs.elevation_m[c];
        if (z > 100.0f && z < 400.0f)
            inputs.lakesurface_m[c] = z + 40.0f;
        if (mesh::dot(p, focus) > 0.995)
            inputs.icethickness_m[c] = 200.0f;
        if (mesh::dot(p, focus) > 0.99)
            inputs.glaciation[c] = 1.0f;
    }
    return inputs;
}

// Synthetic fine fields of a cube window, pure functions of the cell's key: the start elevation G of
// the base elevation, a channel on every 29th key (discharge 2,000 m3/s, bed 15 m below), the sea below 0 m.
struct CubeFine
{
    std::vector<float> elevation;
    std::vector<uint8_t> sea;
    landscape::FineChannels channels;
};

CubeFine cubefine(const RegionalGlobalInputs& inputs, const mesh::CellGrid& window, int k)
{
    CubeFine fine;
    fine.elevation = landscape::upsample(inputs.elevation_m, *inputs.grid, window);
    const std::size_t count = window.cellcount();
    fine.sea.assign(count, 0);
    fine.channels.channel.assign(count, 0);
    fine.channels.direction.assign(count, 0);
    fine.channels.bed_m.assign(count, 0.0f);
    fine.channels.discharge_m3s.assign(count, 0.0f);
    fine.channels.globalcell.assign(count, -1);
    for (mesh::CellId c = 0; c < count; c++)
    {
        const mesh::LevelKey key = window.globalkey(c);
        fine.sea[c] = fine.elevation[c] < 0.0f ? 1 : 0;
        if (fine.sea[c] == 0 && key % 29 == 0)
        {
            fine.channels.channel[c] = 1;
            fine.channels.discharge_m3s[c] = 2000.0f;
            fine.channels.bed_m[c] = fine.elevation[c] - 15.0f;
            fine.channels.globalcell[c] =
                static_cast<std::int32_t>(inputs.grid->localofkey(mesh::cubeparentkey(window.facecells(), k, key)));
        }
    }
    return fine;
}
}

TEST_F(Landscape, RegionalCubeWaterAndFormsOverlap)
{
    // The cube forms of the water (shore water with R4's indicators, the flood of the loop), the glacial
    // fields, the valley floors, the troughs, the cirques and the ice above the ELA on two windows of 7
    // rings around cells 4 apart across a face edge: on the cells of the level-0 cells within 3 rings of
    // both centres every field is bit-identical.
    constexpr int k = 4;
    const RegionalGlobalInputs inputs = lakeworld();
    RegionalParameters parameters = coastparameters();
    parameters.shores.shorenoise = 0.3;
    parameters.shores.terrain_m = 50.0;
    parameters.valleys.floor_km = 3.0;
    parameters.valleys.reference_m3s = 1000.0;
    parameters.valleys.floormax_km = 30.0;
    parameters.valleys.sideslope = 0.2;
    parameters.valleys.reachcells = 3;
    parameters.troughs.width_km = 20.0;
    parameters.troughs.reference_km2 = 1000.0;
    parameters.troughs.widthexponent = 0.5;
    parameters.troughs.widthmax_km = 40.0;
    parameters.troughs.wall_m = 400.0;
    parameters.troughs.exponent = 2.0;
    parameters.troughs.glaciationmin = 0.5;
    parameters.cirques.radius_km = 15.0;
    parameters.cirques.candidates = 2;
    parameters.cirques.band_m = 3000.0;
    parameters.ice.thicknessfraction = 0.5;
    const mesh::LevelKey left = mesh::cubekey(cuben, 0, cuben - 3, 20), right = mesh::cubekey(cuben, 0, cuben - 1, 22);
    struct Result
    {
        std::shared_ptr<const mesh::CellGrid> window;
        std::vector<float> elevation, lake, ice;
    };
    const auto run = [&](mesh::LevelKey centre) {
        Result result;
        result.window = windowof(disc(centre, 7), k);
        const mesh::CellGrid& window = *result.window;
        CubeFine fine = cubefine(inputs, window, k);
        const landscape::FineWarp warp = landscape::regionalwarp(window, cuben, inputs.seed, parameters);
        const landscape::FineGlacial glacial = landscape::regionalglacialfields(inputs, window, parameters);
        landscape::carvevalleys(inputs, window, parameters, fine.channels, fine.sea, fine.elevation, 0.0, 0.0);
        landscape::carvetroughs(inputs, window, parameters, glacial, fine.channels, fine.sea, fine.elevation, 0.0, 0.0);
        landscape::carvecirques(inputs, window, parameters, glacial, fine.channels, fine.sea, fine.elevation, 0.0, 0.0);
        landscape::LakeInheritance lakes;
        lakes.flood = true;
        landscape::FineWater water =
            landscape::regionalwater(inputs, window, warp, parameters, fine.sea, fine.channels, fine.elevation, false, lakes);
        landscape::iceabovetheela(inputs, window, parameters, glacial, fine.sea, fine.elevation, water.icethickness_m);
        result.elevation = fine.elevation;
        result.lake = water.lakesurface_m;
        result.ice = water.icethickness_m;
        return result;
    };
    const Result a = run(left), b = run(right);
    const std::vector<mesh::LevelKey> innera = disc(left, 3), innerb = disc(right, 3);
    int shared = 0, elevations = 0, lakes = 0, ices = 0, wet = 0, iced = 0;
    for (mesh::CellId c = 0; c < a.window->cellcount(); c++)
    {
        const mesh::LevelKey key = a.window->globalkey(c);
        const mesh::LevelKey parent = mesh::cubeparentkey(cuben * k, k, key);
        if (!std::binary_search(innera.begin(), innera.end(), parent) || !std::binary_search(innerb.begin(), innerb.end(), parent))
            continue;
        const mesh::CellId other = b.window->localofkey(key);
        ASSERT_NE(other, mesh::nocell);
        shared++;
        elevations += floatbits(a.elevation[c]) != floatbits(b.elevation[other]) ? 1 : 0;
        lakes += floatbits(a.lake[c]) != floatbits(b.lake[other]) ? 1 : 0;
        ices += floatbits(a.ice[c]) != floatbits(b.ice[other]) ? 1 : 0;
        wet += a.lake[c] > 0.0f ? 1 : 0;
        iced += a.ice[c] > 0.0f ? 1 : 0;
    }
    std::printf("[cube water] %d shared cells (%d wet, %d iced): %d elevations, %d lake surfaces, %d ice thicknesses differ\n",
                shared, wet, iced, elevations, lakes, ices);
    EXPECT_GT(shared, 0);
    EXPECT_EQ(elevations, 0);
    EXPECT_EQ(lakes, 0);
    EXPECT_EQ(ices, 0);
}

namespace
{
// The cube world's drainage network: steepest descent on the world's table (slot codes) to a strictly
// lower neighbour, the sea and sinks without a direction, the areas accumulated downstream.
landscape::GlobalNetwork cubenetwork(const RegionalGlobalInputs& inputs)
{
    const mesh::CellGrid& grid = *inputs.grid;
    const std::size_t size = grid.cellcount();
    landscape::GlobalNetwork network;
    network.direction.assign(size, 0);
    network.sea = inputs.network.sea;
    network.filledelevation_m = inputs.elevation_m;
    network.drainagearea_km2.assign(size, 0.0f);
    network.discharge_m3s.assign(size, 0.0f);
    network.januaryflow.assign(size, 0);
    network.julyflow.assign(size, 0);
    std::vector<std::size_t> receiver(size, size);
    for (mesh::CellId c = 0; c < size; c++)
    {
        if (network.sea[c] != 0)
            continue;
        double steepest = 0.0;
        for (int slot = 1; slot <= mesh::slotcount; slot++)
        {
            const std::int32_t other = grid.neighbour(c, slot);
            if (other < 0)
                continue;
            const double drop = (static_cast<double>(inputs.elevation_m[c]) - inputs.elevation_m[static_cast<std::size_t>(other)]) /
                                grid.distance_m(c, slot);
            if (drop > steepest)
            {
                steepest = drop;
                network.direction[c] = static_cast<uint8_t>(slot);
                receiver[c] = static_cast<std::size_t>(other);
            }
        }
    }
    // Areas in order of decreasing elevation (receivers are strictly lower).
    std::vector<std::size_t> order(size);
    for (std::size_t c = 0; c < size; c++)
        order[c] = c;
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        return inputs.elevation_m[a] > inputs.elevation_m[b] || (inputs.elevation_m[a] == inputs.elevation_m[b] && a < b);
    });
    std::vector<double> area(size, 0.0);
    for (const std::size_t c : order)
    {
        area[c] += grid.area_m2(static_cast<mesh::CellId>(c)) / 1.0e6;
        if (receiver[c] < size)
            area[receiver[c]] += area[c];
    }
    for (std::size_t c = 0; c < size; c++)
    {
        network.drainagearea_km2[c] = static_cast<float>(area[c]);
        network.discharge_m3s[c] = static_cast<float>(0.01 * area[c]);
        network.januaryflow[c] = static_cast<int32_t>(0.01 * area[c]);
        network.julyflow[c] = static_cast<int32_t>(0.005 * area[c]);
    }
    return network;
}
}

TEST_F(Landscape, RegionalCubeChannelsOverlap)
{
    // rasterisechannels on two cube windows of 7 rings around cells 4 apart across a face edge (k 8,
    // the R3 geometry with roughness, transits and springs, ending at the fine coast, beds on the start
    // surface): on the cells of the level-0 cells within 3 rings of both centres every channel field is
    // bit-identical, and every painted cell's direction steps to a neighbour on the fine table.
    // Two worlds: the coast world (short links into many seas) and a funnel (one sea cell at the
    // centre of face 4, the land rising with the angle from it, so long rivers share their outlet and
    // diagonal links take transit blocks).
    constexpr int k = 8;
    for (const bool funnel : { false, true })
    {
    RegionalGlobalInputs inputs = cubeworld();
    if (funnel)
    {
        // The outlet about 50 degrees north-east of the windows, so the links there run diagonally.
        const mesh::CellId leftcell = static_cast<mesh::CellId>(mesh::cubekey(cuben, 0, cuben - 3, 20));
        const mesh::Vec3 towards = mesh::normalise(inputs.grid->centre(leftcell) +
                                                   0.8 * (inputs.grid->east(leftcell) + inputs.grid->north(leftcell)));
        const mesh::CellId outletcell = inputs.grid->locate(towards);
        const mesh::Vec3 focus = inputs.grid->centre(outletcell);
        for (mesh::CellId c = 0; c < inputs.grid->cellcount(); c++)
        {
            const mesh::Vec3 p = inputs.grid->centre(c);
            inputs.elevation_m[c] = static_cast<float>(3000.0 * std::acos(std::clamp(mesh::dot(p, focus), -1.0, 1.0)) +
                                                       40.0 * std::sin(17.0 * p.x + 11.0 * p.y));
            inputs.network.sea[c] = c == outletcell ? 1 : 0;
        }
        inputs.elevation_m[outletcell] = -10.0f;
    }
    inputs.network = cubenetwork(inputs);
    landscape::ChannelParameters parameters;
    parameters.pinarea_km2 = 20000.0;
    parameters.seed = 77;
    parameters.meanderamplitude = 0.2;
    parameters.junctionjitter = 0.25;
    parameters.edgemargin = 0.1875;
    parameters.springspercell = 2;
    parameters.diagonaltransit = 1;
    parameters.meanderroughness = 0.08;
    const mesh::LevelKey left = mesh::cubekey(cuben, 0, cuben - 3, 20), right = mesh::cubekey(cuben, 0, cuben - 1, 22);
    struct Result
    {
        std::shared_ptr<const mesh::CellGrid> window;
        landscape::FineChannels channels;
    };
    const auto run = [&](mesh::LevelKey centre) {
        Result result;
        result.window = windowof(disc(centre, 7), k);
        const std::vector<float> start = landscape::upsample(inputs.elevation_m, *inputs.grid, *result.window);
        std::vector<uint8_t> sea(start.size());
        for (std::size_t cell = 0; cell < sea.size(); cell++)
            sea[cell] = start[cell] < 0.0f ? 1 : 0;
        result.channels = landscape::rasterisechannels(*inputs.grid, *result.window, inputs.network, parameters, 0.0f, start, sea);
        return result;
    };
    const Result a = run(left), b = run(right);
    const std::vector<mesh::LevelKey> innera = disc(left, 3), innerb = disc(right, 3);
    int shared = 0, channel = 0, differ = 0, badsteps = 0;
    for (mesh::CellId c = 0; c < a.window->cellcount(); c++)
    {
        const landscape::FineChannels& ca = a.channels;
        if (ca.channel[c] != 0 && ca.direction[c] != 0 && a.window->neighbour(c, ca.direction[c]) >= 0 &&
            ca.channel[static_cast<std::size_t>(a.window->neighbour(c, ca.direction[c]))] == 0)
            badsteps++;
        const mesh::LevelKey key = a.window->globalkey(c);
        const mesh::LevelKey parent = mesh::cubeparentkey(cuben * k, k, key);
        if (!std::binary_search(innera.begin(), innera.end(), parent) || !std::binary_search(innerb.begin(), innerb.end(), parent))
            continue;
        const mesh::CellId o = b.window->localofkey(key);
        ASSERT_NE(o, mesh::nocell);
        const landscape::FineChannels& cb = b.channels;
        shared++;
        channel += ca.channel[c];
        differ += ca.channel[c] != cb.channel[o] || ca.direction[c] != cb.direction[o] ||
                          floatbits(ca.bed_m[c]) != floatbits(cb.bed_m[o]) || ca.discharge_m3s[c] != cb.discharge_m3s[o] ||
                          ca.globalcell[c] != cb.globalcell[o]
                      ? 1 : 0;
    }
    std::printf("[cube channels] %d shared cells, %d channel cells, %d differ; diagonal links %d (%d through the corner), "
                "merged %d, dry tails %d, steps onto unpainted cells %d\n",
                shared, channel, differ, a.channels.diagonallinks, a.channels.cornerlinks, a.channels.mergedpieces,
                a.channels.drytails, badsteps);
    EXPECT_GT(shared, 0);
    EXPECT_GT(channel, 0);
    EXPECT_EQ(differ, 0);
    EXPECT_EQ(badsteps, 0);
    if (funnel)
        EXPECT_LT(a.channels.cornerlinks, a.channels.diagonallinks);
    }
}

TEST_F(Landscape, RegionalCubeRefinedCoastKeepsTheParentShares)
{
    // refinedcoast at level 2 of a cube window (parents the level-1 cells of a disc, factor 4): the sea
    // it returns lies inside the open cells, and every parent sea cell's children hold at most half land
    // where the displacement moved them (the guard), else the undisplaced sea stays.
    constexpr int k1 = 4, factor = 4;
    const RegionalGlobalInputs inputs = cubeworld();
    const RegionalParameters parameters = coastparameters();
    const mesh::LevelKey centre = mesh::cubekey(cuben, 4, 3, 3);
    const auto parent = windowof(disc(centre, 6), k1);
    std::vector<float> parentelevation = landscape::upsample(inputs.elevation_m, *inputs.grid, *parent);
    std::vector<uint8_t> parentsea(parent->cellcount());
    for (std::size_t cell = 0; cell < parentsea.size(); cell++)
        parentsea[cell] = parentelevation[cell] <= 0.0f ? 1 : 0;
    // The level: the children of the parent cells within 3 level-0 rings of the centre.
    std::vector<mesh::LevelKey> levelparents;
    for (const mesh::LevelKey zero : disc(centre, 3))
        for (const mesh::LevelKey one : mesh::cubechildkeys(cuben, k1, zero))
            levelparents.push_back(one);
    std::sort(levelparents.begin(), levelparents.end());
    std::vector<mesh::LevelKey> keys;
    for (const mesh::LevelKey one : levelparents)
        for (const mesh::LevelKey two : mesh::cubechildkeys(cuben * k1, factor, one))
            keys.push_back(two);
    const auto level = mesh::cubewindowgrid(cuben * k1 * factor, radius_m, keys);
    std::vector<float> surface = landscape::upsample(inputs.elevation_m, *inputs.grid, *level);
    const landscape::FineWarp warp = landscape::regionalwarp(*level, cuben, inputs.seed, parameters);
    landscape::PocketCount filled, left;
    const std::vector<uint8_t> sea = landscape::refinedcoast(inputs, *level, factor, *parent, parentsea, parentelevation, warp,
                                                             parameters, nullptr, surface, nullptr, &filled, &left);
    std::vector<float> below(level->cellcount(), -1.0f);
    const std::vector<uint8_t> open = landscape::regionalsea(inputs, *level, warp, parameters, below, nullptr, true);
    int seacells = 0;
    for (std::size_t cell = 0; cell < sea.size(); cell++)
    {
        seacells += sea[cell];
        if (sea[cell] != 0)
        {
            EXPECT_NE(open[cell], 0) << "cell " << cell;
            EXPECT_LE(surface[cell], 0.0f) << "cell " << cell;
        }
    }
    std::printf("[cube coast] level 2: %d sea of %zu cells, pockets %d filled, %d left\n", seacells, sea.size(),
                filled.pockets, left.pockets);
    EXPECT_GT(seacells, 0);
}
