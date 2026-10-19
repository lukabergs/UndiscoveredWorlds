// Hydrology on a cube world (docs/architecture/sphere-mesh.md sections 3, 5 and 12; T8 W2-H):
// on a cube planet of N 16 with a synthetic sloping surface the routing, the rivers and the lake
// levels run; a river crosses a face edge; the drainage basins tile the land; and the river
// directions show no preferred face-axis direction by latitude band beyond D8's own bias.

#include <cmath>
#include <cstdint>
#include <iostream>
#include <memory>
#include <vector>

#include <gtest/gtest.h>

#include "core/climate_seasons.hpp"
#include "core/mesh/cell_grid.hpp"
#include "simulations/hydrology/drainage_basins.hpp"
#include "simulations/hydrology/flow_queries.hpp"
#include "simulations/hydrology/flow_routing.hpp"
#include "simulations/hydrology/lake_levels.hpp"
#include "simulations/hydrology/rivers.hpp"
#include "simulations/hydrology/tides.hpp"
#include "world/planet.hpp"

namespace
{
constexpr int facecells = 16;
constexpr double pi = 3.14159265358979323846;
constexpr float sealevel = 31000.0f;

double latitudedeg(const mesh::CellGrid& grid, mesh::CellId c) { return grid.latitude_rad(c) * 180.0 / pi; }

// A cube planet sloping down from the north pole to a southern ocean (below 30 S): the
// elevation is a function of latitude only, so the true downhill is due south everywhere.
// With `bowl`, a closed bowl 3 km deep (a Gaussian of 1500 km) is cut into the land at 30 N, 20 E
// and a ripple added, so the routing has a depression to fill and the lake balance a lake.
std::unique_ptr<planet> slopingcubeworld(bool bowl)
{
    auto world = std::make_unique<planet>();
    world->setgrid(mesh::GridKind::cube, 4 * facecells);
    world->setseed(1);
    world->setmaxelevation(48000);
    world->setsealevel(static_cast<int>(sealevel));
    world->setriverfactor(15.0f);
    world->seterosionduration(1.0f);   // the eroded rule: the sea is the routing's ocean
    const mesh::CellGrid& grid = world->grid();
    const mesh::Vec3 bowlcentre = mesh::unitvector(30.0 * pi / 180.0, 20.0 * pi / 180.0);
    for (mesh::CellId c = 0; c < world->cellcount(); c++)
    {
        const double latitude = grid.latitude_rad(c);
        double height = sealevel + 4000.0 * (std::sin(latitude) + 0.5);   // 0 m at 30 S, 6 km at the pole
        if (bowl)
        {
            const mesh::Vec3 p = grid.centre(c);
            const double distance_km = std::acos(std::min(1.0, mesh::dot(p, bowlcentre))) * grid.radius_m() / 1000.0;
            height -= 3000.0 * std::exp(-0.5 * (distance_km / 1500.0) * (distance_km / 1500.0));
            height += 40.0 * std::sin(7.0 * p.x) * std::sin(5.0 * p.y);
        }
        world->setnomfloat(c, static_cast<float>(height));
        for (int season = 0; season < CLIMATESEASONCOUNT; season++)
        {
            world->setseasonalrainfloat(season, c, 120.0f);   // mm/month
            world->setseasonaltempfloat(season, c, 4.0f);     // cool: little evaporation
        }
    }
    return world;
}

void runhydrology(planet& world)
{
    applyflowroutingstage(world);
    applyriverstage(world);
    applylakelevelstage(world);
    physical_layers::generate_drainage_basins(world);
}

int faceof(mesh::CellId c) { return static_cast<int>(c) / (facecells * facecells); }
}

TEST(HydrologyCube, RoutingRiversAndLakeLevelsRun)
{
    const std::unique_ptr<planet> world = slopingcubeworld(true);
    runhydrology(*world);
    const planet& view = *world;
    const mesh::CellGrid& grid = view.grid();
    const std::vector<std::uint8_t> ocean = flowrouting::oceancells(view);

    int land = 0, rivers = 0, lakecells = 0, unrouted = 0, nonfinite = 0;
    double landarea = 0.0, drainedarea = 0.0;
    for (mesh::CellId c = 0; c < view.cellcount(); c++)
    {
        if (!std::isfinite(view.discharge()[c]) || !std::isfinite(view.drainagearea()[c]) ||
            !std::isfinite(view.lakesurfacefloat(c)))
            nonfinite++;
        if (ocean[c] != 0)
            continue;
        land++;
        landarea += grid.area_m2(c);
        const int code = view.flowdirection()[c];
        if (code < 1 || code > mesh::slotcount || grid.neighbour(c, code) < 0)
            unrouted++;
        else if (ocean[static_cast<mesh::CellId>(grid.neighbour(c, code))] != 0)
            drainedarea += static_cast<double>(view.drainagearea()[c]) * 1.0e6;   // an outlet's catchment
        if (view.riverjan(c) + view.riverjul(c) > 0)
            rivers++;
        if (view.lakesurfacefloat(c) != 0.0f)
            lakecells++;
    }
    std::cout << "land " << land << ", rivers " << rivers << ", lake cells " << lakecells << "\n";
    EXPECT_EQ(nonfinite, 0);
    EXPECT_GT(land, 0);
    EXPECT_EQ(unrouted, 0) << "every land cell has a receiver";
    EXPECT_GT(rivers, 0);
    EXPECT_GT(lakecells, 0) << "the bowl holds a lake";
    // Every land cell drains to the sea through an outlet (lakes overflow at their spill):
    // the outlets' catchments sum to the land area.
    EXPECT_NEAR(drainedarea, landarea, 1.0e-4 * landarea);
}

TEST(HydrologyCube, RiverCrossesFaceEdge)
{
    const std::unique_ptr<planet> world = slopingcubeworld(false);
    runhydrology(*world);
    const planet& view = *world;
    const mesh::CellGrid& grid = view.grid();
    int crossings = 0, inflows = 0;
    for (mesh::CellId c = 0; c < view.cellcount(); c++)
    {
        const int code = view.riverdir(c);
        if (view.riverjan(c) + view.riverjul(c) == 0 || code < 1 || code > mesh::slotcount)
            continue;
        const std::int32_t next = grid.neighbour(c, code);
        if (next < 0)
            continue;
        const mesh::CellId n = static_cast<mesh::CellId>(next);
        if (faceof(n) != faceof(c) && view.riverjan(n) + view.riverjul(n) > 0)
        {
            crossings++;
            // The receiver sees the inflow through the reverse slot (the frames rotate across the edge).
            if (checkwaterinflow(*world, n) == 1)
                inflows++;
        }
    }
    std::cout << "rivers crossing a face edge: " << crossings << "\n";
    EXPECT_GT(crossings, 0);
    EXPECT_EQ(inflows, crossings);
}

// Tides on the cube: every sea cell next to land gets a range of at least twice the lunar
// strength, and casts across the southern ocean measure their fetch in km.
TEST(HydrologyCube, TidesReachEveryCoast)
{
    const std::unique_ptr<planet> world = slopingcubeworld(false);
    world->setlunar(1.0f);
    createtidalmap(*world);
    const planet& view = *world;
    const mesh::CellGrid& grid = view.grid();
    int coastal = 0, missing = 0, highest = 0;
    for (mesh::CellId c = 0; c < view.cellcount(); c++)
    {
        if (view.sea(c) != 1)
            continue;
        bool coast = false;
        for (int slot = 1; slot <= mesh::slotcount; slot++)
        {
            const std::int32_t neighbour = grid.neighbour(c, slot);
            if (neighbour >= 0 && view.sea(static_cast<mesh::CellId>(neighbour)) == 0)
                coast = true;
        }
        if (!coast)
            continue;
        coastal++;
        if (view.tide(c) < 2)
            missing++;
        highest = std::max(highest, view.tide(c));
    }
    std::cout << "coastal sea cells " << coastal << ", highest tide " << highest << "\n";
    EXPECT_GT(coastal, 0);
    EXPECT_EQ(missing, 0);
    EXPECT_GT(highest, 2);
}

TEST(HydrologyCube, BasinsTileTheLand)
{
    for (const bool bowl : { false, true })
    {
        const std::unique_ptr<planet> world = slopingcubeworld(bowl);
        runhydrology(*world);
        const planet& view = *world;
        const std::vector<std::uint8_t> ocean = flowrouting::oceancells(view);
        int unclassified = 0, classifiedsea = 0, endorheic = 0;
        for (mesh::CellId c = 0; c < view.cellcount(); c++)
        {
            const BasinClass basin = view.basinclass(c);
            if (ocean[c] != 0)
            {
                if (basin != BasinClass::none)
                    classifiedsea++;
                continue;
            }
            if (basin == BasinClass::none)
                unclassified++;
            if (basin == BasinClass::endorheic)
                endorheic++;
        }
        EXPECT_EQ(unclassified, 0) << "bowl " << bowl;
        EXPECT_EQ(classifiedsea, 0) << "bowl " << bowl;
        if (!bowl)
            EXPECT_EQ(endorheic, 0) << "a plain slope drains to the sea";
    }
}

namespace
{
std::uint64_t splitmix(std::uint64_t& state)
{
    std::uint64_t z = (state += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

double uniform(std::uint64_t& state) { return static_cast<double>(splitmix(state) >> 11) * (1.0 / 9007199254740992.0); }

// A direction-free surface: eight plane waves of random direction, wavenumber and phase on the
// unit sphere, the sea where their sum is below zero (about half the planet).
void noisesurface(planet& world, std::uint64_t seed)
{
    std::uint64_t state = seed;
    mesh::Vec3 direction[8];
    double wavenumber[8], phase[8];
    for (int wave = 0; wave < 8; wave++)
    {
        const double z = 2.0 * uniform(state) - 1.0;
        const double longitude = 2.0 * pi * uniform(state);
        const double r = std::sqrt(1.0 - z * z);
        direction[wave] = { r * std::cos(longitude), r * std::sin(longitude), z };
        wavenumber[wave] = 3.0 + 6.0 * uniform(state);
        phase[wave] = 2.0 * pi * uniform(state);
    }
    const mesh::CellGrid& grid = world.grid();
    for (mesh::CellId c = 0; c < world.cellcount(); c++)
    {
        const mesh::Vec3 p = grid.centre(c);
        double height = 0.0;
        for (int wave = 0; wave < 8; wave++)
            height += std::sin(wavenumber[wave] * mesh::dot(direction[wave], p) + phase[wave]) / wavenumber[wave];
        world.setnomfloat(c, static_cast<float>(sealevel + 3000.0 * height));
    }
}

// Over many direction-free surfaces, the mean of cos(2 b) of the flow bearings b of the land
// cells that drain downhill (outside depressions), per band of |latitude| 0-30, 30-60, 60-90:
// +1 for flows all north-south, -1 for flows all east-west, 0 without a preference.
std::vector<double> meridionalpreference(planet& world, int surfaces, std::vector<double>& counts)
{
    const mesh::CellGrid& grid = world.grid();
    std::vector<double> sum(3, 0.0);
    counts.assign(3, 0.0);
    for (int surface = 0; surface < surfaces; surface++)
    {
        noisesurface(world, 1000 + static_cast<std::uint64_t>(surface));
        applyflowroutingstage(world);
        const planet& view = world;
        const std::vector<std::uint8_t> ocean = flowrouting::oceancells(view);
        for (mesh::CellId c = 0; c < view.cellcount(); c++)
        {
            const int code = view.flowdirection()[c];
            if (ocean[c] != 0 || view.depressionid()[c] != 0 || code < 1 || code > mesh::slotcount)
                continue;
            const std::int32_t next = grid.neighbour(c, code);
            if (next < 0 || ocean[static_cast<mesh::CellId>(next)] != 0 ||
                !(view.nomfloat(static_cast<mesh::CellId>(next)) < view.nomfloat(c)))
                continue;
            const double latitude = std::fabs(latitudedeg(grid, c));
            const int band = std::min(2, static_cast<int>(latitude / 30.0));
            sum[band] += std::cos(2.0 * grid.bearing_rad(c, code));
            counts[band] += 1.0;
        }
    }
    for (int band = 0; band < 3; band++)
        sum[band] = counts[band] > 0.0 ? sum[band] / counts[band] : 0.0;
    return sum;
}
}

// The river-direction histogram by latitude band on direction-free surfaces: a preference
// for the face axes, or on the lat-lon grid for the short east-west spacing of the polar rows
// (its polar D8 bias), shows as a mean cos(2 b) away from 0. Measured (W2-H, 300 surfaces, seeds
// 1000-1299; seeds 5000-5299 alike within 0.02): the lat-lon grid of the same export width has
// -0.43 in the polar band and -0.07 to -0.09 at 30-60; the cube within 0.05 in every band. The
// cube's residue is D8 on its cells, which are 0.71 as wide as long at face-edge midpoints
// (zonal at 30-60, meridional near the equator), not a polar or face-axis preference; it is
// larger than W2-L's 2.5 % face/edge/corner spread of eroded volumes, a different statistic.
// Bounds: within 0.05 in every band, and the polar band within a tenth of the lat-lon one.
TEST(HydrologyCube, RiverDirectionsHaveNoFaceAxisBias)
{
    constexpr int surfaces = 300;
    const std::unique_ptr<planet> cube = slopingcubeworld(false);
    std::vector<double> cubecounts, latloncounts;
    const std::vector<double> cubepreference = meridionalpreference(*cube, surfaces, cubecounts);

    auto latlon = std::make_unique<planet>();
    latlon->setgrid(mesh::GridKind::latlon, 4 * facecells);
    latlon->setseed(1);
    latlon->setmaxelevation(48000);
    latlon->setsealevel(static_cast<int>(sealevel));
    latlon->seterosionduration(1.0f);
    for (mesh::CellId c = 0; c < latlon->cellcount(); c++)
        for (int season = 0; season < CLIMATESEASONCOUNT; season++)
        {
            latlon->setseasonalrainfloat(season, c, 120.0f);
            latlon->setseasonaltempfloat(season, c, 4.0f);
        }
    const std::vector<double> latlonpreference = meridionalpreference(*latlon, surfaces, latloncounts);

    for (int band = 0; band < 3; band++)
    {
        std::cout << "|latitude| " << 30 * band << "-" << 30 * band + 30 << ": cube mean cos 2b " << cubepreference[band]
                  << " (" << cubecounts[band] << " cells), lat-lon " << latlonpreference[band] << " ("
                  << latloncounts[band] << " cells)\n";
        ASSERT_GT(cubecounts[band], 1000.0) << "band " << band;
        EXPECT_LE(std::fabs(cubepreference[band]), 0.05) << "band " << band;
    }
    EXPECT_LE(std::fabs(cubepreference[2]), 0.1 * std::fabs(latlonpreference[2]));
}
