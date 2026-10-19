#include "core/deterministic_random.hpp"
#include "core/mesh/cell_grid.hpp"
#include "core/mesh/cube_geometry.hpp"
#include "core/mesh/mesh_geometry.hpp"
#include "core/mesh/mesh_remap.hpp"
#include "core/parallel_rows.hpp"
#include "simulations/landscape/sphere_metric.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <set>
#include <vector>
#include <gtest/gtest.h>

namespace
{
constexpr double pi = 3.14159265358979323846;
constexpr double radius = 6371000.0;
const double spherearea = 4.0 * pi * radius * radius;

double seconds(std::chrono::steady_clock::time_point start)
{
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

double smoothfield(const mesh::Vec3& p) { return p.x + p.z * p.z + 0.5 * p.y * p.z; }

std::vector<double> sampledfield(const mesh::CellGrid& grid)
{
    std::vector<double> values(grid.cellcount());
    for (mesh::CellId c = 0; c < grid.cellcount(); c++)
        values[c] = smoothfield(grid.centre(c));
    return values;
}

bool areneighbours(const mesh::CellGrid& grid, mesh::CellId a, mesh::CellId b)
{
    for (int k = 1; k <= mesh::slotcount; k++)
        if (grid.neighbour(a, k) == static_cast<std::int32_t>(b))
            return true;
    return false;
}

int faceof(const mesh::CellGrid& grid, mesh::CellId c)
{
    const int n = grid.facecells();
    return static_cast<int>(grid.globalkey(c) / (static_cast<mesh::LevelKey>(n) * static_cast<mesh::LevelKey>(n)));
}
}

TEST(Mesh, AreasSumToTheSphere)
{
    for (int n : { 1, 2, 32, 64 })
    {
        const auto grid = mesh::makecubegrid(n, radius);
        double total = 0.0;
        for (mesh::CellId c = 0; c < grid->cellcount(); c++)
            total += grid->area_m2(c);
        EXPECT_NEAR(total / spherearea, 1.0, 1.0e-12) << "cube N " << n;
        EXPECT_NEAR(grid->totalarea_m2() / spherearea, 1.0, 1.0e-12) << "cube N " << n;
    }
    const auto latlon = mesh::makelatlongrid(128, radius);
    double total = 0.0;
    for (mesh::CellId c = 0; c < latlon->cellcount(); c++)
        total += latlon->area_m2(c);
    EXPECT_NEAR(total / spherearea, 1.0, 1.0e-12);
    EXPECT_NEAR(latlon->totalarea_m2() / spherearea, 1.0, 1.0e-12);
}

TEST(Mesh, NeighboursAreSymmetricAndCornersHaveSeven)
{
    for (int n : { 1, 2, 3, 32 })
    {
        const auto grid = mesh::makecubegrid(n, radius);
        int seven = 0, eight = 0, other = 0;
        // Every distance against the cell's smallest cardinal spacing.
        double worstcardinal = 0.0, worstdiagonal = 0.0;
        for (mesh::CellId c = 0; c < grid->cellcount(); c++)
        {
            int count = 0;
            double cardinal = grid->distance_m(c, 1);
            for (int k = 3; k <= mesh::slotcount; k += 2)
                cardinal = std::min(cardinal, grid->distance_m(c, k));
            for (int k = 1; k <= mesh::slotcount; k++)
            {
                double& worst = k % 2 == 1 ? worstcardinal : worstdiagonal;
                worst = std::max(worst, grid->distance_m(c, k) / cardinal);
            }
            for (int k = 1; k <= mesh::slotcount; k++)
            {
                const std::int32_t m = grid->neighbour(c, k);
                if (m < 0)
                {
                    EXPECT_EQ(grid->reverse(c, k), 0);
                    continue;
                }
                count++;
                const int back = grid->reverse(c, k);
                ASSERT_GE(back, 1) << "N " << n << " cell " << c << " slot " << k;
                EXPECT_EQ(grid->neighbour(static_cast<mesh::CellId>(m), back), static_cast<std::int32_t>(c))
                    << "N " << n << " cell " << c << " slot " << k;
                EXPECT_GT(grid->distance_m(c, k), 0.0);
                EXPECT_EQ(grid->facelength_m(c, k) > 0.0, k % 2 == 1);
            }
            if (count == 7)
                seven++;
            else if (count == 8)
                eight++;
            else
                other++;
            if (n == 1)
                EXPECT_EQ(count, 4);
        }
        std::cout << "[Mesh] cube N " << n << " largest distance / smallest cardinal: cardinal " << worstcardinal
                  << ", diagonal " << worstdiagonal << "\n";
        EXPECT_LE(worstcardinal, 1.6) << "N " << n;
        // Near a cube corner the cells are rhombi with a 120-degree angle at the corner: the long diagonal is
        // about sqrt(3) times the shorter side, and more where the sides differ.
        EXPECT_LE(worstdiagonal, 2.0) << "N " << n;
        if (n >= 2)
        {
            EXPECT_EQ(seven, 24) << "N " << n;
            EXPECT_EQ(other, 0) << "N " << n;
            EXPECT_EQ(static_cast<std::size_t>(eight), grid->cellcount() - 24) << "N " << n;
        }
    }
    const auto latlon = mesh::makelatlongrid(128, radius);
    for (mesh::CellId c = 0; c < latlon->cellcount(); c++)
        for (int k = 1; k <= mesh::slotcount; k++)
        {
            const std::int32_t m = latlon->neighbour(c, k);
            if (m < 0)
                continue;
            EXPECT_EQ(latlon->neighbour(static_cast<mesh::CellId>(m), latlon->reverse(c, k)), static_cast<std::int32_t>(c));
        }
}

TEST(Mesh, EdgeAndDiagonalNeighboursAcrossFaces)
{
    // The uw-cube-1 table of sphere-mesh.md section 2: centre c, right u, up v per face.
    const mesh::Vec3 table[mesh::cubefacecount][3] = {
        { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } },
        { { 0, 1, 0 }, { -1, 0, 0 }, { 0, 0, 1 } },
        { { -1, 0, 0 }, { 0, -1, 0 }, { 0, 0, 1 } },
        { { 0, -1, 0 }, { 1, 0, 0 }, { 0, 0, 1 } },
        { { 0, 0, 1 }, { 0, 1, 0 }, { -1, 0, 0 } },
        { { 0, 0, -1 }, { 0, 1, 0 }, { 1, 0, 0 } },
    };
    for (int f = 0; f < mesh::cubefacecount; f++)
    {
        const mesh::CubeFace& face = mesh::cubeface(f);
        const mesh::Vec3 actual[3] = { face.centre, face.right, face.up };
        for (int axis = 0; axis < 3; axis++)
        {
            EXPECT_EQ(actual[axis].x, table[f][axis].x) << "face " << f << " axis " << axis;
            EXPECT_EQ(actual[axis].y, table[f][axis].y) << "face " << f << " axis " << axis;
            EXPECT_EQ(actual[axis].z, table[f][axis].z) << "face " << f << " axis " << axis;
        }
    }
    // Face 0's top row meets face 4's bottom row and face 0's bottom row face 5's top row, with the same u (same i);
    // face 0's east column meets face 1's west column.
    const int n = 8;
    const auto grid = mesh::makecubegrid(n, radius);
    for (int i = 0; i < n; i++)
    {
        EXPECT_EQ(grid->neighbour(mesh::cubeid(n, 0, i, 0), 1), static_cast<std::int32_t>(mesh::cubeid(n, 4, i, n - 1)));
        EXPECT_EQ(grid->neighbour(mesh::cubeid(n, 0, i, n - 1), 5), static_cast<std::int32_t>(mesh::cubeid(n, 5, i, 0)));
        EXPECT_EQ(grid->neighbour(mesh::cubeid(n, 5, i, 0), 1), static_cast<std::int32_t>(mesh::cubeid(n, 0, i, n - 1)));
        EXPECT_EQ(grid->neighbour(mesh::cubeid(n, 0, n - 1, i), 3), static_cast<std::int32_t>(mesh::cubeid(n, 1, 0, i)));
    }
    for (int j = 1; j < n; j++)
        EXPECT_EQ(grid->neighbour(mesh::cubeid(n, 0, n - 1, j), 2), static_cast<std::int32_t>(mesh::cubeid(n, 1, 0, j - 1)));
    EXPECT_EQ(grid->neighbour(mesh::cubeid(n, 0, n - 1, 0), 2), -1);
}

TEST(Mesh, LocateFindsEveryCentre)
{
    for (int n : { 1, 2, 3, 32, 64 })
    {
        const auto grid = mesh::makecubegrid(n, radius);
        for (mesh::CellId c = 0; c < grid->cellcount(); c++)
            ASSERT_EQ(grid->locate(grid->centre(c)), c) << "cube N " << n;
    }
    const auto latlon = mesh::makelatlongrid(128, radius);
    for (mesh::CellId c = 0; c < latlon->cellcount(); c++)
        ASSERT_EQ(latlon->locate(latlon->centre(c)), c);
    // Longitudes run east from -180 degrees, row 0 is the north.
    EXPECT_NEAR(latlon->longitude_rad(0), -pi + pi / 128.0, 1.0e-15);
    EXPECT_GT(latlon->latitude_rad(0), latlon->latitude_rad(latlon->cellof(0, 1)));
}

TEST(Mesh, LatLonTablesReproduceSphereMetric)
{
    for (int w : { 128, 256 })
    {
        const auto grid = mesh::makelatlongrid(w, radius);
        const landscape::SphereMetric metric(w, w / 2, radius);
        ASSERT_EQ(grid->cellcount(), metric.size());
        EXPECT_EQ(grid->totalarea_m2(), metric.totalarea_m2());
        for (mesh::CellId c = 0; c < grid->cellcount(); c++)
        {
            ASSERT_EQ(grid->area_m2(c), metric.area_m2(c));
            ASSERT_EQ(grid->latitude_rad(c), metric.rowlatitude_rad[static_cast<std::size_t>(metric.row(c))]);
            ASSERT_EQ(grid->column(c), metric.column(c));
            ASSERT_EQ(grid->row(c), metric.row(c));
            for (int k = 1; k <= mesh::slotcount; k++)
            {
                const std::size_t expected = metric.neighbour(c, k);
                const std::int32_t actual = grid->neighbour(c, k);
                if (expected == c)
                    ASSERT_EQ(actual, -1);
                else
                    ASSERT_EQ(actual, static_cast<std::int32_t>(expected));
                ASSERT_EQ(grid->distance_m(c, k), metric.distance_m(c, k));
                ASSERT_EQ(grid->facelength_m(c, k), metric.facelength_m(c, k));
                if (actual >= 0)
                    ASSERT_EQ(grid->reverse(c, k), landscape::oppositedirection(k));
            }
        }
    }
}

TEST(Mesh, CubeAreaRatioAndBuildTime)
{
    const auto grid = mesh::makecubegrid(64, radius);
    double smallest = grid->area_m2(0), largest = grid->area_m2(0);
    for (mesh::CellId c = 0; c < grid->cellcount(); c++)
    {
        smallest = std::min(smallest, grid->area_m2(c));
        largest = std::max(largest, grid->area_m2(c));
    }
    const double ratio = largest / smallest;
    std::cout << "[Mesh] cube N 64 area largest / smallest " << ratio << "\n";
    EXPECT_NEAR(ratio, 1.40, 0.01);

    const auto start = std::chrono::steady_clock::now();
    const auto big = mesh::makecubegrid(128, radius);
    const double elapsed = seconds(start);
    std::cout << "[Mesh] cube N 128 tables built in " << elapsed << " s\n";
    ::testing::Test::RecordProperty("cube128_build_s", std::to_string(elapsed));
    EXPECT_EQ(big->cellcount(), 6u * 128u * 128u);
}

TEST(Mesh, RemapKeepsConstantsTotalsAndSmoothFields)
{
    const auto cube64 = mesh::makecubegrid(64, radius);
    const auto latlon256 = mesh::makelatlongrid(256, radius);
    const auto latlon512 = mesh::makelatlongrid(512, radius);
    const auto latlon2048 = mesh::makelatlongrid(2048, radius);
    struct Case
    {
        const char* name;
        const mesh::CellGrid* source;
        const mesh::CellGrid* destination;
    };
    const Case cases[] = {
        { "cube64->latlon256", cube64.get(), latlon256.get() },
        { "latlon256->cube64", latlon256.get(), cube64.get() },
        { "latlon2048->cube64", latlon2048.get(), cube64.get() },
        { "cube64->latlon512", cube64.get(), latlon512.get() },
        { "latlon512->latlon256", latlon512.get(), latlon256.get() },
        { "latlon256->latlon512", latlon256.get(), latlon512.get() },
    };
    for (const Case& test : cases)
    {
        const mesh::CellGrid& source = *test.source;
        const mesh::CellGrid& destination = *test.destination;
        const auto start = std::chrono::steady_clock::now();
        const mesh::RemapWeights weights = mesh::makeremap(source, destination);
        const double elapsed = seconds(start);
        std::cout << "[Mesh] " << test.name << " weights in " << elapsed << " s, " << weights.source.size()
                  << " entries, " << weights.emptyrows << " empty rows\n";
        ::testing::Test::RecordProperty(std::string(test.name) + "_weights_s", std::to_string(elapsed));

        // Column sums are the source areas.
        std::vector<double> columns(source.cellcount(), 0.0);
        for (std::size_t e = 0; e < weights.source.size(); e++)
            columns[weights.source[e]] += weights.area_m2[e];
        double worstcolumn = 0.0;
        for (mesh::CellId s = 0; s < source.cellcount(); s++)
            worstcolumn = std::max(worstcolumn, std::fabs(columns[s] / source.area_m2(s) - 1.0));
        EXPECT_LT(worstcolumn, 1.0e-12) << test.name;
        // Rows ascending by source.
        for (std::size_t d = 0; d < weights.destinationcount; d++)
            for (std::size_t e = weights.rowstart[d] + 1; e < weights.rowstart[d + 1]; e++)
                ASSERT_LT(weights.source[e - 1], weights.source[e]) << test.name;

        const std::vector<double> constant(source.cellcount(), 3.7);
        const auto mean = mesh::remapmean(weights, source, destination, constant);
        for (double value : mean)
            ASSERT_EQ(value, 3.7) << test.name;

        std::vector<double> mass(source.cellcount());
        double sourcetotal = 0.0;
        for (mesh::CellId s = 0; s < source.cellcount(); s++)
        {
            mass[s] = 2.0 + smoothfield(source.centre(s));
            sourcetotal += mass[s];
        }
        const auto moved = mesh::remaptotal(weights, mass);
        double destinationtotal = 0.0;
        for (double value : moved)
            destinationtotal += value;
        EXPECT_NEAR(destinationtotal / sourcetotal, 1.0, 1.0e-12) << test.name;

        const auto smooth = mesh::remapmean(weights, source, destination, sampledfield(source));
        double sum = 0.0, worst = 0.0;
        for (mesh::CellId d = 0; d < destination.cellcount(); d++)
        {
            const double error = std::fabs(smooth[d] - smoothfield(destination.centre(d)));
            sum += error;
            worst = std::max(worst, error);
        }
        const double meanerror = sum / destination.cellcount();
        std::cout << "[Mesh] " << test.name << " smooth field mean |error| " << meanerror << ", max " << worst << "\n";
        EXPECT_LT(meanerror, 1.0e-2) << test.name;
        EXPECT_LT(worst, 5.0e-2) << test.name;
    }
}

namespace
{
std::uint64_t weightshash(const mesh::RemapWeights& weights)
{
    std::uint64_t hash = 1469598103934665603ull;  // FNV-1a over the bytes of the three arrays
    const auto bytes = [&hash](const void* data, std::size_t size) {
        const unsigned char* p = static_cast<const unsigned char*>(data);
        for (std::size_t n = 0; n < size; n++)
            hash = (hash ^ p[n]) * 1099511628211ull;
    };
    bytes(weights.rowstart.data(), weights.rowstart.size() * sizeof(std::size_t));
    bytes(weights.source.data(), weights.source.size() * sizeof(mesh::CellId));
    bytes(weights.area_m2.data(), weights.area_m2.size() * sizeof(double));
    return hash;
}
}

TEST(Mesh, RemapWeightsDoNotDependOnTheThreadCount)
{
    // parallelforrows runs a call on the calling thread alone while
    // rowworkers::executing is set (its nested-call path): the one-thread build.
    const auto cube = mesh::makecubegrid(64, radius);
    const auto latlon = mesh::makelatlongrid(256, radius);
    for (const auto& pair : { std::make_pair(cube.get(), latlon.get()), std::make_pair(latlon.get(), cube.get()) })
    {
        rowworkers::executing = true;
        const mesh::RemapWeights single = mesh::makeremap(*pair.first, *pair.second);
        rowworkers::executing = false;
        const mesh::RemapWeights parallel = mesh::makeremap(*pair.first, *pair.second);
        std::cout << "[Mesh] weights hash, 1 thread " << std::hex << weightshash(single) << std::dec << ", "
                  << rowworkers::concurrency() << " threads " << std::hex << weightshash(parallel) << std::dec << "\n";
        EXPECT_EQ(weightshash(single), weightshash(parallel));
        EXPECT_EQ(single.rowstart, parallel.rowstart);
        EXPECT_EQ(single.source, parallel.source);
        EXPECT_EQ(single.area_m2, parallel.area_m2);
    }
}

TEST(Mesh, RemapCubeToFineLatLonTiming)
{
    const auto cube = mesh::makecubegrid(64, radius);
    const auto latlon = mesh::makelatlongrid(2048, radius);
    const auto start = std::chrono::steady_clock::now();
    const mesh::RemapWeights weights = mesh::makeremap(*cube, *latlon);
    const double elapsed = seconds(start);
    std::cout << "[Mesh] cube64->latlon2048 weights in " << elapsed << " s, " << weights.source.size() << " entries, "
              << weights.emptyrows << " empty rows\n";
    ::testing::Test::RecordProperty("cube64->latlon2048_weights_s", std::to_string(elapsed));
    const auto mean = mesh::remapmean(weights, *cube, *latlon, std::vector<double>(cube->cellcount(), 3.7));
    for (double value : mean)
        ASSERT_EQ(value, 3.7);
}

TEST(Mesh, MajorityAndNearestCodes)
{
    const auto cube = mesh::makecubegrid(16, radius);
    const auto latlon = mesh::makelatlongrid(32, radius);
    std::vector<int> codes(latlon->cellcount());
    for (mesh::CellId c = 0; c < latlon->cellcount(); c++)
        codes[c] = latlon->centre(c).z > 0.0 ? 2 : 1;
    const auto weights = mesh::makeremap(*latlon, *cube);
    const auto majority = mesh::remapmajority(weights, *latlon, *cube, codes);
    const auto nearest = mesh::remapnearest(*latlon, *cube, codes);
    for (mesh::CellId c = 0; c < cube->cellcount(); c++)
    {
        const double z = cube->centre(c).z;
        if (std::fabs(z) > 0.1)
        {
            EXPECT_EQ(majority[c], z > 0.0 ? 2 : 1);
            EXPECT_EQ(nearest[c], z > 0.0 ? 2 : 1);
        }
    }
}

TEST(Mesh, IntegerTotalsRoundToNearest)
{
    const auto cube = mesh::makecubegrid(16, radius);
    const auto latlon = mesh::makelatlongrid(32, radius);
    const auto weights = mesh::makeremap(*cube, *latlon);
    std::vector<int> counts(cube->cellcount());
    std::vector<double> exact(cube->cellcount());
    long long sourcetotal = 0;
    for (mesh::CellId c = 0; c < cube->cellcount(); c++)
    {
        counts[c] = 1000 + static_cast<int>(c % 7);
        exact[c] = counts[c];
        sourcetotal += counts[c];
    }
    const std::vector<int> moved = mesh::remaptotal(weights, counts);
    const std::vector<double> reference = mesh::remaptotal(weights, exact);
    long long movedtotal = 0;
    for (mesh::CellId d = 0; d < latlon->cellcount(); d++)
    {
        // Rounded, not truncated: a share of 999.7 is 1000.
        EXPECT_EQ(moved[d], static_cast<int>(std::llround(reference[d]))) << d;
        movedtotal += moved[d];
    }
    // The integer total is kept up to half a unit per destination cell.
    EXPECT_LE(std::llabs(movedtotal - sourcetotal), static_cast<long long>(latlon->cellcount() / 2));
}

TEST(Mesh, SampleBilinearIsContinuousAcrossAFaceEdgeAndACorner)
{
    const int n = 32;
    const auto grid = mesh::makecubegrid(n, radius);
    const std::vector<double> field = sampledfield(*grid);
    // The largest difference between neighbouring centres: a seam would jump by about that much.
    double cellstep = 0.0;
    for (mesh::CellId c = 0; c < grid->cellcount(); c++)
        for (int k = 1; k <= mesh::slotcount; k += 2)
            cellstep = std::max(cellstep, std::fabs(field[static_cast<std::size_t>(grid->neighbour(c, k))] - field[c]));
    const double sub = (pi / (2.0 * n)) / 50.0;  // a fiftieth of a cell
    struct Path
    {
        const char* name;
        mesh::Vec3 origin, direction;
    };
    const mesh::Vec3 corner{ 1.0, 1.0, 1.0 };
    const Path paths[] = {
        { "face edge 0-1", mesh::unitvector(0.2, pi / 4.0), { -1.0, 1.0, 0.0 } },
        { "face edge 0-4", mesh::unitvector(pi / 4.0 - 0.0, 0.1), { -1.0, 0.0, 1.0 } },
        { "corner 0-1-4 (through)", corner, { 1.0, -1.0, 0.0 } },
        { "corner 0-1-4 (from face 4)", corner, { 1.0, 1.0, -2.0 } },
    };
    for (const Path& path : paths)
    {
        const mesh::Vec3 origin = mesh::normalise(path.origin);
        const mesh::Vec3 tangent = mesh::normalise(path.direction - mesh::dot(path.direction, origin) * origin);
        std::set<int> faces;
        double previous = 0.0, jump = 0.0;
        const int steps = 400;
        for (int s = -steps; s <= steps; s++)
        {
            const double angle = s * sub;
            const mesh::Vec3 p = std::cos(angle) * origin + std::sin(angle) * tangent;
            faces.insert(mesh::cubeposition(p).face);
            const double value = mesh::samplebilinear(*grid, field, p);
            if (s > -steps)
                jump = std::max(jump, std::fabs(value - previous));
            previous = value;
        }
        std::cout << "[Mesh] samplebilinear " << path.name << ": max jump " << jump << " (cell step " << cellstep << ")\n";
        EXPECT_GE(faces.size(), 2u) << path.name;
        EXPECT_LT(jump, 0.1 * cellstep) << path.name;
    }
    // At a cube corner the sample is the mean of the three corner cells.
    const mesh::Vec3 p = mesh::normalise(corner);
    const mesh::BilinearStencil stencil = mesh::bilinearstencil(*grid, p);
    ASSERT_EQ(stencil.count, 3);
    for (int m = 0; m < 3; m++)
        EXPECT_NEAR(stencil.weight[m], 1.0 / 3.0, 1.0e-9);
}

TEST(Mesh, CellsWithinAndAlongArc)
{
    const int n = 32;
    const auto grid = mesh::makecubegrid(n, radius);
    for (const mesh::Vec3& centre : { mesh::normalise({ 1.0, 1.0, 1.0 }), mesh::unitvector(0.1, pi / 4.0) })
    {
        const double reach = 900000.0;
        const auto cells = mesh::cellswithin(*grid, centre, reach);
        std::vector<mesh::CellId> brute;
        for (mesh::CellId c = 0; c < grid->cellcount(); c++)
            if (mesh::greatcircle_m(centre, grid->centre(c), radius) <= reach)
                brute.push_back(c);
        EXPECT_EQ(cells, brute);
        EXPECT_GT(cells.size(), 20u);
    }
    const mesh::Vec3 p = mesh::unitvector(0.15, 0.4), q = mesh::unitvector(0.35, 1.3);
    const auto arc = mesh::cellsalongarc(*grid, p, q);
    ASSERT_GE(arc.size(), 2u);
    EXPECT_EQ(arc.front(), grid->locate(p));
    EXPECT_EQ(arc.back(), grid->locate(q));
    std::set<int> faces;
    for (std::size_t m = 0; m < arc.size(); m++)
    {
        faces.insert(faceof(*grid, arc[m]));
        if (m > 0)
            EXPECT_TRUE(areneighbours(*grid, arc[m - 1], arc[m])) << "arc cells " << m - 1 << ", " << m;
    }
    EXPECT_GE(faces.size(), 2u);
    // Through a cube corner as well.
    const auto cornerarc = mesh::cellsalongarc(*grid, mesh::normalise({ 1.0, 0.7, 1.3 }), mesh::normalise({ 0.8, 1.3, 0.75 }));
    for (std::size_t m = 1; m < cornerarc.size(); m++)
        EXPECT_TRUE(areneighbours(*grid, cornerarc[m - 1], cornerarc[m]));
}

TEST(Mesh, SubgridKeepsInsideNeighboursOnly)
{
    const auto cube = mesh::makecubegrid(8, radius);
    const auto latlon = mesh::makelatlongrid(32, radius);
    for (const mesh::CellGrid* base : { cube.get(), latlon.get() })
    {
        auto ids = mesh::cellswithin(*base, mesh::normalise({ 1.0, 1.0, 0.9 }), 2500000.0);
        std::reverse(ids.begin(), ids.end());
        ids.push_back(ids.front());  // duplicates and order do not matter
        const auto window = mesh::subgrid(*base, ids);
        ASSERT_TRUE(window->issubgrid());
        ASSERT_EQ(window->cellcount(), ids.size() - 1);
        for (mesh::CellId local = 1; local < window->cellcount(); local++)
            EXPECT_LT(window->globalkey(local - 1), window->globalkey(local));
        double total = 0.0;
        for (mesh::CellId local = 0; local < window->cellcount(); local++)
        {
            const mesh::CellId global = static_cast<mesh::CellId>(window->globalkey(local));
            total += window->area_m2(local);
            EXPECT_EQ(window->locate(window->centre(local)), local);
            for (int k = 1; k <= mesh::slotcount; k++)
            {
                const std::int32_t outside = base->neighbour(global, k);
                const std::int32_t inside = window->neighbour(local, k);
                const mesh::CellId expected = outside < 0 ? mesh::nocell : window->localid(static_cast<mesh::CellId>(outside));
                if (expected == mesh::nocell)
                    EXPECT_EQ(inside, -1);
                else
                {
                    EXPECT_EQ(inside, static_cast<std::int32_t>(expected));
                    EXPECT_EQ(window->neighbour(static_cast<mesh::CellId>(inside), window->reverse(local, k)),
                        static_cast<std::int32_t>(local));
                }
            }
        }
        EXPECT_NEAR(window->totalarea_m2(), total, 1.0e-6 * total);
    }
}

TEST(Mesh, TangentFrameAndRandomCellByArea)
{
    const auto grid = mesh::makecubegrid(16, radius);
    const mesh::TangentFrame frame = mesh::tangentframe(mesh::unitvector(0.5, 1.0), radius);
    const mesh::Vec3 east = mesh::fromtangent(frame, 100.0, 0.0);
    EXPECT_NEAR(mesh::greatcircle_m(frame.origin, east, radius), radius * std::atan(100000.0 / radius), 1.0e-6);
    EXPECT_NEAR(mesh::latitudeof(mesh::fromtangent(frame, 0.0, 100.0)) - 0.5, std::atan(100000.0 / radius), 1.0e-12);
    const mesh::EastNorth components = mesh::tangenttoeastnorth(*grid, 5, grid->east(5));
    EXPECT_NEAR(components.east, 1.0, 1.0e-12);
    EXPECT_NEAR(components.north, 0.0, 1.0e-12);

    std::uint64_t state = 1;
    const auto stream = [&state]() {
        state = deterministicmix64(state);
        return mesh::uniformfrombits(state);
    };
    std::vector<int> counts(grid->cellcount(), 0);
    const int draws = 200000;
    std::vector<mesh::CellId> first;
    for (int d = 0; d < draws; d++)
    {
        const mesh::CellId c = mesh::randomcellbyarea(*grid, stream);
        ASSERT_LT(c, grid->cellcount());
        counts[c]++;
        if (d < 10)
            first.push_back(c);
    }
    // Poleward cube cells and face-centre cells differ in area by 1.4; the draw follows the area.
    double corner = 0.0, centre = 0.0;
    for (int f = 0; f < mesh::cubefacecount; f++)
    {
        corner += counts[mesh::cubeid(16, f, 0, 0)] + counts[mesh::cubeid(16, f, 15, 15)];
        centre += counts[mesh::cubeid(16, f, 7, 7)] + counts[mesh::cubeid(16, f, 8, 8)];
    }
    const double expected = grid->area_m2(mesh::cubeid(16, 0, 7, 7)) / grid->area_m2(mesh::cubeid(16, 0, 0, 0));
    EXPECT_NEAR(centre / corner, expected, 0.15);
    state = 1;
    for (int d = 0; d < 10; d++)
        EXPECT_EQ(mesh::randomcellbyarea(*grid, stream), first[static_cast<std::size_t>(d)]);
}
