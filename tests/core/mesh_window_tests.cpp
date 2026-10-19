// Window grids on level keys (docs/architecture/sphere-mesh.md section 8, T8 W2-W1a):
// 64-bit cube geometry, cubewindowgrid (stored and compact), latlonwindowgrid
// against the planet-wide grid's rows (landscape::SphereMetric), parents and children, and the helpers C,
// G, S and R of window_grid.hpp.

#include "core/mesh/cell_grid.hpp"
#include "core/mesh/cube_geometry.hpp"
#include "core/mesh/window_grid.hpp"
#include "simulations/landscape/regional_grid.hpp"
#include "simulations/landscape/sphere_metric.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>
#include <gtest/gtest.h>

namespace
{
constexpr double pi = 3.14159265358979323846;
constexpr double radius = 6371000.0;

std::uint64_t bits(double value)
{
    std::uint64_t result = 0;
    std::memcpy(&result, &value, sizeof result);
    return result;
}

void expectsamevector(const mesh::Vec3& p, const mesh::Vec3& q)
{
    EXPECT_EQ(bits(p.x), bits(q.x));
    EXPECT_EQ(bits(p.y), bits(q.y));
    EXPECT_EQ(bits(p.z), bits(q.z));
}

std::vector<mesh::LevelKey> randomkeys(std::uint64_t seed, std::uint64_t limit, std::size_t count)
{
    std::mt19937_64 random(seed);
    std::vector<mesh::LevelKey> keys(count);
    for (auto& key : keys)
        key = random() % limit;
    return keys;
}

std::vector<mesh::CellId> cellids(const std::vector<mesh::LevelKey>& keys)
{
    std::vector<mesh::CellId> ids;
    for (mesh::LevelKey key : keys)
        ids.push_back(static_cast<mesh::CellId>(key));
    return ids;
}

// Keys within `reach` cells of every face corner, and a few along each edge's middle.
std::vector<mesh::LevelKey> cornerkeys(std::int64_t n, int reach)
{
    std::vector<mesh::LevelKey> keys;
    for (int face = 0; face < mesh::cubefacecount; face++)
        for (int ci = 0; ci < 2; ci++)
            for (int cj = 0; cj < 2; cj++)
                for (int d = 0; d < reach; d++)
                    for (int e = 0; e < reach; e++)
                    {
                        const std::int64_t i = ci == 0 ? d : n - 1 - d, j = cj == 0 ? e : n - 1 - e;
                        keys.push_back(mesh::cubekey(n, face, i, j));
                        keys.push_back(mesh::cubekey(n, face, n / 2 + d, cj == 0 ? e : n - 1 - e));
                    }
    std::sort(keys.begin(), keys.end());
    keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
    return keys;
}

// Every accessor of two grids over the same cells, bit for bit.
void expectsamegrid(const mesh::CellGrid& a, const mesh::CellGrid& b)
{
    ASSERT_EQ(a.cellcount(), b.cellcount());
    EXPECT_EQ(bits(a.totalarea_m2()), bits(b.totalarea_m2()));
    EXPECT_EQ(bits(a.minimumspacing_m()), bits(b.minimumspacing_m()));
    for (mesh::CellId c = 0; c < a.cellcount(); c++)
    {
        EXPECT_EQ(a.globalkey(c), b.globalkey(c));
        EXPECT_EQ(a.globalid(c), b.globalid(c));
        expectsamevector(a.centre(c), b.centre(c));
        EXPECT_EQ(bits(a.latitude_rad(c)), bits(b.latitude_rad(c)));
        EXPECT_EQ(bits(a.longitude_rad(c)), bits(b.longitude_rad(c)));
        EXPECT_EQ(bits(a.area_m2(c)), bits(b.area_m2(c)));
        for (int k = 1; k <= mesh::slotcount; k++)
        {
            EXPECT_EQ(a.neighbour(c, k), b.neighbour(c, k));
            EXPECT_EQ(a.reverse(c, k), b.reverse(c, k));
            EXPECT_EQ(bits(a.distance_m(c, k)), bits(b.distance_m(c, k)));
            EXPECT_EQ(bits(a.facelength_m(c, k)), bits(b.facelength_m(c, k)));
            EXPECT_EQ(bits(a.bearing_rad(c, k)), bits(b.bearing_rad(c, k)));
        }
    }
}

double smoothfield(const mesh::Vec3& p) { return std::sin(2.0 * p.x + 0.5) * std::cos(1.5 * p.y) + p.z * p.z * p.x + 0.3 * p.y; }
}

TEST(MeshWindow, KeysEqualCubeGridIds)
{
    for (int n : { 4, 16, 64 })
    {
        const auto grid = mesh::makecubegrid(n, radius);
        for (mesh::CellId c = 0; c < grid->cellcount(); c++)
        {
            int face = 0, i = 0, j = 0;
            mesh::cubeindices(n, c, face, i, j);
            ASSERT_EQ(mesh::cubekey(n, face, i, j), c);
            const mesh::CubeCell cell = mesh::cubecell(n, c);
            EXPECT_EQ(cell.face, face);
            EXPECT_EQ(cell.column, i);
            EXPECT_EQ(cell.row, j);
            expectsamevector(mesh::cubekeycentre(n, c), grid->centre(c));
            EXPECT_EQ(bits(radius * radius * mesh::cubekeyunitarea(n, c)), bits(grid->area_m2(c)));
            EXPECT_EQ(mesh::cubelocatekey(n, grid->centre(c)), c);
            for (int k = 1; k <= mesh::slotcount; k++)
            {
                const mesh::LevelKey other = mesh::cubekeyneighbour(n, c, k);
                const std::int32_t expected = grid->neighbour(c, k);
                EXPECT_EQ(other == mesh::nokey ? -1 : static_cast<std::int64_t>(other), expected) << "N " << n << " cell " << c << " slot " << k;
                EXPECT_EQ(mesh::cubekeyreverse(n, c, k), grid->reverse(c, k));
                EXPECT_EQ(bits(mesh::cubekeydistance_m(n, c, k, radius)), bits(grid->distance_m(c, k)));
                EXPECT_EQ(bits(mesh::cubekeyfacelength_m(n, c, k, radius)), bits(grid->facelength_m(c, k)));
                EXPECT_EQ(bits(mesh::cubekeybearing_rad(n, c, k)), bits(grid->bearing_rad(c, k)));
            }
        }
    }
}

TEST(MeshWindow, LevelSpacingIsTheCubeGridsMinimum)
{
    for (int n = 1; n <= 64; n++)
        EXPECT_EQ(bits(mesh::cubeminimumspacing_m(n, radius)), bits(mesh::makecubegrid(n, radius)->minimumspacing_m())) << "N " << n;
}

TEST(MeshWindow, RandomSubsetsEqualSubgrid)
{
    for (int n : { 4, 16, 64 })
    {
        const auto base = mesh::makecubegrid(n, radius);
        for (std::uint64_t trial = 0; trial < 3; trial++)
        {
            const auto keys = randomkeys(100 + trial, base->cellcount(), base->cellcount() / 3);
            const std::vector<mesh::CellId> ids = cellids(keys);
            const auto sub = mesh::subgrid(*base, ids);
            for (bool compact : { false, true })
            {
                const auto window = mesh::cubewindowgrid(n, radius, keys, compact);
                EXPECT_EQ(window->iscompact(), compact);
                EXPECT_TRUE(window->issubgrid());
                expectsamegrid(*window, *sub);
                for (mesh::CellId c = 0; c < window->cellcount(); c++)
                {
                    EXPECT_EQ(window->localofkey(window->globalkey(c)), c);
                    EXPECT_EQ(window->localid(window->globalid(c)), c);
                    EXPECT_EQ(window->locate(window->centre(c)), c);
                }
                for (mesh::CellId g = 0; g < base->cellcount(); g++)
                    EXPECT_EQ(window->localofkey(g), sub->localid(g));
            }
        }
    }
}

TEST(MeshWindow, NeighbourReverseIdentityNearCornersAtLevel28)
{
    const std::int64_t n = std::int64_t{ 1 } << 28;
    const double step = mesh::cubestep(n);
    int corners = 0;
    for (mesh::LevelKey key : cornerkeys(n, 4))
    {
        int count = 0;
        for (int k = 1; k <= mesh::slotcount; k++)
        {
            const mesh::LevelKey other = mesh::cubekeyneighbour(n, key, k);
            if (other == mesh::nokey)
            {
                EXPECT_EQ(k % 2, 0) << "only a diagonal may be missing";
                continue;
            }
            count++;
            const int r = mesh::cubekeyreverse(n, key, k);
            ASSERT_NE(r, 0) << "key " << key << " slot " << k;
            EXPECT_EQ(mesh::cubekeyneighbour(n, other, r), key);
            EXPECT_EQ(mesh::cubekeyreverse(n, other, r), k);
            const double angle = mesh::cubekeydistance_m(n, key, k, 1.0);
            EXPECT_GT(angle, 0.5 * step);
            EXPECT_LT(angle, (k % 2 == 1 ? 1.2 : 2.0) * step);
        }
        const mesh::CubeCell cell = mesh::cubecell(n, key);
        const bool corner = (cell.column == 0 || cell.column == n - 1) && (cell.row == 0 || cell.row == n - 1);
        EXPECT_EQ(count, corner ? 7 : 8);
        corners += corner ? 1 : 0;
    }
    EXPECT_EQ(corners, 24);
}

TEST(MeshWindow, ChildrenAreasSumToTheParent)
{
    constexpr int factor = 4;
    double worst = 0.0;
    for (std::int64_t n : { std::int64_t{ 4 }, std::int64_t{ 16 }, std::int64_t{ 64 }, std::int64_t{ 1024 },
                            std::int64_t{ 1 } << 20, std::int64_t{ 1 } << 28 })
    {
        std::vector<mesh::LevelKey> keys = cornerkeys(n, 2);
        const auto extra = randomkeys(7, static_cast<std::uint64_t>(6 * n * n), 50);
        keys.insert(keys.end(), extra.begin(), extra.end());
        for (mesh::LevelKey key : keys)
        {
            const std::vector<mesh::LevelKey> children = mesh::cubechildkeys(n, factor, key);
            ASSERT_EQ(children.size(), static_cast<std::size_t>(factor * factor));
            EXPECT_TRUE(std::is_sorted(children.begin(), children.end()));
            double sum = 0.0;
            for (mesh::LevelKey child : children)
            {
                EXPECT_EQ(mesh::cubeparentkey(n * factor, factor, child), key);
                sum += mesh::cubekeyunitarea(n * factor, child);
            }
            const double parent = mesh::cubekeyunitarea(n, key);
            const double relative = std::fabs(sum - parent) / parent;
            worst = std::max(worst, relative);
            EXPECT_LE(relative, 1e-12) << "N " << n << " key " << key;
        }
    }
    for (std::int64_t n : { 16, 64, 256 })
    {
        double total = 0.0;
        for (mesh::LevelKey key = 0; key < static_cast<mesh::LevelKey>(6 * n * n); key++)
            total += mesh::cubekeyunitarea(n, key);
        EXPECT_NEAR(total, 4.0 * pi, 1e-12) << "N " << n;
    }
    std::cout << "children over parent: worst relative difference " << worst << "\n";
    EXPECT_THROW(mesh::cubechildkeys(mesh::maxlevelfacecells, factor, 0), std::invalid_argument);
    EXPECT_THROW(mesh::cubechildkeys(mesh::maxlevelfacecells / 2, factor, 0), std::invalid_argument);
    EXPECT_NO_THROW(mesh::cubechildkeys(mesh::maxlevelfacecells / factor, factor, 0));
    EXPECT_THROW(mesh::cubeparentkey(6, factor, 0), std::invalid_argument);
}

TEST(MeshWindow, LocateCentreRoundTripAtLevel30)
{
    const std::int64_t n = mesh::maxlevelfacecells;
    std::vector<mesh::LevelKey> keys = cornerkeys(n, 3);
    const auto extra = randomkeys(11, static_cast<std::uint64_t>(6 * n * n), 2000);
    keys.insert(keys.end(), extra.begin(), extra.end());
    for (mesh::LevelKey key : keys)
        EXPECT_EQ(mesh::cubelocatekey(n, mesh::cubekeycentre(n, key)), key) << "key " << key;
    const auto window = mesh::cubewindowgrid(static_cast<int>(n), radius, keys, true);
    for (mesh::CellId c = 0; c < window->cellcount(); c++)
        EXPECT_EQ(window->locate(window->centre(c)), c);
}

TEST(MeshWindow, CompactEqualsStored)
{
    for (std::int64_t n : { std::int64_t{ 4096 }, std::int64_t{ 1 } << 20 })
    {
        // A patch around a cube corner (faces 0, 1 and 4) and one across an edge, plus scattered keys.
        std::vector<mesh::LevelKey> seeds = { mesh::cubekey(n, 0, n - 1, 0), mesh::cubekey(n, 2, n / 2, n - 1) };
        std::vector<mesh::LevelKey> keys = mesh::dilatekeys(n, seeds, 12);
        const auto extra = randomkeys(13, static_cast<std::uint64_t>(6 * n * n), 500);
        keys.insert(keys.end(), extra.begin(), extra.end());
        const auto stored = mesh::cubewindowgrid(static_cast<int>(n), radius, keys, false);
        const auto compact = mesh::cubewindowgrid(static_cast<int>(n), radius, keys, true);
        expectsamegrid(*stored, *compact);
        for (mesh::CellId c = 0; c < compact->cellcount(); c++)
            EXPECT_EQ(compact->localofkey(compact->globalkey(c)), c);
    }
    // Bytes per cell of the vectors each mode keeps (keys, centres, areas, neighbours, reverses; stored
    // adds latitude, longitude and the three per-slot metrics).
    const std::size_t compactbytes = sizeof(mesh::LevelKey) + sizeof(mesh::Vec3) + sizeof(double) +
                                     mesh::slotcount * (sizeof(std::int32_t) + sizeof(std::int8_t));
    const std::size_t storedbytes = compactbytes + 2 * sizeof(double) + 3 * mesh::slotcount * sizeof(double);
    std::cout << "bytes per cell: stored " << storedbytes << ", compact " << compactbytes << "\n";
    EXPECT_LE(compactbytes, 80u);
}

TEST(MeshWindow, RingsWalksAndChart)
{
    const int n = 16;
    const auto base = mesh::makecubegrid(n, radius);
    const std::vector<mesh::LevelKey> seeds = { mesh::cubekey(n, 4, 0, 0), mesh::cubekey(n, 1, n - 1, 7), mesh::cubekey(n, 5, 8, 8) };
    const std::vector<mesh::CellId> seedids = cellids(seeds);
    for (int rings = 0; rings <= 4; rings++)
    {
        const std::vector<int> ring = mesh::ringdistances(*base, seedids, rings);
        std::vector<mesh::LevelKey> expected;
        for (mesh::CellId c = 0; c < base->cellcount(); c++)
            if (ring[c] >= 0)
                expected.push_back(c);
        EXPECT_EQ(mesh::dilatekeys(n, seeds, rings), expected) << "rings " << rings;
    }
    // Straight walks east along a row and north along a column come back after 4 N steps.
    for (int slot : { 3, 1 })
    {
        const mesh::LevelKey start = mesh::cubekey(n, 0, slot == 3 ? 0 : 5, slot == 3 ? 5 : 0);
        const std::vector<mesh::CellId> cells = mesh::straightwalk(*base, static_cast<mesh::CellId>(start), slot, 4 * n);
        ASSERT_EQ(cells.size(), static_cast<std::size_t>(4 * n));
        EXPECT_EQ(cells.back(), start);
        EXPECT_EQ(std::count(cells.begin(), cells.end(), static_cast<mesh::CellId>(start)), 1);
        mesh::WalkStep step{ start, slot };
        for (mesh::CellId cell : cells)
        {
            step = mesh::cubewalkstep(n, step.key, step.slot);
            EXPECT_EQ(step.key, cell);
        }
    }
    // The chart one cell beyond an edge is the edge neighbour; d cells beyond, d - 1 steps further in.
    for (int face = 0; face < mesh::cubefacecount; face++)
        for (int t = 0; t < n; t++)
        {
            EXPECT_EQ(mesh::cubechartkey(n, face, n, t), mesh::cubekeyneighbour(n, mesh::cubekey(n, face, n - 1, t), 3));
            EXPECT_EQ(mesh::cubechartkey(n, face, -1, t), mesh::cubekeyneighbour(n, mesh::cubekey(n, face, 0, t), 7));
            EXPECT_EQ(mesh::cubechartkey(n, face, t, -1), mesh::cubekeyneighbour(n, mesh::cubekey(n, face, t, 0), 1));
            EXPECT_EQ(mesh::cubechartkey(n, face, t, n), mesh::cubekeyneighbour(n, mesh::cubekey(n, face, t, n - 1), 5));
            mesh::WalkStep step{ mesh::cubekey(n, face, n - 1, t), 3 };
            for (int d = 1; d <= 3; d++)
            {
                step = mesh::cubewalkstep(n, step.key, step.slot);
                EXPECT_EQ(mesh::cubechartkey(n, face, n - 1 + d, t), step.key);
            }
        }
    EXPECT_EQ(mesh::cubechartkey(n, 0, n, n), mesh::nokey);
}

TEST(MeshWindow, LatLonWindowEqualsThePlanetWideGrid)
{
    struct Box
    {
        const char* name;
        double west, south, east, north;
    };
    const Box boxes[] = { { "normal", 10.0, 20.0, 30.0, 40.0 }, { "date line", 170.0, -10.0, -170.0, 10.0 }, { "polar", -20.0, 80.0, 20.0, 90.0 } };
    const int worldwidth = 128, scale = 4, halo = 2;
    const int levelcolumns = worldwidth * scale;
    const auto full = mesh::makelatlongrid(levelcolumns, radius);
    for (const Box& b : boxes)
    {
        const landscape::RegionalWindow window = landscape::regionalwindow(worldwidth, scale, halo, b.west, b.south, b.east, b.north);
        // The rows of the planet-wide kW x kH grid (SphereMetric) and the domain's size.
        const landscape::SphereMetric planet(levelcolumns, levelcolumns / 2, radius);
        const int width = window.domainfinewidth(), height = window.domainfineheight();
        mesh::LatLonWindowBox box;
        box.column0 = static_cast<std::int64_t>(window.domaincolumn0) * scale;
        box.columns = window.domaincolumns * scale;
        box.row0 = window.domainrow0 * scale;
        box.rows = window.domainrows * scale;
        const auto grid = mesh::latlonwindowgrid(levelcolumns, radius, box);
        ASSERT_EQ(grid->cellcount(), static_cast<std::size_t>(width) * static_cast<std::size_t>(height)) << b.name;
        EXPECT_EQ(bits(grid->minimumspacing_m()), bits(full->minimumspacing_m())) << b.name;
        EXPECT_EQ(grid->kind(), mesh::GridKind::latlon);
        for (mesh::CellId c = 0; c < grid->cellcount(); c++)
        {
            const int x = static_cast<int>(c % static_cast<mesh::CellId>(width));
            const int y = static_cast<int>(c / static_cast<mesh::CellId>(width));
            const std::size_t planetrow = static_cast<std::size_t>(window.absolutefinerow(y));
            EXPECT_EQ(bits(grid->latitude_rad(c)), bits(planet.rowlatitude_rad[planetrow])) << b.name;
            EXPECT_EQ(bits(grid->area_m2(c)), bits(planet.rowarea_m2[planetrow])) << b.name;
            EXPECT_NEAR(grid->longitude_rad(c), window.finelongitude_rad(x) - pi, 1e-12) << b.name;
            const std::int64_t column = window.absolutefinecolumn(x), row = window.absolutefinerow(y);
            const mesh::LevelKey key = static_cast<mesh::LevelKey>(row) * static_cast<mesh::LevelKey>(levelcolumns) + static_cast<mesh::LevelKey>(column);
            EXPECT_EQ(grid->globalkey(c), key);
            EXPECT_EQ(grid->localofkey(key), c);
            EXPECT_EQ(grid->column(c), column);
            EXPECT_EQ(grid->row(c), row);
            EXPECT_EQ(grid->cellof(static_cast<int>(column), static_cast<int>(row)), c);
            EXPECT_EQ(grid->locate(grid->centre(c)), c);
            const mesh::CellId whole = static_cast<mesh::CellId>(key);
            expectsamevector(grid->centre(c), full->centre(whole));
            EXPECT_EQ(bits(grid->area_m2(c)), bits(full->area_m2(whole)));
            for (int k = 1; k <= mesh::slotcount; k++)
            {
                // The D8 offsets of slot k (1 N, 2 NE, 3 E, ... 8 NW), the neighbour inside the box or none.
                const int dx[8] = { 0, 1, 1, 1, 0, -1, -1, -1 }, dy[8] = { -1, -1, 0, 1, 1, 1, 0, -1 };
                const int nx = x + dx[k - 1], ny = y + dy[k - 1];
                const bool inside = nx >= 0 && nx < width && ny >= 0 && ny < height;
                const std::size_t expected = inside ? static_cast<std::size_t>(ny) * static_cast<std::size_t>(width) + static_cast<std::size_t>(nx) : c;
                const std::int32_t other = grid->neighbour(c, k);
                if (expected == c)
                    EXPECT_EQ(other, -1) << b.name << " cell " << c << " slot " << k;
                else
                {
                    EXPECT_EQ(static_cast<std::size_t>(other), expected) << b.name << " cell " << c << " slot " << k;
                    EXPECT_EQ(grid->globalkey(static_cast<mesh::CellId>(other)), static_cast<mesh::LevelKey>(full->neighbour(whole, k)));
                    EXPECT_EQ(bits(grid->distance_m(c, k)), bits(full->distance_m(whole, k)));
                    EXPECT_EQ(bits(grid->facelength_m(c, k)), bits(full->facelength_m(whole, k)));
                    EXPECT_EQ(bits(grid->bearing_rad(c, k)), bits(full->bearing_rad(whole, k)));
                    EXPECT_EQ(grid->reverse(c, k), full->reverse(whole, k));
                }
                // The distance (also beyond the box's edge): dy north-south, the row's spacing east-west, the
                // diagonal over the two rows' mean spacing, the cell's own row beyond the edge.
                const std::size_t otherrow = static_cast<std::size_t>(window.absolutefinerow(inside ? ny : y));
                const double distance = dx[k - 1] == 0 ? planet.dy_m
                                        : dy[k - 1] == 0 ? planet.rowdx_m[planetrow]
                                                         : std::hypot(0.5 * (planet.rowdx_m[planetrow] + planet.rowdx_m[otherrow]), planet.dy_m);
                EXPECT_EQ(bits(grid->distance_m(c, k)), bits(distance)) << b.name << " cell " << c << " slot " << k;
            }
        }
        // A cell of the level grid outside the box.
        const mesh::LevelKey outside = static_cast<mesh::LevelKey>(box.row0 + box.rows / 2) * static_cast<mesh::LevelKey>(levelcolumns) +
                                       static_cast<mesh::LevelKey>((box.column0 + box.columns + 3) % levelcolumns);
        EXPECT_EQ(grid->localofkey(outside), mesh::nocell) << b.name;
    }
}

namespace
{
struct GhostMeasure
{
    double interior = 0.0;   // largest error more than 3 cells from every face edge
    double edgeband = 0.0;   // within 3 cells of one edge, away from corners
    double cornerband = 0.0; // within 3 cells of two edges
    double edgejump = 0.0;   // largest difference across an edge more than 3 cells from a corner (the two faces' stencils at one point)
    double nearjump = 0.0;   // the same within 3 cells of a corner
    double cornerjump = 0.0; // largest difference between the three faces' stencils at a cube corner
    double weightsum = 0.0;  // largest |sum of weights - 1|
    double constant = 0.0;   // largest |error| of a constant field (exact by applystencil's relative sum)
};

double stencilvalue(std::int64_t n, int face, double a, double b, mesh::GhostRule rule, GhostMeasure& measure)
{
    const mesh::KeyStencil stencil = mesh::cubicstencil(n, face, a, b, rule);
    double sum = 0.0;
    for (double w : stencil.weight)
        sum += w;
    measure.weightsum = std::max(measure.weightsum, std::fabs(sum - 1.0));
    measure.constant = std::max(measure.constant, std::fabs(mesh::applystencil(stencil, [](mesh::LevelKey) { return 7.25; }) - 7.25));
    return mesh::applystencil(stencil, [n](mesh::LevelKey key) { return smoothfield(mesh::cubekeycentre(n, key)); });
}

GhostMeasure measureghosts(std::int64_t n, mesh::GhostRule rule)
{
    GhostMeasure measure;
    const double edge = 0.25 * pi, step = mesh::cubestep(n), band = 3.0 * step;
    std::mt19937_64 random(20261019);
    const auto uniform = [&random](double lo, double hi) { return lo + (hi - lo) * (static_cast<double>(random() >> 11) * 0x1.0p-53); };
    for (int face = 0; face < mesh::cubefacecount; face++)
        for (int sample = 0; sample < 400; sample++)
        {
            // Interior, edge band (one random edge), corner band (one random corner).
            const double ia = uniform(-edge + band, edge - band), ib = uniform(-edge + band, edge - band);
            const double ei = uniform(-edge + band, edge - band), en = uniform(edge - band, edge);
            const int which = sample % 4;
            const double ea = which == 0 ? en : (which == 1 ? -en : ei), eb = which == 2 ? en : (which == 3 ? -en : ei);
            const double ca = (sample % 2 == 0 ? 1.0 : -1.0) * uniform(edge - band, edge);
            const double cb = (sample / 2 % 2 == 0 ? 1.0 : -1.0) * uniform(edge - band, edge);
            const double points[3][2] = { { ia, ib }, { ea, eb }, { ca, cb } };
            double* worst[3] = { &measure.interior, &measure.edgeband, &measure.cornerband };
            for (int m = 0; m < 3; m++)
            {
                const double value = stencilvalue(n, face, points[m][0], points[m][1], rule, measure);
                const double exact = smoothfield(mesh::cubepoint(face, points[m][0], points[m][1]));
                *worst[m] = std::max(*worst[m], std::fabs(value - exact));
            }
            // Across each edge: the point on the edge seen from both faces.
            for (int e = 0; e < 4; e++)
            {
                const double t = uniform(-edge, edge);
                const double a = e == 1 ? edge : (e == 3 ? -edge : t), b = e == 0 ? edge : (e == 2 ? -edge : t);
                const mesh::Vec3 p = mesh::cubepoint(face, a, b);
                const int other = mesh::cubeedgelink(face, e).face;
                const mesh::FacePosition there = mesh::cubepositiononface(other, p);
                const double here = stencilvalue(n, face, a, b, rule, measure);
                const double across = stencilvalue(n, other, there.a, there.b, rule, measure);
                double& jump = std::fabs(t) <= edge - band ? measure.edgejump : measure.nearjump;
                jump = std::max(jump, std::fabs(here - across));
            }
        }
    // The eight cube corners, each from the three faces that meet there.
    for (int face = 0; face < mesh::cubefacecount; face++)
        for (int corner = 0; corner < 4; corner++)
        {
            const double a = corner % 2 == 0 ? edge : -edge, b = corner / 2 == 0 ? edge : -edge;
            const mesh::Vec3 p = mesh::cubepoint(face, a, b);
            const double here = stencilvalue(n, face, a, b, rule, measure);
            for (int other = 0; other < mesh::cubefacecount; other++)
            {
                if (other == face || mesh::dot(p, mesh::cubeface(other).centre) < 0.5)
                    continue;
                const mesh::FacePosition there = mesh::cubepositiononface(other, p);
                measure.cornerjump = std::max(measure.cornerjump, std::fabs(here - stencilvalue(n, other, there.a, there.b, rule, measure)));
            }
        }
    return measure;
}
}

TEST(MeshWindow, CubicStencilGhostsAtTruePositions)
{
    GhostMeasure previous;
    for (std::int64_t n : { 32, 64, 128 })
    {
        const GhostMeasure truth = measureghosts(n, mesh::GhostRule::trueposition);
        const GhostMeasure unfolded = measureghosts(n, mesh::GhostRule::indexunfolded);
        const auto print = [](const char* name, const GhostMeasure& m) {
            std::cout << "  " << name << " interior " << m.interior << ", edge band " << m.edgeband << ", corner band " << m.cornerband
                      << ", edge jump " << m.edgejump << ", edge jump near a corner " << m.nearjump << ", corner jump " << m.cornerjump << "\n";
        };
        std::cout << std::scientific << std::setprecision(2) << "N " << n << " (max |error| of a smooth field; jumps between faces)\n";
        print("true position: ", truth);
        print("index unfolded:", unfolded);
        EXPECT_LT(truth.weightsum, 1e-12);
        EXPECT_LT(unfolded.weightsum, 1e-12);
        EXPECT_EQ(truth.constant, 0.0);
        EXPECT_EQ(unfolded.constant, 0.0);
        EXPECT_EQ(truth.interior, unfolded.interior);  // no ghost reaches the interior
        // Along an edge the true-position ghosts keep the interior's accuracy and the
        // two faces agree to better than it (fourth order); the index-unfolded ones
        // are continuous by symmetry but first-order wrong.
        EXPECT_LT(truth.edgeband, 2.0 * truth.interior);
        EXPECT_LT(truth.edgejump, truth.interior);
        EXPECT_LT(truth.edgeband, 0.1 * unfolded.edgeband);
        // At a corner (beyond two edges: 1-D Catmull-Rom along the third face) the same,
        // converging at third order or better.
        EXPECT_LT(truth.cornerband, 2.0 * truth.interior);
        EXPECT_LT(truth.nearjump, truth.interior);
        EXPECT_LT(truth.cornerjump, truth.interior);
        EXPECT_LT(truth.cornerband, 0.1 * unfolded.cornerband);
        if (previous.cornerjump > 0.0)
        {
            EXPECT_LT(truth.cornerjump, 0.2 * previous.cornerjump);
            EXPECT_LT(truth.cornerband, 0.2 * previous.cornerband);
        }
        previous = truth;
        std::cout << std::defaultfloat;
    }
}
