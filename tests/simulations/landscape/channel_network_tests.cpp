// Tests of landscape/channel_network.hpp (erosion plan Phase E, step E2b), suite Landscape:
// TEST_F(Landscape, ChannelNetwork<What>) on synthetic 16 x 8 worlds with hand-made
// receiver trees, at scales 4 and 5 (an even and an odd centre offset).

#include "landscape_fixture.hpp"

#include "simulations/landscape/channel_network.hpp"
#include "simulations/landscape/regional_grid.hpp"
#include "simulations/landscape/sphere_metric.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <stdexcept>
#include <vector>

namespace
{
using landscape::ChannelParameters;
using landscape::FineChannels;
using landscape::GlobalNetwork;
using landscape::RegionalWindow;

constexpr int worldwidth = 16;
constexpr int worldheight = 8;
constexpr float sealevel = 0.0f;
constexpr float seafilled_m = -100.0f;  // not the sea level: a link into the sea must end at the sea level
constexpr float cellarea_km2 = 100.0f;
constexpr int north = 1, northeast = 2, east = 3, southeast = 4, south = 5, southwest = 6, west = 7, northwest = 8;

std::size_t globalindex(int i, int j) { return static_cast<std::size_t>(j) * worldwidth + static_cast<std::size_t>(i); }

// The receiver tree of the main world. Column 0 is sea. Columns 1 .. 10 drain west:
// rows 2 and 5 are trunks flowing west into the sea; row 0 drains south to row 1,
// which drains south-west onto row 2's trunk (column 1 west into the sea); row 3
// drains north onto row 2; row 4 drains south-east onto row 5 (column 10 south);
// row 6 drains north-west onto row 5 (column 1 west into the sea); row 7 drains north
// to row 6. Columns 11 .. 15 drain east, column 15 across the date line into the sea
// of column 0, except (12, 7), which drains north-east to (13, 6).
int receiverdirection(int i, int j)
{
    if (i >= 11)
        return i == 12 && j == 7 ? northeast : east;
    switch (j)
    {
    case 0: return south;
    case 1: return i >= 2 ? southwest : west;
    case 2: return west;
    case 3: return north;
    case 4: return i <= 9 ? southeast : south;
    case 5: return west;
    case 6: return i >= 2 ? northwest : west;
    default: return north;
    }
}

std::size_t receiverof(std::size_t cell, int direction)
{
    const int i = static_cast<int>(cell % worldwidth);
    const int j = static_cast<int>(cell / worldwidth);
    return globalindex((i + landscape::directionoffsetx(direction) + worldwidth) % worldwidth,
                       j + landscape::directionoffsety(direction));
}

// A world of sea (filled at seafilled_m) with no links, for the hand-built cases.
GlobalNetwork seaworld()
{
    GlobalNetwork network;
    network.width = worldwidth;
    network.height = worldheight;
    const std::size_t cells = static_cast<std::size_t>(worldwidth) * worldheight;
    network.direction.assign(cells, 0);
    network.sea.assign(cells, 1);
    network.drainagearea_km2.assign(cells, 0.0f);
    network.discharge_m3s.assign(cells, 0.0f);
    network.filledelevation_m.assign(cells, seafilled_m);
    network.januaryflow.assign(cells, 0);
    network.julyflow.assign(cells, 0);
    return network;
}

// A land cell (direction 0: a land sink); discharge 0.01 m3/s per km2, legacy flows
// the area and twice the area.
void setland(GlobalNetwork& network, int i, int j, int direction, float filled_m, float area_km2)
{
    const std::size_t cell = globalindex(i, j);
    network.sea[cell] = 0;
    network.direction[cell] = static_cast<uint8_t>(direction);
    network.filledelevation_m[cell] = filled_m;
    network.drainagearea_km2[cell] = area_km2;
    network.discharge_m3s[cell] = 0.01f * area_km2;
    network.januaryflow[cell] = static_cast<int32_t>(area_km2);
    network.julyflow[cell] = 2 * static_cast<int32_t>(area_km2);
}

// The main world: filled elevations stand above the receiver (the sea level for a
// sea receiver) by 5 + 4 ((3 i + 7 j) mod 11) m, except the flat link (5, 2) -> (4, 2);
// areas 100 km2 per land cell accumulated down the tree.
GlobalNetwork syntheticnetwork()
{
    GlobalNetwork network = seaworld();
    const std::size_t cells = static_cast<std::size_t>(worldwidth) * worldheight;
    for (int j = 0; j < worldheight; j++)
        for (int i = 1; i < worldwidth; i++)
        {
            network.sea[globalindex(i, j)] = 0;
            network.direction[globalindex(i, j)] = static_cast<uint8_t>(receiverdirection(i, j));
        }

    std::vector<int> depth(cells, -1);
    std::function<int(std::size_t)> depthof = [&](std::size_t cell) -> int {
        if (network.sea[cell] != 0)
            return 0;
        if (depth[cell] < 0)
            depth[cell] = depthof(receiverof(cell, network.direction[cell])) + 1;
        return depth[cell];
    };
    std::function<float(std::size_t)> filledof = [&](std::size_t cell) -> float {
        if (network.sea[cell] != 0)
            return seafilled_m;
        const std::size_t receiver = receiverof(cell, network.direction[cell]);
        const float base = network.sea[receiver] != 0 ? sealevel : filledof(receiver);
        const int i = static_cast<int>(cell % worldwidth);
        const int j = static_cast<int>(cell / worldwidth);
        const float drop = i == 5 && j == 2 ? 0.0f : static_cast<float>(5 + 4 * ((3 * i + 7 * j) % 11));
        return base + drop;
    };
    std::vector<std::size_t> order;
    for (std::size_t cell = 0; cell < cells; cell++)
        if (network.sea[cell] == 0)
        {
            depthof(cell);
            network.filledelevation_m[cell] = filledof(cell);
            network.drainagearea_km2[cell] = cellarea_km2;
            order.push_back(cell);
        }
    std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return depth[a] > depth[b]; });
    for (std::size_t cell : order)
    {
        const std::size_t receiver = receiverof(cell, network.direction[cell]);
        if (network.sea[receiver] == 0)
            network.drainagearea_km2[receiver] += network.drainagearea_km2[cell];
    }
    for (std::size_t cell : order)
    {
        const float area = network.drainagearea_km2[cell];
        network.discharge_m3s[cell] = 0.01f * area;
        network.januaryflow[cell] = static_cast<int32_t>(area);
        network.julyflow[cell] = 2 * static_cast<int32_t>(area);
    }
    return network;
}

// A window whose kept box is its domain (no halo), built by hand: the accessors of
// regional_grid.cpp are not used here.
RegionalWindow domainwindow(int scale, int column0, int columns, int row0, int rows)
{
    RegionalWindow window;
    window.worldwidth = worldwidth;
    window.worldheight = worldheight;
    window.scale = scale;
    window.halo = 0;
    window.column0 = window.domaincolumn0 = column0;
    window.columns = window.domaincolumns = columns;
    window.row0 = window.domainrow0 = row0;
    window.rows = window.domainrows = rows;
    return window;
}

RegionalWindow worldwindow(int scale) { return domainwindow(scale, 0, worldwidth, 0, worldheight); }

ChannelParameters parameters(double meanderamplitude, std::uint64_t seed, double pinarea_km2 = 0.0)
{
    ChannelParameters result;
    result.pinarea_km2 = pinarea_km2;
    result.meanderamplitude = meanderamplitude;
    result.seed = seed;
    return result;
}

// Absolute fine coordinates of a domain fine cell, and back (false off the domain).
std::int64_t absolutecolumn(const RegionalWindow& window, int x)
{
    const std::int64_t planet = static_cast<std::int64_t>(window.scale) * worldwidth;
    return (static_cast<std::int64_t>(window.domaincolumn0) * window.scale + x) % planet;
}

std::int64_t absoluterow(const RegionalWindow& window, int y)
{
    return static_cast<std::int64_t>(window.domainrow0) * window.scale + y;
}

bool domaincell(const RegionalWindow& window, std::int64_t column, std::int64_t row, std::size_t& index)
{
    const std::int64_t planet = static_cast<std::int64_t>(window.scale) * worldwidth;
    const std::int64_t y = row - static_cast<std::int64_t>(window.domainrow0) * window.scale;
    std::int64_t x = ((column - static_cast<std::int64_t>(window.domaincolumn0) * window.scale) % planet + planet) % planet;
    if (y < 0 || y >= window.domainfineheight() || x >= window.domainfinewidth())
        return false;
    index = static_cast<std::size_t>(y) * static_cast<std::size_t>(window.domainfinewidth()) + static_cast<std::size_t>(x);
    return true;
}

// The centre fine cell of global (i, j): (k i + k / 2, k j + k / 2).
std::size_t centrecell(const RegionalWindow& window, int i, int j)
{
    std::size_t index = 0;
    EXPECT_TRUE(domaincell(window, static_cast<std::int64_t>(window.scale) * i + window.scale / 2,
                           static_cast<std::int64_t>(window.scale) * j + window.scale / 2, index));
    return index;
}

// The global cell holding a domain fine cell.
std::size_t globalof(const RegionalWindow& window, std::size_t index)
{
    const int x = static_cast<int>(index % static_cast<std::size_t>(window.domainfinewidth()));
    const int y = static_cast<int>(index / static_cast<std::size_t>(window.domainfinewidth()));
    return globalindex(static_cast<int>(absolutecolumn(window, x) / window.scale),
                       static_cast<int>(absoluterow(window, y) / window.scale));
}

bool iscentre(const RegionalWindow& window, std::size_t index)
{
    const int x = static_cast<int>(index % static_cast<std::size_t>(window.domainfinewidth()));
    const int y = static_cast<int>(index / static_cast<std::size_t>(window.domainfinewidth()));
    return absolutecolumn(window, x) % window.scale == window.scale / 2 &&
           absoluterow(window, y) % window.scale == window.scale / 2;
}

// Where a chain may end: any fine cell of a sea cell's block, or a land sink's centre.
bool isend(const GlobalNetwork& network, const RegionalWindow& window, std::size_t index)
{
    const std::size_t cell = globalof(window, index);
    if (network.sea[cell] != 0)
        return true;
    return network.direction[cell] == 0 && iscentre(window, index);
}

// The domain fine cell a direction leads to (false off the domain).
bool nextcell(const RegionalWindow& window, std::size_t index, int direction, std::size_t& next)
{
    const int x = static_cast<int>(index % static_cast<std::size_t>(window.domainfinewidth()));
    const int y = static_cast<int>(index / static_cast<std::size_t>(window.domainfinewidth()));
    return domaincell(window, absolutecolumn(window, x) + landscape::directionoffsetx(direction),
                      absoluterow(window, y) + landscape::directionoffsety(direction), next);
}

// The global receiver of a land cell with a direction; the cell itself otherwise.
std::size_t globalreceiver(const GlobalNetwork& network, std::size_t cell)
{
    if (network.sea[cell] != 0 || network.direction[cell] == 0)
        return cell;
    return receiverof(cell, network.direction[cell]);
}

// The global outlet of a cell: the first sea cell or land sink down the tree.
std::size_t globaloutlet(const GlobalNetwork& network, std::size_t cell)
{
    while (globalreceiver(network, cell) != cell)
        cell = globalreceiver(network, cell);
    return cell;
}

// Links adjacent on the global tree: the same receiver, or one's receiver is the
// other's upstream cell.
bool adjacentlinks(const GlobalNetwork& network, std::size_t a, std::size_t b)
{
    const std::size_t ra = globalreceiver(network, a);
    const std::size_t rb = globalreceiver(network, b);
    return ra == rb || ra == b || rb == a;
}

// Walks the directions from every pinned global cell's centre fine cell (the
// window must be the whole world): counts the walks, the changes of owner
// (globalcell) to a link not adjacent on the global tree, and the walks that do not
// end in the starting cell's global outlet (a sea cell's block or a land sink's
// centre) or loop.
struct WalkReport
{
    int walks = 0;
    int nonadjacent = 0;
    int wrongoutlet = 0;
};

WalkReport walkpinned(const GlobalNetwork& network, const RegionalWindow& window, const FineChannels& channels,
                      double pinarea_km2 = 0.0)
{
    WalkReport report;
    const int limit = static_cast<int>(channels.channel.size());
    for (int j = 0; j < worldheight; j++)
        for (int i = 0; i < worldwidth; i++)
        {
            const std::size_t start = globalindex(i, j);
            if (network.sea[start] != 0 || network.direction[start] == 0 ||
                !(network.drainagearea_km2[start] >= pinarea_km2))
                continue;
            report.walks++;
            std::size_t walk = centrecell(window, i, j);
            int steps = 0;
            bool ok = channels.globalcell[walk] == static_cast<int32_t>(start);
            while (ok && channels.direction[walk] != 0 && steps <= limit)
            {
                std::size_t next = 0;
                ok = nextcell(window, walk, channels.direction[walk], next) && channels.channel[next] == 1;
                if (!ok)
                    break;
                const int32_t from = channels.globalcell[walk];
                const int32_t to = channels.globalcell[next];
                if (from != to &&
                    !adjacentlinks(network, static_cast<std::size_t>(from), static_cast<std::size_t>(to)))
                    report.nonadjacent++;
                walk = next;
                steps++;
            }
            const std::size_t outlet = globaloutlet(network, start);
            if (!ok || steps > limit || globalof(window, walk) != outlet || !isend(network, window, walk))
                report.wrongoutlet++;
        }
    return report;
}

// Steps from a channel cell to a non-centre cell owned by another link.
int contestedsteps(const RegionalWindow& window, const FineChannels& channels)
{
    int count = 0;
    for (std::size_t cell = 0; cell < channels.channel.size(); cell++)
    {
        std::size_t next = 0;
        if (channels.channel[cell] == 0 || channels.direction[cell] == 0 ||
            !nextcell(window, cell, channels.direction[cell], next))
            continue;
        if (channels.globalcell[next] != channels.globalcell[cell] && !iscentre(window, next))
            count++;
    }
    return count;
}

struct Case
{
    int scale;
    double meanderamplitude;
    std::uint64_t seed;
};

const Case cases[] = { { 4, 0.3, 1 }, { 5, 0.3, 1 }, { 4, 0.45, 2 }, { 5, 0.45, 2 }, { 8, 0.3, 3 },
                       { 4, 0.0, 1 }, { 5, 0.0, 1 } };

int channelcount(const FineChannels& channels)
{
    return static_cast<int>(std::count(channels.channel.begin(), channels.channel.end(), uint8_t{ 1 }));
}
}

TEST_F(Landscape, ChannelNetworkDirectionsLeadToTheSea)
{
    const GlobalNetwork network = syntheticnetwork();
    int contested = 0;
    for (const Case& c : cases)
    {
        SCOPED_TRACE(testing::Message() << "k " << c.scale << ", meander " << c.meanderamplitude << ", seed " << c.seed);
        const RegionalWindow window = worldwindow(c.scale);
        const FineChannels channels = rasterisechannels(window, network, parameters(c.meanderamplitude, c.seed), sealevel);
        const int count = channelcount(channels);
        EXPECT_GT(count, 0);
        for (std::size_t cell = 0; cell < channels.channel.size(); cell++)
        {
            if (channels.channel[cell] == 0)
            {
                EXPECT_EQ(channels.direction[cell], 0);
                EXPECT_EQ(channels.globalcell[cell], -1);
                continue;
            }
            // One step: a channel neighbour, or the end of a link in a sea cell's block.
            if (channels.direction[cell] == 0)
            {
                EXPECT_TRUE(network.sea[globalof(window, cell)] != 0) << "fine cell " << cell;
                EXPECT_EQ(channels.bed_m[cell], sealevel) << "fine cell " << cell;
                continue;
            }
            EXPECT_EQ(network.sea[globalof(window, cell)], 0) << "fine cell " << cell;
            std::size_t next = 0;
            ASSERT_TRUE(nextcell(window, cell, channels.direction[cell], next)) << "fine cell " << cell;
            EXPECT_EQ(channels.channel[next], 1) << "fine cell " << cell;
            // Following the directions reaches the sea within count steps (no loop).
            std::size_t walk = cell;
            int steps = 0;
            while (channels.direction[walk] != 0 && steps <= count)
            {
                ASSERT_TRUE(nextcell(window, walk, channels.direction[walk], walk));
                ASSERT_EQ(channels.channel[walk], 1);
                steps++;
            }
            EXPECT_LE(steps, count) << "loop from fine cell " << cell;
            EXPECT_TRUE(isend(network, window, walk)) << "from fine cell " << cell;
        }
        // Every owner change along the way is to a link adjacent on the global tree,
        // and every pinned cell's chain ends in its own global outlet.
        const WalkReport report = walkpinned(network, window, channels);
        EXPECT_EQ(report.walks, 120);
        EXPECT_EQ(report.nonadjacent, 0);
        EXPECT_EQ(report.wrongoutlet, 0);
        contested += contestedsteps(window, channels);
    }
    // The ownership rule is exercised: some chains run into another link's cell
    // that is not a centre.
    EXPECT_GT(contested, 0);
}

TEST_F(Landscape, ChannelNetworkBedNeverRises)
{
    const GlobalNetwork network = syntheticnetwork();
    for (const Case& c : cases)
    {
        SCOPED_TRACE(testing::Message() << "k " << c.scale << ", meander " << c.meanderamplitude << ", seed " << c.seed);
        const RegionalWindow window = worldwindow(c.scale);
        const FineChannels channels = rasterisechannels(window, network, parameters(c.meanderamplitude, c.seed), sealevel);
        for (std::size_t cell = 0; cell < channels.channel.size(); cell++)
        {
            if (channels.channel[cell] == 0 || channels.direction[cell] == 0)
                continue;
            std::size_t next = 0;
            ASSERT_TRUE(nextcell(window, cell, channels.direction[cell], next));
            EXPECT_LE(channels.bed_m[next], channels.bed_m[cell]) << "fine cell " << cell;
            EXPECT_GE(channels.bed_m[cell], sealevel);
        }
    }
}

TEST_F(Landscape, ChannelNetworkCentresArePinned)
{
    const GlobalNetwork network = syntheticnetwork();
    for (const Case& c : cases)
    {
        SCOPED_TRACE(testing::Message() << "k " << c.scale << ", meander " << c.meanderamplitude << ", seed " << c.seed);
        const RegionalWindow window = worldwindow(c.scale);
        const FineChannels channels = rasterisechannels(window, network, parameters(c.meanderamplitude, c.seed), sealevel);
        for (int j = 0; j < worldheight; j++)
            for (int i = 0; i < worldwidth; i++)
            {
                const std::size_t cell = globalindex(i, j);
                const std::size_t centre = centrecell(window, i, j);
                if (network.sea[cell] != 0)
                {
                    // A link into the sea ends at its first fine cell in the sea cell's
                    // block, on the block's edge: the centre is not reached.
                    EXPECT_EQ(channels.channel[centre], 0) << "sea " << i << ", " << j;
                    continue;
                }
                EXPECT_EQ(channels.channel[centre], 1) << "land " << i << ", " << j;
                EXPECT_EQ(channels.bed_m[centre], network.filledelevation_m[cell]) << "land " << i << ", " << j;
                EXPECT_EQ(channels.globalcell[centre], static_cast<int32_t>(cell));
                EXPECT_NE(channels.direction[centre], 0);
                EXPECT_EQ(channels.discharge_m3s[centre], network.discharge_m3s[cell]);
                EXPECT_EQ(channels.januaryflow[centre], network.januaryflow[cell]);
                EXPECT_EQ(channels.julyflow[centre], network.julyflow[cell]);
            }
    }
}

TEST_F(Landscape, ChannelNetworkSmallLinksAreAbsent)
{
    const GlobalNetwork network = syntheticnetwork();
    constexpr double pinarea_km2 = 250.0;  // headwaters (100 km2) and cells with one headwater donor (200 km2)
    for (int scale : { 4, 5 })
    {
        const RegionalWindow window = worldwindow(scale);
        const FineChannels channels = rasterisechannels(window, network, parameters(0.3, 1, pinarea_km2), sealevel);
        int pinned = 0;
        int unpinned = 0;
        for (int j = 0; j < worldheight; j++)
            for (int i = 1; i < worldwidth; i++)
            {
                const std::size_t cell = globalindex(i, j);
                const std::size_t centre = centrecell(window, i, j);
                const bool expected = network.drainagearea_km2[cell] >= pinarea_km2;
                (expected ? pinned : unpinned)++;
                EXPECT_EQ(channels.channel[centre], expected ? 1 : 0) << i << ", " << j;
                const bool drawn = std::find(channels.globalcell.begin(), channels.globalcell.end(),
                                             static_cast<int32_t>(cell)) != channels.globalcell.end();
                EXPECT_EQ(drawn, expected) << i << ", " << j;
            }
        // Of the 120 land cells, 41 reach 250 km2: rows 2 (columns 1 .. 9) and 5
        // (columns 1 .. 10), columns 13 and 14 (rows 0 .. 6) and column 15.
        EXPECT_EQ(pinned, 41);
        EXPECT_EQ(unpinned, 79);
        const WalkReport report = walkpinned(network, window, channels, pinarea_km2);
        EXPECT_EQ(report.walks, 41);
        EXPECT_EQ(report.nonadjacent, 0);
        EXPECT_EQ(report.wrongoutlet, 0);
    }
}

TEST_F(Landscape, ChannelNetworkStraightLinkWithoutMeander)
{
    const GlobalNetwork network = syntheticnetwork();
    for (int k : { 4, 5 })
    {
        SCOPED_TRACE(testing::Message() << "k " << k);
        const RegionalWindow window = worldwindow(k);
        const FineChannels channels = rasterisechannels(window, network, parameters(0.0, 1), sealevel);
        // The link (6, 2) -> (5, 2) along row 2: fine row 2 k + k / 2, from fine column
        // 6 k + k / 2 west to 5 k + k / 2; its bed falls by 45 m (the drop 5 + 4 (32 mod 11))
        // linearly over k fine cells.
        const std::size_t upstream = globalindex(6, 2);
        const float top = network.filledelevation_m[upstream];
        const float bottom = network.filledelevation_m[globalindex(5, 2)];
        EXPECT_EQ(top - bottom, 45.0f);
        const std::int64_t row = 2 * k + k / 2;
        int owned = 0;
        for (std::size_t cell = 0; cell < channels.globalcell.size(); cell++)
            owned += channels.globalcell[cell] == static_cast<int32_t>(upstream) ? 1 : 0;
        EXPECT_EQ(owned, k);  // the start centre and the k - 1 cells before the receiver's centre
        for (int m = 0; m < k; m++)
        {
            std::size_t cell = 0;
            ASSERT_TRUE(domaincell(window, 6 * k + k / 2 - m, row, cell));
            EXPECT_EQ(channels.channel[cell], 1) << "step " << m;
            EXPECT_EQ(channels.globalcell[cell], static_cast<int32_t>(upstream)) << "step " << m;
            EXPECT_EQ(channels.direction[cell], west) << "step " << m;
            EXPECT_FLOAT_EQ(channels.bed_m[cell], top - 45.0f * static_cast<float>(m) / static_cast<float>(k))
                << "step " << m;
        }
        std::size_t receiver = 0;
        ASSERT_TRUE(domaincell(window, 5 * k + k / 2, row, receiver));
        EXPECT_EQ(channels.globalcell[receiver], static_cast<int32_t>(globalindex(5, 2)));
        EXPECT_EQ(channels.bed_m[receiver], bottom);
        // The flat link (5, 2) -> (4, 2) keeps its bed level.
        for (int m = 1; m < k; m++)
        {
            std::size_t cell = 0;
            ASSERT_TRUE(domaincell(window, 5 * k + k / 2 - m, row, cell));
            EXPECT_EQ(channels.bed_m[cell], bottom) << "step " << m;
        }
        // The link (1, 2) -> (0, 2) into the sea: from fine column k + k / 2 west, its
        // bed falls linearly to the sea level at fine column k - 1, its first cell in
        // the sea cell's block, over k / 2 + 1 steps; the rest of the block's row holds
        // nothing of it (fine column 0 is the end of (15, 2) -> (0, 2), which enters the
        // block from the west across the date line).
        const std::size_t coast = globalindex(1, 2);
        const float coasttop = network.filledelevation_m[coast];
        const int steps = k / 2 + 1;
        for (int m = 0; m <= steps; m++)
        {
            std::size_t cell = 0;
            ASSERT_TRUE(domaincell(window, k + k / 2 - m, row, cell));
            EXPECT_EQ(channels.globalcell[cell], static_cast<int32_t>(coast)) << "step " << m;
            EXPECT_EQ(channels.direction[cell], m == steps ? 0 : west) << "step " << m;
            EXPECT_FLOAT_EQ(channels.bed_m[cell],
                            coasttop - coasttop * static_cast<float>(m) / static_cast<float>(steps))
                << "step " << m;
        }
        for (int column = 1; column < k - 1; column++)
        {
            std::size_t cell = 0;
            ASSERT_TRUE(domaincell(window, column, row, cell));
            EXPECT_EQ(channels.channel[cell], 0) << "column " << column;
        }
        std::size_t westend = 0;
        ASSERT_TRUE(domaincell(window, 0, row, westend));
        EXPECT_EQ(channels.globalcell[westend], static_cast<int32_t>(globalindex(15, 2)));
        EXPECT_EQ(channels.direction[westend], 0);
    }
}

TEST_F(Landscape, ChannelNetworkLinkAcrossTheDateLine)
{
    const GlobalNetwork network = syntheticnetwork();
    for (int k : { 4, 5 })
    {
        SCOPED_TRACE(testing::Message() << "k " << k);
        const std::int64_t planet = static_cast<std::int64_t>(k) * worldwidth;
        const RegionalWindow across = domainwindow(k, 14, 4, 2, 3);  // columns 14, 15, 0, 1
        // Without a meander: (15, 3) -> (0, 3) is the straight fine row 3 k + k / 2 from
        // column 15 k + k / 2 east across the date line to fine column 0, its first cell
        // in the sea cell's block, where it ends at the sea level.
        for (const RegionalWindow& window : { worldwindow(k), across })
        {
            const FineChannels channels = rasterisechannels(window, network, parameters(0.0, 1), sealevel);
            const std::int64_t row = 3 * k + k / 2;
            for (std::int64_t column = 15 * k + k / 2; column < planet; column++)
            {
                std::size_t cell = 0;
                ASSERT_TRUE(domaincell(window, column, row, cell));
                EXPECT_EQ(channels.globalcell[cell], static_cast<int32_t>(globalindex(15, 3))) << column;
                EXPECT_EQ(channels.direction[cell], east) << column;
            }
            std::size_t mouth = 0;
            ASSERT_TRUE(domaincell(window, 0, row, mouth));
            EXPECT_EQ(channels.globalcell[mouth], static_cast<int32_t>(globalindex(15, 3)));
            EXPECT_EQ(channels.direction[mouth], 0);
            EXPECT_EQ(channels.bed_m[mouth], sealevel);
            EXPECT_EQ(channels.channel[centrecell(window, 0, 3)], 0);
        }
        // With a meander: following the directions from (15, 3)'s centre crosses the date
        // line, one 8-neighbour step at a time, into the sea block of column 0.
        for (const RegionalWindow& window : { worldwindow(k), across })
        {
            const FineChannels channels = rasterisechannels(window, network, parameters(0.45, 2), sealevel);
            std::size_t walk = centrecell(window, 15, 3);
            bool crossed = false;
            int steps = 0;
            while (channels.direction[walk] != 0 && steps < 4 * k)
            {
                const int x = static_cast<int>(walk % static_cast<std::size_t>(window.domainfinewidth()));
                const std::int64_t before = absolutecolumn(window, x);
                ASSERT_TRUE(nextcell(window, walk, channels.direction[walk], walk));
                ASSERT_EQ(channels.channel[walk], 1);
                const int after = static_cast<int>(walk % static_cast<std::size_t>(window.domainfinewidth()));
                crossed = crossed || (before == planet - 1 && absolutecolumn(window, after) == 0);
                steps++;
            }
            EXPECT_TRUE(crossed);
            EXPECT_EQ(globalof(window, walk), globalindex(0, 3));
            const int x = static_cast<int>(walk % static_cast<std::size_t>(window.domainfinewidth()));
            EXPECT_EQ(absolutecolumn(window, x), 0);  // the block's first column
        }
    }
}

TEST_F(Landscape, ChannelNetworkOverlappingWindowsAgree)
{
    // Every fine cell's values come from the links whose paths cross it, and every
    // window draws every link that reaches its domain: windows agree on every shared
    // fine cell, up to their edges (no margin is needed).
    const GlobalNetwork network = syntheticnetwork();
    for (const Case& c : cases)
    {
        SCOPED_TRACE(testing::Message() << "k " << c.scale << ", meander " << c.meanderamplitude << ", seed " << c.seed);
        const ChannelParameters settings = parameters(c.meanderamplitude, c.seed);
        const RegionalWindow windows[] = {
            domainwindow(c.scale, 5, worldwidth, 0, worldheight),  // the whole world from column 5
            domainwindow(c.scale, 13, 8, 1, 5),                    // across the date line
            domainwindow(c.scale, 2, 10, 0, 8),
            domainwindow(c.scale, 9, 5, 3, 4),
            domainwindow(c.scale, 0, 3, 6, 2),
        };
        std::vector<FineChannels> results;
        for (const RegionalWindow& window : windows)
            results.push_back(rasterisechannels(window, network, settings, sealevel));
        for (std::size_t p = 0; p < results.size(); p++)
            for (std::size_t q = 0; q < results.size(); q++)
            {
                if (p == q)
                    continue;
                const RegionalWindow& a = windows[p];
                const RegionalWindow& b = windows[q];
                int shared = 0;
                for (int y = 0; y < a.domainfineheight(); y++)
                    for (int x = 0; x < a.domainfinewidth(); x++)
                    {
                        std::size_t other = 0;
                        if (!domaincell(b, absolutecolumn(a, x), absoluterow(a, y), other))
                            continue;
                        const std::size_t cell =
                            static_cast<std::size_t>(y) * static_cast<std::size_t>(a.domainfinewidth()) + x;
                        shared++;
                        ASSERT_EQ(results[p].channel[cell], results[q].channel[other]) << p << " " << q << " " << cell;
                        ASSERT_EQ(results[p].direction[cell], results[q].direction[other]) << p << " " << q << " " << cell;
                        ASSERT_EQ(results[p].bed_m[cell], results[q].bed_m[other]) << p << " " << q << " " << cell;
                        ASSERT_EQ(results[p].discharge_m3s[cell], results[q].discharge_m3s[other]);
                        ASSERT_EQ(results[p].januaryflow[cell], results[q].januaryflow[other]);
                        ASSERT_EQ(results[p].julyflow[cell], results[q].julyflow[other]);
                        ASSERT_EQ(results[p].globalcell[cell], results[q].globalcell[other]);
                    }
                // Windows 1 (columns 13 .. 4, rows 1 .. 5) and 2 (columns 2 .. 11) share 3 x 5 global cells.
                if ((p == 1 && q == 2) || (p == 2 && q == 1))
                    EXPECT_EQ(shared, 15 * c.scale * c.scale);
                if (q == 0)
                    EXPECT_EQ(shared, a.domainfinewidth() * a.domainfineheight());
            }
    }
}

TEST_F(Landscape, ChannelNetworkRepeatable)
{
    const GlobalNetwork network = syntheticnetwork();
    const RegionalWindow window = domainwindow(5, 13, 8, 1, 5);
    const ChannelParameters settings = parameters(0.3, 7);
    const FineChannels first = rasterisechannels(window, network, settings, sealevel);
    const FineChannels second = rasterisechannels(window, network, settings, sealevel);
    EXPECT_EQ(first.channel, second.channel);
    EXPECT_EQ(first.direction, second.direction);
    EXPECT_EQ(first.bed_m, second.bed_m);
    EXPECT_EQ(first.discharge_m3s, second.discharge_m3s);
    EXPECT_EQ(first.januaryflow, second.januaryflow);
    EXPECT_EQ(first.julyflow, second.julyflow);
    EXPECT_EQ(first.globalcell, second.globalcell);
    // Another seed moves the meanders.
    const FineChannels other = rasterisechannels(window, network, parameters(0.3, 8), sealevel);
    EXPECT_NE(first.globalcell, other.globalcell);
}

TEST_F(Landscape, ChannelNetworkParallelDiagonalsStayApart)
{
    // Two parallel south-east links in different basins, k / sqrt(2) fine cells apart:
    // (5, 3) -> (6, 4), which drains west to the sea at (4, 4), and (6, 3) -> (7, 4),
    // which drains east to the sea at (9, 4). (6, 4) -> (5, 4) passes k / sqrt(2) from
    // the second link's middle as well. Over many seeds the meanders bulge towards
    // each other; no chain may be captured into the other basin. (The cap of E2b's
    // first version, (k - 1) / 2, let such curves come to about 0.7 fine cells, but
    // only with every meander weight at its extreme: of 20000 seeds at each of k = 4
    // and 8, and 256 at each of k = 3, 4, 5, 6, 8, 12 and 16, none shared a cell. The
    // cap 0.35 (k - 2) makes sharing impossible; this test guards it at k = 3, 4, 5,
    // 6, 8 and 16, amplitudes 0.3 and 1.0, seeds 1 .. 64.)
    GlobalNetwork network = seaworld();
    setland(network, 5, 3, southeast, 40.0f, 100.0f);
    setland(network, 6, 4, west, 20.0f, 200.0f);
    setland(network, 5, 4, west, 10.0f, 300.0f);
    setland(network, 6, 3, southeast, 40.0f, 100.0f);
    setland(network, 7, 4, east, 20.0f, 200.0f);
    setland(network, 8, 4, east, 10.0f, 300.0f);
    setland(network, 7, 3, 0, 40.0f, 100.0f);  // a land sink beside the second diagonal: not a sea corner
    for (int k : { 3, 4, 5, 6, 8, 16 })
        for (double meander : { 0.3, 1.0 })
            for (std::uint64_t seed = 1; seed <= 64; seed++)
            {
                SCOPED_TRACE(testing::Message() << "k " << k << ", meander " << meander << ", seed " << seed);
                const RegionalWindow window = worldwindow(k);
                const FineChannels channels = rasterisechannels(window, network, parameters(meander, seed), sealevel);
                const WalkReport report = walkpinned(network, window, channels);
                ASSERT_EQ(report.walks, 6);
                ASSERT_EQ(report.nonadjacent, 0);
                ASSERT_EQ(report.wrongoutlet, 0);
            }
}

TEST_F(Landscape, ChannelNetworkCrossingDiagonalsMeet)
{
    // An X across the 2 x 2 block (6 .. 7, 3 .. 4): (6, 3) -> (7, 4), then east to the
    // sea at (9, 4); (7, 3) -> (6, 4), then west to the sea at (4, 4). The beds where
    // the two diagonals cross (their middles) are 60 m and 25 m. Without a meander the
    // two staircases share fine cell (28, 16) at k = 4, which the lower link takes, so
    // the chain of (6, 3) is captured into the west basin (the accepted exception);
    // at k = 5 they cross between diagonal neighbours without sharing a cell and
    // each reaches its own outlet.
    GlobalNetwork network = seaworld();
    setland(network, 6, 3, southeast, 100.0f, 100.0f);
    setland(network, 7, 4, east, 20.0f, 200.0f);
    setland(network, 8, 4, east, 10.0f, 300.0f);
    setland(network, 7, 3, southwest, 30.0f, 100.0f);
    setland(network, 6, 4, west, 20.0f, 200.0f);
    setland(network, 5, 4, west, 10.0f, 300.0f);
    {
        const RegionalWindow window = worldwindow(4);
        const FineChannels channels = rasterisechannels(window, network, parameters(0.0, 1), sealevel);
        std::size_t crossing = 0;
        ASSERT_TRUE(domaincell(window, 28, 16, crossing));
        EXPECT_EQ(channels.globalcell[crossing], static_cast<int32_t>(globalindex(7, 3)));
        EXPECT_FLOAT_EQ(channels.bed_m[crossing], 25.0f);
        const WalkReport report = walkpinned(network, window, channels);
        EXPECT_EQ(report.walks, 6);
        EXPECT_EQ(report.nonadjacent, 1);   // (6, 3) runs into (7, 3)'s link
        EXPECT_EQ(report.wrongoutlet, 1);   // and ends in the west basin
    }
    {
        const RegionalWindow window = worldwindow(5);
        const FineChannels channels = rasterisechannels(window, network, parameters(0.0, 1), sealevel);
        const WalkReport report = walkpinned(network, window, channels);
        EXPECT_EQ(report.walks, 6);
        EXPECT_EQ(report.nonadjacent, 0);
        EXPECT_EQ(report.wrongoutlet, 0);
    }
    // With meanders the crossing may be shared or not; the beds never rise and every
    // chain still ends in the sea.
    for (int k : { 4, 5 })
        for (std::uint64_t seed = 1; seed <= 16; seed++)
        {
            SCOPED_TRACE(testing::Message() << "k " << k << ", seed " << seed);
            const RegionalWindow window = worldwindow(k);
            const FineChannels channels = rasterisechannels(window, network, parameters(0.3, seed), sealevel);
            const WalkReport report = walkpinned(network, window, channels);
            EXPECT_LE(report.nonadjacent, 1);
            EXPECT_EQ(report.wrongoutlet, report.nonadjacent);
            for (std::size_t cell = 0; cell < channels.channel.size(); cell++)
            {
                std::size_t next = 0;
                if (channels.channel[cell] != 0 && channels.direction[cell] != 0)
                {
                    ASSERT_TRUE(nextcell(window, cell, channels.direction[cell], next));
                    EXPECT_LE(channels.bed_m[next], channels.bed_m[cell]);
                }
            }
        }
}

TEST_F(Landscape, ChannelNetworkGentleTributaryKeepsDescending)
{
    // A gentle tributary and a steep, larger link meet at X = (6, 4), which drains
    // west to the sea: A = (7, 4) -> X, 1000 m above X and 5000 km2; B = (7, 3) -> X
    // (south-west), 1 m above X and 100 km2. Near X their meanders share fine cells,
    // where A's bed is far above B's. Under ownership by the larger area, B's chain
    // would step from its own cell up onto A's; under the lowest bed it does not.
    // A's path alone is drawn by making B a sink, so the cells A would hold are known.
    GlobalNetwork network = seaworld();
    setland(network, 6, 4, west, 20.0f, 5200.0f);
    setland(network, 5, 4, west, 10.0f, 5300.0f);
    setland(network, 7, 4, west, 1020.0f, 5000.0f);
    setland(network, 7, 3, southwest, 21.0f, 100.0f);
    setland(network, 6, 3, 0, 30.0f, 100.0f);  // a land sink beside B's diagonal: not a sea corner
    GlobalNetwork alone = network;
    alone.direction[globalindex(7, 3)] = 0;
    int wouldrise = 0;
    for (int k : { 4, 5, 8 })
        for (std::uint64_t seed = 1; seed <= 64; seed++)
        {
            SCOPED_TRACE(testing::Message() << "k " << k << ", seed " << seed);
            const RegionalWindow window = worldwindow(k);
            const FineChannels both = rasterisechannels(window, network, parameters(0.3, seed), sealevel);
            const FineChannels steep = rasterisechannels(window, alone, parameters(0.3, seed), sealevel);
            const int32_t a = static_cast<int32_t>(globalindex(7, 4));
            const int32_t b = static_cast<int32_t>(globalindex(7, 3));
            for (std::size_t cell = 0; cell < both.channel.size(); cell++)
            {
                if (both.channel[cell] == 0 || both.direction[cell] == 0)
                    continue;
                std::size_t next = 0;
                ASSERT_TRUE(nextcell(window, cell, both.direction[cell], next));
                EXPECT_LE(both.bed_m[next], both.bed_m[cell]);
                // B steps into a cell that A's path also crosses, higher than B's cell.
                if (both.globalcell[cell] == b && steep.globalcell[next] == a && steep.bed_m[next] > both.bed_m[cell])
                    wouldrise++;
            }
            const WalkReport report = walkpinned(network, window, both);
            EXPECT_EQ(report.nonadjacent, 0);
            EXPECT_EQ(report.wrongoutlet, 0);
        }
    EXPECT_GT(wouldrise, 0);
}

TEST_F(Landscape, ChannelNetworkLandSinkEndsItsChain)
{
    // (10, 2) is a land sink 50 m above datum; (11, 2) drains west into it from 80 m.
    GlobalNetwork network = seaworld();
    setland(network, 10, 2, 0, 50.0f, 200.0f);
    setland(network, 11, 2, west, 80.0f, 100.0f);
    for (int k : { 4, 5 })
    {
        SCOPED_TRACE(testing::Message() << "k " << k);
        const RegionalWindow window = worldwindow(k);
        const FineChannels channels = rasterisechannels(window, network, parameters(0.3, 1), sealevel);
        const std::size_t sink = centrecell(window, 10, 2);
        EXPECT_EQ(channels.channel[sink], 1);
        EXPECT_EQ(channels.direction[sink], 0);
        EXPECT_EQ(channels.bed_m[sink], 50.0f);
        EXPECT_EQ(channels.globalcell[sink], static_cast<int32_t>(globalindex(11, 2)));
        const WalkReport report = walkpinned(network, window, channels);
        EXPECT_EQ(report.walks, 1);
        EXPECT_EQ(report.wrongoutlet, 0);
        EXPECT_EQ(channelcount(channels), static_cast<int>(std::count(channels.globalcell.begin(),
                                                                       channels.globalcell.end(),
                                                                       static_cast<int32_t>(globalindex(11, 2)))));
    }
}

TEST_F(Landscape, ChannelNetworkSeaCornersEndLinks)
{
    // Two diagonals along a coast. (5, 3) -> (6, 4) is land to land, but its side
    // neighbour (6, 3) is sea; (8, 3) -> (9, 4) runs into the sea, and its side
    // neighbour (9, 3) is sea as well. A diagonal segment runs through the corner its
    // two side neighbours share, so a meander that bulges towards (6, 3) or (9, 3)
    // enters that sea block, and the link must end there at the sea level.
    GlobalNetwork network = seaworld();
    setland(network, 5, 4, west, 10.0f, 300.0f);
    setland(network, 6, 4, west, 20.0f, 200.0f);
    setland(network, 5, 3, southeast, 40.0f, 100.0f);
    setland(network, 8, 4, east, 10.0f, 100.0f);
    setland(network, 8, 3, southeast, 30.0f, 100.0f);
    const int32_t coastal = static_cast<int32_t>(globalindex(5, 3));
    int cutatcorner = 0;
    for (int k : { 4, 5, 8 })
        for (std::uint64_t seed = 1; seed <= 64; seed++)
        {
            SCOPED_TRACE(testing::Message() << "k " << k << ", seed " << seed);
            const RegionalWindow window = worldwindow(k);
            const FineChannels channels = rasterisechannels(window, network, parameters(0.3, seed), sealevel);
            for (std::size_t cell = 0; cell < channels.channel.size(); cell++)
            {
                if (channels.channel[cell] == 0)
                    continue;
                const bool sea = network.sea[globalof(window, cell)] != 0;
                if (channels.direction[cell] != 0)
                {
                    ASSERT_FALSE(sea) << "fine cell " << cell;  // a cell that leads on lies on land
                    std::size_t next = 0;
                    ASSERT_TRUE(nextcell(window, cell, channels.direction[cell], next));
                    ASSERT_EQ(channels.channel[next], 1);
                    ASSERT_LE(channels.bed_m[next], channels.bed_m[cell]);
                    continue;
                }
                ASSERT_TRUE(sea) << "fine cell " << cell;  // no land sink here: every end is in the sea
                ASSERT_EQ(channels.bed_m[cell], sealevel);
                if (channels.globalcell[cell] == coastal && globalof(window, cell) == globalindex(6, 3))
                    cutatcorner++;
            }
            // Every walk from a pinned centre ends in a sea block (the coastal link's,
            // when cut, in (6, 3)'s rather than its global outlet's).
            for (std::size_t start : { globalindex(5, 3), globalindex(6, 4), globalindex(5, 4), globalindex(8, 3),
                                       globalindex(8, 4) })
            {
                std::size_t walk = centrecell(window, static_cast<int>(start % worldwidth),
                                              static_cast<int>(start / worldwidth));
                int steps = 0;
                while (channels.direction[walk] != 0 && steps < 8 * k)
                {
                    ASSERT_TRUE(nextcell(window, walk, channels.direction[walk], walk));
                    steps++;
                }
                EXPECT_LT(steps, 8 * k);
                EXPECT_NE(network.sea[globalof(window, walk)], 0);
            }
        }
    // The cut at a sea side neighbour's corner happens for some seeds.
    EXPECT_GT(cutatcorner, 0);
}

TEST_F(Landscape, ChannelNetworkFlatChainTiesByOutletDistance)
{
    // A flat chain: X = (6, 4) -> (7, 4) -> (8, 4) east, all at 50 m, then (8, 4) into
    // the sea at (9, 4); the diagonal child (7, 3) -> X (south-west), also at 50 m,
    // has the smaller index (55 against 70), and (6, 3) is a land sink (no link), so
    // the child's side neighbours are land. Near X the child's path and X's outgoing
    // path share cells at the same bed. Ties by index alone would give them to the
    // child, whose direction leads back to X's centre: a loop. Ties by the distance to
    // the outlet give them to X's link (3 links from the sea against the child's 4).
    GlobalNetwork network = seaworld();
    setland(network, 7, 3, southwest, 50.0f, 100.0f);
    setland(network, 6, 4, east, 50.0f, 200.0f);
    setland(network, 7, 4, east, 50.0f, 300.0f);
    setland(network, 8, 4, east, 50.0f, 400.0f);
    setland(network, 6, 3, 0, 50.0f, 100.0f);
    int sharedflat = 0;
    for (int k : { 4, 5, 8 })
        for (std::uint64_t seed = 1; seed <= 64; seed++)
        {
            SCOPED_TRACE(testing::Message() << "k " << k << ", seed " << seed);
            const RegionalWindow window = worldwindow(k);
            const FineChannels channels = rasterisechannels(window, network, parameters(0.3, seed), sealevel);
            for (std::size_t cell = 0; cell < channels.channel.size(); cell++)
            {
                std::size_t next = 0;
                if (channels.channel[cell] == 0 || channels.direction[cell] == 0 ||
                    !nextcell(window, cell, channels.direction[cell], next))
                    continue;
                if (channels.globalcell[next] != channels.globalcell[cell] && !iscentre(window, next) &&
                    channels.bed_m[next] == channels.bed_m[cell])
                    sharedflat++;
            }
            const WalkReport report = walkpinned(network, window, channels);
            ASSERT_EQ(report.walks, 4);
            ASSERT_EQ(report.nonadjacent, 0);
            ASSERT_EQ(report.wrongoutlet, 0);  // includes: no loop
        }
    EXPECT_GT(sharedflat, 0);
}

TEST_F(Landscape, ChannelNetworkRejectsBadInput)
{
    const GlobalNetwork network = syntheticnetwork();
    EXPECT_THROW(rasterisechannels(worldwindow(0), network, parameters(0.3, 1), sealevel), std::invalid_argument);
    EXPECT_THROW(rasterisechannels(worldwindow(4), network, parameters(-0.1, 1), sealevel), std::invalid_argument);
    EXPECT_THROW(rasterisechannels(domainwindow(4, 0, worldwidth + 1, 0, worldheight), network, parameters(0.3, 1),
                                   sealevel),
                 std::invalid_argument);
    EXPECT_THROW(rasterisechannels(domainwindow(4, 0, 4, 6, 3), network, parameters(0.3, 1), sealevel),
                 std::invalid_argument);
    GlobalNetwork shortnetwork = network;
    shortnetwork.julyflow.pop_back();
    EXPECT_THROW(rasterisechannels(worldwindow(4), shortnetwork, parameters(0.3, 1), sealevel), std::invalid_argument);
    // A land receiver above its upstream cell on the filled surface.
    GlobalNetwork rising = seaworld();
    setland(rising, 6, 4, west, 30.0f, 200.0f);
    setland(rising, 7, 4, west, 20.0f, 100.0f);
    EXPECT_THROW(rasterisechannels(worldwindow(4), rising, parameters(0.3, 1), sealevel), std::invalid_argument);
    // A filled elevation that is not finite on a drawn link.
    GlobalNetwork notfinite = seaworld();
    setland(notfinite, 6, 4, west, std::numeric_limits<float>::quiet_NaN(), 100.0f);
    EXPECT_THROW(rasterisechannels(worldwindow(4), notfinite, parameters(0.3, 1), sealevel), std::invalid_argument);
    // The sea level above the filled elevation of a land cell that drains into the sea.
    GlobalNetwork belowsea = seaworld();
    setland(belowsea, 6, 4, west, -5.0f, 100.0f);
    EXPECT_THROW(rasterisechannels(worldwindow(4), belowsea, parameters(0.3, 1), sealevel), std::invalid_argument);
    // A sea level that is not finite.
    EXPECT_THROW(rasterisechannels(worldwindow(4), network, parameters(0.3, 1), std::numeric_limits<float>::infinity()),
                 std::invalid_argument);
    EXPECT_THROW(rasterisechannels(worldwindow(4), network, parameters(0.3, 1), std::numeric_limits<float>::quiet_NaN()),
                 std::invalid_argument);
    // Receivers that form a loop: (6, 4) and (7, 4) drain into each other on a flat.
    GlobalNetwork loop = seaworld();
    setland(loop, 6, 4, east, 10.0f, 100.0f);
    setland(loop, 7, 4, west, 10.0f, 100.0f);
    EXPECT_THROW(rasterisechannels(worldwindow(4), loop, parameters(0.3, 1), sealevel), std::invalid_argument);
    // ... but a bad link out of reach of the domain is not looked at.
    EXPECT_NO_THROW(rasterisechannels(domainwindow(4, 0, 2, 0, 2), notfinite, parameters(0.3, 1), sealevel));
}
