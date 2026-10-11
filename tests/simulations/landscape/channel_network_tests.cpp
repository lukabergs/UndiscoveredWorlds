// Tests of landscape/channel_network.hpp (erosion plan Phase E, step E2b; Phase R, steps
// R1 and R3; multiscale plan S3.10, the transit blocks), suite Landscape: TEST_F(Landscape, ChannelNetwork<What>) on synthetic
// 16 x 8 worlds with hand-made receiver trees, at scales 1 to 16.

#include "landscape_fixture.hpp"

#include "simulations/landscape/channel_network.hpp"
#include "simulations/landscape/regional_grid.hpp"
#include "simulations/landscape/sphere_metric.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <limits>
#include <set>
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

// The geometry of Phase R, R3 (junctions within 0.25 k of the centre, edge margin 0.1875,
// every diagonal through the shared corner): the tests before S3.10 test it as it was.
ChannelParameters parameters(double meanderamplitude, std::uint64_t seed, double pinarea_km2 = 0.0, int springs = 0)
{
    ChannelParameters result;
    result.junctionjitter = 0.25;
    result.edgemargin = 0.1875;
    result.diagonaltransit = 0;
    result.pinarea_km2 = pinarea_km2;
    result.meanderamplitude = meanderamplitude;
    result.seed = seed;
    result.springspercell = springs;
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

int columnof(const RegionalWindow& window, std::size_t index)
{
    return static_cast<int>(index % static_cast<std::size_t>(window.domainfinewidth()));
}

int rowof(const RegionalWindow& window, std::size_t index)
{
    return static_cast<int>(index / static_cast<std::size_t>(window.domainfinewidth()));
}

// The junction fine cell of global (i, j) (channeljunction), which must lie in the domain.
std::size_t junctioncell(const RegionalWindow& window, const ChannelParameters& settings, int i, int j)
{
    const landscape::FinePosition junction = landscape::channeljunction(window, settings, i, j);
    std::size_t index = 0;
    EXPECT_TRUE(domaincell(window, junction.column, junction.row, index));
    return index;
}

// The global cell holding a domain fine cell.
std::size_t globalof(const RegionalWindow& window, std::size_t index)
{
    return globalindex(static_cast<int>(absolutecolumn(window, columnof(window, index)) / window.scale),
                       static_cast<int>(absoluterow(window, rowof(window, index)) / window.scale));
}

bool isjunction(const RegionalWindow& window, const ChannelParameters& settings, std::size_t index)
{
    const std::size_t global = globalof(window, index);
    const landscape::FinePosition junction = landscape::channeljunction(
        window, settings, static_cast<int>(global % worldwidth), static_cast<int>(global / worldwidth));
    return junction.column == absolutecolumn(window, columnof(window, index)) &&
           junction.row == absoluterow(window, rowof(window, index));
}

// Where a chain may end: any fine cell of a sea cell's block, or a land sink's junction.
bool isend(const GlobalNetwork& network, const RegionalWindow& window, const ChannelParameters& settings,
           std::size_t index)
{
    const std::size_t cell = globalof(window, index);
    if (network.sea[cell] != 0)
        return true;
    return network.direction[cell] == 0 && isjunction(window, settings, index);
}

// The domain fine cell a direction leads to (false off the domain).
bool nextcell(const RegionalWindow& window, std::size_t index, int direction, std::size_t& next)
{
    return domaincell(window, absolutecolumn(window, columnof(window, index)) + landscape::directionoffsetx(direction),
                      absoluterow(window, rowof(window, index)) + landscape::directionoffsety(direction), next);
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

// Walks the directions from every pinned global cell's junction (the window must be
// the whole world): counts the walks, the changes of owner (globalcell) to a link not
// adjacent on the global tree, and the walks that do not end in the starting cell's
// global outlet (a sea cell's block or a land sink's junction) or loop.
struct WalkReport
{
    int walks = 0;
    int nonadjacent = 0;
    int wrongoutlet = 0;
};

WalkReport walkpinned(const GlobalNetwork& network, const RegionalWindow& window, const ChannelParameters& settings,
                      const FineChannels& channels)
{
    WalkReport report;
    const int limit = static_cast<int>(channels.channel.size());
    for (int j = 0; j < worldheight; j++)
        for (int i = 0; i < worldwidth; i++)
        {
            const std::size_t start = globalindex(i, j);
            if (network.sea[start] != 0 || network.direction[start] == 0 ||
                !(network.drainagearea_km2[start] >= settings.pinarea_km2))
                continue;
            report.walks++;
            std::size_t walk = junctioncell(window, settings, i, j);
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
            if (!ok || steps > limit || globalof(window, walk) != outlet || !isend(network, window, settings, walk))
                report.wrongoutlet++;
        }
    return report;
}

// Steps from a channel cell to a cell owned by another link that is not a junction.
int contestedsteps(const RegionalWindow& window, const ChannelParameters& settings, const FineChannels& channels)
{
    int count = 0;
    for (std::size_t cell = 0; cell < channels.channel.size(); cell++)
    {
        std::size_t next = 0;
        if (channels.channel[cell] == 0 || channels.direction[cell] == 0 ||
            !nextcell(window, cell, channels.direction[cell], next))
            continue;
        if (channels.globalcell[next] != channels.globalcell[cell] && !isjunction(window, settings, next))
            count++;
    }
    return count;
}

// Channel cells outside the blocks of their owner's link (its upstream cell and its
// receiver; a spring lies in its own cell's block), plus the steps that break the legs'
// order: a link's path is its outflow leg in the upstream block, then its inflow leg in
// the receiver's block, so a step between two cells of the same owner never leads from
// the receiver's block back into the upstream block, and a step from the upstream block
// into the receiver's block crosses to the neighbour block in the link's direction.
int cellsoutsidetheirblocks(const GlobalNetwork& network, const RegionalWindow& window, const FineChannels& channels)
{
    int outside = 0;
    for (std::size_t cell = 0; cell < channels.channel.size(); cell++)
    {
        if (channels.channel[cell] == 0)
            continue;
        const std::size_t owner = static_cast<std::size_t>(channels.globalcell[cell]);
        const std::size_t receiver = globalreceiver(network, owner);
        const std::size_t block = globalof(window, cell);
        if (block != owner && block != receiver)
        {
            outside++;
            continue;
        }
        std::size_t next = 0;
        if (channels.direction[cell] == 0 || !nextcell(window, cell, channels.direction[cell], next) ||
            channels.globalcell[next] != channels.globalcell[cell])
            continue;
        const std::size_t nextblock = globalof(window, next);
        if (block == receiver && nextblock == owner && receiver != owner)
            outside++;  // back from the inflow leg into the outflow leg's block
        if (block == owner && nextblock == receiver && receiver != owner &&
            channels.direction[cell] != network.direction[owner])
            outside++;  // the crossing is not the step from the exit point to the entry point
    }
    return outside;
}

struct Case
{
    int scale;
    double meanderamplitude;
    std::uint64_t seed;
    int springs;
};

const Case cases[] = { { 4, 0.3, 1, 0 },  { 5, 0.3, 1, 0 },  { 4, 0.45, 2, 0 }, { 5, 0.45, 2, 0 },
                       { 8, 0.3, 3, 2 },  { 4, 0.0, 1, 0 },  { 5, 0.0, 1, 0 },  { 1, 0.3, 1, 0 },
                       { 2, 0.3, 1, 1 },  { 3, 0.3, 1, 1 },  { 16, 0.2, 4, 3 } };

ChannelParameters settingsof(const Case& c) { return parameters(c.meanderamplitude, c.seed, 0.0, c.springs); }

int channelcount(const FineChannels& channels)
{
    return static_cast<int>(std::count(channels.channel.begin(), channels.channel.end(), uint8_t{ 1 }));
}

// Checks that every channel cell's bed is at or above the next cell's.
int risingsteps(const RegionalWindow& window, const FineChannels& channels)
{
    int rising = 0;
    for (std::size_t cell = 0; cell < channels.channel.size(); cell++)
    {
        std::size_t next = 0;
        if (channels.channel[cell] == 0 || channels.direction[cell] == 0 ||
            !nextcell(window, cell, channels.direction[cell], next))
            continue;
        rising += channels.bed_m[next] > channels.bed_m[cell] ? 1 : 0;
    }
    return rising;
}
}

TEST_F(Landscape, ChannelNetworkDirectionsLeadToTheSea)
{
    const GlobalNetwork network = syntheticnetwork();
    int contested = 0;
    for (const Case& c : cases)
    {
        SCOPED_TRACE(testing::Message() << "k " << c.scale << ", meander " << c.meanderamplitude << ", seed " << c.seed
                                        << ", springs " << c.springs);
        const RegionalWindow window = worldwindow(c.scale);
        const ChannelParameters settings = settingsof(c);
        const FineChannels channels = rasterisechannels(window, network, settings, sealevel);
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
            EXPECT_TRUE(isend(network, window, settings, walk)) << "from fine cell " << cell;
        }
        // Every owner change along the way is to a link adjacent on the global tree,
        // and every pinned cell's chain ends in its own global outlet.
        const WalkReport report = walkpinned(network, window, settings, channels);
        EXPECT_EQ(report.walks, 120);
        EXPECT_EQ(report.nonadjacent, 0);
        EXPECT_EQ(report.wrongoutlet, 0);
        contested += contestedsteps(window, settings, channels);
    }
    // The ownership rule is exercised: some chains run into another link's cell
    // that is not a junction.
    EXPECT_GT(contested, 0);
}

TEST_F(Landscape, ChannelNetworkBedNeverRises)
{
    const GlobalNetwork network = syntheticnetwork();
    for (const Case& c : cases)
    {
        SCOPED_TRACE(testing::Message() << "k " << c.scale << ", meander " << c.meanderamplitude << ", seed " << c.seed);
        const RegionalWindow window = worldwindow(c.scale);
        const FineChannels channels = rasterisechannels(window, network, settingsof(c), sealevel);
        EXPECT_EQ(risingsteps(window, channels), 0);
        for (std::size_t cell = 0; cell < channels.channel.size(); cell++)
            if (channels.channel[cell] != 0)
                EXPECT_GE(channels.bed_m[cell], sealevel);
    }
}

TEST_F(Landscape, ChannelNetworkPathsStayInTheirBlocks)
{
    // Every channel cell lies in the block of its owner's upstream cell or receiver
    // (springs in their own cell's): a link occupies only its two cells' blocks, so
    // links not adjacent on the global tree never share a fine cell. With and without a
    // start surface and a speckled fine sea.
    const GlobalNetwork network = syntheticnetwork();
    for (const Case& c : cases)
    {
        SCOPED_TRACE(testing::Message() << "k " << c.scale << ", meander " << c.meanderamplitude << ", seed " << c.seed);
        const RegionalWindow window = worldwindow(c.scale);
        const ChannelParameters settings = settingsof(c);
        EXPECT_EQ(cellsoutsidetheirblocks(network, window, rasterisechannels(window, network, settings, sealevel)), 0);
        std::vector<float> start(static_cast<std::size_t>(window.domainfinewidth()) * window.domainfineheight());
        std::vector<uint8_t> finesea(start.size());
        for (std::size_t cell = 0; cell < start.size(); cell++)
        {
            start[cell] = static_cast<float>((cell * 2654435761u) % 400) - 50.0f;
            finesea[cell] = network.sea[globalof(window, cell)] != 0 || (cell * 40503u) % 97 == 0 ? 1 : 0;
        }
        const FineChannels channels = rasterisechannels(window, network, settings, sealevel, start, finesea);
        EXPECT_EQ(cellsoutsidetheirblocks(network, window, channels), 0);
        EXPECT_EQ(risingsteps(window, channels), 0);
    }
}

TEST_F(Landscape, ChannelNetworkJunctionsArePinned)
{
    // A land cell's junction belongs to its own outgoing link, carries its flows and,
    // without a start surface, lies at its filled elevation; with one, at the start
    // surface held between the receiver's filled elevation (the sea level's next float
    // for a sea receiver) and its own. A sea cell's junction is not reached when it lies
    // off its block's edge (a link into the sea ends at its entry point).
    const GlobalNetwork network = syntheticnetwork();
    const float abovesea = std::nextafter(sealevel, std::numeric_limits<float>::infinity());
    for (const Case& c : cases)
    {
        SCOPED_TRACE(testing::Message() << "k " << c.scale << ", meander " << c.meanderamplitude << ", seed " << c.seed);
        const RegionalWindow window = worldwindow(c.scale);
        const ChannelParameters settings = settingsof(c);
        std::vector<float> start(static_cast<std::size_t>(window.domainfinewidth()) * window.domainfineheight());
        for (std::size_t cell = 0; cell < start.size(); cell++)
            start[cell] = static_cast<float>((cell * 2654435761u) % 300) - 20.0f;
        for (const bool withstart : { false, true })
        {
            const FineChannels channels = withstart ? rasterisechannels(window, network, settings, sealevel, start)
                                                    : rasterisechannels(window, network, settings, sealevel);
            int clamped = 0, followed = 0;
            for (int j = 0; j < worldheight; j++)
                for (int i = 0; i < worldwidth; i++)
                {
                    const std::size_t cell = globalindex(i, j);
                    const std::size_t junction = junctioncell(window, settings, i, j);
                    if (network.sea[cell] != 0)
                    {
                        const int x = columnof(window, junction) % c.scale, y = rowof(window, junction) % c.scale;
                        if (x > 0 && y > 0 && x < c.scale - 1 && y < c.scale - 1)
                            EXPECT_EQ(channels.channel[junction], 0) << "sea " << i << ", " << j;
                        continue;
                    }
                    EXPECT_EQ(channels.channel[junction], 1) << "land " << i << ", " << j;
                    EXPECT_EQ(channels.globalcell[junction], static_cast<int32_t>(cell));
                    EXPECT_NE(channels.direction[junction], 0);
                    EXPECT_EQ(channels.discharge_m3s[junction], network.discharge_m3s[cell]);
                    EXPECT_EQ(channels.januaryflow[junction], network.januaryflow[cell]);
                    EXPECT_EQ(channels.julyflow[junction], network.julyflow[cell]);
                    const float filled = network.filledelevation_m[cell];
                    if (!withstart)
                    {
                        EXPECT_EQ(channels.bed_m[junction], filled) << "land " << i << ", " << j;
                        continue;
                    }
                    const std::size_t receiver = globalreceiver(network, cell);
                    const float base = network.sea[receiver] != 0 ? sealevel : network.filledelevation_m[receiver];
                    const float lo = std::max(base, std::min(filled, abovesea));
                    const float expected = std::max(lo, std::min(filled, start[junction]));
                    EXPECT_EQ(channels.bed_m[junction], expected) << "land " << i << ", " << j;
                    clamped += start[junction] < lo ? 1 : 0;
                    followed += start[junction] >= lo && start[junction] < filled ? 1 : 0;
                }
            if (withstart)
            {
                EXPECT_GT(clamped, 0);
                EXPECT_GT(followed, 0);
            }
        }
    }
}

TEST_F(Landscape, ChannelNetworkEdgePointSharedByBothCells)
{
    // The point where a cardinal link crosses an edge is drawn from the edge, so it is
    // the same whichever of the two cells drains into the other: (5, 3) east into (6, 3)
    // against (6, 3) west into (5, 3) (a vertical edge; also across the date line, (15, 3)
    // east into (0, 3) against (0, 3) west into (15, 3)), and (5, 3) south into (5, 4)
    // against (5, 4) north into (5, 3) (a horizontal edge). It lies between
    // edgemargin x k and (1 - edgemargin) x k along the edge. A diagonal link steps
    // through the shared corner: (5, 3) south-east into (6, 4).
    for (int k : { 4, 5, 8, 16 })
    {
        std::set<std::int64_t> offsets;
        for (std::uint64_t seed = 1; seed <= 32; seed++)
        {
            SCOPED_TRACE(testing::Message() << "k " << k << ", seed " << seed);
            const RegionalWindow window = worldwindow(k);
            const ChannelParameters settings = parameters(0.3, seed);
            // The crossing: the cell of `column` (or `row`) owned by `owner` whose
            // direction is `direction`, returned as its row (or column).
            auto crossing = [&](const FineChannels& channels, bool vertical, std::int64_t line, std::size_t owner,
                                int direction) {
                std::int64_t found = -1;
                int count = 0;
                for (std::size_t cell = 0; cell < channels.channel.size(); cell++)
                {
                    const std::int64_t x = absolutecolumn(window, columnof(window, cell));
                    const std::int64_t y = absoluterow(window, rowof(window, cell));
                    if (channels.channel[cell] == 0 || channels.globalcell[cell] != static_cast<int32_t>(owner) ||
                        channels.direction[cell] != direction || (vertical ? x : y) != line)
                        continue;
                    found = vertical ? y : x;
                    count++;
                }
                EXPECT_EQ(count, 1);
                return found;
            };
            GlobalNetwork eastward = seaworld();
            setland(eastward, 5, 3, east, 20.0f, 100.0f);
            setland(eastward, 6, 3, east, 10.0f, 200.0f);
            GlobalNetwork westward = seaworld();
            setland(westward, 6, 3, west, 20.0f, 100.0f);
            setland(westward, 5, 3, west, 10.0f, 200.0f);
            const std::int64_t a = crossing(rasterisechannels(window, eastward, settings, sealevel), true, 6 * k - 1,
                                            globalindex(5, 3), east);
            const std::int64_t b = crossing(rasterisechannels(window, westward, settings, sealevel), true, 6 * k,
                                            globalindex(6, 3), west);
            EXPECT_EQ(a, b);
            const std::int64_t offset = a - 3 * k;
            EXPECT_GE(offset, static_cast<std::int64_t>(std::floor(settings.edgemargin * k)));
            EXPECT_LT(offset, k - static_cast<std::int64_t>(std::floor(settings.edgemargin * k)));
            offsets.insert(offset);

            // Across the date line: (15, 3) east into (0, 3) against (0, 3) west into (15, 3).
            GlobalNetwork eastacross = seaworld();
            setland(eastacross, 15, 3, east, 20.0f, 100.0f);
            setland(eastacross, 0, 3, east, 10.0f, 200.0f);
            GlobalNetwork westacross = seaworld();
            setland(westacross, 0, 3, west, 20.0f, 100.0f);
            setland(westacross, 15, 3, west, 10.0f, 200.0f);
            const std::int64_t c = crossing(rasterisechannels(window, eastacross, settings, sealevel), true,
                                            static_cast<std::int64_t>(worldwidth) * k - 1, globalindex(15, 3), east);
            const std::int64_t d =
                crossing(rasterisechannels(window, westacross, settings, sealevel), true, 0, globalindex(0, 3), west);
            EXPECT_EQ(c, d);
            EXPECT_GE(c - 3 * k, static_cast<std::int64_t>(std::floor(settings.edgemargin * k)));
            EXPECT_LT(c - 3 * k, k - static_cast<std::int64_t>(std::floor(settings.edgemargin * k)));

            GlobalNetwork southward = seaworld();
            setland(southward, 5, 3, south, 20.0f, 100.0f);
            setland(southward, 5, 4, south, 10.0f, 200.0f);
            GlobalNetwork northward = seaworld();
            setland(northward, 5, 4, north, 20.0f, 100.0f);
            setland(northward, 5, 3, north, 10.0f, 200.0f);
            EXPECT_EQ(crossing(rasterisechannels(window, southward, settings, sealevel), false, 4 * k - 1,
                               globalindex(5, 3), south),
                      crossing(rasterisechannels(window, northward, settings, sealevel), false, 4 * k, globalindex(5, 4),
                               north));

            GlobalNetwork diagonal = seaworld();
            setland(diagonal, 5, 3, southeast, 20.0f, 100.0f);
            setland(diagonal, 6, 4, east, 10.0f, 200.0f);
            const FineChannels channels = rasterisechannels(window, diagonal, settings, sealevel);
            std::size_t corner = 0, across = 0;
            ASSERT_TRUE(domaincell(window, 6 * k - 1, 4 * k - 1, corner));
            ASSERT_TRUE(domaincell(window, 6 * k, 4 * k, across));
            EXPECT_EQ(channels.globalcell[corner], static_cast<int32_t>(globalindex(5, 3)));
            EXPECT_EQ(channels.direction[corner], southeast);
            EXPECT_EQ(channels.channel[across], 1);
        }
        if (k >= 8)
            EXPECT_GT(offsets.size(), 2u) << "k " << k;
    }
}

TEST_F(Landscape, ChannelNetworkSmallLinksAreAbsent)
{
    const GlobalNetwork network = syntheticnetwork();
    constexpr double pinarea_km2 = 250.0;  // headwaters (100 km2) and cells with one headwater donor (200 km2)
    for (int scale : { 4, 5 })
    {
        const RegionalWindow window = worldwindow(scale);
        const ChannelParameters settings = parameters(0.3, 1, pinarea_km2);
        const FineChannels channels = rasterisechannels(window, network, settings, sealevel);
        int pinned = 0;
        int unpinned = 0;
        for (int j = 0; j < worldheight; j++)
            for (int i = 1; i < worldwidth; i++)
            {
                const std::size_t cell = globalindex(i, j);
                const bool expected = network.drainagearea_km2[cell] >= pinarea_km2;
                (expected ? pinned : unpinned)++;
                if (expected)
                    EXPECT_EQ(channels.globalcell[junctioncell(window, settings, i, j)], static_cast<int32_t>(cell));
                const bool drawn = std::find(channels.globalcell.begin(), channels.globalcell.end(),
                                             static_cast<int32_t>(cell)) != channels.globalcell.end();
                EXPECT_EQ(drawn, expected) << i << ", " << j;
            }
        // Of the 120 land cells, 41 reach 250 km2: rows 2 (columns 1 .. 9) and 5
        // (columns 1 .. 10), columns 13 and 14 (rows 0 .. 6) and column 15.
        EXPECT_EQ(pinned, 41);
        EXPECT_EQ(unpinned, 79);
        const WalkReport report = walkpinned(network, window, settings, channels);
        EXPECT_EQ(report.walks, 41);
        EXPECT_EQ(report.nonadjacent, 0);
        EXPECT_EQ(report.wrongoutlet, 0);
    }
}

TEST_F(Landscape, ChannelNetworkLinkAcrossTheDateLine)
{
    // (15, 3) -> (0, 3) runs east across the date line into the sea of column 0: from
    // (15, 3)'s junction, one 8-neighbour step at a time, to fine column 0, its entry
    // point and first cell in the sea cell's block, where it ends at the sea level; in
    // the whole world and in a window across the date line.
    const GlobalNetwork network = syntheticnetwork();
    for (int k : { 4, 5, 16 })
        for (double meander : { 0.0, 0.45 })
        {
            SCOPED_TRACE(testing::Message() << "k " << k << ", meander " << meander);
            const std::int64_t planet = static_cast<std::int64_t>(k) * worldwidth;
            const ChannelParameters settings = parameters(meander, 2);
            for (const RegionalWindow& window : { worldwindow(k), domainwindow(k, 14, 4, 2, 3) })
            {
                const FineChannels channels = rasterisechannels(window, network, settings, sealevel);
                std::size_t walk = junctioncell(window, settings, 15, 3);
                bool crossed = false;
                int steps = 0;
                while (channels.direction[walk] != 0 && steps < 4 * k)
                {
                    const std::int64_t before = absolutecolumn(window, columnof(window, walk));
                    ASSERT_TRUE(nextcell(window, walk, channels.direction[walk], walk));
                    ASSERT_EQ(channels.channel[walk], 1);
                    crossed = crossed || (before == planet - 1 && absolutecolumn(window, columnof(window, walk)) == 0);
                    steps++;
                }
                EXPECT_TRUE(crossed);
                EXPECT_EQ(globalof(window, walk), globalindex(0, 3));
                EXPECT_EQ(absolutecolumn(window, columnof(window, walk)), 0);  // the block's first column
                EXPECT_EQ(channels.globalcell[walk], static_cast<int32_t>(globalindex(15, 3)));
                EXPECT_EQ(channels.bed_m[walk], sealevel);
            }
        }
}

TEST_F(Landscape, ChannelNetworkOverlappingWindowsAgree)
{
    // Every fine cell's values come from the pieces whose paths cross it, and every
    // window draws every link that reaches its domain: windows agree on every shared
    // fine cell, up to their edges (no margin is needed without a start surface).
    const GlobalNetwork network = syntheticnetwork();
    for (const Case& c : cases)
    {
        SCOPED_TRACE(testing::Message() << "k " << c.scale << ", meander " << c.meanderamplitude << ", seed " << c.seed
                                        << ", springs " << c.springs);
        const ChannelParameters settings = settingsof(c);
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
    const ChannelParameters settings = parameters(0.3, 7, 0.0, 2);
    const FineChannels first = rasterisechannels(window, network, settings, sealevel);
    const FineChannels second = rasterisechannels(window, network, settings, sealevel);
    EXPECT_EQ(first.channel, second.channel);
    EXPECT_EQ(first.direction, second.direction);
    EXPECT_EQ(first.bed_m, second.bed_m);
    EXPECT_EQ(first.discharge_m3s, second.discharge_m3s);
    EXPECT_EQ(first.januaryflow, second.januaryflow);
    EXPECT_EQ(first.julyflow, second.julyflow);
    EXPECT_EQ(first.globalcell, second.globalcell);
    // Another seed moves the junctions, edge points and meanders.
    const FineChannels other = rasterisechannels(window, network, parameters(0.3, 8, 0.0, 2), sealevel);
    EXPECT_NE(first.globalcell, other.globalcell);
}

TEST_F(Landscape, ChannelNetworkCrossingDiagonalsOnAFlat)
{
    // An X across the 2 x 2 block (6 .. 7, 3 .. 4) on a flat (every cell at 50 m, the
    // routing of a flood order): (6, 3) -> (7, 4), then east to the sea at (9, 4);
    // (7, 3) -> (6, 4), then west to the sea at (4, 4). The two diagonals pass through
    // different corner cells, so they share no fine cell and each chain reaches its own
    // outlet, with or without a start surface, at every scale and seed; the beds never
    // rise and no direction loops.
    GlobalNetwork network = seaworld();
    setland(network, 6, 3, southeast, 50.0f, 100.0f);
    setland(network, 7, 4, east, 50.0f, 200.0f);
    setland(network, 8, 4, east, 50.0f, 300.0f);
    setland(network, 7, 3, southwest, 50.0f, 100.0f);
    setland(network, 6, 4, west, 50.0f, 200.0f);
    setland(network, 5, 4, west, 50.0f, 300.0f);
    for (int k : { 1, 2, 3, 4, 5, 8, 16 })
        for (std::uint64_t seed = 1; seed <= 16; seed++)
            for (const bool withstart : { false, true })
            {
                SCOPED_TRACE(testing::Message() << "k " << k << ", seed " << seed << (withstart ? ", start" : ""));
                const RegionalWindow window = worldwindow(k);
                const ChannelParameters settings = parameters(0.3, seed, 0.0, static_cast<int>(seed % 3));
                std::vector<float> start;
                if (withstart)
                {
                    start.resize(static_cast<std::size_t>(window.domainfinewidth()) * window.domainfineheight());
                    for (std::size_t cell = 0; cell < start.size(); cell++)
                        start[cell] = static_cast<float>((cell * 2654435761u + seed) % 120);
                }
                const FineChannels channels = rasterisechannels(window, network, settings, sealevel, start);
                const WalkReport report = walkpinned(network, window, settings, channels);
                ASSERT_EQ(report.walks, 6);
                ASSERT_EQ(report.nonadjacent, 0);
                ASSERT_EQ(report.wrongoutlet, 0);  // includes: no loop
                ASSERT_EQ(risingsteps(window, channels), 0);
                ASSERT_EQ(cellsoutsidetheirblocks(network, window, channels), 0);
            }
}

TEST_F(Landscape, ChannelNetworkGentleTributaryKeepsDescending)
{
    // A gentle tributary and a steep, larger link meet at X = (6, 4), which drains
    // west to the sea: A = (7, 4) -> X, 1000 m above X and 5000 km2; B = (7, 3) -> X
    // (south-west), 1 m above X and 100 km2. Near X's junction their inflow legs share
    // fine cells, where A's bed is far above B's. Under ownership by the larger area,
    // B's chain would step from its own cell up onto A's; under the lowest bed it does
    // not. A's path alone is drawn by making B a sink, so the cells A would hold are known.
    GlobalNetwork network = seaworld();
    setland(network, 6, 4, west, 20.0f, 5200.0f);
    setland(network, 5, 4, west, 10.0f, 5300.0f);
    setland(network, 7, 4, west, 1020.0f, 5000.0f);
    setland(network, 7, 3, southwest, 21.0f, 100.0f);
    setland(network, 6, 3, 0, 30.0f, 100.0f);  // a land sink beside B's diagonal
    GlobalNetwork alone = network;
    alone.direction[globalindex(7, 3)] = 0;
    int wouldrise = 0;
    for (int k : { 4, 5, 8 })
        for (std::uint64_t seed = 1; seed <= 64; seed++)
        {
            SCOPED_TRACE(testing::Message() << "k " << k << ", seed " << seed);
            const RegionalWindow window = worldwindow(k);
            const ChannelParameters settings = parameters(0.3, seed);
            const FineChannels both = rasterisechannels(window, network, settings, sealevel);
            const FineChannels steep = rasterisechannels(window, alone, settings, sealevel);
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
            const WalkReport report = walkpinned(network, window, settings, both);
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
        const ChannelParameters settings = parameters(0.3, 1);
        const FineChannels channels = rasterisechannels(window, network, settings, sealevel);
        const std::size_t sink = junctioncell(window, settings, 10, 2);
        EXPECT_EQ(channels.channel[sink], 1);
        EXPECT_EQ(channels.direction[sink], 0);
        EXPECT_EQ(channels.bed_m[sink], 50.0f);
        EXPECT_EQ(channels.globalcell[sink], static_cast<int32_t>(globalindex(11, 2)));
        const WalkReport report = walkpinned(network, window, settings, channels);
        EXPECT_EQ(report.walks, 1);
        EXPECT_EQ(report.wrongoutlet, 0);
        EXPECT_EQ(channelcount(channels), static_cast<int>(std::count(channels.globalcell.begin(),
                                                                       channels.globalcell.end(),
                                                                       static_cast<int32_t>(globalindex(11, 2)))));
    }
}

TEST_F(Landscape, ChannelNetworkDiagonalsAlongACoastStayOnLand)
{
    // Two diagonals along a coast. (5, 3) -> (6, 4) is land to land, but its side
    // neighbour (6, 3) is sea; (8, 3) -> (9, 4) runs into the sea, and its side
    // neighbour (9, 3) is sea as well. A link occupies only its two cells' blocks (a
    // diagonal steps through the shared corner), so no link enters a side neighbour's
    // sea block: every chain ends in its own outlet, in the sea, at the sea level.
    GlobalNetwork network = seaworld();
    setland(network, 5, 4, west, 10.0f, 300.0f);
    setland(network, 6, 4, west, 20.0f, 200.0f);
    setland(network, 5, 3, southeast, 40.0f, 100.0f);
    setland(network, 8, 4, east, 10.0f, 100.0f);
    setland(network, 8, 3, southeast, 30.0f, 100.0f);
    for (int k : { 2, 3, 4, 5, 8, 16 })
        for (std::uint64_t seed = 1; seed <= 64; seed++)
        {
            SCOPED_TRACE(testing::Message() << "k " << k << ", seed " << seed);
            const RegionalWindow window = worldwindow(k);
            const ChannelParameters settings = parameters(0.45, seed);
            const FineChannels channels = rasterisechannels(window, network, settings, sealevel);
            for (std::size_t cell = 0; cell < channels.channel.size(); cell++)
            {
                if (channels.channel[cell] == 0)
                    continue;
                const std::size_t block = globalof(window, cell);
                ASSERT_TRUE(block != globalindex(6, 3) && block != globalindex(9, 3)) << "fine cell " << cell;
                if (channels.direction[cell] == 0)
                {
                    ASSERT_NE(network.sea[block], 0) << "fine cell " << cell;
                    ASSERT_EQ(channels.bed_m[cell], sealevel);
                }
            }
            const WalkReport report = walkpinned(network, window, settings, channels);
            ASSERT_EQ(report.walks, 5);
            ASSERT_EQ(report.nonadjacent, 0);
            ASSERT_EQ(report.wrongoutlet, 0);
        }
}

TEST_F(Landscape, ChannelNetworkFlatChainTiesByOutletDistance)
{
    // A flat chain: X = (6, 4) -> (7, 4) -> (8, 4) east, all at 50 m, then (8, 4) into
    // the sea at (9, 4); the diagonal child (7, 3) -> X (south-west), also at 50 m,
    // has the smaller index (55 against 70), and (6, 3) is a land sink (no link). In X's
    // block the child's inflow leg and X's outflow leg share cells at the same bed. Ties
    // by index alone would give them to the child, whose direction leads back to X's
    // junction: a loop. Ties by the distance to the outlet give them to X's link (3
    // links from the sea against the child's 4).
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
            const ChannelParameters settings = parameters(0.3, seed);
            const FineChannels channels = rasterisechannels(window, network, settings, sealevel);
            for (std::size_t cell = 0; cell < channels.channel.size(); cell++)
            {
                std::size_t next = 0;
                if (channels.channel[cell] == 0 || channels.direction[cell] == 0 ||
                    !nextcell(window, cell, channels.direction[cell], next))
                    continue;
                if (channels.globalcell[next] != channels.globalcell[cell] && !isjunction(window, settings, next) &&
                    channels.bed_m[next] == channels.bed_m[cell])
                    sharedflat++;
            }
            const WalkReport report = walkpinned(network, window, settings, channels);
            ASSERT_EQ(report.walks, 4);
            ASSERT_EQ(report.nonadjacent, 0);
            ASSERT_EQ(report.wrongoutlet, 0);  // includes: no loop
        }
    EXPECT_GT(sharedflat, 0);
}

TEST_F(Landscape, ChannelNetworkSprings)
{
    // A source cell (no drawn link into it) gets springspercell springs: legs from
    // hashed points of its block to its junction, owned by the cell, carrying its flows,
    // inside its block, at or above its junction's bed; they add channel cells to the
    // source blocks only, and every chain still ends in the sea without a rising bed.
    const GlobalNetwork network = syntheticnetwork();
    std::vector<uint8_t> source(network.sea.size(), 0);
    for (std::size_t cell = 0; cell < network.sea.size(); cell++)
        source[cell] = network.sea[cell] == 0 ? 1 : 0;
    for (std::size_t cell = 0; cell < network.sea.size(); cell++)
        if (network.sea[cell] == 0)
            source[globalreceiver(network, cell)] = 0;
    for (int k : { 3, 8, 16 })
        for (const bool withstart : { false, true })
        {
            SCOPED_TRACE(testing::Message() << "k " << k << (withstart ? ", start" : ""));
            const RegionalWindow window = worldwindow(k);
            std::vector<float> start;
            if (withstart)
            {
                start.resize(static_cast<std::size_t>(window.domainfinewidth()) * window.domainfineheight());
                for (std::size_t cell = 0; cell < start.size(); cell++)
                    start[cell] = static_cast<float>((cell * 2654435761u) % 500) - 20.0f;
            }
            const FineChannels plain = rasterisechannels(window, network, parameters(0.3, 5), sealevel, start);
            const ChannelParameters settings = parameters(0.3, 5, 0.0, 3);
            const FineChannels springs = rasterisechannels(window, network, settings, sealevel, start);
            int added = 0;
            for (std::size_t cell = 0; cell < springs.channel.size(); cell++)
            {
                if (springs.channel[cell] == plain.channel[cell])
                    continue;
                ASSERT_EQ(plain.channel[cell], 0);
                added++;
                const std::size_t block = globalof(window, cell);
                ASSERT_EQ(springs.globalcell[cell], static_cast<int32_t>(block));
                ASSERT_EQ(source[block], 1);
                ASSERT_EQ(springs.discharge_m3s[cell], network.discharge_m3s[block]);
                ASSERT_GE(springs.bed_m[cell],
                          springs.bed_m[junctioncell(window, settings, static_cast<int>(block % worldwidth),
                                                     static_cast<int>(block / worldwidth))]);
            }
            EXPECT_GT(added, 0);
            EXPECT_EQ(risingsteps(window, springs), 0);
            const WalkReport report = walkpinned(network, window, settings, springs);
            EXPECT_EQ(report.walks, 120);
            EXPECT_EQ(report.nonadjacent, 0);
            EXPECT_EQ(report.wrongoutlet, 0);
        }
}

TEST_F(Landscape, ChannelNetworkMeanderFollowsTheDischarge)
{
    // The midpoint displacement scales with min(1, (Q / reference)^exponent): a huge
    // reference draws the legs of a zero amplitude, a tiny one those of exponent 0.
    const GlobalNetwork network = syntheticnetwork();
    const RegionalWindow window = worldwindow(8);
    ChannelParameters straight = parameters(0.0, 3);
    ChannelParameters full = parameters(0.3, 3);
    ChannelParameters calm = full;
    calm.meanderdischargeexponent = 1.0;
    calm.meanderreference_m3s = 1e30;
    ChannelParameters wild = full;
    wild.meanderdischargeexponent = 1.0;
    wild.meanderreference_m3s = 1e-6;
    EXPECT_EQ(rasterisechannels(window, network, calm, sealevel).globalcell,
              rasterisechannels(window, network, straight, sealevel).globalcell);
    EXPECT_EQ(rasterisechannels(window, network, wild, sealevel).globalcell,
              rasterisechannels(window, network, full, sealevel).globalcell);
    EXPECT_NE(rasterisechannels(window, network, full, sealevel).globalcell,
              rasterisechannels(window, network, straight, sealevel).globalcell);
}

TEST_F(Landscape, ChannelNetworkStepDirections)
{
    // Reported, not asserted beyond the diagonals being used: the histogram of the
    // channel steps on the main world at k = 16 (its links run mostly west and east).
    const GlobalNetwork network = syntheticnetwork();
    const RegionalWindow window = worldwindow(16);
    const FineChannels channels = rasterisechannels(window, network, parameters(0.2, 1, 0.0, 1), sealevel);
    int counts[9] = {};
    for (std::size_t cell = 0; cell < channels.channel.size(); cell++)
        counts[channels.direction[cell]] += channels.channel[cell] != 0 ? 1 : 0;
    std::printf("[channels] k 16 step directions N %d NE %d E %d SE %d S %d SW %d W %d NW %d\n", counts[1], counts[2],
                counts[3], counts[4], counts[5], counts[6], counts[7], counts[8]);
    EXPECT_GT(counts[northeast] + counts[southeast] + counts[southwest] + counts[northwest], 0);
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
    // Bad geometry parameters.
    for (int bad = 0; bad < 7; bad++)
    {
        ChannelParameters settings = parameters(0.3, 1);
        switch (bad)
        {
        case 0: settings.meanderdischargeexponent = -1.0; break;
        case 1: settings.meanderreference_m3s = 0.0; break;
        case 2: settings.junctionjitter = -0.1; break;
        case 3: settings.edgemargin = 0.6; break;
        case 4: settings.edgemargin = std::numeric_limits<double>::quiet_NaN(); break;
        case 5: settings.diagonaltransit = 2; break;
        default: settings.springspercell = -1; break;
        }
        EXPECT_THROW(rasterisechannels(worldwindow(4), network, settings, sealevel), std::invalid_argument) << bad;
    }
    // S3.11: a negative or non-finite roughness, a persistence outside 0 .. 1.
    for (int bad = 0; bad < 5; bad++)
    {
        ChannelParameters settings = parameters(0.3, 1);
        switch (bad)
        {
        case 0: settings.meanderroughness = -0.01; break;
        case 1: settings.meanderroughness = std::numeric_limits<double>::infinity(); break;
        case 2: settings.meanderroughness = std::numeric_limits<double>::quiet_NaN(); break;
        case 3: settings.meanderroughnesspersistence = -0.1; break;
        default: settings.meanderroughnesspersistence = 1.5; break;
        }
        EXPECT_THROW(rasterisechannels(worldwindow(4), network, settings, sealevel), std::invalid_argument) << bad;
    }
    // A land receiver above its upstream cell on the filled surface.
    GlobalNetwork rising = seaworld();
    setland(rising, 6, 4, west, 30.0f, 200.0f);
    setland(rising, 7, 4, west, 20.0f, 100.0f);
    EXPECT_THROW(rasterisechannels(worldwindow(4), rising, parameters(0.3, 1), sealevel), std::invalid_argument);
    // ... also where a fine sea cuts the link before that receiver (the whole of the
    // receiver's block is fine sea), and a link cut at the fine coast whose upstream
    // cell lies below the sea level.
    {
        const RegionalWindow window = worldwindow(4);
        std::vector<uint8_t> finesea(static_cast<std::size_t>(window.domainfinewidth()) * window.domainfineheight(), 0);
        for (std::size_t cell = 0; cell < finesea.size(); cell++)
            finesea[cell] = globalof(window, cell) == globalindex(6, 4) ? 1 : 0;
        EXPECT_THROW(rasterisechannels(window, rising, parameters(0.3, 1), sealevel, {}, finesea), std::invalid_argument);
        GlobalNetwork low = seaworld();
        setland(low, 7, 4, west, -5.0f, 100.0f);
        setland(low, 6, 4, 0, -10.0f, 200.0f);
        EXPECT_NO_THROW(rasterisechannels(window, low, parameters(0.3, 1), sealevel));
        EXPECT_THROW(rasterisechannels(window, low, parameters(0.3, 1), sealevel, {}, finesea), std::invalid_argument);
    }
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

namespace
{
// Phase R, R1: domain fields of a window as functions of the absolute fine cell (x, y).
std::vector<float> finefield(const RegionalWindow& window, const std::function<float(std::int64_t, std::int64_t)>& value)
{
    std::vector<float> field(static_cast<std::size_t>(window.domainfinewidth()) * window.domainfineheight());
    for (int y = 0; y < window.domainfineheight(); y++)
        for (int x = 0; x < window.domainfinewidth(); x++)
            field[static_cast<std::size_t>(y) * window.domainfinewidth() + x] =
                value(absolutecolumn(window, x), absoluterow(window, y));
    return field;
}

std::vector<uint8_t> finemask(const RegionalWindow& window, const std::function<bool(std::int64_t, std::int64_t)>& test)
{
    const std::vector<float> field =
        finefield(window, [&](std::int64_t x, std::int64_t y) { return test(x, y) ? 1.0f : 0.0f; });
    std::vector<uint8_t> mask(field.size());
    for (std::size_t cell = 0; cell < field.size(); cell++)
        mask[cell] = field[cell] != 0.0f ? 1 : 0;
    return mask;
}

// The global sea blocks as a fine mask (where links end without a fine sea).
std::vector<uint8_t> seablocks(const GlobalNetwork& network, const RegionalWindow& window)
{
    return finemask(window, [&](std::int64_t x, std::int64_t y) {
        return network.sea[globalindex(static_cast<int>(x / window.scale), static_cast<int>(y / window.scale))] != 0;
    });
}

// A start surface of white noise in [low, high) from a hash of the absolute fine cell.
double hashunit(std::int64_t x, std::int64_t y, std::uint64_t salt = 0)
{
    std::uint64_t h = static_cast<std::uint64_t>(y * 1000003 + x) + 0x9e3779b97f4a7c15ULL + salt;
    h = (h ^ (h >> 30)) * 0xbf58476d1ce4e5b9ULL;
    h = (h ^ (h >> 27)) * 0x94d049bb133111ebULL;
    h ^= h >> 31;
    return static_cast<double>(h >> 11) / 9007199254740992.0;
}

std::vector<float> noisesurface(const RegionalWindow& window, float low, float high)
{
    return finefield(window, [&](std::int64_t x, std::int64_t y) {
        return low + (high - low) * static_cast<float>(hashunit(x, y));
    });
}

// Checks every channel cell against a fine sea mask: a cell that leads on is not sea,
// its bed above the sea level (every land cell of these worlds stands above it) and at
// or above the next cell's; a link's last cell is sea at the sea level or a land
// sink's junction; following the directions from any channel cell ends there without a
// loop. Returns the number of channel cells.
int checkfinecoast(const GlobalNetwork& network, const RegionalWindow& window, const ChannelParameters& settings,
                   const FineChannels& channels, const std::vector<uint8_t>& finesea)
{
    const int count = channelcount(channels);
    for (std::size_t cell = 0; cell < channels.channel.size(); cell++)
    {
        if (channels.channel[cell] == 0)
            continue;
        if (channels.direction[cell] == 0)
        {
            const std::size_t global = globalof(window, cell);
            const bool sink = network.sea[global] == 0 && network.direction[global] == 0 && isjunction(window, settings, cell);
            EXPECT_TRUE(finesea[cell] != 0 || sink) << "fine cell " << cell;
            if (finesea[cell] != 0)
                EXPECT_EQ(channels.bed_m[cell], sealevel) << "fine cell " << cell;
            continue;
        }
        EXPECT_EQ(finesea[cell], 0) << "a channel runs on through the sea at fine cell " << cell;
        EXPECT_GT(channels.bed_m[cell], sealevel) << "fine cell " << cell;
        std::size_t next = 0;
        EXPECT_TRUE(nextcell(window, cell, channels.direction[cell], next));
        EXPECT_LE(channels.bed_m[next], channels.bed_m[cell]) << "fine cell " << cell;
        std::size_t walk = cell;
        int steps = 0;
        while (channels.direction[walk] != 0 && steps <= count && nextcell(window, walk, channels.direction[walk], walk))
            steps++;
        EXPECT_LE(steps, count) << "loop from fine cell " << cell;
        EXPECT_EQ(channels.direction[walk], 0) << "from fine cell " << cell;
    }
    return count;
}

// The cells of the walk from a cell's junction along the directions to the end.
std::vector<std::size_t> walkfrom(const RegionalWindow& window, const FineChannels& channels, std::size_t start)
{
    std::vector<std::size_t> walk{ start };
    while (channels.direction[walk.back()] != 0 && walk.size() <= channels.channel.size())
    {
        std::size_t next = 0;
        EXPECT_TRUE(nextcell(window, walk.back(), channels.direction[walk.back()], next));
        walk.push_back(next);
    }
    return walk;
}
}

TEST_F(Landscape, ChannelNetworkEndsAtTheFineCoast)
{
    // The main world with a fine sea that reaches beyond the global sea column 0: the
    // west half of column 1's blocks (fine columns k .. k + k / 2 - 1). Links end at
    // their first fine sea cell, so some end inside column 1's land blocks (the fine
    // coast) at the sea level, and no channel cell that leads on is sea.
    const GlobalNetwork network = syntheticnetwork();
    for (const Case& c : cases)
    {
        SCOPED_TRACE(testing::Message() << "k " << c.scale << ", meander " << c.meanderamplitude << ", seed " << c.seed);
        const int k = c.scale;
        if (k < 2)
            continue;  // no fine coast inside a block of one cell
        const RegionalWindow window = worldwindow(k);
        const ChannelParameters settings = settingsof(c);
        const std::vector<uint8_t> finesea = finemask(window, [&](std::int64_t x, std::int64_t) { return x < k + k / 2; });
        for (const bool withstart : { false, true })
        {
            SCOPED_TRACE(withstart ? "with a start surface" : "without a start surface");
            const std::vector<float> start = withstart ? noisesurface(window, -50.0f, 400.0f) : std::vector<float>{};
            const FineChannels channels = rasterisechannels(window, network, settings, sealevel, start, finesea);
            EXPECT_GT(checkfinecoast(network, window, settings, channels, finesea), 0);
            int coastends = 0;
            for (std::size_t cell = 0; cell < channels.channel.size(); cell++)
                coastends += channels.channel[cell] != 0 && channels.direction[cell] == 0 &&
                                     absolutecolumn(window, columnof(window, cell)) >= k
                                 ? 1
                                 : 0;
            EXPECT_GT(coastends, 0);
        }
    }
}

TEST_F(Landscape, ChannelNetworkBedFollowsTheStartSurface)
{
    // One link at k = 8: (6, 4) at 100 m drains west into the land sink (5, 4) at 50 m,
    // under a start surface of white noise from 30 to 130 m. Along the walk from
    // (6, 4)'s junction: the junction bed b(6, 4) = max(50, min(100, s)), then the running
    // minimum of the start surface held at or above b(5, 4), and b(5, 4) at the sink's
    // junction. Then with (5, 4) draining west into the sea at (4, 4) and the start
    // surface below the sea level in its block: its beds are held at the next float above
    // the sea level and its link ends at the sea level in the sea block.
    const int k = 8;
    const RegionalWindow window = worldwindow(k);
    const float abovesea = std::nextafter(sealevel, std::numeric_limits<float>::infinity());
    for (std::uint64_t seed = 1; seed <= 16; seed++)
    {
        SCOPED_TRACE(testing::Message() << "seed " << seed);
        const ChannelParameters settings = parameters(0.3, seed);
        GlobalNetwork network = seaworld();
        setland(network, 6, 4, west, 100.0f, 200.0f);
        setland(network, 5, 4, 0, 50.0f, 300.0f);
        const std::vector<float> start = finefield(window, [&](std::int64_t x, std::int64_t y) {
            return 30.0f + 100.0f * static_cast<float>(hashunit(x, y, seed));
        });
        const FineChannels channels = rasterisechannels(window, network, settings, sealevel, start);
        const std::size_t upper = junctioncell(window, settings, 6, 4);
        const std::size_t lower = junctioncell(window, settings, 5, 4);
        const float top = std::max(50.0f, std::min(100.0f, start[upper]));
        const float end = std::max(std::min(50.0f, abovesea), std::min(50.0f, start[lower]));
        const std::vector<std::size_t> walk = walkfrom(window, channels, upper);
        ASSERT_EQ(walk.back(), lower);
        float running = top;
        for (std::size_t m = 0; m < walk.size(); m++)
        {
            running = std::min(running, start[walk[m]]);
            const float expected = m == 0 ? top : (m + 1 == walk.size() ? end : std::max(running, end));
            EXPECT_EQ(channels.bed_m[walk[m]], expected) << "step " << m;
            EXPECT_EQ(channels.globalcell[walk[m]], static_cast<int32_t>(globalindex(6, 4))) << "step " << m;
        }

        setland(network, 5, 4, west, 50.0f, 300.0f);
        const std::vector<float> low = finefield(window, [&](std::int64_t x, std::int64_t y) {
            return x < 6 * k ? -20.0f : 30.0f + 100.0f * static_cast<float>(hashunit(x, y, seed));
        });
        const FineChannels coastal = rasterisechannels(window, network, settings, sealevel, low);
        const std::vector<std::size_t> down = walkfrom(window, coastal, lower);
        EXPECT_EQ(coastal.bed_m[down.front()], abovesea);
        for (std::size_t m = 1; m + 1 < down.size(); m++)
            EXPECT_EQ(coastal.bed_m[down[m]], abovesea) << "step " << m;
        EXPECT_EQ(coastal.bed_m[down.back()], sealevel);
        EXPECT_EQ(globalof(window, down.back()), globalindex(4, 4));
        EXPECT_EQ(absolutecolumn(window, columnof(window, down.back())), 5 * k - 1);  // the sea block's east column
    }
}

TEST_F(Landscape, ChannelNetworkStartSurfaceBedNeverRises)
{
    // Under a noisy start surface (white noise from 50 m below the sea level to 150 m
    // above, so it dips below receivers and the sea level everywhere) the beds never
    // rise along the directions, stay at or above the sea level and at most the
    // upstream cell's filled elevation (springs: at most their top), and every pinned
    // chain still ends in its own outlet.
    const GlobalNetwork network = syntheticnetwork();
    for (const Case& c : cases)
    {
        SCOPED_TRACE(testing::Message() << "k " << c.scale << ", meander " << c.meanderamplitude << ", seed " << c.seed);
        const RegionalWindow window = worldwindow(c.scale);
        const ChannelParameters settings = settingsof(c);
        const std::vector<float> start = noisesurface(window, -50.0f, 150.0f);
        const FineChannels channels = rasterisechannels(window, network, settings, sealevel, start);
        int flat = 0;
        for (std::size_t cell = 0; cell < channels.channel.size(); cell++)
        {
            if (channels.channel[cell] == 0)
                continue;
            EXPECT_GE(channels.bed_m[cell], sealevel);
            if (c.springs == 0)
                EXPECT_LE(channels.bed_m[cell],
                          network.filledelevation_m[static_cast<std::size_t>(channels.globalcell[cell])]);
            std::size_t next = 0;
            if (channels.direction[cell] != 0 && nextcell(window, cell, channels.direction[cell], next))
                flat += channels.bed_m[next] == channels.bed_m[cell] ? 1 : 0;
        }
        if (c.scale > 2)
            EXPECT_GT(flat, 0);  // a running minimum is flat between the start surface's new lows
        EXPECT_GT(checkfinecoast(network, window, settings, channels, seablocks(network, window)), 0);
        const WalkReport report = walkpinned(network, window, settings, channels);
        EXPECT_EQ(report.walks, 120);
        EXPECT_EQ(report.nonadjacent, 0);
        EXPECT_EQ(report.wrongoutlet, 0);
    }
}

TEST_F(Landscape, ChannelNetworkRejectsWrongSizedSurfaces)
{
    const GlobalNetwork network = syntheticnetwork();
    const RegionalWindow window = worldwindow(4);
    const std::size_t cells = static_cast<std::size_t>(window.domainfinewidth()) * window.domainfineheight();
    EXPECT_THROW(rasterisechannels(window, network, parameters(0.3, 1), sealevel, std::vector<float>(cells - 1, 1.0f)),
                 std::invalid_argument);
    EXPECT_THROW(rasterisechannels(window, network, parameters(0.3, 1), sealevel, {}, std::vector<uint8_t>(cells + 1, 0)),
                 std::invalid_argument);
    EXPECT_NO_THROW(rasterisechannels(window, network, parameters(0.3, 1), sealevel, std::vector<float>(cells, 1.0f),
                                      std::vector<uint8_t>(cells, 0)));
}

TEST_F(Landscape, ChannelNetworkWindowsAgreeAwayFromTheirEdges)
{
    // The seam claim of R1 at channel level: given the same start surface and fine sea
    // (functions of the absolute fine cell: noise from 50 m below the sea level to 400 m
    // above, the fine sea of ChannelNetworkEndsAtTheFineCoast plus 3 % of the fine cells
    // as sea specks that cut links anywhere), windows agree on every channel field at the
    // fine cells at least one global cell inside both domains: a link lies in its two
    // cells' blocks (channel_network.hpp), so what a cell's pieces read lies in both
    // domains. Cells nearer an edge may differ (counted, not asserted).
    const GlobalNetwork network = syntheticnetwork();
    int compared = 0, comparedchannels = 0, edgediffer = 0;
    for (const Case& c : cases)
    {
        SCOPED_TRACE(testing::Message() << "k " << c.scale << ", meander " << c.meanderamplitude << ", seed " << c.seed);
        const int k = c.scale;
        const ChannelParameters settings = settingsof(c);
        const RegionalWindow windows[] = {
            domainwindow(k, 5, worldwidth, 0, worldheight),
            domainwindow(k, 13, 8, 1, 5),
            domainwindow(k, 2, 10, 0, 8),
            domainwindow(k, 9, 5, 3, 4),
            domainwindow(k, 0, 3, 6, 2),
        };
        std::vector<FineChannels> results;
        for (const RegionalWindow& window : windows)
        {
            const std::vector<uint8_t> finesea = finemask(window, [&](std::int64_t x, std::int64_t y) {
                return x < k + k / 2 || hashunit(x, y, 77) < 0.03;
            });
            results.push_back(rasterisechannels(window, network, settings, sealevel, noisesurface(window, -50.0f, 400.0f),
                                                finesea));
        }
        auto inside = [&](const RegionalWindow& window, int x, int y) {
            return x >= k && y >= k && x < window.domainfinewidth() - k && y < window.domainfineheight() - k;
        };
        for (std::size_t p = 0; p < results.size(); p++)
            for (std::size_t q = 0; q < results.size(); q++)
            {
                if (p == q)
                    continue;
                const RegionalWindow& a = windows[p];
                const RegionalWindow& b = windows[q];
                for (int y = 0; y < a.domainfineheight(); y++)
                    for (int x = 0; x < a.domainfinewidth(); x++)
                    {
                        std::size_t other = 0;
                        if (!domaincell(b, absolutecolumn(a, x), absoluterow(a, y), other))
                            continue;
                        const std::size_t cell =
                            static_cast<std::size_t>(y) * static_cast<std::size_t>(a.domainfinewidth()) + x;
                        const int bx = columnof(b, other);
                        const int by = rowof(b, other);
                        const FineChannels& ra = results[p];
                        const FineChannels& rb = results[q];
                        const bool same = ra.channel[cell] == rb.channel[other] &&
                                          ra.direction[cell] == rb.direction[other] && ra.bed_m[cell] == rb.bed_m[other] &&
                                          ra.discharge_m3s[cell] == rb.discharge_m3s[other] &&
                                          ra.januaryflow[cell] == rb.januaryflow[other] &&
                                          ra.julyflow[cell] == rb.julyflow[other] &&
                                          ra.globalcell[cell] == rb.globalcell[other];
                        if (!inside(a, x, y) || !inside(b, bx, by))
                        {
                            edgediffer += same ? 0 : 1;
                            continue;
                        }
                        compared++;
                        comparedchannels += ra.channel[cell];
                        ASSERT_TRUE(same) << p << " " << q << " fine cell " << cell;
                    }
            }
    }
    EXPECT_GT(compared, 0);
    EXPECT_GT(comparedchannels, 0);
    std::printf("[channels] R1 seams: %d interior cells compared (%d on channels), %d cells near an edge differ\n",
                compared, comparedchannels, edgediffer);
}

TEST_F(Landscape, ChannelNetworkCutBeforeALandReceiver)
{
    // At k = 8, (6, 4) at 100 m drains west to the land cell (5, 4) at 50 m, which drains
    // into the sea at (4, 4). A fine sea along the east column of (5, 4)'s block (fine
    // column 47, the column of (6, 4)'s entry point there) cuts (6, 4)'s link at its
    // entry point before its land receiver: the end becomes the sea level and the floor
    // the next float above it, so under a start surface 10 m below the sea level the
    // bed is the junction's 50 m (the receiver's filled elevation holds it), then the
    // next float above the sea level, then the sea level at column 47.
    GlobalNetwork network = seaworld();
    setland(network, 6, 4, west, 100.0f, 200.0f);
    setland(network, 5, 4, west, 50.0f, 300.0f);
    const int k = 8;
    const RegionalWindow window = worldwindow(k);
    const std::vector<float> start = finefield(window, [&](std::int64_t, std::int64_t) { return -10.0f; });
    const std::vector<uint8_t> blocks = seablocks(network, window);
    const std::vector<uint8_t> finesea = finemask(window, [&](std::int64_t x, std::int64_t y) {
        return blocks[static_cast<std::size_t>(y) * window.domainfinewidth() + static_cast<std::size_t>(x)] != 0 ||
               (y / k == 4 && x == 47);
    });
    const float abovesea = std::nextafter(sealevel, std::numeric_limits<float>::infinity());
    for (std::uint64_t seed = 1; seed <= 8; seed++)
    {
        SCOPED_TRACE(testing::Message() << "seed " << seed);
        const ChannelParameters settings = parameters(0.3, seed);
        const FineChannels channels = rasterisechannels(window, network, settings, sealevel, start, finesea);
        const std::vector<std::size_t> walk = walkfrom(window, channels, junctioncell(window, settings, 6, 4));
        ASSERT_GE(walk.size(), 2u);
        EXPECT_EQ(channels.bed_m[walk.front()], 50.0f);
        for (std::size_t m = 1; m + 1 < walk.size(); m++)
        {
            EXPECT_EQ(channels.bed_m[walk[m]], abovesea) << "step " << m;
            EXPECT_GE(absolutecolumn(window, columnof(window, walk[m])), 6 * k) << "step " << m;
        }
        EXPECT_EQ(absolutecolumn(window, columnof(window, walk.back())), 47);
        EXPECT_EQ(channels.bed_m[walk.back()], sealevel);
        EXPECT_EQ(channels.globalcell[walk.back()], static_cast<int32_t>(globalindex(6, 4)));
        // Nothing of (6, 4)'s link lies west of its entry point.
        for (std::size_t cell = 0; cell < channels.channel.size(); cell++)
            if (channels.globalcell[cell] == static_cast<int32_t>(globalindex(6, 4)))
                EXPECT_GE(absolutecolumn(window, columnof(window, cell)), 47);
        // The receiver's junction: the start surface held at the next float above the sea level.
        EXPECT_EQ(channels.bed_m[junctioncell(window, settings, 5, 4)], abovesea);
        checkfinecoast(network, window, settings, channels, finesea);
    }
}

TEST_F(Landscape, ChannelNetworkNaNInTheStartSurfaceLowersNothing)
{
    // A NaN in the start surface is skipped by the running minimum and by the junction
    // bed (no throw): the beds equal those with the NaN replaced by a value above every
    // link. (6, 4) at 100 m drains west to (5, 4) at 50 m, which drains into the sea;
    // k = 8, the start surface NaN (or 1000 m) on every other fine cell, 40 to 90 m elsewhere.
    GlobalNetwork network = seaworld();
    setland(network, 6, 4, west, 100.0f, 200.0f);
    setland(network, 5, 4, west, 50.0f, 300.0f);
    const int k = 8;
    const RegionalWindow window = worldwindow(k);
    auto surface = [&](float other) {
        return finefield(window, [&](std::int64_t x, std::int64_t y) {
            return (x + y) % 2 == 0 ? other : 40.0f + 50.0f * static_cast<float>(hashunit(x, y));
        });
    };
    for (std::uint64_t seed = 1; seed <= 8; seed++)
    {
        const ChannelParameters settings = parameters(0.3, seed);
        const FineChannels withnan =
            rasterisechannels(window, network, settings, sealevel, surface(std::numeric_limits<float>::quiet_NaN()));
        const FineChannels without = rasterisechannels(window, network, settings, sealevel, surface(1000.0f));
        EXPECT_EQ(withnan.channel, without.channel) << "seed " << seed;
        EXPECT_EQ(withnan.direction, without.direction) << "seed " << seed;
        EXPECT_EQ(withnan.globalcell, without.globalcell) << "seed " << seed;
        EXPECT_EQ(withnan.bed_m, without.bed_m) << "seed " << seed;
        // A spring's top comes from the finite start values only: no NaN bed.
        ChannelParameters withsprings = settings;
        withsprings.springspercell = 2;
        for (float bed : rasterisechannels(window, network, withsprings, sealevel,
                                           surface(std::numeric_limits<float>::quiet_NaN()))
                             .bed_m)
            EXPECT_FALSE(std::isnan(bed));
    }
}

// Multiscale plan S3.10: diagonal links through a transit block (channel_network.hpp).
namespace
{
// S3.10's geometry (the ChannelParameters defaults: junctions anywhere in their block, edge
// margin 0, diagonals through a transit block where one qualifies).
ChannelParameters transitparameters(double meanderamplitude, std::uint64_t seed, double pinarea_km2 = 0.0,
                                    int springs = 0)
{
    ChannelParameters result;
    result.pinarea_km2 = pinarea_km2;
    result.meanderamplitude = meanderamplitude;
    result.seed = seed;
    result.springspercell = springs;
    return result;
}

// The global cell (i + ox, j + oy), columns wrapped (the row must lie on the grid).
std::size_t offsetglobal(std::size_t cell, int ox, int oy)
{
    const int i = static_cast<int>(cell % worldwidth);
    const int j = static_cast<int>(cell / worldwidth);
    return globalindex((i + ox + worldwidth) % worldwidth, j + oy);
}

// Whether `neighbour` qualifies as the transit block of the diagonal link from `owner`
// (channel_network.hpp): land, filled between the receiver and the owner, neither the
// receiver's receiver nor a donor of the owner, draining to the receiver's outlet.
// With the default transittolerance_m a side cell missing that range by at most the tolerance
// may be chosen when neither lies inside it, so this admits a miss up to the tolerance (a link
// goes through its corner exactly when both side cells miss by more).
bool qualifiedtransit(const GlobalNetwork& network, std::size_t owner, std::size_t neighbour)
{
    const std::size_t receiver = globalreceiver(network, owner);
    const double filled = network.filledelevation_m[neighbour];
    const double miss = std::max({ 0.0, static_cast<double>(network.filledelevation_m[receiver]) - filled,
                                   filled - static_cast<double>(network.filledelevation_m[owner]) });
    return neighbour != owner && neighbour != receiver && network.sea[neighbour] == 0 &&
           miss <= ChannelParameters{}.transittolerance_m && globalreceiver(network, receiver) != neighbour &&
           globalreceiver(network, neighbour) != owner && globaloutlet(network, neighbour) == globaloutlet(network, receiver);
}

// The diagonal links of a network (land cells with a diagonal direction) and those with no
// qualifying transit block (drawn through their corner).
void diagonalcounts(const GlobalNetwork& network, int& diagonals, int& corners)
{
    diagonals = corners = 0;
    for (std::size_t cell = 0; cell < network.sea.size(); cell++)
    {
        const int direction = network.direction[cell];
        if (network.sea[cell] != 0 || direction == 0)
            continue;
        const int ox = landscape::directionoffsetx(direction);
        const int oy = landscape::directionoffsety(direction);
        if (ox == 0 || oy == 0)
            continue;
        diagonals++;
        corners += qualifiedtransit(network, cell, offsetglobal(cell, ox, 0)) ||
                           qualifiedtransit(network, cell, offsetglobal(cell, 0, oy))
                       ? 0
                       : 1;
    }
}

// Whether `block` may hold cells of `owner`'s link (or springs): the owner's block, its
// receiver's or, for a diagonal link, a qualifying transit block (`transit` set).
bool linkblock(const GlobalNetwork& network, std::size_t owner, std::size_t block, bool& transit)
{
    transit = false;
    const std::size_t receiver = globalreceiver(network, owner);
    if (block == owner || block == receiver)
        return true;
    const int direction = network.direction[owner];
    const int ox = landscape::directionoffsetx(direction);
    const int oy = landscape::directionoffsety(direction);
    if (receiver == owner || ox == 0 || oy == 0)
        return false;
    if (block != offsetglobal(owner, ox, 0) && block != offsetglobal(owner, 0, oy))
        return false;
    transit = qualifiedtransit(network, owner, block);
    return transit;
}

// Channel cells outside their link's blocks, cells in a transit block, and the steps along
// one owner's path from one block into another by an owner with cells in a transit block:
// through a corner (diagonal) or across an edge (4-connected).
struct TransitReport
{
    int outside = 0;
    int transitcells = 0;
    int cornersteps = 0;
    int edgesteps = 0;
};

TransitReport transitreport(const GlobalNetwork& network, const RegionalWindow& window, const FineChannels& channels)
{
    TransitReport report;
    std::set<int32_t> transitowners;
    for (std::size_t cell = 0; cell < channels.channel.size(); cell++)
    {
        if (channels.channel[cell] == 0)
            continue;
        bool transit = false;
        if (!linkblock(network, static_cast<std::size_t>(channels.globalcell[cell]), globalof(window, cell), transit))
            report.outside++;
        else if (transit)
        {
            report.transitcells++;
            transitowners.insert(channels.globalcell[cell]);
        }
    }
    for (std::size_t cell = 0; cell < channels.channel.size(); cell++)
    {
        std::size_t next = 0;
        if (channels.channel[cell] == 0 || channels.direction[cell] == 0 ||
            !nextcell(window, cell, channels.direction[cell], next) ||
            channels.globalcell[next] != channels.globalcell[cell] || globalof(window, next) == globalof(window, cell) ||
            transitowners.count(channels.globalcell[cell]) == 0)
            continue;
        const bool diagonal = channels.direction[cell] % 2 == 0;  // 2, 4, 6, 8: NE, SE, SW, NW
        (diagonal ? report.cornersteps : report.edgesteps)++;
    }
    return report;
}

// Walks the directions from every pinned global cell's junction (the window must be the
// whole world, no fine sea): the junction belongs to a piece of the cell's outlet (its own
// link, or a transit leg with a lower bed there: `takenjunctions`), every owner on the way
// drains to that outlet (`foreign` counts the others), and the walk ends in it (a sea
// block or a land sink's junction) without a loop (`wrongoutlet` counts the others).
struct OutletReport
{
    int walks = 0;
    int foreign = 0;
    int wrongoutlet = 0;
    int takenjunctions = 0;
    int dry = 0;  // S3.11: walks that step onto an unpainted cell (a dry tail), not counted as wrong
};

OutletReport walkoutlets(const GlobalNetwork& network, const RegionalWindow& window, const ChannelParameters& settings,
                         const FineChannels& channels)
{
    OutletReport report;
    const int limit = static_cast<int>(channels.channel.size());
    for (int j = 0; j < worldheight; j++)
        for (int i = 0; i < worldwidth; i++)
        {
            const std::size_t start = globalindex(i, j);
            if (network.sea[start] != 0 || network.direction[start] == 0 ||
                !(network.drainagearea_km2[start] >= settings.pinarea_km2))
                continue;
            report.walks++;
            const std::size_t outlet = globaloutlet(network, start);
            std::size_t walk = junctioncell(window, settings, i, j);
            bool ok = channels.channel[walk] == 1 &&
                      globaloutlet(network, static_cast<std::size_t>(channels.globalcell[walk])) == outlet;
            report.takenjunctions += ok && channels.globalcell[walk] != static_cast<int32_t>(start) ? 1 : 0;
            int steps = 0;
            bool dry = false;
            while (ok && channels.direction[walk] != 0 && steps <= limit)
            {
                std::size_t next = 0;
                ok = nextcell(window, walk, channels.direction[walk], next);
                dry = ok && channels.channel[next] == 0;  // S3.11: a dry tail leaves the pinned network
                if (!ok || dry)
                    break;
                if (globaloutlet(network, static_cast<std::size_t>(channels.globalcell[next])) != outlet)
                    report.foreign++;
                walk = next;
                steps++;
            }
            report.dry += dry ? 1 : 0;
            if (!dry && (!ok || steps > limit || globalof(window, walk) != outlet || !isend(network, window, settings, walk)))
                report.wrongoutlet++;
        }
    return report;
}

// Channel cells whose chain does not end at an end cell (a sea block or a land sink's
// junction) within the channel count's steps: a loop or a broken chain (no fine sea). A chain
// that steps onto an unpainted cell (a dry tail of the confluences, multiscale plan S3.11) ends
// there and is not counted.
int badchains(const GlobalNetwork& network, const RegionalWindow& window, const ChannelParameters& settings,
              const FineChannels& channels)
{
    const int count = channelcount(channels);
    int bad = 0;
    for (std::size_t cell = 0; cell < channels.channel.size(); cell++)
    {
        if (channels.channel[cell] == 0)
            continue;
        std::size_t walk = cell;
        int steps = 0;
        bool ok = true, dry = false;
        while (ok && !dry && channels.direction[walk] != 0 && steps <= count)
        {
            ok = nextcell(window, walk, channels.direction[walk], walk);
            dry = ok && channels.channel[walk] == 0;
            steps++;
        }
        bad += !ok || steps > count || (!dry && !isend(network, window, settings, walk)) ? 1 : 0;
    }
    return bad;
}

// A flat world: every land cell (columns 1 .. 15) filled at 50 m and draining to the previous
// column (west, north-west or south-west by a hash of the cell and `seed`), column 0 sea.
// Many diagonals, whose side cells drain to the same outlet or to a neighbouring one; on the
// flat every bed is 50 m without a start surface, so ownership goes by rank alone and
// transit legs cross other cells' junctions at equal beds.
GlobalNetwork flatnetwork(std::uint64_t seed)
{
    GlobalNetwork network = seaworld();
    for (int j = 0; j < worldheight; j++)
        for (int i = 1; i < worldwidth; i++)
        {
            std::uint64_t hash = static_cast<std::uint64_t>(i) * 73856093u ^ static_cast<std::uint64_t>(j) * 19349663u ^
                                 seed * 83492791u;
            hash ^= hash >> 13;
            hash *= 0x9E3779B97F4A7C15ull;
            hash ^= hash >> 29;
            int direction = hash % 3 == 0 ? west : (hash % 3 == 1 ? northwest : southwest);
            if ((j == 0 && direction == northwest) || (j == worldheight - 1 && direction == southwest))
                direction = west;
            setland(network, i, j, direction, 50.0f, 100.0f);
        }
    return network;
}

std::vector<float> hashedstart(const RegionalWindow& window, std::uint64_t seed)
{
    std::vector<float> start(static_cast<std::size_t>(window.domainfinewidth()) * window.domainfineheight());
    for (std::size_t cell = 0; cell < start.size(); cell++)
        start[cell] = static_cast<float>((cell * 2654435761u + seed) % 400) - 50.0f;
    return start;
}
}

TEST_F(Landscape, ChannelNetworkTransitBlocks)
{
    // A diagonal link's cells lie in its upstream cell's block, its receiver's and a transit
    // block that is land and drains to the receiver's outlet only (never a sea side cell or
    // one of another basin), and a link with cells in a transit block crosses each block
    // edge by a 4-connected step, never through a corner. The main world's 28
    // diagonal links go through a qualifying side cell where one exists, else through their
    // corner (with the default tolerance every one has one; with none, some keep their
    // corner); with diagonaltransit 0 all 28 go through their corner. With and without a start
    // surface and a speckled fine sea.
    const GlobalNetwork network = syntheticnetwork();
    int diagonals = 0, cornerlinks = 0;
    diagonalcounts(network, diagonals, cornerlinks);
    EXPECT_EQ(diagonals, 28);
    int transitcells = 0, edgesteps = 0;
    for (const Case& c : cases)
    {
        SCOPED_TRACE(testing::Message() << "k " << c.scale << ", meander " << c.meanderamplitude << ", seed " << c.seed
                                        << ", springs " << c.springs);
        const RegionalWindow window = worldwindow(c.scale);
        const ChannelParameters settings = transitparameters(c.meanderamplitude, c.seed, 0.0, c.springs);
        const std::vector<float> start = hashedstart(window, c.seed);
        std::vector<uint8_t> finesea(start.size());
        for (std::size_t cell = 0; cell < start.size(); cell++)
            finesea[cell] = network.sea[globalof(window, cell)] != 0 || (cell * 40503u) % 97 == 0 ? 1 : 0;
        for (const bool surfaces : { false, true })
        {
            const FineChannels channels = surfaces ? rasterisechannels(window, network, settings, sealevel, start, finesea)
                                                   : rasterisechannels(window, network, settings, sealevel);
            const TransitReport report = transitreport(network, window, channels);
            EXPECT_EQ(report.outside, 0) << (surfaces ? "with surfaces" : "");
            EXPECT_EQ(report.cornersteps, 0) << (surfaces ? "with surfaces" : "");
            EXPECT_EQ(risingsteps(window, channels), 0) << (surfaces ? "with surfaces" : "");
            EXPECT_EQ(channels.diagonallinks, diagonals);
            EXPECT_EQ(channels.cornerlinks, cornerlinks);
            transitcells += report.transitcells;
            edgesteps += report.edgesteps;
        }
    }
    EXPECT_GT(transitcells, 0);
    EXPECT_GT(edgesteps, 0);
    const FineChannels corners = rasterisechannels(worldwindow(8), network, parameters(0.3, 1), sealevel);
    EXPECT_EQ(corners.diagonallinks, 28);
    EXPECT_EQ(corners.cornerlinks, 28);
    ChannelParameters strict = transitparameters(0.3, 1);
    strict.transittolerance_m = 0.0;
    const FineChannels strictcorners = rasterisechannels(worldwindow(8), network, strict, sealevel);
    EXPECT_EQ(strictcorners.diagonallinks, 28);
    EXPECT_GT(strictcorners.cornerlinks, cornerlinks);
    EXPECT_LT(strictcorners.cornerlinks, 28);
}

TEST_F(Landscape, ChannelNetworkTransitCornerFallback)
{
    // A diagonal with no qualifying side cell keeps the shared corner: on the X of
    // ChannelNetworkCrossingDiagonalsOnAFlat each diagonal's side cells drain to the other
    // basin's outlet, so both links are drawn exactly as with diagonaltransit 0 (the same
    // junctions and edge margins) and counted as corner links. Along a coast (the network of
    // ChannelNetworkDiagonalsAlongACoastStayOnLand) a sea side cell is never a transit block:
    // (8, 3) -> (9, 4) crosses its land side cell (8, 4), which drains into the receiver;
    // (5, 3) -> (6, 4) keeps its corner, since its land side cell (5, 4) is the receiver's
    // receiver and lies below it; no channel cell lies in the sea side blocks.
    GlobalNetwork cross = seaworld();
    setland(cross, 6, 3, southeast, 50.0f, 100.0f);
    setland(cross, 7, 4, east, 50.0f, 200.0f);
    setland(cross, 8, 4, east, 50.0f, 300.0f);
    setland(cross, 7, 3, southwest, 50.0f, 100.0f);
    setland(cross, 6, 4, west, 50.0f, 200.0f);
    setland(cross, 5, 4, west, 50.0f, 300.0f);
    GlobalNetwork coast = seaworld();
    setland(coast, 5, 4, west, 10.0f, 300.0f);
    setland(coast, 6, 4, west, 20.0f, 200.0f);
    setland(coast, 5, 3, southeast, 40.0f, 100.0f);
    setland(coast, 8, 4, east, 10.0f, 100.0f);
    setland(coast, 8, 3, southeast, 30.0f, 100.0f);
    int coasttransit = 0;
    for (int k : { 1, 2, 4, 5, 8, 16 })
        for (std::uint64_t seed = 1; seed <= 16; seed++)
        {
            SCOPED_TRACE(testing::Message() << "k " << k << ", seed " << seed);
            const RegionalWindow window = worldwindow(k);
            const ChannelParameters on = transitparameters(0.3, seed, 0.0, static_cast<int>(seed % 3));
            ChannelParameters off = on;
            off.diagonaltransit = 0;
            const FineChannels a = rasterisechannels(window, cross, on, sealevel);
            const FineChannels b = rasterisechannels(window, cross, off, sealevel);
            EXPECT_EQ(a.channel, b.channel);
            EXPECT_EQ(a.direction, b.direction);
            EXPECT_EQ(a.bed_m, b.bed_m);
            EXPECT_EQ(a.globalcell, b.globalcell);
            EXPECT_EQ(a.diagonallinks, 2);
            EXPECT_EQ(a.cornerlinks, 2);
            const OutletReport crossreport = walkoutlets(cross, window, on, a);
            EXPECT_EQ(crossreport.walks, 6);
            EXPECT_EQ(crossreport.wrongoutlet, 0);

            const FineChannels shore = rasterisechannels(window, coast, on, sealevel);
            EXPECT_EQ(shore.diagonallinks, 2);
            EXPECT_EQ(shore.cornerlinks, 1);
            for (std::size_t cell = 0; cell < shore.channel.size(); cell++)
            {
                if (shore.channel[cell] == 0)
                    continue;
                const std::size_t block = globalof(window, cell);
                ASSERT_TRUE(block != globalindex(6, 3) && block != globalindex(9, 3)) << "fine cell " << cell;
                coasttransit += (block == globalindex(5, 4) && shore.globalcell[cell] == static_cast<int32_t>(globalindex(5, 3))) ||
                                        (block == globalindex(8, 4) &&
                                         shore.globalcell[cell] == static_cast<int32_t>(globalindex(8, 3)))
                                    ? 1
                                    : 0;
            }
            const OutletReport shorereport = walkoutlets(coast, window, on, shore);
            EXPECT_EQ(shorereport.walks, 5);
            EXPECT_EQ(shorereport.foreign, 0);
            EXPECT_EQ(shorereport.wrongoutlet, 0);
        }
    EXPECT_GT(coasttransit, 0);
}

TEST_F(Landscape, ChannelNetworkTransitNoLoop)
{
    // With transit legs the bed never rises along the directions, no direction loops, every
    // chain ends at a sea block or a land sink's junction, and every owner a walk from a
    // junction meets drains to that junction's outlet (no capture across a divide): on the
    // main world and on flat worlds of many diagonals, where every bed is equal without a
    // start surface and transit legs take other cells' junctions by rank (the case a
    // junction-first rule would turn into a loop), with and without a start surface.
    int taken = 0, corners = 0, diagonals = 0;
    for (std::uint64_t world = 0; world <= 4; world++)
    {
        const GlobalNetwork network = world == 0 ? syntheticnetwork() : flatnetwork(world);
        for (const Case& c : cases)
            for (const bool withstart : { false, true })
            {
                SCOPED_TRACE(testing::Message() << "world " << world << ", k " << c.scale << ", meander "
                                                << c.meanderamplitude << ", seed " << c.seed << ", springs " << c.springs
                                                << (withstart ? ", start" : ""));
                const RegionalWindow window = worldwindow(c.scale);
                const ChannelParameters settings = transitparameters(c.meanderamplitude, c.seed, 0.0, c.springs);
                const FineChannels channels =
                    withstart ? rasterisechannels(window, network, settings, sealevel, hashedstart(window, c.seed + world))
                              : rasterisechannels(window, network, settings, sealevel);
                ASSERT_EQ(risingsteps(window, channels), 0);
                ASSERT_EQ(badchains(network, window, settings, channels), 0);
                const OutletReport report = walkoutlets(network, window, settings, channels);
                ASSERT_EQ(report.walks, 120);
                ASSERT_EQ(report.foreign, 0);
                ASSERT_EQ(report.wrongoutlet, 0);
                const TransitReport blocks = transitreport(network, window, channels);
                ASSERT_EQ(blocks.outside, 0);
                ASSERT_EQ(blocks.cornersteps, 0);
                taken += report.takenjunctions;
                corners += channels.cornerlinks;
                diagonals += channels.diagonallinks;
            }
    }
    std::printf("[channels] transit worlds: %d diagonal links, %d through a corner, %d junctions taken by a transit\n",
                diagonals, corners, taken);
    EXPECT_GT(taken, 0);
    EXPECT_GT(corners, 0);
    EXPECT_LT(corners, diagonals);
}

TEST_F(Landscape, ChannelNetworkTransitWindowsAgree)
{
    // Every window draws the links within one global cell of its domain, and a transit link
    // lies in the 2 x 2 blocks around its corner: windows agree with the whole world on
    // every fine cell, including windows whose edge column or row is a transit block of a
    // link whose upstream cell lies outside them (counted: `edgetransit`).
    int edgetransit = 0;
    for (std::uint64_t world = 0; world <= 2; world++)
    {
        const GlobalNetwork network = world == 0 ? syntheticnetwork() : flatnetwork(world);
        for (const Case& c : cases)
        {
            SCOPED_TRACE(testing::Message() << "world " << world << ", k " << c.scale << ", meander "
                                            << c.meanderamplitude << ", seed " << c.seed << ", springs " << c.springs);
            const ChannelParameters settings = transitparameters(c.meanderamplitude, c.seed, 0.0, c.springs);
            const RegionalWindow whole = worldwindow(c.scale);
            const FineChannels reference = rasterisechannels(whole, network, settings, sealevel);
            const RegionalWindow windows[] = {
                domainwindow(c.scale, 2, 4, 4, 2),   // row 4's diagonal from (1, 4) through (2, 4)
                domainwindow(c.scale, 3, 6, 1, 3),   // rows 1 .. 3: row 1's diagonals into row 2
                domainwindow(c.scale, 4, 5, 5, 3),   // rows 5 .. 7: row 6's diagonals, row 4 outside
                domainwindow(c.scale, 13, 8, 1, 5),  // across the date line
                domainwindow(c.scale, 9, 5, 3, 4),
            };
            for (const RegionalWindow& window : windows)
            {
                const FineChannels channels = rasterisechannels(window, network, settings, sealevel);
                for (int y = 0; y < window.domainfineheight(); y++)
                    for (int x = 0; x < window.domainfinewidth(); x++)
                    {
                        std::size_t other = 0;
                        ASSERT_TRUE(domaincell(whole, absolutecolumn(window, x), absoluterow(window, y), other));
                        const std::size_t cell =
                            static_cast<std::size_t>(y) * static_cast<std::size_t>(window.domainfinewidth()) + x;
                        ASSERT_EQ(channels.channel[cell], reference.channel[other]) << cell;
                        ASSERT_EQ(channels.direction[cell], reference.direction[other]) << cell;
                        ASSERT_EQ(channels.bed_m[cell], reference.bed_m[other]) << cell;
                        ASSERT_EQ(channels.discharge_m3s[cell], reference.discharge_m3s[other]) << cell;
                        ASSERT_EQ(channels.globalcell[cell], reference.globalcell[other]) << cell;
                        if (channels.channel[cell] == 0)
                            continue;
                        const std::size_t owner = static_cast<std::size_t>(channels.globalcell[cell]);
                        const int ownerrow = static_cast<int>(owner / worldwidth);
                        const int ownercolumn =
                            ((static_cast<int>(owner % worldwidth) - window.domaincolumn0) % worldwidth + worldwidth) %
                            worldwidth;
                        const bool outside = ownerrow < window.domainrow0 ||
                                             ownerrow >= window.domainrow0 + window.domainrows ||
                                             ownercolumn >= window.domaincolumns;
                        bool transit = false;
                        if (outside && linkblock(network, owner, globalof(window, cell), transit) && transit)
                            edgetransit++;
                    }
            }
        }
    }
    EXPECT_GT(edgetransit, 0);
}

TEST_F(Landscape, ChannelNetworkTransitWindowsAgreeWithSurfaces)
{
    // ChannelNetworkWindowsAgreeAwayFromTheirEdges with transit legs: given the same start
    // surface and fine sea (functions of the absolute fine cell, sea specks cutting links
    // anywhere), windows agree on every channel field at the fine cells at least two global
    // cells inside both domains (one until S3.10): a transit link lies in the 2 x 2 blocks around
    // its corner and its junction beds, running minimum and fine-sea cut read only cells there,
    // and with the confluences (S3.11) a piece's prefix end reads the owners along its path, whose
    // beds read the start surface one global cell beyond (channel_network.hpp, Seams); including
    // transit cells whose link starts outside the smaller window (counted).
    int compared = 0, comparedtransit = 0;
    for (std::uint64_t world = 0; world <= 2; world++)
    {
        const GlobalNetwork network = world == 0 ? syntheticnetwork() : flatnetwork(world);
        for (const Case& c : cases)
        {
            SCOPED_TRACE(testing::Message() << "world " << world << ", k " << c.scale << ", meander "
                                            << c.meanderamplitude << ", seed " << c.seed);
            const int k = c.scale;
            const ChannelParameters settings = transitparameters(c.meanderamplitude, c.seed, 0.0, c.springs);
            const RegionalWindow windows[] = {
                worldwindow(k),
                domainwindow(k, 1, 6, 3, 4),
                domainwindow(k, 3, 7, 0, 4),
                domainwindow(k, 13, 8, 1, 5),
                domainwindow(k, 9, 5, 3, 4),
            };
            std::vector<FineChannels> results;
            for (const RegionalWindow& window : windows)
            {
                const std::vector<uint8_t> finesea = finemask(window, [&](std::int64_t x, std::int64_t y) {
                    return x < k + k / 2 || hashunit(x, y, 77) < 0.03;
                });
                results.push_back(rasterisechannels(window, network, settings, sealevel,
                                                    noisesurface(window, -50.0f, 400.0f), finesea));
            }
            const RegionalWindow& whole = windows[0];
            for (std::size_t q = 1; q < results.size(); q++)
            {
                const RegionalWindow& window = windows[q];
                for (int y = 2 * k; y < window.domainfineheight() - 2 * k; y++)
                    for (int x = 2 * k; x < window.domainfinewidth() - 2 * k; x++)
                    {
                        std::size_t other = 0;
                        ASSERT_TRUE(domaincell(whole, absolutecolumn(window, x), absoluterow(window, y), other));
                        const std::size_t cell =
                            static_cast<std::size_t>(y) * static_cast<std::size_t>(window.domainfinewidth()) + x;
                        const FineChannels& a = results[q];
                        const FineChannels& b = results[0];
                        compared++;
                        ASSERT_EQ(a.channel[cell], b.channel[other]) << q << " fine cell " << cell;
                        ASSERT_EQ(a.direction[cell], b.direction[other]) << q << " fine cell " << cell;
                        ASSERT_EQ(a.bed_m[cell], b.bed_m[other]) << q << " fine cell " << cell;
                        ASSERT_EQ(a.discharge_m3s[cell], b.discharge_m3s[other]) << q << " fine cell " << cell;
                        ASSERT_EQ(a.globalcell[cell], b.globalcell[other]) << q << " fine cell " << cell;
                        if (a.channel[cell] == 0)
                            continue;
                        bool transit = false;
                        if (linkblock(network, static_cast<std::size_t>(a.globalcell[cell]), globalof(window, cell),
                                      transit) &&
                            transit)
                            comparedtransit++;
                    }
            }
        }
    }
    EXPECT_GT(compared, 0);
    EXPECT_GT(comparedtransit, 0);
}
// Multiscale plan S3.11: confluences and leg roughness (channel_network.hpp).
namespace
{
// What the confluence rule leaves visible on a whole-world window without springs.
struct ConfluenceReport
{
    int openprefixes = 0;  // links whose walk from their junction, while they own the cells, steps onto an
                           // unpainted cell: a prefix that does not reach its first lost cell
    int orphans = 0;       // painted cells, not a land cell's junction, that no painted 8-neighbour steps onto
                           // (a suffix that runs on past a cell another piece owns)
    int loops = 0;         // channel cells whose walk along the directions does not stop
    int drysteps = 0;      // channel cells whose direction leads to an unpainted cell (dry tails)
};

ConfluenceReport confluencereport(const GlobalNetwork& network, const RegionalWindow& window,
                                  const ChannelParameters& settings, const FineChannels& channels)
{
    ConfluenceReport report;
    std::vector<uint8_t> entered(channels.channel.size(), 0), junction(channels.channel.size(), 0);
    for (int j = 0; j < worldheight; j++)
        for (int i = 0; i < worldwidth; i++)
            junction[junctioncell(window, settings, i, j)] = 1;
    const int count = channelcount(channels);
    for (std::size_t cell = 0; cell < channels.channel.size(); cell++)
    {
        std::size_t next = 0;
        if (channels.channel[cell] == 0 || channels.direction[cell] == 0 ||
            !nextcell(window, cell, channels.direction[cell], next))
            continue;
        if (channels.channel[next] != 0)
            entered[next] = 1;
        else
            report.drysteps++;
        std::size_t walk = cell;
        int steps = 0;
        while (channels.channel[walk] != 0 && channels.direction[walk] != 0 && steps <= count &&
               nextcell(window, walk, channels.direction[walk], walk))
            steps++;
        report.loops += steps > count ? 1 : 0;
    }
    for (std::size_t cell = 0; cell < channels.channel.size(); cell++)
        report.orphans += channels.channel[cell] != 0 && entered[cell] == 0 && junction[cell] == 0 ? 1 : 0;
    for (int j = 0; j < worldheight; j++)
        for (int i = 0; i < worldwidth; i++)
        {
            const std::size_t start = globalindex(i, j);
            if (network.sea[start] != 0 || network.direction[start] == 0 ||
                !(network.drainagearea_km2[start] >= settings.pinarea_km2))
                continue;
            std::size_t walk = junctioncell(window, settings, i, j);
            int steps = 0;
            while (channels.channel[walk] != 0 && channels.globalcell[walk] == static_cast<int32_t>(start) &&
                   channels.direction[walk] != 0 && steps <= count &&
                   nextcell(window, walk, channels.direction[walk], walk))
                steps++;
            report.openprefixes += channels.channel[walk] == 0 ? 1 : 0;
        }
    return report;
}
}

TEST_F(Landscape, ChannelNetworkConfluences)
{
    // A piece is painted on its prefix (to its first lost cell, included) and its suffix (from
    // another piece's first lost cell after its prefix): every link's walk from its junction
    // reaches its first lost cell painted, the beds never rise along the directions and no
    // direction loops; every owner a walk meets drains to the walk's outlet, and links stay in
    // their blocks; on the main world (R3's geometry and S3.10's) and on flat worlds of many
    // diagonals, with and without a start surface. Pieces do merge (mergedpieces); chains leaving
    // the pinned network (drytails) and painted cells nothing steps onto (a suffix that runs on past
    // a cell another piece owns) are counted, not asserted.
    int merged = 0, dry = 0, drysteps = 0, orphans = 0, runs = 0;
    for (std::uint64_t world = 0; world <= 2; world++)
    {
        const GlobalNetwork network = world == 0 ? syntheticnetwork() : flatnetwork(world);
        for (const Case& c : cases)
            for (const bool transit : { false, true })
                for (const bool withstart : { false, true })
                {
                    if (world != 0 && !transit)
                        continue;
                    SCOPED_TRACE(testing::Message() << "world " << world << ", k " << c.scale << ", meander "
                                                    << c.meanderamplitude << ", seed " << c.seed
                                                    << (transit ? ", transit" : ", R3") << (withstart ? ", start" : ""));
                    const RegionalWindow window = worldwindow(c.scale);
                    const ChannelParameters settings =
                        transit ? transitparameters(c.meanderamplitude, c.seed) : parameters(c.meanderamplitude, c.seed);
                    const FineChannels channels =
                        withstart ? rasterisechannels(window, network, settings, sealevel, hashedstart(window, c.seed + world))
                                  : rasterisechannels(window, network, settings, sealevel);
                    const ConfluenceReport report = confluencereport(network, window, settings, channels);
                    ASSERT_EQ(report.openprefixes, 0);
                    ASSERT_EQ(report.loops, 0);
                    ASSERT_EQ(risingsteps(window, channels), 0);
                    ASSERT_EQ(walkoutlets(network, window, settings, channels).foreign, 0);
                    ASSERT_EQ(transitreport(network, window, channels).outside, 0);
                    merged += channels.mergedpieces;
                    dry += channels.drytails;
                    drysteps += report.drysteps;
                    orphans += report.orphans;
                    runs++;
                }
    }
    std::printf("[channels] confluences: %d pieces merged, %d dry tails (%d dry steps), %d unentered suffix cells in %d "
                "rasterisations\n",
                merged, dry, drysteps, orphans, runs);
    EXPECT_GT(merged, 0);
}

TEST_F(Landscape, ChannelNetworkRoughLegs)
{
    // Leg roughness at k 32 and 64 with a large roughness (0.2, persistence 0.7): the confluence
    // report stays clean, the beds never rise, every channel cell lies in its link's blocks with
    // 4-connected block crossings, and the legs move (the owners differ from the smooth arcs').
    // Overlapping windows agree on every shared cell.
    for (std::uint64_t world = 0; world <= 1; world++)
    {
        const GlobalNetwork network = world == 0 ? syntheticnetwork() : flatnetwork(1);
        for (int k : { 32, 64 })
            for (std::uint64_t seed = 1; seed <= 2; seed++)
            {
                SCOPED_TRACE(testing::Message() << "world " << world << ", k " << k << ", seed " << seed);
                const RegionalWindow window = worldwindow(k);
                ChannelParameters rough = transitparameters(0.2, seed);
                rough.meanderroughness = 0.2;
                rough.meanderroughnesspersistence = 0.7;
                ChannelParameters smooth = rough;
                smooth.meanderroughness = 0.0;
                const FineChannels channels = rasterisechannels(window, network, rough, sealevel);
                const ConfluenceReport report = confluencereport(network, window, rough, channels);
                ASSERT_EQ(report.openprefixes, 0);
                ASSERT_EQ(report.loops, 0);
                ASSERT_EQ(risingsteps(window, channels), 0);
                ASSERT_EQ(walkoutlets(network, window, rough, channels).foreign, 0);
                const TransitReport blocks = transitreport(network, window, channels);
                ASSERT_EQ(blocks.outside, 0);
                ASSERT_EQ(blocks.cornersteps, 0);
                EXPECT_NE(channels.globalcell, rasterisechannels(window, network, smooth, sealevel).globalcell);
                const RegionalWindow part = domainwindow(k, 3, 6, 1, 3);
                const FineChannels windowed = rasterisechannels(part, network, rough, sealevel);
                for (int y = 0; y < part.domainfineheight(); y++)
                    for (int x = 0; x < part.domainfinewidth(); x++)
                    {
                        std::size_t other = 0;
                        ASSERT_TRUE(domaincell(window, absolutecolumn(part, x), absoluterow(part, y), other));
                        const std::size_t cell =
                            static_cast<std::size_t>(y) * static_cast<std::size_t>(part.domainfinewidth()) + x;
                        ASSERT_EQ(windowed.channel[cell], channels.channel[other]) << cell;
                        ASSERT_EQ(windowed.globalcell[cell], channels.globalcell[other]) << cell;
                        ASSERT_EQ(windowed.direction[cell], channels.direction[other]) << cell;
                        ASSERT_EQ(windowed.bed_m[cell], channels.bed_m[other]) << cell;
                    }
            }
    }
}

TEST_F(Landscape, ChannelNetworkRoughnessOffKeepsTheArcs)
{
    // meanderroughness 0 draws S3.10's smooth arcs: the walk of a single link (no other piece, so
    // no confluence) from (6, 4) west, or (6, 3) south-west through the corner, into a land sink,
    // with S3.10's geometry, pinned as its length and an FNV-1a hash of its absolute fine cells.
    // The values were printed by this code and match S3.10 because, with the confluences switched
    // off by a temporary hook, roughness 0 reproduced the S3.10 binary's regional bundle of
    // window B at k 16 (erosion run 9) file for file (2026-10-11). Roughness 0.05 moves the k 64 legs.
    GlobalNetwork cardinal = seaworld();
    setland(cardinal, 6, 4, west, 100.0f, 200.0f);
    setland(cardinal, 5, 4, 0, 50.0f, 300.0f);
    GlobalNetwork diagonal = seaworld();
    setland(diagonal, 6, 3, southwest, 100.0f, 200.0f);
    setland(diagonal, 5, 4, 0, 50.0f, 300.0f);
    setland(diagonal, 5, 3, 0, 80.0f, 100.0f);
    struct Pinned
    {
        int world;
        int k;
        std::uint64_t seed;
        std::size_t cells;
        std::uint64_t hash;
    };
    const Pinned pinned[] = {
        { 0, 4, 1, 4, 11359039419027424115ull },  { 0, 4, 2, 7, 1674562202105952313ull },
        { 0, 4, 3, 4, 3720837384589586742ull },   { 0, 16, 1, 15, 12990339599241022070ull },
        { 0, 16, 2, 23, 13287087815768985187ull }, { 0, 16, 3, 13, 15329570031918160423ull },
        { 0, 64, 1, 68, 13916128713388234474ull }, { 0, 64, 2, 91, 18007080842466672190ull },
        { 0, 64, 3, 52, 2286038539570859934ull },  { 1, 4, 1, 6, 10692345346872103989ull },
        { 1, 4, 2, 7, 1888782723492141791ull },    { 1, 4, 3, 6, 9773813615215116780ull },
        { 1, 16, 1, 22, 2647325632442881105ull },  { 1, 16, 2, 26, 6896447059521687915ull },
        { 1, 16, 3, 24, 4645429708154544421ull },  { 1, 64, 1, 88, 606380696182006003ull },
        { 1, 64, 2, 106, 5046840530215962682ull }, { 1, 64, 3, 95, 1148829669704141986ull },
    };
    int moved = 0;
    for (const Pinned& p : pinned)
    {
        SCOPED_TRACE(testing::Message() << "world " << p.world << ", k " << p.k << ", seed " << p.seed);
        const RegionalWindow window = worldwindow(p.k);
        ChannelParameters settings = transitparameters(0.2, p.seed);
        auto walkhash = [&](double roughness, std::size_t& cells) {
            settings.meanderroughness = roughness;
            const FineChannels channels =
                rasterisechannels(window, p.world == 0 ? cardinal : diagonal, settings, sealevel);
            const std::vector<std::size_t> walk =
                walkfrom(window, channels, junctioncell(window, settings, 6, p.world == 0 ? 4 : 3));
            std::uint64_t hash = 1469598103934665603ull;
            for (const std::size_t cell : walk)
            {
                hash = (hash ^ static_cast<std::uint64_t>(absolutecolumn(window, columnof(window, cell)))) *
                       1099511628211ull;
                hash = (hash ^ static_cast<std::uint64_t>(absoluterow(window, rowof(window, cell)))) * 1099511628211ull;
            }
            cells = walk.size();
            return hash;
        };
        std::size_t cells = 0;
        EXPECT_EQ(walkhash(0.0, cells), p.hash);
        EXPECT_EQ(cells, p.cells);
        if (p.k == 64)
            moved += walkhash(0.05, cells) != p.hash ? 1 : 0;
    }
    EXPECT_GT(moved, 0);
}
