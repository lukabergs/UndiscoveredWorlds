#include "simulations/landscape/channel_network.hpp"

// Pinned channels of the regional pass (erosion plan step E2b,
// docs/architecture/erosion-hydrology-plan.md, Phase E and section 6, "Seams of the
// solve"). Every quantity below is computed in absolute fine coordinates from one
// global link (its two global cells and the seed) and only mapped to the window's
// domain when a fine cell is written, so a link is the same fine path in every
// window, and a fine cell's values depend only on the links whose paths cross it.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

#include "core/deterministic_random.hpp"
#include "simulations/landscape/sphere_metric.hpp"

namespace landscape
{
namespace
{
constexpr double pi = 3.14159265358979323846;
// Largest value of sin(pi t) sin(2 pi t) = 2 sin^2(pi t) cos(pi t) on 0 .. 1, at
// cos(pi t) = 1 / sqrt(3): 4 / (3 sqrt(3)).
constexpr double sshapemaximum = 0.76980035891950105;
// Curve samples per fine cell of link length. The curve's speed is at most
// L (1 + 2 pi) fine cells per unit t (the bump's slope is bounded by 2 pi D and
// D < L), so consecutive samples lie less than one fine cell apart and fall in the
// same or 8-neighbouring fine cells.
constexpr double samplesperfinecell = 16.0;
// Global cells around the domain whose outgoing links are drawn. A path stays
// within its two global cells' blocks widened by its largest sideways displacement,
// which is below half a global cell (see linkpath), so only a link with an end in a
// global cell next to the domain can reach it; its upstream cell is then at most
// two global cells away.
constexpr int linkreach = 2;

// An unwrapped absolute fine cell: column k i .. on the planet-wide fine grid
// (outside 0 .. kW - 1 for a link across the date line; wrapped when written) and
// absolute fine row.
struct FineCell
{
    std::int64_t x = 0;
    std::int64_t y = 0;
};

bool samecell(const FineCell& a, const FineCell& b) { return a.x == b.x && a.y == b.y; }

bool neighbours8(const FineCell& a, const FineCell& b)
{
    return !samecell(a, b) && std::llabs(a.x - b.x) <= 1 && std::llabs(a.y - b.y) <= 1;
}

int directioncode(std::int64_t ox, std::int64_t oy)
{
    for (int direction = 1; direction <= directioncount; direction++)
        if (directionoffsetx(direction) == ox && directionoffsety(direction) == oy)
            return direction;
    return 0;
}

// The largest sideways displacement D of any link at scale k, in fine cells:
// 0.35 (k - 2), at least 0 (a link is straight at k <= 2). Derivation, with |d| <= D
// along every curve (the bound linkpath's shape gives):
// - Two links not adjacent on the global tree share no end (a cell has one outgoing
//   link), so their straight segments either cross (two diagonals of one 2 x 2 block
//   of global cells, an X) or lie at least k / sqrt(2) apart: the nearest such
//   segments are a lattice point (a link's end, centre to centre the lattice has
//   spacing k) and a diagonal through the middle of a neighbouring block side
//   (parallel diagonals, or a link ending beside a diagonal's middle), k / sqrt(2);
//   a lattice point off a cardinal segment's line lies at least k from it, and one
//   on the line of a segment that does not end there lies at least k from it.
// - Two curves each within D of their segments are then at least k / sqrt(2) - 2 D
//   apart, and two fine cells (unit squares) can be the same only for points less
//   than sqrt(2) apart. D <= (k - 2) / (2 sqrt(2)) = 0.3536 (k - 2) keeps the curves
//   sqrt(2) apart, so non-adjacent links share no fine cell and one never captures
//   the other's chain; 0.35 (k - 2) keeps that with a margin against rounding. The
//   same bound serves cardinal and diagonal links, because a cardinal link's end
//   can lie k / sqrt(2) from a diagonal's middle; only the X is left (header).
// - A global cell's centre middle point lies at least k / sqrt(2) from any link's
//   segment that does not end at it, and k / sqrt(2) - 0.35 (k - 2) > sqrt(2) / 2, so
//   no path enters another cell's centre fine cell.
// - D < k / 2: a path stays within half a global cell of its two blocks (linkreach).
// The derivation is for k >= 3. For k <= 2 the cap is 0 and every link is straight:
// at k = 2 non-adjacent segments lie at least sqrt(2) apart, and two points sqrt(2)
// apart never fall in one half-open fine cell; at k = 1 a link is its two centre
// cells, which non-adjacent links do not share.
double meandercap(int k) { return k > 2 ? 0.35 * (k - 2) : 0.0; }

// A value in [-1, 1) from the top 53 bits of a hash.
double signedunit(std::uint64_t hash)
{
    return static_cast<double>(hash >> 11) * (2.0 / 9007199254740992.0) - 1.0;
}

// The fine path of the link from global cell (i, j) to its neighbour (i + ox, j + oy)
// (unwrapped column), as an 8-connected chain of distinct fine cells, each step to
// an 8-neighbour, from the upstream centre fine cell to the receiver's.
//
// Centre fine cell of global (i, j): absolute fine (k i + k / 2, k j + k / 2), integer
// division (the block's middle cell for odd k, the cell south-east of the block's
// middle point for even k). The curve runs between the two centre cells' middle
// points, displaced sideways along the segment's left normal by
// d(t) = D (a1 sin^2(pi t) + a2 sin(pi t) sin(2 pi t) / sshapemaximum) / 2, |d| <= D:
// one bump and one S, both zero with a zero slope at the ends, a1 and a2 in [-1, 1)
// from the link's hash. D = min(meanderamplitude L, meandercap(k)) with L the segment
// length in fine cells (see meandercap for what the cap guarantees). The curve is
// sampled every
// 1 / samplesperfinecell of a fine cell or finer; each sample's fine cell is
// appended to the chain, a cell already on the chain cuts the chain back to it (the
// loop is erased), and a cell that neighbours the cell before the chain's last
// replaces that last cell (an L-shaped corner becomes a diagonal step). The
// arithmetic is in the upstream block's own coordinates, so it does not depend on
// the cell's absolute position beyond the integer offset added at the end.
std::vector<FineCell> linkpath(int k, int i, int j, int ox, int oy, double meanderamplitude, std::uint64_t hash)
{
    const double centre = static_cast<double>(k / 2) + 0.5;
    const double vx = static_cast<double>(ox) * k;
    const double vy = static_cast<double>(oy) * k;
    const double length = std::hypot(vx, vy);
    const double displacement = std::min(meanderamplitude * length, meandercap(k));
    const double nx = -vy / length;
    const double ny = vx / length;
    const double a1 = signedunit(hash);
    const double a2 = signedunit(deterministicmix64(hash));

    const FineCell first{ k / 2, k / 2 };
    const FineCell last{ k / 2 + static_cast<std::int64_t>(ox) * k, k / 2 + static_cast<std::int64_t>(oy) * k };
    const int samples = static_cast<int>(std::ceil(samplesperfinecell * length));

    std::vector<FineCell> chain{ first };
    for (int s = 1; s <= samples; s++)
    {
        FineCell cell = last;
        if (s < samples)
        {
            const double t = static_cast<double>(s) / samples;
            const double sine = std::sin(pi * t);
            const double shape = 0.5 * (a1 * sine * sine + a2 * sine * std::sin(2.0 * pi * t) / sshapemaximum);
            const double d = displacement * shape;
            cell.x = static_cast<std::int64_t>(std::floor(centre + t * vx + d * nx));
            cell.y = static_cast<std::int64_t>(std::floor(centre + t * vy + d * ny));
        }
        if (samecell(cell, chain.back()))
            continue;
        const auto found =
            std::find_if(chain.begin(), chain.end(), [&](const FineCell& other) { return samecell(other, cell); });
        if (found != chain.end())
        {
            chain.erase(found + 1, chain.end());
            continue;
        }
        while (chain.size() >= 2 && neighbours8(chain[chain.size() - 2], cell))
            chain.pop_back();
        chain.push_back(cell);
    }

    const std::int64_t originx = static_cast<std::int64_t>(k) * i;
    const std::int64_t originy = static_cast<std::int64_t>(k) * j;
    for (FineCell& cell : chain)
    {
        cell.x += originx;
        cell.y += originy;
    }
    return chain;
}

void requiresize(std::size_t size, std::size_t expected, const char* name)
{
    if (size != expected)
        throw std::invalid_argument(std::string("rasterisechannels: GlobalNetwork.") + name +
                                    " does not have width x height values");
}

// The receiver of a global cell on the drainage tree, or the cell itself for an
// outlet: a sea cell, a land cell without a direction (a sink), or one whose
// direction points beyond a pole.
std::size_t globalreceiver(const GlobalNetwork& network, std::size_t cell)
{
    const int direction = network.direction[cell];
    if (network.sea[cell] != 0 || direction < 1 || direction > directioncount)
        return cell;
    const int i = static_cast<int>(cell % static_cast<std::size_t>(network.width));
    const int j = static_cast<int>(cell / static_cast<std::size_t>(network.width)) + directionoffsety(direction);
    if (j < 0 || j >= network.height)
        return cell;
    return static_cast<std::size_t>(j) * static_cast<std::size_t>(network.width) +
           static_cast<std::size_t>((i + directionoffsetx(direction) + network.width) % network.width);
}

// The number of links from a global cell down the tree to its outlet (0 for an
// outlet), memoised in `distance` (-1 unknown). Receivers that form a loop throw.
int32_t outletdistance(const GlobalNetwork& network, std::size_t cell, std::vector<int32_t>& distance)
{
    std::vector<std::size_t> pending;
    std::size_t current = cell;
    while (distance[current] < 0)
    {
        const std::size_t next = globalreceiver(network, current);
        if (next == current)
        {
            distance[current] = 0;
            break;
        }
        pending.push_back(current);
        if (pending.size() > distance.size())
            throw std::invalid_argument("rasterisechannels: the global receivers form a loop");
        current = next;
    }
    int32_t value = distance[current];
    for (auto it = pending.rbegin(); it != pending.rend(); ++it)
        distance[*it] = ++value;
    return distance[cell];
}
}

// Ownership of a fine cell on several links (header changes of E2b: the contract
// gave it to the larger drainage area, which can make the bed rise where a gentle
// link enters a cell that a steeper, larger link holds higher up, with no way back
// down to the gentle link's own centre): a centre fine cell belongs to its own global
// cell's outgoing link; any other fine cell to the link with the lowest bed there,
// equal beds to the link whose upstream cell is fewer links from its outlet on the
// global tree, then to the smaller absolute index of the upstream global cell. Each
// cell's direction is its owner's step to the owner's next cell. A chain that
// crosses a cell taken by another link therefore runs into it and on along the
// owner's path (its own later cells stay its own where it owns them):
// - the bed never rises: the next cell's owner holds it at or below this link's bed
//   there, which is at or below this cell's (each link's bed descends: the inputs are
//   checked so that no receiver stands above its upstream cell);
// - no loop, without any assumption on the drainage areas: along the directions the
//   bed descends, or stays equal while the owner stays the same (its path position
//   grows) or changes to a link that ranks strictly higher in (outlet distance,
//   index), a fixed order of the links; the centre rule agrees with it, because the
//   only links that reach a centre (meandercap) are the cell's outgoing link, at the
//   cell's filled elevation, and the links into the cell, which end there at the same
//   bed and start one link further from the outlet;
// - every chain ends at a link's last cell, where the direction is 0: its first fine
//   cell in any sea block, or a land sink's centre. A link cut at a sea side
//   neighbour's corner drops to the sea level there, which is a descent (the sea
//   level is checked to be at or below the upstream cell).
// Links adjacent on the global tree (the same receiver, or one's receiver the other's
// upstream cell) may share fine cells, and a chain then runs on along the other link
// to the same outlet; links that are not adjacent share none (meandercap), except two
// diagonals crossing in an X, whose crossing cells, if their staircases share any,
// go to the lower bed, so that the higher link's chain continues into the lower
// link's basin (accepted: the X exists only where the global routing crosses its own
// flow on a flat).
FineChannels rasterisechannels(const RegionalWindow& window, const GlobalNetwork& network,
                               const ChannelParameters& parameters, float sealevel)
{
    const int k = window.scale;
    const int width = network.width;
    const int height = network.height;
    if (k < 1)
        throw std::invalid_argument("rasterisechannels: the scale must be at least 1");
    if (width < 1 || height < 1 || width != window.worldwidth || height != window.worldheight)
        throw std::invalid_argument("rasterisechannels: the network's grid must be the window's world grid");
    if (window.domaincolumns < 1 || window.domaincolumns > width || window.domainrows < 1 || window.domainrow0 < 0 ||
        window.domainrow0 + window.domainrows > height)
        throw std::invalid_argument("rasterisechannels: the domain must lie on the world grid");
    if (!(parameters.meanderamplitude >= 0.0) || !std::isfinite(parameters.meanderamplitude))
        throw std::invalid_argument("rasterisechannels: meanderamplitude must be finite and non-negative");
    if (!std::isfinite(sealevel))
        throw std::invalid_argument("rasterisechannels: the sea level must be finite");
    const std::size_t cells = static_cast<std::size_t>(width) * static_cast<std::size_t>(height);
    requiresize(network.direction.size(), cells, "direction");
    requiresize(network.sea.size(), cells, "sea");
    requiresize(network.drainagearea_km2.size(), cells, "drainagearea_km2");
    requiresize(network.discharge_m3s.size(), cells, "discharge_m3s");
    requiresize(network.filledelevation_m.size(), cells, "filledelevation_m");
    requiresize(network.januaryflow.size(), cells, "januaryflow");
    requiresize(network.julyflow.size(), cells, "julyflow");

    const int domainwidth = window.domaincolumns * k;
    const int domainheight = window.domainrows * k;
    const std::size_t finecells = static_cast<std::size_t>(domainwidth) * static_cast<std::size_t>(domainheight);
    FineChannels out;
    out.channel.assign(finecells, 0);
    out.direction.assign(finecells, 0);
    out.bed_m.assign(finecells, 0.0f);
    out.discharge_m3s.assign(finecells, 0.0f);
    out.januaryflow.assign(finecells, 0);
    out.julyflow.assign(finecells, 0);
    out.globalcell.assign(finecells, -1);
    std::vector<uint8_t> owncentre(finecells, 0);
    std::vector<int32_t> ownerdistance(finecells, 0);
    std::vector<int32_t> distance(cells, -1);

    // Domain fine cell of an unwrapped absolute fine cell: the absolute column
    // wrapped to 0 .. kW - 1, then counted eastwards from the domain's first fine
    // column (domaincolumn0 k) modulo kW; the row counted from domainrow0 k.
    const std::int64_t planetfinewidth = static_cast<std::int64_t>(k) * width;
    const std::int64_t domainx0 = static_cast<std::int64_t>(k) * window.domaincolumn0;
    const std::int64_t domainy0 = static_cast<std::int64_t>(k) * window.domainrow0;
    auto domainindex = [&](const FineCell& cell, std::size_t& index) {
        const std::int64_t y = cell.y - domainy0;
        if (y < 0 || y >= domainheight)
            return false;
        std::int64_t x = (cell.x - domainx0) % planetfinewidth;
        if (x < 0)
            x += planetfinewidth;
        if (x >= domainwidth)
            return false;
        index = static_cast<std::size_t>(y) * static_cast<std::size_t>(domainwidth) + static_cast<std::size_t>(x);
        return true;
    };

    // Upstream global cells within linkreach of the domain, each column once.
    std::vector<int> columns;
    if (window.domaincolumns + 2 * linkreach >= width)
    {
        for (int i = 0; i < width; i++)
            columns.push_back(i);
    }
    else
    {
        for (int m = -linkreach; m < window.domaincolumns + linkreach; m++)
            columns.push_back(((window.domaincolumn0 + m) % width + width) % width);
    }
    const int rowfirst = std::max(0, window.domainrow0 - linkreach);
    const int rowlast = std::min(height - 1, window.domainrow0 + window.domainrows - 1 + linkreach);

    for (int j = rowfirst; j <= rowlast; j++)
        for (int i : columns)
        {
            const std::size_t upstream = static_cast<std::size_t>(j) * width + i;
            const int direction = network.direction[upstream];
            if (network.sea[upstream] != 0 || direction < 1 || direction > directioncount ||
                !(network.drainagearea_km2[upstream] >= parameters.pinarea_km2))
                continue;
            const int ox = directionoffsetx(direction);
            const int oy = directionoffsety(direction);
            if (j + oy < 0 || j + oy >= height)
                continue;
            const std::size_t receiver = static_cast<std::size_t>(j + oy) * width + (i + ox + width) % width;

            const std::uint64_t hash = deterministiccombine(parameters.seed, static_cast<std::int64_t>(upstream),
                                                            static_cast<std::int64_t>(receiver));
            std::vector<FineCell> path = linkpath(k, i, j, ox, oy, parameters.meanderamplitude, hash);
            // A link ends at its first fine cell inside any global sea cell's k x k
            // block: the sea receiver's, or a sea side neighbour's whose corner a
            // diagonal or a meander crosses (a diagonal segment runs through the
            // corner its two side neighbours share). The rest of the path is not drawn.
            const auto insea = std::find_if(path.begin(), path.end(), [&](const FineCell& cell) {
                const std::int64_t column = cell.x >= 0 ? cell.x / k : (cell.x - (k - 1)) / k;  // floor
                const std::int64_t row = cell.y >= 0 ? cell.y / k : (cell.y - (k - 1)) / k;
                if (row < 0 || row >= height)
                    return false;
                const std::int64_t wrapped = ((column % width) + width) % width;
                return network.sea[static_cast<std::size_t>(row) * width + static_cast<std::size_t>(wrapped)] != 0;
            });
            const bool intosea = insea != path.end();
            if (intosea)
                path.erase(insea + 1, path.end());

            // Bed: linear in the path length (fine cell edges, sqrt(2) per diagonal
            // step) from the upstream filled elevation to the receiver's, or to the
            // sea level, exactly, for a link that ends in the sea. The global routing's
            // receivers lie on the filled surface and land stands at or above the sea
            // level, so an end above the upstream cell, or a filled elevation that is
            // not a number, is bad input.
            const float top = network.filledelevation_m[upstream];
            const float end = intosea ? sealevel : network.filledelevation_m[receiver];
            if (!std::isfinite(top) || !std::isfinite(end))
                throw std::invalid_argument("rasterisechannels: a drawn link's filled elevation is not finite (cell " +
                                            std::to_string(upstream) + ")");
            if (end > top)
                throw std::invalid_argument("rasterisechannels: a receiver (or the sea level, for a link that ends "
                                            "in the sea) stands above its upstream cell's filled elevation (cell " +
                                            std::to_string(upstream) + ")");
            std::vector<double> along(path.size(), 0.0);
            for (std::size_t m = 1; m < path.size(); m++)
            {
                const bool diagonal = path[m].x != path[m - 1].x && path[m].y != path[m - 1].y;
                along[m] = along[m - 1] + (diagonal ? std::sqrt(2.0) : 1.0);
            }
            const double total = along.back();
            const int32_t rank = outletdistance(network, upstream, distance);
            const int32_t upstreamindex = static_cast<int32_t>(upstream);

            for (std::size_t m = 0; m < path.size(); m++)
            {
                std::size_t index = 0;
                if (!domainindex(path[m], index))
                    continue;
                const bool lastcell = m + 1 == path.size();
                const float bed =
                    m == 0 ? top
                           : (lastcell ? end : static_cast<float>(top + (static_cast<double>(end) - top) * along[m] / total));
                const uint8_t own = m == 0 ? 1 : 0;
                if (out.channel[index] != 0)
                {
                    bool better = false;
                    if (own != owncentre[index])
                        better = own != 0;
                    else if (bed != out.bed_m[index])
                        better = bed < out.bed_m[index];
                    else if (rank != ownerdistance[index])
                        better = rank < ownerdistance[index];
                    else
                        better = upstreamindex < out.globalcell[index];
                    if (!better)
                        continue;
                }
                out.channel[index] = 1;
                out.direction[index] = lastcell ? 0
                                                : static_cast<uint8_t>(directioncode(path[m + 1].x - path[m].x,
                                                                                     path[m + 1].y - path[m].y));
                out.bed_m[index] = bed;
                out.discharge_m3s[index] = network.discharge_m3s[upstream];
                out.januaryflow[index] = network.januaryflow[upstream];
                out.julyflow[index] = network.julyflow[upstream];
                out.globalcell[index] = upstreamindex;
                owncentre[index] = own;
                ownerdistance[index] = rank;
            }
        }
    return out;
}
}
