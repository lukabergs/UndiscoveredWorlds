#include "simulations/landscape/channel_network.hpp"

// Pinned channels of the regional pass (erosion plan step E2b,
// docs/architecture/erosion-hydrology-plan.md, Phase E and section 6, "Seams of the
// solve"; Phase R, step R1: the fine coast and the start surface's running minimum;
// step R3: the legacy geometry of junctions, edge points and Catmull-Rom legs). Every
// quantity below is computed in a global cell's own block coordinates from that cell,
// its edges, its neighbours on the global network (and their outlets on the global tree,
// S3.10's transit blocks) and the seed, and only mapped to the
// window's domain when a fine cell is written, so a link is the same fine path in every
// window, and a fine cell's values depend only on the pieces whose paths cross it.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "core/deterministic_random.hpp"
#include "simulations/landscape/refinement_levels.hpp"
#include "simulations/landscape/regional_pass.hpp"
#include "simulations/landscape/regional_tributaries.hpp"
#include "simulations/landscape/sphere_metric.hpp"

namespace landscape
{
namespace
{
constexpr double pi = 3.14159265358979323846;
// Global cells around the domain whose outgoing links are drawn: a link lies inside
// its two cells' blocks (a diagonal through a transit block, S3.10, inside the 2 x 2
// blocks around its corner, mutual 8-neighbours), so only a link whose upstream cell is
// within one cell of the domain reaches it. The confluences (S3.11) read the first lost cells of
// the pieces crossing a painted piece's path, and those read the owners on their own paths: the
// links and springs of cells within three cells of the domain are drawn (outside it they compete
// for ownership and end prefixes, never painted).
constexpr int linkreach = 3;
// Hash salts of the geometry's pieces.
constexpr std::int64_t junctionsalt = 0x4a554e43;
constexpr std::int64_t edgesalt = 0x45444745;
constexpr std::int64_t legsalt = 0x4c454753;
constexpr std::int64_t springsalt = 0x53505247;
constexpr std::int64_t transitsalt = 0x5452414e;
// The edge hashes' first value: 0 and 1 a cardinal link's vertical and horizontal edge, 2 and 3 a
// transit link's first edge (upstream | transit block) and second (transit block | receiver).
constexpr std::int64_t transitfirstedge = 2;
constexpr std::int64_t transitsecondedge = 3;
// A spring starts within this angle of the direction opposite the junction's tangent
// (upstream of the junction), between springnear and springfar x k from it.
constexpr double springspread = 2.0 * pi / 3.0;
constexpr double springnear = 0.25;
constexpr double springfar = 0.45;
constexpr std::size_t nocell = static_cast<std::size_t>(-1);

// A fine cell: block-local (0 .. k - 1) while a leg is drawn, then an unwrapped
// absolute fine cell (column outside 0 .. kW - 1 for a link across the date line;
// wrapped when written).
struct FineCell
{
    std::int64_t x = 0;
    std::int64_t y = 0;
};

struct Point
{
    double x = 0.0;
    double y = 0.0;
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

// A value in [-1, 1) and one in [0, 1) from the top 53 bits of a hash.
double signedunit(std::uint64_t hash)
{
    return static_cast<double>(hash >> 11) * (2.0 / 9007199254740992.0) - 1.0;
}

double unitinterval(std::uint64_t hash) { return static_cast<double>(hash >> 11) * (1.0 / 9007199254740992.0); }

// The block-local fine index of a coordinate: its floor, clamped to 0 .. k - 1.
std::int64_t blockindex(double value, int k)
{
    const double floored = std::floor(value);
    if (!(floored >= 0.0))
        return 0;
    if (floored > static_cast<double>(k - 1))
        return k - 1;
    return static_cast<std::int64_t>(floored);
}

Point middle(const FineCell& cell) { return { static_cast<double>(cell.x) + 0.5, static_cast<double>(cell.y) + 0.5 }; }

double length(const Point& p) { return std::hypot(p.x, p.y); }

Point unitvector(const Point& p)
{
    const double size = length(p);
    return size > 1e-9 ? Point{ p.x / size, p.y / size } : Point{};
}

bool iszero(const Point& p) { return p.x == 0.0 && p.y == 0.0; }

// The junction of a global cell (row-major index `cell`) in its block's coordinates.
FineCell junctionoffset(int k, std::size_t cell, const ChannelParameters& parameters)
{
    const std::uint64_t hash =
        deterministiccombine(parameters.seed, junctionsalt, static_cast<std::int64_t>(cell));
    const double reach = parameters.junctionjitter * k;
    const double centre = 0.5 * k;
    return { blockindex(centre + reach * signedunit(hash), k),
             blockindex(centre + reach * signedunit(deterministicmix64(hash)), k) };
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

// The outlet of a global cell on the drainage tree (where globalreceiver stops: a sea
// cell, a land sink, a direction beyond a pole), memoised in `outlet` (nocell unknown).
// Receivers that form a loop throw.
std::size_t outletcell(const GlobalNetwork& network, std::size_t cell, std::vector<std::size_t>& outlet)
{
    std::vector<std::size_t> pending;
    std::size_t current = cell;
    while (outlet[current] == nocell)
    {
        const std::size_t next = globalreceiver(network, current);
        if (next == current)
        {
            outlet[current] = current;
            break;
        }
        pending.push_back(current);
        if (pending.size() > outlet.size())
            throw std::invalid_argument("rasterisechannels: the global receivers form a loop");
        current = next;
    }
    for (const std::size_t waiting : pending)
        outlet[waiting] = outlet[current];
    return outlet[cell];
}

// The global cell (i + ox, j + oy), columns wrapped (the row must lie on the grid).
std::size_t offsetcell(const GlobalNetwork& network, std::size_t cell, int ox, int oy)
{
    const int width = network.width;
    const int i = static_cast<int>(cell % static_cast<std::size_t>(width));
    const int j = static_cast<int>(cell / static_cast<std::size_t>(width));
    return static_cast<std::size_t>(j + oy) * static_cast<std::size_t>(width) +
           static_cast<std::size_t>(((i + ox) % width + width) % width);
}

// The transit block of a diagonal link from `upstream` to `receiver` (offsets ox, oy, both
// non-zero; multiscale plan S3.10): 1 the orthogonal neighbour across the link's column step
// (upstream + (ox, 0)), 2 the one across its row step (upstream + (0, oy)), 0 the shared corner.
// With diagonaltransit on, a neighbour qualifies when it is land, its filled elevation lies
// between the receiver's and the upstream cell's (inclusive), it is neither the receiver's
// receiver nor a donor of the upstream cell, and its outlet on the global tree is the
// receiver's, so every piece in a block drains to that block's outlet; both qualify: a hash
// of the link; neither (a neck between basins, an X on a flat, a block above or below the
// link) or diagonaltransit 0: the corner.
int transitside(const GlobalNetwork& network, std::size_t upstream, std::size_t receiver, int ox, int oy,
                const ChannelParameters& parameters, std::vector<std::size_t>& outlet)
{
    if (parameters.diagonaltransit == 0)
        return 0;
    const std::size_t target = outletcell(network, receiver, outlet);
    const float top = network.filledelevation_m[upstream];
    const float bottom = network.filledelevation_m[receiver];
    // The miss of a candidate: how far its filled elevation lies outside [bottom, top] (0
    // between the link's ends: no pinned trench through a higher block, no levee or dam across
    // a lower one); infinity for one that is the upstream cell or the receiver, sea, the
    // receiver's receiver or the upstream cell's donor (a hairpin), or of another outlet.
    const double infinity = std::numeric_limits<double>::infinity();
    auto miss = [&](std::size_t neighbour) {
        if (neighbour == upstream || neighbour == receiver || network.sea[neighbour] != 0 ||
            globalreceiver(network, receiver) == neighbour || globalreceiver(network, neighbour) == upstream ||
            outletcell(network, neighbour, outlet) != target)
            return infinity;
        const double filled = network.filledelevation_m[neighbour];
        const double value = std::max({ 0.0, static_cast<double>(bottom) - filled, filled - static_cast<double>(top) });
        return std::isnan(value) ? infinity : value;
    };
    const double first = miss(offsetcell(network, upstream, ox, 0));
    const double second = miss(offsetcell(network, upstream, 0, oy));
    // Between the ends (miss 0) qualifies; when neither does, the smaller miss within
    // transittolerance_m; equal misses: a hash of the link; none: the corner.
    const double limit = std::min(first, second) > 0.0 ? parameters.transittolerance_m : 0.0;
    const bool firstok = first <= limit, secondok = second <= limit;
    if (firstok && secondok && first == second)
        return (deterministiccombine(parameters.seed, transitsalt, static_cast<std::int64_t>(upstream),
                                     static_cast<std::int64_t>(receiver)) & 1u) != 0 ? 2 : 1;
    if (firstok && (!secondok || first < second))
        return 1;
    return secondok ? 2 : 0;
}

// The exit point of a link in its upstream cell's block and the entry point in its
// receiver's block (block coordinates): a cardinal link across their shared edge at an
// offset along the edge drawn from the edge's two cells in a canonical order (west then
// east, north then south), so the two cells agree whichever drains into the other, between
// edgemargin x k and (1 - edgemargin) x k; a diagonal link through the shared corner, or
// (transitside, S3.10) through the transit block: the exit on the edge between the upstream
// cell and the transit block, the entry on the edge between the transit block and the
// receiver, each at an offset drawn from the link (edgesalt, upstream, receiver) between
// the same margins, with the transit leg's two ends across those edges.
struct LinkEnds
{
    FineCell exit;
    FineCell entry;
    int transit = 0;          // transitside: 0 a cardinal link or a diagonal through the corner
    FineCell transitentry;    // in the transit block, across the first edge from `exit`
    FineCell transitexit;     // in the transit block, across the second edge from `entry`
};

std::int64_t edgeoffset(int k, std::uint64_t hash, const ChannelParameters& parameters)
{
    return blockindex(parameters.edgemargin * k + (1.0 - 2.0 * parameters.edgemargin) * k * unitinterval(hash), k);
}

LinkEnds linkends(int k, const GlobalNetwork& network, std::size_t upstream, std::size_t receiver, int ox, int oy,
                  const ChannelParameters& parameters, std::vector<std::size_t>& outlet)
{
    const std::int64_t last = k - 1;
    LinkEnds ends;
    if (ox != 0 && oy != 0)
    {
        ends.transit = transitside(network, upstream, receiver, ox, oy, parameters, outlet);
        if (ends.transit == 0)
        {
            ends.exit = { ox > 0 ? last : 0, oy > 0 ? last : 0 };
            ends.entry = { ox > 0 ? 0 : last, oy > 0 ? 0 : last };
            return ends;
        }
        const std::int64_t first = edgeoffset(
            k,
            deterministiccombine(parameters.seed, edgesalt, transitfirstedge, static_cast<std::int64_t>(upstream),
                                 static_cast<std::int64_t>(receiver)),
            parameters);
        const std::int64_t second = edgeoffset(
            k,
            deterministiccombine(parameters.seed, edgesalt, transitsecondedge, static_cast<std::int64_t>(upstream),
                                 static_cast<std::int64_t>(receiver)),
            parameters);
        if (ends.transit == 1)
        {
            // The transit block lies across the upstream cell's east or west edge, the
            // receiver across the transit block's south or north edge.
            ends.exit = { ox > 0 ? last : 0, first };
            ends.transitentry = { ox > 0 ? 0 : last, first };
            ends.transitexit = { second, oy > 0 ? last : 0 };
            ends.entry = { second, oy > 0 ? 0 : last };
        }
        else
        {
            // Across the upstream cell's south or north edge, then the transit block's east or west.
            ends.exit = { first, oy > 0 ? last : 0 };
            ends.transitentry = { first, oy > 0 ? 0 : last };
            ends.transitexit = { ox > 0 ? last : 0, second };
            ends.entry = { ox > 0 ? 0 : last, second };
        }
        return ends;
    }
    const bool vertical = ox != 0;  // the shared edge runs north to south
    const bool forward = vertical ? ox > 0 : oy > 0;
    const std::size_t first = forward ? upstream : receiver;
    const std::size_t second = forward ? receiver : upstream;
    const std::uint64_t hash = deterministiccombine(parameters.seed, edgesalt, vertical ? 0 : 1,
                                                    static_cast<std::int64_t>(first), static_cast<std::int64_t>(second));
    const std::int64_t offset = edgeoffset(k, hash, parameters);
    if (vertical)
    {
        ends.exit = { forward ? last : 0, offset };
        ends.entry = { forward ? 0 : last, offset };
    }
    else
    {
        ends.exit = { offset, forward ? last : 0 };
        ends.entry = { offset, forward ? 0 : last };
    }
    return ends;
}

// Whether a global cell's link is drawn (land, a receiver on the grid, the drainage
// area at least pinarea_km2), and its receiver.
bool drawnlink(const GlobalNetwork& network, std::size_t cell, double pinarea_km2, std::size_t& receiver)
{
    receiver = globalreceiver(network, cell);
    return receiver != cell && network.drainagearea_km2[cell] >= pinarea_km2;
}

// A global cell's junction (block coordinates), its main upstream cell (the drawn link
// into it with the largest drainage area, ties to the smaller index; nocell for a
// source) and the unit tangent of its main stem at the junction: from the main entry
// point towards the exit point, or along the one of them that exists (zero if neither
// does or the points coincide).
struct CellGeometry
{
    FineCell junction;
    std::size_t mainupstream = nocell;
    Point tangent;
};

CellGeometry cellgeometry(int k, const GlobalNetwork& network, std::size_t cell, const ChannelParameters& parameters,
                          std::vector<std::size_t>& outlet)
{
    CellGeometry geometry;
    geometry.junction = junctionoffset(k, cell, parameters);
    const int width = network.width;
    const int i = static_cast<int>(cell % static_cast<std::size_t>(width));
    const int j = static_cast<int>(cell / static_cast<std::size_t>(width));
    float best = 0.0f;
    for (int direction = 1; direction <= directioncount; direction++)
    {
        const int nj = j + directionoffsety(direction);
        if (nj < 0 || nj >= network.height)
            continue;
        const int ni = ((i + directionoffsetx(direction)) % width + width) % width;
        const std::size_t neighbour = static_cast<std::size_t>(nj) * static_cast<std::size_t>(width) + ni;
        std::size_t receiver = nocell;
        if (neighbour == cell || !drawnlink(network, neighbour, parameters.pinarea_km2, receiver) || receiver != cell)
            continue;
        const float area = network.drainagearea_km2[neighbour];
        if (geometry.mainupstream == nocell || area > best || (area == best && neighbour < geometry.mainupstream))
        {
            best = area;
            geometry.mainupstream = neighbour;
        }
    }

    const Point junction = middle(geometry.junction);
    bool hasexit = false, hasentry = false;
    Point exit, entry;
    const std::size_t receiver = globalreceiver(network, cell);
    if (receiver != cell)
    {
        const int direction = network.direction[cell];
        exit = middle(linkends(k, network, cell, receiver, directionoffsetx(direction), directionoffsety(direction),
                               parameters, outlet)
                          .exit);
        hasexit = true;
    }
    if (geometry.mainupstream != nocell)
    {
        const int direction = network.direction[geometry.mainupstream];
        entry = middle(linkends(k, network, geometry.mainupstream, cell, directionoffsetx(direction),
                                directionoffsety(direction), parameters, outlet)
                           .entry);
        hasentry = true;
    }
    Point direction;
    if (hasexit && hasentry)
        direction = { exit.x - entry.x, exit.y - entry.y };
    if (length(direction) < 0.5)
    {
        if (hasexit)
            direction = { exit.x - junction.x, exit.y - junction.y };
        else if (hasentry)
            direction = { junction.x - entry.x, junction.y - entry.y };
    }
    geometry.tangent = unitvector(direction);
    return geometry;
}

// The largest sideways displacement of a leg's midpoint, in fine cells.
double meanderdisplacement(const ChannelParameters& parameters, int k, float discharge_m3s)
{
    double scale = 1.0;
    if (parameters.meanderdischargeexponent > 0.0)
    {
        const double ratio = std::max(0.0, static_cast<double>(discharge_m3s)) / parameters.meanderreference_m3s;
        scale = std::min(1.0, std::pow(ratio, parameters.meanderdischargeexponent));
    }
    return parameters.meanderamplitude * k * scale;
}

// The fine cells of a leg inside one k x k block, from `from` to `to` (block
// coordinates), as an 8-connected chain of distinct cells, each step to an
// 8-neighbour. The curve is a Catmull-Rom spline (cubic Hermite segments) through the
// two cells' middle points P0, P2 and a midpoint P1: the chord's middle displaced along
// its left normal by +-(0.5 to 1) x `displacement` (sign and size from the hash),
// clamped inside the block. The tangent at P1 is (P2 - P0) / 2; at an end, a given unit
// direction times the adjacent segment's chord, or the chord itself (a zero direction).
// Each segment is sampled at steps under a quarter of a fine cell (the Hermite curve's
// speed is at most 1.5 |B - A| + |mA| + |mB|), each sample's fine cell clamped into the
// block; a cell already on the chain cuts the chain back to it (the loop is erased), and
// a cell that neighbours the cell before the chain's last replaces that last cell (an
// L-shaped corner becomes a diagonal step). The chain starts at `from` and ends at `to`.
// Multiscale plan S3.11 (header, "Roughness"): each sample moves along the curve's left
// normal by the octaves' sum d(t), t = s / 2 over the first segment and 1/2 + s / 2 over the
// second (s the segment's parameter); the sample bound grows by the octaves' largest speed,
// sum a_i (pi / 2) 2^i per unit s, and a sample cell that is not an 8-neighbour of the chain's
// last (the normal turning fast) is reached by diagonal-first steps. With `roughness` 0 there
// are no octaves and the chain is the smooth arc's.
std::vector<FineCell> legpath(int k, const FineCell& from, const FineCell& to, double displacement,
                              std::uint64_t hash, const Point& startdirection, const Point& enddirection,
                              double roughness, double persistence)
{
    std::vector<FineCell> chain{ from };
    if (samecell(from, to))
        return chain;
    const Point p0 = middle(from);
    const Point p2 = middle(to);
    const double cx = p2.x - p0.x;
    const double cy = p2.y - p0.y;
    const double chord = std::hypot(cx, cy);
    const double side = (hash & 1u) != 0 ? 1.0 : -1.0;
    const double offset = side * displacement * (0.5 + 0.5 * unitinterval(deterministicmix64(hash)));
    const double low = 0.5;
    const double high = static_cast<double>(k) - 0.5;
    const Point p1{ std::clamp(0.5 * (p0.x + p2.x) - cy / chord * offset, low, high),
                    std::clamp(0.5 * (p0.y + p2.y) + cx / chord * offset, low, high) };
    const double first = std::hypot(p1.x - p0.x, p1.y - p0.y);
    const double second = std::hypot(p2.x - p1.x, p2.y - p1.y);
    const Point m1{ 0.5 * cx, 0.5 * cy };
    const Point m0 = iszero(startdirection) ? Point{ p1.x - p0.x, p1.y - p0.y }
                                            : Point{ startdirection.x * first, startdirection.y * first };
    const Point m2 = iszero(enddirection) ? Point{ p2.x - p1.x, p2.y - p1.y }
                                          : Point{ enddirection.x * second, enddirection.y * second };

    // S3.11: the octaves (2^i intervals, amplitude a_i, node values; nodes 0 and 2^i are 0).
    struct Octave
    {
        int intervals = 0;
        double amplitude = 0.0;
        std::vector<double> nodes;
    };
    std::vector<Octave> octaves;
    double roughspeed = 0.0;
    if (roughness > 0.0)
    {
        double amplitude = roughness * k;
        for (int i = 2; amplitude > 0.0 && chord / std::ldexp(1.0, i) >= 2.0; i++)
        {
            Octave octave;
            octave.intervals = 1 << i;
            octave.amplitude = amplitude;
            octave.nodes.assign(static_cast<std::size_t>(octave.intervals) + 1, 0.0);
            for (int node = 1; node < octave.intervals; node++)
                octave.nodes[static_cast<std::size_t>(node)] = signedunit(deterministiccombine(hash, i, node));
            roughspeed += amplitude * 0.5 * pi * octave.intervals;
            octaves.push_back(std::move(octave));
            amplitude *= persistence;
        }
    }
    auto rough = [&](double t) {
        double sum = 0.0;
        for (const Octave& octave : octaves)
        {
            const double x = t * octave.intervals;
            const int node = std::clamp(static_cast<int>(std::floor(x)), 0, octave.intervals - 1);
            const double w = 0.5 - 0.5 * std::cos(pi * (x - node));
            sum += octave.amplitude * ((1.0 - w) * octave.nodes[static_cast<std::size_t>(node)] +
                                       w * octave.nodes[static_cast<std::size_t>(node) + 1]);
        }
        return sum;
    };

    auto place = [&](const FineCell& cell) {
        if (samecell(cell, chain.back()))
            return;
        const auto found =
            std::find_if(chain.begin(), chain.end(), [&](const FineCell& other) { return samecell(other, cell); });
        if (found != chain.end())
        {
            chain.erase(found + 1, chain.end());
            return;
        }
        while (chain.size() >= 2 && neighbours8(chain[chain.size() - 2], cell))
            chain.pop_back();
        chain.push_back(cell);
    };
    // A cell beyond the chain's last cell's 8-neighbours (only with roughness: without it the
    // samples step under a quarter of a fine cell) is reached by diagonal-first steps; each
    // placed step becomes the chain's last cell, one step nearer.
    auto add = [&](const FineCell& cell) {
        while (!samecell(cell, chain.back()) && !neighbours8(cell, chain.back()))
        {
            const FineCell& last = chain.back();
            place({ last.x + (cell.x > last.x ? 1 : (cell.x < last.x ? -1 : 0)),
                    last.y + (cell.y > last.y ? 1 : (cell.y < last.y ? -1 : 0)) });
        }
        place(cell);
    };
    auto segment = [&](const Point& a, const Point& b, const Point& ma, const Point& mb, bool lastsegment) {
        const double bound = 1.5 * std::hypot(b.x - a.x, b.y - a.y) + length(ma) + length(mb) + roughspeed;
        const int samples = static_cast<int>(std::ceil(4.0 * bound)) + 1;
        for (int s = 1; s <= samples; s++)
        {
            if (lastsegment && s == samples)
            {
                add(to);
                continue;
            }
            const double t = static_cast<double>(s) / samples;
            const double t2 = t * t;
            const double t3 = t2 * t;
            const double h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
            const double h10 = t3 - 2.0 * t2 + t;
            const double h01 = -2.0 * t3 + 3.0 * t2;
            const double h11 = t3 - t2;
            double x = h00 * a.x + h10 * ma.x + h01 * b.x + h11 * mb.x;
            double y = h00 * a.y + h10 * ma.y + h01 * b.y + h11 * mb.y;
            if (!octaves.empty())
            {
                // The curve's derivative (Hermite basis derivatives), its left normal.
                const double d00 = 6.0 * t2 - 6.0 * t;
                const double d10 = 3.0 * t2 - 4.0 * t + 1.0;
                const double d01 = -d00;
                const double d11 = 3.0 * t2 - 2.0 * t;
                const double dx = d00 * a.x + d10 * ma.x + d01 * b.x + d11 * mb.x;
                const double dy = d00 * a.y + d10 * ma.y + d01 * b.y + d11 * mb.y;
                const double speed = std::hypot(dx, dy);
                const Point normal = speed > 1e-9 ? Point{ -dy / speed, dx / speed } : Point{ -cy / chord, cx / chord };
                const double d = rough((lastsegment ? 0.5 : 0.0) + 0.5 * t);
                x += d * normal.x;
                y += d * normal.y;
            }
            add({ blockindex(x, k), blockindex(y, k) });
        }
    };
    segment(p0, p1, m0, m1, false);
    segment(p1, p2, m1, m2, true);
    return chain;
}

void requiresize(std::size_t size, std::size_t expected, const char* name)
{
    if (size != expected)
        throw std::invalid_argument(std::string("rasterisechannels: GlobalNetwork.") + name +
                                    " does not have width x height values");
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

FinePosition channeljunction(const RegionalWindow& window, const ChannelParameters& parameters, int i, int j)
{
    const int k = window.scale;
    const std::size_t cell = static_cast<std::size_t>(j) * static_cast<std::size_t>(window.worldwidth) +
                             static_cast<std::size_t>(i);
    const FineCell offset = junctionoffset(k, cell, parameters);
    return { static_cast<std::int64_t>(k) * i + offset.x, static_cast<std::int64_t>(k) * j + offset.y };
}

// Ownership of a fine cell on several pieces (links and springs): the piece with the
// lowest bed there, equal beds to the piece whose upstream cell is fewer links from its
// outlet on the global tree, then to the smaller absolute index of the upstream global
// cell, then to the link before its springs and to the lower spring number (a strict
// order of the pieces at a cell, so the result does not depend on the drawing order).
// Each cell's direction is its owner's step to the owner's next cell. A chain that
// crosses a cell taken by another piece therefore runs into it and on along the owner's
// path:
// - the bed never rises: the next cell's owner holds it at or below this piece's bed
//   there, which is at or below this cell's, because each piece's bed never rises along
//   its own path: the junction beds satisfy b(c) <= F(c) <= lo(u) <= b(u) for a link
//   from u to a land receiver c (F the filled elevation, checked so that no receiver
//   stands above its upstream cell), and the sea level <= min(F(u), next float above
//   it) <= b(u) for a link that ends in the sea (F(u) is checked to be at or above the
//   sea level), so end <= top; the linear profile runs from top to end; the R1 bed is
//   top at the first cell, then max(running, floor), where the running minimum starts at
//   top and never rises and end <= floor <= top, so it never rises and stays in
//   [end, top], then end at the last cell; a spring is the same with its own top, which
//   is at least its end;
// - no loop, without any assumption on the drainage areas: along the directions the
//   bed descends, or stays equal while the owner stays the same (its path position
//   grows) or changes to a piece that ranks strictly higher in (outlet distance, index,
//   piece), a fixed order of the pieces; equal beds are common under the R1 bed (a
//   running minimum is flat between the start surface's new lows, and the floor is
//   flat), and each such step is one of these two. A junction needs no rule of its own
//   (until S3.10 it went to its cell's outgoing link first, which gave the same owners):
//   without transit legs the only pieces that reach it are its cell's outgoing link (at
//   b(u), the first cell), the links into the cell and the cell's springs, which end
//   there (a chain holds each cell once, so the junction is their last cell) at the
//   same bed b(u), and rank lower (one link further from the outlet, or the same link's
//   springs); where a fine sea makes the junction sea, the outgoing link ends at its
//   first cell and every piece that reaches it ends there too, all at the sea level. A
//   transit leg crossing a junction takes it where its bed is lower (then the cell's
//   outgoing link starts on the transit's path, which reaches the same outlet); a rule
//   that gave the junction to the outgoing link first would let a lower transit cell
//   step up onto it and back (a loop);
// - every chain ends at a piece's last cell, where the direction is 0: its first fine
//   cell in a sea block (no fine sea) or its first fine sea cell (R1), or a land sink's
//   junction. A link cut at a fine coast before its receiver drops to the sea level
//   there, which is a descent.
// A block holds only its cell's outflow leg and springs, the inflow legs of the links
// into the cell and the transit legs through it (S3.10), and every one of those pieces
// drains to the block's outlet on the global tree (a transit block's outlet is its
// receiver's, transitside), so a chain that runs onto another piece reaches the same
// outlet: no capture across a divide. Without transits only links adjacent on the global
// tree (the same receiver, or one's receiver the other's upstream cell) share fine cells,
// and two diagonals crossing in an X (flats routed by flood order) pass through different
// corner cells.
FineChannels rasterisechannels(const RegionalWindow& window, const GlobalNetwork& network,
                               const ChannelParameters& parameters, float sealevel,
                               const std::vector<float>& start_m, const std::vector<uint8_t>& finesea,
                               const std::vector<float>& relief_m)
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
    if (!(parameters.meanderdischargeexponent >= 0.0) || !std::isfinite(parameters.meanderdischargeexponent) ||
        !(parameters.meanderreference_m3s > 0.0) || !std::isfinite(parameters.meanderreference_m3s))
        throw std::invalid_argument("rasterisechannels: the meander's discharge exponent must be finite and "
                                    "non-negative, its reference discharge positive and finite");
    if (!(parameters.meanderroughness >= 0.0) || !std::isfinite(parameters.meanderroughness) ||
        !(parameters.meanderroughnesspersistence >= 0.0) || !(parameters.meanderroughnesspersistence <= 1.0))
        throw std::invalid_argument("rasterisechannels: meanderroughness must be finite and non-negative, "
                                    "meanderroughnesspersistence in 0 .. 1");
    if (!(parameters.meanderroughnessreliefgain >= 0.0) || !std::isfinite(parameters.meanderroughnessreliefgain) ||
        !(parameters.meanderroughnessrelief_m > 0.0) || !std::isfinite(parameters.meanderroughnessrelief_m))
        throw std::invalid_argument("rasterisechannels: meanderroughnessreliefgain must be finite and non-negative, "
                                    "meanderroughnessrelief_m positive and finite");
    if (!(parameters.junctionjitter >= 0.0) || !std::isfinite(parameters.junctionjitter) ||
        !(parameters.edgemargin >= 0.0) || !(parameters.edgemargin <= 0.5))
        throw std::invalid_argument("rasterisechannels: junctionjitter must be finite and non-negative, edgemargin "
                                    "in 0 .. 0.5");
    if (parameters.springspercell < 0)
        throw std::invalid_argument("rasterisechannels: springspercell must not be negative");
    if (parameters.diagonaltransit != 0 && parameters.diagonaltransit != 1)
        throw std::invalid_argument("rasterisechannels: diagonaltransit must be 0 or 1");
    if (!(parameters.transittolerance_m >= 0.0))
        throw std::invalid_argument("rasterisechannels: transittolerance_m must not be negative");
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
    if (!relief_m.empty() && relief_m.size() != cells)
        throw std::invalid_argument("rasterisechannels: the relief must be empty or hold the world grid's cells");
    // S3.14: the roughness of a piece's legs, meanderroughness x (1 + meanderroughnessreliefgain x
    // min(1, R / meanderroughnessrelief_m)), R the local relief of its upstream global cell (0 without
    // relief_m, for a non-finite value or below 0); roughness 0 stays 0, gain 0 the plain roughness.
    auto roughnessof = [&](std::size_t upstream) {
        if (parameters.meanderroughness == 0.0 || relief_m.empty())
            return parameters.meanderroughness;
        const double relief = relief_m[upstream];
        const double share =
            std::isfinite(relief) && relief > 0.0 ? std::min(1.0, relief / parameters.meanderroughnessrelief_m) : 0.0;
        return parameters.meanderroughness * (1.0 + parameters.meanderroughnessreliefgain * share);
    };

    const int domainwidth = window.domaincolumns * k;
    const int domainheight = window.domainrows * k;
    const std::size_t finecells = static_cast<std::size_t>(domainwidth) * static_cast<std::size_t>(domainheight);
    if (!start_m.empty() && start_m.size() != finecells)
        throw std::invalid_argument("rasterisechannels: the start surface must be empty or hold the domain's fine cells");
    if (!finesea.empty() && finesea.size() != finecells)
        throw std::invalid_argument("rasterisechannels: the fine sea must be empty or hold the domain's fine cells");
    const bool usestart = !start_m.empty();
    const bool usesea = !finesea.empty();
    const float abovesea = std::nextafter(sealevel, std::numeric_limits<float>::infinity());
    FineChannels out;
    out.channel.assign(finecells, 0);
    out.direction.assign(finecells, 0);
    out.bed_m.assign(finecells, 0.0f);
    out.discharge_m3s.assign(finecells, 0.0f);
    out.januaryflow.assign(finecells, 0);
    out.julyflow.assign(finecells, 0);
    out.globalcell.assign(finecells, -1);
    std::vector<int32_t> distance(cells, -1);
    std::vector<std::size_t> outlet(cells, nocell);

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
    // Whether an unwrapped absolute fine cell lies in a global sea cell's block.
    auto inseablock = [&](const FineCell& cell) {
        const std::int64_t column = cell.x >= 0 ? cell.x / k : (cell.x - (k - 1)) / k;  // floor
        const std::int64_t row = cell.y >= 0 ? cell.y / k : (cell.y - (k - 1)) / k;
        if (row < 0 || row >= height)
            return false;
        const std::int64_t wrapped = ((column % width) + width) % width;
        return network.sea[static_cast<std::size_t>(row) * width + static_cast<std::size_t>(wrapped)] != 0;
    };
    // A piece ends at its first fine sea cell inside the domain (R1) or, without a fine
    // sea or outside the domain, at its first fine cell in a global sea block; the rest
    // of the path is erased. Returns whether it was cut.
    auto cutatsea = [&](std::vector<FineCell>& path) {
        const auto insea = std::find_if(path.begin(), path.end(), [&](const FineCell& cell) {
            std::size_t index = 0;
            if (usesea && domainindex(cell, index))
                return finesea[index] != 0;
            return inseablock(cell);
        });
        if (insea == path.end())
            return false;
        path.erase(insea + 1, path.end());
        return true;
    };
    // The junction bed of a global cell (header): the start surface at its junction,
    // held between lo and the filled elevation.
    auto junctionbed = [&](std::size_t cell) {
        const float filled = network.filledelevation_m[cell];
        if (!usestart)
            return filled;
        const FineCell offset = junctionoffset(k, cell, parameters);
        const FineCell junction{ static_cast<std::int64_t>(k) * static_cast<std::int64_t>(cell % width) + offset.x,
                                 static_cast<std::int64_t>(k) * static_cast<std::int64_t>(cell / width) + offset.y };
        std::size_t index = 0;
        if (!domainindex(junction, index) || std::isnan(start_m[index]))
            return filled;
        const std::size_t receiver = globalreceiver(network, cell);
        float base = -std::numeric_limits<float>::infinity();
        if (receiver != cell)
            base = network.sea[receiver] != 0 ? sealevel : network.filledelevation_m[receiver];
        if (!std::isfinite(base))
            base = -std::numeric_limits<float>::infinity();
        const float lo = std::max(base, std::min(filled, abovesea));
        return std::max(lo, std::min(filled, start_m[index]));
    };
    // The pieces (links and springs) as drawn, before the confluences: per path cell its domain
    // index (nocell outside the domain), its ownership slot (the domain index, or one of the cells
    // outside the domain that a drawn piece crosses, numbered from finecells on), its step to the
    // next cell (0 at the last) and its bed; `end`, `lost` and `suffix` are the confluences' (below).
    struct Piece
    {
        std::vector<std::size_t> index;
        std::vector<std::size_t> slot;
        std::vector<uint8_t> direction;
        std::vector<float> bed;
        std::size_t upstream = 0;
        int32_t rank = 0;
        int32_t number = 0;  // 0 the link, 1, 2, ... its springs
        bool kept = false;   // the upstream cell lies in the kept box
        std::size_t end = 0;         // the prefix's last cell (its first lost cell, else the path's last)
        bool lost = false;           // the prefix ends on a cell another piece owns
        std::size_t suffix = nocell; // the suffix's first cell (nocell: none)
    };
    std::vector<Piece> pieces;
    std::unordered_map<std::int64_t, std::size_t> outsideslots;  // wrapped absolute fine cell -> slot
    auto startpiece = [&](const std::vector<FineCell>& path, std::size_t upstream, int32_t rank, int32_t number,
                          bool kept) -> Piece& {
        Piece& piece = pieces.emplace_back();
        piece.index.assign(path.size(), nocell);
        piece.slot.assign(path.size(), nocell);
        piece.direction.assign(path.size(), 0);
        piece.bed.assign(path.size(), 0.0f);
        for (std::size_t m = 0; m < path.size(); m++)
        {
            std::size_t index = 0;
            if (domainindex(path[m], index))
                piece.index[m] = piece.slot[m] = index;
            else
            {
                const std::int64_t column = ((path[m].x % planetfinewidth) + planetfinewidth) % planetfinewidth;
                const auto inserted =
                    outsideslots.try_emplace(path[m].y * planetfinewidth + column, finecells + outsideslots.size());
                piece.slot[m] = inserted.first->second;
            }
            if (m + 1 < path.size())
                piece.direction[m] =
                    static_cast<uint8_t>(directioncode(path[m + 1].x - path[m].x, path[m + 1].y - path[m].y));
        }
        piece.upstream = upstream;
        piece.rank = rank;
        piece.number = number;
        piece.kept = kept;
        piece.end = path.size() - 1;
        return piece;
    };
    // Writes a channel cell's values (its owner's).
    auto paint = [&](std::size_t index, uint8_t direction, float bed, std::size_t upstream) {
        out.channel[index] = 1;
        out.direction[index] = direction;
        out.bed_m[index] = bed;
        out.discharge_m3s[index] = network.discharge_m3s[upstream];
        out.januaryflow[index] = network.januaryflow[upstream];
        out.julyflow[index] = network.julyflow[upstream];
        out.globalcell[index] = static_cast<int32_t>(upstream);
    };

    // Whether a global cell lies in the window's kept box (the link counts).
    auto keptcell = [&](int i, int j) {
        const int column = ((i - window.column0) % width + width) % width;
        return column < window.columns && j >= window.row0 && j < window.row0 + window.rows;
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
            std::size_t receiver = nocell;
            if (!drawnlink(network, upstream, parameters.pinarea_km2, receiver))
                continue;
            const int direction = network.direction[upstream];
            const int ox = directionoffsetx(direction);
            const int oy = directionoffsety(direction);

            // The two legs (three with a transit block), each in its own cell's block
            // coordinates, then moved to the blocks' unwrapped absolute positions.
            const CellGeometry from = cellgeometry(k, network, upstream, parameters, outlet);
            const CellGeometry to = cellgeometry(k, network, receiver, parameters, outlet);
            const LinkEnds ends = linkends(k, network, upstream, receiver, ox, oy, parameters, outlet);
            if (ox != 0 && oy != 0 && keptcell(i, j))
            {
                out.diagonallinks++;
                out.cornerlinks += ends.transit == 0 ? 1 : 0;
            }
            const double displacement = meanderdisplacement(parameters, k, network.discharge_m3s[upstream]);
            const Point flow = unitvector({ static_cast<double>(ox), static_cast<double>(oy) });
            const std::uint64_t outhash = deterministiccombine(parameters.seed, legsalt, static_cast<std::int64_t>(upstream),
                                                               static_cast<std::int64_t>(receiver), 0);
            const std::uint64_t inhash = deterministiccombine(parameters.seed, legsalt, static_cast<std::int64_t>(upstream),
                                                              static_cast<std::int64_t>(receiver), 1);
            const double roughness = roughnessof(upstream);
            const double persistence = parameters.meanderroughnesspersistence;
            std::vector<FineCell> path = legpath(k, from.junction, ends.exit, displacement, outhash, from.tangent, flow,
                                                 roughness, persistence);
            const std::vector<FineCell> inflow = legpath(k, ends.entry, to.junction, displacement, inhash, flow,
                                                         to.mainupstream == upstream ? to.tangent : Point{}, roughness,
                                                         persistence);
            const std::int64_t originx = static_cast<std::int64_t>(k) * i;
            const std::int64_t originy = static_cast<std::int64_t>(k) * j;
            for (FineCell& cell : path)
            {
                cell.x += originx;
                cell.y += originy;
            }
            if (ends.transit != 0)
            {
                // S3.10: the transit leg across the transit block, from the first edge
                // point to the second, along the link's direction at both ends.
                const std::uint64_t transithash = deterministiccombine(
                    parameters.seed, legsalt, static_cast<std::int64_t>(upstream), static_cast<std::int64_t>(receiver), 2);
                const std::int64_t tx = ends.transit == 1 ? ox : 0;
                const std::int64_t ty = ends.transit == 1 ? 0 : oy;
                for (const FineCell& cell :
                     legpath(k, ends.transitentry, ends.transitexit, displacement, transithash, flow, flow, roughness,
                             persistence))
                    path.push_back({ cell.x + originx + static_cast<std::int64_t>(k) * tx,
                                     cell.y + originy + static_cast<std::int64_t>(k) * ty });
            }
            for (const FineCell& cell : inflow)
                path.push_back({ cell.x + originx + static_cast<std::int64_t>(k) * ox,
                                 cell.y + originy + static_cast<std::int64_t>(k) * oy });

            const bool cut = cutatsea(path);
            // A link into a sea receiver ends at the sea level even where the fine sea
            // leaves the receiver's junction land (without a fine sea the cut always
            // happens).
            const bool intosea = cut || network.sea[receiver] != 0;

            // Bed (header). The global routing's receivers lie on the filled surface and
            // land stands at or above the sea level, so an end above the upstream cell,
            // or a filled elevation that is not a number, is bad input. The receiver is
            // checked whether or not the link was cut at the fine coast (junctionbed reads
            // its filled elevation as the upstream junction's lower bound), and the sea
            // level too when the link was cut, so the check does not depend on the window.
            const float filled = network.filledelevation_m[upstream];
            const float receiverend = network.sea[receiver] != 0 ? sealevel : network.filledelevation_m[receiver];
            if (!std::isfinite(filled) || !std::isfinite(receiverend))
                throw std::invalid_argument("rasterisechannels: a drawn link's filled elevation is not finite (cell " +
                                            std::to_string(upstream) + ")");
            if (receiverend > filled || (intosea && sealevel > filled))
                throw std::invalid_argument("rasterisechannels: a receiver (or the sea level, for a link that ends "
                                            "in the sea) stands above its upstream cell's filled elevation (cell " +
                                            std::to_string(upstream) + ")");
            const float top = junctionbed(upstream);
            const float end = intosea ? sealevel : junctionbed(receiver);
            std::vector<double> along(path.size(), 0.0);
            for (std::size_t m = 1; m < path.size(); m++)
            {
                const bool diagonal = path[m].x != path[m - 1].x && path[m].y != path[m - 1].y;
                along[m] = along[m - 1] + (diagonal ? std::sqrt(2.0) : 1.0);
            }
            const double total = along.back();
            const int32_t rank = outletdistance(network, upstream, distance);
            // The floor of a cell before the last, which is a land cell where a fine
            // sea is given: the next float above the sea level, but never above top (a
            // coastal cell exactly at the sea level keeps its level) and never below the
            // end value. end <= floor <= top, so the bed stays in [end, top] and the
            // clamp keeps the running minimum's order.
            const float bedfloor = std::max(end, std::min(top, abovesea));
            float running = top;

            Piece& piece = startpiece(path, upstream, rank, 0, keptcell(i, j));
            for (std::size_t m = 0; m < path.size(); m++)
            {
                // The bed outside the domain too (as this window knows it: no start surface
                // there), where it competes for ownership but is not painted.
                const std::size_t index = piece.index[m];
                if (usestart && index != nocell)
                    running = std::min(running, start_m[index]);  // a NaN start value lowers nothing
                const bool lastcell = m + 1 == path.size();
                float bed = 0.0f;
                if (lastcell)
                    bed = end;  // a path of one cell is cut at its first cell: a fine sea junction
                else if (m == 0)
                    bed = top;
                else if (usestart)
                    bed = std::max(running, bedfloor);
                else
                    bed = static_cast<float>(top + (static_cast<double>(end) - top) * along[m] / total);
                piece.bed[m] = bed;
            }
        }

    // Springs of the source cells within linkreach of the domain (they lie in their own cell's block).
    if (parameters.springspercell > 0)
        for (int j = rowfirst; j <= rowlast; j++)
            for (int i : columns)
            {
                const std::size_t cell = static_cast<std::size_t>(j) * width + i;
                std::size_t receiver = nocell;
                if (!drawnlink(network, cell, parameters.pinarea_km2, receiver))
                    continue;
                const CellGeometry geometry = cellgeometry(k, network, cell, parameters, outlet);
                if (geometry.mainupstream != nocell)
                    continue;  // not a source
                const int32_t rank = outletdistance(network, cell, distance);
                const float junction = junctionbed(cell);
                const double displacement = meanderdisplacement(parameters, k, network.discharge_m3s[cell]);
                Point upstreamward{ -geometry.tangent.x, -geometry.tangent.y };
                if (iszero(upstreamward))
                    upstreamward = unitvector({ -static_cast<double>(directionoffsetx(network.direction[cell])),
                                                -static_cast<double>(directionoffsety(network.direction[cell])) });
                const double heading = std::atan2(upstreamward.y, upstreamward.x);
                const Point centre = middle(geometry.junction);
                const double margin = std::min(0.5 * k, std::max(0.5, parameters.edgemargin * k));
                for (int spring = 1; spring <= parameters.springspercell; spring++)
                {
                    const std::uint64_t hash =
                        deterministiccombine(parameters.seed, springsalt, static_cast<std::int64_t>(cell), spring);
                    const double angle = heading + springspread * signedunit(hash);
                    const std::uint64_t radiushash = deterministicmix64(hash);
                    const double radius = k * (springnear + (springfar - springnear) * unitinterval(radiushash));
                    const FineCell start{
                        blockindex(std::clamp(centre.x + radius * std::cos(angle), margin, k - margin), k),
                        blockindex(std::clamp(centre.y + radius * std::sin(angle), margin, k - margin), k) };
                    std::vector<FineCell> path =
                        legpath(k, start, geometry.junction, displacement, deterministicmix64(radiushash), Point{},
                                Point{}, roughnessof(cell), parameters.meanderroughnesspersistence);
                    if (path.size() < 2)
                        continue;
                    for (FineCell& fine : path)
                    {
                        fine.x += static_cast<std::int64_t>(k) * i;
                        fine.y += static_cast<std::int64_t>(k) * j;
                    }
                    const bool cut = cutatsea(path);
                    const float end = cut ? sealevel : junction;
                    float top = junction;
                    if (usestart)
                        for (const FineCell& fine : path)
                        {
                            std::size_t index = 0;
                            if (domainindex(fine, index) && std::isfinite(start_m[index]))
                                top = std::max(top, start_m[index]);
                        }
                    const float bedfloor = std::max(end, std::min(top, abovesea));
                    float running = top;
                    Piece& piece = startpiece(path, cell, rank, spring, keptcell(i, j));
                    for (std::size_t m = 0; m < path.size(); m++)
                    {
                        const std::size_t index = piece.index[m];
                        if (usestart && index != nocell)
                            running = std::min(running, start_m[index]);
                        const bool lastcell = m + 1 == path.size();
                        piece.bed[m] = lastcell ? end : (usestart ? std::max(running, bedfloor) : top);
                    }
                }
            }

    // Ownership on the full paths (the ownership rules above), on every cell a drawn piece
    // crosses (outside the domain among the drawn pieces, with the beds this window knows; not
    // painted there): the lowest bed, then the outlet distance, the upstream cell's index and the
    // piece number (a strict order, so the owners do not depend on the order of the pieces).
    const std::size_t slots = finecells + outsideslots.size();
    auto beats = [&](const Piece& a, std::size_t ma, const Piece& b, std::size_t mb) {
        if (a.bed[ma] != b.bed[mb])
            return a.bed[ma] < b.bed[mb];
        if (a.rank != b.rank)
            return a.rank < b.rank;
        if (a.upstream != b.upstream)
            return a.upstream < b.upstream;
        return a.number < b.number;
    };
    std::vector<int32_t> owner(slots, -1);
    std::vector<std::uint32_t> ownerat(slots, 0);
    for (std::size_t p = 0; p < pieces.size(); p++)
    {
        const Piece& piece = pieces[p];
        for (std::size_t m = 0; m < piece.slot.size(); m++)
        {
            const std::size_t slot = piece.slot[m];
            if (owner[slot] < 0 || beats(piece, m, pieces[static_cast<std::size_t>(owner[slot])], ownerat[slot]))
            {
                owner[slot] = static_cast<int32_t>(p);
                ownerat[slot] = static_cast<std::uint32_t>(m);
            }
        }
    }
    // S3.11, the confluences (header), without iteration: a piece's prefix runs from its first cell
    // to its first lost cell (included; its last cell if it loses none), and its suffix from the
    // earliest of its own cells after its prefix that is another piece's first lost cell to its
    // last cell (a prefix end on its own prefix needs none: the water runs on along the prefix,
    // and a suffix from there would paint the piece past its first lost cell again). Painted
    // cells carry their owner's values; a path cell in no prefix or suffix is not painted.
    for (std::size_t p = 0; p < pieces.size(); p++)
    {
        Piece& piece = pieces[p];
        piece.end = piece.slot.size() - 1;
        piece.lost = false;
        for (std::size_t m = 0; m < piece.slot.size(); m++)
            if (owner[piece.slot[m]] != static_cast<int32_t>(p))
            {
                piece.end = m;
                piece.lost = true;
                break;
            }
    }
    for (const Piece& piece : pieces)
        if (piece.lost)
        {
            const std::size_t slot = piece.slot[piece.end];
            Piece& taker = pieces[static_cast<std::size_t>(owner[slot])];
            if (taker.lost && ownerat[slot] > taker.end)
                taker.suffix = std::min<std::size_t>(taker.suffix, ownerat[slot]);
        }
    std::vector<uint8_t> painted(slots, 0);
    for (const Piece& piece : pieces)
    {
        for (std::size_t m = 0; m <= piece.end; m++)
            painted[piece.slot[m]] = 1;
        for (std::size_t m = piece.suffix; m < piece.slot.size(); m++)
            painted[piece.slot[m]] = 1;
    }
    // S3.14, the painted network connected (header): a painted cell's chain runs on along its owner's
    // path. Where the owner's next cell is not painted (a suffix ran over a cell that a third piece
    // owns past that piece's own prefix end), the chain is followed over the unpainted cells, each
    // painted, until it meets a painted cell or reaches a piece's last cell. The directions hold no
    // loop, so every walk ends; cells are only added, so a painting without dry tails is unchanged,
    // and the result is the set of cells the painted cells' chains reach, whatever the order.
    for (std::size_t slot = 0; slot < slots; slot++)
    {
        if (painted[slot] == 0)
            continue;
        for (std::size_t at = slot;;)
        {
            const Piece& piece = pieces[static_cast<std::size_t>(owner[at])];
            const std::size_t m = ownerat[at];
            if (m + 1 >= piece.slot.size() || painted[piece.slot[m + 1]] != 0)
                break;
            at = piece.slot[m + 1];
            painted[at] = 1;
        }
    }
    for (std::size_t index = 0; index < finecells; index++)
    {
        if (owner[index] < 0 || painted[index] == 0)
            continue;
        const Piece& piece = pieces[static_cast<std::size_t>(owner[index])];
        const std::size_t m = ownerat[index];
        paint(index, piece.direction[m], piece.bed[m], piece.upstream);
    }
    // The counters (kept pieces): a prefix that ends before the piece's last cell; a cell the piece
    // owns and paints whose next path cell is not painted (its chain leaves the pinned network).
    for (std::size_t p = 0; p < pieces.size(); p++)
    {
        const Piece& piece = pieces[p];
        if (!piece.kept)
            continue;
        out.mergedpieces += piece.end + 1 < piece.slot.size() ? 1 : 0;
        for (std::size_t m = 0; m + 1 < piece.slot.size(); m++)
            if (painted[piece.slot[m]] != 0 && owner[piece.slot[m]] == static_cast<int32_t>(p) &&
                painted[piece.slot[m + 1]] == 0)
            {
                out.drytails++;
                break;
            }
    }

    // S3.9 (pinsinks): the junction of every land sink inside the domain that no link reached, so that
    // every land block holds a pinned cell (its outflow link's junction, else this one).
    if (parameters.pinsinks != 0)
        for (int jj = 0; jj < window.domainrows; jj++)
            for (int ii = 0; ii < window.domaincolumns; ii++)
            {
                const int j = window.domainrow0 + jj;
                const int i = (window.domaincolumn0 + ii) % width;
                const std::size_t cell = static_cast<std::size_t>(j) * width + i;
                if (network.sea[cell] != 0 || globalreceiver(network, cell) != cell ||
                    !std::isfinite(network.filledelevation_m[cell]))
                    continue;
                const FineCell offset = junctionoffset(k, cell, parameters);
                const FineCell junction{ static_cast<std::int64_t>(k) * i + offset.x, static_cast<std::int64_t>(k) * j + offset.y };
                std::size_t index = 0;
                if (!domainindex(junction, index) || out.channel[index] != 0 || (usesea && finesea[index] != 0))
                    continue;
                paint(index, 0, junctionbed(cell), cell);
                const int x = static_cast<int>(index % static_cast<std::size_t>(domainwidth));
                const int y = static_cast<int>(index / static_cast<std::size_t>(domainwidth));
                const int x0 = window.keptfinecolumnoffset(), y0 = window.keptfinerowoffset();
                if (x >= x0 && x < x0 + window.finewidth() && y >= y0 && y < y0 + window.fineheight())
                    out.sinkjunctions++;
            }
    return out;
}

void channelslevel(const LevelContext& context, LevelState& state)
{
    // Multiscale plan S3.0: regional_pass.cpp's step 3 moved here unchanged.
    // 3. The pinned channels: they end at their first fine sea cell, and their beds
    // follow the start surface's running minimum (channel_network.hpp). Beds are
    // pinned on land; a link's last cell in the sea keeps its own (seabed) elevation,
    // since the solve pins every sea cell at the sea level anyway. R9 (regional_valleys.hpp):
    // with the incision active the bed law reads the network's incised filled elevations
    // and the start surface lowered by the upsampled incision, so the channels cut below
    // the land by the world's uplift over the elapsed time; inactive, the inputs as they are.
    const RegionalGlobalInputs& inputs = context.inputs;
    const bool incised = context.incision.active;
    const std::vector<uint8_t>& sea = state.sea;
    std::vector<float>& surface = state.surface_m;
    const std::size_t count = surface.size();
    // S3.9: with bounded catchments in the refinement loop every land sink's junction is pinned.
    ChannelParameters channelparameters = context.parameters.channels;
    channelparameters.pinsinks = context.refinement != nullptr && context.parameters.relief.catchmentreach > 0 ? 1 : 0;
    // S3.14: the legs' roughness grows with the relief the regional noise amplitude reads
    // (regional_surface.cpp's noiserelief: the world's relief amplitude with noise_relief_field 1,
    // else the local relief, the largest absolute difference between a global cell's elevation and
    // its 8 neighbours', columns wrapped, no rows beyond the poles), not needed without roughness or gain.
    std::vector<float> relief;
    if (channelparameters.meanderroughness > 0.0 && channelparameters.meanderroughnessreliefgain > 0.0)
    {
        if (context.parameters.noiserelieffield == 1 && !inputs.reliefamplitude_m.empty())
            relief = inputs.reliefamplitude_m;
        else
        {
            const int width = inputs.width, height = inputs.height;
            const std::vector<float>& elevation = inputs.elevation_m;
            relief.assign(elevation.size(), 0.0f);
            for (int j = 0; j < height; j++)
                for (int i = 0; i < width; i++)
                {
                    const float own = elevation[static_cast<std::size_t>(j) * width + i];
                    float largest = 0.0f;
                    for (int dj = -1; dj <= 1; dj++)
                    {
                        const int row = j + dj;
                        if (row < 0 || row >= height)
                            continue;
                        for (int di = -1; di <= 1; di++)
                        {
                            const int column = ((i + di) % width + width) % width;
                            largest = std::max(largest,
                                               std::fabs(elevation[static_cast<std::size_t>(row) * width + column] - own));
                        }
                    }
                    relief[static_cast<std::size_t>(j) * width + i] = largest;
                }
        }
    }
    state.channels = rasterisechannels(context.window, incised ? state.incisednetwork : inputs.network,
                                       channelparameters, inputs.sealevel,
                                       incised ? state.bedsurface_m : surface, sea, relief);
    tributarieslevel(context, state);
    const FineChannels& channels = state.channels;
    for (std::size_t cell = 0; cell < count; cell++)
        if (channels.channel[cell] != 0 && sea[cell] == 0)
            surface[cell] = channels.bed_m[cell];
    state.bed_m = channels.bed_m;
    state.sinkjunctions = channels.sinkjunctions;
    state.diagonallinks = channels.diagonallinks;
    state.cornerlinks = channels.cornerlinks;
    state.mergedpieces = channels.mergedpieces;
    state.drytails = channels.drytails;
    recordlevelstep(state, LevelStep::channels);
}
}
