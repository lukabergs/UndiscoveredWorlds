#include "simulations/landscape/analytical_relief.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <functional>
#include <limits>
#include <map>
#include <queue>
#include <stdexcept>
#include <tuple>
#include <utility>

#include "simulations/landscape/sphere_metric.hpp"

// Erosion plan step E2c (docs/architecture/erosion-hydrology-plan.md, Phase E and
// section 6, "The solve"). With n = 1 and no uplift the stream power law
// dz/dt = -K A^m |grad z| is an advection of the start surface up the flow paths at
// the speed K Q^m (Tzathas et al. 2024, Eqs. 1, 5-9 with u = 0): z(x, t) = z0(D)
// with D the point whose response time is tau(x) - t, or the boundary's elevation
// when tau(x) <= t. Each step below is a pass over the cells in the flood's order
// (receivers before donors) and reads only the cell, its receiver path and its
// donors (in direction-code order), so the result does not depend on how ties in
// that order fall elsewhere in the patch.

namespace landscape
{
namespace
{
constexpr double cubicmetrespercubickilometre = 1.0e9;

struct Routing
{
    std::vector<int32_t> receiver;   // -1 on pinned cells
    std::vector<uint8_t> direction;  // 1..8 towards the receiver, 0 on pinned cells
    std::vector<int32_t> order;      // every cell once, each after its receiver
};

// The direction code from a cell to one of its neighbours, the lowest on a tie; 0 when none.
uint8_t directionto(const PatchMetric& metric, std::size_t cell, std::size_t target)
{
    for (int direction = 1; direction <= directioncount; direction++)
    {
        const std::size_t next = metric.neighbour(cell, direction);
        if (next != cell && next == target)
            return static_cast<uint8_t>(direction);
    }
    return 0;
}

// A diagonal step whose two corner cells (the cells orthogonally adjacent to both
// ends) are pinned passes between two cells of an 8-connected pinned line: refused,
// so that catchments stay on their side of the channel network.
bool crossespinnedcorner(const PatchMetric& metric, const std::vector<uint8_t>& pinned, std::size_t cell,
                         int direction)
{
    const int ox = directionoffsetx(direction), oy = directionoffsety(direction);
    if (ox == 0 || oy == 0)
        return false;
    const int x = static_cast<int>(cell % static_cast<std::size_t>(metric.width));
    const int y = static_cast<int>(cell / static_cast<std::size_t>(metric.width));
    return pinned[metric.index(x + ox, y)] != 0 && pinned[metric.index(x, y + oy)] != 0;
}

// S3.9's bounded catchments (AnalyticalReliefParameters::catchmentreach >= 1): per cell its
// catchment's label (the pinned cell it drains to; itself on a pinned cell) and the pinned cells with
// the unreached ones added (labelcatchments).
struct Catchments
{
    std::vector<int32_t> label;
    std::vector<uint8_t> pinned;
};

// The tie key of the bounded solve's queues: a bijective 64-bit mix (splitmix64's finaliser) of the
// absolute key, so that equal levels are taken in an order without a direction (the row-major key
// takes the north-west first, which stretches catchments south-east of their outlets).
std::uint64_t hashedkey(std::uint64_t value)
{
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ull;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebull;
    return value ^ (value >> 31);
}

// The labels (header, AnalyticalReliefParameters::catchmentreach r): one Priority-Flood per global
// cell g of the patch, seeded at the pinned cells inside g's block, over the unpinned cells whose
// global cell lies within r of g (Chebyshev, columns wrapped), as route's flood (levels, pit queue,
// the corner rule) with the hashed tie key; it gives each cell it reaches a filled level, the seed it
// came from and its parent. Each unpinned cell takes the seed with the smallest (filled level, hashed
// key of the seed) over the floods reaching it, so its label is a function of the cells within 2r
// global cells. A labelled cell whose parent chain in its winning flood leaves its label (a catchment
// cut into pieces by the floods' boxes), and a cell no flood reaches, is pinned at its start.
Catchments labelcatchments(const PatchMetric& metric, const std::vector<float>& surface, const AnalyticalReliefInput& input,
                           int reach)
{
    const std::size_t count = metric.size();
    const std::vector<uint8_t>& pinned = input.pinned;
    const int width = input.globalwidth;
    const auto keyof = [&](std::size_t cell) {
        return hashedkey(input.key.empty() ? static_cast<std::uint64_t>(cell) : input.key[cell]);
    };
    // The global cells (row, column) of the patch, each with its cells, in a canonical order.
    std::map<std::pair<int, int>, std::size_t> groupof;
    std::vector<std::vector<int32_t>> groups;
    std::vector<std::pair<int, int>> coordinates;
    for (std::size_t cell = 0; cell < count; cell++)
    {
        const std::pair<int, int> coordinate{ input.globalrow[cell], input.globalcolumn[cell] };
        const auto found = groupof.emplace(coordinate, groups.size());
        if (found.second)
        {
            groups.emplace_back();
            coordinates.push_back(coordinate);
        }
        groups[found.first->second].push_back(static_cast<int32_t>(cell));
    }

    Catchments catchments;
    catchments.pinned = pinned;
    catchments.label.assign(count, -1);
    std::vector<float> bestlevel(count, std::numeric_limits<float>::infinity());
    std::vector<std::uint64_t> besttie(count, 0);
    std::vector<int32_t> besthops(count, 0);
    std::vector<int32_t> bestparent(count, -1);
    for (std::size_t cell = 0; cell < count; cell++)
        if (pinned[cell] != 0)
            catchments.label[cell] = static_cast<int32_t>(cell);

    std::vector<int32_t> localof(count, -1);
    std::vector<uint8_t> boxgroup(groups.size(), 0);
    using Entry = std::tuple<float, std::uint64_t, int32_t>;
    for (std::size_t g = 0; g < groups.size(); g++)
    {
        bool seeded = false;
        for (int32_t cell : groups[g])
            seeded = seeded || pinned[static_cast<std::size_t>(cell)] != 0;
        if (!seeded)
            continue;
        // The box: the global cells within the reach, each once (columns wrap).
        std::vector<int32_t> box;
        std::vector<std::size_t> members;
        for (int dr = -reach; dr <= reach; dr++)
            for (int dc = -reach; dc <= reach; dc++)
            {
                const int column = ((coordinates[g].second + dc) % width + width) % width;
                const auto found = groupof.find({ coordinates[g].first + dr, column });
                if (found == groupof.end() || boxgroup[found->second] != 0)
                    continue;
                boxgroup[found->second] = 1;
                members.push_back(found->second);
                for (int32_t cell : groups[found->second])
                {
                    localof[static_cast<std::size_t>(cell)] = static_cast<int32_t>(box.size());
                    box.push_back(cell);
                }
            }
        const std::size_t size = box.size();
        std::vector<float> filled(size);
        std::vector<uint8_t> closed(size, 0);
        std::vector<int32_t> seed(size, -1), parent(size, -1), hops(size, 0);
        std::priority_queue<Entry, std::vector<Entry>, std::greater<Entry>> queue;
        for (int32_t cell : groups[g])
        {
            const std::size_t index = static_cast<std::size_t>(cell);
            if (pinned[index] == 0)
                continue;
            const int32_t local = localof[index];
            closed[static_cast<std::size_t>(local)] = 1;
            seed[static_cast<std::size_t>(local)] = cell;
            filled[static_cast<std::size_t>(local)] = surface[index];
            queue.push(Entry{ surface[index], keyof(index), local });
        }
        std::vector<int32_t> pit;
        std::size_t pithead = 0;
        while (pithead < pit.size() || !queue.empty())
        {
            std::size_t local;
            if (pithead < pit.size())
                local = static_cast<std::size_t>(pit[pithead++]);
            else
            {
                local = static_cast<std::size_t>(std::get<2>(queue.top()));
                queue.pop();
                pit.clear();
                pithead = 0;
            }
            const std::size_t cell = static_cast<std::size_t>(box[local]);
            for (int direction = 1; direction <= directioncount; direction++)
            {
                const std::size_t next = metric.neighbour(cell, direction);
                if (next == cell || pinned[next] != 0 || localof[next] < 0)
                    continue;
                const std::size_t nextlocal = static_cast<std::size_t>(localof[next]);
                if (closed[nextlocal] != 0 || crossespinnedcorner(metric, pinned, cell, direction))
                    continue;
                closed[nextlocal] = 1;
                seed[nextlocal] = seed[local];
                hops[nextlocal] = hops[local] + 1;
                parent[nextlocal] = static_cast<int32_t>(cell);
                if (surface[next] < filled[local])
                {
                    filled[nextlocal] = filled[local];
                    pit.push_back(static_cast<int32_t>(nextlocal));
                }
                else
                {
                    filled[nextlocal] = surface[next];
                    queue.push(Entry{ surface[next], keyof(next), static_cast<int32_t>(nextlocal) });
                }
            }
        }
        // Each reached unpinned cell keeps the better of its labels.
        for (std::size_t local = 0; local < size; local++)
        {
            const std::size_t cell = static_cast<std::size_t>(box[local]);
            if (pinned[cell] != 0 || closed[local] == 0)
                continue;
            const std::uint64_t tie = keyof(static_cast<std::size_t>(seed[local]));
            const bool better =
                catchments.label[cell] < 0 || filled[local] < bestlevel[cell] ||
                (filled[local] == bestlevel[cell] &&
                 (hops[local] < besthops[cell] || (hops[local] == besthops[cell] && tie < besttie[cell])));
            if (better)
            {
                bestlevel[cell] = filled[local];
                besthops[cell] = hops[local];
                besttie[cell] = tie;
                bestparent[cell] = parent[local];
                catchments.label[cell] = seed[local];
            }
        }
        for (int32_t cell : box)
            localof[static_cast<std::size_t>(cell)] = -1;
        for (std::size_t member : members)
            boxgroup[member] = 0;
    }

    // Connected: the parent chain in the winning flood keeps the label down to the seed (the chain's
    // cells won the same flood, since a seed lies in one block only). 0 unknown, 1 yes, 2 no.
    std::vector<uint8_t> connected(count, 0);
    std::vector<int32_t> path;
    for (std::size_t start = 0; start < count; start++)
    {
        if (pinned[start] != 0 || connected[start] != 0)
            continue;
        path.clear();
        std::size_t cell = start;
        uint8_t result = 2;
        for (;;)
        {
            if (catchments.label[cell] < 0)
                break;
            path.push_back(static_cast<int32_t>(cell));
            const std::size_t next = static_cast<std::size_t>(bestparent[cell]);
            if (catchments.label[next] != catchments.label[cell])
                break;
            if (pinned[next] != 0)
            {
                result = 1;
                break;
            }
            if (connected[next] != 0)
            {
                result = connected[next];
                break;
            }
            cell = next;
        }
        if (path.empty())
            connected[start] = 2;
        for (int32_t member : path)
            connected[static_cast<std::size_t>(member)] = result;
    }
    for (std::size_t cell = 0; cell < count; cell++)
        if (pinned[cell] == 0 && connected[cell] != 1)
        {
            catchments.pinned[cell] = 1;
            catchments.label[cell] = static_cast<int32_t>(cell);
        }
    return catchments;
}

// S3.12 (multiscale plan §3, the addendum of 2026-10-11): route's flood receivers resolved as
// resolveflatcells (flow_graph.cpp) resolves the global grid's (Barnes, Lehman and Mulla 2014), with
// metric distances and hashed ties. A step is allowed as in route: not a diagonal between two pinned
// corner cells of the input, and with catchments inside one label. A flat is a connected set (allowed
// steps) of the unpinned cells that take a flood receiver (`flooded`) at one filled elevation L; its low
// edges are the cells an allowed step away from it that are pinned or drain by steepest descent, at a
// filled elevation at most L (route's flood makes it exactly L). Two Dijkstra distances in metres over
// each flat's cells (steps of distance_m; the heap ties by the hashed key): towards, from the low edges
// (a flat cell next to one starts at its step to the nearest), and away, from the flat cells next to
// higher terrain (a neighbour outside the flat higher on the filled surface; on the surface across a
// pinned corner or a label, where the filled surface is another catchment's). The increment
// 2 towards + (height - away), height the flat's largest away (the term 0 in a flat without higher
// terrain), decreases towards the low edges, whose value is 0; a flat cell drains by steepest descent of
// the increment over the step's distance to a cell of its flat or a low edge, a tie to the neighbour
// with the smaller hashed key, so no raster order prints a direction. The level sets are octagons (an
// 8-neighbour Dijkstra), not the Chebyshev squares of breadth-first distances; on them a cardinal and a
// diagonal step off the axes drop the towards term equally per metre, so the away term or rounding
// decides there, and in a long flat Barnes' 2 : 1 weights send the cells diagonally onto the flat's
// axis (a diagonal step drops 3 / sqrt 2 per cell size against 2). Every reached flat cell has a strictly lower
// neighbour: the predecessor p on its towards path (the low edge, for a seed) has towards(c) =
// towards(p) + d and |away(c) - away(p)| <= d (a shortest-path distance), so the increment drops by at
// least d. No loop: a resolved link lowers the increment strictly inside one flat, so its chain ends at
// a low edge, which is pinned or drains by steepest descent to a strictly lower filled elevation, and no
// resolved link points to an unreached cell; route's argument covers the other links. A flat cell the
// towards distance does not reach keeps its flood receiver (none should: the flood chain of a flat cell
// stays at L inside the flat until it meets a low edge). The resolution reads the flat and the cells
// next to it, which the flood's sub-basin or label holds, so the seam property is unchanged. Returns
// whether a receiver changed.
bool resolveflats(const PatchMetric& metric, const std::vector<float>& surface, const std::vector<float>& filled,
                  const std::vector<uint8_t>& inputpinned, const std::vector<std::uint64_t>& key,
                  const Catchments* catchments, const std::vector<uint8_t>& flooded, Routing& routing)
{
    const std::size_t count = metric.size();
    const int width = metric.width, height = metric.height;
    const auto tieof = [&](std::size_t cell) {
        return hashedkey(key.empty() ? static_cast<std::uint64_t>(cell) : key[cell]);
    };
    // The neighbour offsets and, per row, the steps' distances (PatchMetric::distance_m's values).
    std::ptrdiff_t offset[directioncount + 1] = {};
    for (int direction = 1; direction <= directioncount; direction++)
        offset[direction] = directionoffsetx(direction)
                            + static_cast<std::ptrdiff_t>(directionoffsety(direction)) * static_cast<std::ptrdiff_t>(width);
    std::vector<double> step(static_cast<std::size_t>(height) * (directioncount + 1), 0.0);
    for (int y = 0; y < height; y++)
        for (int direction = 1; direction <= directioncount; direction++)
        {
            const int x = directionoffsetx(direction) < 0 ? width - 1 : 0;  // a column with the neighbour inside
            step[static_cast<std::size_t>(y) * (directioncount + 1) + static_cast<std::size_t>(direction)]
                = metric.distance_m(metric.index(x, y), direction);
        }
    const auto distance = [&](std::size_t cell, int direction) {
        return step[(cell / static_cast<std::size_t>(width)) * (directioncount + 1) + static_cast<std::size_t>(direction)];
    };
    const auto neighbourof = [&](std::size_t cell, int direction) {
        return static_cast<std::size_t>(static_cast<std::ptrdiff_t>(cell) + offset[direction]);
    };

    // The flats: their cells (members) contiguous, flat by flat; per member the flat, the directions
    // whose neighbour lies inside the patch and those of an allowed step.
    std::vector<int32_t> slot(count, -1);
    std::vector<int32_t> members, flatof, flatstart;
    std::vector<uint8_t> insidemask, allowedmask;
    for (std::size_t seed = 0; seed < count; seed++)
    {
        if (flooded[seed] == 0 || slot[seed] >= 0)
            continue;
        const int32_t id = static_cast<int32_t>(flatstart.size());
        flatstart.push_back(static_cast<int32_t>(members.size()));
        slot[seed] = static_cast<int32_t>(members.size());
        members.push_back(static_cast<int32_t>(seed));
        for (std::size_t head = static_cast<std::size_t>(flatstart.back()); head < members.size(); head++)
        {
            const std::size_t cell = static_cast<std::size_t>(members[head]);
            const int x = static_cast<int>(cell % static_cast<std::size_t>(width));
            const int y = static_cast<int>(cell / static_cast<std::size_t>(width));
            uint8_t inside = 0, allowed = 0;
            for (int direction = 1; direction <= directioncount; direction++)
            {
                const int ox = directionoffsetx(direction), oy = directionoffsety(direction);
                if (x + ox < 0 || x + ox >= width || y + oy < 0 || y + oy >= height)
                    continue;
                const uint8_t bit = static_cast<uint8_t>(1u << (direction - 1));
                inside |= bit;
                const std::size_t next = neighbourof(cell, direction);
                // crossespinnedcorner and the label, as route's steps
                if (ox != 0 && oy != 0 && inputpinned[cell + static_cast<std::size_t>(static_cast<std::ptrdiff_t>(ox))] != 0
                    && inputpinned[static_cast<std::size_t>(static_cast<std::ptrdiff_t>(cell) + offset[oy < 0 ? 1 : 5])] != 0)
                    continue;
                if (catchments != nullptr && catchments->label[next] != catchments->label[cell])
                    continue;
                allowed |= bit;
                if (flooded[next] == 0 || slot[next] >= 0 || filled[next] != filled[cell])
                    continue;
                slot[next] = static_cast<int32_t>(members.size());
                members.push_back(static_cast<int32_t>(next));
            }
            insidemask.push_back(inside);
            allowedmask.push_back(allowed);
            flatof.push_back(id);
        }
    }
    if (members.empty())
        return false;
    flatstart.push_back(static_cast<int32_t>(members.size()));
    const std::size_t size = members.size();
    const auto inflat = [&](std::size_t next, int32_t id) {
        return slot[next] >= 0 && flatof[static_cast<std::size_t>(slot[next])] == id;
    };
    // A low edge: outside the flats (pinned or draining by steepest descent), not higher.
    const auto lowedge = [&](std::size_t cell, std::size_t next) { return slot[next] < 0 && filled[next] <= filled[cell]; };

    const double infinity = std::numeric_limits<double>::infinity();
    using Entry = std::tuple<double, std::uint64_t, int32_t>;  // (distance, hashed key, member)
    std::vector<Entry> heap;
    // Dijkstra over one flat's members from the seeds in `heap` (their distances already set).
    const auto spread = [&](std::vector<double>& reached, int32_t id) {
        std::make_heap(heap.begin(), heap.end(), std::greater<Entry>());
        while (!heap.empty())
        {
            std::pop_heap(heap.begin(), heap.end(), std::greater<Entry>());
            const double at = std::get<0>(heap.back());
            const std::size_t member = static_cast<std::size_t>(std::get<2>(heap.back()));
            heap.pop_back();
            if (at > reached[member])
                continue;  // superseded
            const std::size_t cell = static_cast<std::size_t>(members[member]);
            for (int direction = 1; direction <= directioncount; direction++)
            {
                if ((allowedmask[member] & (1u << (direction - 1))) == 0)
                    continue;
                const std::size_t next = neighbourof(cell, direction);
                if (!inflat(next, id))
                    continue;
                const std::size_t other = static_cast<std::size_t>(slot[next]);
                const double candidate = at + distance(cell, direction);
                if (candidate < reached[other])
                {
                    reached[other] = candidate;
                    heap.push_back(Entry{ candidate, tieof(next), static_cast<int32_t>(other) });
                    std::push_heap(heap.begin(), heap.end(), std::greater<Entry>());
                }
            }
        }
    };

    std::vector<double> towards(size, infinity), away(size, infinity);
    bool changed = false;
    for (std::size_t flat = 0; flat + 1 < flatstart.size(); flat++)
    {
        const int32_t id = static_cast<int32_t>(flat);
        const std::size_t first = static_cast<std::size_t>(flatstart[flat]), last = static_cast<std::size_t>(flatstart[flat + 1]);
        // Towards the low edges: a member next to one starts at its step to the nearest.
        heap.clear();
        for (std::size_t member = first; member < last; member++)
        {
            const std::size_t cell = static_cast<std::size_t>(members[member]);
            for (int direction = 1; direction <= directioncount; direction++)
            {
                if ((allowedmask[member] & (1u << (direction - 1))) == 0)
                    continue;
                const std::size_t next = neighbourof(cell, direction);
                if (lowedge(cell, next))
                    towards[member] = std::min(towards[member], distance(cell, direction));
            }
            if (towards[member] < infinity)
                heap.push_back(Entry{ towards[member], tieof(cell), static_cast<int32_t>(member) });
        }
        spread(towards, id);
        // Away from higher terrain: 0 on the members next to it.
        heap.clear();
        for (std::size_t member = first; member < last; member++)
        {
            const std::size_t cell = static_cast<std::size_t>(members[member]);
            for (int direction = 1; direction <= directioncount; direction++)
            {
                const unsigned bit = 1u << (direction - 1);
                if ((insidemask[member] & bit) == 0)
                    continue;
                const std::size_t next = neighbourof(cell, direction);
                if (inflat(next, id))
                    continue;
                const float level = (allowedmask[member] & bit) != 0 ? filled[next] : surface[next];
                if (level > filled[cell])
                {
                    away[member] = 0.0;
                    heap.push_back(Entry{ 0.0, tieof(cell), static_cast<int32_t>(member) });
                    break;
                }
            }
        }
        spread(away, id);
        double top = 0.0;  // the flat's height
        for (std::size_t member = first; member < last; member++)
            if (away[member] < infinity)
                top = std::max(top, away[member]);
        // The increment, in place of towards (infinite where towards is).
        for (std::size_t member = first; member < last; member++)
            if (towards[member] < infinity)
                towards[member] = 2.0 * towards[member] + (away[member] < infinity ? top - away[member] : 0.0);
        const std::vector<double>& increment = towards;

        for (std::size_t member = first; member < last; member++)
        {
            const double own = increment[member];
            if (!(own < infinity))
                continue;  // not reached from a low edge: keeps its flood receiver
            const std::size_t cell = static_cast<std::size_t>(members[member]);
            int best = 0;
            double steepest = 0.0;
            std::uint64_t besttie = 0;
            for (int direction = 1; direction <= directioncount; direction++)
            {
                if ((allowedmask[member] & (1u << (direction - 1))) == 0)
                    continue;
                const std::size_t next = neighbourof(cell, direction);
                double value;
                if (inflat(next, id))
                    value = increment[static_cast<std::size_t>(slot[next])];  // infinite where unreached
                else if (lowedge(cell, next))
                    value = 0.0;
                else
                    continue;
                if (!(value < own))
                    continue;
                const double gradient = (own - value) / distance(cell, direction);
                const std::uint64_t tie = tieof(next);
                if (gradient > steepest || (best != 0 && gradient == steepest && tie < besttie))
                {
                    steepest = gradient;
                    best = direction;
                    besttie = tie;
                }
            }
            if (best == 0)
                continue;
            const int32_t receiver = static_cast<int32_t>(neighbourof(cell, best));
            changed = changed || routing.receiver[cell] != receiver;
            routing.receiver[cell] = receiver;
            routing.direction[cell] = static_cast<uint8_t>(best);
        }
    }
    return changed;
}

// Receivers on `surface`. A Priority-Flood seeded at the pinned cells (Barnes et al.
// 2014, Algorithm 2: a min-queue on (filled elevation, key) and a first-in first-out
// pit queue, as depressions.cpp, which ties by the local index instead) gives every
// cell a filled elevation, the neighbour it was reached from and the order in which
// it was taken; the key (the absolute fine index) makes that order the same in every
// window, and neighbours are visited in direction-code order. A cell drains
// by steepest descent (drop over metric distance, the lowest code on a tie) to a
// neighbour lower on both the surface and the filled surface; with none, by the flat
// resolution (S3.12, resolveflats: the flood's pits and flats drain along a metric
// gradient towards their low edges, ties by the hashed key), and to its flood receiver
// where that does not reach. The plain rule (steepest descent to any lower neighbour) loops:
// a rim cell descending into a pit whose flood receiver is that rim cell. Here a
// steepest-descent link lowers the filled elevation strictly and a flood link never
// raises it, while flood links alone form a tree towards the pinned cells and a resolved
// link ends at a low edge (resolveflats), so no loop exists; the flood takes cells in
// non-decreasing filled elevation, so its order lists every receiver before its donors,
// and when the resolution changes a receiver the order is rebuilt breadth-first from the
// pinned cells over the final receivers' donors (the solve reads only that property:
// its sums and walks are the same in any such order). On a surface without depressions the
// filled surface is the surface and every cell with a lower neighbour drains by
// steepest descent. Neither the flood nor the steepest descent takes a diagonal step
// between two pinned corner cells (crossespinnedcorner, on the input's pinned cells), so
// no water crosses an 8-connected channel; every cell is still reached, through orthogonal steps.
// With `catchments` (S3.9, fixed labels from labelcatchments): the flood steps only into cells of
// the seed's label and ties by the hashed key; the steepest descent takes only a neighbour of the
// cell's own label.
Routing route(const PatchMetric& metric, const std::vector<float>& surface, const std::vector<uint8_t>& inputpinned,
              const std::vector<std::uint64_t>& key, const Catchments* catchments)
{
    const std::size_t count = metric.size();
    Routing routing;
    routing.receiver.assign(count, -1);
    routing.direction.assign(count, 0);
    routing.order.reserve(count);
    const std::vector<uint8_t>& pinned = catchments != nullptr ? catchments->pinned : inputpinned;

    std::vector<float> filled = surface;
    std::vector<int32_t> floodreceiver(count, -1);
    const auto keyof = [&](std::size_t cell) {
        const std::uint64_t value = key.empty() ? static_cast<std::uint64_t>(cell) : key[cell];
        return catchments != nullptr ? hashedkey(value) : value;
    };
    // (level, key, cell): the cell decides only between equal keys, which the contract excludes.
    using Entry = std::tuple<float, std::uint64_t, int32_t>;
    std::priority_queue<Entry, std::vector<Entry>, std::greater<Entry>> queue;
    std::vector<uint8_t> closed(count, 0);
    for (std::size_t cell = 0; cell < count; cell++)
    {
        if (pinned[cell] == 0)
            continue;
        closed[cell] = 1;
        queue.push(Entry{ surface[cell], keyof(cell), static_cast<int32_t>(cell) });
    }
    std::vector<int32_t> pit;
    std::size_t pithead = 0;
    while (pithead < pit.size() || !queue.empty())
    {
        std::size_t cell;
        if (pithead < pit.size())
            cell = static_cast<std::size_t>(pit[pithead++]);
        else
        {
            cell = static_cast<std::size_t>(std::get<2>(queue.top()));
            queue.pop();
            pit.clear();
            pithead = 0;
        }
        routing.order.push_back(static_cast<int32_t>(cell));
        for (int direction = 1; direction <= directioncount; direction++)
        {
            const std::size_t next = metric.neighbour(cell, direction);
            if (next == cell || closed[next] != 0 || crossespinnedcorner(metric, inputpinned, cell, direction))
                continue;
            if (catchments != nullptr && catchments->label[next] != catchments->label[cell])
                continue;
            closed[next] = 1;
            floodreceiver[next] = static_cast<int32_t>(cell);
            if (surface[next] < filled[cell])
            {
                filled[next] = filled[cell];
                pit.push_back(static_cast<int32_t>(next));
            }
            else
                queue.push(Entry{ surface[next], keyof(next), static_cast<int32_t>(next) });
        }
    }
    std::vector<uint8_t> flooded(count, 0);  // takes a flood receiver: no lower neighbour on both surfaces
    for (std::size_t cell = 0; cell < count; cell++)
    {
        if (pinned[cell] != 0)
            continue;
        int best = 0;
        double steepest = 0.0;
        for (int direction = 1; direction <= directioncount; direction++)
        {
            const std::size_t next = metric.neighbour(cell, direction);
            if (next == cell || !(surface[next] < surface[cell]) || !(filled[next] < filled[cell])
                || crossespinnedcorner(metric, inputpinned, cell, direction))
                continue;
            if (catchments != nullptr && catchments->label[next] != catchments->label[cell])
                continue;
            const double gradient = (static_cast<double>(surface[cell]) - static_cast<double>(surface[next]))
                                    / metric.distance_m(cell, direction);
            if (gradient > steepest)
            {
                steepest = gradient;
                best = direction;
            }
        }
        if (best != 0)
        {
            routing.receiver[cell] = static_cast<int32_t>(metric.neighbour(cell, best));
            routing.direction[cell] = static_cast<uint8_t>(best);
        }
        else
        {
            routing.receiver[cell] = floodreceiver[cell];
            routing.direction[cell] = directionto(metric, cell, static_cast<std::size_t>(floodreceiver[cell]));
            flooded[cell] = 1;
        }
    }
    if (!resolveflats(metric, surface, filled, inputpinned, key, catchments, flooded, routing))
        return routing;
    // The order on the final receivers: the pinned cells, then breadth-first each listed cell's
    // donors (in increasing index, from a count of donors per cell); every cell is listed, since
    // every chain ends at a pinned cell.
    std::vector<int32_t> firstdonor(count + 1, 0);
    for (std::size_t cell = 0; cell < count; cell++)
        if (routing.receiver[cell] >= 0)
            firstdonor[static_cast<std::size_t>(routing.receiver[cell]) + 1]++;
    for (std::size_t cell = 0; cell < count; cell++)
        firstdonor[cell + 1] += firstdonor[cell];
    std::vector<int32_t> donors(static_cast<std::size_t>(firstdonor[count]));
    {
        std::vector<int32_t> cursor(firstdonor.begin(), firstdonor.end() - 1);
        for (std::size_t cell = 0; cell < count; cell++)
            if (routing.receiver[cell] >= 0)
                donors[static_cast<std::size_t>(cursor[static_cast<std::size_t>(routing.receiver[cell])]++)]
                    = static_cast<int32_t>(cell);
    }
    std::vector<int32_t>& order = routing.order;
    order.clear();
    for (std::size_t cell = 0; cell < count; cell++)
        if (pinned[cell] != 0)
            order.push_back(static_cast<int32_t>(cell));
    for (std::size_t head = 0; head < order.size(); head++)
    {
        const std::size_t cell = static_cast<std::size_t>(order[head]);
        order.insert(order.end(), donors.begin() + firstdonor[cell], donors.begin() + firstdonor[cell + 1]);
    }
    if (order.size() != count)
        throw std::logic_error("analytical relief: a receiver chain does not reach a pinned cell");
    return routing;
}

// The solution on fixed receivers: discharge, response time, the start surface
// advected up the flow paths, then the threshold slope. `pinned` is the input's, with S3.9's
// unreached cells added; `outlet` receives each cell's pinned cell.
std::vector<float> solve(const PatchMetric& metric, const AnalyticalReliefInput& input, const std::vector<uint8_t>& pinned,
                         const AnalyticalReliefParameters& parameters, const Routing& routing,
                         std::vector<double>& discharge, std::vector<int32_t>& outlet)
{
    const std::size_t count = metric.size();
    const std::vector<float>& start = input.start_m;
    const std::size_t width = static_cast<std::size_t>(metric.width);

    // Discharge: runoff x cell area, summed over each unpinned catchment, donors last
    // first; a cell adds its donors in direction-code order (a canonical sum).
    discharge.assign(count, 0.0);
    for (std::size_t k = routing.order.size(); k-- > 0;)
    {
        const std::size_t cell = static_cast<std::size_t>(routing.order[k]);
        if (pinned[cell] != 0)
            continue;
        double water = static_cast<double>(input.runoff_m_yr[cell]) * metric.rowarea_m2[cell / width];
        for (int direction = 1; direction <= directioncount; direction++)
        {
            const std::size_t next = metric.neighbour(cell, direction);
            if (next != cell && routing.receiver[next] == static_cast<int32_t>(cell))
                water += discharge[next];
        }
        discharge[cell] = water;
    }

    // Response time (Myr) down to the first pinned cell: tau = tau_receiver + L / (K Q^m),
    // Q in km^3/yr as stream_power.cpp (a negative sum counts as 0); a link without
    // flow (Q^m = 0) never responds (infinite time).
    const double infinity = std::numeric_limits<double>::infinity();
    std::vector<double> tau(count, 0.0);
    outlet.assign(count, -1);
    for (int32_t entry : routing.order)
    {
        const std::size_t cell = static_cast<std::size_t>(entry);
        if (pinned[cell] != 0)
        {
            outlet[cell] = entry;
            continue;
        }
        const std::size_t receiver = static_cast<std::size_t>(routing.receiver[cell]);
        const double length = metric.distance_m(cell, routing.direction[cell]);
        const double water = std::max(0.0, discharge[cell]) / cubicmetrespercubickilometre;
        const double speed = parameters.erodibility * std::pow(water, parameters.areaexponent);
        tau[cell] = tau[receiver] + (speed > 0.0 ? length / speed : infinity);
        outlet[cell] = outlet[receiver];
    }

    // The characteristic solution (Tzathas et al. 2024, Eqs. 8 and 9 with u = 0): with
    // tau <= t the pinned cell's elevation (with uplift_m_myr, plus U tau); otherwise the start elevation of the point
    // at response time s = tau - t downstream, linear in tau along the link a -> r
    // with tau(r) <= s < tau(a) (the paper's locally constant 1 / a). The walk covers
    // the links within time t of the cell. Neither value may exceed the cell's own
    // start: with u = 0 the erosion term K A^m |grad z| never raises a cell, while a
    // receiver path that climbs on the start surface (a depression's flood links, or
    // receivers taken from a later surface) would carry a higher z0 up to it.
    const double elapsed = parameters.elapsed_myr;
    const std::vector<float>& uplift = input.uplift_m_myr;
    std::vector<double> advected(count, 0.0);
    for (std::size_t cell = 0; cell < count; cell++)
    {
        const double own = static_cast<double>(start[cell]);
        if (pinned[cell] != 0)
        {
            advected[cell] = own;
            continue;
        }
        if (tau[cell] <= elapsed)
        {
            if (uplift.empty())
                advected[cell] = std::min(own, static_cast<double>(start[static_cast<std::size_t>(outlet[cell])]));
            else
                // The uplift term (R9, analytical_relief.hpp): the outlet's elevation plus
                // U tau, the cell's own U for the rate along its path.
                advected[cell] = std::min(own, static_cast<double>(start[static_cast<std::size_t>(outlet[cell])])
                                                   + std::max(0.0, static_cast<double>(uplift[cell])) * tau[cell]);
            continue;
        }
        const double target = tau[cell] - elapsed;
        std::size_t upper = cell;
        for (;;)
        {
            const std::size_t lower = static_cast<std::size_t>(routing.receiver[upper]);
            if (tau[lower] <= target)
            {
                const double zupper = static_cast<double>(start[upper]);
                const double zlower = static_cast<double>(start[lower]);
                if (target >= tau[upper])
                    advected[cell] = zupper;  // t = 0, or a link that never responds
                else
                    advected[cell] = zlower + (zupper - zlower) * (target - tau[lower]) / (tau[upper] - tau[lower]);
                break;
            }
            upper = lower;
        }
        advected[cell] = std::min(own, advected[cell]);
    }

    // Threshold slope (Tzathas et al. 2024, section 6): from the pinned cells upwards,
    // a cell above its receiver + S x L is lowered to it. The stored float is rounded
    // down when rounding to nearest would put it above the limit.
    std::vector<float> elevation(count, 0.0f);
    for (int32_t entry : routing.order)
    {
        const std::size_t cell = static_cast<std::size_t>(entry);
        if (pinned[cell] != 0)
        {
            elevation[cell] = start[cell];
            continue;
        }
        const std::size_t receiver = static_cast<std::size_t>(routing.receiver[cell]);
        const double limit = static_cast<double>(elevation[receiver])
                             + parameters.thresholdslope * metric.distance_m(cell, routing.direction[cell]);
        const double value = std::min(advected[cell], limit);
        float stored = static_cast<float>(value);
        if (static_cast<double>(stored) > limit)
            stored = std::nextafter(stored, -std::numeric_limits<float>::infinity());
        elevation[cell] = stored;
    }
    return elevation;
}
}

AnalyticalReliefResult solveanalyticalrelief(const PatchMetric& metric, const AnalyticalReliefInput& input,
                                             const AnalyticalReliefParameters& parameters)
{
    const std::size_t count = metric.size();
    const std::size_t rows = static_cast<std::size_t>(metric.height < 0 ? 0 : metric.height);
    if (metric.width < 0 || metric.height < 0 || metric.rowdx_m.size() != rows || metric.rowarea_m2.size() != rows
        || input.start_m.size() != count || input.pinned.size() != count || input.runoff_m_yr.size() != count
        || (!input.key.empty() && input.key.size() != count)
        || (!input.uplift_m_myr.empty() && input.uplift_m_myr.size() != count)
        || (parameters.catchmentreach > 0
            && (input.globalcolumn.size() != count || input.globalrow.size() != count || input.globalwidth < 1)))
        throw std::invalid_argument("solveanalyticalrelief: input sizes differ from the metric's");
    if (!(parameters.erodibility > 0.0) || !(parameters.areaexponent >= 0.0) || !(parameters.elapsed_myr >= 0.0)
        || !std::isfinite(parameters.elapsed_myr) || !(parameters.thresholdslope > 0.0) || parameters.maxiterations < 1
        || !(parameters.damping > 0.0) || !(parameters.damping <= 1.0) || parameters.catchmentreach < 0)
        throw std::invalid_argument("solveanalyticalrelief: a parameter is out of range");
    if (!input.key.empty())
    {
        std::vector<std::uint64_t> sorted = input.key;
        std::sort(sorted.begin(), sorted.end());
        if (std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end())
            throw std::invalid_argument("solveanalyticalrelief: two cells share a key");
    }
    bool anypinned = false;
    for (uint8_t flag : input.pinned)
        anypinned = anypinned || flag != 0;
    if (count > 0 && !anypinned)
        throw std::invalid_argument("solveanalyticalrelief: no cell is pinned");

    AnalyticalReliefResult result;
    result.elevation_m = input.start_m;
    result.direction.assign(count, 0);
    result.discharge_m3_yr.assign(count, 0.0);
    // The routing surface: the start, then the exponential moving average of the
    // solutions (Tzathas et al. 2024, section 5.1), held in double and read as float.
    std::vector<double> average(input.start_m.begin(), input.start_m.end());
    std::vector<float> routed = input.start_m;
    const double damping = parameters.damping;
    // Always maxiterations passes: stopping when no receiver changes would make the
    // pass count, and so every cell's result, depend on cells anywhere in the window.
    // S3.9: the bounded catchments, labelled once on the start surface.
    Catchments catchments;
    const Catchments* bounded = nullptr;
    if (parameters.catchmentreach > 0)
    {
        catchments = labelcatchments(metric, input.start_m, input, parameters.catchmentreach);
        bounded = &catchments;
    }
    Routing previous;
    for (int pass = 1; pass <= parameters.maxiterations; pass++)
    {
        Routing routing = route(metric, routed, input.pinned, input.key, bounded);
        result.iterations = pass;
        result.converged = pass > 1 && routing.receiver == previous.receiver;
        result.changedreceivers = 0;
        if (pass > 1)
            for (std::size_t cell = 0; cell < count; cell++)
                result.changedreceivers += routing.receiver[cell] != previous.receiver[cell] ? 1 : 0;
        // Every pass solves from the start surface on the receivers of the damped surface.
        result.elevation_m = solve(metric, input, bounded != nullptr ? catchments.pinned : input.pinned, parameters,
                                   routing, result.discharge_m3_yr, result.outlet);
        result.direction = routing.direction;
        previous = std::move(routing);
        if (pass == parameters.maxiterations)
            break;
        for (std::size_t cell = 0; cell < count; cell++)
        {
            average[cell] = damping * static_cast<double>(result.elevation_m[cell]) + (1.0 - damping) * average[cell];
            routed[cell] = static_cast<float>(average[cell]);
        }
    }
    if (bounded != nullptr)
    {
        result.unreached.assign(count, 0);
        for (std::size_t cell = 0; cell < count; cell++)
            if (catchments.pinned[cell] != input.pinned[cell])
            {
                result.unreached[cell] = 1;
                result.unreachedcells++;
            }
    }
    return result;
}
}
