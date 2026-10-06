#pragma once

// Drainage tree of a terrain on a SphereMetric: D8 receivers, the Braun-Willett
// order and the accumulated area and discharge
// (docs/architecture/erosion-hydrology-plan.md, section 3, Kernels).

#include <cstddef>
#include <cstdint>
#include <vector>

#include "simulations/landscape/sphere_metric.hpp"

namespace landscape
{
struct Depressions;

struct FlowGraph
{
    // Per cell. A land cell's receiver is the neighbour of steepest descent (drop
    // over centre distance) among the eight, a sea neighbour counting with
    // elevation = sea level. A land cell whose steepest descent is a sea neighbour
    // is an outlet: receiver = the cell itself, direction = the code towards that
    // sea neighbour, outletsea = that sea neighbour (so every land cell has a
    // direction 1..8, which the drained check of the evaluation relies on). A land
    // cell with no lower neighbour, or flooded by a depression (depressionid != 0),
    // takes the depression kernel's flood receiver (towards its spill), so every
    // land cell has a receiver and every chain ends at an outlet; a flood receiver
    // that is a sea cell makes the cell an outlet. A sea cell's receiver is its
    // lowest sea neighbour (the seabed's steepest descent), itself when none is
    // lower; sea cells have direction 0. A land cell the flood never reached (a
    // region without a sea cell) is its own receiver with direction 0 and a root
    // of `order`, like an outlet without a sea cell; the evaluation's drained
    // check reports such cells, as intended.
    std::vector<int32_t> receiver;        // index of the receiver; the cell itself for outlets and seabed pits
    std::vector<uint8_t> direction;       // 1..8 towards the receiver, or towards outletsea for an outlet; 0 for sea cells
    std::vector<float> distance_m;        // centre to centre, to outletsea for an outlet; 0 for a seabed pit
    std::vector<int32_t> outletsea;       // outlets: the sea cell the water enters; -1 elsewhere
    std::vector<uint8_t> land;            // 1 land, 0 sea (the input mask)
    // The Braun-Willett stack over the land cells: outlets first, every cell after
    // its receiver. Accumulation and the implicit solve run over it; accumulation
    // reversed (upstream cells first).
    std::vector<int32_t> order;
    // Accumulated over the tree: a cell includes every cell upstream of it. Sea cells: 0.
    std::vector<double> area_m2;          // drainage area
    std::vector<double> discharge_m3_yr;  // sum over the upstream cells of runoff x area

    std::size_t size() const { return receiver.size(); }
};

// How receivers are chosen (erosion plan, Phase C, step C3). The defaults are the
// deterministic D8 described above.
// - randomised: a land cell that drains to a lower neighbour (the steepest-descent
//   case above) picks, among the two neighbours whose directions bracket its local
//   gradient direction, one with probabilities that make the mean displacement
//   parallel to the gradient (the intent of Rho8, Fairfield and Leymarie 1991, on
//   the latitude-longitude metric, whose neighbour directions are not 45 degrees
//   apart away from the equator). The gradient g is D-infinity's (Tarboton 1997):
//   the steepest of the eight facets formed by the cell, a cardinal and the
//   adjacent diagonal neighbour, with the metric's distances; with p1, p2 the
//   offsets to the facet's two neighbours in metres, the second is drawn with
//   probability |p1 x g| / (|p1 x g| + |p2 x g|), so (1 - f) p1 + f p2 is parallel
//   to g. The receiver is lower on the filled surface or a sea neighbour (an
//   outlet, as above); when neither drawn neighbour is, the D8 choice stays. The
//   draw is a pure function of `seed` and the cell index
//   (core/deterministic_random.hpp), never the shared random stream, so a routing
//   repeated on the same terrain is identical.
// - resolveflats: the cells that follow their flood receiver (a flooded depression
//   or a flat without a lower neighbour) drain towards their flat's low edges and
//   away from higher terrain (Barnes, Lehman and Mulla 2014) instead of following
//   the flood order. Raised cells whose flood-receiver chains leave the fill at the
//   same cell form a flat whose only low edge is that cell, so a depression's cells
//   drain through its spill; a connected region of cells at their own terrain and
//   one filled elevation is a flat whose low edges are the adjacent land cells at
//   that elevation that drain downhill (a region without one keeps the flood
//   order). Depressions are not changed.
struct RoutingOptions
{
    bool randomised = false;
    std::uint64_t seed = 0;
    bool resolveflats = false;
};

// Receivers, directions, distances, outlets and the order; area and discharge are
// zero until accumulateflow. `sea` marks sea cells (1); `elevation` is the terrain
// of every cell; a sea neighbour of a land cell counts with elevation `sealevel`;
// `depressions` is the Priority-Flood result on the same inputs (depressions.hpp)
// and supplies the receivers of flooded and pit cells. Ties in the steepest descent
// go to the lowest direction code, so the result is deterministic.
FlowGraph buildflowgraph(const SphereMetric& metric, const std::vector<float>& elevation,
                         const std::vector<uint8_t>& sea, float sealevel, const Depressions& depressions,
                         const RoutingOptions& options = {});

// Area and discharge accumulated along `order` reversed. `runoff_m_yr` is the
// runoff depth per cell per year (sea cells ignored): discharge = sum over the
// cell and its upstream cells of runoff x area. Serial, so bit-identical for any thread count.
void accumulateflow(const SphereMetric& metric, const std::vector<float>& runoff_m_yr, FlowGraph& graph);
}
