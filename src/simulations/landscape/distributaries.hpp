#pragma once

// Distributaries on prograded land (docs/architecture/erosion-hydrology-plan.md,
// step C5b): the river flow that enters land built by marine deposition splits over
// the lower neighbours, as distributaries do. Planet-free; hydrology.rivers passes
// the world's fields and writes the result into the legacy delta fields.
//
// The flows are the legacy units of the rivers (January and July), carried as
// doubles. A land cell off the prograded land whose flow_direction receiver is a
// prograded land cell hands its flows to it. Every prograded land cell splits what
// it carries over its neighbours that are lower on `surface` (the routing's filled
// elevation): prograded land at its surface, sea cells at the sea level, in
// proportion to the drop to each (a multiple-flow split). A cell with no such
// neighbour passes everything to its flow_direction receiver when that is prograded
// land no higher than the cell, or sea. Cells are processed in a topological order
// of these links (Kahn's algorithm; ready cells leave a first-in first-out queue
// seeded in increasing cell index), so no flow reaches a cell after it was split; on
// the filled surface the links have no cycle, and whatever enters the prograded land
// reaches the sea unless a cell has no link at all (its flow is then reported as
// stranded).
//
// Output, on prograded land only (sea cells carry nothing): a cell with a link onward
// (and processed: not left out by a cycle) is a delta cell where the mean of its
// rounded January and July flows reaches `minimumflow`; a stranded cell gets nothing. Its
// flows are what it carries and its direction (1..8, sphere_metric's codes) points
// to the neighbour that delivered its largest share, upstream: the layout of the
// legacy delta fields, whose branches form one tree rooted at the apex.

#include <cstdint>
#include <vector>

#include "simulations/landscape/sphere_metric.hpp"

namespace landscape
{
struct Distributaries
{
    std::vector<int16_t> direction;  // 0, or the code towards the largest supplier
    std::vector<int32_t> january;    // legacy units; 0 where direction is 0
    std::vector<int32_t> july;
    double inflowjanuary = 0.0;      // totals entering the prograded land
    double inflowjuly = 0.0;
    double seajanuary = 0.0;         // totals delivered to sea cells
    double seajuly = 0.0;
    double strandedjanuary = 0.0;    // totals left on cells without a link
    double strandedjuly = 0.0;
};

// Every vector has metric.size() elements in the grid's cell order: `sea` and
// `prograded` 0/1, `surface` in m on the datum of `sealevel`, `flowdirection` the
// routing's codes (0 none), the flows in legacy units. Throws std::invalid_argument
// on a size mismatch.
Distributaries splitdistributaries(const SphereMetric& metric, const std::vector<uint8_t>& sea,
                                   const std::vector<uint8_t>& prograded, const std::vector<float>& surface,
                                   double sealevel, const std::vector<uint8_t>& flowdirection,
                                   const std::vector<int32_t>& januaryflow, const std::vector<int32_t>& julyflow,
                                   double minimumflow);
}
