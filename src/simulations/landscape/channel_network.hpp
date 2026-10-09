#pragma once

// Pinned channels of the regional pass (docs/architecture/erosion-hydrology-plan.md,
// Phase E, step E2b, section 6, "Seams of the solve", and Phase R, steps R1 and R3):
// the global drainage network drawn onto a window's fine grid after the legacy
// regional generator's rivers. Every global land cell whose drainage area reaches
// pinarea_km2 contributes one link to its receiver, drawn as two legs: the cell's
// outflow leg from its junction (a fine cell near its block's centre, a hash of the
// cell) to its exit point on the edge it shares with the receiver, then the
// receiver's inflow leg from the entry point (the fine cell across the edge) to the
// receiver's junction. A cardinal link crosses the shared edge at an offset drawn from
// a hash of the edge (both cells draw the same), a diagonal link through the shared
// corner. Each leg is a Catmull-Rom curve through its two ends and a midpoint displaced
// sideways (a hash of the leg), rasterised 8-connected inside its own global cell's
// k x k block, so a link occupies only its two cells' blocks. Every value is a pure
// function of the global cells, the seed and the absolute fine position, so a link that
// two windows share is drawn identically in both. Planet-free.

#include <cstdint>
#include <vector>

#include "simulations/landscape/regional_grid.hpp"

namespace landscape
{
// The global drainage network, W x H row-major (the planet's hydrosphere fields).
struct GlobalNetwork
{
    int width = 0;
    int height = 0;
    std::vector<uint8_t> direction;          // flow_direction: 1..8 towards the receiver, 0 for sea and sinks
    std::vector<uint8_t> sea;                // 1 sea
    std::vector<float> drainagearea_km2;
    std::vector<float> discharge_m3s;        // annual mean
    std::vector<float> filledelevation_m;    // the Priority-Flood surface, m above datum
    std::vector<int32_t> januaryflow;        // river_flow_january, legacy units
    std::vector<int32_t> julyflow;           // river_flow_july, legacy units
};

struct ChannelParameters
{
    double pinarea_km2 = 0.0;          // links of cells with a smaller drainage area are not pinned
    std::uint64_t seed = 0;            // derived from the world seed with a fixed salt
    // The geometry of Phase R, R3 (the legacy scheme, plan section 6).
    double meanderamplitude = 0.2;     // a leg's midpoint is displaced sideways by 0.5 to 1 times
                                       // meanderamplitude x k x min(1, (Q / meanderreference_m3s)^
                                       // meanderdischargeexponent) fine cells (Q the upstream cell's
                                       // discharge; exponent 0: no dependence), kept inside the block
    double meanderdischargeexponent = 0.0;
    double meanderreference_m3s = 1000.0;
    double junctionjitter = 0.25;      // the junction lies within this fraction of k of the block's centre
    double edgemargin = 0.1875;        // edge points stay this fraction of k away from the cell's corners
    int springspercell = 0;            // springs from hashed points into a source cell's junction (0: none)
};

// The pinned channels on the window's domain grid (domainfinewidth() x
// domainfineheight(), row-major).
//
// Geometry. The junction of global cell (i, j) is channeljunction(); a link's path is
// its outflow leg in the upstream cell's block followed by the receiver's inflow leg in
// the receiver's block, each an 8-connected chain of distinct fine cells (loops erased),
// so consecutive path cells are 8-neighbours and the step across the shared edge goes
// from the exit point to the entry point. The curve's tangent at the edge point is the
// link's direction, so the two legs join smoothly; at a junction, the outflow leg and
// the inflow leg of the main upstream cell (the largest drainage area among the drawn
// links into the cell, ties to the smaller index) share one tangent (from the main
// entry point towards the exit point), and other inflows join at their own angle. A
// source cell (a drawn link and no drawn link into it) also gets springspercell
// springs: legs from hashed points of its block (away from the block's edge, upstream
// of the junction) to its junction, carrying the cell's discharge and flows.
//
// Ownership. A junction belongs to its own cell's outgoing link; any other fine cell
// on several pieces (links and springs) to the piece with the lower bed there, equal
// beds to the link whose upstream cell is fewer links from its outlet on the global
// tree, then to the smaller absolute index of the upstream global cell, then to the
// link before its springs (spring 1 before spring 2 ...), so that the bed never rises
// along the directions and the directions hold no loop whatever the drainage areas
// (.cpp). A block holds only its cell's outflow leg, its springs and the inflow legs of
// the links into it, so only links adjacent on the global tree (the same receiver, or
// one's receiver the other's upstream cell) share fine cells, and a chain that runs onto
// another link reaches the same outlet; no capture across a divide can arise.
//
// Ends. A link ends at its first fine cell inside a global sea cell's block (its sea
// receiver's entry point) without `finesea`; with it, see below. A land sink's
// junction ends the links into it. The rest of a path is not drawn.
//
// Beds. The junction bed of a land cell with junction J is b = max(lo, min(F, s(J)))
// with F its filled elevation, s(J) the start surface at J (F where `start_m` is empty,
// J lies outside the domain or s(J) is not a number) and lo = max(e, min(F, next float
// above the sea level)), e the receiver's filled elevation, the sea level for a sea
// receiver, and minus infinity for a land sink. So b lies in [e, F] (the junction
// follows the start surface where the drop to the receiver allows: no pillar), and
// every link's receiver end lies at or below its upstream end: b(receiver) <= F(receiver)
// <= b(upstream). A link from u to c runs from top = b(u) at its first cell (u's
// junction) to end = b(c) at its last (c's junction), or to the sea level for a link
// that ends in the sea. Without `start_m` the bed between them falls linearly with the
// path length (fine cell edges, sqrt(2) per diagonal step). A spring's top is the
// largest of b(u) and the start values along it, its end b(u) (the sea level if it ends
// in the sea), and it is flat at b(u) without `start_m`.
//
// Throws std::invalid_argument for a scale below 1, a network not on the window's
// world grid, a domain off it, a negative or non-finite meanderamplitude,
// meanderdischargeexponent or junctionjitter, a meanderreference_m3s that is not
// positive and finite, an edgemargin outside 0 .. 0.5, a negative springspercell, a
// non-finite sea level, receivers that form a loop downstream of a link within reach of
// the domain (loops elsewhere go unseen), and, on a link within reach of the domain, a
// filled elevation that is not finite or a receiver (the sea level for a sea receiver;
// checked whether or not the link is cut at the fine coast) or, for a link that ends in
// the sea, the sea level above the upstream cell's filled elevation.
struct FineChannels
{
    std::vector<uint8_t> channel;        // 1 on a pinned channel cell
    std::vector<uint8_t> direction;      // 1..8 towards the owner piece's next cell (0 at a link's last cell:
                                         // its first fine sea cell, or a land sink's junction)
    std::vector<float> bed_m;            // pinned bed elevation, m above datum
    std::vector<float> discharge_m3s;    // the upstream global cell's discharge
    std::vector<int32_t> januaryflow;    // the upstream global cell's legacy flows
    std::vector<int32_t> julyflow;
    std::vector<int32_t> globalcell;     // the upstream global cell's index, -1 off the network
};

// Phase R, R1 (empty vectors keep the rules without them): `start_m` is the regional
// start surface and `finesea` the fine sea of the domain (domainfinewidth() x
// domainfineheight(), regional_surface.hpp), decided before the channels; each is used
// on its own.
// - With `finesea`, a link (or spring) ends at its first fine sea cell inside the domain
//   (the fine coast, not the global sea block), whose bed is the sea level exactly and
//   whose direction is 0, so no channel cell is pinned through the fine sea and every
//   channel cell before a link's last is fine land. A path cell outside the domain,
//   where the fine sea is unknown, ends the link in a global sea block as without it.
//   A link into a sea receiver ends at the sea level even if the fine sea leaves the
//   receiver's junction land. A land global cell whose own junction is fine sea has a
//   link of that one cell, at the sea level.
// - With `start_m`, the bed of a link's path cell m (not its first or last) is
//   max(r_m, f): r_m the minimum of top and of the start surface over the path cells
//   0 .. m that lie inside the domain (a NaN and cells outside it lower nothing),
//   f = max(end, min(top, next float above the sea level)). So the bed lies in
//   [end, top], follows the start surface's lows (no straight-profile trenches), never
//   rises along a link or along the directions (ownership, .cpp), and stands above the
//   sea level on every cell before the last unless top is the sea level itself. Springs
//   follow the same law with their own top and end.
// Seams: a channel cell's values depend on the pieces crossing it and on `start_m` and
// `finesea` at the cells of those pieces' paths (and their two junctions). A link lies
// inside its two global cells' blocks, so those cells lie at most one global cell nearer
// the domain's edge than the channel cell. The regional pass's start surface and fine
// sea agree between windows from 2 x meanpasses global cells inside the domain's edge
// (each block-mean pass spoils the residuals of the outer two), so its channels agree
// from 2 x meanpasses + 1 cells inside it.
// Throws std::invalid_argument also for a non-empty `start_m` or `finesea` whose size
// is not the domain's.
FineChannels rasterisechannels(const RegionalWindow& window, const GlobalNetwork& network,
                               const ChannelParameters& parameters, float sealevel,
                               const std::vector<float>& start_m = {}, const std::vector<uint8_t>& finesea = {});

// The junction of global cell (i, j), 0 <= i < worldwidth and 0 <= j < worldheight, at
// the window's scale: its absolute fine column k i + x and row k j + y, with x and y in
// 0 .. k - 1 drawn from a hash of the cell and the seed within junctionjitter x k of the
// block's centre (rasterisechannels' junctions; exposed for tests and probes).
struct FinePosition
{
    std::int64_t column = 0;
    std::int64_t row = 0;
};
FinePosition channeljunction(const RegionalWindow& window, const ChannelParameters& parameters, int i, int j);
}
