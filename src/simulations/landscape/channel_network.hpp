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
// corner or (multiscale plan S3.10, diagonaltransit) through a transit block: one of the
// two cells orthogonal to both, crossed by a third leg between two hashed edge points.
// Each leg is a Catmull-Rom curve through its two ends and a midpoint displaced
// sideways (a hash of the leg), roughened by hashed value-noise octaves along its normal
// (multiscale plan S3.11), rasterised 8-connected inside its own global cell's
// k x k block, so a link occupies only its two cells' blocks (and its transit block's).
// Every value is a pure
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
    // Multiscale plan S3.11: a leg's roughness, octaves i = 2, 3, ... of value noise along its
    // normal with amplitude meanderroughness x k x meanderroughnesspersistence^(i - 2) fine cells
    // (0: the smooth arc of S3.10), times (S3.14) 1 + meanderroughnessreliefgain x min(1, R /
    // meanderroughnessrelief_m), R the local relief of the piece's upstream global cell
    // (rasterisechannels' relief_m; gain 0: the plain roughness).
    double meanderroughness = 0.08;
    double meanderroughnesspersistence = 0.5;
    double meanderroughnessreliefgain = 1.0;
    double meanderroughnessrelief_m = 500.0;
    double junctionjitter = 0.5;       // the junction lies within this fraction of k of the block's centre
                                       // (0.5: anywhere in the block; R3: 0.25)
    double edgemargin = 0.0;           // edge points stay this fraction of k away from the cell's corners (R3: 0.1875;
                                       // 0: every fine cell of the edge equally likely at every k)
    int springspercell = 0;            // springs from hashed points into a source cell's junction (0: none)
    // Multiscale plan S3.10: 1 draws a diagonal link through a transit block (an orthogonal
    // neighbour of both cells that is land, filled between the link's two cells, not the receiver's
    // receiver nor the upstream cell's donor, and drains to the receiver's outlet on the global tree;
    // the corner where none does), 0 through the shared corner (R3).
    int diagonaltransit = 1;
    // S3.10: when neither side cell's filled elevation lies between the link's two cells', the one
    // missing that range by less is the transit block if the miss is at most this (m; 0: none).
    double transittolerance_m = 100.0;
    // Multiscale plan S3.9 (set by channelslevel in the refinement loop with bounded catchments,
    // regional.catchment_reach_global_cells >= 1; 0: none): 1 pins the junction of every land sink of
    // the domain (a land cell without a receiver on the grid) that no link reached and that is not fine
    // sea, as a link's last cell (direction 0) with the sink's junction bed, discharge, flows and global
    // cell, so that every land block holds a pinned cell.
    int pinsinks = 0;
};

// The pinned channels on the window's domain grid (domainfinewidth() x
// domainfineheight(), row-major).
//
// Geometry. The junction of global cell (i, j) is channeljunction(); a link's path is
// its outflow leg in the upstream cell's block followed by the receiver's inflow leg in
// the receiver's block, each an 8-connected chain of distinct fine cells (loops erased),
// so consecutive path cells are 8-neighbours and the step across the shared edge goes
// from the exit point to the entry point. A diagonal link u -> c with a transit block n
// (S3.10: n = u + (ox, 0) or u + (0, oy), land, filled between c and u, neither c's
// receiver nor u's donor, with the receiver's outlet on the global tree; both qualify: a
// hash of the link) runs u's outflow leg to a point on the edge
// u | n, a transit leg across n's block from the fine cell across that edge to a point on
// the edge n | c, then c's inflow leg from the fine cell across it, each edge point at an
// offset drawn from the link between edgemargin x k and (1 - edgemargin) x k, so both
// steps across an edge are 4-connected; without one (neither neighbour qualifies, or
// diagonaltransit 0) it steps diagonally through the shared corner. Roughness (S3.11): each
// sample of a leg's curve at its parameter t (0 .. 1/2 over the first Hermite segment, 1/2 .. 1
// over the second) moves along the curve's left normal (the chord's where the curve stands
// still) by d(t), the sum over octaves i = 2, 3, ... of a_i v_i(t): v_i a cosine-interpolated
// value noise with 2^i intervals over t, node values in [-1, 1) from a hash of (the leg's
// hash, i, node), 0 at t = 0 and t = 1, so legs still meet at their edge points and
// junctions; a_i = meanderroughness x k x meanderroughnesspersistence^(i - 2) x (S3.14)
// (1 + meanderroughnessreliefgain x min(1, R / meanderroughnessrelief_m)), R the local relief of
// the piece's upstream global cell (relief_m; 0 without it), the same for every leg of a link and
// for a cell's springs; octaves while chord / 2^i is at least 2 fine cells. A leg at a finer k is
// the same curve with more octaves; samples stay clamped inside the block, a gap between two samples' cells is filled
// by diagonal-first steps and the chain's loop erasure removes self-crossings. The curve's tangent
// at an edge point is the link's direction, so the legs join smoothly; at a junction, the outflow leg and
// the inflow leg of the main upstream cell (the largest drainage area among the drawn
// links into the cell, ties to the smaller index) share one tangent (from the main
// entry point towards the exit point), and other inflows join at their own angle. A
// source cell (a drawn link and no drawn link into it) also gets springspercell
// springs: legs from hashed points of its block (away from the block's edge, upstream
// of the junction) to its junction, carrying the cell's discharge and flows.
//
// Ownership. A fine cell on several pieces (links and springs) belongs to the piece
// with the lower bed there, equal beds to the link whose upstream cell is fewer links
// from its outlet on the global tree, then to the smaller absolute index of the upstream
// global cell, then to the link before its springs (spring 1 before spring 2 ...), so
// that the bed never rises along the directions and the directions hold no loop whatever
// the drainage areas (.cpp). A junction is reached at the same bed by its cell's outgoing
// link, the links into it and its springs, so it belongs to the outgoing link; a transit
// leg with a lower bed there takes it. A block holds only its cell's outflow leg, its
// springs, the inflow legs of the links into it and the transit legs through it, all
// draining to the block's outlet on the global tree, so a chain that runs onto another
// piece reaches the same outlet; no capture across a divide can arise.
//
// Confluences (multiscale plan S3.11), without iteration. Ownership and the directions are
// decided on the full paths as above (every cell a drawn piece crosses, also outside the domain,
// where only the drawn pieces compete with the beds this window knows; painted only inside it).
// A piece is painted on its prefix, from its first cell up to and including its first lost cell
// (the first cell of its path another piece owns, where the chain continues on the owner; its
// last cell if it loses none), and on its suffix, from the earliest of its cells after its own
// prefix end that is another piece's first lost cell to its last cell (a prefix end on its
// prefix needs no suffix: the water runs on along the prefix, and a suffix from there would
// paint the piece past its first lost cell again). A painted cell carries its owner's values; a path
// cell in no prefix or suffix is not painted (channel 0), so the crossing a piece made past its
// first lost cell (S3.10) is gone. The network is then closed (S3.14): a suffix can run over a
// cell that a third piece owns past that piece's own prefix end, and the chain there steps on
// along the third piece's path onto a cell no prefix or suffix paints (S3.11's dry tail, where
// the water left the pinned network); every painted cell's chain is followed along the
// directions over such cells, each painted, until it meets a painted cell or a piece's last
// cell. So every painted chain runs on painted cells to its last cell (a sea cell, a land
// sink's junction) or off the domain; a painting without dry tails is unchanged, and drytails
// counts what is left (0). A junction is reached by its
// cell's outflow leg (its first cell); an inflow leg that meets another piece before the
// junction merges there, and from there to the junction the painted discharge is the owner's
// alone (recorded, not corrected). Reach: a cell's painting depends on its owner's first lost
// cell and on the first lost cells of the pieces crossing the owner's path, each decided by the
// ownership on that piece's own cells: pieces within two global cells of the cell's block, whose
// owners involve pieces starting up to three cells away; a cell the closure paints depends also
// on the suffix painted upstream on its owner's path, that is on one more piece crossing it
// (measured on windows A, B, C at k 16 and A at k 64 of erosion run 9: every walk stays within
// one global cell, at most 22 fine cells and one change of owner, and the two-window tests hold
// with the same reach). So the links and springs of the global
// cells up to three cells outside the domain are drawn: they compete for ownership and end
// prefixes there and are never painted (about 7 % more time at level 2 of window B at k 16).
//
// Ends. A link ends at its first fine cell inside a global sea cell's block (its sea
// receiver's entry point) without `finesea`; with it, see below. A land sink's
// junction ends the links into it (with pinsinks it is pinned without them too). The rest
// of a path is not drawn.
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
// positive and finite, a negative or non-finite meanderroughness, a
// meanderroughnesspersistence outside 0 .. 1, an edgemargin outside 0 .. 0.5, a negative springspercell, a
// diagonaltransit other than 0 and 1, a non-finite sea level, receivers that form a loop
// downstream of a link within reach of the domain or of a transit candidate (loops
// elsewhere go unseen), and, on a link within reach of the domain, a
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
    std::vector<int32_t> globalcell;     // the upstream global cell's index, -1 off the network and on the
                                         // tributaries (regional_tributaries.hpp: a channel cell with -1 is one)
    int sinkjunctions = 0;               // ChannelParameters::pinsinks: the land-sink junctions it pinned in the kept box
    int diagonallinks = 0;               // S3.10: drawn diagonal links whose upstream cell lies in the kept box
    int cornerlinks = 0;                 // of them, those through the shared corner (no transit block, or
                                         // diagonaltransit 0)
    int mergedpieces = 0;                // S3.11: pieces (links and springs) whose upstream cell lies in the kept
                                         // box and whose prefix ends before their last cell (a confluence)
    int drytails = 0;                    // S3.11: of the same pieces, those that own a painted cell whose next path
                                         // cell is not painted (the chain leaves the pinned network there); 0 since
                                         // S3.14's closure
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
// inside its two global cells' blocks (with a transit block, inside the 2 x 2 blocks
// around its corner, mutual 8-neighbours), so those cells lie at most one global cell
// nearer the domain's edge than the channel cell. The regional pass's start surface and fine
// sea agree between windows from 2 x meanpasses global cells inside the domain's edge
// (each block-mean pass spoils the residuals of the outer two), so its channels agree from
// 2 x meanpasses + 2 cells inside it (2 x meanpasses + 1 until S3.10): a prefix end reads the
// owners along the piece's path, whose beds read the start surface one global cell beyond
// (unchanged by S3.14's closure, whose walks stay within a global cell of a painted cell).
// regionalhalo and refinementhalo cover it (channelreachglobalcells, refinement_grid.hpp).
// Without a start surface and fine sea every window agrees to its edge.
// `relief_m` (multiscale plan S3.14; empty: 0 everywhere) is the local relief of the global cells
// (W x H, m: the regional noise amplitude's), which scales the legs' roughness (ChannelParameters).
// Throws std::invalid_argument also for a non-empty `start_m` or `finesea` whose size
// is not the domain's, a non-empty `relief_m` whose size is not the world grid's, a negative or
// non-finite meanderroughnessreliefgain and a meanderroughnessrelief_m that is not positive and finite.
FineChannels rasterisechannels(const RegionalWindow& window, const GlobalNetwork& network,
                               const ChannelParameters& parameters, float sealevel,
                               const std::vector<float>& start_m = {}, const std::vector<uint8_t>& finesea = {},
                               const std::vector<float>& relief_m = {});

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

struct LevelContext;  // regional_pass.hpp (the hook interface, multiscale plan S3.0)
struct LevelState;

// LevelStep::channels: rasterisechannels on the level's window with the incision's network and
// bed surface when the incision is active (state.incisednetwork, state.bedsurface_m;
// incisionlevel), else the inputs' network and state.surface_m, ending at the fine coast
// (state.sea), with pinsinks set in the refinement loop when parameters.relief.catchmentreach >= 1
// (S3.9), the legs' roughness scaled by the local relief of the inputs' elevation (S3.14: the
// regional noise amplitude's relief, the world's relief amplitude with regional.noise_relief_field 1);
// then the tributaries (tributarieslevel, regional_tributaries.hpp); then every
// land channel cell's bed pinned into state.surface_m. Writes state.channels and state.bed_m
// (a copy of the channels' beds).
void channelslevel(const LevelContext& context, LevelState& state);
}
