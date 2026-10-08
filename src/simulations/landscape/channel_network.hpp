#pragma once

// Pinned channels of the regional pass (docs/architecture/erosion-hydrology-plan.md,
// Phase E, step E2b, and section 6, "Seams of the solve"): the global drainage
// network drawn onto a window's fine grid. Every global land cell whose drainage
// area reaches pinarea_km2 contributes one link, a path of fine cells from its
// centre fine cell (k i + k / 2, k j + k / 2) to its receiver's centre fine cell, or
// to its first fine cell inside any global sea cell's block. The path is a pure
// function of the two global cells and the seed (a meander displaced sideways by
// noise that is zero at both ends), so a link that two windows share is drawn
// identically in both; links meet at the centre fine cells, except where they end in
// the sea. Planet-free.

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
    double meanderamplitude = 0.3;     // largest sideways displacement as a fraction of the link length,
                                       // capped at 0.35 (k - 2) fine cells (0 for k <= 2): two links'
                                       // segments not adjacent on the global tree lie at least
                                       // k / sqrt(2) apart (cardinal or diagonal) unless they cross,
                                       // and the cap keeps their curves sqrt(2) apart, so they share
                                       // no fine cell; no path reaches another cell's centre fine cell
    std::uint64_t seed = 0;            // derived from the world seed with a fixed salt
    // Phase R, R3 (the legacy scheme, plan section 6): a link enters and leaves each
    // global cell at points on the shared edges (a hash of the edge; a diagonal through
    // the shared corner), bends through a junction near the cell's centre (a hash of
    // the cell) and follows Catmull-Rom legs through jittered midpoints, inside its
    // own cell. Read by R3; Phase E's geometry ignores them.
    double meanderdischargeexponent = 0.0;  // midpoint displacement x min(1, (Q / meanderreference_m3s)^exponent)
    double meanderreference_m3s = 1000.0;
    double junctionjitter = 0.25;      // the junction lies within this fraction of k of the cell's centre
    double edgemargin = 0.1875;        // edge points stay this fraction of k away from the cell's corners
    int springspercell = 0;            // springs from hashed points into the cell's own river (source cells)
};

// The pinned channels on the window's domain grid (domainfinewidth() x
// domainfineheight(), row-major). A global cell's centre fine cell belongs to its
// own outgoing link; any other fine cell on two links belongs to the link with the
// lower bed there (equal beds: the link whose upstream cell is fewer links from its
// outlet on the global tree, then the smaller absolute index of the upstream global
// cell), so that the bed never rises along the directions (ownership by drainage
// area cannot keep that) and the directions hold no loop whatever the drainage
// areas. Only links adjacent on the global tree (the same receiver, or one's
// receiver the other's upstream cell) share fine cells, so a chain that runs onto
// another link reaches the same outlet; the one accepted exception is two
// diagonals crossing in an X (possible on flats routed by flood order): where their
// staircases share a cell, the lower bed takes it and the higher link's chain
// continues into the lower link's basin. A link ends at its first fine cell inside
// any global sea cell's k x k block (its sea receiver's, or the corner of a sea side
// neighbour that a diagonal or a meander crosses), whose bed is the sea level
// exactly; the rest of its path is not drawn. Along a link the bed descends linearly
// with the path length from the upstream cell's filled elevation to the receiver's,
// or to the sea level for a link that ends in the sea. Every link whose path
// reaches the domain is drawn there, up to the domain's edge (a path stays within
// half a global cell of its two global cells), so a fine cell's values are the same
// in every window that holds it. Throws std::invalid_argument for a scale below 1, a network not on the
// window's world grid, a domain off it, a negative or non-finite meanderamplitude,
// a non-finite sea level, receivers that form a loop downstream of a link within
// reach of the domain (loops elsewhere go unseen), and, on a link within reach of
// the domain, a filled elevation that is not finite or a receiver (or, for a link
// that ends in the sea, the sea level) above the upstream cell's filled elevation.
struct FineChannels
{
    std::vector<uint8_t> channel;        // 1 on a pinned channel cell
    std::vector<uint8_t> direction;      // 1..8 towards the owner link's next cell (0 at a link's last cell:
                                         // its first fine cell in a sea block, or a land sink's centre)
    std::vector<float> bed_m;            // pinned bed elevation, m above datum
    std::vector<float> discharge_m3s;    // the upstream global cell's discharge
    std::vector<int32_t> januaryflow;    // the upstream global cell's legacy flows
    std::vector<int32_t> julyflow;
    std::vector<int32_t> globalcell;     // the upstream global cell's index, -1 off the network
};

// Phase R, R1 (contract of R0; empty vectors keep the Phase E behaviour above):
// `start_m` is the regional start surface and `finesea` the fine sea of the domain
// (domainfinewidth() x domainfineheight(), regional_surface.hpp), decided before the
// channels. With them, a link ends at its first fine sea cell (the fine coast, not
// the global sea block) and no channel cell is pinned through the fine sea; the bed
// follows the start surface's running minimum downstream, held between the
// receiver's end value and the upstream cell's filled elevation (no straight-profile
// trenches), never rising along the directions, and at least the next float above
// the sea level on a land cell.
FineChannels rasterisechannels(const RegionalWindow& window, const GlobalNetwork& network,
                               const ChannelParameters& parameters, float sealevel,
                               const std::vector<float>& start_m = {}, const std::vector<uint8_t>& finesea = {});
}
