#pragma once

// Pinned channels of the regional pass (docs/architecture/erosion-hydrology-plan.md,
// Phase E, step E2b, and section 6, "Seams of the solve"): the global drainage
// network drawn onto a window's fine grid. Every global land cell whose drainage
// area reaches pinarea_km2 contributes one link, a path of fine cells from its
// centre fine cell to its receiver's centre fine cell. The path is a pure function
// of the two global cells and the seed (a meander displaced sideways by noise that
// is zero at both ends), so a link that two windows share is drawn identically in
// both; links meet at the centre fine cells. Planet-free.

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
    double meanderamplitude = 0.3;     // largest sideways displacement as a fraction of the link length
    std::uint64_t seed = 0;            // derived from the world seed with a fixed salt
};

// The pinned channels on the window's domain grid (domainfinewidth() x
// domainfineheight(), row-major). A fine cell on two links belongs to the link of
// the larger drainage area (equal areas: the smaller absolute index of the upstream
// global cell). Along a link the bed descends linearly with the path length from
// the upstream cell's filled elevation to the receiver's (never rising: a receiver
// is never higher on the filled surface); a link into a sea cell ends at the sea
// cell's centre at the sea level. Links whose two global cells both lie outside the
// domain are not drawn; a link that leaves the domain is drawn up to its edge.
struct FineChannels
{
    std::vector<uint8_t> channel;        // 1 on a pinned channel cell
    std::vector<uint8_t> direction;      // 1..8 towards the next cell of the link (0 at a link's last cell in the sea)
    std::vector<float> bed_m;            // pinned bed elevation, m above datum
    std::vector<float> discharge_m3s;    // the upstream global cell's discharge
    std::vector<int32_t> januaryflow;    // the upstream global cell's legacy flows
    std::vector<int32_t> julyflow;
    std::vector<int32_t> globalcell;     // the upstream global cell's index, -1 off the network
};

FineChannels rasterisechannels(const RegionalWindow& window, const GlobalNetwork& network,
                               const ChannelParameters& parameters, float sealevel);
}
