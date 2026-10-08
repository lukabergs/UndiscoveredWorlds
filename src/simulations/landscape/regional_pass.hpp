#pragma once

// The regional pass (docs/architecture/erosion-hydrology-plan.md, Phase E, step E3,
// and section 6, "Decided for Phase E"): global fields and a window in, the window's
// fine fields out. Steps: the start surface (the global filled elevation upsampled,
// band-limited noise added with an amplitude from the local relief and tapered near
// the sea level, block means corrected); the pinned channels; the analytical relief
// on the domain; sea, lakes, rivers and climate on the fine grid by local rules;
// the halo discarded. Every step is a pure function of the global fields, the seed
// and absolute position within the halo's reach, so two windows give the same
// values on their overlap (test Landscape.RegionalOverlap). Planet-free: the planet
// adapter is src/simulations/regional/regional_inputs.hpp.

#include <cstdint>
#include <vector>

#include "simulations/landscape/analytical_relief.hpp"
#include "simulations/landscape/channel_network.hpp"
#include "simulations/landscape/regional_grid.hpp"

namespace landscape
{
// The global fields the pass reads, W x H row-major, north row first, float getters
// only (no whole-metre values).
struct RegionalGlobalInputs
{
    int width = 0;
    int height = 0;
    double radius_m = 0.0;               // planetbody.radius
    float sealevel = 0.0f;               // geosphere.sea_level, m above datum
    std::uint64_t seed = 0;              // the world seed
    float lapserate_c_per_km = 0.0f;     // atmosphere.lapse_rate
    std::vector<float> elevation_m;      // geosphere.base_elevation (the land surface of procedural worlds)
    GlobalNetwork network;               // the routing and river fields (sea mask included)
    std::vector<float> runoff_mm_yr;     // hydrosphere.runoff
    std::vector<float> lakesurface_m;    // hydrosphere.lake_surface, 0 where no lake
    std::vector<float> januarytemperature_c, julytemperature_c;
    std::vector<float> januaryprecipitation_mm, julyprecipitation_mm;  // mm/month
    std::vector<float> icethickness_m;   // cryosphere.ice_thickness (present-day ice)
    std::vector<int16_t> climate;        // the global climate class
};

struct RegionalParameters
{
    double noiseamplitude_m = 20.0;      // noise amplitude where the global relief is zero
    double noiserelieffraction = 0.25;   // added amplitude per metre of local global relief
    double noisepersistence = 0.6;
    double coastband_m = 50.0;           // the noise tapers to zero within this height of the sea level
    int meanpasses = 2;                  // block-mean correction passes (upsampling.hpp)
    ChannelParameters channels;
    AnalyticalReliefParameters relief;
};

// The fine fields of the kept window (finewidth() x fineheight(), row-major, north
// row first); the bundle writer stores them with the window.
struct RegionalFields
{
    RegionalWindow window;
    std::vector<float> elevation_m;
    std::vector<uint8_t> sea;
    std::vector<uint8_t> riverdirection;     // 1..8 on channel cells, 0 elsewhere
    std::vector<float> discharge_m3s;        // pinned channels: the global discharge; elsewhere the fine accumulation
    std::vector<int32_t> januaryflow, julyflow;  // legacy units, on pinned channels
    std::vector<float> lakesurface_m;        // 0 where no lake water
    std::vector<float> januarytemperature_c, julytemperature_c;
    std::vector<float> januaryprecipitation_mm, julyprecipitation_mm;
    std::vector<float> icethickness_m;
    std::vector<int16_t> climate;
};

RegionalFields runregionalpass(const RegionalGlobalInputs& inputs, const RegionalWindow& window,
                               const RegionalParameters& parameters);
}
