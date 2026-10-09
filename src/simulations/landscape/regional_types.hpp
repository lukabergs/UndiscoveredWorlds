#pragma once

// Inputs and parameters of the regional pass (docs/architecture/erosion-hydrology-plan.md,
// Phase E, step E3, Phase R, step R0, and Phase R2, step R6): shared by regional_pass.hpp,
// which composes the pass, and the parts it calls (regional_surface.hpp, regional_coast.hpp,
// regional_edges.hpp, regional_water.hpp, regional_valleys.hpp).
// Planet-free: the planet adapter is src/simulations/regional/regional_inputs.hpp.

#include <cstdint>
#include <vector>

#include "simulations/landscape/analytical_relief.hpp"
#include "simulations/landscape/channel_network.hpp"

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

    // Phase R2 (R6): the world's erosion and tectonic fields read by R7 (the coast's
    // ruggedness and the island sink, regional_coast.hpp) and R9 (the incision,
    // regional_valleys.hpp). Each is W x H or empty; an empty field means 0 everywhere
    // (the tests' synthetic worlds), so the pass runs without them.
    std::vector<float> upliftrate_m_myr;   // geosphere.rock_uplift_rate, signed (negative is subsidence)
    std::vector<float> glaciation;         // cryosphere.glaciation, 0..1
    std::vector<float> denudation_m;       // geosphere.denudation
    std::vector<uint8_t> island;           // geosphere.island: a one-cell island of the world
    std::vector<uint8_t> mountainisland;   // geosphere.mountain_island
    double erosionduration_myr = 0.0;      // geosphere.erosion_duration (0: unknown; the denudation rate is then 0)
};

// The domain warp (Phase R, R2): a smooth vector noise that displaces the sampling
// positions of the global fields, so that the global grid's rows, columns and
// diagonals do not print through coasts, lake shores and ice edges. amplitude 0: no
// warp (FineWarp stays empty, the Phase E results exactly).
struct WarpParameters
{
    double amplitude_cells = 0.0;        // largest displacement, global cells
    double wavelength_cells = 2.0;       // the warp noise's first wavelength, global cells
    double fade_m = 0.0;                 // the elevation's warp fades to 0 this far from the sea level (0: no fade)
};

// The coast (Phase R, R2): a scale-proportional fractal noise kept across the coast,
// the sea side pushed down so the coast is steep, and faithfulness by each global
// cell's land share instead of its centre fine cell. noiseamplitude_m 0, seadrop_m 0
// and landshare 0 keep the Phase E coast (the taper of coastband_m and the centre rule).
// Phase R2, R7 (regional_coast.hpp): the one-sided noise, the amplitude and the push
// from the world's glaciation and denudation, and the island sink; all off at 0.
struct CoastParameters
{
    double noiseamplitude_m = 0.0;       // amplitude at the global-cell scale (octaves down to two fine cells),
                                         // plus noiserelieffraction x the local global relief
    double noisepersistence = 0.5;
    double noiserelieffraction = 0.0;    // added amplitude per metre of local global relief
    double seadrop_m = 0.0;              // the decision surface's push below the sea level on the sea side and lift
                                         // towards global cells that may hold no sea (with landshare > 0 only)
    double shorefloor_m = 20.0;          // the elevation's shore: a fine cell at least this deep or high (as far as
                                         // the decision is), capped at its global cell's distance from the sea level
    int sharepasses = 8;                 // the cap of the land-share offset's Jacobi passes
    int signedpasses = 12;               // the sign-keeping block-mean passes after the land-share mask
    double landshare = 0.0;              // a global land cell's block is at least this share land, a sea cell's
                                         // at most 1 - landshare (0: the centre rule of Phase E)
    // R7. onesided 1: inside a global sea cell the coast noise may lower the decision
    // surface but never lift it above the sea level (the legacy rcoastsquare/rcoastdiamond
    // take an offset at a sea point only if the result stays at or below the sea level), so
    // noise never makes sea into land there; 0: the two-sided noise of R2.
    int onesided = 0;
    double glacialamplitude_m = 0.0;     // added coast-noise amplitude x the global cell's glaciation (0..1)
    double glacialdropfraction = 0.0;    // the sea-side push is multiplied by 1 - this x glaciation (the legacy
                                         // sinks glacial tiles' sea points to 1/12 of the depth: skerries and fjords)
    double denudationamplitude_m = 0.0;  // added amplitude x min(1, denudation / denudationreference_m)
    double denudationreference_m = 1000.0;
    double islandmargin = 0.0;           // the island sink's ring as a fraction of k (the legacy removesmallislands:
                                         // 8 of 16 fine cells); 0: no sink (regional_coast.hpp)
};

// The sea gate and the straight-run breaker (Phase R2, R8, regional_edges.hpp): both off
// at the defaults, which keep the 3 x 3 gate of Phase E and no breaker.
struct EdgeParameters
{
    int smoothgate = 0;                  // 0: a fine cell may be sea only where its global cell's 3 x 3 holds a sea
                                         // cell; 1: R8's smooth gate (the interpolated global sea indicator)
    int straightmax = 0;                 // straight shore runs longer than this many fine cells are broken (0: never)
    int chunkmin = 2;                    // the chunks' sizes in fine cells (the legacy removestraights: 2 to 7)
    int chunkmax = 7;
};

// Lake shores and ice edges (Phase R, R4): the zero contour of a signed per-cell
// indicator interpolated at the warped positions plus fractal noise (regional_water.hpp;
// the noise's persistence is coast.noisepersistence). shorenoise 0 keeps the Phase E
// rules (lake water on every fine bed below the global surface, ice masked by the global
// cells).
struct ShoreParameters
{
    double shorenoise = 0.0;             // noise amplitude as a fraction of the indicator's scale (0: Phase E rules)
    double terrain_m = 50.0;             // the terrain term's scale D: clamp((z - surface) / D, -w, 1) joins the sum
                                         // (w the lake cells' mask interpolated at the warped position, clamped to
                                         // 0 .. 1), water only where z < surface + D; beds lowered by less than 2 D,
                                         // shores within D below the water raised by at most 2 D (regional_water.hpp;
                                         // 0: no terrain term, no bound)
};

// The incision of the pinned channels and the valley floors (Phase R2, R9,
// regional_valleys.hpp and the uplift term of analytical_relief.hpp): all off at the
// defaults (no incision, no uplift term, no valleys: the Phase R pass exactly).
struct ValleyParameters
{
    double upliftscale = 0.0;            // the incision rate takes upliftscale x max(0, the world's rock uplift rate)
    double denudationratescale = 0.0;    // ... or denudationratescale x denudation / erosion duration, the larger
    double floor_km = 0.0;               // a valley floor's half-width is floor_km x sqrt(Q / reference_m3s) km on
                                         // the ground, Q the channel's discharge (0: no valleys); in fine cells
                                         // by the fine north-south spacing, so the width does not change with k
    double reference_m3s = 1000.0;
    double floormax_km = 40.0;           // the half-width's cap, km (and at most k fine cells)
    double sideslope = 0.0;              // rise over run of the valley sides beyond the floor (0: no sides)
    int reachcells = 0;                  // the farthest fine cell a channel shapes; at most k (the halo rule)
};

struct RegionalParameters
{
    double noiseamplitude_m = 20.0;      // noise amplitude where the global relief is zero
    double noiserelieffraction = 0.25;   // added amplitude per metre of local global relief
    double noisepersistence = 0.6;
    double coastband_m = 50.0;           // the noise tapers to zero within this height of the sea level
    int meanpasses = 2;                  // block-mean correction passes (upsampling.hpp)
    WarpParameters warp;
    CoastParameters coast;
    EdgeParameters edges;
    ShoreParameters shores;
    ValleyParameters valleys;
    ChannelParameters channels;
    AnalyticalReliefParameters relief;
};
}
