#pragma once

// Inputs and parameters of the regional pass (docs/architecture/erosion-hydrology-plan.md,
// Phase E, step E3, and Phase R, step R0): shared by regional_pass.hpp, which composes
// the pass, and the parts it calls (regional_surface.hpp, regional_water.hpp).
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
struct CoastParameters
{
    double noiseamplitude_m = 0.0;       // amplitude at the global-cell scale (octaves down to two fine cells)
    double noisepersistence = 0.5;
    double seadrop_m = 0.0;              // how far the sea side is pushed below the sea level
    double landshare = 0.0;              // a global land cell's block is at least this share land, a sea cell's
                                         // at most 1 - landshare (0: the centre rule of Phase E)
};

// Lake shores and ice edges (Phase R, R4): the zero contour of a signed per-cell
// indicator interpolated at the warped positions plus fractal noise. shorenoise 0
// keeps the Phase E rules (lake water on every fine bed below the global surface,
// ice masked by the global cells).
struct ShoreParameters
{
    double shorenoise = 0.0;             // noise amplitude as a fraction of the indicator's scale (0: Phase E rules)
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
    ShoreParameters shores;
    ChannelParameters channels;
    AnalyticalReliefParameters relief;
};
}
