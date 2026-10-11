#pragma once

// Inputs and parameters of the regional pass (docs/architecture/erosion-hydrology-plan.md,
// Phase E, step E3, Phase R, step R0, Phase R2, step R6, and Phase E5, step E5.0): shared by
// regional_pass.hpp, which composes the pass, and the parts it calls (regional_surface.hpp,
// regional_coast.hpp, regional_edges.hpp, regional_water.hpp, regional_valleys.hpp,
// regional_glacial.hpp, regional_troughs.hpp, regional_cirques.hpp, regional_ice.hpp).
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

    // Phase E5 (the glacial forms, regional_glacial.hpp): the erosion run's
    // equilibrium-line altitude (cryosphere.equilibrium_line_altitude, m above datum, the
    // run's climate cooled by landscape.glacial.cooling_c), W x H or empty (no glacial form
    // is cut and no ice added by the ELA rule without it; the adapter leaves it empty when
    // every value is zero, the world having had no glacial term); and the height by which the
    // present climate's ELA lies above the run's: cooling_c x 1000 / the lapse rate in C/km
    // (the adapter computes it; 0 when unknown).
    std::vector<float> equilibriumline_m;
    double presentelaoffset_m = 0.0;

    // Multiscale plan Phase S3 (contracts of S3.0): the world's sub-grid statistics from S3.6's
    // diagnoser (geosphere.relief_amplitude, m; geosphere.valley_spacing, km), W x H or empty
    // (absent: S3.4's noise falls back to the local relief, S3.3's density to
    // regional.tributary_spacing_km). The adapter fills them from S3.6 on.
    std::vector<float> reliefamplitude_m;
    std::vector<float> valleyspacing_km;
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
    // Multiscale plan S3.15 (regional_coast.hpp; with landshare > 0 only): the shoreline's position
    // displacement and the sea pockets; 0 for warpcells and pocketreach keeps S3.14's coast bit for bit.
    double warpcells = 0.0;              // the displacement's first-octave amplitude in global cells (0: none)
    double warppersistence = 0.7;        // octave o of the displacement weighs warppersistence^o
    int pocketreach = 0;                 // fine sea in a global land cell is kept only if it connects through fine sea
                                         // to a global sea cell's fine sea within this many global cells (0: kept)
};

// The sea gate and the straight-run breaker (Phase R2, R8, regional_edges.hpp): both off
// at the defaults, which keep the 3 x 3 gate of Phase E and no breaker.
struct EdgeParameters
{
    int smoothgate = 0;                  // 0: a fine cell may be sea only where its global cell's 3 x 3 holds a sea
                                         // cell; n >= 1: R8's smooth gate, a sea cell within n global cells of the
                                         // warped position (regional_edges.hpp, seagate and inlandcells)
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
    // Multiscale plan S3.5 (contract of S3.0, not read yet): the sides' reach on the ground, km,
    // tested against the level's band like the floor; 0 keeps reachcells.
    double reach_km = 0.0;
};

// Tributaries per level (multiscale plan S3.3, regional_tributaries.hpp): grown from the level's
// channels, hashed by the absolute key, Hack-law lengths in [2 level cells, one parent cell),
// junction angles, Strahler order as the growth priority; the density target is spacing_km until
// S3.6's valley-spacing field. enabled 0: no tributaries (the S2 pass exactly); enabled 1 also
// widens the halo rules by tributaryreach (regional_tributaries.hpp).
struct TributaryParameters
{
    int enabled = 0;
    double spacing_km = 0.13;            // the target valley spacing on the ground (Perron et al. 2009: 30 to 320 m)
    double hackexponent = 0.57;          // h in L = C A^h (0.55 to 0.7, Tarboton 1996)
    double hackcoefficient = 1.27;       // C in km per km2^h
    double junctionangle_deg = 72.0;     // the mean angle between a tributary and its receiving channel's upstream
    double anglespread = 0.2;            // the angle is the mean x (1 + anglespread x a hashed value in [-1, 1))
    double nodemargin = 0.25;            // an inherited node lies floor(nodemargin x factor) cells inside its parent cell
    double lengthminfraction = 0.75;     // the target length is the capped Hack length x [lengthminfraction, 1)
};

// Glacial troughs (Phase E5, step E5.1, regional_troughs.hpp): a land channel cell of a
// glaciated cell is the axis of a trough whose cross-section is a power law of the
// distance (exponent 2: a parabola, the U shape), its half-width from the ice catchment of
// the channel's global cell (regional_glacial.hpp, glacialcatchment). width_km 0: off (the
// Phase R2 pass exactly).
struct TroughParameters
{
    double width_km = 0.0;               // half-width on the ground at the reference ice catchment (0: no troughs)
    double reference_km2 = 1000.0;       // the ice catchment at which the half-width is width_km
    double widthexponent = 0.5;          // half-width = width_km x (A_i / reference_km2)^widthexponent
    double widthmax_km = 40.0;           // the half-width's cap, km (and at most k fine cells)
    double wall_m = 1000.0;              // the cross-section's height above the bed at the half-width
    double exponent = 2.0;               // target = bed + wall_m x (d / W)^exponent, d the ground distance
    double glaciationmin = 0.5;          // a channel cell is a trough axis where the fine glaciation is at least this
};

// Cirques (Phase E5, step E5.2, regional_cirques.hpp): hashed candidate points in every
// glaciated global cell's block become parabolic bowls where the fine surface lies above
// the run's ELA by less than band_m (above that is the plateau ice, not a cirque's
// headwall). radius_km 0: off.
struct CirqueParameters
{
    double radius_km = 0.0;              // the bowl's radius on the ground (0: no cirques); at most k fine cells
    double depth_m = 300.0;              // the floor below the centre's surface, unless the run's ELA at the centre is higher
    double band_m = 1500.0;              // candidates more than this above the run's ELA are not cirques
    int candidates = 2;                  // hashed candidate points per glaciated global cell
    double glaciationmin = 0.5;          // the global cell's glaciation at least this
};

// Ice above the ELA (Phase E5, step E5.3, regional_ice.hpp): a fine land cell whose
// surface lies above the present climate's ELA carries ice, added to the present-day ice
// of regional_water.hpp (never removed). thicknessfraction 0: off.
struct IceParameters
{
    double thicknessfraction = 0.0;      // thickness = this x (surface - present ELA), at least 1 m (0: no ELA ice)
    double thicknessmax_m = 1000.0;      // the thickness's cap
    double edgenoise_m = 100.0;          // the threshold's noise: ice where surface >= ELA + edgenoise_m x noise (own salt)
};

struct RegionalParameters
{
    double noiseamplitude_m = 20.0;      // noise amplitude where the global relief is zero
    double noiserelieffraction = 0.25;   // added amplitude per metre of local global relief
    double noisepersistence = 0.6;
    double coastband_m = 50.0;           // the noise tapers to zero within this height of the sea level
    int meanpasses = 2;                  // block-mean correction passes (upsampling.hpp)
    // Multiscale plan S3.4 (contract of S3.0, not read yet): 1 takes the noise amplitude per level
    // from inputs.reliefamplitude_m (S3.6's field; localrelief without it); 0: the local relief.
    int noiserelieffield = 0;
    WarpParameters warp;
    CoastParameters coast;
    EdgeParameters edges;
    ShoreParameters shores;
    ValleyParameters valleys;
    TributaryParameters tributaries;
    TroughParameters troughs;
    CirqueParameters cirques;
    IceParameters ice;
    ChannelParameters channels;
    AnalyticalReliefParameters relief;
};
}
