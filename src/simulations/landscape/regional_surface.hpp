#pragma once

// The start surface and the fine sea of the regional pass (docs/architecture/erosion-hydrology-plan.md,
// Phase E, step E3, and Phase R: R0 moved them here from regional_pass.cpp
// unchanged; R2 owns this file: the warp, the coast noise and the land-share rule).
// Every value is a pure function of the global fields, the seed and the fine cell's
// absolute position within the halo's reach, so two windows agree on their overlap.
// Planet-free.

#include <cstdint>
#include <vector>

#include "simulations/landscape/regional_grid.hpp"
#include "simulations/landscape/regional_types.hpp"
#include "simulations/landscape/upsampling.hpp"

namespace landscape
{
// The domain warp of the window (WarpParameters): empty when amplitude_cells is 0.
// R2: a smooth vector noise on the sphere (spherical_noise.hpp, its own seed salt),
// at most amplitude_cells global cells long, a pure function of the fine cell's
// absolute position and the seed.
FineWarp regionalwarp(const RegionalWindow& window, std::uint64_t seed, const RegionalParameters& parameters);

struct StartSurface
{
    std::vector<float> base_m;    // the global base elevation upsampled at the (faded, R2) warped positions
    // The global base elevation upsampled at the unwarped fine centres: the reference of
    // the lapse-rate correction, because the temperatures (and the runoff) are upsampled
    // without the warp (plan section 6, 'Decided for Phase R'), so the correction takes
    // the fine elevation's whole difference from the elevation the climate was
    // interpolated with, the warp's displacement included.
    std::vector<float> lapsebase_m;
    std::vector<float> start_m;   // base_m plus the regional noise, block means corrected
    double blockmeanerror_m = 0.0;  // the largest |block mean - global value| over the kept window's
                                    // global cells (RegionalFields::blockmeanerror_m)
};

// The start surface on the window's domain (domainfinewidth() x domainfineheight()):
// the base elevation upsampled, noise of the band from twice the global north-south
// spacing (exclusive) down to twice the fine one with an amplitude from the local
// global relief, tapered within coastband_m plus the amplitude of the sea level
// (Phase E; R2 replaces the taper by CoastParameters), block means corrected
// meanpasses times.
StartSurface regionalstartsurface(const RegionalGlobalInputs& inputs, const RegionalWindow& window,
                                  const FineWarp& warp, const RegionalParameters& parameters);

// The fine sea of the domain, by local rules only (no flood fill): a fine cell is
// sea when its surface is at or below the sea level and its global cell or one of
// the eight neighbours is sea, and a fine cell of a global lake cell is never sea
// unless its global cell is. Phase E's faithfulness (the centre rule, kept while
// coast.landshare is 0): each global cell's centre fine cell follows the global
// mask, so `surface_m` is changed there (a sea centre lowered to the sea level, a
// land centre counted as sea raised to the next float above it). R2 replaces the
// centre rule by the land share (CoastParameters::landshare) and may change
// `surface_m` near the coast (the sea side pushed down). R1 calls it before the
// channels are drawn, so the channels end at the fine coast.
std::vector<uint8_t> regionalsea(const RegionalGlobalInputs& inputs, const RegionalWindow& window,
                                 const FineWarp& warp, const RegionalParameters& parameters,
                                 std::vector<float>& surface_m);
}
