#pragma once

// The start surface and the fine sea of the regional pass (docs/architecture/erosion-hydrology-plan.md,
// Phase E, step E3, and Phase R: R0 moved them here from regional_pass.cpp
// unchanged; R2 added the warp, the coast noise and the land-share rule).
// Every value is a pure function of the global fields, the seed and the fine cell's
// absolute position within the halo's reach, so two windows agree on their overlap.
// Planet-free.

#include <cstdint>
#include <optional>
#include <vector>

#include "simulations/landscape/regional_grid.hpp"
#include "simulations/landscape/regional_types.hpp"
#include "simulations/landscape/spherical_noise.hpp"
#include "simulations/landscape/upsampling.hpp"

namespace landscape
{
// A point-evaluated fractal noise on the planet-wide fine grid of a scale k (R2's coast
// noise; R4 uses it for the lake shores and ice edges with a salt of its own). Gradient
// noise on the unit sphere (spherical_noise.hpp) at the fine cell's centre, a pure
// function of the fine cell's absolute position, the seed and the salt: octave o has the
// wavelength of one global north-south spacing (pi / H on the unit sphere) divided by
// 2^o, down to the last octave of at least two fine spacings (floor(log2 k) octaves; none
// for k < 2, where the value is 0), weight persistence^o, octave seeds from
// deterministiccombine(seed, salt). Normalised by FractalNoiseField::normalisation over
// those octaves: the standard deviation over the sphere is about 1 (octavedeviation,
// 0.268 per octave; 1.016 measured, Landscape.RegionalSurfaceNoise); |value| is bounded by
// (sum of the weights) / normalisation if an octave stays within +-1 (6.1 at k = 16 with
// persistence 0.5); the largest of 35,055 cells of a 64 x 32 world at k = 16 was 3.46.
class RegionalNoise
{
public:
    RegionalNoise(int worldwidth, int worldheight, int scale, std::uint64_t seed, std::uint64_t salt, double persistence);
    // At the fine cell (I, J) of the kW x kH grid (I taken modulo kW, 0 <= J < kH).
    double at(std::int64_t finecolumn, std::int64_t finerow) const;
    bool active() const { return field_.has_value(); }

private:
    int finewidth_ = 0;
    int fineheight_ = 0;
    std::optional<FractalNoiseField> field_;
};

// The domain warp of the window (WarpParameters): empty when amplitude_cells is 0.
// Two fractal gradient noises on the unit sphere (salts of their own) with the first
// wavelength wavelength_cells global north-south spacings, octaves halved down to one
// spacing, persistence 0.5, each normalised to a standard deviation of about 1; the
// vector w of the two is shortened to amplitude_cells x tanh(|w| / 2) w / |w|, so the
// displacement is smooth and shorter than amplitude_cells (about half of it typically).
// A pure function of the fine cell's absolute position and the seed.
FineWarp regionalwarp(const RegionalWindow& window, std::uint64_t seed, const RegionalParameters& parameters);

struct StartSurface
{
    std::vector<float> base_m;    // the global base elevation upsampled at the faded warped positions
    // The global base elevation upsampled at the unwarped fine centres: the reference of
    // the lapse-rate correction, because the temperatures (and the runoff) are upsampled
    // without the warp (plan section 6, 'Decided for Phase R'), so the correction takes
    // the fine elevation's whole difference from the elevation the climate was
    // interpolated with, the warp's displacement included.
    std::vector<float> lapsebase_m;
    std::vector<float> start_m;   // base_m plus the regional noise, block means corrected
    double blockmeanerror_m = 0.0;  // the largest |block mean - global value| over the kept window's
                                    // global cells (RegionalFields::blockmeanerror_m)
    // The land-share rule's diagnostics: the offset passes after which no block of the
    // domain missed its share (-1: some still missed at the cap), and the domain's blocks
    // that miss it in the final mask (0 unless the cap was reached).
    int sharepasses = -1;
    int sharemisses = 0;
};

// The start surface on the window's domain (domainfinewidth() x domainfineheight()).
// `warp` must be regionalwarp's (the land-share rule evaluates the same warp beyond the
// domain). Per fine cell:
// - the elevation's warp is the full warp times max(0, 1 - |unwarped base - sea level| /
//   warp.fade_m) (fade_m 0: no fade), so the coast takes nearly all of it and inland
//   valleys stay under the pinned channels; base_m is the base elevation at that position;
// - Phase E's noise: the band from twice the global north-south spacing (exclusive) down
//   to twice the fine one, amplitude noiseamplitude_m + noiserelieffraction x the local
//   global relief, tapered within coastband_m plus the amplitude of the sea level;
// - the coast noise (coast.noiseamplitude_m > 0): RegionalNoise (persistence
//   coast.noisepersistence) times a = coast.noiseamplitude_m + coast.noiserelieffraction x
//   the local global relief, times the coastal weight (the global cells whose 3 x 3
//   neighbourhood holds sea and land, interpolated at the faded warped position and
//   clamped to 0 .. 1). With coast.landshare > 0 only the decision surface (below) takes
//   it, so the coastline is fractal while coastal land stays as rough as inland land;
//   without the land-share rule the elevation takes it no larger than Phase E's noise at
//   that height (its amplitude with its taper at the sea level);
// then the block means corrected meanpasses times (Phase E's correction).
// With coast.landshare > 0 the fine sea is decided here, on a decision surface: the
// corrected surface with the full coast noise, its sea side pushed down (a height h below
// the sea level, before the coast noise, sinks by seadrop_m x the coastal weight x
// min(1, -h / a)) and lifted by seadrop_m x the interpolated share of global cells that may
// hold no sea (so sea does not run into them along a block edge), plus a continuous offset
// (bilinear between global cell centres, weights from the fine cell's position inside its
// global cell) that brings every global cell's block to at least `landshare` land (a land
// cell) or at most 1 - landshare land (a sea cell): Jacobi passes until no block misses, at
// most coast.sharepasses (StartSurface::sharepasses and sharemisses report it); sea only
// where Phase E allowed it (regionalsea); diagonal pinches (the sea cells of a 2 x 2 square
// whose land touches only diagonally) made land in one sweep, then one-cell holes (a sea
// cell whose four neighbours are land) filled, the shares counted after both. The elevation
// then follows the mask: a cell the mask moves across the sea level is at least as deep (a
// sea cell) or as high (a land cell) as the decision there, capped at coast.shorefloor_m and
// at the distance of the cell's own global cell from the sea level (so every block can
// still reach its mean); the other cells take that floor only in a band along the shoreline
// (full while the decision is within the floor of the sea level, none beyond twice it), so
// shelves keep their shape; the rest is clamped (a sea cell above the sea level to it, a
// land cell at or below it to the next float above it), and the block means are restored by
// coast.signedpasses passes that never cross the sea level and never bring a cell nearer it
// than it then is or than its floor, so no zero-depth sheet forms. All of it is computed
// on the domain widened by r x sharepasses + 2 + 2 x signedpasses global cells (r, a share
// pass's reach: 1 for k >= 3, 2 for k <= 2), so the domain's values depend on the window
// only within 2 x meanpasses global cells of its edge, as in Phase E: the halo rule
// halo >= 2 x meanpasses + 2 holds.
StartSurface regionalstartsurface(const RegionalGlobalInputs& inputs, const RegionalWindow& window,
                                  const FineWarp& warp, const RegionalParameters& parameters);

// The fine sea of the domain, by local rules only (no flood fill): a fine cell is
// sea when its surface is at or below the sea level and its global cell or one of
// the eight neighbours is sea, and a fine cell of a global lake cell is never sea
// unless its global cell is. With coast.landshare 0 (Phase E's faithfulness, the centre
// rule): each global cell's centre fine cell follows the global mask, so `surface_m` is
// changed there (a sea centre lowered to the sea level, a land centre counted as sea
// raised to the next float above it). With coast.landshare > 0 the start surface already
// carries the land-share mask in its sign (regionalstartsurface) and `surface_m` is not
// changed. R1 calls it before the channels are drawn, so the channels end at the fine coast.
std::vector<uint8_t> regionalsea(const RegionalGlobalInputs& inputs, const RegionalWindow& window,
                                 const FineWarp& warp, const RegionalParameters& parameters,
                                 std::vector<float>& surface_m);
}
