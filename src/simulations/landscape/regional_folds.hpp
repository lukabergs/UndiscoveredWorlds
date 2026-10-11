#pragma once

// The fold fabric of the regional pass (docs/architecture/multiscale-plan.md, section 3, "S5.2:
// the fold fabric (design, 2026-10-11)" with its amendment, step S5.2, and "Phacelle noise for the
// fold fabric", step S5.3b): fold-train relief across the shortening direction the tectonics
// exports (S5.1), its wavelength from the sediment cover (Biot), its amplitude from the strain.
// Called by regional_pass.cpp after the solve and before the valley floors (LevelStep::foldfabric),
// so the floors, troughs and cirques cut the folded land; around the lakes step it takes the fabric
// out of the surface the water decides on and puts it back on dry land (unfoldforlakes,
// refoldafterlakes). Planet-free: the adapter (src/simulations/regional/regional_inputs.cpp) fills
// RegionalGlobalInputs::foldstrain, folddirection_rad and sediment_m.
//
// Contract (S5.2, S5.3b). Inert, and no field read, while folds.amplitudemax_m <= 0,
// folds.splatspacing_km <= 0 or inputs.foldstrain is empty (erosion run 9 has no shortening: its
// windows are bit-identical to the pass without the unit). Per global cell: the strain eps
// (foldstrain), the wavelength L = foldwavelength_km of the sediment cover (the cover 0 without
// sediment_m) and the across-strike axis theta (folddirection_rad, radians in [0, pi) from east
// toward south; 0 without the field). Per domain cell (land: not sea, not a channel), at its centre
// p in absolute global coordinates:
// - the cell's eps, L and height above the sea level are the global fields (eps, L, max(0, base
//   elevation - sea level)) interpolated by upsample, as the noise amplitude's relief is; the
//   cell's amplitude A = foldamplitude_m (the weight smoothstep(eps; strainmin, 2 strainmin)
//   included); a cell with A = 0 is not changed;
// - the gates (S5.3b): no fabric where inputs.glaciation (cryosphere.glaciation, 0 to 1) upsampled
//   and clamped to 0 .. 1 is above folds.icemax (icemax 1, or no field: no gate), nor in the block of
//   a global lake cell (lake surface above 0, not sea: regional_water.cpp's lake cells);
// - the cell is drawn only when its L lies in the level's octave band [longest / k, longest /
//   parentscale) with longest = 2 pi R / H (twice the global north-south cell: the band of
//   addlevelnoise, [2 c_l, 2 c_(l-1)) in level cells); parentscale 0 (level 1, the single pass):
//   no upper end; finer levels inherit the fabric by the prolongation;
// - the pattern (S5.3b, Phacelle's layout and weight): one splat per cell of an absolute lattice,
//   the global grid subdivided by s = 2^j (j >= 0 with the lattice cell pi R / (H s) nearest
//   splatspacing_km on the ground, in the ratio: j = round(log2(the global cell / the spacing)),
//   0 to 20), the splat at (column + fx, row + fy) / s with (fx, fy) the two halves of
//   deterministiccombine(seed, foldsalt, column, row) over 2^32; each splat's axis d = (cos theta,
//   sin theta) (east, south) with 2 theta the angle of the doubled-angle field
//   w(eps) (cos 2 theta, sin 2 theta) of the global cells interpolated by bicubic at the splat
//   (the unweighted field where the weighted one vanishes there), and its wavelength the global L
//   interpolated by bicubic there, clamped to the min and max; z(p) = sum over the splats within
//   r = 2 lattice cells on the ground (the lattice rows floor(v s) - 2 .. + 2 and the columns
//   within 2 / the row's east-west ratio lattice cells: every splat of positive weight, 5 x 5 at
//   the equator, wider towards the poles) of g_i exp(i 2 pi <d_i, p - c_i> / L_i), the phase
//   origin the splat's own position c_i, g_i = max(0, exp(-|p - c_i|^2 / sigma^2) - exp(-r^2 /
//   sigma^2)) (sigma = kernelwidth_km, or L_i with kernelwidth_km 0), ground distances in km with
//   the column offsets (wrapped) scaled by the cell's row ratio (eastwestratio); crests run at
//   theta + pi / 2. A pure function of the global fields, the seed and absolute position: every
//   level and window evaluates the same splats;
// - the fabric F = A sin(arg z) (phasor 1: constant amplitude, fold trains) or A Re z / sum g
//   (phasor 0), times the water gap 1 - max over the channel cells within gapfraction x L (ground
//   distance, rows at most k, columns columnextent) of 1 - smoothstep(d / (gapfraction x L));
// - the cell's elevation e becomes e + F, never lowered below S = the sea level + coast.shorefloor_m
//   (at least the next float above the sea level): max(e + F, S); a cell with e <= S (dry land at the
//   shore or below the sea level: a filled pocket, a basin) is not changed, since it could only rise
//   and the projection would have to sink its parent's other children (S5.2's first judging run:
//   24 forced parents and a land cell 751 m below the sea level at window A's bay). Channel and sea
//   cells are never changed; nothing else is read. The refinement loop's landbelowsea counts what the
//   projection may still push below the sea level.
// The water gap reads the level's channels within gapfraction x L: under half a parent cell at
// level >= 2 (the band), at most one global cell at level 1 (the cap k), within formreach.

#include <cstdint>
#include <vector>

#include "simulations/landscape/channel_network.hpp"
#include "simulations/landscape/regional_grid.hpp"
#include "simulations/landscape/regional_types.hpp"

namespace landscape
{
// True when the fold fabric acts: folds.amplitudemax_m > 0, folds.splatspacing_km > 0 and
// inputs.foldstrain holds width x height values.
bool foldfabricon(const RegionalGlobalInputs& inputs, const RegionalParameters& parameters);

// L in km: clamp(wavelengthbase_km + wavelengthpercover_km x sediment_m / 1000 (0 when !withsediment), the smaller
// and the larger of wavelengthmin_km and wavelengthmax_km).
double foldwavelength_km(const FoldParameters& folds, double sediment_m, bool withsediment);

// A in m: min(relieffraction x (L / pi) sqrt(eps) (L in m), amplitudemax_m, reliefelevationfraction x
// max(0, heightabovesea_m)) x smoothstep(eps; strainmin, 2 strainmin); 0 when amplitudemax_m, eps or
// L is not positive.
double foldamplitude_m(const FoldParameters& folds, double strain, double wavelength_km, double heightabovesea_m);

// The contract above on the window's domain (domainfinewidth() x domainfineheight(); channels.channel,
// sea and elevation_m of that size, otherwise nothing happens). parentscale: the parent level's
// scale k_(l-1) (the band's upper end), 0 for none. floorclamped (optional): set to one value per domain
// cell, 1 where the fabric was raised to the shore floor.
void addfoldfabric(const RegionalGlobalInputs& inputs, const RegionalWindow& window, const RegionalParameters& parameters,
                   const FineChannels& channels, const std::vector<uint8_t>& sea, std::vector<float>& elevation_m,
                   int parentscale = 0, std::vector<uint8_t>* floorclamped = nullptr);

struct LevelContext;  // regional_pass.hpp (the hook interface, multiscale plan S3.0)
struct LevelState;

// LevelStep::foldfabric: with foldfabricon, addfoldfabric on state.elevation_m with the level's band
// (parentscale k_(l-1) at level >= 2, none at level 1 and in the single pass), then channelcellmeans
// with channelcellcap_m (regional_valleys.hpp). The loop's diagnostics: LevelDiagnostics::foldcells,
// foldamplitude_p50_m and foldamplitude_p99_m (the kept cells the fabric changed, before the channel
// cells' means), and foldfloorclamped (the kept cells the shore floor raised). S5.3b: keeps the
// elevation before the fabric in state.unfolded_m and after it in state.folded_m for the lakes step.
// Off: nothing, and the step is not recorded.
void foldfabriclevel(const LevelContext& context, LevelState& state);

// S5.3b, the lakes step on the surface without the fabric (regional_pass.cpp calls it just before
// lakeslevel; nothing without state.folded_m). Every cell that is neither sea nor a channel and whose
// value the fabric changed takes its un-folded value u: the elevation before the fabric where the forms
// after it (valley floors, troughs, cirques) left the folded value, else the smaller of the present
// value and the one before the fabric (a form cut below the folded surface); state.folded_m keeps the
// present (folded) value and state.unfolded_m u. Other cells are left.
void unfoldforlakes(const LevelContext& context, LevelState& state);

// S5.3b, after lakeslevel: the fabric back on the dry land the lakes left as it was. A cell that
// unfoldforlakes changed and that is dry with the elevation u takes its folded value back, never at
// or below the highest surface of its wet 8-neighbours (the next float above it), and keeps u where
// u itself lies at or below that surface; a wet cell, or one the lakes changed (a shore raise, a
// rim), keeps the lakes' value: no fabric under the water or on its shore. Then state.folded_m and
// state.unfolded_m are emptied.
void refoldafterlakes(const LevelContext& context, LevelState& state);
}
