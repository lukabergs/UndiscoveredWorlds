#pragma once

// The fold fabric of the regional pass (docs/architecture/multiscale-plan.md, section 3, "S5.2: the fold
// fabric (design, 2026-10-11)" with its amendment, step S5.2; "Phacelle noise for the fold fabric", step
// S5.3b; "S5.3c: the coherent fold phase", step S5.3c): fold-train relief across the shortening direction
// the tectonics exports (S5.1), its wavelength from the sediment cover (Biot), its amplitude from the
// strain, its crests the level sets of one phase field. Called by regional_pass.cpp after the solve and
// before the valley floors (LevelStep::foldfabric), so the floors, troughs and cirques cut the folded land;
// around the lakes step it takes the fabric out of the surface the water decides on and puts it back on dry
// land (unfoldforlakes, refoldafterlakes). Planet-free: the adapter
// (src/simulations/regional/regional_inputs.cpp) fills RegionalGlobalInputs::foldstrain, folddirection_rad
// and sediment_m.
//
// Contract (S5.2, S5.3b, S5.3c). Inert, and no field read, while folds.amplitudemax_m <= 0 or
// inputs.foldstrain is empty (erosion run 9 has no shortening: its windows are bit-identical to the pass
// without the unit). Per global cell: the strain eps (foldstrain), the wavelength L = foldwavelength_km of
// the sediment cover (the cover 0 without sediment_m) and the across-strike axis theta (folddirection_rad,
// radians in [0, pi) from east toward south; 0 without the field).
//
// The phase field (S5.3c, foldphasefield; a pure function of the global fields and the seed, computed per
// call, never stored):
// - the mask: the global cells with smoothstep(eps; strainmin, 2 strainmin) > 0;
// - the axis a = (cos theta, sin theta) (east, south) smoothed (round two): per mask cell, the doubled-angle
//   vectors w(eps) (cos 2 theta, sin 2 theta) of the mask cells within ceil(folds.axissmoothing_cells) rows and
//   columns summed with the weight exp(-d^2 / (2 sigma^2)), sigma = axissmoothing_cells and d the ground distance in
//   global north-south cells (a column offset times the mean of the two rows' east-west cells over the north-south
//   cell), the angle halved; a vanishing sum, or sigma <= 0, keeps the cell's own axis;
// - the orientation, by continuity alone (round two: no vergence; it returns as a profile asymmetry in S5.5):
//   each component's first cell (row-major) takes the sign s that points its axis north (east for an
//   east-west axis), and the other mask cells take theirs from the nearest oriented cell (a Dijkstra flood
//   over the mask's 8-neighbours by ground distance, the sign that keeps the axis within 90 degrees of the
//   predecessor's oriented axis, equal distances to the predecessor of the smaller
//   deterministiccombine(seed, salt, cell));
// - the wavenumber k = s (2 pi / L) a per mask cell (rad per km);
// - Phi, the least-squares potential of k over the mask: the edges between 4-neighbouring mask cells
//   (periodic east-west, none across the poles) ask Phi_n - Phi_c = the two cells' mean k times the step
//   on the ground, weighted by the finite-volume Laplacian on the sphere (dy / dx_row across a row,
//   dx_edge / dy between rows, dx = the column's 2 pi R / W times the cosine of the latitude, dy = pi R / H);
//   natural (Neumann) conditions at the mask's edge. Started from k integrated along a breadth-first tree
//   from each component's first cell, solved by over-relaxed Gauss-Seidel (factor 1.97) in double over the
//   mask cells in row-major order until the largest change of a sweep is below 0.01 rad (at most 20,000
//   sweeps: a tolerance reached deterministically, the sweep count reported), then each component's mean
//   set to 0;
// - two rings of global cells around the mask take an extrapolated phase (the mean over their defined
//   8-neighbours n of Phi_n + k_n . the step from n, ring by ring) and the mean k of those neighbours;
// - at a position (u, v) in global cell units (cell (i, j) spans [i, i + 1) x [j, j + 1)) the phase is
//   the bicubic Hermite interpolation over the 2 x 2 surrounding cell centres of Phi with the slopes
//   k_east dx_row and k_south dy and the cross term (the mean of d/dt (k_east dx_row) and d/ds (k_south dy)
//   by central differences between cells, one-sided at the defined set's edge) (foldphase_rad); the polar
//   half rows carry the edge row's curve by its slope; undefined where a corner is undefined. The misfit there
//   is m = |grad phase - k~| / |k~| with grad phase the interpolant's gradient on the ground (d/ds over the
//   column's east-west width at v, d/dt over dy) and k~ the bilinear interpolation of the corners' k.
//
// Per domain cell (land: not sea, not a channel), at its centre p in absolute global coordinates:
// - the cell's eps, L and height above the sea level are the global fields (eps, L, max(0, base
//   elevation - sea level)) interpolated by upsample, as the noise amplitude's relief is; the
//   cell's amplitude A = foldamplitude_m (the weight smoothstep(eps; strainmin, 2 strainmin)
//   included); a cell with A = 0 is not changed;
// - the gates (S5.3b): no fabric where inputs.glaciation (cryosphere.glaciation, 0 to 1) upsampled
//   and clamped to 0 .. 1 is above folds.icemax (icemax 1, or no field: no gate), nor in the block of
//   a global lake cell (lake surface above 0, not sea: regional_water.cpp's lake cells);
// - the cell is drawn only at the coarsest level whose cell is at most L / 4 (round two: at least 4 cells
//   per wavelength): the level's north-south cell pi R / (H k) <= L / 4 and the parent's pi R / (H
//   parentscale) > L / 4 (parentscale 0, level 1 and the single pass: no parent condition); finer levels
//   inherit the fabric by the prolongation, and a cell whose L / 4 is below the finest level's cell is
//   never drawn;
// - the fabric F = A sin(phase(p)) (no fabric where the phase is undefined), times the misfit fade
//   1 - smoothstep(m; misfitfadestart, misfitfadeend) (round two: where the phase cannot follow k the folds
//   end instead of bending), times the water gap 1 - max over the channel cells within gapfraction x L
//   (ground distance, rows at most k, columns columnextent) of 1 - smoothstep(d / (gapfraction x L));
// - the cell's elevation e becomes e + F, never lowered below S = the sea level + coast.shorefloor_m
//   (at least the next float above the sea level): max(e + F, S); a cell with e <= S (dry land at the
//   shore or below the sea level: a filled pocket, a basin) is not changed, since it could only rise
//   and the projection would have to sink its parent's other children (S5.2's first judging run:
//   24 forced parents and a land cell 751 m below the sea level at window A's bay). Channel and sea
//   cells are never changed; nothing else is read. The refinement loop's landbelowsea counts what the
//   projection may still push below the sea level.
// The water gap reads the level's channels within gapfraction x L, at most k rows (the cap), within
// formreach.

#include <cstdint>
#include <vector>

#include "simulations/landscape/channel_network.hpp"
#include "simulations/landscape/regional_grid.hpp"
#include "simulations/landscape/regional_types.hpp"

namespace landscape
{
// True when the fold fabric acts: folds.amplitudemax_m > 0 and inputs.foldstrain holds width x height values.
bool foldfabricon(const RegionalGlobalInputs& inputs, const RegionalParameters& parameters);

// L in km: clamp(wavelengthbase_km + wavelengthpercover_km x sediment_m / 1000 (0 when !withsediment), the smaller
// and the larger of wavelengthmin_km and wavelengthmax_km).
double foldwavelength_km(const FoldParameters& folds, double sediment_m, bool withsediment);

// A in m: min(relieffraction x (L / pi) sqrt(eps) (L in m), amplitudemax_m, reliefelevationfraction x
// max(0, heightabovesea_m)) x smoothstep(eps; strainmin, 2 strainmin); 0 when amplitudemax_m, eps or
// L is not positive.
double foldamplitude_m(const FoldParameters& folds, double strain, double wavelength_km, double heightabovesea_m);

// The phase field of the contract above on the global grid (S5.3c), W x H row-major, with what decided it.
struct FoldPhaseField
{
    int width = 0, height = 0;
    double dy_km = 0.0;                    // the global north-south cell on the ground
    double column_km = 0.0;                // a global column's east-west width at the equator
    std::vector<double> rowdx_km;          // the east-west cell of each row at its centre
    std::vector<double> phase_rad;         // Phi (0 where undefined)
    std::vector<double> wavenumbereast;    // k, rad per km (east, south; 0 where undefined)
    std::vector<double> wavenumbersouth;
    std::vector<double> cross_rad;         // the Hermite cross term, rad per cell squared
    std::vector<uint8_t> state;            // 0 undefined, 1 the mask (solved), 2 an extension ring
    std::vector<int8_t> orientation;       // the axis' sign on the mask (+1, -1), 0 elsewhere
    int maskcells = 0;
    int components = 0;                    // the mask's 4-connected components
    int orientedbyflood = 0;               // by the flood
    int orientedbynorth = 0;               // the components' first cells
    int sweeps = 0;                        // Gauss-Seidel sweeps
    double lastchange_rad = 0.0;           // the largest change of the last sweep
    double misfit = 0.0;                   // sqrt(sum w (dPhi - target)^2 / sum w target^2) over the edges
    int reversaledges = 0;                 // edges whose two cells' k point more than 90 degrees apart
    double seconds = 0.0;                  // wall time of the computation
};

// The phase field of `inputs` (empty, maskcells 0, without foldstrain or without a mask cell).
FoldPhaseField foldphasefield(const RegionalGlobalInputs& inputs, const FoldParameters& folds);

// The phase at (u, v) in global cell units (u wraps); false where it is undefined. misfit (optional): m there.
bool foldphase_rad(const FoldPhaseField& field, double u, double v, double& phase, double* misfit = nullptr);

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
