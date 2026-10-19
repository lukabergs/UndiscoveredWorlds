#include "simulations/landscape/regional_surface.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>

#include "core/deterministic_random.hpp"
#include "simulations/landscape/level_blocks.hpp"
#include "simulations/landscape/refinement_levels.hpp"
#include "simulations/landscape/regional_coast.hpp"
#include "simulations/landscape/regional_edges.hpp"
#include "simulations/landscape/regional_pass.hpp"

// Erosion plan steps E3, R0, R2 and R6 (docs/architecture/erosion-hydrology-plan.md, Phase E,
// Phase R and Phase R2): the Phase E start surface and sea moved from regional_pass.cpp by
// R0; the warp, the coast noise and the land-share rule of R2; R6 moved the sea gate to
// regional_edges.cpp (R8's) and the coast factors and the island sink to regional_coast.cpp
// (R7's), called from here with the R2 behaviour at their defaults.

namespace landscape
{
namespace
{
constexpr double pi = 3.14159265358979323846;
// Salt of the regional noise (deterministiccombine(seed, salt)).
constexpr std::uint64_t noisesalt = 0x5e6a;
// Salts of the warp's two components and of the coast noise (R2).
constexpr std::uint64_t warpxsalt = 0x3a7d;
constexpr std::uint64_t warpysalt = 0x3a7e;
constexpr std::uint64_t coastsalt = 0xc0a5;
// The form of the warp noise (not a tuning value): octave weights halve, as the legacy's
// diamond-square offsets and the plan's coast noise do.
constexpr double warppersistence = 0.5;

std::size_t at(int x, int y, int width)
{
    return static_cast<std::size_t>(y) * static_cast<std::size_t>(width) + static_cast<std::size_t>(x);
}

int wrapcolumn(int column, int width) { return ((column % width) + width) % width; }

std::int64_t wrapfine(std::int64_t column, std::int64_t width) { return ((column % width) + width) % width; }

// True when global cell (i, j) or one of its 8 neighbours satisfies `test` (columns
// wrap, rows beyond the poles do not exist).
template<typename Test>
bool aroundcell(int width, int height, int i, int j, Test test)
{
    for (int dj = -1; dj <= 1; dj++)
    {
        const int row = j + dj;
        if (row < 0 || row >= height)
            continue;
        for (int di = -1; di <= 1; di++)
            if (test(at(wrapcolumn(i + di, width), row, width)))
                return true;
    }
    return false;
}

// The local global relief: the largest absolute difference between a cell and its 8
// neighbours (columns wrap, rows beyond the poles do not exist).
std::vector<float> localrelief(const std::vector<float>& elevation, int width, int height)
{
    std::vector<float> relief(elevation.size(), 0.0f);
    for (int j = 0; j < height; j++)
        for (int i = 0; i < width; i++)
        {
            const float own = elevation[at(i, j, width)];
            float largest = 0.0f;
            aroundcell(width, height, i, j, [&](std::size_t other) {
                largest = std::max(largest, std::fabs(elevation[other] - own));
                return false;
            });
            relief[at(i, j, width)] = largest;
        }
    return relief;
}

// Whether Phase E's noise takes its amplitude from the world's relief amplitude (multiscale plan
// S3.4: regional.noise_relief_field 1 and inputs.reliefamplitude_m present).
bool relieffield(const RegionalGlobalInputs& inputs, const RegionalParameters& parameters)
{
    return parameters.noiserelieffield == 1 && !inputs.reliefamplitude_m.empty();
}

// The relief Phase E's noise amplitude grows with (noise_amplitude_m + noise_relief_fraction x it):
// the world's relief amplitude where relieffield, else the local global relief.
std::vector<float> noiserelief(const RegionalGlobalInputs& inputs, const RegionalParameters& parameters)
{
    if (relieffield(inputs, parameters))
        return inputs.reliefamplitude_m;
    return localrelief(inputs.elevation_m, inputs.width, inputs.height);
}

// The coastal indicator: 1 where the global cell's 3 x 3 neighbourhood holds a cell of
// the other kind (sea or land), else 0.
std::vector<float> coastalcells(const std::vector<uint8_t>& sea, int width, int height)
{
    std::vector<float> coastal(sea.size(), 0.0f);
    for (int j = 0; j < height; j++)
        for (int i = 0; i < width; i++)
        {
            const uint8_t own = sea[at(i, j, width)] != 0 ? 1 : 0;
            coastal[at(i, j, width)] =
                aroundcell(width, height, i, j, [&](std::size_t other) { return (sea[other] != 0 ? 1 : 0) != own; })
                    ? 1.0f
                    : 0.0f;
        }
    return coastal;
}

// A fine cell centre's latitude and longitude from its absolute indices: the expressions
// of RegionalWindow::finelatitude_rad and finelongitude_rad, so the values agree bit for bit.
double finelatitude(std::int64_t finerow, int fineheight)
{
    const int j = static_cast<int>(finerow);
    const double north = pi / 2.0 - pi * j / fineheight;
    const double south = pi / 2.0 - pi * (j + 1) / fineheight;
    return 0.5 * (north + south);
}

double finelongitude(std::int64_t finecolumn, int finewidth)
{
    return 2.0 * pi * (static_cast<double>(finecolumn) + 0.5) / finewidth;
}

// A fine cell centre in global cell units: (I + 0.5) / k (upsampling.cpp's expression).
double globalunits(std::int64_t fine, int scale) { return (static_cast<double>(fine) + 0.5) / scale; }

// The warp noise (regional_surface.hpp, regionalwarp): positions on the unit sphere in
// units of the first wavelength, wavelength_cells x pi / H.
class WarpField
{
public:
    WarpField(int worldwidth, int worldheight, int scale, std::uint64_t seed, const WarpParameters& parameters)
        : finewidth_(worldwidth * scale), fineheight_(worldheight * scale), amplitude_(parameters.amplitude_cells)
    {
        if (!(amplitude_ > 0.0) || !(parameters.wavelength_cells > 0.0))
            return;
        FractalNoise noise;
        noise.radius_km = 1.0;
        noise.wavelength_km = parameters.wavelength_cells * pi / worldheight;
        noise.persistence = warppersistence;
        noise.octaves = 1;
        while (parameters.wavelength_cells / std::ldexp(1.0, noise.octaves) >= 1.0)
            noise.octaves++;
        noise.seed = deterministiccombine(seed, warpxsalt);
        x_.emplace(noise);
        noise.seed = deterministiccombine(seed, warpysalt);
        y_.emplace(noise);
    }

    bool active() const { return x_.has_value(); }

    // The displacement in global cell units at the fine cell (I, J).
    void at(std::int64_t finecolumn, std::int64_t finerow, double& dx, double& dy) const
    {
        dx = dy = 0.0;
        if (!active())
            return;
        double px = 0.0, py = 0.0, pz = 0.0;
        unitvector(finelatitude(finerow, fineheight_), finelongitude(finecolumn, finewidth_), px, py, pz);
        const double u = x_->sum(px, py, pz) / x_->normalisation();
        const double v = y_->sum(px, py, pz) / y_->normalisation();
        const double length = std::sqrt(u * u + v * v);
        const double factor = length > 0.0 ? std::tanh(0.5 * length) / length : 0.5;
        dx = amplitude_ * factor * u;
        dy = amplitude_ * factor * v;
    }

private:
    int finewidth_ = 0;
    int fineheight_ = 0;
    double amplitude_ = 0.0;
    std::optional<FractalNoiseField> x_, y_;
};

// What the start surface and the decision of one fine cell are made of.
struct CellSample
{
    float lapsebase = 0.0f;  // the base elevation at the unwarped centre
    float base = 0.0f;       // the base elevation at the faded warped centre
    float surface = 0.0f;    // base + Phase E noise + the elevation's coast noise, before the block means
    float noise = 0.0f;      // the elevation's coast noise: the full one tapered to 0 at one amplitude from the sea level
    float fullnoise = 0.0f;  // the full coast noise (the decision's)
    float weight = 0.0f;     // the coastal weight (0 .. 1)
    float amplitude = 0.0f;  // the coast noise's amplitude before the coastal weight
    float inland = 0.0f;     // the global cells that may hold no sea, interpolated (0 .. 1; the land-share rule)
    float dropfactor = 1.0f; // multiplies the sea-side push (R7's coast factors; 1 without them)
};

// The per-cell part of the start surface, a pure function of the global fields, the seed
// and the fine cell's absolute position (and its full warp).
class Surfacer
{
public:
    Surfacer(const RegionalGlobalInputs& inputs, int scale, const RegionalParameters& parameters)
        : inputs_(inputs), parameters_(parameters), scale_(scale),
          relief_(localrelief(inputs.elevation_m, inputs.width, inputs.height)),
          noiserelief_(relieffield(inputs, parameters) ? inputs.reliefamplitude_m : std::vector<float>{}),
          coastal_(coastalcells(inputs.network.sea, inputs.width, inputs.height)),
          inland_(inlandcells(inputs, parameters)),
          factors_(coastfactors(inputs, parameters)),
          coast_(inputs.width, inputs.height, scale, inputs.seed, coastsalt, parameters.coast.noisepersistence)
    {
        // Phase E's noise (erosion plan E3): the band from twice the global north-south
        // spacing (exclusive) down to twice the fine one.
        FractalNoise noise;
        noise.seed = deterministiccombine(inputs.seed, noisesalt);
        noise.radius_km = inputs.radius_m / 1000.0;
        noise.longest_km = 2.0 * pi * inputs.radius_m / inputs.height / 1000.0;
        noise.wavelength_km = noise.longest_km;
        // Twice the fine north-south spacing, written as longest / k so that the octave
        // whose wavelength equals it (k a power of two) is not lost to rounding.
        noise.shortest_km = noise.longest_km / scale;
        noise.octaves = 1;
        while ((1 << noise.octaves) <= scale)
            noise.octaves++;
        noise.persistence = parameters.noisepersistence;
        band_.emplace(noise);
        coastnoise_ = parameters.coast.noiseamplitude_m > 0.0 && coast_.active();
    }

    const std::vector<float>& relief() const { return relief_; }

    CellSample sample(std::int64_t finecolumn, std::int64_t finerow, double warpx, double warpy) const
    {
        const int width = inputs_.width, height = inputs_.height;
        const double sealevel = static_cast<double>(inputs_.sealevel);
        const double gx = globalunits(finecolumn, scale_), gy = globalunits(finerow, scale_);
        CellSample cell;
        cell.lapsebase = static_cast<float>(bicubic(inputs_.elevation_m, width, height, gx, gy));
        double fade = 1.0;
        if (parameters_.warp.fade_m > 0.0)
            fade = std::clamp(1.0 - std::fabs(static_cast<double>(cell.lapsebase) - sealevel) / parameters_.warp.fade_m,
                              0.0, 1.0);
        const double wx = gx + fade * warpx, wy = gy + fade * warpy;
        cell.base = warpx == 0.0 && warpy == 0.0 ? cell.lapsebase
                                                 : static_cast<float>(bicubic(inputs_.elevation_m, width, height, wx, wy));
        const float relief = static_cast<float>(bicubic(relief_, width, height, gx, gy));
        // S3.4: Phase E's noise grows with the world's relief amplitude when it is used (the coast
        // noise keeps the local global relief).
        const float noiserelief =
            noiserelief_.empty() ? relief : static_cast<float>(bicubic(noiserelief_, width, height, gx, gy));

        // Phase E's noise, tapered at the sea level (regional_pass.cpp's expression of E3).
        float banded = cell.base;
        if (band_->resolved() > 0)
        {
            double px = 0.0, py = 0.0, pz = 0.0;
            unitvector(finelatitude(finerow, height * scale_), finelongitude(finecolumn, width * scale_), px, py, pz);
            const double amplitude =
                parameters_.noiseamplitude_m + parameters_.noiserelieffraction * std::max(0.0, static_cast<double>(noiserelief));
            const double height_m = std::fabs(static_cast<double>(cell.base) - sealevel);
            const double band = parameters_.coastband_m + amplitude;
            const double taper = band > 0.0 ? std::clamp(height_m / band, 0.0, 1.0) : 1.0;
            banded = static_cast<float>(static_cast<double>(cell.base) +
                                        taper * amplitude * band_->sum(px, py, pz) / band_->normalisation());
        }
        cell.surface = banded;

        // The coast noise on the coastal global cells: in full on the decision surface (not
        // tapered at the sea level: it sets the coastline); on the elevation, only without the
        // land-share rule (where it sets the coastline itself), no larger than Phase E's noise
        // at that height (its amplitude with its taper at the sea level).
        if (coastnoise_)
        {
            const double weight = std::clamp(bicubic(coastal_, width, height, wx, wy), 0.0, 1.0);
            double amplitude = parameters_.coast.noiseamplitude_m +
                               parameters_.coast.noiserelieffraction * std::max(0.0, static_cast<double>(relief));
            // R7's factors (regional_coast.hpp): the world's glaciation and denudation add to the
            // amplitude and weaken the sea-side push; not read while inactive (bit-identity).
            if (factors_.active)
            {
                amplitude += std::max(0.0, bicubic(factors_.amplitude_m, width, height, wx, wy));
                cell.dropfactor = static_cast<float>(std::clamp(bicubic(factors_.dropfactor, width, height, wx, wy), 0.0, 1.0));
            }
            if (weight > 0.0 && amplitude > 0.0)
            {
                const double inlandamplitude = std::max(
                    0.0, parameters_.noiseamplitude_m + parameters_.noiserelieffraction * std::max(0.0, static_cast<double>(noiserelief)));
                const double height_m = std::fabs(static_cast<double>(cell.base) - sealevel);
                const double band = parameters_.coastband_m + inlandamplitude;
                const double inlandtaper = band > 0.0 ? std::clamp(height_m / band, 0.0, 1.0) : 1.0;
                const double elevationamplitude = std::min(amplitude, inlandtaper * inlandamplitude);
                const double full = weight * amplitude * coast_.at(finecolumn, finerow);
                cell.fullnoise = static_cast<float>(full);
                // With the land-share rule the coastline is the decision's, and the elevation takes none.
                cell.noise = parameters_.coast.landshare > 0.0 ? 0.0f : static_cast<float>((elevationamplitude / amplitude) * full);
                cell.weight = static_cast<float>(weight);
                cell.amplitude = static_cast<float>(amplitude);
                cell.surface = static_cast<float>(static_cast<double>(banded) + static_cast<double>(cell.noise));
            }
        }
        if (parameters_.coast.landshare > 0.0)
        {
            cell.inland = static_cast<float>(std::clamp(bicubic(inland_, width, height, wx, wy), 0.0, 1.0));
        }
        return cell;
    }

private:
    const RegionalGlobalInputs& inputs_;
    const RegionalParameters& parameters_;
    int scale_ = 1;
    std::vector<float> relief_;
    std::vector<float> noiserelief_;  // the world's relief amplitude when relieffield, else empty (relief_ is used)
    std::vector<float> coastal_;
    std::vector<float> inland_;
    CoastFactors factors_;
    RegionalNoise coast_;
    std::optional<FractalNoiseField> band_;
    bool coastnoise_ = false;
};

// The elevation made to agree with the mask on the gated cells (regional_edges.hpp,
// seagate): a sea cell above the sea level lowered to it, a land cell at or below it
// raised to the next float above it.
void clamptomask(const RegionalGlobalInputs& inputs, const RegionalWindow& window, const std::vector<uint8_t>& gate,
                 const std::vector<uint8_t>& sea, std::vector<float>& surface)
{
    const int domainwidth = window.domainfinewidth(), domainheight = window.domainfineheight();
    const float abovesea = std::nextafter(inputs.sealevel, std::numeric_limits<float>::infinity());
    for (int y = 0; y < domainheight; y++)
    {
        for (int x = 0; x < domainwidth; x++)
        {
            const std::size_t cell = at(x, y, domainwidth);
            if (gate[cell] == 0)
                continue;
            if (sea[cell] != 0)
                surface[cell] = std::min(surface[cell], inputs.sealevel);
            else if (surface[cell] <= inputs.sealevel)
                surface[cell] = abovesea;
        }
    }
}

// One pass of the block-mean correction that keeps every fine cell on its side of the
// sea level (the land-share rule; Phase E's additive correction, upsampling.hpp, would
// lift a shallow bank between deep cells above the sea level, or sink a low coast below
// it, and a clamp after it converges slowly). The correction is Phase E's: the residual
// of every domain block (global value minus block mean) interpolated by `bicubic` at the
// fine centres, the residuals beyond the domain's edge clamped as in correctblockmeans
// (its reach, two global cells per pass); but a cell of height h above the sea level
// moves by it at most until signedlimit of h is gone (h' = max(h + c, (1 - signedlimit) h)
// on land, min(...) at sea), so no sign changes and the rest is left to the next pass.
constexpr double signedlimit = 0.9;

// A fine value moved by a correction under that limit.
float signedmove(float value, double correction, double sealevel)
{
    const double h = static_cast<double>(value) - sealevel;
    const double moved = h > 0.0 ? std::max(h + correction, (1.0 - signedlimit) * h)
                                 : std::min(h + correction, (1.0 - signedlimit) * h);
    return static_cast<float>(sealevel + moved);
}

void correctsignedblockmeans(const RegionalGlobalInputs& inputs, const RegionalWindow& window, std::vector<float>& fine)
{
    const int width = inputs.width, height = inputs.height;
    const int k = window.scale;
    const int finewidth = window.domainfinewidth(), fineheight = window.domainfineheight();
    const int dcolumns = window.domaincolumns, drows = window.domainrows;
    const double sealevel = static_cast<double>(inputs.sealevel);
    // The domain's blocks over the per-face quadtree (level_blocks.hpp; keys only, the radius nominal),
    // in the domain's row-major order of global cells, as correctblockmeans takes them.
    const WindowBlocks blocks = windowblocks(*regionalwindowgrid(window, 1.0), k);
    if (blocks.count() != static_cast<std::size_t>(dcolumns) * static_cast<std::size_t>(drows))
        throw std::logic_error("correctsignedblockmeans: the domain's blocks are not its global cells");
    std::vector<double> residual(blocks.count());
    for (std::size_t b = 0; b < blocks.count(); b++)
        // Summed in ascending key order inside the block (its row-major order), as correctblockmeans sums.
        residual[b] = static_cast<double>(inputs.elevation_m[blocks.parents[b]]) - blockmean(blocks, b, fine);

    // On the global grid for `bicubic`, the two cells beyond the domain's edges taking the
    // nearest domain cell's value (correctblockmeans' ring).
    std::vector<float> field(static_cast<std::size_t>(width) * static_cast<std::size_t>(height), 0.0f);
    const int rowfirst = std::max(0, window.domainrow0 - 2);
    const int rowlast = std::min(height - 1, window.domainrow0 + drows + 1);
    const bool wholeturn = dcolumns >= width;
    for (int row = rowfirst; row <= rowlast; row++)
    {
        const int dr = std::max(0, std::min(drows - 1, row - window.domainrow0));
        for (int rel = -2; rel < dcolumns + 2; rel++)
        {
            const bool outside = rel < 0 || rel >= dcolumns;
            if (outside && wholeturn)
                continue;
            const int column = wrapcolumn(window.domaincolumn0 + rel, width);
            if (outside && wrapcolumn(column - window.domaincolumn0, width) < dcolumns)
                continue;
            const int dc = std::max(0, std::min(dcolumns - 1, rel));
            field[at(column, row, width)] = static_cast<float>(residual[at(dc, dr, dcolumns)]);
        }
    }
    for (int y = 0; y < fineheight; y++)
    {
        const double gy = globalunits(window.absolutefinerow(y), k);
        for (int x = 0; x < finewidth; x++)
        {
            const double gx = globalunits(window.absolutefinecolumn(x), k);
            const double correction = bicubic(field, width, height, gx, gy);
            float& value = fine[at(x, y, finewidth)];
            value = signedmove(value, correction, sealevel);
        }
    }
}

// The start surface under the land-share rule (regional_surface.hpp). Everything is
// computed on a box, the domain widened by `margin` global cells (rows stop at the poles;
// the whole turn when the widened domain would reach around the planet), described as a
// RegionalWindow whose domain is the box: the samples (Surfacer, with the canonical warp),
// Phase E's block-mean correction, the decision surface, the land-share offset, the mask,
// the clamp to it and the sign-keeping passes. A value that depends on the box's edge
// (the residual clamp of the block means, the offset's clamp) is at most
// 2 meanpasses + r sharepasses + 2 + 2 signedpasses global cells from it, where r is a share
// pass's reach: the offset's bilinear field reaches one global cell and the mask, after
// the pinch sweep and the hole fill, reads two fine cells beyond a block, which stay in the
// neighbour's half towards it for k >= 3 (r = 1) and reach the next one for k <= 2 (r = 2).
// With margin = r sharepasses + 2 + 2 signedpasses the domain's values depend on the window
// only within 2 meanpasses of the domain's edge: Phase E's reach, and the halo rule
// (halo >= 2 meanpasses + 2) stays sufficient. The passes stop early when no block of the
// box misses its share; a pass that changes nothing within the box is the same as no pass,
// and changes beyond the box travel one reach per pass, so stopping early keeps the result.
// R6: the cells that may be sea come from the gate (regional_edges.hpp, seagate: the 3 x 3
// rule per global cell until R8), the mask passes through the island sink
// (regional_coast.hpp) and the straight-run breaker (regional_edges.hpp), both no-ops at
// their defaults, and the coast factors scale the push (Surfacer).
void landsharesurface(const RegionalGlobalInputs& inputs, const RegionalWindow& window, const RegionalParameters& parameters,
                      const Surfacer& surfacer, StartSurface& start, bool lakegate)
{
    const int width = inputs.width, height = inputs.height;
    const int k = window.scale;
    const int sharecount = std::max(0, parameters.coast.sharepasses);
    const int signedcount = std::max(0, parameters.coast.signedpasses);
    // R6: a mask hook (the island sink's ring of m fine cells, regional_coast.hpp; the
    // straight-run breaker, which reads its block and one cell beyond, regional_edges.hpp)
    // widens what a block's share reads beyond the block: the final pinch sweep and hole fill
    // read two cells into the neighbouring block, whose hook read that whole block and its
    // ring on the first-cleaned mask, so e = k + hookreach + 2 fine cells (the second clean's
    // two cells lie inside the neighbour, already counted by k), and a share pass's reach is
    // one block (the offset's tent) plus the blocks e exceeds half a block by: r = 1 +
    // ceil(max(0, e - k / 2) / k) (the R6 and R7 reviews' derivation: r = 2 at k 16 with m 2,
    // 3 at k 4 with m 1). Without a hook the R2 rule stands.
    const int hookreach = std::max(parameters.coast.islandmargin > 0.0 ? std::max(1, static_cast<int>(std::lround(parameters.coast.islandmargin * k))) : 0,
                                   parameters.edges.straightmax > 0 ? 1 : 0);
    const int sharereach = hookreach > 0 ? 1 + (std::max(0, k + hookreach + 2 - k / 2) + k - 1) / k : (k <= 2 ? 2 : 1);
    // S3.15: the shoreline's displacement reads the decision surface within warpreach global cells
    // (its bound and the bilinear stencil, ShorelineWarp::reachcells) before the share passes, and
    // the pocket fill the final mask within coast.pocketreach global cells after them, so the
    // margin grows by both (0 and 0 without them: the R2 box).
    const ShorelineWarp shorewarp(width, height, k, 1, inputs.seed, parameters.coast);
    const int warpreach = (shorewarp.reachcells() + k - 1) / k;
    const int pocketreach = std::max(0, parameters.coast.pocketreach);
    const int margin = sharereach * sharecount + 2 + 2 * signedcount + warpreach + pocketreach;
    const double sealevel = static_cast<double>(inputs.sealevel);

    RegionalWindow box = window;
    const bool wholeturn = window.domaincolumns + 2 * margin >= width;
    box.domaincolumns = wholeturn ? width : window.domaincolumns + 2 * margin;
    box.domaincolumn0 = wholeturn ? window.domaincolumn0 : wrapcolumn(window.domaincolumn0 - margin, width);
    box.domainrow0 = std::max(0, window.domainrow0 - margin);
    box.domainrows = std::min(height - 1, window.domainrow0 + window.domainrows - 1 + margin) - box.domainrow0 + 1;
    const int columns = box.domaincolumns, rows = box.domainrows;
    const int finewidth = box.domainfinewidth(), fineheight = box.domainfineheight();
    const std::size_t finecount = static_cast<std::size_t>(finewidth) * static_cast<std::size_t>(fineheight);
    const int domainx = wholeturn ? 0 : margin * k, domainy = (window.domainrow0 - box.domainrow0) * k;
    const int domainwidth = window.domainfinewidth(), domainheight = window.domainfineheight();
    auto globalof = [&](int bx, int by) { return at(wrapcolumn(box.domaincolumn0 + bx, width), box.domainrow0 + by, width); };

    // The samples and Phase E's block means on the box.
    const FineWarp warp = regionalwarp(box, inputs.seed, parameters);
    // S3.4: in the refinement loop the gate's lake exclusion is the lake cells' warped half-mask.
    const std::vector<float> lakeweight = lakegate ? lakecellweight(inputs, box, warp) : std::vector<float>{};
    const std::vector<uint8_t> gate = seagate(inputs, box, warp, parameters, lakegate ? &lakeweight : nullptr);
    // A block may hold sea when any of its cells may (with the 3 x 3 gate: its global cell's
    // nearseacells value); the land-share rule and the shore apply to those blocks only.
    std::vector<uint8_t> blockgate(static_cast<std::size_t>(columns) * static_cast<std::size_t>(rows), 0);
    for (int y = 0; y < fineheight; y++)
        for (int x = 0; x < finewidth; x++)
            if (gate[at(x, y, finewidth)] != 0)
                blockgate[at(x / k, y / k, columns)] = 1;
    const bool onesided = parameters.coast.onesided != 0;
    const bool maskhooks = parameters.coast.islandmargin > 0.0 || parameters.edges.straightmax > 0;
    std::vector<float> surface(finecount), noise(finecount), fullnoise(finecount), weight(finecount),
        amplitude(finecount), inland(finecount), dropfactor(finecount);
    const std::size_t domaincount = static_cast<std::size_t>(domainwidth) * static_cast<std::size_t>(domainheight);
    start.base_m.assign(domaincount, 0.0f);
    start.lapsebase_m.assign(domaincount, 0.0f);
    for (int y = 0; y < fineheight; y++)
        for (int x = 0; x < finewidth; x++)
        {
            const std::size_t cell = at(x, y, finewidth);
            const CellSample sample = surfacer.sample(box.absolutefinecolumn(x), box.absolutefinerow(y),
                                                      warp.empty() ? 0.0 : warp.dx[cell], warp.empty() ? 0.0 : warp.dy[cell]);
            surface[cell] = sample.surface;
            noise[cell] = sample.noise;
            fullnoise[cell] = sample.fullnoise;
            weight[cell] = sample.weight;
            amplitude[cell] = sample.amplitude;
            inland[cell] = sample.inland;
            dropfactor[cell] = sample.dropfactor;
            const int dx = x - domainx, dy = y - domainy;
            if (dx >= 0 && dx < domainwidth && dy >= 0 && dy < domainheight)
            {
                start.base_m[at(dx, dy, domainwidth)] = sample.base;
                start.lapsebase_m[at(dx, dy, domainwidth)] = sample.lapsebase;
            }
        }
    correctblockmeans(inputs.elevation_m, width, height, box, parameters.meanpasses, surface);

    // The decision surface: the corrected surface with the full coast noise in place of the
    // elevation's tapered one, its sea side pushed down (a height h below the sea level
    // before the coast noise sinks by seadrop x coastal weight x min(1, -h / amplitude), so
    // the push fades with the weight) and lifted by seadrop x inland. R7: the push is
    // multiplied by the coast factors' drop factor (1 without them), and with coast.onesided
    // the noise inside a global sea cell may not lift the decision above the sea level: where
    // it would, the cell takes the decision without the noise (the legacy's rule).
    std::vector<double> decision(finecount);
    // S3.15 (b): the core sea of the pocket rule, the cells of global sea cells whose decision without
    // the coast noise (and before the displacement below) lies at or below the sea level (the one-sided
    // rule's `base`).
    std::vector<uint8_t> core(pocketreach > 0 ? finecount : 0, 0);
    for (int y = 0; y < fineheight; y++)
        for (int x = 0; x < finewidth; x++)
        {
            const std::size_t cell = at(x, y, finewidth);
            const double value = static_cast<double>(surface[cell]);
            const double h = value - static_cast<double>(noise[cell]) - sealevel;
            const double local = static_cast<double>(amplitude[cell]);
            double push = h < 0.0 && local > 0.0
                              ? parameters.coast.seadrop_m * static_cast<double>(weight[cell]) * std::min(1.0, -h / local)
                              : 0.0;
            if (dropfactor[cell] != 1.0f)
                push *= static_cast<double>(dropfactor[cell]);
            // R2's expression, kept textually (the R6 review: a re-associated sum rounds differently).
            decision[cell] = value - static_cast<double>(noise[cell]) + static_cast<double>(fullnoise[cell]) - push +
                             parameters.coast.seadrop_m * static_cast<double>(inland[cell]);
            if (onesided && inputs.network.sea[globalof(x / k, y / k)] != 0)
            {
                const double base = value - static_cast<double>(noise[cell]) - push + parameters.coast.seadrop_m * static_cast<double>(inland[cell]);
                if (base <= sealevel && decision[cell] > sealevel)
                    decision[cell] = base;
            }
            if (!core.empty() && inputs.network.sea[globalof(x / k, y / k)] != 0)
            {
                const double base = value - static_cast<double>(noise[cell]) - push + parameters.coast.seadrop_m * static_cast<double>(inland[cell]);
                core[cell] = base <= sealevel ? 1 : 0;
            }
        }
    // S3.15 (a): the shoreline's position, the decision sampled at the displaced position (the
    // elevation keeps its own warp and noise; the share rule, the passes and the shore read this one).
    if (shorewarp.active())
        decision = displacedsample(decision, box, shorewarp);

    // The offset: one value per block, interpolated bilinearly between block centres
    // (clamped at the box's edges, wrapped around a whole turn), so it is continuous
    // across block edges. Each Jacobi pass gives every block that misses its share the
    // increment that would meet it if its neighbours' offsets stayed (its own weight at a
    // cell is the bilinear tent, 0.25 to 1).
    std::vector<double> offset(static_cast<std::size_t>(columns) * static_cast<std::size_t>(rows), 0.0);
    std::vector<double> field(finecount, 0.0);
    auto blockcolumn = [&](int b) { return wholeturn ? wrapcolumn(b, columns) : std::clamp(b, 0, columns - 1); };
    // The bilinear weights from the fine cell's position inside its global cell (o = the
    // absolute fine index mod k, the same in every window): (o + 0.5) / k - 0.5 is the offset
    // from the cell's centre in global cells; below 0 the cell lies between the previous
    // centre and its own.
    std::vector<int> lower(static_cast<std::size_t>(k));
    std::vector<double> fraction(static_cast<std::size_t>(k));
    for (int o = 0; o < k; o++)
    {
        const double t = (o + 0.5) / k - 0.5;
        lower[static_cast<std::size_t>(o)] = t < 0.0 ? -1 : 0;
        fraction[static_cast<std::size_t>(o)] = t < 0.0 ? t + 1.0 : t;
    }
    auto interpolate = [&]() {
        for (int y = 0; y < fineheight; y++)
        {
            const std::size_t oy = static_cast<std::size_t>(y % k);
            const int rowbase = y / k + lower[oy];
            const double ty = fraction[oy];
            const int r0 = std::clamp(rowbase, 0, rows - 1), r1 = std::clamp(rowbase + 1, 0, rows - 1);
            for (int x = 0; x < finewidth; x++)
            {
                const std::size_t ox = static_cast<std::size_t>(x % k);
                const int columnbase = x / k + lower[ox];
                const double tx = fraction[ox];
                const int c0 = blockcolumn(columnbase), c1 = blockcolumn(columnbase + 1);
                const double north = offset[at(c0, r0, columns)] + tx * (offset[at(c1, r0, columns)] - offset[at(c0, r0, columns)]);
                const double south = offset[at(c0, r1, columns)] + tx * (offset[at(c1, r1, columns)] - offset[at(c0, r1, columns)]);
                field[at(x, y, finewidth)] = north + ty * (south - north);
            }
        }
    };
    // The mask on the box (1 = sea), then the one-cell holes filled (a sea cell whose four
    // neighbours are land; a neighbour beyond the box or a pole does not count). The
    // shares are counted on the filled mask.
    std::vector<uint8_t> raw(finecount, 0), mask(finecount, 0);
    // The pinch sweep and the hole fill, raw -> mask (below); R6: run again after the mask
    // hooks (the island sink and the straight-run breaker) when either is on, on the mask
    // they left, so the shares are counted on a mask without pinches or holes.
    auto cleanmask = [&]() {
        // Diagonal pinches first: in a 2 x 2 square whose land and sea cells touch only
        // diagonally, the sea cells become land (all squares at once, `mask` as the scratch).
        std::fill(mask.begin(), mask.end(), 0);
        for (int y = 0; y + 1 < fineheight; y++)
            for (int x = 0; x < finewidth; x++)
            {
                const int east = wholeturn ? wrapcolumn(x + 1, finewidth) : x + 1;
                if (east >= finewidth)
                    continue;
                const uint8_t a = raw[at(x, y, finewidth)], b = raw[at(east, y, finewidth)];
                const uint8_t c = raw[at(x, y + 1, finewidth)], d = raw[at(east, y + 1, finewidth)];
                if (a == d && b == c && a != b)
                {
                    mask[at(a != 0 ? x : east, y, finewidth)] = 1;
                    mask[at(a != 0 ? east : x, y + 1, finewidth)] = 1;
                }
            }
        for (std::size_t cell = 0; cell < finecount; cell++)
            if (mask[cell] != 0)
                raw[cell] = 0;
        // Then the one-cell holes.
        for (int y = 0; y < fineheight; y++)
            for (int x = 0; x < finewidth; x++)
            {
                if (raw[at(x, y, finewidth)] == 0)
                {
                    mask[at(x, y, finewidth)] = 0;
                    continue;
                }
                bool hole = true;
                const int west = wholeturn ? wrapcolumn(x - 1, finewidth) : x - 1;
                const int east = wholeturn ? wrapcolumn(x + 1, finewidth) : x + 1;
                const int nx[4] = {west, east, x, x};
                const int ny[4] = {y, y, y - 1, y + 1};
                for (int n = 0; n < 4 && hole; n++)
                    if (nx[n] >= 0 && nx[n] < finewidth && ny[n] >= 0 && ny[n] < fineheight)
                        hole = raw[at(nx[n], ny[n], finewidth)] == 0;
                mask[at(x, y, finewidth)] = hole ? 0 : 1;
            }
    };
    auto makemask = [&]() {
        interpolate();
        for (int y = 0; y < fineheight; y++)
            for (int x = 0; x < finewidth; x++)
                raw[at(x, y, finewidth)] =
                    gate[at(x, y, finewidth)] != 0 && decision[at(x, y, finewidth)] + field[at(x, y, finewidth)] <= sealevel ? 1 : 0;
        cleanmask();
        if (maskhooks)
        {
            sinksmallislands(inputs, box, parameters, gate, mask);
            breakstraightruns(inputs, box, parameters, gate, mask);
            raw = mask;
            cleanmask();
        }
    };
    const int blockcells = k * k;
    const double share = parameters.coast.landshare;
    const int leastland = static_cast<int>(std::ceil(share * blockcells - 1.0e-9));
    const int mostland = static_cast<int>(std::floor((1.0 - share) * blockcells + 1.0e-9));
    std::vector<double> increment(offset.size(), 0.0);
    std::vector<double> keys;
    keys.reserve(static_cast<std::size_t>(blockcells));
    const int domainblock0 = domainx / k, domainblockrow0 = domainy / k;
    for (int pass = 0;; pass++)
    {
        // The mask of the offsets so far, verified: every block that misses its share counts
        // (those of the domain apart); none missing in the box, or the cap reached, ends it.
        makemask();
        std::fill(increment.begin(), increment.end(), 0.0);
        int misses = 0, domainmisses = 0;
        for (int by = 0; by < rows; by++)
            for (int bx = 0; bx < columns; bx++)
            {
                const std::size_t global = globalof(bx, by);
                if (blockgate[at(bx, by, columns)] == 0)
                    continue;
                const bool seacell = inputs.network.sea[global] != 0;
                int land = 0;
                for (int y = by * k; y < (by + 1) * k; y++)
                    for (int x = bx * k; x < (bx + 1) * k; x++)
                        land += mask[at(x, y, finewidth)] == 0 ? 1 : 0;
                const int missing = seacell ? land - mostland : leastland - land;
                if (missing <= 0)
                    continue;
                misses++;
                if (bx >= domainblock0 && bx < domainblock0 + window.domaincolumns && by >= domainblockrow0 &&
                    by < domainblockrow0 + window.domainrows)
                    domainmisses++;
                // The offset that turns `missing` cells, cheapest first: a cell's cost is its
                // distance from the sea level divided by the block's own tent weight there.
                keys.clear();
                for (int y = by * k; y < (by + 1) * k; y++)
                {
                    const double wy = 1.0 - std::fabs((y - by * k + 0.5) / k - 0.5);
                    for (int x = bx * k; x < (bx + 1) * k; x++)
                    {
                        const double wx = 1.0 - std::fabs((x - bx * k + 0.5) / k - 0.5);
                        const double value = decision[at(x, y, finewidth)] + field[at(x, y, finewidth)] - sealevel;
                        // Filled holes (sea by value, land by the fill) are not candidates.
                        if (seacell ? value > 0.0 : mask[at(x, y, finewidth)] != 0)
                            keys.push_back(std::fabs(value) / (wx * wy));
                    }
                }
                if (keys.empty())
                    continue;
                // One cell beyond the share, so that the neighbours' increments of the same pass
                // rarely leave it a cell short.
                const std::size_t pick = std::min(keys.size(), static_cast<std::size_t>(missing) + 1) - 1;
                std::nth_element(keys.begin(), keys.begin() + static_cast<std::ptrdiff_t>(pick), keys.end());
                const double needed = keys[pick] * 1.001 + 0.001;
                increment[at(bx, by, columns)] = seacell ? -needed : needed;
            }
        if (domainmisses == 0 && start.sharepasses < 0)
            start.sharepasses = pass;
        start.sharemisses = domainmisses;
        if (misses == 0 || pass == sharecount)
            break;
        for (std::size_t b = 0; b < offset.size(); b++)
            offset[b] += increment[b];
    }

    // S3.15 (b): the sea pockets of the final mask (regional_coast.hpp, seapockets, anchored on the core
    // sea) filled after the shares are counted, a global sea cell's only while its block keeps at most
    // 1 - landshare land (a global land cell's fill only adds land), and a filled 4-connected component
    // leaves no pinch or one-cell hole behind. The pocket's basin (pocketbasins: the land at or below the
    // sea level around it within the reach) takes the shore floor in full (below). A basin bounded by the
    // sea level + coast.shorefloor_m instead spread over the low plains and raised them to 20 m: on window A
    // at k 64 a 67 m block-mean error at global cell (189, 46), a 96 m border step at level 1 and 494 m at
    // level 2.
    std::vector<uint8_t> basin;
    if (pocketreach > 0)
    {
        std::vector<uint8_t> filled = seapockets(inputs, box, pocketreach, mask, core);
        for (int by = 0; by < rows; by++)
            for (int bx = 0; bx < columns; bx++)
            {
                if (inputs.network.sea[globalof(bx, by)] == 0)
                    continue;
                int land = 0, fill = 0;
                for (int y = by * k; y < (by + 1) * k; y++)
                    for (int x = bx * k; x < (bx + 1) * k; x++)
                    {
                        land += mask[at(x, y, finewidth)] == 0 ? 1 : 0;
                        fill += filled[at(x, y, finewidth)] != 0 ? 1 : 0;
                    }
                if (fill > 0 && land + fill > mostland)
                    for (int y = by * k; y < (by + 1) * k; y++)
                        for (int x = bx * k; x < (bx + 1) * k; x++)
                            filled[at(x, y, finewidth)] = 0;
            }
        for (std::size_t cell = 0; cell < finecount; cell++)
            if (filled[cell] != 0)
                mask[cell] = 0;
        const int keptx = domainx + window.keptfinecolumnoffset(), kepty = domainy + window.keptfinerowoffset();
        const PocketCount done = countpockets(filled, finewidth, keptx, kepty, window.finewidth(), window.fineheight());
        const PocketCount left = countpockets(seapockets(inputs, box, pocketreach, mask, core), finewidth, keptx, kepty,
                                              window.finewidth(), window.fineheight());
        start.pocketsfilled = done.pockets;
        start.pocketcellsfilled = done.pocketcells;
        start.pockets = left.pockets;
        start.pocketcells = left.pocketcells;
        basin = pocketbasins(inputs, box, pocketreach, mask, filled, surface,
                             std::nextafter(inputs.sealevel, std::numeric_limits<float>::infinity()));
        // In global land cells only: a pocket filled in a global sea cell takes the ordinary shore (the clamp
        // above the sea level), so the sea block's mean stays within the sign-keeping passes' reach.
        for (int y = 0; y < fineheight; y++)
            for (int x = 0; x < finewidth; x++)
                if (inputs.network.sea[globalof(x / k, y / k)] != 0)
                    basin[at(x, y, finewidth)] = 0;
    }

    // The elevation made to agree with the mask and given a shore, then its block means
    // corrected without crossing the sea level, the shore kept after every pass; the domain's
    // part kept. The shore: a sea cell is at least as deep as the decision surface there and a
    // land cell at least as high, both capped at coast.shorefloor_m and at the distance of the
    // cell's own global cell from the sea level (a pure function of the global cell, so a
    // block whose global value lies near the sea level can still reach its mean), so the depth
    // and the height grow from the shoreline and no zero-depth sheet stays, while the land the
    // decision adds is low and as smooth as inland land (a deeper shore, the push of up to
    // seadrop_m, broke shallow shelves' block means: 76 m on window A; an uncapped 20 m floor
    // 29 m on window C; a cap interpolated between global cells left 569 zero-depth cells on
    // A). Then the clamp (clamptomask).
    std::vector<float> shoredepth(finecount, 0.0f), shoreheight(finecount, 0.0f), shorefloor(finecount, 0.0f);
    for (int y = 0; y < fineheight; y++)
        for (int x = 0; x < finewidth; x++)
        {
            const std::size_t cell = at(x, y, finewidth);
            const std::size_t global = globalof(x / k, y / k);
            // S3.15: a filled pocket's basin, gated or not, raised to the shore floor in full.
            if (!basin.empty() && basin[cell] != 0 && mask[cell] == 0)
            {
                shoreheight[cell] = static_cast<float>(std::max(
                    0.0, std::min(parameters.coast.shorefloor_m, std::fabs(static_cast<double>(inputs.elevation_m[global]) - sealevel))));
                continue;
            }
            if (gate[cell] == 0)
                continue;
            const double value = decision[cell] + field[cell] - sealevel;
            const double floor = std::max(
                0.0, std::min(parameters.coast.shorefloor_m, std::fabs(static_cast<double>(inputs.elevation_m[global]) - sealevel)));
            shorefloor[cell] = static_cast<float>(floor);
            // A cell the mask moves across the sea level takes the floor in full; the others
            // only in a band along the shoreline (full while the decision's distance from the
            // sea level is within the floor, fading to nothing at twice it), so the shelf and
            // the land beyond keep the corrected surface (a floor everywhere flattened shelves
            // into block-shaped plateaus).
            const double reach = std::fabs(value);
            const double band = floor > 0.0 ? std::clamp(2.0 - reach / floor, 0.0, 1.0) : 0.0;
            const bool crossed = mask[cell] != 0 ? surface[cell] > inputs.sealevel : surface[cell] <= inputs.sealevel;
            const double fade = crossed ? 1.0 : band;
            // R7: a sea cell whose decision is land (only the mask hooks make one: the island
            // sink) takes the floor in full, so sunk land leaves no zero-depth sheet.
            if (mask[cell] != 0)
                shoredepth[cell] = static_cast<float>(fade * (value > 0.0 ? floor : std::max(0.0, std::min(-value, floor))));
            else
                shoreheight[cell] = static_cast<float>(fade * std::max(0.0, std::min(value, floor)));
        }
    auto keepshore = [&]() {
        for (std::size_t cell = 0; cell < finecount; cell++)
        {
            if (mask[cell] != 0 && shoredepth[cell] > 0.0f)
                surface[cell] = std::min(surface[cell], static_cast<float>(sealevel - static_cast<double>(shoredepth[cell])));
            else if (mask[cell] == 0 && shoreheight[cell] > 0.0f)
                surface[cell] = std::max(surface[cell], static_cast<float>(sealevel + static_cast<double>(shoreheight[cell])));
        }
        clamptomask(inputs, box, gate, mask, surface);
    };
    keepshore();
    // The sign-keeping passes may not bring a cell nearer the sea level than it is now or
    // than its floor, whichever is less: repeated passes would otherwise lift shallow sea
    // towards a zero-depth sheet (and lower low land towards a sea-level one).
    for (std::size_t cell = 0; cell < finecount; cell++)
    {
        const double h = static_cast<double>(surface[cell]) - sealevel;
        const double keep = std::min(std::fabs(h), static_cast<double>(shorefloor[cell]));
        if (mask[cell] != 0)
            shoredepth[cell] = std::max(shoredepth[cell], static_cast<float>(keep));
        else if (h > 0.0)
            shoreheight[cell] = std::max(shoreheight[cell], static_cast<float>(keep));
    }
    for (int pass = 0; pass < signedcount; pass++)
    {
        correctsignedblockmeans(inputs, box, surface);
        keepshore();
    }
    start.start_m.assign(domaincount, 0.0f);
    for (int y = 0; y < domainheight; y++)
        for (int x = 0; x < domainwidth; x++)
            start.start_m[at(x, y, domainwidth)] = surface[at(x + domainx, y + domainy, finewidth)];
}
}

RegionalNoise::RegionalNoise(int worldwidth, int worldheight, int scale, std::uint64_t seed, std::uint64_t salt,
                             double persistence)
    : finewidth_(worldwidth * scale), fineheight_(worldheight * scale)
{
    int octaves = 0;
    while ((2 << octaves) <= scale)
        octaves++;
    if (octaves < 1 || worldheight < 1)
        return;
    FractalNoise noise;
    noise.seed = deterministiccombine(seed, salt);
    noise.radius_km = 1.0;
    noise.wavelength_km = pi / worldheight;
    noise.octaves = octaves;
    noise.persistence = persistence;
    field_.emplace(noise);
}

double RegionalNoise::at(std::int64_t finecolumn, std::int64_t finerow) const
{
    if (!field_)
        return 0.0;
    double px = 0.0, py = 0.0, pz = 0.0;
    unitvector(finelatitude(finerow, fineheight_), finelongitude(wrapfine(finecolumn, finewidth_), finewidth_), px, py, pz);
    return field_->sum(px, py, pz) / field_->normalisation();
}

FineWarp regionalwarp(const RegionalWindow& window, std::uint64_t seed, const RegionalParameters& parameters)
{
    const WarpField field(window.worldwidth, window.worldheight, window.scale, seed, parameters.warp);
    FineWarp warp;
    if (!field.active())
        return warp;
    const int domainwidth = window.domainfinewidth(), domainheight = window.domainfineheight();
    const std::size_t count = static_cast<std::size_t>(domainwidth) * static_cast<std::size_t>(domainheight);
    warp.dx.resize(count);
    warp.dy.resize(count);
    for (int y = 0; y < domainheight; y++)
        for (int x = 0; x < domainwidth; x++)
        {
            const std::size_t cell = at(x, y, domainwidth);
            field.at(window.absolutefinecolumn(x), window.absolutefinerow(y), warp.dx[cell], warp.dy[cell]);
        }
    return warp;
}

StartSurface regionalstartsurface(const RegionalGlobalInputs& inputs, const RegionalWindow& window,
                                  const FineWarp& warp, const RegionalParameters& parameters, bool lakegate)
{
    const int width = inputs.width, height = inputs.height;
    const int k = window.scale;
    const int domainwidth = window.domainfinewidth(), domainheight = window.domainfineheight();
    const std::size_t count = static_cast<std::size_t>(domainwidth) * static_cast<std::size_t>(domainheight);
    const bool coastnoise = parameters.coast.noiseamplitude_m > 0.0;
    const bool landshare = parameters.coast.landshare > 0.0;

    StartSurface result;
    if (warp.empty() && !coastnoise && !landshare)
    {
        // Phase E (erosion plan E3): the base elevation upsampled at the fine centres, plus
        // noise of the band from twice the global north-south spacing (exclusive) down to
        // twice the fine one, with an amplitude that grows with the local global relief and
        // tapers to zero within coastband_m of the sea level; then the block means.
        const double sealevel = static_cast<double>(inputs.sealevel);
        result.base_m = upsample(inputs.elevation_m, *inputs.grid, *regionalwindowgrid(window, inputs.radius_m));
        result.lapsebase_m = result.base_m;
        const std::vector<float>& base = result.base_m;
        const std::vector<float> relief = upsample(noiserelief(inputs, parameters), *inputs.grid, *regionalwindowgrid(window, inputs.radius_m));
        FractalNoise noise;
        noise.seed = deterministiccombine(inputs.seed, noisesalt);
        noise.radius_km = inputs.radius_m / 1000.0;
        noise.longest_km = 2.0 * pi * inputs.radius_m / height / 1000.0;
        noise.wavelength_km = noise.longest_km;
        noise.shortest_km = noise.longest_km / k;
        noise.octaves = 1;
        while ((1 << noise.octaves) <= k)
            noise.octaves++;
        noise.persistence = parameters.noisepersistence;
        const FractalNoiseField field(noise);

        std::vector<float>& surface = result.start_m;
        surface = base;
        if (field.resolved() > 0)
        {
            for (int y = 0; y < domainheight; y++)
            {
                const double latitude = window.finelatitude_rad(y);
                for (int x = 0; x < domainwidth; x++)
                {
                    const std::size_t cell = at(x, y, domainwidth);
                    double px = 0.0, py = 0.0, pz = 0.0;
                    unitvector(latitude, window.finelongitude_rad(x), px, py, pz);
                    const double amplitude =
                        parameters.noiseamplitude_m +
                        parameters.noiserelieffraction * std::max(0.0, static_cast<double>(relief[cell]));
                    const double height_m = std::fabs(static_cast<double>(base[cell]) - sealevel);
                    // The taper scales with the local amplitude, so a strong noise cannot
                    // lift a seabed or sink a lowland across the sea level far from the coast.
                    const double band = parameters.coastband_m + amplitude;
                    const double taper = band > 0.0 ? std::clamp(height_m / band, 0.0, 1.0) : 1.0;
                    surface[cell] = static_cast<float>(static_cast<double>(base[cell]) +
                                                       taper * amplitude * field.sum(px, py, pz) / field.normalisation());
                }
            }
        }
        correctblockmeans(inputs.elevation_m, width, height, window, parameters.meanpasses, surface);
    }
    else
    {
        // Phase R (R2): the same per cell (Surfacer, which reproduces the expressions above),
        // with the faded warp and the coast noise; with the land-share rule the mask is
        // decided first and kept through the block-mean passes.
        const Surfacer surfacer(inputs, k, parameters);
        if (!landshare)
        {
            result.base_m.resize(count);
            result.lapsebase_m.resize(count);
            result.start_m.resize(count);
            for (int y = 0; y < domainheight; y++)
                for (int x = 0; x < domainwidth; x++)
                {
                    const std::size_t cell = at(x, y, domainwidth);
                    const CellSample sample =
                        surfacer.sample(window.absolutefinecolumn(x), window.absolutefinerow(y),
                                        warp.empty() ? 0.0 : warp.dx[cell], warp.empty() ? 0.0 : warp.dy[cell]);
                    result.base_m[cell] = sample.base;
                    result.lapsebase_m[cell] = sample.lapsebase;
                    result.start_m[cell] = sample.surface;
                }
            correctblockmeans(inputs.elevation_m, width, height, window, parameters.meanpasses, result.start_m);
        }
        else
        {
            // The warp is evaluated again on the widened box (regionalwarp's values, a pure
            // function of position, so `warp` is not read).
            landsharesurface(inputs, window, parameters, surfacer, result, lakegate);
        }
    }

    {
        const std::vector<float>& surface = result.start_m;
        const int x0 = window.keptfinecolumnoffset(), y0 = window.keptfinerowoffset();
        double largest = 0.0;
        for (int r = 0; r < window.rows; r++)
            for (int c = 0; c < window.columns; c++)
            {
                double sum = 0.0;
                for (int y = y0 + r * k; y < y0 + (r + 1) * k; y++)
                    for (int x = x0 + c * k; x < x0 + (c + 1) * k; x++)
                        sum += static_cast<double>(surface[at(x, y, domainwidth)]);
                const float global = inputs.elevation_m[at(wrapcolumn(window.column0 + c, width), window.row0 + r, width)];
                const double error = std::fabs(sum / (static_cast<double>(k) * k) - static_cast<double>(global));
                if (error > largest)
                {
                    largest = error;
                    result.maxblockerrorcolumn = wrapcolumn(window.column0 + c, width);
                    result.maxblockerrorrow = window.row0 + r;
                }
            }
        result.blockmeanerror_m = largest;
    }
    return result;
}

void correctsignedblockmeans(const RegionalGlobalInputs& inputs, const mesh::CellGrid& window, std::vector<float>& fine)
{
    if (inputs.grid == nullptr)
        throw std::invalid_argument("correctsignedblockmeans: the inputs need their level-0 grid");
    const mesh::CellGrid& levelzero = *inputs.grid;
    if (levelzero.kind() != mesh::GridKind::cube || window.kind() != mesh::GridKind::cube || levelzero.facecells() < 1 ||
        window.facecells() % levelzero.facecells() != 0 || fine.size() != window.cellcount() ||
        inputs.elevation_m.size() != levelzero.cellcount())
        throw std::invalid_argument("correctsignedblockmeans: a cube level 0 and a cube window refining it");
    const WindowBlocks blocks = windowblocks(window, window.facecells() / levelzero.facecells());
    std::vector<double> residual(blocks.count());
    for (std::size_t b = 0; b < blocks.count(); b++)
        residual[b] = static_cast<double>(static_cast<float>(
            static_cast<double>(inputs.elevation_m[levelzero.localofkey(blocks.parents[b])]) - blockmean(blocks, b, fine)));
    const std::vector<double> correction = gatherblocks(blocks, residual, window);
    const double sealevel = static_cast<double>(inputs.sealevel);
    for (std::size_t cell = 0; cell < fine.size(); cell++)
        fine[cell] = signedmove(fine[cell], correction[cell], sealevel);
}

double levelnoiseshare(int scale, int parentscale, double persistence)
{
    if (parentscale < 1 || scale < parentscale)
        throw std::invalid_argument("levelnoiseshare: the scale must be at least its parent's, which is at least 1");
    // The octaves o = 0 .. floor(log2 scale) of addlevelnoise's set, weight persistence^o; in the band
    // when parentscale < 2^o <= scale (the wavelength longest / 2^o in [longest / scale, longest / parentscale)).
    double all = 0.0, band = 0.0;
    for (int o = 0; (std::int64_t{1} << o) <= scale; o++)
    {
        const double weight = std::pow(persistence, o);
        all += weight * weight;
        if ((std::int64_t{1} << o) > parentscale)
            band += weight * weight;
    }
    return all > 0.0 ? std::sqrt(band / all) : 0.0;
}

void addlevelnoise(const RegionalGlobalInputs& inputs, const RegionalWindow& window, int parentscale,
                   const RegionalParameters& parameters, std::vector<float>& surface_m)
{
    const int width = inputs.width, height = inputs.height;
    const int k = window.scale;
    const int domainwidth = window.domainfinewidth(), domainheight = window.domainfineheight();
    if (parentscale < 1 || k <= parentscale || k % parentscale != 0)
        throw std::invalid_argument("addlevelnoise: the level's scale must be a larger multiple of its parent's");
    if (surface_m.size() != static_cast<std::size_t>(domainwidth) * static_cast<std::size_t>(domainheight))
        throw std::invalid_argument("addlevelnoise: the surface must cover the window's domain");
    const double sealevel = static_cast<double>(inputs.sealevel);
    // Phase E's noise (the first branch of regionalstartsurface), its band cut at the parent's scale.
    FractalNoise noise;
    noise.seed = deterministiccombine(inputs.seed, noisesalt);
    noise.radius_km = inputs.radius_m / 1000.0;
    const double longest = 2.0 * pi * inputs.radius_m / height / 1000.0;
    noise.wavelength_km = longest;
    noise.longest_km = longest / parentscale;
    noise.shortest_km = longest / k;
    noise.octaves = 1;
    while ((1 << noise.octaves) <= k)
        noise.octaves++;
    noise.persistence = parameters.noisepersistence;
    const FractalNoiseField field(noise);
    if (field.resolved() == 0)
        return;
    const std::vector<float> relief = upsample(noiserelief(inputs, parameters), *inputs.grid, *regionalwindowgrid(window, inputs.radius_m));
    for (int y = 0; y < domainheight; y++)
    {
        const double latitude = window.finelatitude_rad(y);
        for (int x = 0; x < domainwidth; x++)
        {
            const std::size_t cell = at(x, y, domainwidth);
            double px = 0.0, py = 0.0, pz = 0.0;
            unitvector(latitude, window.finelongitude_rad(x), px, py, pz);
            const double amplitude =
                parameters.noiseamplitude_m + parameters.noiserelieffraction * std::max(0.0, static_cast<double>(relief[cell]));
            const double base = static_cast<double>(surface_m[cell]);
            const double band = parameters.coastband_m + amplitude;
            const double taper = band > 0.0 ? std::clamp(std::fabs(base - sealevel) / band, 0.0, 1.0) : 1.0;
            surface_m[cell] = static_cast<float>(base + taper * amplitude * field.sum(px, py, pz) / field.normalisation());
        }
    }
}

std::vector<uint8_t> regionalsea(const RegionalGlobalInputs& inputs, const RegionalWindow& window,
                                 const FineWarp& warp, const RegionalParameters& parameters,
                                 std::vector<float>& surface_m, const std::vector<uint8_t>* parentgate, bool lakegate)
{
    const int k = window.scale;
    const int domainwidth = window.domainfinewidth(), domainheight = window.domainfineheight();
    const std::size_t count = static_cast<std::size_t>(domainwidth) * static_cast<std::size_t>(domainheight);
    const bool centrerule = !(parameters.coast.landshare > 0.0);
    if (parentgate != nullptr && parentgate->size() != count)
        throw std::invalid_argument("regionalsea: the parent gate must hold one value per domain cell");

    // At or below the sea level, with a global sea cell among the global cell and its
    // 8 neighbours. Then (Phase E, centre rule) each global cell's centre fine cell
    // follows the global mask: a sea centre above the sea level is lowered to the sea
    // level, a land centre counted as sea is raised to the next float above it. Without
    // `lakegate` a fine cell of a global lake cell (lake surface above zero) is never sea unless
    // its global cell is sea: a lake below the sea level beside the ocean stays a lake. With it
    // (S3.4, the loop) a lake cell's fringe cells, where the warped lake weight is below 0.5, may
    // be sea when the gate's rule (a global sea cell in the 3 x 3, or the smooth gate's reach)
    // allows it; cells where the weight is at least 0.5 never are. R8: the global rule is the
    // gate's (regional_edges.hpp, seagate; the 3 x 3 rule at edges.smoothgate 0), so the sea
    // agrees with the start surface's mask.
    // S3.4: in the refinement loop the gate's lake exclusion follows the lake cells' warped half-mask
    // (regional_edges.hpp, seagate's lakeweight), so a lake beside the sea meets it off the block edge.
    const std::vector<float> lakeweight = lakegate ? lakecellweight(inputs, window, warp) : std::vector<float>{};
    const std::vector<uint8_t> gate = seagate(inputs, window, warp, parameters, lakegate ? &lakeweight : nullptr);
    const float abovesea = std::nextafter(inputs.sealevel, std::numeric_limits<float>::infinity());
    std::vector<uint8_t> sea(count, 0);
    for (int y = 0; y < domainheight; y++)
    {
        const bool centrerow = window.absolutefinerow(y) % k == k / 2;
        for (int x = 0; x < domainwidth; x++)
        {
            const std::size_t cell = at(x, y, domainwidth);
            const std::size_t global = window.globalcell(x, y);
            if (centrerule && centrerow && window.absolutefinecolumn(x) % k == k / 2)
            {
                if (inputs.network.sea[global] != 0)
                    surface_m[cell] = std::min(surface_m[cell], inputs.sealevel);
                else if (surface_m[cell] <= inputs.sealevel && gate[cell] != 0)
                    surface_m[cell] = abovesea;
            }
            // S2.8: at level >= 2 the sea is inherited from the parent level within one parent
            // cell; every other cell is land at its own height (the surface is not written).
            const bool open = gate[cell] != 0 && (parentgate == nullptr || (*parentgate)[cell] != 0);
            sea[cell] = surface_m[cell] <= inputs.sealevel && open ? 1 : 0;
        }
    }
    return sea;
}

void noiselevel(const LevelContext& context, LevelState& state)
{
    // Multiscale plan S3.0: the start surface of regional_pass.cpp's step 1 (level 1, the
    // single pass) and the loop's octave band (level >= 2), moved here unchanged.
    if (context.parent == nullptr)
    {
        state.start = regionalstartsurface(context.inputs, context.window, context.warp, context.parameters,
                                           context.levelwindow != nullptr);
        // S3.15: level 1's pockets are filled by the land-share rule (coastlevel fills the deeper levels').
        if (state.diagnostics != nullptr)
        {
            state.diagnostics->pockets = state.start.pockets;
            state.diagnostics->pocketcells = state.start.pocketcells;
            state.diagnostics->pocketsfilled = state.start.pocketsfilled;
            state.diagnostics->pocketcellsfilled = state.start.pocketcellsfilled;
            state.diagnostics->maxblockerrorcolumn = state.start.maxblockerrorcolumn;
            state.diagnostics->maxblockerrorrow = state.start.maxblockerrorrow;
        }
    }
    else
        addlevelnoise(context.inputs, context.window, static_cast<int>(levelscale(context.factor, context.level - 1)),
                      context.parameters, state.start.start_m);
    // The level's noise band per cell (S3.4): the amplitude of the noise above (before the taper,
    // at the unwarped cell centres) times the share of the level's octaves in it.
    const RegionalGlobalInputs& inputs = context.inputs;
    const RegionalParameters& parameters = context.parameters;
    const int parentscale = context.parent == nullptr ? 1 : static_cast<int>(levelscale(context.factor, context.level - 1));
    const double share = levelnoiseshare(context.window.scale, parentscale, parameters.noisepersistence);
    state.noiseband_m = upsample(noiserelief(inputs, parameters), *inputs.grid, context.grid);  // G on the window grid
    for (float& band : state.noiseband_m)
        band = static_cast<float>(share * (parameters.noiseamplitude_m +
                                           parameters.noiserelieffraction * std::max(0.0, static_cast<double>(band))));
    recordlevelstep(state, LevelStep::noiseoctave);
}

namespace
{
// S3.15 at level >= 2 (regional_surface.hpp, coastlevel): the inherited sea with the level's band of
// the shoreline displacement, the parent-share guard and the pocket fill.
void refinedcoast(const LevelContext& context, LevelState& state)
{
    const RegionalGlobalInputs& inputs = context.inputs;
    const RegionalWindow& window = context.window;
    const RegionalParameters& parameters = context.parameters;
    const LevelWindow& level = *context.levelwindow;
    const ParentFields& parent = *context.parent;
    const int k = window.scale, factor = context.factor, width = inputs.width;
    const int domainwidth = window.domainfinewidth(), domainheight = window.domainfineheight();
    const std::size_t count = static_cast<std::size_t>(domainwidth) * static_cast<std::size_t>(domainheight);
    const double sealevel = static_cast<double>(inputs.sealevel);
    const float abovesea = std::nextafter(inputs.sealevel, std::numeric_limits<float>::infinity());
    std::vector<float>& surface = state.surface_m;

    // The cells the sea may hold: regionalsea's gates (its result on a surface below the sea level
    // everywhere; with the land-share rule it writes no surface), and its sea on the level's surface.
    std::vector<float> below(count, std::nextafter(inputs.sealevel, -std::numeric_limits<float>::infinity()));
    const std::vector<uint8_t> open =
        regionalsea(inputs, window, context.warp, parameters, below, context.details.parentseagate, true);
    std::vector<uint8_t> sea(count, 0);
    for (std::size_t cell = 0; cell < count; cell++)
        sea[cell] = surface[cell] <= inputs.sealevel && open[cell] != 0 ? 1 : 0;

    // The parent's own sea per parent cell of the domain (the domain is whole parent cells), and its land
    // at the projection's floor (at most four float steps above the sea level: the level-1 clamp and the
    // projection leave land there), whose children cannot all be land (they may not sink to the sea level,
    // so their mean cannot come down to the parent's: the synthetic coast windows of the tests forced one).
    const int parentcolumns = domainwidth / factor, parentrows = domainheight / factor;
    const std::size_t parents = static_cast<std::size_t>(parentcolumns) * static_cast<std::size_t>(parentrows);
    std::vector<uint8_t> parentsea(parents), parentfloor(parents);
    float floorlevel = inputs.sealevel;
    for (int step = 0; step < 4; step++)
        floorlevel = std::nextafter(floorlevel, std::numeric_limits<float>::infinity());
    for (int py = 0; py < parentrows; py++)
        for (int px = 0; px < parentcolumns; px++)
        {
            const std::int64_t pc = level.absolutecolumn(px * factor) / factor, pr = level.absoluterow(py * factor) / factor;
            const std::size_t index = parent.window.index(parent.window.localcolumn(pc), static_cast<int>(parent.window.localrow(pr)));
            parentsea[at(px, py, parentcolumns)] = parent.sea[index];
            parentfloor[at(px, py, parentcolumns)] = parent.sea[index] == 0 && parent.elevation_m[index] <= floorlevel ? 1 : 0;
        }

    // (a) The level's band of the displacement: a cell is sea where the surface at its displaced
    // position lies at or below the sea level (and the gates allow it). Per parent cell the result is
    // kept only while the parent's share rule holds for its children (a parent sea cell at most
    // 1 - landshare land, a land cell at least landshare) and a land parent at the floor that held sea
    // keeps some, else the undisplaced sea stays. A cell that changes side takes the displaced surface
    // (a sea cell at most the sea level, a land cell at least the next float above it; the least move
    // instead left window A's level-3 border step at 13 m and raised its land below the sea level).
    const int children = factor * factor;
    const double share = parameters.coast.landshare;
    const int leastland = static_cast<int>(std::ceil(share * children - 1.0e-9));
    const int mostland = static_cast<int>(std::floor((1.0 - share) * children + 1.0e-9));
    const std::vector<uint8_t> undisplaced = sea;
    const std::vector<float> original = surface;
    const ShorelineWarp warp(width, inputs.height, k, k / factor, inputs.seed, parameters.coast);
    if (warp.active())
    {
        const std::vector<double> displaced = displacedsample(surface, window, warp);
        std::vector<uint8_t> moved(static_cast<std::size_t>(children));
        for (int py = 0; py < parentrows; py++)
            for (int px = 0; px < parentcolumns; px++)
            {
                const bool seaparent = parentsea[at(px, py, parentcolumns)] != 0;
                int land = 0, seabefore = 0;
                for (int y = py * factor; y < (py + 1) * factor; y++)
                    for (int x = px * factor; x < (px + 1) * factor; x++)
                    {
                        const std::size_t cell = at(x, y, domainwidth);
                        const uint8_t value = open[cell] != 0 && displaced[cell] <= sealevel ? 1 : 0;
                        moved[at(x - px * factor, y - py * factor, factor)] = value;
                        land += value == 0 ? 1 : 0;
                        seabefore += sea[cell] != 0 ? 1 : 0;
                    }
                if (seaparent ? land > mostland
                              : land < leastland || (parentfloor[at(px, py, parentcolumns)] != 0 && seabefore > 0 && land == children))
                    continue;
                for (int y = py * factor; y < (py + 1) * factor; y++)
                    for (int x = px * factor; x < (px + 1) * factor; x++)
                    {
                        const std::size_t cell = at(x, y, domainwidth);
                        const uint8_t value = moved[at(x - px * factor, y - py * factor, factor)];
                        if (value == sea[cell])
                            continue;
                        surface[cell] = value != 0 ? std::min(inputs.sealevel, static_cast<float>(displaced[cell]))
                                                   : std::max(abovesea, static_cast<float>(displaced[cell]));
                        sea[cell] = value;
                    }
            }
    }

    // (b) The sea pockets (seapockets, anchored on the core sea: the undisplaced sea of the global sea
    // cells, the level's decision without the displacement; the interpolated global elevation as the core
    // left 103 enclosed components on window A at k 64 against 47) filled per parent cell under the same
    // rule: a parent sea cell's pocket cells only while its
    // children stay at most 1 - landshare land, a land parent's at the floor only while some sea is left
    // among its children (the test's forced parent was a land parent at the next float above the sea
    // level whose one sea child was a pocket). The filled pockets' basins (pocketbasins: the land at or
    // below the sea level around them within the reach) are raised to the sea level + coast.shorefloor_m.
    if (parameters.coast.pocketreach > 0)
    {
        const int reach = parameters.coast.pocketreach;
        std::vector<uint8_t> core(count, 0);
        for (int y = 0; y < domainheight; y++)
            for (int x = 0; x < domainwidth; x++)
            {
                const std::size_t cell = at(x, y, domainwidth);
                core[cell] = undisplaced[cell] != 0 && inputs.network.sea[window.globalcell(x, y)] != 0 ? 1 : 0;
            }
        std::vector<uint8_t> filled = seapockets(inputs, window, reach, sea, core);
        // A pocket the share rule keeps in a parent sea cell is most often sea the displacement cut off:
        // around such parents (and their 8 neighbours) the displacement is undone and the pockets found again.
        if (warp.active())
        {
            std::vector<uint8_t> undo(parents, 0);
            bool any = false;
            for (int py = 0; py < parentrows; py++)
                for (int px = 0; px < parentcolumns; px++)
                {
                    if (parentsea[at(px, py, parentcolumns)] == 0)
                        continue;
                    int land = 0, fill = 0;
                    for (int y = py * factor; y < (py + 1) * factor; y++)
                        for (int x = px * factor; x < (px + 1) * factor; x++)
                        {
                            land += sea[at(x, y, domainwidth)] == 0 ? 1 : 0;
                            fill += filled[at(x, y, domainwidth)] != 0 ? 1 : 0;
                        }
                    if (fill == 0 || land + fill <= mostland)
                        continue;
                    any = true;
                    for (int dy = -1; dy <= 1; dy++)
                        for (int dx = -1; dx <= 1; dx++)
                            if (py + dy >= 0 && py + dy < parentrows && px + dx >= 0 && px + dx < parentcolumns)
                                undo[at(px + dx, py + dy, parentcolumns)] = 1;
                }
            if (any)
            {
                for (int y = 0; y < domainheight; y++)
                    for (int x = 0; x < domainwidth; x++)
                        if (undo[at(x / factor, y / factor, parentcolumns)] != 0)
                        {
                            const std::size_t cell = at(x, y, domainwidth);
                            sea[cell] = undisplaced[cell];
                            surface[cell] = original[cell];
                        }
                filled = seapockets(inputs, window, reach, sea, core);
            }
        }
        for (int py = 0; py < parentrows; py++)
            for (int px = 0; px < parentcolumns; px++)
            {
                int land = 0, fill = 0;
                for (int y = py * factor; y < (py + 1) * factor; y++)
                    for (int x = px * factor; x < (px + 1) * factor; x++)
                    {
                        land += sea[at(x, y, domainwidth)] == 0 ? 1 : 0;
                        fill += filled[at(x, y, domainwidth)] != 0 ? 1 : 0;
                    }
                const bool keep = fill > 0 && (parentsea[at(px, py, parentcolumns)] != 0 ? land + fill > mostland
                                                                                          : parentfloor[at(px, py, parentcolumns)] != 0 && land + fill == children);
                for (int y = py * factor; y < (py + 1) * factor; y++)
                    for (int x = px * factor; x < (px + 1) * factor; x++)
                    {
                        const std::size_t cell = at(x, y, domainwidth);
                        if (keep)
                            filled[cell] = 0;
                        else if (filled[cell] != 0)
                            sea[cell] = 0;
                    }
            }
        const float raised = std::max(abovesea, static_cast<float>(sealevel + parameters.coast.shorefloor_m));
        const std::vector<uint8_t> basin = pocketbasins(inputs, window, reach, sea, filled, surface, abovesea);
        // In every parent (a land parent at the floor whose children end all land above it is the
        // projection's to force, counted by its forcedparents: a filled bed left hundreds of metres below
        // the sea level as dry land made a 411 m border step on window A at level 2).
        for (std::size_t cell = 0; cell < count; cell++)
            if (basin[cell] != 0)
                surface[cell] = std::max(surface[cell], raised);
        if (state.diagnostics != nullptr)
        {
            const int x0 = window.keptfinecolumnoffset(), y0 = window.keptfinerowoffset();
            const PocketCount done = countpockets(filled, domainwidth, x0, y0, window.finewidth(), window.fineheight());
            const PocketCount left = countpockets(seapockets(inputs, window, reach, sea, core), domainwidth, x0, y0,
                                                  window.finewidth(), window.fineheight());
            state.diagnostics->pocketsfilled = done.pockets;
            state.diagnostics->pocketcellsfilled = done.pocketcells;
            state.diagnostics->pockets = left.pockets;
            state.diagnostics->pocketcells = left.pocketcells;
        }
    }
    state.sea = std::move(sea);
}
}

void coastlevel(const LevelContext& context, LevelState& state)
{
    // 2. The sea by local rules on the start surface, before the channels (R1), so the
    // channels end at the fine coast and never change the mask (regional_pass.cpp's step 2,
    // moved here unchanged in S3.0).
    const CoastParameters& coast = context.parameters.coast;
    if (context.parent != nullptr && context.levelwindow != nullptr && coast.landshare > 0.0 &&
        (coast.warpcells > 0.0 || coast.pocketreach > 0))
        refinedcoast(context, state);
    else
        state.sea = regionalsea(context.inputs, context.window, context.warp, context.parameters, state.surface_m,
                                context.details.parentseagate, context.levelwindow != nullptr);
    recordlevelstep(state, LevelStep::coast);
}
}
