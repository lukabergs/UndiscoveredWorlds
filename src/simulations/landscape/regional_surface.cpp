#include "simulations/landscape/regional_surface.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>

#include "core/deterministic_random.hpp"
#include "simulations/landscape/regional_coast.hpp"
#include "simulations/landscape/regional_edges.hpp"

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
          coastal_(coastalcells(inputs.network.sea, inputs.width, inputs.height)),
          inland_(inlandcells(inputs)),
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

        // Phase E's noise, tapered at the sea level (regional_pass.cpp's expression of E3).
        float banded = cell.base;
        if (band_->resolved() > 0)
        {
            double px = 0.0, py = 0.0, pz = 0.0;
            unitvector(finelatitude(finerow, height * scale_), finelongitude(finecolumn, width * scale_), px, py, pz);
            const double amplitude =
                parameters_.noiseamplitude_m + parameters_.noiserelieffraction * std::max(0.0, static_cast<double>(relief));
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
                    0.0, parameters_.noiseamplitude_m + parameters_.noiserelieffraction * std::max(0.0, static_cast<double>(relief)));
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

void correctsignedblockmeans(const RegionalGlobalInputs& inputs, const RegionalWindow& window, std::vector<float>& fine)
{
    const int width = inputs.width, height = inputs.height;
    const int k = window.scale;
    const int finewidth = window.domainfinewidth(), fineheight = window.domainfineheight();
    const int dcolumns = window.domaincolumns, drows = window.domainrows;
    const double sealevel = static_cast<double>(inputs.sealevel);
    std::vector<double> residual(static_cast<std::size_t>(dcolumns) * static_cast<std::size_t>(drows));
    for (int dr = 0; dr < drows; dr++)
        for (int dc = 0; dc < dcolumns; dc++)
        {
            // Summed in row-major order inside the block, as correctblockmeans sums.
            double sum = 0.0;
            for (int y = dr * k; y < (dr + 1) * k; y++)
                for (int x = dc * k; x < (dc + 1) * k; x++)
                    sum += static_cast<double>(fine[at(x, y, finewidth)]);
            const std::size_t global = at(wrapcolumn(window.domaincolumn0 + dc, width), window.domainrow0 + dr, width);
            residual[at(dc, dr, dcolumns)] =
                static_cast<double>(inputs.elevation_m[global]) - sum / (static_cast<double>(k) * k);
        }

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
            const double h = static_cast<double>(value) - sealevel;
            const double moved = h > 0.0 ? std::max(h + correction, (1.0 - signedlimit) * h)
                                         : std::min(h + correction, (1.0 - signedlimit) * h);
            value = static_cast<float>(sealevel + moved);
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
                      const Surfacer& surfacer, StartSurface& start)
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
    const int margin = sharereach * sharecount + 2 + 2 * signedcount;
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
    const std::vector<uint8_t> gate = seagate(inputs, box, warp, parameters);
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
        }

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
            if (mask[cell] != 0)
                shoredepth[cell] = static_cast<float>(fade * std::max(0.0, std::min(-value, floor)));
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
                                  const FineWarp& warp, const RegionalParameters& parameters)
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
        result.base_m = upsample(inputs.elevation_m, width, height, window);
        result.lapsebase_m = result.base_m;
        const std::vector<float>& base = result.base_m;
        const std::vector<float> relief = upsample(localrelief(inputs.elevation_m, width, height), width, height, window);
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
            landsharesurface(inputs, window, parameters, surfacer, result);
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
                largest = std::max(largest, std::fabs(sum / (static_cast<double>(k) * k) - static_cast<double>(global)));
            }
        result.blockmeanerror_m = largest;
    }
    return result;
}

std::vector<uint8_t> regionalsea(const RegionalGlobalInputs& inputs, const RegionalWindow& window,
                                 const FineWarp& warp, const RegionalParameters& parameters,
                                 std::vector<float>& surface_m)
{
    (void)warp;
    const int width = inputs.width;
    const int k = window.scale;
    const int domainwidth = window.domainfinewidth(), domainheight = window.domainfineheight();
    const std::size_t count = static_cast<std::size_t>(domainwidth) * static_cast<std::size_t>(domainheight);
    const bool centrerule = !(parameters.coast.landshare > 0.0);

    // At or below the sea level, with a global sea cell among the global cell and its
    // 8 neighbours. Then (Phase E, centre rule) each global cell's centre fine cell
    // follows the global mask: a sea centre above the sea level is lowered to the sea
    // level, a land centre counted as sea is raised to the next float above it. A fine
    // cell of a global lake cell (lake surface above zero) is never sea unless its global
    // cell is sea: a lake below the sea level beside the ocean stays a lake.
    const std::vector<uint8_t> nearsea = nearseacells(inputs);
    const float abovesea = std::nextafter(inputs.sealevel, std::numeric_limits<float>::infinity());
    std::vector<uint8_t> sea(count, 0);
    for (int y = 0; y < domainheight; y++)
    {
        const int row = window.globalrow(y);
        const bool centrerow = window.absolutefinerow(y) % k == k / 2;
        for (int x = 0; x < domainwidth; x++)
        {
            const std::size_t cell = at(x, y, domainwidth);
            const std::size_t global = at(window.globalcolumn(x), row, width);
            if (centrerule && centrerow && window.absolutefinecolumn(x) % k == k / 2)
            {
                if (inputs.network.sea[global] != 0)
                    surface_m[cell] = std::min(surface_m[cell], inputs.sealevel);
                else if (surface_m[cell] <= inputs.sealevel && nearsea[global] != 0)
                    surface_m[cell] = abovesea;
            }
            sea[cell] = surface_m[cell] <= inputs.sealevel && nearsea[global] != 0 ? 1 : 0;
        }
    }
    return sea;
}
}
