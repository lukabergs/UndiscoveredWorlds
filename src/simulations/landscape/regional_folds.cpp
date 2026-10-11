#include "simulations/landscape/regional_folds.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <unordered_map>

#include "core/deterministic_random.hpp"
#include "simulations/landscape/refinement_grid.hpp"
#include "simulations/landscape/refinement_levels.hpp"
#include "simulations/landscape/regional_pass.hpp"
#include "simulations/landscape/regional_valleys.hpp"
#include "simulations/landscape/upsampling.hpp"

// Multiscale plan S5.2 and S5.3b (regional_folds.hpp): the fold fabric, phasor-normalised stripes
// from one splat per cell of a jittered absolute lattice (Phacelle's layout and compact weight),
// across the exported shortening direction interpolated at each splat.

namespace landscape
{
namespace
{
constexpr double pi = 3.14159265358979323846;
// Salt of the splats' draws (deterministiccombine(seed, salt, lattice column, lattice row)).
constexpr std::int64_t foldsalt = 0x464f4c44;  // "FOLD"
// The compact weight's cutoff in lattice cells (Phacelle: the splats of the cells within two).
constexpr double cutoffcells = 2.0;

std::size_t at(int x, int y, int width)
{
    return static_cast<std::size_t>(y) * static_cast<std::size_t>(width) + static_cast<std::size_t>(x);
}

double smoothstep(double value, double lower, double upper)
{
    if (!(upper > lower))
        return value >= lower ? 1.0 : 0.0;
    const double t = std::clamp((value - lower) / (upper - lower), 0.0, 1.0);
    return t * t * (3.0 - 2.0 * t);
}

double strainweight(const FoldParameters& folds, double strain)
{
    return smoothstep(strain, folds.strainmin, 2.0 * folds.strainmin);
}

// A splat: its centre in absolute global coordinates (columns, rows; global cell (i, j) spans
// [i, i + 1) x [j, j + 1)), the phase origin; its wavelength and its across-strike unit axis (east,
// south), interpolated at the centre.
struct Splat
{
    double u = 0.0, v = 0.0;
    double wavelength_km = 0.0;
    double east = 1.0, south = 0.0;
};

// The global fields a splat interpolates: the strain-weighted and the plain doubled-angle axis and the
// wavelength, one value per global cell.
struct SplatFields
{
    std::vector<float> weightedcos, weightedsin, plaincos, plainsin, wavelength_km;
};

// The splats of the absolute lattice, drawn once per lattice cell and kept for the call.
class SplatLattice
{
public:
    SplatLattice(const RegionalGlobalInputs& inputs, const FoldParameters& folds, std::int64_t subdivision,
                 const SplatFields& fields)
        : inputs_(inputs), folds_(folds), subdivision_(subdivision), fields_(fields)
    {
        columns_ = static_cast<std::int64_t>(inputs.width) * subdivision;
        rows_ = static_cast<std::int64_t>(inputs.height) * subdivision;
    }

    std::int64_t columns() const { return columns_; }
    std::int64_t rows() const { return rows_; }

    // The splat of lattice cell (column, row), column wrapped to 0 .. columns() - 1, row in 0 .. rows() - 1.
    const Splat& splat(std::int64_t column, std::int64_t row)
    {
        const std::int64_t key = row * columns_ + column;
        const auto found = splats_.find(key);
        if (found != splats_.end())
            return found->second;
        Splat& splat = splats_[key];
        const std::uint64_t hash = deterministiccombine(inputs_.seed, foldsalt, column, row);
        const double fx = static_cast<double>(hash & 0xffffffffull) / 4294967296.0;
        const double fy = static_cast<double>(hash >> 32) / 4294967296.0;
        const double scale = static_cast<double>(subdivision_);
        splat.u = (static_cast<double>(column) + fx) / scale;
        splat.v = (static_cast<double>(row) + fy) / scale;
        const int width = inputs_.width, height = inputs_.height;
        // The axis on doubled angles, interpolated at the splat: the strain-weighted field where it does
        // not vanish there, else the plain one (east without either).
        double c = bicubic(fields_.weightedcos, width, height, splat.u, splat.v);
        double s = bicubic(fields_.weightedsin, width, height, splat.u, splat.v);
        if (!(std::hypot(c, s) > 1.0e-9))
        {
            c = bicubic(fields_.plaincos, width, height, splat.u, splat.v);
            s = bicubic(fields_.plainsin, width, height, splat.u, splat.v);
        }
        const double theta = std::hypot(c, s) > 0.0 ? 0.5 * std::atan2(s, c) : 0.0;
        splat.east = std::cos(theta);
        splat.south = std::sin(theta);
        const double lower = std::min(folds_.wavelengthmin_km, folds_.wavelengthmax_km);
        const double upper = std::max(folds_.wavelengthmin_km, folds_.wavelengthmax_km);
        const double wavelength = bicubic(fields_.wavelength_km, width, height, splat.u, splat.v);
        splat.wavelength_km = std::clamp(std::isfinite(wavelength) ? wavelength : lower, lower, upper);
        return splat;
    }

private:
    const RegionalGlobalInputs& inputs_;
    const FoldParameters& folds_;
    std::int64_t subdivision_ = 1;
    std::int64_t columns_ = 0, rows_ = 0;
    const SplatFields& fields_;
    std::unordered_map<std::int64_t, Splat> splats_;
};

// The nearest-rank percentile of `values` (sorted in place; 0 when empty).
double nearestrank(std::vector<double>& values, double fraction)
{
    if (values.empty())
        return 0.0;
    const std::size_t rank = static_cast<std::size_t>(std::ceil(fraction * static_cast<double>(values.size())));
    const std::size_t index = std::min(values.size() - 1, rank > 0 ? rank - 1 : 0);
    std::nth_element(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(index), values.end());
    return values[index];
}
}

bool foldfabricon(const RegionalGlobalInputs& inputs, const RegionalParameters& parameters)
{
    const std::size_t size = static_cast<std::size_t>(inputs.width) * static_cast<std::size_t>(inputs.height);
    return parameters.folds.amplitudemax_m > 0.0 && parameters.folds.splatspacing_km > 0.0 && size > 0 &&
           inputs.foldstrain.size() == size;
}

double foldwavelength_km(const FoldParameters& folds, double sediment_m, bool withsediment)
{
    const double lower = std::min(folds.wavelengthmin_km, folds.wavelengthmax_km);
    const double upper = std::max(folds.wavelengthmin_km, folds.wavelengthmax_km);
    const double cover_km = withsediment ? std::max(0.0, sediment_m) / 1000.0 : 0.0;
    const double raw = folds.wavelengthbase_km + folds.wavelengthpercover_km * cover_km;
    return std::clamp(std::isfinite(raw) ? raw : 0.0, lower, upper);
}

double foldamplitude_m(const FoldParameters& folds, double strain, double wavelength_km, double heightabovesea_m)
{
    if (!(folds.amplitudemax_m > 0.0) || !(strain > 0.0) || !(wavelength_km > 0.0))
        return 0.0;
    const double structural = folds.relieffraction * wavelength_km * 1000.0 / pi * std::sqrt(strain);
    const double elevation = folds.reliefelevationfraction * std::max(0.0, heightabovesea_m);
    const double amplitude = std::min(std::min(structural, folds.amplitudemax_m), elevation);
    return std::isfinite(amplitude) ? std::max(0.0, amplitude) * strainweight(folds, strain) : 0.0;
}

void addfoldfabric(const RegionalGlobalInputs& inputs, const RegionalWindow& window, const RegionalParameters& parameters,
                   const FineChannels& channels, const std::vector<uint8_t>& sea, std::vector<float>& elevation_m,
                   int parentscale, std::vector<uint8_t>* floorclamped)
{
    if (!foldfabricon(inputs, parameters))
        return;
    const FoldParameters& folds = parameters.folds;
    const int width = window.domainfinewidth(), height = window.domainfineheight();
    const std::size_t count = static_cast<std::size_t>(width) * static_cast<std::size_t>(height);
    if (count == 0 || channels.channel.size() != count || sea.size() != count || elevation_m.size() != count)
        return;
    if (floorclamped != nullptr)
        floorclamped->assign(count, 0);
    const int globalwidth = inputs.width, globalheight = inputs.height;
    const std::size_t globalsize = static_cast<std::size_t>(globalwidth) * static_cast<std::size_t>(globalheight);
    const int k = window.scale;
    const PatchMetric metric = finepatchmetric(window, inputs.radius_m);
    const double cell_km = metric.dy_m / 1000.0;
    // The global north-south cell on the ground and the octave band of the level (addlevelnoise's).
    const double globalcell_km = pi * inputs.radius_m / static_cast<double>(globalheight) / 1000.0;
    if (!(cell_km > 0.0) || !(globalcell_km > 0.0))
        return;
    const double longest_km = 2.0 * globalcell_km;
    const double bandlower = longest_km / static_cast<double>(k);
    const double bandupper = parentscale > 0 ? longest_km / static_cast<double>(parentscale)
                                             : std::numeric_limits<double>::infinity();

    // The global fields: the wavelength and the doubled-angle axis per global cell, strain-weighted and plain.
    const bool withsediment = inputs.sediment_m.size() == globalsize;
    const bool withdirection = inputs.folddirection_rad.size() == globalsize;
    SplatFields fields;
    fields.weightedcos.assign(globalsize, 0.0f);
    fields.weightedsin.assign(globalsize, 0.0f);
    fields.plaincos.assign(globalsize, 1.0f);
    fields.plainsin.assign(globalsize, 0.0f);
    fields.wavelength_km.resize(globalsize);
    std::vector<float> heightfield(globalsize);
    const double sealevel = static_cast<double>(inputs.sealevel);
    for (std::size_t cell = 0; cell < globalsize; cell++)
    {
        fields.wavelength_km[cell] = static_cast<float>(
            foldwavelength_km(folds, withsediment ? static_cast<double>(inputs.sediment_m[cell]) : 0.0, withsediment));
        heightfield[cell] = static_cast<float>(std::max(0.0, static_cast<double>(inputs.elevation_m[cell]) - sealevel));
        const double angle = withdirection ? 2.0 * static_cast<double>(inputs.folddirection_rad[cell]) : 0.0;
        const double weight = strainweight(folds, static_cast<double>(inputs.foldstrain[cell]));
        const double c = std::cos(angle), s = std::sin(angle);
        fields.plaincos[cell] = static_cast<float>(c);
        fields.plainsin[cell] = static_cast<float>(s);
        fields.weightedcos[cell] = static_cast<float>(weight * c);
        fields.weightedsin[cell] = static_cast<float>(weight * s);
    }

    // The lattice: the global grid subdivided by the power of two whose cell is nearest the spacing.
    const double ratiolog = std::log2(globalcell_km / folds.splatspacing_km);
    const int exponent = std::clamp(static_cast<int>(std::lround(std::isfinite(ratiolog) ? ratiolog : 0.0)), 0, 20);
    const std::int64_t subdivision = std::int64_t{1} << exponent;
    const double lattice_km = globalcell_km / static_cast<double>(subdivision);
    const double cutoff_km = cutoffcells * lattice_km;
    SplatLattice lattice(inputs, folds, subdivision, fields);
    const double worldcolumns = static_cast<double>(globalwidth);
    const double scale = static_cast<double>(subdivision);

    // The cell fields, interpolated as the noise amplitude's relief is.
    const std::vector<float> strain = upsample(inputs.foldstrain, globalwidth, globalheight, window);
    const std::vector<float> cellwavelength = upsample(fields.wavelength_km, globalwidth, globalheight, window);
    const std::vector<float> above = upsample(heightfield, globalwidth, globalheight, window);
    const bool withice = inputs.glaciation.size() == globalsize;
    const std::vector<float> glaciation = withice ? upsample(inputs.glaciation, globalwidth, globalheight, window)
                                                  : std::vector<float>{};
    const bool withlakes = inputs.lakesurface_m.size() == globalsize && inputs.network.sea.size() == globalsize;
    const double wavelengthlower = std::min(folds.wavelengthmin_km, folds.wavelengthmax_km);
    const double wavelengthupper = std::max(folds.wavelengthmin_km, folds.wavelengthmax_km);
    const float floorlevel = std::max(std::nextafter(inputs.sealevel, std::numeric_limits<float>::infinity()),
                                      static_cast<float>(sealevel + std::max(0.0, parameters.coast.shorefloor_m)));

    // Every cell reads the elevation before the fabric only at itself, so the order does not matter.
    for (int y = 0; y < height; y++)
    {
        const double ratio = eastwestratio(metric, y);
        const int globalrow = window.globalrow(y);
        const double v = (static_cast<double>(window.absolutefinerow(y)) + 0.5) / static_cast<double>(k);
        const std::int64_t latticerow = static_cast<std::int64_t>(std::floor(v * scale));
        // The lattice columns within the cutoff on the ground: 2 / ratio lattice cells either side.
        const double columnreach = ratio > 0.0 ? cutoffcells / ratio : std::numeric_limits<double>::infinity();
        for (int x = 0; x < width; x++)
        {
            const std::size_t cell = at(x, y, width);
            if (sea[cell] != 0 || channels.channel[cell] != 0)
                continue;
            // Dry land at or below the shore floor (a filled sea pocket, a basin below the sea level)
            // takes no fabric: it could only rise, and the projection must then sink the parent's
            // other children (the first judging run forced parents there).
            if (!(elevation_m[cell] > floorlevel))
                continue;
            // The gates (S5.3b): under the global ice cover, and in a global lake cell's block.
            if (withice && std::clamp(static_cast<double>(glaciation[cell]), 0.0, 1.0) > folds.icemax)
                continue;
            if (withlakes)
            {
                const std::size_t global = at(window.globalcolumn(x), globalrow, globalwidth);
                if (inputs.lakesurface_m[global] > 0.0f && inputs.network.sea[global] == 0)
                    continue;
            }
            const double epsilon = std::max(0.0, static_cast<double>(strain[cell]));
            if (!(strainweight(folds, epsilon) > 0.0))
                continue;
            const double wavelengthcell = std::clamp(static_cast<double>(cellwavelength[cell]), wavelengthlower, wavelengthupper);
            if (!(wavelengthcell >= bandlower) || !(wavelengthcell < bandupper))
                continue;  // drawn at another level of the refinement
            const double amplitude = foldamplitude_m(folds, epsilon, wavelengthcell, static_cast<double>(above[cell]));
            if (!(amplitude > 0.0))
                continue;

            // The water gap: the largest falloff over the channel cells within gapfraction x L.
            double gap = 1.0;
            const double gapcells = folds.gapfraction * wavelengthcell / cell_km;
            if (gapcells > 0.0)
            {
                const int rows = std::min(k, static_cast<int>(std::floor(gapcells)));
                const int columns = columnextent(gapcells, ratio, k);
                double strongest = 0.0;
                for (int dy = -rows; dy <= rows; dy++)
                {
                    const int ny = y + dy;
                    if (ny < 0 || ny >= height)
                        continue;
                    for (int dx = -columns; dx <= columns; dx++)
                    {
                        const int nx = x + dx;
                        if (nx < 0 || nx >= width || channels.channel[at(nx, ny, width)] == 0)
                            continue;
                        const double east = static_cast<double>(dx) * ratio;
                        const double distance = std::sqrt(east * east + static_cast<double>(dy * dy));
                        if (distance >= gapcells)
                            continue;
                        strongest = std::max(strongest, 1.0 - smoothstep(distance / gapcells, 0.0, 1.0));
                    }
                }
                gap = 1.0 - strongest;
            }
            if (!(gap > 0.0))
                continue;

            // The phasor sum over the splats within the cutoff of the cell's centre.
            const double u = (static_cast<double>(window.absolutefinecolumn(x)) + 0.5) / static_cast<double>(k);
            const double latticeu = u * scale;
            std::int64_t firstcolumn = 0, lastcolumn = lattice.columns() - 1;
            if (2.0 * columnreach + 2.0 < static_cast<double>(lattice.columns()))
            {
                firstcolumn = static_cast<std::int64_t>(std::floor(latticeu - columnreach));
                lastcolumn = static_cast<std::int64_t>(std::floor(latticeu + columnreach));
            }
            double real = 0.0, imaginary = 0.0, weights = 0.0;
            for (std::int64_t row = latticerow - 2; row <= latticerow + 2; row++)
            {
                if (row < 0 || row >= lattice.rows())
                    continue;
                for (std::int64_t column = firstcolumn; column <= lastcolumn; column++)
                {
                    const Splat& splat = lattice.splat(wraplevel(column, lattice.columns()), row);
                    if (!(splat.wavelength_km > 0.0))
                        continue;
                    double du = u - splat.u;
                    du -= worldcolumns * std::round(du / worldcolumns);
                    const double east = du * globalcell_km * ratio;
                    const double south = (v - splat.v) * globalcell_km;
                    const double distance2 = east * east + south * south;
                    if (!(distance2 < cutoff_km * cutoff_km))
                        continue;
                    const double width_km = folds.kernelwidth_km > 0.0 ? folds.kernelwidth_km : splat.wavelength_km;
                    const double sigma2 = width_km * width_km;
                    const double weight = std::exp(-distance2 / sigma2) - std::exp(-cutoff_km * cutoff_km / sigma2);
                    if (!(weight > 0.0))
                        continue;
                    const double phase = 2.0 * pi * (splat.east * east + splat.south * south) / splat.wavelength_km;
                    real += weight * std::cos(phase);
                    imaginary += weight * std::sin(phase);
                    weights += weight;
                }
            }
            double unit = 0.0;
            if (folds.phasor != 0)
            {
                const double magnitude = std::hypot(real, imaginary);
                unit = magnitude > 0.0 ? imaginary / magnitude : 0.0;  // sin(arg z)
            }
            else
                unit = weights > 0.0 ? real / weights : 0.0;
            const double fabric = amplitude * gap * unit;
            if (fabric == 0.0 || !std::isfinite(fabric))
                continue;
            const float before = elevation_m[cell];
            const float value = static_cast<float>(static_cast<double>(before) + fabric);
            // Never lowered below the shore floor above the sea level.
            if (value < floorlevel && floorclamped != nullptr)
                (*floorclamped)[cell] = 1;
            elevation_m[cell] = std::max(value, floorlevel);
        }
    }
}

void foldfabriclevel(const LevelContext& context, LevelState& state)
{
    if (!foldfabricon(context.inputs, context.parameters))
        return;
    const RegionalWindow& window = context.window;
    const int parentscale = context.parent != nullptr ? static_cast<int>(levelscale(context.factor, context.level - 1)) : 0;
    const double cap = channelcellcap_m(context);
    const bool record = state.diagnostics != nullptr;
    std::vector<float> before = state.elevation_m;
    std::vector<uint8_t> floorclamped;
    addfoldfabric(context.inputs, window, context.parameters, state.channels, state.sea, state.elevation_m, parentscale,
                  record ? &floorclamped : nullptr);
    if (record)
    {
        const int domainwidth = window.domainfinewidth();
        const int x0 = window.keptfinecolumnoffset(), y0 = window.keptfinerowoffset();
        std::vector<double> changes;
        for (int y = y0; y < y0 + window.fineheight(); y++)
            for (int x = x0; x < x0 + window.finewidth(); x++)
            {
                const std::size_t cell = at(x, y, domainwidth);
                if (!floorclamped.empty() && floorclamped[cell] != 0)
                    state.diagnostics->foldfloorclamped++;
                if (before[cell] != state.elevation_m[cell])
                    changes.push_back(std::fabs(static_cast<double>(state.elevation_m[cell]) - static_cast<double>(before[cell])));
            }
        state.diagnostics->foldcells = static_cast<int>(changes.size());
        state.diagnostics->foldamplitude_p50_m = nearestrank(changes, 0.5);
        state.diagnostics->foldamplitude_p99_m = nearestrank(changes, 0.99);
    }
    channelcellmeans(window, state.channels, state.sea, cap, before, state.elevation_m);
    // S5.3b: the surfaces before and after the fabric, for the lakes step (unfoldforlakes).
    state.folded_m = state.elevation_m;
    state.unfolded_m = std::move(before);
    recordlevelstep(state, LevelStep::foldfabric);
}

void unfoldforlakes(const LevelContext& context, LevelState& state)
{
    (void)context;
    std::vector<float>& elevation = state.elevation_m;
    const std::size_t count = elevation.size();
    if (state.folded_m.size() != count || state.unfolded_m.size() != count || state.sea.size() != count ||
        state.channels.channel.size() != count)
    {
        state.folded_m.clear();
        state.unfolded_m.clear();
        return;
    }
    for (std::size_t cell = 0; cell < count; cell++)
    {
        const float folded = state.folded_m[cell], prior = state.unfolded_m[cell];
        if (folded == prior)
            continue;
        if (state.sea[cell] != 0 || state.channels.channel[cell] != 0)
        {
            state.unfolded_m[cell] = folded;  // a channel cell's mean: left, and not refolded
            continue;
        }
        const float present = elevation[cell];
        const float unfolded = present == folded ? prior : std::min(present, prior);
        state.folded_m[cell] = present;
        state.unfolded_m[cell] = unfolded;
        elevation[cell] = unfolded;
    }
}

void refoldafterlakes(const LevelContext& context, LevelState& state)
{
    std::vector<float>& elevation = state.elevation_m;
    const std::size_t count = elevation.size();
    const int width = context.window.domainfinewidth(), height = context.window.domainfineheight();
    if (state.folded_m.size() == count && state.unfolded_m.size() == count &&
        count == static_cast<std::size_t>(width) * static_cast<std::size_t>(height))
    {
        const std::vector<float>& surface = state.water.lakesurface_m;
        const bool withwater = surface.size() == count;
        for (int y = 0; y < height; y++)
            for (int x = 0; x < width; x++)
            {
                const std::size_t cell = at(x, y, width);
                const float folded = state.folded_m[cell], unfolded = state.unfolded_m[cell];
                if (folded == unfolded || elevation[cell] != unfolded)
                    continue;  // not folded, or changed by the lakes (a lake bed, a shore raise, a rim)
                float level = 0.0f;
                if (withwater)
                {
                    if (surface[cell] > 0.0f)
                        continue;  // under the water
                    for (int dy = -1; dy <= 1; dy++)
                        for (int dx = -1; dx <= 1; dx++)
                        {
                            const int nx = x + dx, ny = y + dy;
                            if ((dx != 0 || dy != 0) && nx >= 0 && nx < width && ny >= 0 && ny < height)
                                level = std::max(level, surface[at(nx, ny, width)]);
                        }
                }
                if (level > 0.0f && !(unfolded > level))
                    continue;  // already at or below the water beside it: no fabric on the shore
                elevation[cell] =
                    level > 0.0f && !(folded > level) ? std::nextafter(level, std::numeric_limits<float>::infinity()) : folded;
            }
    }
    state.folded_m.clear();
    state.folded_m.shrink_to_fit();
    state.unfolded_m.clear();
    state.unfolded_m.shrink_to_fit();
}
}
