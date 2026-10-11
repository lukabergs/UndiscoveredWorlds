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

// Multiscale plan S5.2 (regional_folds.hpp): the fold fabric, phasor-normalised Gabor stripes from
// kernels hashed per absolute tile, across the exported shortening direction.

namespace landscape
{
namespace
{
constexpr double pi = 3.14159265358979323846;
// Salt of the fold kernels' draws (deterministiccombine(seed, salt, tile column, tile row, n)).
constexpr std::int64_t foldsalt = 0x464f4c44;  // "FOLD"

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

// A Gabor kernel: its centre in absolute global coordinates (columns, rows; global cell (i, j)
// spans [i, i + 1) x [j, j + 1)), its wavelength and its across-strike unit axis (east, south).
struct Kernel
{
    double u = 0.0, v = 0.0;
    double originu = 0.0, originv = 0.0;  // the phase origin: the centre of the global cell that owns it
    double wavelength_km = 0.0;
    double east = 1.0, south = 0.0;
};

// The kernels of the absolute lattice, drawn once per tile and kept for the call.
class KernelLattice
{
public:
    KernelLattice(const RegionalGlobalInputs& inputs, const FoldParameters& folds, std::int64_t subdivision,
                  const std::vector<double>& wavelength_km, const std::vector<double>& axiseast,
                  const std::vector<double>& axissouth)
        : inputs_(inputs), folds_(folds), subdivision_(subdivision), wavelength_(wavelength_km), axiseast_(axiseast),
          axissouth_(axissouth)
    {
        columns_ = static_cast<std::int64_t>(inputs.width) * subdivision;
        rows_ = static_cast<std::int64_t>(inputs.height) * subdivision;
    }

    std::int64_t columns() const { return columns_; }
    std::int64_t rows() const { return rows_; }

    // The kernels of tile (column, row), column wrapped to 0 .. columns() - 1, row in 0 .. rows() - 1.
    const std::vector<Kernel>& tile(std::int64_t column, std::int64_t row)
    {
        const std::int64_t key = row * columns_ + column;
        const auto found = tiles_.find(key);
        if (found != tiles_.end())
            return found->second;
        std::vector<Kernel>& kernels = tiles_[key];
        const int count = std::max(0, folds_.kernelspertile);
        kernels.reserve(static_cast<std::size_t>(count));
        const double scale = static_cast<double>(subdivision_);
        for (int n = 0; n < count; n++)
        {
            const std::uint64_t hash = deterministiccombine(inputs_.seed, foldsalt, column, row, n);
            const double fx = static_cast<double>(hash & 0xffffffffull) / 4294967296.0;
            const double fy = static_cast<double>(hash >> 32) / 4294967296.0;
            Kernel kernel;
            kernel.u = (static_cast<double>(column) + fx) / scale;
            kernel.v = (static_cast<double>(row) + fy) / scale;
            const int i = std::clamp(static_cast<int>(std::floor(kernel.u)), 0, inputs_.width - 1);
            const int j = std::clamp(static_cast<int>(std::floor(kernel.v)), 0, inputs_.height - 1);
            // The kernels of one global cell share its d, L and phase origin, so they add coherently and
            // interfere only where d or L change, at the cell boundaries (S5.2's review).
            kernel.originu = static_cast<double>(i) + 0.5;
            kernel.originv = static_cast<double>(j) + 0.5;
            const std::size_t global = at(i, j, inputs_.width);
            kernel.wavelength_km = wavelength_[global];
            kernel.east = axiseast_[global];
            kernel.south = axissouth_[global];
            kernels.push_back(kernel);
        }
        return kernels;
    }

private:
    const RegionalGlobalInputs& inputs_;
    const FoldParameters& folds_;
    std::int64_t subdivision_ = 1;
    std::int64_t columns_ = 0, rows_ = 0;
    const std::vector<double>& wavelength_;
    const std::vector<double>& axiseast_;
    const std::vector<double>& axissouth_;
    std::unordered_map<std::int64_t, std::vector<Kernel>> tiles_;
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
    return parameters.folds.amplitudemax_m > 0.0 && parameters.folds.kernelspertile >= 1 && size > 0 &&
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

    // The global fields: the wavelength and the strain-weighted doubled-angle axis per global cell.
    const bool withsediment = inputs.sediment_m.size() == globalsize;
    const bool withdirection = inputs.folddirection_rad.size() == globalsize;
    std::vector<double> wavelength(globalsize), axiseast(globalsize, 1.0), axissouth(globalsize, 0.0);
    std::vector<float> wavelengthfield(globalsize), heightfield(globalsize);
    const double sealevel = static_cast<double>(inputs.sealevel);
    for (std::size_t cell = 0; cell < globalsize; cell++)
    {
        wavelength[cell] = foldwavelength_km(folds, withsediment ? static_cast<double>(inputs.sediment_m[cell]) : 0.0,
                                             withsediment);
        wavelengthfield[cell] = static_cast<float>(wavelength[cell]);
        heightfield[cell] = static_cast<float>(std::max(0.0, static_cast<double>(inputs.elevation_m[cell]) - sealevel));
    }
    if (withdirection)
    {
        for (int j = 0; j < globalheight; j++)
            for (int i = 0; i < globalwidth; i++)
            {
                double weighted[2] = {0.0, 0.0}, plain[2] = {0.0, 0.0}, weights = 0.0;
                for (int dj = -1; dj <= 1; dj++)
                {
                    const int nj = j + dj;
                    if (nj < 0 || nj >= globalheight)
                        continue;
                    for (int di = -1; di <= 1; di++)
                    {
                        const std::size_t other = at((i + di + globalwidth) % globalwidth, nj, globalwidth);
                        const double angle = 2.0 * static_cast<double>(inputs.folddirection_rad[other]);
                        const double weight = strainweight(folds, static_cast<double>(inputs.foldstrain[other]));
                        const double c = std::cos(angle), s = std::sin(angle);
                        weighted[0] += weight * c;
                        weighted[1] += weight * s;
                        plain[0] += c;
                        plain[1] += s;
                        weights += weight;
                    }
                }
                const double* sum = weights > 0.0 ? weighted : plain;
                const double theta = sum[0] == 0.0 && sum[1] == 0.0
                                         ? static_cast<double>(inputs.folddirection_rad[at(i, j, globalwidth)])
                                         : 0.5 * std::atan2(sum[1], sum[0]);
                axiseast[at(i, j, globalwidth)] = std::cos(theta);
                axissouth[at(i, j, globalwidth)] = std::sin(theta);
            }
    }

    // The lattice: the largest subdivision whose tile is at least r = 2 max(L_max, the kernel width) on the
    // ground (the envelope's width is the kernel width, or the kernel's own L at most L_max with width 0).
    const double kernelwidth = folds.kernelwidth_km > 0.0 ? folds.kernelwidth_km : 0.0;
    const double reach_km = 2.0 * std::max(kernelwidth, std::max(folds.wavelengthmin_km, folds.wavelengthmax_km));
    if (!(reach_km > 0.0))
        return;
    std::int64_t subdivision = 1;
    while (subdivision < (std::int64_t{1} << 20) && globalcell_km / static_cast<double>(2 * subdivision) >= reach_km)
        subdivision *= 2;
    const double tile_km = globalcell_km / static_cast<double>(subdivision);
    KernelLattice lattice(inputs, folds, subdivision, wavelength, axiseast, axissouth);
    const std::int64_t rowreach = std::max<std::int64_t>(1, static_cast<std::int64_t>(std::ceil(reach_km / tile_km)));
    const std::int64_t columnreachmax = std::max<std::int64_t>(0, (lattice.columns() - 1) / 2);

    // The cell fields, interpolated as the noise amplitude's relief is.
    const std::vector<float> strain = upsample(inputs.foldstrain, globalwidth, globalheight, window);
    const std::vector<float> cellwavelength = upsample(wavelengthfield, globalwidth, globalheight, window);
    const std::vector<float> above = upsample(heightfield, globalwidth, globalheight, window);
    const double wavelengthlower = std::min(folds.wavelengthmin_km, folds.wavelengthmax_km);
    const double wavelengthupper = std::max(folds.wavelengthmin_km, folds.wavelengthmax_km);
    const float floorlevel = std::max(std::nextafter(inputs.sealevel, std::numeric_limits<float>::infinity()),
                                      static_cast<float>(sealevel + std::max(0.0, parameters.coast.shorefloor_m)));
    const double worldcolumns = static_cast<double>(globalwidth);

    // Every cell reads the elevation before the fabric only at itself, so the order does not matter.
    for (int y = 0; y < height; y++)
    {
        const double ratio = eastwestratio(metric, y);
        const double v = (static_cast<double>(window.absolutefinerow(y)) + 0.5) / static_cast<double>(k);
        const std::int64_t tilerow = static_cast<std::int64_t>(std::floor(v * static_cast<double>(subdivision)));
        const std::int64_t columnreach =
            ratio > 0.0 ? std::min(columnreachmax, std::max<std::int64_t>(1, static_cast<std::int64_t>(
                                                                                  std::ceil(reach_km / (tile_km * ratio)))))
                        : columnreachmax;
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

            // The phasor sum over the kernels of the tiles within r of the cell's centre.
            const double u = (static_cast<double>(window.absolutefinecolumn(x)) + 0.5) / static_cast<double>(k);
            const std::int64_t tilecolumn = static_cast<std::int64_t>(std::floor(u * static_cast<double>(subdivision)));
            double real = 0.0, imaginary = 0.0, weights = 0.0;
            for (std::int64_t row = tilerow - rowreach; row <= tilerow + rowreach; row++)
            {
                if (row < 0 || row >= lattice.rows())
                    continue;
                for (std::int64_t column = tilecolumn - columnreach; column <= tilecolumn + columnreach; column++)
                {
                    for (const Kernel& kernel : lattice.tile(wraplevel(column, lattice.columns()), row))
                    {
                        if (!(kernel.wavelength_km > 0.0))
                            continue;
                        double du = u - kernel.u;
                        du -= worldcolumns * std::round(du / worldcolumns);
                        const double east = du * globalcell_km * ratio;
                        const double south = (v - kernel.v) * globalcell_km;
                        const double distance2 = east * east + south * south;
                        const double envelope = kernelwidth > 0.0 ? kernelwidth : kernel.wavelength_km;
                        const double weight = std::exp(-distance2 / (2.0 * envelope * envelope));
                        double originu = u - kernel.originu;
                        originu -= worldcolumns * std::round(originu / worldcolumns);
                        const double origineast = originu * globalcell_km * ratio;
                        const double originsouth = (v - kernel.originv) * globalcell_km;
                        const double phase =
                            2.0 * pi * (kernel.east * origineast + kernel.south * originsouth) / kernel.wavelength_km;
                        real += weight * std::cos(phase);
                        imaginary += weight * std::sin(phase);
                        weights += weight;
                    }
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
    std::vector<float> before;
    if (record || cap > 0.0)
        before = state.elevation_m;
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
    recordlevelstep(state, LevelStep::foldfabric);
}
}
