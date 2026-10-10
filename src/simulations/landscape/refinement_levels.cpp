#include "simulations/landscape/refinement_levels.hpp"

#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>

#include "simulations/landscape/regional_surface.hpp"
#include "simulations/landscape/regional_valleys.hpp"
#include "simulations/landscape/upsampling.hpp"

// The cell size and the level assignment are step S0.4; the prolongation, the parent-mean
// correction and the loop are step S2.3 (docs/architecture/multiscale-plan.md, section 3,
// "S2.3's level placement").

namespace landscape
{
namespace
{
constexpr double pi = 3.14159265358979323846;

std::size_t at(int x, int y, int width)
{
    return static_cast<std::size_t>(y) * static_cast<std::size_t>(width) + static_cast<std::size_t>(x);
}

// The interpolation position of a child at offset o (0 .. factor - 1) inside its parent:
// t = (o + 0.5) / factor - 0.5 parent cells from the parent's centre; the stencil's second
// point is the parent itself (lower 0) or its predecessor (lower -1, t < 0), at fraction
// t or t + 1 between it and the next. A function of the offset only, so every window
// interpolates a child alike.
struct ChildOffsets
{
    std::vector<int> lower;
    std::vector<double> fraction;
};

ChildOffsets childoffsets(int factor)
{
    ChildOffsets offsets;
    offsets.lower.resize(static_cast<std::size_t>(factor));
    offsets.fraction.resize(static_cast<std::size_t>(factor));
    for (int o = 0; o < factor; o++)
    {
        const double t = (o + 0.5) / factor - 0.5;
        offsets.lower[static_cast<std::size_t>(o)] = t < 0.0 ? -1 : 0;
        offsets.fraction[static_cast<std::size_t>(o)] = t < 0.0 ? t + 1.0 : t;
    }
    return offsets;
}

void checkpair(const LevelWindow& parent, const LevelWindow& level, const char* name)
{
    if (level.level < 1 || level.factor < 2 || parent.level != level.level - 1 || parent.factor != level.factor ||
        parent.worldwidth != level.worldwidth || parent.worldheight != level.worldheight)
        throw std::invalid_argument(std::string(name) + ": the windows are not a level and its parent");
}

// The parent values over the parents of a level window's domain widened by `ring` parent
// cells on every side: (domainparentcolumns() + 2 ring) x (domainparentrows() + 2 ring),
// row-major; patch row r holds the parent row domainrow0 / factor - ring + r clamped to
// the parent level's grid, patch column c the parent column domaincolumn0 / factor - ring + c
// wrapped. Every value read lies in the parent window's kept box, or this throws.
struct ParentPatch
{
    int columns = 0;
    int rows = 0;
    std::vector<float> values;
};

ParentPatch parentpatch(const LevelWindow& parent, const std::vector<float>& parentfield, const LevelWindow& level,
                        int ring, const char* name)
{
    checkpair(parent, level, name);
    if (parentfield.size() != parent.domainsize())
        throw std::invalid_argument(std::string(name) + ": the parent field must hold the parent window's domain");
    const int factor = level.factor;
    const std::int64_t parentwidth = parent.levelwidth(), parentheight = parent.levelheight();
    const std::int64_t column0 = level.domaincolumn0 / factor - ring;
    const std::int64_t row0 = level.domainrow0 / factor - ring;
    ParentPatch patch;
    patch.columns = level.domainparentcolumns() + 2 * ring;
    patch.rows = level.domainparentrows() + 2 * ring;
    patch.values.resize(static_cast<std::size_t>(patch.columns) * static_cast<std::size_t>(patch.rows));
    const int keptx = parent.keptcolumnoffset(), kepty = parent.keptrowoffset();
    for (int r = 0; r < patch.rows; r++)
    {
        const std::int64_t row = std::clamp<std::int64_t>(row0 + r, 0, parentheight - 1);
        const std::int64_t localrow = parent.localrow(row);
        if (localrow < kepty || localrow >= kepty + parent.rows)
            throw std::invalid_argument(std::string(name) + ": a parent row read lies outside the parent's kept box");
        for (int c = 0; c < patch.columns; c++)
        {
            const int localcolumn = parent.localcolumn(wraplevel(column0 + c, parentwidth));
            if (localcolumn < keptx || localcolumn >= keptx + parent.columns)
                throw std::invalid_argument(std::string(name) + ": a parent column read lies outside the parent's kept box");
            patch.values[at(c, r, patch.columns)] = parentfield[parent.index(localcolumn, static_cast<int>(localrow))];
        }
    }
    return patch;
}

// The level's projected field averaged over every global cell of its kept box, against
// the global elevation: the largest absolute difference.
double globalmeanerror(const LevelWindow& window, const std::vector<float>& field, const RegionalGlobalInputs& inputs)
{
    const std::int64_t scale = window.scale();
    const int s = static_cast<int>(scale);
    const int x0 = window.keptcolumnoffset(), y0 = window.keptrowoffset();
    const double cells = static_cast<double>(s) * static_cast<double>(s);
    double largest = 0.0;
    for (int r = 0; r < window.rows / s; r++)
        for (int c = 0; c < window.columns / s; c++)
        {
            double sum = 0.0;
            for (int y = y0 + r * s; y < y0 + (r + 1) * s; y++)
                for (int x = x0 + c * s; x < x0 + (c + 1) * s; x++)
                    sum += static_cast<double>(field[window.index(x, y)]);
            const int column = static_cast<int>(wraplevel(window.column0 / scale + c, inputs.width));
            const int row = static_cast<int>(window.row0 / scale + r);
            const double global = static_cast<double>(inputs.elevation_m[at(column, row, inputs.width)]);
            largest = std::max(largest, std::fabs(sum / cells - global));
        }
    return largest;
}
}

double levelcellsize_km(double radius_m, int worldwidth, int factor, int level)
{
    if (!std::isfinite(radius_m) || !(radius_m > 0.0))
        throw std::invalid_argument("levelcellsize_km: the radius must be positive and finite");
    if (worldwidth < 4 || worldwidth % 2 != 0)
        throw std::invalid_argument("levelcellsize_km: the world width must be even and at least 4");
    if (factor < 2 || level < 0)
        throw std::invalid_argument("levelcellsize_km: the factor must be at least 2 and the level not negative");
    const double rows = static_cast<double>(worldwidth / 2) * static_cast<double>(levelscale(factor, level));
    return pi * radius_m / rows / 1000.0;
}

int detaillevel(double size_km, double detailcells, double radius_m, int worldwidth, int factor, int levels)
{
    if (!std::isfinite(size_km) || !(size_km > 0.0))
        throw std::invalid_argument("detaillevel: the detail's size must be positive and finite");
    if (!std::isfinite(detailcells) || !(detailcells > 0.0))
        throw std::invalid_argument("detaillevel: the cells per detail must be positive and finite");
    if (levels < 1)
        throw std::invalid_argument("detaillevel: at least one level must run");
    const double cell_km = size_km / detailcells;
    for (int level = 1; level <= levels; level++)
        if (levelcellsize_km(radius_m, worldwidth, factor, level) <= cell_km)
            return level;
    return levels;
}

std::vector<float> prolong(const LevelWindow& parent, const std::vector<float>& parentfield, const LevelWindow& level)
{
    const ParentPatch patch = parentpatch(parent, parentfield, level, prolongationreach, "prolong");
    const int factor = level.factor;
    const ChildOffsets offsets = childoffsets(factor);
    std::vector<float> fine(level.domainsize());
    for (int y = 0; y < level.domainrows; y++)
    {
        const std::size_t oy = static_cast<std::size_t>(y % factor);
        const int rowbase = y / factor + prolongationreach + offsets.lower[oy];
        const double ty = offsets.fraction[oy];
        for (int x = 0; x < level.domaincolumns; x++)
        {
            const std::size_t ox = static_cast<std::size_t>(x % factor);
            const int columnbase = x / factor + prolongationreach + offsets.lower[ox];
            const double tx = offsets.fraction[ox];
            double rows[4];
            for (int b = 0; b < 4; b++)
            {
                const std::size_t base = static_cast<std::size_t>(rowbase - 1 + b) * static_cast<std::size_t>(patch.columns);
                const auto value = [&](int a) {
                    return static_cast<double>(patch.values[base + static_cast<std::size_t>(columnbase - 1 + a)]);
                };
                rows[b] = catmullrom(value(0), value(1), value(2), value(3), tx);
            }
            fine[level.index(x, y)] = static_cast<float>(catmullrom(rows[0], rows[1], rows[2], rows[3], ty));
        }
    }
    return fine;
}

void correctparentmeans(const LevelWindow& parent, const std::vector<float>& parentfield, const LevelWindow& level,
                        int passes, std::vector<float>& fine)
{
    const ParentPatch patch = parentpatch(parent, parentfield, level, 0, "correctparentmeans");
    if (passes < 0 || fine.size() != level.domainsize())
        throw std::invalid_argument("correctparentmeans: passes must not be negative and fine must cover the domain");
    const int factor = level.factor;
    const int columns = patch.columns, rows = patch.rows;
    const bool wholeturn = static_cast<std::int64_t>(level.domaincolumns) >= level.levelwidth();
    const ChildOffsets offsets = childoffsets(factor);
    const double children = static_cast<double>(factor) * factor;
    std::vector<float> residual(static_cast<std::size_t>(columns) * static_cast<std::size_t>(rows));
    const auto column = [&](int c) { return wholeturn ? ((c % columns) + columns) % columns : std::clamp(c, 0, columns - 1); };
    for (int pass = 0; pass < passes; pass++)
    {
        // The residual of every domain parent: the mean of its children, summed in row-major
        // order (the absolute order inside a parent), minus the parent value.
        for (int pr = 0; pr < rows; pr++)
            for (int pc = 0; pc < columns; pc++)
            {
                double sum = 0.0;
                for (int y = pr * factor; y < (pr + 1) * factor; y++)
                    for (int x = pc * factor; x < (pc + 1) * factor; x++)
                        sum += static_cast<double>(fine[level.index(x, y)]);
                residual[at(pc, pr, columns)] =
                    static_cast<float>(sum / children - static_cast<double>(patch.values[at(pc, pr, columns)]));
            }
        // Interpolated at the child centres (the two parents beyond the domain's edge take
        // the nearest domain parent's residual) and subtracted.
        for (int y = 0; y < level.domainrows; y++)
        {
            const std::size_t oy = static_cast<std::size_t>(y % factor);
            const int rowbase = y / factor + offsets.lower[oy];
            const double ty = offsets.fraction[oy];
            for (int x = 0; x < level.domaincolumns; x++)
            {
                const std::size_t ox = static_cast<std::size_t>(x % factor);
                const int columnbase = x / factor + offsets.lower[ox];
                const double tx = offsets.fraction[ox];
                double interpolated[4];
                for (int b = 0; b < 4; b++)
                {
                    const std::size_t base = static_cast<std::size_t>(std::clamp(rowbase - 1 + b, 0, rows - 1)) *
                                             static_cast<std::size_t>(columns);
                    const auto value = [&](int a) {
                        return static_cast<double>(residual[base + static_cast<std::size_t>(column(columnbase - 1 + a))]);
                    };
                    interpolated[b] = catmullrom(value(0), value(1), value(2), value(3), tx);
                }
                float& cell = fine[level.index(x, y)];
                cell = static_cast<float>(static_cast<double>(cell) -
                                          catmullrom(interpolated[0], interpolated[1], interpolated[2], interpolated[3], ty));
            }
        }
    }
}

RefinementResult runrefinement(const RegionalGlobalInputs& inputs, const LevelTarget& target,
                               const RegionalParameters& parameters, const RefinementParameters& refinement)
{
    const int factor = refinement.factor, levels = target.level;
    if (levels < 1 || factor < 2)
        throw std::invalid_argument("runrefinement: the level must be at least 1 and the factor at least 2");
    {
        RegionalWindow check = target.window;
        check.scale = 1;  // the target's scale is not read
        checkregionalinputs(inputs, check);
    }
    if (refinement.smoothpasses < 0 || !std::isfinite(refinement.smoothrelaxation) || !(refinement.smoothrelaxation > 0.0))
        throw std::invalid_argument("runrefinement: refinement.smooth_passes must not be negative and "
                                    "refinement.smooth_relaxation must be positive and finite");
    if (parameters.meanpasses < 0 || target.window.halo < refinementhalo(parameters.meanpasses, refinement.smoothpasses))
        throw std::invalid_argument("runrefinement: the halo must be at least 2 x meanpasses + 2 + 2 x smooth_passes + "
                                    "2 (the lake flood) global cells");
    std::int64_t finest = 1;
    for (int level = 0; level < levels; level++)
    {
        finest *= factor;
        if (finest > INT_MAX)
            throw std::invalid_argument("runrefinement: factor^level must fit an int");
    }

    const std::vector<LevelWindow> chain = halochain(target.window, factor, levels, prolongationreach);
    const GlobalIncision incision = globalincision(inputs, parameters);
    // Every detail unit runs at every level when its size key is above 0; each form is
    // carved where its size in level cells lies in the band (RegionalDetails).
    RegionalDetails details;
    details.valleys = parameters.valleys.floor_km > 0.0;
    details.troughs = parameters.troughs.width_km > 0.0;
    details.cirques = parameters.cirques.radius_km > 0.0;
    details.factor = factor;
    details.detailcells = refinement.detailcells;
    details.globallakedepth = true;
    if (!std::isfinite(refinement.detailcells) || !(refinement.detailcells > 0.0))
        throw std::invalid_argument("runrefinement: refinement.detail_cells must be positive and finite");

    // Level 0: the global elevation on the chain's level-0 window.
    std::vector<float> parentfield(chain[0].domainsize());
    for (int y = 0; y < chain[0].domainrows; y++)
        for (int x = 0; x < chain[0].domaincolumns; x++)
            parentfield[chain[0].index(x, y)] =
                inputs.elevation_m[at(static_cast<int>(chain[0].absolutecolumn(x)), static_cast<int>(chain[0].absoluterow(y)),
                                      inputs.width)];

    // The parent level's sea (1 or 0) and lake surface on its domain (S2.8); level 0's are the
    // global sea and the global lake cells' surfaces (read for the diagnostics only: level 1
    // keeps the single pass's sea and its lakes have no parent core).
    std::vector<float> parentsea(chain[0].domainsize()), parentlake(chain[0].domainsize());
    for (int y = 0; y < chain[0].domainrows; y++)
        for (int x = 0; x < chain[0].domaincolumns; x++)
        {
            const std::size_t global = at(static_cast<int>(chain[0].absolutecolumn(x)),
                                          static_cast<int>(chain[0].absoluterow(y)), inputs.width);
            const bool sea = inputs.network.sea[global] != 0;
            parentsea[chain[0].index(x, y)] = sea ? 1.0f : 0.0f;
            parentlake[chain[0].index(x, y)] = !sea && inputs.lakesurface_m[global] > 0.0f ? inputs.lakesurface_m[global] : 0.0f;
        }
    details.lakes.flood = true;

    RefinementResult result;
    RegionalLevel state;
    RegionalWindow window;
    int sharepasses = -1, sharemisses = 0;
    for (int l = 1; l <= levels; l++)
    {
        const auto started = std::chrono::steady_clock::now();
        const LevelWindow& level = chain[static_cast<std::size_t>(l)];
        const LevelWindow& parent = chain[static_cast<std::size_t>(l - 1)];
        window = regionalwindowof(level);
        LevelDiagnostics diagnostics;
        diagnostics.level = l;
        diagnostics.cellsize_km = levelcellsize_km(inputs.radius_m, inputs.width, factor, l);
        diagnostics.window = level;

        // The start surface: level 1 the single-scale pass's (its warp, noise, coast decision
        // and block means against the global grid); deeper levels the parent prolonged, the
        // level's octave band added and the parent means corrected. The warp of a deeper level
        // is only the sea gate's and the lakes' position warp.
        const FineWarp warp = regionalwarp(window, inputs.seed, parameters);
        StartSurface start;
        if (l == 1)
        {
            start = regionalstartsurface(inputs, window, warp, parameters);
            sharepasses = start.sharepasses;
            sharemisses = start.sharemisses;
        }
        else
        {
            start.start_m = prolong(parent, parentfield, level);
            addlevelnoise(inputs, window, static_cast<int>(levelscale(factor, l - 1)), parameters, start.start_m);
            correctparentmeans(parent, parentfield, level, parameters.meanpasses, start.start_m);
            start.blockmeanerror_m = levelparentmeanerror(level, start.start_m, parent, parentfield);
            start.lapsebase_m = upsample(inputs.elevation_m, inputs.width, inputs.height, window);
        }
        diagnostics.blockmeanerror_m = start.blockmeanerror_m;
        diagnostics.steps = {LevelStep::noiseoctave, LevelStep::presmoothing, LevelStep::coast};
        if (incision.active)
            diagnostics.steps.push_back(LevelStep::incision);
        diagnostics.steps.push_back(LevelStep::channels);
        diagnostics.steps.push_back(LevelStep::solve);
        if (details.valleys)
            diagnostics.steps.push_back(LevelStep::valleyfloors);
        if (details.troughs || details.cirques)
            diagnostics.steps.push_back(LevelStep::glacialforms);
        diagnostics.steps.push_back(LevelStep::lakes);

        // The parent's masks per domain cell (S2.8): its parent and the parent's 8 neighbours
        // read from the parent's kept box; given to the sea and the lakes at level >= 2.
        const std::size_t count = level.domainsize();
        const ParentPatch seapatch = parentpatch(parent, parentsea, level, 1, "runrefinement (parent sea)");
        const ParentPatch lakepatch = parentpatch(parent, parentlake, level, 1, "runrefinement (parent lakes)");
        std::vector<uint8_t> seagate(count, 0), ownseaparent(count, 0), lakereach(count, 0);
        std::vector<float> lakeparent(count, 0.0f);
        for (int y = 0; y < level.domainrows; y++)
            for (int x = 0; x < level.domaincolumns; x++)
            {
                const std::size_t cell = level.index(x, y);
                const int pc = x / factor + 1, pr = y / factor + 1;
                bool sea = false, wet = false;
                for (int dr = -1; dr <= 1; dr++)
                    for (int dc = -1; dc <= 1; dc++)
                    {
                        const std::size_t neighbour = at(pc + dc, pr + dr, seapatch.columns);
                        sea = sea || seapatch.values[neighbour] != 0.0f;
                        wet = wet || lakepatch.values[neighbour] > 0.0f;
                    }
                seagate[cell] = sea ? 1 : 0;
                lakereach[cell] = wet ? 1 : 0;
                ownseaparent[cell] = seapatch.values[at(pc, pr, seapatch.columns)] != 0.0f ? 1 : 0;
                lakeparent[cell] = lakepatch.values[at(pc, pr, lakepatch.columns)];
            }
        details.parentseagate = l >= 2 ? &seagate : nullptr;
        details.lakes.parentsurface_m = l >= 2 ? &lakeparent : nullptr;
        details.lakes.parentreach = l >= 2 ? &lakereach : nullptr;
        state = regionallevelsteps(inputs, window, parameters, warp, std::move(start), incision, details);
        details.parentseagate = nullptr;
        details.lakes.parentsurface_m = nullptr;
        details.lakes.parentreach = nullptr;
        diagnostics.reliefconverged = state.reliefconverged;
        diagnostics.changedreceivers = state.changedreceivers;

        // The projection by bounds: sea and wet cells may sink but not rise past their water,
        // dry land may rise but not sink to the sea level or the water beside it (its own
        // global lake cell's surface included), land channel beds take the channel rate.
        std::vector<float>& elevation = state.elevation_m;
        std::vector<float>& lakesurface = state.water.lakesurface_m;
        std::vector<uint8_t> wet(count, 0);
        for (std::size_t cell = 0; cell < count; cell++)
            wet[cell] = lakesurface[cell] > 0.0f ? 1 : 0;
        std::vector<float> ownlake(count, 0.0f);
        for (int y = 0; y < level.domainrows; y++)
            for (int x = 0; x < level.domaincolumns; x++)
            {
                const std::size_t global = at(level.globalcolumn(x), level.globalrow(y), inputs.width);
                if (inputs.lakesurface_m[global] > 0.0f && inputs.network.sea[global] == 0)
                    ownlake[level.index(x, y)] = inputs.lakesurface_m[global];
            }
        ProjectionBounds bounds = projectionbounds(level, state.channels.channel, state.sea, lakesurface, elevation,
                                                   inputs.sealevel, refinement.projectionweightchannel, ownlake);
        // The smoothing before the projection (S2.7) with the same bounds; the projection's
        // bounds are then rebuilt on the smoothed field (the same floors and ceilings).
        diagnostics.smoothing = smoothparentmeans(level, elevation, bounds, parent, parentfield, refinement.smoothpasses,
                                                  refinement.smoothrelaxation);
        if (refinement.smoothpasses > 0)
        {
            diagnostics.steps.push_back(LevelStep::smoothing);
            bounds = projectionbounds(level, state.channels.channel, state.sea, lakesurface, elevation, inputs.sealevel,
                                      refinement.projectionweightchannel, ownlake);
        }
        diagnostics.projection = projectlevel(level, elevation, bounds, parent, parentfield);
        diagnostics.steps.push_back(LevelStep::projection);

        // After the projection: wet non-channel cells at or above their surface dried, dry
        // cells at or below a wet neighbour's surface counted (kept box only).
        const int keptx0 = level.keptcolumnoffset(), kepty0 = level.keptrowoffset();
        const auto kept = [&](int x, int y) {
            return x >= keptx0 && x < keptx0 + level.columns && y >= kepty0 && y < kepty0 + level.rows;
        };
        for (int y = 0; y < level.domainrows; y++)
            for (int x = 0; x < level.domaincolumns; x++)
            {
                const std::size_t cell = level.index(x, y);
                if (wet[cell] == 0 || state.channels.channel[cell] != 0 || elevation[cell] < lakesurface[cell])
                    continue;
                lakesurface[cell] = 0.0f;
                if (kept(x, y))
                    diagnostics.driedcells++;
            }
        for (int y = kepty0; y < kepty0 + level.rows; y++)
            for (int x = keptx0; x < keptx0 + level.columns; x++)
            {
                const std::size_t cell = level.index(x, y);
                if (lakesurface[cell] > 0.0f || state.sea[cell] != 0 || state.channels.channel[cell] != 0)
                    continue;
                bool below = false;
                for (int dy = -1; dy <= 1 && !below; dy++)
                    for (int dx = -1; dx <= 1 && !below; dx++)
                    {
                        const int nx = x + dx, ny = y + dy;
                        if ((dx == 0 && dy == 0) || nx < 0 || ny < 0 || nx >= level.domaincolumns || ny >= level.domainrows)
                            continue;
                        const float surface = lakesurface[level.index(nx, ny)];
                        below = surface > 0.0f && elevation[cell] <= surface;
                    }
                diagnostics.shoresbelowlake += below ? 1 : 0;
            }
        float largest = 0.0f;
        for (int y = kepty0; y < kepty0 + level.rows; y++)
            for (int x = keptx0; x < keptx0 + level.columns; x++)
            {
                const std::size_t cell = level.index(x, y);
                largest = std::max(largest, std::fabs(elevation[cell]));
                if (lakesurface[cell] <= 0.0f && state.sea[cell] == 0 && elevation[cell] <= inputs.sealevel)
                {
                    diagnostics.landbelowsea++;
                    const double depth = static_cast<double>(inputs.sealevel) - static_cast<double>(elevation[cell]);
                    if (depth > diagnostics.maxlandbelowsea_m || diagnostics.maxlandbelowseacolumn < 0)
                    {
                        diagnostics.maxlandbelowsea_m = depth;
                        diagnostics.maxlandbelowseacolumn = level.absolutecolumn(x);
                        diagnostics.maxlandbelowsearow = level.absoluterow(y);
                    }
                }
                if (!state.water.rimraise_m.empty() && state.water.rimraise_m[cell] > 0.0f)
                {
                    diagnostics.rimraised++;
                    diagnostics.maxrimraise_m = std::max(diagnostics.maxrimraise_m, static_cast<double>(state.water.rimraise_m[cell]));
                }
                if (state.sea[cell] != 0 && ownseaparent[cell] == 0)
                {
                    diagnostics.seainlandparents++;
                    diagnostics.seabeyondparent += seagate[cell] == 0 ? 1 : 0;
                }
                if (lakesurface[cell] > 0.0f || state.sea[cell] != 0 || state.channels.channel[cell] != 0)
                    continue;
                const int offsets[4][2] = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}};
                bool below = false;
                for (const auto& offset : offsets)
                {
                    const int nx = x + offset[0], ny = y + offset[1];
                    if (nx < 0 || ny < 0 || nx >= level.domaincolumns || ny >= level.domainrows)
                        continue;
                    const float surface = lakesurface[level.index(nx, ny)];
                    below = below || (surface > 0.0f && elevation[cell] < surface);
                }
                diagnostics.drybelowlake += below ? 1 : 0;
            }
        diagnostics.meanbound_m = largest > 0.0f ? std::ldexp(1.0, std::ilogb(largest) - 24) : 0.0;
        diagnostics.globalmeanerror_m = globalmeanerror(level, elevation, inputs);

        if (l < levels)
        {
            std::vector<float> keptfield(level.keptsize());
            for (int y = 0; y < level.rows; y++)
                for (int x = 0; x < level.columns; x++)
                    keptfield[at(x, y, level.columns)] = elevation[level.index(keptx0 + x, kepty0 + y)];
            result.coarserelevation_m.push_back(std::move(keptfield));
            std::vector<uint8_t> keptsea(level.keptsize());
            std::vector<float> keptlake(level.keptsize());
            for (int y = 0; y < level.rows; y++)
                for (int x = 0; x < level.columns; x++)
                {
                    keptsea[at(x, y, level.columns)] = state.sea[level.index(keptx0 + x, kepty0 + y)];
                    keptlake[at(x, y, level.columns)] = lakesurface[level.index(keptx0 + x, kepty0 + y)];
                }
            result.coarsersea.push_back(std::move(keptsea));
            result.coarserlakesurface_m.push_back(std::move(keptlake));
            parentsea.assign(count, 0.0f);
            for (std::size_t cell = 0; cell < count; cell++)
                parentsea[cell] = state.sea[cell] != 0 ? 1.0f : 0.0f;
            parentlake = lakesurface;
            parentfield = std::move(elevation);
        }
        else
            diagnostics.steps.push_back(LevelStep::ice);
        diagnostics.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
        result.levels.push_back(std::move(diagnostics));
    }

    // The ice above the ELA, the rivers and the climate at the final level.
    const auto started = std::chrono::steady_clock::now();
    result.fields = finishregionallevel(inputs, window, parameters, state);
    result.levels.back().seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
    result.fields.window = target.window;
    result.fields.window.scale = static_cast<int>(finest);
    result.fields.sharepasses = sharepasses;
    result.fields.sharemisses = sharemisses;
    result.chain = chain;
    return result;
}

int resamplesourcelevel(double metresperpixel, double radius_m, int worldwidth, int factor, int maxlevel)
{
    if (!std::isfinite(metresperpixel) || !(metresperpixel > 0.0) || !std::isfinite(radius_m) || !(radius_m > 0.0))
        throw std::invalid_argument("resamplesourcelevel: metres per pixel and the radius must be positive and finite");
    if (worldwidth < 4 || worldwidth % 2 != 0 || factor < 2 || maxlevel < 0)
        throw std::invalid_argument("resamplesourcelevel: an even world width of at least 4, a factor of at least 2 "
                                    "and a level not negative");
    // resampleto's expression: pi R / levelheight.
    for (int level = 0; level <= maxlevel; level++)
    {
        const std::int64_t rows = static_cast<std::int64_t>(worldwidth / 2) * levelscale(factor, level);
        if (pi * radius_m / static_cast<double>(rows) <= metresperpixel * (1.0 + 1.0e-9))
            return level;
    }
    return -1;
}

ResampledRegionalFields resampleregionalfields(const RefinementResult& result, double radius_m, double metresperpixel)
{
    const std::vector<LevelWindow>& chain = result.chain;
    if (chain.empty())
        throw std::invalid_argument("resampleregionalfields: the chain is empty");
    const LevelWindow& level = chain.back();
    const RegionalFields& source = result.fields;
    const std::size_t kept = level.keptsize();
    for (const std::size_t size : {source.elevation_m.size(), source.sea.size(), source.riverdirection.size(),
                                   source.discharge_m3s.size(), source.januaryflow.size(), source.julyflow.size(),
                                   source.lakesurface_m.size(), source.januarytemperature_c.size(),
                                   source.julytemperature_c.size(), source.januaryprecipitation_mm.size(),
                                   source.julyprecipitation_mm.size(), source.icethickness_m.size(), source.climate.size()})
        if (size != kept)
            throw std::invalid_argument("resampleregionalfields: a field is not the final level's kept size");

    // resampleto reads the source level's kept box only; the coarser levels' fields are
    // placeholders of their domain sizes, and the final level's kept field is placed in its
    // domain.
    std::vector<std::vector<float>> stack(chain.size());
    for (std::size_t l = 0; l < chain.size(); l++)
        stack[l].assign(chain[l].domainsize(), 0.0f);
    const int keptx0 = level.keptcolumnoffset(), kepty0 = level.keptrowoffset();
    const auto average = [&](const std::vector<float>& field)
    {
        std::vector<float>& domain = stack.back();
        for (int y = 0; y < level.rows; y++)
            for (int x = 0; x < level.columns; x++)
                domain[level.index(keptx0 + x, kepty0 + y)] = field[at(x, y, level.columns)];
        return resampleto(chain, stack, radius_m, metresperpixel);
    };

    ResampledRegionalFields resampled;
    ResampledField elevation = average(source.elevation_m);
    if (elevation.report.sourcelevel != level.level)
        throw std::invalid_argument("resampleregionalfields: the source level of " + std::to_string(metresperpixel) +
                                    " m per pixel is " + std::to_string(elevation.report.sourcelevel) +
                                    ", not the final level " + std::to_string(level.level));
    resampled.width = elevation.width;
    resampled.height = elevation.height;
    resampled.report = elevation.report;
    RegionalFields& fields = resampled.fields;
    fields.window = source.window;
    fields.blockmeanerror_m = source.blockmeanerror_m;
    fields.reliefconverged = source.reliefconverged;
    fields.changedreceivers = source.changedreceivers;
    fields.meanlowering_m = source.meanlowering_m;
    fields.maxlowering_m = source.maxlowering_m;
    fields.sharepasses = source.sharepasses;
    fields.sharemisses = source.sharemisses;
    fields.elevation_m = std::move(elevation.values);
    fields.januarytemperature_c = average(source.januarytemperature_c).values;
    fields.julytemperature_c = average(source.julytemperature_c).values;
    fields.januaryprecipitation_mm = average(source.januaryprecipitation_mm).values;
    fields.julyprecipitation_mm = average(source.julyprecipitation_mm).values;
    fields.icethickness_m = average(source.icethickness_m).values;

    // The discrete fields on resampleto's geometry: pixel (p, q) covers [p e, (p + 1) e) x
    // [q e, (q + 1) e) of the kept box in source cells, e = round(n) when n is whole, else n.
    const double n = resampled.report.pixelcells;
    const double whole = std::round(n);
    const double edge = std::fabs(n - whole) <= 1.0e-9 * n ? whole : n;
    const std::size_t pixels = static_cast<std::size_t>(resampled.width) * static_cast<std::size_t>(resampled.height);
    fields.sea.assign(pixels, 0);
    fields.riverdirection.assign(pixels, 0);
    fields.discharge_m3s.assign(pixels, 0.0f);
    fields.januaryflow.assign(pixels, 0);
    fields.julyflow.assign(pixels, 0);
    fields.lakesurface_m.assign(pixels, 0.0f);
    fields.climate.assign(pixels, 0);
    const int columns = level.columns, rows = level.rows;
    for (int q = 0; q < resampled.height; q++)
        for (int p = 0; p < resampled.width; p++)
        {
            const double west = p * edge, east = (p + 1) * edge, north = q * edge, south = (q + 1) * edge;
            const int ifirst = static_cast<int>(std::floor(west));
            const int ilast = std::min(columns - 1, static_cast<int>(std::ceil(east)) - 1);
            const int jfirst = static_cast<int>(std::floor(north));
            const int jlast = std::min(rows - 1, static_cast<int>(std::ceil(south)) - 1);
            double total = 0.0, seaweight = 0.0, wetweight = 0.0, surface = 0.0;
            float discharge = 0.0f;
            std::size_t largest = at(ifirst, jfirst, columns);
            bool found = false;
            int32_t january = 0, july = 0;
            for (int j = jfirst; j <= jlast; j++)
            {
                const double oy = std::min(static_cast<double>(j + 1), south) - std::max(static_cast<double>(j), north);
                if (oy <= 0.0)
                    continue;
                for (int i = ifirst; i <= ilast; i++)
                {
                    const double ox = std::min(static_cast<double>(i + 1), east) - std::max(static_cast<double>(i), west);
                    if (ox <= 0.0)
                        continue;
                    const double weight = ox * oy;
                    const std::size_t cell = at(i, j, columns);
                    total += weight;
                    if (source.sea[cell] != 0)
                        seaweight += weight;
                    else if (source.lakesurface_m[cell] > 0.0f)
                    {
                        wetweight += weight;
                        surface += weight * static_cast<double>(source.lakesurface_m[cell]);
                    }
                    if (!found || source.discharge_m3s[cell] > discharge)
                    {
                        discharge = source.discharge_m3s[cell];
                        largest = cell;
                    }
                    january = found ? std::max(january, source.januaryflow[cell]) : source.januaryflow[cell];
                    july = found ? std::max(july, source.julyflow[cell]) : source.julyflow[cell];
                    found = true;
                }
            }
            const std::size_t pixel = at(p, q, resampled.width);
            const bool sea = 2.0 * seaweight >= total;
            fields.sea[pixel] = sea ? 1 : 0;
            if (!sea && wetweight > 0.0 && 2.0 * wetweight >= total)
                fields.lakesurface_m[pixel] = static_cast<float>(surface / wetweight);
            fields.discharge_m3s[pixel] = discharge;
            fields.riverdirection[pixel] = source.riverdirection[largest];
            fields.januaryflow[pixel] = january;
            fields.julyflow[pixel] = july;
            const int ci = std::min(columns - 1, static_cast<int>(std::floor((p + 0.5) * edge)));
            const int cj = std::min(rows - 1, static_cast<int>(std::floor((q + 0.5) * edge)));
            fields.climate[pixel] = source.climate[at(ci, cj, columns)];
        }
    return resampled;
}
}
