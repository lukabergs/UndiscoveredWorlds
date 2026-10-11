// Tests of landscape/regional_folds.hpp (multiscale plan step S5.2: the fold fabric), suite
// Landscape: TEST_F(Landscape, RegionalFolds<What>) on a hand-built bench on window A's domain at
// k 8 and the synthetic world of regional_test_world.hpp with withfoldfields: off by default and
// at zero strain bit-identical, the amplitude and wavelength laws, a uniform field's wavelength,
// bearing and amplitude, the water gap of a channel line, the overlaps A/B and C/D with the folds
// on, and the refinement loop's level consistency and overlap with the folds drawn at two levels.

#include "landscape_fixture.hpp"

#include "regional_test_world.hpp"

#include "simulations/landscape/refinement_levels.hpp"
#include "simulations/landscape/regional_folds.hpp"

#include <cstring>
#include <limits>
#include <string>
#include <utility>

namespace
{
// The synthetic fold fields: a strain of 0 to 0.4 varying with longitude and latitude (so the
// threshold 0.1 cuts belts), an axis wandering over the map (kept in [0, pi)), a sediment cover of
// 1 to 9 km.
RegionalGlobalInputs withfoldfields(RegionalGlobalInputs inputs)
{
    const std::size_t size = static_cast<std::size_t>(inputs.width) * static_cast<std::size_t>(inputs.height);
    inputs.foldstrain.assign(size, 0.0f);
    inputs.folddirection_rad.assign(size, 0.0f);
    inputs.sediment_m.assign(size, 0.0f);
    for (int j = 0; j < inputs.height; j++)
        for (int i = 0; i < inputs.width; i++)
        {
            const std::size_t cell = at(i, j, inputs.width);
            const double phi = latitude_rad(j, inputs.height), lambda = longitude_rad(i, inputs.width);
            inputs.foldstrain[cell] = static_cast<float>(0.4 * (0.5 + 0.5 * std::sin(3.0 * lambda + 2.0 * phi)));
            double theta = 0.3 + 0.5 * std::sin(2.0 * lambda) + 0.4 * std::cos(3.0 * phi);
            theta = std::fmod(theta + 2.0 * pi, pi);
            inputs.folddirection_rad[cell] = static_cast<float>(theta);
            inputs.sediment_m[cell] = static_cast<float>(5000.0 + 4000.0 * std::sin(lambda + 2.0 * phi));
        }
    return inputs;
}

// The folds on the synthetic worlds, whose fine cells are 156 to 313 km at k 4: wavelengths of
// 150 x the cover in km clamped to 700 to 1,350 km (2.2 to 8.6 fine cells), amplitudes up to 300 m.
RegionalParameters foldparameters()
{
    RegionalParameters parameters;
    parameters.folds.amplitudemax_m = 300.0;
    parameters.folds.wavelengthmin_km = 700.0;
    parameters.folds.wavelengthmax_km = 1350.0;
    parameters.folds.wavelengthpercover_km = 150.0;
    parameters.folds.wavelengthbase_km = 0.0;
    parameters.folds.kernelwidth_km = 1350.0;
    return parameters;
}

bool samebits(const std::vector<float>& a, const std::vector<float>& b)
{
    return a.size() == b.size() && (a.empty() || std::memcmp(a.data(), b.data(), a.size() * sizeof(float)) == 0);
}

void expectsamefields(const RegionalFields& a, const RegionalFields& b, const std::string& label)
{
    EXPECT_TRUE(samebits(a.elevation_m, b.elevation_m)) << label << ": elevation";
    EXPECT_TRUE(samebits(a.discharge_m3s, b.discharge_m3s)) << label << ": discharge";
    EXPECT_TRUE(samebits(a.lakesurface_m, b.lakesurface_m)) << label << ": lakes";
    EXPECT_TRUE(samebits(a.januarytemperature_c, b.januarytemperature_c)) << label << ": temperature";
    EXPECT_TRUE(samebits(a.icethickness_m, b.icethickness_m)) << label << ": ice";
    EXPECT_EQ(a.sea, b.sea) << label << ": sea";
    EXPECT_EQ(a.riverdirection, b.riverdirection) << label << ": rivers";
}

// A hand-built bench on the domain of window A at k 8 (156 km cells): land at 1,000 m, no sea, no
// channel; a uniform strain, axis and cover; the global elevation 2,000 m everywhere (the
// elevation term 0.5 x 2,000 m does not bind the 100 m cap).
struct FoldBench
{
    RegionalGlobalInputs inputs = syntheticworld(32);
    RegionalWindow window = windowa(8);
    int width = window.domainfinewidth(), height = window.domainfineheight();
    std::size_t count = static_cast<std::size_t>(width) * static_cast<std::size_t>(height);
    landscape::FineChannels channels;
    std::vector<uint8_t> sea;
    std::vector<float> elevation;
    double cell_km = 0.0;
    RegionalParameters parameters;

    FoldBench(double wavelength_cells, float strain, float direction_rad)
    {
        const std::size_t size = static_cast<std::size_t>(inputs.width) * static_cast<std::size_t>(inputs.height);
        inputs.elevation_m.assign(size, 2000.0f);
        inputs.foldstrain.assign(size, strain);
        inputs.folddirection_rad.assign(size, direction_rad);
        inputs.sediment_m.clear();
        channels.channel.assign(count, 0);
        sea.assign(count, 0);
        elevation.assign(count, 1000.0f);
        cell_km = landscape::finepatchmetric(window, inputs.radius_m).dy_m / 1000.0;
        parameters.folds.amplitudemax_m = 100.0;
        parameters.folds.wavelengthmin_km = wavelength_cells * cell_km;
        parameters.folds.wavelengthmax_km = wavelength_cells * cell_km;
        parameters.folds.gapfraction = 0.0;
        parameters.folds.kernelwidth_km = 0.0;
    }

    std::size_t cell(int x, int y) const { return at(x, y, width); }
    void fold() { landscape::addfoldfabric(inputs, window, parameters, channels, sea, elevation); }
    double fabric(int x, int y) const { return static_cast<double>(elevation[cell(x, y)]) - 1000.0; }
    // The domain rows within 15 degrees of the equator.
    std::vector<int> equatorialrows() const
    {
        std::vector<int> rows;
        for (int y = 0; y < height; y++)
            if (std::fabs(window.finelatitude_rad(y)) < pi / 12.0)
                rows.push_back(y);
        return rows;
    }
};

// The number of kept cells whose elevation differs between two passes, and of them those in the
// global cells `other`'s kept box shares with a's window.
std::pair<int, int> foldchangedcells(const RegionalFields& a, const RegionalFields& b, const RegionalWindow& other)
{
    const RegionalWindow& window = a.window;
    const int k = window.scale, width = window.finewidth();
    int changed = 0, onoverlap = 0;
    for (std::size_t cell = 0; cell < a.elevation_m.size(); cell++)
    {
        if (a.elevation_m[cell] == b.elevation_m[cell])
            continue;
        changed++;
        const int x = static_cast<int>(cell % static_cast<std::size_t>(width));
        const int y = static_cast<int>(cell / static_cast<std::size_t>(width));
        const int i = (window.column0 + x / k) % window.worldwidth, j = window.row0 + y / k;
        const int column = ((i - other.column0) % window.worldwidth + window.worldwidth) % window.worldwidth;
        if (column < other.columns && j >= other.row0 && j < other.row0 + other.rows)
            onoverlap++;
    }
    return {changed, onoverlap};
}
}

TEST_F(Landscape, RegionalFoldsAreOffByDefault)
{
    // At the defaults (amplitude cap 0) the fields change nothing; with the folds on, an empty strain
    // field and a strain of zero everywhere change nothing either: every field bit-identical.
    const RegionalGlobalInputs plain = syntheticworld(32);
    const RegionalGlobalInputs folded = withfoldfields(plain);
    const RegionalWindow window = windowa(4);
    const RegionalFields reference = landscape::runregionalpass(plain, window, RegionalParameters{});
    expectsamefields(landscape::runregionalpass(folded, window, RegionalParameters{}), reference, "off by default");
    EXPECT_FALSE(landscape::foldfabricon(folded, RegionalParameters{}));
    const RegionalParameters on = foldparameters();
    EXPECT_FALSE(landscape::foldfabricon(plain, on));
    expectsamefields(landscape::runregionalpass(plain, window, on), reference, "no strain field");
    RegionalGlobalInputs zero = folded;
    std::fill(zero.foldstrain.begin(), zero.foldstrain.end(), 0.0f);
    EXPECT_TRUE(landscape::foldfabricon(zero, on));
    expectsamefields(landscape::runregionalpass(zero, window, on), reference, "zero strain");
    // And the folds on with the strain act.
    const RegionalFields active = landscape::runregionalpass(folded, window, on);
    EXPECT_GT(foldchangedcells(active, reference, window).first, 0);
}

TEST_F(Landscape, RegionalFoldsAmplitudeAndWavelengthLaws)
{
    landscape::FoldParameters folds;
    folds.amplitudemax_m = 1000.0;
    // L = clamp(10 + 2 x the cover in km, 5, 20); the base 10 km without the cover.
    EXPECT_NEAR(landscape::foldwavelength_km(folds, 4000.0, true), 18.0, 1.0e-12);
    EXPECT_NEAR(landscape::foldwavelength_km(folds, 1000.0, true), 12.0, 1.0e-12);
    EXPECT_NEAR(landscape::foldwavelength_km(folds, 20000.0, true), 20.0, 1.0e-12);
    EXPECT_NEAR(landscape::foldwavelength_km(folds, 0.0, true), 10.0, 1.0e-12);
    EXPECT_NEAR(landscape::foldwavelength_km(folds, 4000.0, false), 10.0, 1.0e-12);
    // Without the base (the design's law): clamp(2 x the cover, 5, 20).
    folds.wavelengthbase_km = 0.0;
    EXPECT_NEAR(landscape::foldwavelength_km(folds, 4000.0, true), 8.0, 1.0e-12);
    EXPECT_NEAR(landscape::foldwavelength_km(folds, 1000.0, true), 5.0, 1.0e-12);
    folds.wavelengthbase_km = 10.0;
    // The Zagros: L 12 km, eps 0.3: 0.5 x 12,000 / pi x sqrt(0.3) = 1,046 m, capped at 1,000.
    EXPECT_NEAR(landscape::foldamplitude_m(folds, 0.3, 12.0, 5000.0), 1000.0, 1.0e-9);
    EXPECT_NEAR(landscape::foldamplitude_m(folds, 0.3, 8.0, 5000.0), 0.5 * 8000.0 / pi * std::sqrt(0.3), 1.0e-9);
    // The elevation term: 300 m above the sea level stands 150 m.
    EXPECT_NEAR(landscape::foldamplitude_m(folds, 0.3, 12.0, 300.0), 150.0, 1.0e-9);
    EXPECT_EQ(landscape::foldamplitude_m(folds, 0.3, 12.0, -20.0), 0.0);
    // The weight smoothstep(eps; 0.1, 0.2): 0.5 at 0.15, 0 at 0.1 and below, 1 from 0.2.
    EXPECT_NEAR(landscape::foldamplitude_m(folds, 0.15, 12.0, 300.0), 75.0, 1.0e-9);
    EXPECT_EQ(landscape::foldamplitude_m(folds, 0.1, 12.0, 300.0), 0.0);
    EXPECT_NEAR(landscape::foldamplitude_m(folds, 0.2, 12.0, 300.0), 150.0, 1.0e-9);
    folds.amplitudemax_m = 0.0;
    EXPECT_EQ(landscape::foldamplitude_m(folds, 0.3, 12.0, 5000.0), 0.0);
}

TEST_F(Landscape, RegionalFoldsUniformStripes)
{
    // A uniform strain of 0.3 with the axis east (direction 0) and L = 8 fine cells: the crests run
    // north-south, so near the equator (cells about square) the row autocorrelation peaks at about
    // 8 cells and the fabric hardly changes along a column; the amplitude is A = min(0.5 (L / pi)
    // sqrt(0.3), 100, 0.5 x 2,000) = 100 m: the largest |fabric| within 10 percent of it, none above.
    // At the kernel width 0 (the design's form: each kernel's envelope its own L) and at twice L.
    const std::vector<int> rows = FoldBench(8.0, 0.3f, 0.0f).equatorialrows();
    ASSERT_GT(rows.size(), 4u);
    for (const double widthcells : {0.0, 16.0})
    {
        SCOPED_TRACE("kernel width " + std::to_string(widthcells) + " cells");
        FoldBench bench(8.0, 0.3f, 0.0f);
        bench.parameters.folds.kernelwidth_km = widthcells * bench.cell_km;
        const double amplitude = landscape::foldamplitude_m(bench.parameters.folds, 0.3, 8.0 * bench.cell_km, 2000.0);
        ASSERT_NEAR(amplitude, 100.0, 1.0e-9);
        bench.fold();
        double largest = 0.0;
        for (int y = 0; y < bench.height; y++)
            for (int x = 0; x < bench.width; x++)
                largest = std::max(largest, std::fabs(bench.fabric(x, y)));
        EXPECT_LE(largest, amplitude * (1.0 + 1.0e-6));
        EXPECT_GE(largest, 0.9 * amplitude);

        std::vector<double> correlation(25, 0.0);
        double alongrow = 0.0, alongcolumn = 0.0;
        int rowpairs = 0, columnpairs = 0;
        for (int y : rows)
        {
            for (int lag = 1; lag < static_cast<int>(correlation.size()); lag++)
            {
                double sum = 0.0;
                int pairs = 0;
                for (int x = 0; x + lag < bench.width; x++)
                {
                    sum += bench.fabric(x, y) * bench.fabric(x + lag, y);
                    pairs++;
                }
                correlation[static_cast<std::size_t>(lag)] += sum / pairs;
            }
            for (int x = 0; x + 1 < bench.width; x++)
            {
                alongrow += std::fabs(bench.fabric(x + 1, y) - bench.fabric(x, y));
                rowpairs++;
            }
            if (y + 1 < bench.height)
                for (int x = 0; x < bench.width; x++)
                {
                    alongcolumn += std::fabs(bench.fabric(x, y + 1) - bench.fabric(x, y));
                    columnpairs++;
                }
        }
        // The first peak after the first trough of the autocorrelation.
        int trough = 1;
        while (trough + 1 < static_cast<int>(correlation.size()) && correlation[static_cast<std::size_t>(trough + 1)] <
                                                                        correlation[static_cast<std::size_t>(trough)])
            trough++;
        // (the first local maximum: a coherent sinusoid repeats its peak at 2 L, 3 L with equal height)
        int peak = trough;
        while (peak + 1 < static_cast<int>(correlation.size()) &&
               correlation[static_cast<std::size_t>(peak + 1)] > correlation[static_cast<std::size_t>(peak)])
            peak++;
        alongrow /= rowpairs;
        alongcolumn /= columnpairs;
        std::printf("[folds] uniform stripes, kernel width %.0f cells: largest |fabric| %.3f m of A %.1f m; row autocorrelation trough at %d, peak at "
                    "%d cells (L 8); mean |step| along rows %.2f m, along columns %.2f m\n",
                    widthcells, largest, amplitude, trough, peak, alongrow, alongcolumn);
        EXPECT_NEAR(peak, 8, 1);
        EXPECT_NEAR(trough, 4, 1);
        EXPECT_LT(alongcolumn, 0.3 * alongrow);

        // The axis north-south (direction pi / 2): the crests run east-west, the steps swap.
        FoldBench turned(8.0, 0.3f, static_cast<float>(pi / 2.0));
        turned.parameters.folds.kernelwidth_km = widthcells * turned.cell_km;
        turned.fold();
        double turnedrow = 0.0, turnedcolumn = 0.0;
        for (int y : rows)
        {
            for (int x = 0; x + 1 < turned.width; x++)
                turnedrow += std::fabs(turned.fabric(x + 1, y) - turned.fabric(x, y));
            if (y + 1 < turned.height)
                for (int x = 0; x < turned.width; x++)
                    turnedcolumn += std::fabs(turned.fabric(x, y + 1) - turned.fabric(x, y));
        }
        EXPECT_LT(turnedrow / rowpairs, 0.3 * turnedcolumn / columnpairs);
    }

    // Sea, channels and the shore floor: a sea cell and a channel cell keep their values; land at
    // 30 m above the sea level is never lowered below the shore floor of 20 m; land at 10 m (below
    // the floor) is not changed.
    FoldBench kept(8.0, 0.3f, 0.0f);
    const int x0 = kept.width / 2, y0 = rows[rows.size() / 2];
    kept.sea[kept.cell(x0, y0)] = 1;
    kept.channels.channel[kept.cell(x0 + 1, y0)] = 1;
    for (int x = 0; x < kept.width; x++)
    {
        kept.elevation[kept.cell(x, y0 + 1)] = 30.0f;
        kept.elevation[kept.cell(x, y0 + 2)] = 10.0f;
    }
    kept.fold();
    for (int x = 0; x < kept.width; x++)
        EXPECT_EQ(kept.elevation[kept.cell(x, y0 + 2)], 10.0f) << "x = " << x;
    EXPECT_EQ(kept.elevation[kept.cell(x0, y0)], 1000.0f);
    EXPECT_EQ(kept.elevation[kept.cell(x0 + 1, y0)], 1000.0f);
    int floored = 0;
    for (int x = 0; x < kept.width; x++)
    {
        const float value = kept.elevation[kept.cell(x, y0 + 1)];
        EXPECT_GE(value, 20.0f) << "x = " << x;
        floored += value == 20.0f ? 1 : 0;
    }
    EXPECT_GT(floored, 0);
}

TEST_F(Landscape, RegionalFoldsChannelGap)
{
    // An east-west channel line across north-south crests (L = 10 cells, gap 0.25 L = 2.5 cells): the
    // channel cells keep their values, the fabric of the cells within 2.5 cells of the line is
    // reduced (rows 1 and 2 away), and the cells farther away are unchanged.
    FoldBench open(10.0, 0.3f, 0.0f);
    FoldBench gapped(10.0, 0.3f, 0.0f);
    gapped.parameters.folds.gapfraction = 0.25;
    open.parameters.folds.gapfraction = 0.25;
    const std::vector<int> rows = gapped.equatorialrows();
    const int y0 = rows[rows.size() / 2];
    for (int x = 0; x < gapped.width; x++)
        gapped.channels.channel[gapped.cell(x, y0)] = 1;
    open.fold();
    gapped.fold();
    int reduced = 0, compared = 0;
    for (int x = 0; x < gapped.width; x++)
    {
        EXPECT_EQ(gapped.elevation[gapped.cell(x, y0)], 1000.0f) << "channel cell x = " << x;
        for (int d : {-2, -1, 1, 2})
        {
            const double with = std::fabs(gapped.fabric(x, y0 + d)), without = std::fabs(open.fabric(x, y0 + d));
            if (without < 1.0)
                continue;
            compared++;
            reduced += with < without ? 1 : 0;
            EXPECT_LE(with, without) << "x = " << x << ", d = " << d;
        }
        for (int d : {-4, -3, 3, 4})
            EXPECT_EQ(gapped.elevation[gapped.cell(x, y0 + d)], open.elevation[open.cell(x, y0 + d)]) << "x = " << x << ", d = " << d;
    }
    std::printf("[folds] channel gap: %d of %d cells within 0.25 L of the line reduced\n", reduced, compared);
    EXPECT_GT(compared, 0);
    EXPECT_EQ(reduced, compared);
}

TEST_F(Landscape, RegionalOverlapFolds)
{
    // The folds on at the smallest halo the pass allows: the overlaps stay bit-identical (compareoverlap
    // asserts it), and the folds act on the shared cells.
    const RegionalParameters parameters = foldparameters();
    const RegionalGlobalInputs columns = withfoldfields(syntheticworld(32));
    for (int scale : {4, 5})
    {
        SCOPED_TRACE("k = " + std::to_string(scale));
        const OverlapCounts counts = compareoverlap(columns, windowa(scale), windowb(scale), "folds A/B", parameters);
        EXPECT_EQ(counts.shared, 4 * 5 * scale * scale);
        EXPECT_GT(counts.channel, 0);
        const RegionalFields on = landscape::runregionalpass(columns, windowb(scale), parameters);
        const RegionalFields off = landscape::runregionalpass(columns, windowb(scale), RegionalParameters{});
        const auto [changed, onoverlap] = foldchangedcells(on, off, windowa(scale));
        std::printf("[regional] folds B, k %d: %d of %zu kept cells changed, %d on the overlap with A\n", scale, changed,
                    on.elevation_m.size(), onoverlap);
        EXPECT_GT(changed, 0);
        EXPECT_GT(onoverlap, 0);
    }
    const RegionalGlobalInputs rows = withfoldfields(syntheticworld(64));
    const OverlapCounts counts = compareoverlap(rows, windowc(4), windowd(4), "folds C/D", parameters);
    EXPECT_GT(counts.sea, 0);
    EXPECT_GT(counts.channel, 0);
    const RegionalFields on = landscape::runregionalpass(rows, windowd(4), parameters);
    const RegionalFields off = landscape::runregionalpass(rows, windowd(4), RegionalParameters{});
    const auto [changed, onoverlap] = foldchangedcells(on, off, windowc(4));
    std::printf("[regional] folds D, k 4: %d of %zu kept cells changed, %d on the overlap with C\n", changed,
                on.elevation_m.size(), onoverlap);
    EXPECT_GT(changed, 0);
    EXPECT_GT(onoverlap, 0);
}

TEST_F(Landscape, RegionalFoldsRefinementConsistency)
{
    // The folds in the refinement loop on the 256-column synthetic world (levels 1 and 2 of 39.1 and
    // 9.77 km): wavelengths of 15 x the cover in km clamped to 20 to 120 km, so level 1 draws L >= 78 km
    // and level 2 the rest. Level 1 of two levels is the one-level result bit for bit; the k = 16 result
    // box-averaged equals both within the level-2 bound; windows P and Q (Q across the date line) agree
    // bit for bit on their overlap at both levels; both levels draw folds.
    const RegionalGlobalInputs inputs = withfoldfields(syntheticworld(256));
    RegionalParameters parameters;
    parameters.folds.amplitudemax_m = 500.0;
    parameters.folds.wavelengthmin_km = 20.0;
    parameters.folds.wavelengthmax_km = 120.0;
    parameters.folds.wavelengthpercover_km = 15.0;
    parameters.folds.wavelengthbase_km = 0.0;
    parameters.folds.kernelwidth_km = 120.0;
    const int halo = landscape::refinementhalo(2, 0);
    const RegionalWindow windowp = landscape::regionalwindow(256, 16, halo, 160.0, -15.0, 180.0, 10.0);
    const RegionalWindow windowq = landscape::regionalwindow(256, 16, halo, 170.0, -10.0, -170.0, 15.0);
    const auto run = [&](const RegionalWindow& window, int level) {
        landscape::LevelTarget target;
        target.window = window;
        target.level = level;
        landscape::RefinementResult result = landscape::runrefinement(inputs, target, parameters, landscape::RefinementParameters{});
        for (const landscape::LevelDiagnostics& diagnostics : result.levels)
            std::printf("[folds] L %d, level %d: %d fold cells, |fabric| p50 %.1f m p99 %.1f m; forced %d, border step "
                        "max %.1f m, land below the sea level %d\n",
                        level, diagnostics.level, diagnostics.foldcells, diagnostics.foldamplitude_p50_m,
                        diagnostics.foldamplitude_p99_m, diagnostics.projection.forcedparents,
                        diagnostics.projection.maxborderstep_m, diagnostics.landbelowsea);
        return result;
    };
    const landscape::RefinementResult one = run(windowp, 1);
    const landscape::RefinementResult two = run(windowp, 2);
    ASSERT_EQ(two.levels.size(), 2u);
    ASSERT_EQ(two.coarserelevation_m.size(), 1u);
    EXPECT_GT(two.levels[0].foldcells, 0);
    EXPECT_GT(two.levels[1].foldcells, 0);
    for (const landscape::LevelDiagnostics& diagnostics : two.levels)
        EXPECT_NE(std::find(diagnostics.steps.begin(), diagnostics.steps.end(), landscape::LevelStep::foldfabric),
                  diagnostics.steps.end());
    // The cells two level windows' kept boxes share, and how many of them differ in their bits.
    const auto overlap = [](const landscape::LevelWindow& a, const std::vector<float>& fa, const landscape::LevelWindow& b,
                            const std::vector<float>& fb) {
        const std::int64_t width = a.levelwidth();
        int shared = 0, differing = 0;
        for (int y = 0; y < a.rows; y++)
            for (int x = 0; x < a.columns; x++)
            {
                const std::int64_t column = landscape::wraplevel(a.column0 + x, width), row = a.row0 + y;
                const std::int64_t xb = landscape::wraplevel(column - b.column0, width), yb = row - b.row0;
                if (xb >= b.columns || yb < 0 || yb >= b.rows)
                    continue;
                shared++;
                const float va = fa[at(x, y, a.columns)], vb = fb[at(static_cast<int>(xb), static_cast<int>(yb), b.columns)];
                differing += std::memcmp(&va, &vb, sizeof(float)) != 0 ? 1 : 0;
            }
        return std::make_pair(shared, differing);
    };
    const auto [levelone, levelonediffer] =
        overlap(one.levels[0].window, one.fields.elevation_m, two.levels[0].window, two.coarserelevation_m[0]);
    std::printf("[folds] level 1 of two levels against one level: %d shared cells, %d differ\n", levelone, levelonediffer);
    EXPECT_EQ(static_cast<std::size_t>(levelone), one.levels[0].window.keptsize());
    EXPECT_EQ(levelonediffer, 0) << "level 1 of two levels";

    const int columns = windowp.columns * 4, rows = windowp.rows * 4;
    const int finewidth = two.fields.window.finewidth();
    ASSERT_EQ(finewidth, 4 * columns);
    double largest = 0.0;
    for (int y = 0; y < rows; y++)
        for (int x = 0; x < columns; x++)
        {
            double sum = 0.0;
            for (int fy = 4 * y; fy < 4 * y + 4; fy++)
                for (int fx = 4 * x; fx < 4 * x + 4; fx++)
                    sum += static_cast<double>(two.fields.elevation_m[at(fx, fy, finewidth)]);
            largest = std::max(largest, std::fabs(sum / 16.0 - static_cast<double>(one.fields.elevation_m[at(x, y, columns)])));
        }
    std::printf("[folds] consistency P: k = 16 box-averaged against k = 4 %.3g m (bound %.3g m)\n", largest,
                two.levels[1].meanbound_m);
    EXPECT_LE(largest, two.levels[1].meanbound_m);

    // The overlap of P and Q at both levels.
    const landscape::RefinementResult q = run(windowq, 2);
    ASSERT_EQ(q.levels.size(), 2u);
    for (int level = 1; level <= 2; level++)
    {
        const landscape::LevelWindow& a = two.levels[static_cast<std::size_t>(level - 1)].window;
        const landscape::LevelWindow& b = q.levels[static_cast<std::size_t>(level - 1)].window;
        const std::vector<float>& fa = level == 1 ? two.coarserelevation_m[0] : two.fields.elevation_m;
        const std::vector<float>& fb = level == 1 ? q.coarserelevation_m[0] : q.fields.elevation_m;
        const auto [shared, differing] = overlap(a, fa, b, fb);
        std::printf("[folds] P/Q, level %d: %d shared cells, %d differ\n", level, shared, differing);
        EXPECT_GT(shared, 0);
        EXPECT_EQ(differing, 0) << "level " << level;
    }
}
