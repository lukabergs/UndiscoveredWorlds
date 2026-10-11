// Tests of landscape/regional_folds.hpp (multiscale plan steps S5.2 and S5.3b: the fold fabric),
// suite Landscape: TEST_F(Landscape, RegionalFolds<What>) on a hand-built bench on window A's domain
// at k 8 and the synthetic world of regional_test_world.hpp with withfoldfields: off by default and
// at zero strain bit-identical, the amplitude and wavelength laws, a uniform field's wavelength,
// bearing and amplitude, the water gap of a channel line, a direction turning by 90 degrees over 8
// global cells (S5.3b), the ice and lake gates, the lakes decided on the surface without the fabric,
// the overlaps A/B and C/D with the folds on (phasor, Gabor, the ice gate acting), and the refinement
// loop's level consistency and overlap with the folds drawn at two levels.

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
    parameters.folds.splatspacing_km = 625.0;  // half the shortest wavelength: a 625 km lattice on 32 and 64 columns
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
// channel, no lake cell; a uniform strain, axis and cover; the global elevation 2,000 m everywhere
// (the elevation term 0.5 x 2,000 m does not bind the 100 m cap); the splat lattice at half the
// wavelength.
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
        inputs.lakesurface_m.assign(size, 0.0f);
        channels.channel.assign(count, 0);
        sea.assign(count, 0);
        elevation.assign(count, 1000.0f);
        cell_km = landscape::finepatchmetric(window, inputs.radius_m).dy_m / 1000.0;
        parameters.folds.amplitudemax_m = 100.0;
        parameters.folds.wavelengthmin_km = wavelength_cells * cell_km;
        parameters.folds.wavelengthmax_km = wavelength_cells * cell_km;
        parameters.folds.gapfraction = 0.0;
        parameters.folds.kernelwidth_km = 0.0;
        parameters.folds.splatspacing_km = 0.5 * wavelength_cells * cell_km;
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
    // S5.3b's relief fraction 0.15: A = 0.15 x 10,000 / pi x sqrt(eps) = 477 m x sqrt(eps) at L 10 km,
    // under the cap: 151 m at eps 0.1 (the weight then 0), 261 m at 0.3, 427 m at 0.8.
    EXPECT_NEAR(folds.relieffraction, 0.15, 1.0e-12);
    EXPECT_NEAR(landscape::foldamplitude_m(folds, 0.3, 10.0, 5000.0), 1500.0 / pi * std::sqrt(0.3), 1.0e-9);
    EXPECT_NEAR(landscape::foldamplitude_m(folds, 0.8, 10.0, 5000.0), 1500.0 / pi * std::sqrt(0.8), 1.0e-9);
    // The laws below at S5.2's fraction 0.5.
    folds.relieffraction = 0.5;
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
    // asserts it), and the folds act on the shared cells. Variants: the phasor form, the Gabor form
    // (phasor 0), and the phasor form with the glaciation of withglacialfields gating the fabric (S5.3b).
    for (const int variant : {0, 1, 2})
    {
        SCOPED_TRACE(variant == 0 ? "phasor" : variant == 1 ? "Gabor" : "ice gate");
        RegionalParameters parameters = foldparameters();
        if (variant == 1)
            parameters.folds.phasor = 0;
        const auto world = [variant](int width) {
            RegionalGlobalInputs inputs = withfoldfields(syntheticworld(width));
            if (variant == 2)
            {
                const RegionalGlobalInputs glacial = withglacialfields(inputs);
                inputs.glaciation = glacial.glaciation;
            }
            return inputs;
        };
        const RegionalGlobalInputs columns = world(32);
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
        const RegionalGlobalInputs rows = world(64);
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
}

TEST_F(Landscape, RegionalFoldsTurningDirection)
{
    // S5.3b: the across-strike axis turns by 90 degrees over 8 global cells (east of global column 130,
    // from 0, crests north-south, to pi / 2, crests east-west) on the 256-column world near the equator,
    // at k 64 (2.44 km cells) with L = 64 cells (one global cell), the lattice at L / 2 and sigma = 2 L
    // (the defaults' proportions: 10 km, 9.77 km, 20 km). Each splat interpolates the axis at its own
    // position, so the crests turn smoothly: between a cell and its 8-neighbour nearest the local strike
    // (the strike's angle as the field's, at most 22.5 degrees off it: an ideal sinusoid steps by at most
    // 2 pi x 0.54 / 64 A = 0.053 A there) the fabric steps by at most 0.1 A at p90, and the turn adds at
    // most 0.1 A to the uniform axis' p99. Not at every cell: with each splat's phase origin its own
    // position, the phasor sum has isolated zeros (the fold trains' dislocations), around which the
    // phase winds by 2 pi; the uniform axis' own along-strike steps (p99 about 0.17 A, the largest about
    // 1.8 A) are theirs.
    RegionalGlobalInputs inputs = syntheticworld(256);
    const RegionalWindow window = landscape::regionalwindow(256, 64, 1, 0.0, -2.8125, 16.875, 2.8125);
    const int width = window.domainfinewidth(), height = window.domainfineheight();
    const std::size_t count = static_cast<std::size_t>(width) * static_cast<std::size_t>(height);
    const std::size_t size = static_cast<std::size_t>(inputs.width) * static_cast<std::size_t>(inputs.height);
    const double cell_km = landscape::finepatchmetric(window, inputs.radius_m).dy_m / 1000.0;
    inputs.elevation_m.assign(size, 2000.0f);
    inputs.foldstrain.assign(size, 0.3f);
    inputs.sediment_m.clear();
    inputs.lakesurface_m.assign(size, 0.0f);
    const double first = 130.0;  // the turn's first global column (window columns 128 to 139)
    const auto axis = [first](double column) { return 0.5 * pi * std::clamp((column - first) / 8.0, 0.0, 1.0); };
    RegionalParameters parameters;
    parameters.folds.amplitudemax_m = 100.0;
    parameters.folds.wavelengthmin_km = 64.0 * cell_km;
    parameters.folds.wavelengthmax_km = 64.0 * cell_km;
    parameters.folds.gapfraction = 0.0;
    parameters.folds.splatspacing_km = 32.0 * cell_km;
    parameters.folds.kernelwidth_km = 128.0 * cell_km;
    landscape::FineChannels channels;
    channels.channel.assign(count, 0);
    const std::vector<uint8_t> sea(count, 0);

    struct Steps
    {
        double p50 = 0.0, p90 = 0.0, p99 = 0.0, largest = 0.0;
        int cells = 0;
    };
    const auto measure = [&](bool turning) {
        inputs.folddirection_rad.assign(size, 0.0f);
        if (turning)
            for (int j = 0; j < inputs.height; j++)
                for (int i = 0; i < inputs.width; i++)
                    inputs.folddirection_rad[at(i, j, inputs.width)] = static_cast<float>(axis(static_cast<double>(i)));
        std::vector<float> elevation(count, 1000.0f);
        landscape::addfoldfabric(inputs, window, parameters, channels, sea, elevation);
        std::vector<double> steps;
        const int x0 = window.keptfinecolumnoffset(), y0 = window.keptfinerowoffset();
        for (int y = y0; y < y0 + window.fineheight(); y++)
            for (int x = x0; x < x0 + window.finewidth(); x++)
            {
                // The strike (crests at the axis + pi / 2) at the cell's centre, in global columns.
                const double u = (static_cast<double>(window.absolutefinecolumn(x)) + 0.5) / 64.0 - 0.5;
                const double theta = turning ? axis(u) : 0.0;
                const double east = -std::sin(theta), south = std::cos(theta);
                int bestx = 0, besty = 0;
                double best = -2.0;
                for (int dy = -1; dy <= 1; dy++)
                    for (int dx = -1; dx <= 1; dx++)
                    {
                        if (dx == 0 && dy == 0)
                            continue;
                        const double alignment = (dx * east + dy * south) / std::hypot(dx, dy);
                        if (alignment > best)
                        {
                            best = alignment;
                            bestx = dx;
                            besty = dy;
                        }
                    }
                const int nx = x + bestx, ny = y + besty;
                if (nx < 0 || nx >= width || ny < 0 || ny >= height)
                    continue;
                steps.push_back(std::fabs(static_cast<double>(elevation[at(nx, ny, width)]) -
                                          static_cast<double>(elevation[at(x, y, width)])));
            }
        Steps result;
        std::sort(steps.begin(), steps.end());
        const auto rank = [&](double fraction) { return steps[static_cast<std::size_t>(fraction * (steps.size() - 1))]; };
        result.cells = static_cast<int>(steps.size());
        result.p50 = rank(0.5);
        result.p90 = rank(0.9);
        result.p99 = rank(0.99);
        result.largest = steps.back();
        return result;
    };
    const Steps uniform = measure(false);
    const Steps turning = measure(true);
    for (const auto& [label, steps] : {std::make_pair("uniform", uniform), std::make_pair("turning", turning)})
        std::printf("[folds] along-strike steps, %s axis: %d cells, p50 %.2f, p90 %.2f, p99 %.2f, max %.2f m of A 100 m\n",
                    label, steps.cells, steps.p50, steps.p90, steps.p99, steps.largest);
    EXPECT_LE(uniform.p90, 10.0);
    EXPECT_LE(turning.p90, 10.0);
    EXPECT_LE(turning.p99, uniform.p99 + 10.0);
}

TEST_F(Landscape, RegionalFoldsIceGate)
{
    // S5.3b: no fabric where the glaciation upsampled is above folds.icemax (0.5): a 3 x 3 group of
    // global cells at glaciation 1 leaves its cells above 0.5 unchanged; every other cell is folded as
    // without the field; icemax 1 is no gate.
    FoldBench plain(8.0, 0.3f, 0.0f);
    plain.fold();
    FoldBench iced(8.0, 0.3f, 0.0f);
    const std::size_t size = static_cast<std::size_t>(iced.inputs.width) * static_cast<std::size_t>(iced.inputs.height);
    iced.inputs.glaciation.assign(size, 0.0f);
    const int x0 = iced.width / 2, y0 = iced.equatorialrows()[iced.equatorialrows().size() / 2];
    const int i0 = iced.window.globalcolumn(x0), j0 = iced.window.globalrow(y0);
    for (int dj = -1; dj <= 1; dj++)
        for (int di = -1; di <= 1; di++)
            iced.inputs.glaciation[at((i0 + di + iced.inputs.width) % iced.inputs.width, j0 + dj, iced.inputs.width)] = 1.0f;
    const std::vector<float> cover =
        landscape::upsample(iced.inputs.glaciation, iced.inputs.width, iced.inputs.height, iced.window);
    FoldBench open = iced;
    open.parameters.folds.icemax = 1.0;
    iced.fold();
    open.fold();
    EXPECT_TRUE(samebits(open.elevation, plain.elevation)) << "icemax 1";
    int gated = 0, folded = 0;
    for (std::size_t cell = 0; cell < iced.count; cell++)
    {
        if (cover[cell] > 0.5f)
        {
            gated++;
            EXPECT_EQ(iced.elevation[cell], 1000.0f) << "cell " << cell;
        }
        else
        {
            EXPECT_EQ(iced.elevation[cell], plain.elevation[cell]) << "cell " << cell;
            folded += iced.elevation[cell] != 1000.0f ? 1 : 0;
        }
    }
    std::printf("[folds] ice gate: %d cells under the cover above 0.5 unchanged, %d others folded\n", gated, folded);
    EXPECT_GT(gated, 64);
    EXPECT_GT(folded, 0);
}

TEST_F(Landscape, RegionalFoldsLakeGate)
{
    // S5.3b: no fabric in the block of a global lake cell (lake surface above 0, not sea); the cells of
    // the other blocks are folded as without the lake.
    FoldBench plain(8.0, 0.3f, 0.0f);
    plain.fold();
    FoldBench lake(8.0, 0.3f, 0.0f);
    const int x0 = lake.width / 2, y0 = lake.equatorialrows()[lake.equatorialrows().size() / 2];
    const int i0 = lake.window.globalcolumn(x0), j0 = lake.window.globalrow(y0);
    const std::size_t global = at(i0, j0, lake.inputs.width);
    lake.inputs.lakesurface_m[global] = 500.0f;
    lake.inputs.network.sea[global] = 0;
    lake.fold();
    int inside = 0, changedinside = 0, foldedwithout = 0, outside = 0;
    for (int y = 0; y < lake.height; y++)
        for (int x = 0; x < lake.width; x++)
        {
            const std::size_t cell = lake.cell(x, y);
            if (lake.window.globalcolumn(x) == i0 && lake.window.globalrow(y) == j0)
            {
                inside++;
                changedinside += lake.elevation[cell] != 1000.0f ? 1 : 0;
                foldedwithout += plain.elevation[cell] != 1000.0f ? 1 : 0;
            }
            else
            {
                EXPECT_EQ(lake.elevation[cell], plain.elevation[cell]) << "x = " << x << ", y = " << y;
                outside += lake.elevation[cell] != 1000.0f ? 1 : 0;
            }
        }
    std::printf("[folds] lake gate: %d cells of the lake block, %d changed (%d without the lake); %d folded outside\n",
                inside, changedinside, foldedwithout, outside);
    EXPECT_EQ(inside, 64);
    EXPECT_EQ(changedinside, 0);
    EXPECT_GT(foldedwithout, 0);
    EXPECT_GT(outside, 0);
}

TEST_F(Landscape, RegionalFoldsLakesOnTheUnfoldedSurface)
{
    // S5.3b: the lakes step decides on the surface without the fabric, so in the single pass the lake
    // water of the synthetic world is the same with the folds on as with them off, and no folded dry cell
    // stands at or below the surface of a wet 8-neighbour.
    const RegionalGlobalInputs inputs = withfoldfields(syntheticworld(32));
    RegionalParameters off = RegionalParameters{};
    RegionalParameters on = foldparameters();
    for (RegionalParameters* parameters : {&off, &on})
        parameters->shores.shorenoise = 0.5;  // R4's lakes: the indicator, the terrain term and the shores
    for (int scale : {4, 8})
    {
        SCOPED_TRACE("k = " + std::to_string(scale));
        const RegionalWindow window = windowa(scale);
        const RegionalFields without = landscape::runregionalpass(inputs, window, off);
        const RegionalFields with = landscape::runregionalpass(inputs, window, on);
        int wet = 0, differing = 0;
        for (std::size_t cell = 0; cell < with.lakesurface_m.size(); cell++)
        {
            wet += with.lakesurface_m[cell] > 0.0f ? 1 : 0;
            differing += with.lakesurface_m[cell] != without.lakesurface_m[cell] ? 1 : 0;
        }
        const auto [changed, unused] = foldchangedcells(with, without, window);
        (void)unused;
        // A folded dry cell never at or below the surface of a wet 8-neighbour.
        const int width = window.finewidth(), height = window.fineheight();
        int shorecells = 0, sunk = 0;
        for (int y = 0; y < height; y++)
            for (int x = 0; x < width; x++)
            {
                const std::size_t cell = at(x, y, width);
                if (with.lakesurface_m[cell] > 0.0f || with.elevation_m[cell] == without.elevation_m[cell])
                    continue;
                float level = 0.0f;
                for (int dy = -1; dy <= 1; dy++)
                    for (int dx = -1; dx <= 1; dx++)
                        if (x + dx >= 0 && x + dx < width && y + dy >= 0 && y + dy < height)
                            level = std::max(level, with.lakesurface_m[at(x + dx, y + dy, width)]);
                if (level > 0.0f)
                {
                    shorecells++;
                    sunk += with.elevation_m[cell] > level ? 0 : 1;
                }
            }
        std::printf("[folds] lakes on the unfolded surface, k %d: %d wet cells, %d differ; %d cells folded, %d of them "
                    "beside the water, %d at or below it\n",
                    scale, wet, differing, changed, shorecells, sunk);
        EXPECT_GT(wet, 0);
        EXPECT_GT(changed, 0);
        EXPECT_EQ(differing, 0);
        EXPECT_EQ(sunk, 0);
    }
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
    parameters.folds.splatspacing_km = 40.0;  // the 156 km global cell by 4: a 39 km lattice
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
