// Tests of landscape/regional_folds.hpp (multiscale plan steps S5.2, S5.3b and S5.3c: the fold fabric),
// suite Landscape: TEST_F(Landscape, RegionalFolds<What>) on a hand-built bench on window A's domain
// at k 8, the synthetic world of regional_test_world.hpp with withfoldfields, and the phase field
// (foldphasefield, foldphase_rad) on the 256-column synthetic world sampled at k 64 near the equator:
// off by default and at zero strain bit-identical, the amplitude and wavelength laws, a uniform field's
// wavelength, bearing and amplitude, the water gap of a channel line, the phase of a uniform field, of a
// belt turning by 90 degrees, of a wavelength ramp and of a mask with a hole (S5.3c), the orientation
// flood (by continuity alone), the ice and lake gates, the lakes decided on the surface without the fabric, the overlaps A/B
// and C/D with the folds on (with and without the ice gate acting), and the refinement loop's level
// consistency and overlap with the folds drawn at two levels.

#include "landscape_fixture.hpp"

#include "regional_test_world.hpp"

#include "simulations/landscape/refinement_levels.hpp"
#include "simulations/landscape/regional_folds.hpp"

#include <algorithm>
#include <cmath>
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
// 300 x the cover in km clamped to 1,300 to 2,600 km (at least 4 cells of 313 km per wavelength: drawn
// in the single pass at k 4, 5 and 8), amplitudes up to 300 m.
RegionalParameters foldparameters()
{
    RegionalParameters parameters;
    parameters.folds.amplitudemax_m = 300.0;
    parameters.folds.wavelengthmin_km = 1300.0;
    parameters.folds.wavelengthmax_km = 2600.0;
    parameters.folds.wavelengthpercover_km = 300.0;
    parameters.folds.wavelengthbase_km = 0.0;
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
// channel, no lake cell; the global elevation 2,000 m everywhere (the elevation term 0.5 x 2,000 m does
// not bind the 100 m cap). Crests east-west: the axis north-south (pi / 2) and a strain of 0.3 over the
// whole planet (the phase is the latitude's, Phi = k y). Crests north-south: the axis east (0) on a band
// of global columns from 4 west of the domain to 2 east of it and the four global rows within 22.5
// degrees of the equator, the strain falling eastward from 0.25 + 0.03 (n - 1) to 0.25, so every cell
// but the band's two westernmost verges east (a band over every column would wind around the planet,
// which a single-valued phase cannot follow, and one from pole to pole would make its crests meridians:
// the converging meridians leave no uniform wavelength to north-south crests).
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

    FoldBench(double wavelength_cells, bool crestseastwest)
    {
        const std::size_t size = static_cast<std::size_t>(inputs.width) * static_cast<std::size_t>(inputs.height);
        inputs.elevation_m.assign(size, 2000.0f);
        inputs.sediment_m.clear();
        inputs.lakesurface_m.assign(size, 0.0f);
        if (crestseastwest)
        {
            inputs.foldstrain.assign(size, 0.3f);
            inputs.folddirection_rad.assign(size, static_cast<float>(pi / 2.0));
        }
        else
        {
            inputs.foldstrain.assign(size, 0.0f);
            inputs.folddirection_rad.assign(size, 0.0f);
            const int band = window.domaincolumns + 6;
            for (int offset = 0; offset < band; offset++)
            {
                const int column = ((window.domaincolumn0 - 4 + offset) % inputs.width + inputs.width) % inputs.width;
                for (int j = inputs.height / 2 - 2; j < inputs.height / 2 + 2; j++)
                    inputs.foldstrain[at(column, j, inputs.width)] = static_cast<float>(0.25 + 0.03 * (band - 1 - offset));
            }
        }
        channels.channel.assign(count, 0);
        sea.assign(count, 0);
        elevation.assign(count, 1000.0f);
        cell_km = landscape::finepatchmetric(window, inputs.radius_m).dy_m / 1000.0;
        parameters.folds.amplitudemax_m = 100.0;
        parameters.folds.wavelengthmin_km = wavelength_cells * cell_km;
        parameters.folds.wavelengthmax_km = wavelength_cells * cell_km;
        parameters.folds.gapfraction = 0.0;
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

// The 256-column bench of the phase field (156 km global cells, 2.44 km fine cells at k 64): a uniform global
// elevation of 2,000 m, no lakes, no sediment unless a test sets it; the strain and the axis set by the test.
RegionalGlobalInputs phaseworld()
{
    RegionalGlobalInputs inputs = syntheticworld(256);
    const std::size_t size = static_cast<std::size_t>(inputs.width) * static_cast<std::size_t>(inputs.height);
    inputs.elevation_m.assign(size, 2000.0f);
    inputs.foldstrain.assign(size, 0.0f);
    inputs.folddirection_rad.assign(size, 0.0f);
    inputs.sediment_m.clear();
    inputs.lakesurface_m.assign(size, 0.0f);
    return inputs;
}

// Samples of a phase field in global cell units with the ground metric of its world.
struct PhaseProbe
{
    const landscape::FoldPhaseField& field;
    double column_km = 0.0;

    explicit PhaseProbe(const landscape::FoldPhaseField& f)
        : field(f), column_km(2.0 * field.dy_km * static_cast<double>(field.height) / static_cast<double>(field.width))
    {
    }

    // The ground width of a global column at v (km).
    double eastkm(double v) const { return column_km * std::cos(0.5 * pi - pi * v / static_cast<double>(field.height)); }
    bool phase(double u, double v, double& value) const { return landscape::foldphase_rad(field, u, v, value); }
    // The local wavelength 2 pi / |grad Phi| (km) by central differences over h global cells; 0 where undefined.
    double wavelength(double u, double v, double h = 1.0 / 256.0) const
    {
        double west = 0.0, east = 0.0, north = 0.0, south = 0.0;
        if (!phase(u - h, v, west) || !phase(u + h, v, east) || !phase(u, v - h, north) || !phase(u, v + h, south))
            return 0.0;
        const double gradeast = (east - west) / (2.0 * h * eastkm(v)), gradsouth = (south - north) / (2.0 * h * field.dy_km);
        return 2.0 * pi / std::hypot(gradeast, gradsouth);
    }
    // |sin Phi(p + d (east, south)) - sin Phi(p)| for a ground step of d km along the unit (east, south); -1 where undefined.
    double step(double u, double v, double east, double south, double distance_km) const
    {
        double here = 0.0, there = 0.0;
        if (!phase(u, v, here) || !phase(u + distance_km * east / eastkm(v), v + distance_km * south / field.dy_km, there))
            return -1.0;
        return std::fabs(std::sin(there) - std::sin(here));
    }
};
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
    EXPECT_EQ(landscape::foldphasefield(zero, on.folds).maskcells, 0);
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
    // L = 8 fine cells and A = min(0.15 (L / pi) sqrt(eps), 100, 0.5 x 2,000) = 100 m. Crests north-south
    // (FoldBench's band) and east-west (the whole planet): across the crests the autocorrelation near the
    // equator peaks at about 8 cells, along them the fabric hardly changes; the largest |fabric| is within 10
    // percent of A, none above; every cell's fabric is A sin of the phase field at its centre.
    for (const bool eastwest : {false, true})
    {
        SCOPED_TRACE(eastwest ? "crests east-west" : "crests north-south");
        FoldBench bench(8.0, eastwest);
        const std::vector<int> rows = bench.equatorialrows();
        ASSERT_GT(rows.size(), 4u);
        const double amplitude = landscape::foldamplitude_m(bench.parameters.folds, 0.3, 8.0 * bench.cell_km, 2000.0);
        ASSERT_NEAR(amplitude, 100.0, 1.0e-9);
        bench.fold();
        const landscape::FoldPhaseField field = landscape::foldphasefield(bench.inputs, bench.parameters.folds);
        const std::vector<float> strain = landscape::upsample(bench.inputs.foldstrain, bench.inputs.width, bench.inputs.height, bench.window);
        double largest = 0.0, mismatch = 0.0;
        for (int y = 0; y < bench.height; y++)
            for (int x = 0; x < bench.width; x++)
            {
                if (bench.inputs.foldstrain[at(bench.window.globalcolumn(x), bench.window.globalrow(y), bench.inputs.width)] <= 0.0f)
                    continue;  // the band's edge
                largest = std::max(largest, std::fabs(bench.fabric(x, y)));
                double phase = 0.0;
                const double u = (static_cast<double>(bench.window.absolutefinecolumn(x)) + 0.5) / 8.0;
                const double v = (static_cast<double>(bench.window.absolutefinerow(y)) + 0.5) / 8.0;
                const double own = landscape::foldamplitude_m(bench.parameters.folds, strain[bench.cell(x, y)], 8.0 * bench.cell_km, 2000.0);
                if (landscape::foldphase_rad(field, u, v, phase))
                    mismatch = std::max(mismatch, std::fabs(bench.fabric(x, y) - own * std::sin(phase)));
            }
        EXPECT_LE(largest, amplitude * (1.0 + 1.0e-6));
        EXPECT_GE(largest, 0.9 * amplitude);
        EXPECT_LE(mismatch, 1.0e-3);

        // Across and along the crests near the equator.
        std::vector<double> correlation(25, 0.0);
        double alongcrest = 0.0, acrosscrest = 0.0;
        int alongpairs = 0, acrosspairs = 0, lines = 0;
        // Lines across the crests: rows (north-south crests) or the columns of the equatorial rows' span.
        const int linecount = eastwest ? bench.width : static_cast<int>(rows.size());
        const int linelength = eastwest ? static_cast<int>(rows.size()) : bench.width;
        const auto value = [&](int line, int position) {
            return eastwest ? bench.fabric(line, rows[static_cast<std::size_t>(position)]) : bench.fabric(position, rows[static_cast<std::size_t>(line)]);
        };
        // The east-west crests' lines run over the equatorial rows only: lengthen them to the domain's rows within 30
        // degrees of the equator for the autocorrelation.
        std::vector<int> longrows;
        for (int y = 0; y < bench.height; y++)
            if (std::fabs(bench.window.finelatitude_rad(y)) < pi / 6.0)
                longrows.push_back(y);
        for (int line = 0; line < linecount; line++)
        {
            if (!eastwest && bench.inputs.foldstrain[at(bench.window.globalcolumn(0), bench.window.globalrow(rows[static_cast<std::size_t>(line)]), bench.inputs.width)] <= 0.0f)
                continue;
            lines++;
            for (int lag = 1; lag < static_cast<int>(correlation.size()); lag++)
            {
                double sum = 0.0;
                int pairs = 0;
                if (eastwest)
                    for (std::size_t p = 0; p + static_cast<std::size_t>(lag) < longrows.size(); p++)
                    {
                        sum += bench.fabric(line, longrows[p]) * bench.fabric(line, longrows[p + static_cast<std::size_t>(lag)]);
                        pairs++;
                    }
                else
                    for (int x = 0; x + lag < bench.width; x++)
                    {
                        sum += value(line, x) * value(line, x + lag);
                        pairs++;
                    }
                correlation[static_cast<std::size_t>(lag)] += sum / pairs;
            }
            for (int p = 0; p + 1 < linelength; p++)
            {
                acrosscrest += std::fabs(value(line, p + 1) - value(line, p));
                acrosspairs++;
            }
        }
        // Along the crests: between neighbours along columns (north-south crests) or rows (east-west).
        for (int y : rows)
            for (int x = 0; x + 1 < bench.width; x++)
            {
                if (eastwest)
                    alongcrest += std::fabs(bench.fabric(x + 1, y) - bench.fabric(x, y));
                else if (y + 1 < bench.height)
                    alongcrest += std::fabs(bench.fabric(x, y + 1) - bench.fabric(x, y));
                alongpairs++;
            }
        ASSERT_GT(lines, 0);
        int trough = 1;
        while (trough + 1 < static_cast<int>(correlation.size()) && correlation[static_cast<std::size_t>(trough + 1)] <
                                                                        correlation[static_cast<std::size_t>(trough)])
            trough++;
        int peak = trough;
        while (peak + 1 < static_cast<int>(correlation.size()) &&
               correlation[static_cast<std::size_t>(peak + 1)] > correlation[static_cast<std::size_t>(peak)])
            peak++;
        alongcrest /= alongpairs;
        acrosscrest /= acrosspairs;
        std::printf("[folds] uniform stripes, %s: largest |fabric| %.3f m of A %.1f m (A sin phase within %.2g m); "
                    "autocorrelation across the crests: trough at %d, peak at %d cells (L 8); mean |step| across %.2f m, along %.2f m\n",
                    eastwest ? "crests east-west" : "crests north-south", largest, amplitude, mismatch, trough, peak, acrosscrest,
                    alongcrest);
        EXPECT_NEAR(peak, 8, 1);
        EXPECT_NEAR(trough, 4, 1);
        EXPECT_LT(alongcrest, 0.3 * acrosscrest);
    }

    // Sea, channels and the shore floor: a sea cell and a channel cell keep their values; land at
    // 30 m above the sea level is never lowered below the shore floor of 20 m; land at 10 m (below
    // the floor) is not changed.
    FoldBench kept(8.0, false);
    const std::vector<int> rows = kept.equatorialrows();
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
    FoldBench open(10.0, false);
    FoldBench gapped(10.0, false);
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
    // asserts it), and the folds act on the shared cells. Variants: the folds alone, and with the
    // glaciation of withglacialfields gating the fabric (S5.3b).
    for (const int variant : {0, 1})
    {
        SCOPED_TRACE(variant == 0 ? "folds" : "ice gate");
        const RegionalParameters parameters = foldparameters();
        const auto world = [variant](int width) {
            RegionalGlobalInputs inputs = withfoldfields(syntheticworld(width));
            if (variant == 1)
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

TEST_F(Landscape, RegionalFoldsPhaseUniform)
{
    // S5.3c: a uniform field, L = 10 km (4.1 fine cells at k 64) and a strain of 0.3. The axis north-south
    // over the whole planet (crests east-west): the phase is k times the distance south, so on the k 64
    // centres of 4 x 12 global cells at the equator the stripes are straight (along-strike steps of one fine
    // cell below 1e-3 A) and the local wavelength is L within 1 percent. The axis east on a band of 32 x 24
    // global cells about the equator (crests north-south; the strain falling eastward orients it, the band's
    // west edge verging west): reported over the same box, the sphere's meridians converging (a uniform
    // wavelength north-south crests cannot keep).
    for (const bool eastwest : {true, false})
    {
        SCOPED_TRACE(eastwest ? "crests east-west" : "crests north-south");
        RegionalGlobalInputs inputs = phaseworld();
        landscape::FoldParameters folds;
        folds.amplitudemax_m = 100.0;
        if (eastwest)
        {
            std::fill(inputs.foldstrain.begin(), inputs.foldstrain.end(), 0.3f);
            std::fill(inputs.folddirection_rad.begin(), inputs.folddirection_rad.end(), static_cast<float>(pi / 2.0));
        }
        else
            for (int j = 52; j < 76; j++)
                for (int i = 116; i < 148; i++)
                    inputs.foldstrain[at(i, j, inputs.width)] = static_cast<float>(0.25 + 0.025 * (147 - i));
        const landscape::FoldPhaseField field = landscape::foldphasefield(inputs, folds);
        const PhaseProbe probe(field);
        const double fine_km = field.dy_km / 64.0;
        double steplargest = 0.0, wavelengthlow = 1.0e9, wavelengthhigh = 0.0;
        int samples = 0;
        for (int y = 0; y < 4 * 64; y++)
            for (int x = 0; x < 12 * 64; x++)
            {
                const double u = 128.0 + (x + 0.5) / 64.0, v = 62.0 + (y + 0.5) / 64.0;
                const double step = eastwest ? probe.step(u, v, 1.0, 0.0, fine_km) : probe.step(u, v, 0.0, 1.0, fine_km);
                const double wavelength = probe.wavelength(u, v);
                ASSERT_GE(step, 0.0);
                ASSERT_GT(wavelength, 0.0);
                steplargest = std::max(steplargest, step);
                wavelengthlow = std::min(wavelengthlow, wavelength);
                wavelengthhigh = std::max(wavelengthhigh, wavelength);
                samples++;
            }
        std::printf("[folds] phase, uniform field, %s: %d samples; along-strike steps at most %.2e A; local wavelength %.4f to %.4f "
                    "km (L 10); mask %d, oriented by the flood %d, north %d; %d sweeps, last change %.2g rad, misfit "
                    "%.4f, %d reversal edges, %.3f s\n",
                    eastwest ? "crests east-west" : "crests north-south", samples, steplargest, wavelengthlow, wavelengthhigh,
                    field.maskcells, field.orientedbyflood, field.orientedbynorth, field.sweeps,
                    field.lastchange_rad, field.misfit, field.reversaledges, field.seconds);
        if (eastwest)
        {
            EXPECT_LT(steplargest, 1.0e-3);
            EXPECT_NEAR(wavelengthlow, 10.0, 0.1);
            EXPECT_NEAR(wavelengthhigh, 10.0, 0.1);
        }
        else
        {
            EXPECT_LT(steplargest, 0.02);
            EXPECT_NEAR(wavelengthlow, 10.0, 0.5);
            EXPECT_NEAR(wavelengthhigh, 10.0, 0.5);
        }
    }
}

TEST_F(Landscape, RegionalFoldsPhaseTurning)
{
    // S5.3c: a belt turning by 90 degrees over about 8 global cells, an orocline: the global cells whose centres
    // lie 3 to 8 cells (on the ground) from a centre at (128, 60), from 10 degrees north of east to 10 degrees
    // west of south of it; the axis across the belt (radial), the strain falling outward from 0.5 to 0.25, oriented
    // by the flood (no reversal); L one global cell (64 fine cells). The arc-shaped belt (an orocline) is the turning
    // test the orchestrator confirmed on 2026-10-11: an axis that varies only east-west with a uniform L has a curl
    // no single-valued phase can follow over a tall mask (the solve fans the crests out). Over the k 64 positions 5 to 7 cells from the centre between east and
    // south: no step above
    // 0.02 A along the strike (one fine cell along the circle), the local wavelength within 10 percent of L.
    RegionalGlobalInputs inputs = phaseworld();
    const double uc = 128.0, vc = 60.0;
    const auto local = [&](double u, double v, double& x, double& y) {
        x = (u - uc) * std::cos(0.5 * pi - pi * v / static_cast<double>(inputs.height));
        y = v - vc;
    };
    for (int j = 0; j < inputs.height; j++)
        for (int i = 0; i < inputs.width; i++)
        {
            double x = 0.0, y = 0.0;
            local(i + 0.5, j + 0.5, x, y);
            const double r = std::hypot(x, y), angle = std::atan2(y, x);
            if (r < 3.0 || r > 8.0 || angle < -pi / 18.0 || angle > pi / 2.0 + pi / 18.0)
                continue;
            inputs.foldstrain[at(i, j, inputs.width)] = static_cast<float>(0.25 + 0.05 * (8.0 - r));
            inputs.folddirection_rad[at(i, j, inputs.width)] = static_cast<float>(std::fmod(angle + pi, pi));
        }
    landscape::FoldParameters folds;
    folds.amplitudemax_m = 100.0;
    const double length_km = pi * inputs.radius_m / static_cast<double>(inputs.height) / 1000.0;
    folds.wavelengthmin_km = length_km;
    folds.wavelengthmax_km = length_km;
    const landscape::FoldPhaseField field = landscape::foldphasefield(inputs, folds);
    EXPECT_EQ(field.reversaledges, 0);
    const PhaseProbe probe(field);
    const double fine_km = field.dy_km / 64.0;
    double steplargest = 0.0, wavelengthlow = 1.0e18, wavelengthhigh = 0.0;
    int samples = 0;
    for (int y = 0; y < 12 * 64; y++)
        for (int x = 0; x < 12 * 64; x++)
        {
            const double u = uc - 2.0 + (x + 0.5) / 64.0, v = vc - 2.0 + (y + 0.5) / 64.0;
            double ex = 0.0, ey = 0.0;
            local(u, v, ex, ey);
            const double r = std::hypot(ex, ey), angle = std::atan2(ey, ex);
            if (r < 5.0 || r > 7.0 || angle < 0.0 || angle > pi / 2.0)
                continue;
            const double step = probe.step(u, v, -std::sin(angle), std::cos(angle), fine_km);
            const double wavelength = probe.wavelength(u, v);
            ASSERT_GE(step, 0.0);
            ASSERT_GT(wavelength, 0.0);
            steplargest = std::max(steplargest, step);
            wavelengthlow = std::min(wavelengthlow, wavelength);
            wavelengthhigh = std::max(wavelengthhigh, wavelength);
            samples++;
        }
    std::printf("[folds] phase, turning belt: %d samples; along-strike steps at most %.4f A; local wavelength %.1f to %.1f km (L %.1f); "
                "mask %d, oriented by the flood %d; %d sweeps, misfit %.4f, %d reversal edges\n",
                samples, steplargest, wavelengthlow, wavelengthhigh, length_km, field.maskcells,
                field.orientedbyflood, field.sweeps, field.misfit, field.reversaledges);
    EXPECT_GT(samples, 10000);
    EXPECT_LE(steplargest, 0.02);
    EXPECT_GE(wavelengthlow, 0.9 * length_km);
    EXPECT_LE(wavelengthhigh, 1.1 * length_km);
}

TEST_F(Landscape, RegionalFoldsPhaseWavelengthRamp)
{
    // S5.3c: the axis east on a band of 32 x 24 global cells about the equator, the strain falling eastward (the
    // band verges east), L from 10 to 20 km across it (the sediment cover 0 to 5 km over global columns 124 to 140:
    // L = 10 + 2 x the cover in km). On every k 64 row of 4 global rows at the equator, between global columns
    // 124.5 and 139.5: the phase rises by less than pi from one fine centre to the next (no jump), the local
    // wavelength never falls (monotone within 0.1 percent), and every row crosses as many crests (no dislocation).
    RegionalGlobalInputs inputs = phaseworld();
    inputs.sediment_m.assign(inputs.foldstrain.size(), 0.0f);
    for (int j = 52; j < 76; j++)
        for (int i = 116; i < 148; i++)
        {
            inputs.foldstrain[at(i, j, inputs.width)] = static_cast<float>(0.25 + 0.025 * (147 - i));
            inputs.sediment_m[at(i, j, inputs.width)] = static_cast<float>(5000.0 * std::clamp((i + 0.5 - 124.5) / 16.0, 0.0, 1.0));
        }
    landscape::FoldParameters folds;
    folds.amplitudemax_m = 100.0;
    const landscape::FoldPhaseField field = landscape::foldphasefield(inputs, folds);
    const PhaseProbe probe(field);
    int jumps = 0, falls = 0, rows = 0, crestslow = 1 << 30, crestshigh = 0;
    double first = 0.0, last = 0.0;
    for (int y = 0; y < 4 * 64; y++)
    {
        const double v = 62.0 + (y + 0.5) / 64.0;
        double previous = 0.0, previouswavelength = 0.0;
        int crests = 0;
        for (int x = 0; x <= 15 * 64; x++)
        {
            const double u = 124.5 + x / 64.0;
            double phase = 0.0;
            ASSERT_TRUE(probe.phase(u, v, phase));
            const double wavelength = probe.wavelength(u, v);
            if (x > 0)
            {
                jumps += !(phase > previous && phase - previous < pi) ? 1 : 0;
                falls += wavelength < previouswavelength * (1.0 - 1.0e-3) ? 1 : 0;
                crests += static_cast<int>(std::floor((phase - pi / 2.0) / (2.0 * pi)) - std::floor((previous - pi / 2.0) / (2.0 * pi)));
            }
            if (y == 128 && x == 0)
                first = wavelength;
            if (y == 128 && x == 15 * 64)
                last = wavelength;
            previous = phase;
            previouswavelength = wavelength;
        }
        rows++;
        crestslow = std::min(crestslow, crests);
        crestshigh = std::max(crestshigh, crests);
    }
    std::printf("[folds] phase, wavelength ramp: %d rows, %d to %d crests per row, %d jumps, %d falls of the local wavelength; "
                "L %.2f km at column 124.5, %.2f km at 139.5; %d sweeps, misfit %.4f, %d reversal edges\n",
                rows, crestslow, crestshigh, jumps, falls, first, last, field.sweeps, field.misfit, field.reversaledges);
    EXPECT_EQ(jumps, 0);
    EXPECT_EQ(falls, 0);
    EXPECT_EQ(crestslow, crestshigh);
    EXPECT_NEAR(first, 10.0, 0.5);
    EXPECT_NEAR(last, 20.0, 1.0);
}

TEST_F(Landscape, RegionalFoldsPhaseAroundAHole)
{
    // S5.3c: a disk of strain 0.3 (the global cells within 8 cells on the ground of (128, 64)) with a hole of 2.5
    // cells at its centre, the axis radial, L one global cell; the flood orients every cell from the first (north) and must carry one sign around the hole: the crests are rings
    // about it. Over the k 64 positions 3.5 to 7 cells from the centre at every angle: no step above 0.02 A
    // along the strike, the local wavelength within 10 percent of L, and the phase rising by 3.5 wavelengths within
    // 10 percent along every radius from 3.5 to 7 cells.
    RegionalGlobalInputs inputs = phaseworld();
    const double uc = 128.0, vc = 64.0;
    const auto local = [&](double u, double v, double& x, double& y) {
        x = (u - uc) * std::cos(0.5 * pi - pi * v / static_cast<double>(inputs.height));
        y = v - vc;
    };
    for (int j = 0; j < inputs.height; j++)
        for (int i = 0; i < inputs.width; i++)
        {
            double x = 0.0, y = 0.0;
            local(i + 0.5, j + 0.5, x, y);
            const double r = std::hypot(x, y);
            if (r < 2.5 || r > 8.0)
                continue;
            inputs.foldstrain[at(i, j, inputs.width)] = 0.3f;
            inputs.folddirection_rad[at(i, j, inputs.width)] = static_cast<float>(std::fmod(std::atan2(y, x) + pi, pi));
        }
    landscape::FoldParameters folds;
    folds.amplitudemax_m = 100.0;
    const double length_km = pi * inputs.radius_m / static_cast<double>(inputs.height) / 1000.0;
    folds.wavelengthmin_km = length_km;
    folds.wavelengthmax_km = length_km;
    const landscape::FoldPhaseField field = landscape::foldphasefield(inputs, folds);
    EXPECT_EQ(field.orientedbynorth, 1);
    EXPECT_EQ(field.components, 1);
    // One sign around the hole: every mask cell's oriented axis points the same way along the radius.
    int outward = 0, inward = 0;
    for (int j = 0; j < inputs.height; j++)
        for (int i = 0; i < inputs.width; i++)
        {
            const std::size_t cell = at(i, j, inputs.width);
            if (field.state[cell] != 1)
                continue;
            double x = 0.0, y = 0.0;
            local(i + 0.5, j + 0.5, x, y);
            const double radial = field.wavenumbereast[cell] * x + field.wavenumbersouth[cell] * y;
            outward += radial > 0.0 ? 1 : 0;
            inward += radial < 0.0 ? 1 : 0;
        }
    EXPECT_TRUE(outward == 0 || inward == 0) << outward << " outward, " << inward << " inward";
    const PhaseProbe probe(field);
    const double fine_km = field.dy_km / 64.0;
    double steplargest = 0.0, wavelengthlow = 1.0e18, wavelengthhigh = 0.0;
    int samples = 0;
    for (int y = 0; y < 16 * 64; y++)
        for (int x = 0; x < 16 * 64; x++)
        {
            const double u = uc - 8.0 + (x + 0.5) / 64.0, v = vc - 8.0 + (y + 0.5) / 64.0;
            double ex = 0.0, ey = 0.0;
            local(u, v, ex, ey);
            const double r = std::hypot(ex, ey), angle = std::atan2(ey, ex);
            if (r < 3.5 || r > 7.0)
                continue;
            const double step = probe.step(u, v, -std::sin(angle), std::cos(angle), fine_km);
            const double wavelength = probe.wavelength(u, v);
            ASSERT_GE(step, 0.0);
            ASSERT_GT(wavelength, 0.0);
            steplargest = std::max(steplargest, step);
            wavelengthlow = std::min(wavelengthlow, wavelength);
            wavelengthhigh = std::max(wavelengthhigh, wavelength);
            samples++;
        }
    // The phase's rise from r 3.5 to 7 along 64 radii: 3.5 wavelengths on every one (no crest ends between them).
    double riselow = 1.0e18, risehigh = 0.0;
    for (int ray = 0; ray < 64; ray++)
    {
        const double angle = 2.0 * pi * (ray + 0.5) / 64.0;
        double ends[2] = {0.0, 0.0};
        for (int end = 0; end < 2; end++)
        {
            const double r = end == 0 ? 3.5 : 7.0;
            const double v = vc + r * std::sin(angle);
            const double u = uc + r * std::cos(angle) / std::cos(0.5 * pi - pi * v / static_cast<double>(inputs.height));
            ASSERT_TRUE(probe.phase(u, v, ends[end]));
        }
        const double rise = std::fabs(ends[1] - ends[0]) / (2.0 * pi);
        riselow = std::min(riselow, rise);
        risehigh = std::max(risehigh, rise);
    }
    std::printf("[folds] phase around a hole: %d samples; along-strike steps at most %.4f A; local wavelength %.1f to %.1f km (L %.1f); "
                "%.3f to %.3f wavelengths from r 3.5 to 7 (3.5); %d sweeps, misfit %.4f, %d reversal edges\n",
                samples, steplargest, wavelengthlow, wavelengthhigh, length_km, riselow, risehigh, field.sweeps, field.misfit,
                field.reversaledges);
    EXPECT_LE(steplargest, 0.02);
    EXPECT_GE(wavelengthlow, 0.9 * length_km);
    EXPECT_LE(wavelengthhigh, 1.1 * length_km);
    EXPECT_NEAR(riselow, 3.5, 0.35);
    EXPECT_NEAR(risehigh, 3.5, 0.35);
    EXPECT_EQ(field.reversaledges, 0);
}

TEST_F(Landscape, RegionalFoldsOrientation)
{
    // S5.3c (round two: by continuity alone). A belt over every column, rows 56 to 71 (16 rows about the equator),
    // the strain rising from 0.2 at its edges by 0.15 per row to a plateau of 0.8, the axis wandering about north-south
    // (pi / 2 + 0.5 sin(4 pi i / W)): the flood from the first cell, pointed north, carries one sign over the whole
    // belt (every oriented axis points north, no reversal edge). The fabric of two overlapping windows (P, Q across
    // the date line) at k 16 is bit-identical on their overlap.
    RegionalGlobalInputs inputs = phaseworld();
    const float profile[16] = {0.2f, 0.35f, 0.5f, 0.65f, 0.8f, 0.8f, 0.8f, 0.8f, 0.8f, 0.8f, 0.8f, 0.8f, 0.65f, 0.5f, 0.35f, 0.2f};
    for (int j = 56; j < 72; j++)
        for (int i = 0; i < inputs.width; i++)
        {
            inputs.foldstrain[at(i, j, inputs.width)] = profile[j - 56];
            inputs.folddirection_rad[at(i, j, inputs.width)] =
                static_cast<float>(pi / 2.0 + 0.5 * std::sin(4.0 * pi * i / static_cast<double>(inputs.width)));
        }
    RegionalParameters parameters;
    parameters.folds.amplitudemax_m = 100.0;
    parameters.folds.wavelengthmin_km = 40.0;  // drawn in the single pass at k 16 (9.77 km cells: L / 4 = 10 km)
    parameters.folds.wavelengthmax_km = 40.0;
    const landscape::FoldPhaseField field = landscape::foldphasefield(inputs, parameters.folds);
    int wrong = 0;
    for (int j = 56; j < 72; j++)
        for (int i = 0; i < inputs.width; i++)
            wrong += field.wavenumbersouth[at(i, j, inputs.width)] < 0.0 ? 0 : 1;
    std::printf("[folds] orientation: oriented by the flood %d, north %d of %d; %d not pointing north; %d reversal edges\n",
                field.orientedbyflood, field.orientedbynorth, field.maskcells, wrong, field.reversaledges);
    EXPECT_EQ(field.orientedbyflood, field.maskcells - 1);
    EXPECT_EQ(field.orientedbynorth, 1);
    EXPECT_EQ(wrong, 0);
    EXPECT_EQ(field.reversaledges, 0);

    // Two windows over the belt, the second across the date line, sharing global columns 0 to 3 of rows 60 to 67.
    const RegionalWindow windowp = landscape::regionalwindow(256, 16, 1, -180.0, -5.625, -168.75, 5.625);
    const RegionalWindow windowq = landscape::regionalwindow(256, 16, 1, 174.375, -5.625, -174.375, 5.625);
    const auto fold = [&](const RegionalWindow& window) {
        const std::size_t count = static_cast<std::size_t>(window.domainfinewidth()) * static_cast<std::size_t>(window.domainfineheight());
        landscape::FineChannels channels;
        channels.channel.assign(count, 0);
        const std::vector<uint8_t> sea(count, 0);
        std::vector<float> elevation(count, 1000.0f);
        landscape::addfoldfabric(inputs, window, parameters, channels, sea, elevation);
        return elevation;
    };
    const std::vector<float> p = fold(windowp), q = fold(windowq);
    int shared = 0, differing = 0, folded = 0;
    for (int y = 0; y < windowp.domainfineheight(); y++)
        for (int x = 0; x < windowp.domainfinewidth(); x++)
        {
            const std::int64_t column = windowp.absolutefinecolumn(x), row = windowp.absolutefinerow(y);
            for (int yq = 0; yq < windowq.domainfineheight(); yq++)
            {
                if (windowq.absolutefinerow(yq) != row)
                    continue;
                for (int xq = 0; xq < windowq.domainfinewidth(); xq++)
                    if (windowq.absolutefinecolumn(xq) == column)
                    {
                        shared++;
                        const float a = p[at(x, y, windowp.domainfinewidth())], b = q[at(xq, yq, windowq.domainfinewidth())];
                        differing += std::memcmp(&a, &b, sizeof(float)) != 0 ? 1 : 0;
                        folded += a != 1000.0f ? 1 : 0;
                    }
            }
        }
    std::printf("[folds] orientation, windows P and Q at k 16: %d shared domain cells, %d folded, %d differ\n", shared, folded, differing);
    EXPECT_GT(shared, 0);
    EXPECT_GT(folded, 0);
    EXPECT_EQ(differing, 0);
}

TEST_F(Landscape, RegionalFoldsIceGate)
{
    // S5.3b: no fabric where the glaciation upsampled is above folds.icemax (0.5): a 3 x 3 group of
    // global cells at glaciation 1 leaves its cells above 0.5 unchanged; every other cell is folded as
    // without the field; icemax 1 is no gate.
    FoldBench plain(8.0, true);
    plain.fold();
    FoldBench iced(8.0, true);
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
    FoldBench plain(8.0, true);
    plain.fold();
    FoldBench lake(8.0, true);
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
    // 9.77 km): wavelengths of 30 x the cover in km clamped to 40 to 300 km, so level 1 draws L >= 156 km
    // (4 cells per wavelength) and level 2 the rest. Level 1 of two levels is the one-level result bit for bit; the k = 16 result
    // box-averaged equals both within the level-2 bound; windows P and Q (Q across the date line) agree
    // bit for bit on their overlap at both levels; both levels draw folds.
    const RegionalGlobalInputs inputs = withfoldfields(syntheticworld(256));
    RegionalParameters parameters;
    parameters.folds.amplitudemax_m = 500.0;
    parameters.folds.wavelengthmin_km = 40.0;
    parameters.folds.wavelengthmax_km = 300.0;
    parameters.folds.wavelengthpercover_km = 30.0;
    parameters.folds.wavelengthbase_km = 0.0;
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
