// Tests of landscape/regional_cirques.hpp (erosion plan Phase E5: the E5.0 contract, which
// step E5.2 implements: cirques, the parabolic bowls on glaciated high ground), suite
// Landscape: TEST_F(Landscape, RegionalCirques<What>) on the synthetic world of
// regional_test_world.hpp with withglacialfields, and the overlap variants with the
// cirques on (RegionalOverlapCirques*), on a world wider than the widened box.

#include "landscape_fixture.hpp"

#include "regional_test_world.hpp"

#include "simulations/landscape/regional_cirques.hpp"
#include "simulations/landscape/regional_glacial.hpp"

#include <set>
#include <utility>

namespace
{
// One call of the hook on a window's domain with hand-built fine fields: a flat surface,
// the run's ELA constant, no channel, no sea, every global land cell glaciated (1).
struct CirqueCase
{
    RegionalGlobalInputs inputs;
    RegionalWindow window;
    RegionalParameters parameters;
    landscape::FineGlacial glacial;
    landscape::FineChannels channels;
    std::vector<uint8_t> sea;
    std::vector<float> elevation;

    std::size_t count() const { return elevation.size(); }
    int width() const { return window.domainfinewidth(); }
    double cell_m() const { return landscape::finepatchmetric(window, inputs.radius_m).dy_m; }
    // The radius in kilometres that is `cells` fine cells on the ground.
    void radiuscells(double cells) { parameters.cirques.radius_km = cells * cell_m() / 1000.0; }
    // A flat surface: the fine elevation and the global base elevation (the candidate's
    // band test and floor read the base upsampled at the centre, the review of E5.0).
    void surface(float value)
    {
        elevation.assign(count(), value);
        inputs.elevation_m.assign(inputs.elevation_m.size(), value);
    }
    std::vector<float> carve() const
    {
        std::vector<float> out = elevation;
        landscape::carvecirques(inputs, window, parameters, glacial, channels, sea, out);
        return out;
    }
};

CirqueCase cirquecase(const RegionalWindow& window, float surface_m, float ela_m)
{
    CirqueCase c;
    c.inputs = syntheticworld(window.worldwidth);
    const std::size_t size = static_cast<std::size_t>(c.inputs.width) * static_cast<std::size_t>(c.inputs.height);
    c.inputs.glaciation.assign(size, 0.0f);
    for (std::size_t cell = 0; cell < size; cell++)
        c.inputs.glaciation[cell] = c.inputs.network.sea[cell] != 0 ? 0.0f : 1.0f;
    c.window = window;
    const std::size_t count = static_cast<std::size_t>(window.domainfinewidth()) * static_cast<std::size_t>(window.domainfineheight());
    c.glacial.active = true;
    c.glacial.runela_m.assign(count, ela_m);
    c.glacial.presentela_m.assign(count, ela_m + 500.0f);
    c.glacial.glaciation.assign(count, 1.0f);
    c.channels.channel.assign(count, 0);
    c.sea.assign(count, 0);
    c.elevation.assign(count, surface_m);
    c.inputs.elevation_m.assign(size, surface_m);  // the base the candidates read; network.sea stays the synthetic world's
    c.parameters.cirques.radius_km = 0.0;
    c.parameters.cirques.depth_m = 300.0;
    c.parameters.cirques.band_m = 1500.0;
    c.parameters.cirques.candidates = 2;
    c.parameters.cirques.glaciationmin = 0.5;
    return c;
}

// The candidates' fine cells: on a flat surface 3000 m over an ELA of 2000 m every
// candidate is a cirque, and with a radius under one cell each lowers its centre alone
// (to 2700 m).
std::vector<std::size_t> centres(CirqueCase c)
{
    c.surface(3000.0f);
    c.glacial.runela_m.assign(c.count(), 2000.0f);
    c.radiuscells(0.25);
    const std::vector<float> out = c.carve();
    std::vector<std::size_t> found;
    for (std::size_t cell = 0; cell < out.size(); cell++)
        if (out[cell] != 3000.0f)
        {
            EXPECT_EQ(out[cell], 2700.0f);
            found.push_back(cell);
        }
    return found;
}

int changed(const std::vector<float>& a, const std::vector<float>& b)
{
    int count = 0;
    for (std::size_t cell = 0; cell < a.size(); cell++)
        count += a[cell] != b[cell] ? 1 : 0;
    return count;
}

int changedcirques(const RegionalFields& a, const RegionalFields& b) { return changed(a.elevation_m, b.elevation_m); }

// The changed kept cells of `a`'s window inside the global cells `other`'s kept box
// shares with it: the overlap test shows a bowl reached the shared cells (E5.2's review).
int changedonoverlap(const RegionalFields& a, const RegionalFields& b, const RegionalWindow& other)
{
    const RegionalWindow& window = a.window;
    const int k = window.scale, width = window.finewidth();
    int count = 0;
    for (std::size_t cell = 0; cell < a.elevation_m.size(); cell++)
    {
        if (a.elevation_m[cell] == b.elevation_m[cell])
            continue;
        const int x = static_cast<int>(cell % static_cast<std::size_t>(width));
        const int y = static_cast<int>(cell / static_cast<std::size_t>(width));
        const int i = (window.column0 + x / k) % window.worldwidth, j = window.row0 + y / k;
        const int column = ((i - other.column0) % window.worldwidth + window.worldwidth) % window.worldwidth;
        if (column < other.columns && j >= other.row0 && j < other.row0 + other.rows)
            count++;
    }
    return count;
}

bool indomain(const RegionalWindow& window, int i, int j)
{
    const int column = ((i - window.domaincolumn0) % window.worldwidth + window.worldwidth) % window.worldwidth;
    return column < window.domaincolumns && j >= window.domainrow0 && j < window.domainrow0 + window.domainrows;
}

// The parameters of the overlap variants: bowls of about 2 to 4 fine cells.
RegionalParameters cirqueparameters()
{
    RegionalParameters parameters;
    parameters.cirques.radius_km = 600.0;  // 1.9 cells at k 4 and 2.4 at k 5 on W 32, 3.8 at k 4 on W 64
    return parameters;
}
}

TEST_F(Landscape, RegionalCirquesAreOffByDefault)
{
    // At the defaults (radius_km 0) the hook changes nothing, with or without the fields.
    const RegionalGlobalInputs inputs = withglacialfields(syntheticworld(32));
    const RegionalWindow window = windowa(4);
    const std::size_t count = static_cast<std::size_t>(window.domainfinewidth()) * static_cast<std::size_t>(window.domainfineheight());
    RegionalParameters parameters;
    parameters.troughs.width_km = 10.0;  // another feature on: the fields are active, the cirques still off
    const landscape::FineGlacial glacial = landscape::regionalglacialfields(inputs, window, parameters);
    ASSERT_TRUE(glacial.active);
    landscape::FineChannels channels;
    channels.channel.assign(count, 0);
    const std::vector<uint8_t> sea(count, 0);
    std::vector<float> elevation(count, 3000.0f);
    landscape::carvecirques(inputs, window, parameters, glacial, channels, sea, elevation);
    for (std::size_t cell = 0; cell < count; cell++)
        ASSERT_EQ(elevation[cell], 3000.0f);
}

TEST_F(Landscape, RegionalCirquesInertWithoutFields)
{
    // The radius on but the glacial fields inactive (no ELA field): nothing changes.
    CirqueCase c = cirquecase(windowb(4), 3000.0f, 2000.0f);
    c.radiuscells(2.5);
    ASSERT_GT(changed(c.carve(), c.elevation), 0);
    c.glacial.active = false;
    EXPECT_EQ(changed(c.carve(), c.elevation), 0);
}

TEST_F(Landscape, RegionalCirquesPlacement)
{
    // The candidates are a pure function of the global cell, the index and the seed: two
    // windows whose domains share global cells place the same fine cells there; at most
    // `candidates` per block, only in glaciated land blocks, and a third candidate adds
    // cirques without moving the first two.
    for (int scale : {4, 5})
    {
        SCOPED_TRACE("k = " + std::to_string(scale));
        const CirqueCase a = cirquecase(windowa(scale), 3000.0f, 2000.0f);
        const CirqueCase b = cirquecase(windowb(scale), 3000.0f, 2000.0f);
        std::set<std::pair<std::int64_t, std::int64_t>> placeda, placedb;
        std::vector<int> perblock(static_cast<std::size_t>(a.inputs.width) * static_cast<std::size_t>(a.inputs.height), 0);
        for (const CirqueCase* c : {&a, &b})
        {
            auto& placed = c == &a ? placeda : placedb;
            for (std::size_t cell : centres(*c))
            {
                const int x = static_cast<int>(cell % static_cast<std::size_t>(c->width()));
                const int y = static_cast<int>(cell / static_cast<std::size_t>(c->width()));
                const int i = c->window.globalcolumn(x), j = c->window.globalrow(y);
                const std::size_t global = at(i, j, c->inputs.width);
                EXPECT_EQ(c->inputs.network.sea[global], 0);
                if (!indomain(a.window, i, j) || !indomain(b.window, i, j))
                    continue;
                placed.insert({c->window.absolutefinecolumn(x), c->window.absolutefinerow(y)});
                if (c == &a)
                    perblock[global]++;
            }
        }
        EXPECT_GT(placeda.size(), 50u);
        EXPECT_EQ(placeda, placedb);
        for (int value : perblock)
            EXPECT_LE(value, 2);

        CirqueCase three = b;
        three.parameters.cirques.candidates = 3;
        const std::vector<std::size_t> two = centres(b), more = centres(three);
        const std::set<std::size_t> moreset(more.begin(), more.end());
        for (std::size_t cell : two)
            EXPECT_EQ(moreset.count(cell), 1u);
        EXPECT_GT(more.size(), two.size());
    }
}

TEST_F(Landscape, RegionalCirquesBandAndEla)
{
    // A flat surface, centres alone (radius a quarter cell), the run's ELA 3500 m, the band
    // 1500 m, the depth 300 m.
    CirqueCase c = cirquecase(windowb(4), 0.0f, 3500.0f);
    const std::vector<std::size_t> found = centres(c);
    ASSERT_GT(found.size(), 50u);
    c.radiuscells(0.25);

    // Just under the band's top (4999 m; the upsampled base of a constant field may differ
    // from it by rounding, so the edge itself is not asserted): a cirque, its floor 300 m
    // down (4699 m).
    c.surface(4999.0f);
    std::vector<float> out = c.carve();
    EXPECT_EQ(changed(out, c.elevation), static_cast<int>(found.size()));
    for (std::size_t cell : found)
        EXPECT_EQ(out[cell], 4699.0f);

    // Above the band (5002 m): plateau ice, no cirque.
    c.surface(5002.0f);
    EXPECT_EQ(changed(c.carve(), c.elevation), 0);

    // 100 m above the ELA (3600 m): the floor stops at the ELA (3500 m), not 300 m down.
    c.surface(3600.0f);
    out = c.carve();
    EXPECT_EQ(changed(out, c.elevation), static_cast<int>(found.size()));
    for (std::size_t cell : found)
        EXPECT_EQ(out[cell], 3500.0f);

    // Below the ELA: with a bowl of 2.5 cells, a base 2 m under the ELA (3498 m) under fine
    // ground at 4000 m cuts nothing, although the bowls would lower the ground around the
    // centres; a base 2 m over it (3502 m) does, and the floor is the ELA (max(3500, 3202)),
    // which the centre cells take (d = 0).
    c.radiuscells(2.5);
    c.elevation.assign(c.count(), 4000.0f);
    c.inputs.elevation_m.assign(c.inputs.elevation_m.size(), 3498.0f);
    EXPECT_EQ(changed(c.carve(), c.elevation), 0);
    c.inputs.elevation_m.assign(c.inputs.elevation_m.size(), 3502.0f);
    out = c.carve();
    EXPECT_GT(changed(out, c.elevation), static_cast<int>(found.size()));
    for (std::size_t cell : found)
        EXPECT_EQ(out[cell], 3500.0f);
}

TEST_F(Landscape, RegionalCirquesBowlProfile)
{
    // One cirque: only global cell (2, 8) glaciated, one candidate; a flat surface 3000 m
    // over an ELA of 2000 m, depth 300 m, so the floor is 2700 m; radius R = 2.5 fine cells:
    // target 2700 + 300 (d / 2.5)^2: 2700 at d 0, 2748 at d 1, 2796 at d sqrt 2, 2892 at
    // d 2, 2940 at d sqrt 5; 3000 (kept) beyond 2.5.
    CirqueCase c = cirquecase(windowb(4), 3000.0f, 2000.0f);
    const std::size_t global = at(2, 8, c.inputs.width);
    ASSERT_EQ(c.inputs.network.sea[global], 0);
    std::fill(c.inputs.glaciation.begin(), c.inputs.glaciation.end(), 0.0f);
    c.inputs.glaciation[global] = 1.0f;
    c.parameters.cirques.candidates = 1;
    c.radiuscells(2.5);
    const std::vector<float> out = c.carve();
    const std::size_t centre = static_cast<std::size_t>(std::min_element(out.begin(), out.end()) - out.begin());
    const int cx = static_cast<int>(centre % static_cast<std::size_t>(c.width()));
    const int cy = static_cast<int>(centre / static_cast<std::size_t>(c.width()));
    EXPECT_EQ(c.window.globalcolumn(cx), 2);
    EXPECT_EQ(c.window.globalrow(cy), 8);
    int lowered = 0;
    for (int dy = -3; dy <= 3; dy++)
        for (int dx = -3; dx <= 3; dx++)
        {
            const int d2 = dx * dx + dy * dy;
            const float value = out[at(cx + dx, cy + dy, c.width())];
            const double expected = d2 <= 6 ? 2700.0 + 300.0 * d2 / 6.25 : 3000.0;
            EXPECT_NEAR(value, expected, 1.0e-3) << "dx " << dx << " dy " << dy;
            lowered += value < 3000.0f ? 1 : 0;
        }
    EXPECT_EQ(lowered, 21);  // the cells with d^2 <= 6: 1 + 4 + 4 + 4 + 8
    EXPECT_EQ(changed(out, c.elevation), 21);
}

TEST_F(Landscape, RegionalCirquesKeepChannelsAndSea)
{
    // Bowls of 2.5 cells on a flat surface 3000 m over an ELA of 2000 m: channel and sea
    // cells are never lowered, and a candidate on a channel or sea cell is not a cirque.
    CirqueCase c = cirquecase(windowb(4), 3000.0f, 2000.0f);
    const std::vector<std::size_t> found = centres(c);
    c.radiuscells(2.5);
    for (int y = 0; y < c.window.domainfineheight(); y++)
        for (int x = 0; x < c.width(); x++)
        {
            if (x % 3 == 0)
                c.channels.channel[at(x, y, c.width())] = 1;
            if (y % 5 == 0)
                c.sea[at(x, y, c.width())] = 1;
        }
    const std::vector<float> out = c.carve();
    EXPECT_GT(changed(out, c.elevation), 0);
    for (std::size_t cell = 0; cell < c.count(); cell++)
        if (c.channels.channel[cell] != 0 || c.sea[cell] != 0)
            EXPECT_EQ(out[cell], 3000.0f);

    CirqueCase onchannels = cirquecase(windowb(4), 3000.0f, 2000.0f);
    onchannels.radiuscells(2.5);
    for (std::size_t cell : found)
        onchannels.channels.channel[cell] = 1;
    EXPECT_EQ(changed(onchannels.carve(), onchannels.elevation), 0);
    CirqueCase onsea = cirquecase(windowb(4), 3000.0f, 2000.0f);
    onsea.radiuscells(2.5);
    for (std::size_t cell : found)
        onsea.sea[cell] = 1;
    EXPECT_EQ(changed(onsea.carve(), onsea.elevation), 0);
}

TEST_F(Landscape, RegionalCirquesLowerOnly)
{
    // An uneven surface (1500 to 4500 m) over an ELA of 2000 m, bowls of 3 cells, three
    // candidates: no cell is raised, some are lowered, none to or below the sea level.
    CirqueCase c = cirquecase(windowb(4), 3000.0f, 2000.0f);  // the base 3000 m: every candidate in the band
    c.parameters.cirques.candidates = 3;
    c.radiuscells(3.0);
    for (int y = 0; y < c.window.domainfineheight(); y++)
        for (int x = 0; x < c.width(); x++)
            c.elevation[at(x, y, c.width())] =
                static_cast<float>(3000.0 + 1000.0 * std::sin(0.7 * x) * std::cos(0.5 * y) + 500.0 * std::sin(0.13 * (x + 2 * y)));
    const std::vector<float> out = c.carve();
    int lowered = 0;
    for (std::size_t cell = 0; cell < c.count(); cell++)
    {
        EXPECT_LE(out[cell], c.elevation[cell]);
        EXPECT_GT(out[cell], c.inputs.sealevel);
        lowered += out[cell] < c.elevation[cell] ? 1 : 0;
    }
    EXPECT_GT(lowered, 100);
}

TEST_F(Landscape, RegionalCirquesSeaLevelFloor)
{
    // The sea level at 2900 m, a flat surface 3000 m over an ELA of 2000 m, depth 300 m:
    // the floor (2700 m) would lie under the sea, so it is the next float above the sea
    // level; nothing reaches the sea level.
    CirqueCase c = cirquecase(windowb(4), 3000.0f, 2000.0f);
    const std::vector<std::size_t> found = centres(c);
    c.inputs.sealevel = 2900.0f;
    c.radiuscells(2.5);
    const std::vector<float> out = c.carve();
    const float floorlevel = std::nextafter(2900.0f, 1.0e9f);
    for (std::size_t cell : found)
        EXPECT_EQ(out[cell], floorlevel);
    for (float value : out)
        EXPECT_GE(value, floorlevel);
}

TEST_F(Landscape, RegionalOverlapCirques)
{
    // The cirques at the smallest halo the pass allows: the overlaps stay bit-identical
    // (compareoverlap asserts it), and the cirques act.
    const RegionalParameters parameters = cirqueparameters();
    const RegionalGlobalInputs columns = withglacialfields(syntheticworld(32));
    for (int scale : {4, 5})
    {
        SCOPED_TRACE("k = " + std::to_string(scale));
        const OverlapCounts counts = compareoverlap(columns, windowa(scale), windowb(scale), "cirques A/B", parameters);
        EXPECT_EQ(counts.shared, 4 * 5 * scale * scale);
        const RegionalFields on = landscape::runregionalpass(columns, windowb(scale), parameters);
        const RegionalFields off = landscape::runregionalpass(columns, windowb(scale), RegionalParameters{});
        const int cut = changedcirques(on, off);
        const int onoverlap = changedonoverlap(on, off, windowa(scale));
        std::printf("[regional] cirques B, k %d: %d cells lowered, %d on the overlap with A\n", scale, cut, onoverlap);
        EXPECT_GT(cut, 0);
        EXPECT_GT(onoverlap, 0);
    }
    const RegionalGlobalInputs rows = withglacialfields(syntheticworld(64));
    const OverlapCounts counts = compareoverlap(rows, windowc(4), windowd(4), "cirques C/D", parameters);
    EXPECT_GT(counts.sea, 0);
    EXPECT_GT(counts.channel, 0);
    const RegionalFields on = landscape::runregionalpass(rows, windowd(4), parameters);
    const RegionalFields off = landscape::runregionalpass(rows, windowd(4), RegionalParameters{});
    const int cut = changedcirques(on, off);
    const int onoverlap = changedonoverlap(on, off, windowc(4));
    std::printf("[regional] cirques D, k 4: %d cells lowered, %d on the overlap with C\n", cut, onoverlap);
    EXPECT_GT(cut, 0);
    EXPECT_GT(onoverlap, 0);
}
