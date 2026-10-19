// Distributed deformation on the cube (W2-E, E3): the taper, the shares of the closing at
// convergent continental contacts, and the persistent subcells.
#include "gtest/gtest.h"
#include "cube_fixture.hpp"
#include "simulations/geology/plate_tectonics/tectonic_deformation.hpp"
#include "simulations/geology/plate_tectonics/tectonic_columns.hpp"
#include "core/mesh/cube_geometry.hpp"
#include "core/mesh/mesh_geometry.hpp"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

namespace {
using cubefixture::Band;
using cubefixture::CellId;
using cubefixture::Vec3;
using tectonics::columns::Model;
constexpr double pi = cubefixture::pi;

double along(const Band& band, CellId c, const Vec3& v) { return mesh::dot(v, cubefixture::east(*band.grid, c)); }
}

TEST(TectonicDeformation, DistanceCurveHasAnInlandStrainPeakAndCompactSupport) {
    using tectonics::deformation::taper;
    EXPECT_DOUBLE_EQ(taper(0,300000),1);
    EXPECT_DOUBLE_EQ(taper(600000,300000),0);
    EXPECT_DOUBLE_EQ(taper(900000,300000),0);
    const double shore=taper(0,300000)-taper(10000,300000);
    const double inland=taper(295000,300000)-taper(305000,300000);
    EXPECT_GT(inland,shore*10);
    EXPECT_NEAR(taper(150000,300000),taper(300000,600000),1e-15);
}

// Two continental halves (split by y = 0) on a 50 km band: sliding along the boundary, a
// disabled profile and the diverging contact leave the motion unchanged; the convergent
// contact spreads shortening into both plates.
TEST(TectonicDeformation, TransformDivergenceAndDisabledProfileLeaveMaterialMotionUnchanged) {
    using namespace tectonics;
    const Band band(16, cubefixture::radiusfor_km(16, 50000));
    const columns::Metric metric(*band.grid);
    const auto owners = cubefixture::hemispheres(*band.grid);
    const std::vector<double> continent(owners.size(), 1);
    const double turn = .1 * pi / 32;
    // Slip along the boundary: the halves turn about y, the boundary plane's normal.
    const auto slip = Model::displacement(*band.grid, cubefixture::plates(owners, {{0, turn, 0}, {0, -turn, 0}}));
    const auto shear = deformation::distribute(metric, owners, continent, slip, 1, {});
    for (size_t i = 0; i < slip.size(); ++i) EXPECT_LT(mesh::length(shear[i] - slip[i]), 1e-6 * band.spacing_m());
    const auto motion = Model::displacement(*band.grid, cubefixture::plates(owners, {band.eastward(.1), band.eastward(-.1)}));
    const auto disabled = deformation::distribute(metric, owners, continent, motion, 1, {0, .2, .2});
    for (size_t i = 0; i < motion.size(); ++i) EXPECT_EQ(mesh::length(disabled[i] - motion[i]), 0);
    // The contact at 180 degrees parts, the one at 0 closes; only the latter spreads shortening.
    const auto active = deformation::distribute(metric, owners, continent, motion, 1, {200, 0, .2});
    const CellId ridge = band.cell(56, 8), west = band.cell(23, 8), east = band.cell(24, 8);
    EXPECT_EQ(mesh::length(active[ridge] - motion[ridge]), 0);
    EXPECT_LT(along(band, west, active[west]), along(band, west, motion[west]));
    EXPECT_GT(along(band, east, active[east]), along(band, east, motion[east]));
}

// The width's seeded variation is repeatable, and a common translation added to both plates
// (a rotation of the whole sphere) changes nothing but itself.
TEST(TectonicDeformation, SeededWidthIsRepeatableAndDoesNotChangeWithACommonMotion) {
    using namespace tectonics;
    const Band band(16, cubefixture::radiusfor_km(16, 50000));
    const columns::Metric metric(*band.grid);
    const auto owners = cubefixture::hemispheres(*band.grid);
    const std::vector<double> continent(owners.size(), 1);
    const auto v = Model::displacement(*band.grid, cubefixture::plates(owners, {band.eastward(.1), band.eastward(-.1)}));
    const Vec3 common{.05 * pi / 32, .02 * pi / 32, .03 * pi / 32};
    const auto offset = Model::displacement(*band.grid, cubefixture::plates(owners,
        {band.eastward(.1) + common, band.eastward(-.1) + common}));
    const auto a = deformation::distribute(metric, owners, continent, v, 1, {});
    const auto b = deformation::distribute(metric, owners, continent, v, 1, {});
    const auto moved = deformation::distribute(metric, owners, continent, offset, 1, {});
    double along_strike = 0, largest = 0;
    for (size_t i = 0; i < v.size(); ++i) {
        EXPECT_EQ(mesh::length(a[i] - b[i]), 0);
        // The corrections are the same: they read the motion less its common rotation.
        EXPECT_LT(mesh::length((moved[i] - offset[i]) - (a[i] - v[i])), 2e-3 * band.spacing_m());
        largest = std::max(largest, mesh::length(a[i] - v[i]));
    }
    // Along the boundary the shortening varies with the seeded width.
    for (int y = 2; y < 14; ++y) {
        const CellId c = band.cell(22, y), d = band.cell(22, 8);
        along_strike = std::max(along_strike, std::abs(mesh::length(a[c] - v[c]) - mesh::length(a[d] - v[d])));
    }
    EXPECT_GT(largest, 0);
    EXPECT_GT(along_strike, 1e-4 * band.spacing_m());
}

// An ocean plate closing on a continent builds an inland belt: the shortening peaks a few
// cells inland of the contact, the elevation rises there above the rigid run's, and no
// continental mass is created.
TEST(TectonicDeformation, OceanConvergenceBuildsAnInlandBeltWithoutCreatingContinentalMass) {
    using namespace tectonics;
    const Band band(32, cubefixture::radiusfor_km(32, 50000));
    const columns::Metric metric(*band.grid);
    // The continent: band columns 32..95 (longitude -45..135), owner 1; ocean elsewhere, owner 0.
    std::vector<uint32_t> owners(band.cellcount(), 0);
    std::vector<double> continent(band.cellcount(), 0);
    std::vector<columns::Cell> initial(band.cellcount(), columns::seed(1000, false, 100));
    for (int y = 0; y < band.rows(); ++y) for (int x = 32; x < 96; ++x) {
        const CellId c = band.cell(x, y);
        owners[c] = 1; continent[c] = 1; initial[c] = columns::seed(1000, true, 100);
    }
    // Closing at 0.02 column per update from both sides (the offsets of a uniform east-west
    // motion on the band, zero on the polar faces).
    std::vector<Vec3> rigid(band.cellcount());
    for (int y = 0; y < band.rows(); ++y) for (int x = 0; x < band.columns(); ++x) {
        const CellId c = band.cell(x, y);
        rigid[c] = ((owners[c] ? -.02 : .02) * band.spacing_m()) * cubefixture::east(*band.grid, c);
    }
    const auto v = deformation::distribute(metric, owners, continent, rigid, 1, {300, 0, .2});
    const auto motion = [&](const std::vector<Vec3>& d) {
        std::vector<Vec3> offset = d;
        return cubefixture::offsets(*band.grid, offset, owners);
    };
    Model a, b; a.initialize(band.grid, initial); b.initialize(band.grid, initial);
    double mass = 0;
    for (CellId c = 0; c < initial.size(); ++c) mass += initial[c].continent * band.grid->area_m2(c);
    for (int t = 0; t < 100; ++t) { a.advance(motion(v), .25, 0, 1); b.advance(motion(rigid), .25, 0, 1); }
    int peak = 32; double final_mass = 0, control_peak = 0;
    for (CellId c = 0; c < initial.size(); ++c) final_mass += a.cells()[c].continent * band.grid->area_m2(c);
    const int row = band.rows() / 2;
    for (int x = 32; x < 64; ++x) {
        if (columns::elevation_m(a.cells()[band.cell(x, row)]) > columns::elevation_m(a.cells()[band.cell(peak, row)])) peak = x;
        control_peak = std::max(control_peak, columns::elevation_m(b.cells()[band.cell(x, row)]));
    }
    RecordProperty("peak_cell_from_initial_contact", peak - 32);
    EXPECT_GT(peak, 34); EXPECT_LT(peak, 46);
    EXPECT_GT(columns::elevation_m(a.cells()[band.cell(peak, row)]), 1100);
    EXPECT_LT(columns::elevation_m(a.cells()[band.cell(32, row)]), control_peak);
    EXPECT_NEAR(final_mass, mass, mass * 1e-12);
    EXPECT_LT(std::abs(a.budget().relative_mass_residual), 1e-12);
}

// Persistent subcells keep relief of four cells' wavelength through 200 updates of rigid
// rotation (without them the remap at the cells' scale smears it); age and range stay.
TEST(TectonicColumns, PersistentSubcellsRetainFourCellTerrainWithoutNoiseOrExtraMass) {
    using namespace tectonics::columns;
    const Band band(16, cubefixture::radiusfor_km(16, 1e7));  // wide cells: no root relaxation
    const auto initial = cubefixture::layout(band, seed(1000, true, 100), [](int x, int) {
        return seed(1000 + 100 * std::sin(2 * pi * x / 4), true, 100); });
    double retained[2]{};
    for (int mode = 0; mode < 2; ++mode) {
        Model model; model.initialize(band.grid, initial, mode ? 2 : 1);
        const auto motion = cubefixture::rotating(*band.grid, band.eastward(.16));
        for (int t = 0; t < 200; ++t) model.advance(motion, 1, 0, 0);
        double signal = 0, norm = 0;
        for (int y = 0; y < band.rows(); ++y) for (int x = 0; x < band.columns(); ++x) {
            const auto& c = model.cells()[band.cell(x, y)];
            const double wave = std::sin(2 * pi * (x - 200 * .16) / 4);
            signal += (elevation_m(c) - 1000) * wave;
            norm += 100 * wave * wave;
            EXPECT_NEAR(age_myr(c), 300, 1e-8);
            EXPECT_GE(elevation_m(c), 899.9);
            EXPECT_LE(elevation_m(c), 1100.1);
        }
        retained[mode] = signal / norm;
        EXPECT_LT(std::abs(model.budget().relative_mass_residual), 1e-11);
    }
    RecordProperty("native_signal_retained", retained[0]); RecordProperty("subcell_signal_retained", retained[1]);
    EXPECT_LT(retained[0], .2); EXPECT_GT(retained[1], retained[0] + .2);
}

// Distributed deformation on the cube (W2-E, E2B): an ocean plate on face 0 closes on a
// continent beyond its east edge; the upper plate's share of the motion keeps its direction
// on the far side of the face edge (each cell's own frame).
TEST(TectonicDeformationCube, ContactCorrectionsKeepTheirDirectionAcrossAFaceEdge) {
    constexpr int n = 32;
    constexpr double earth_km = 6371.0;
    const auto grid = tectonics::grid::TectonicGrid::cube(n, earth_km);
    const tectonics::columns::Metric metric(*grid);
    const auto& mesh = *grid->mesh();
    const size_t count = grid->cellcount();
    std::vector<uint32_t> owners(count, 1);
    std::vector<double> continent(count, 1);
    std::vector<mesh::Vec3> rigid(count);
    const mesh::Vec3 right = mesh::cubeface(0).right;
    const auto alongright = [&](CellId c) { const mesh::Vec3 r = mesh.centre(c); return mesh::normalise(right - mesh::dot(right, r) * r); };
    for (CellId c = 0; c < count; ++c) {
        int face = 0, column = 0, row = 0;
        mesh::cubeindices(n, c, face, column, row);
        if (face != 0 || column >= n - 1) continue;
        owners[c] = 0; continent[c] = 0;
        rigid[c] = 30000.0 * alongright(c);
    }
    tectonics::deformation::Parameters p;
    p.peak_distance_km = 500;
    p.distance_variation = 0;
    const auto result = tectonics::deformation::distribute(metric, owners, continent, rigid, 1, p);
    int east_face = 0, east_column = 0, east_row = 0;
    mesh::cubeindices(n, grid->neighbour(mesh::cubeid(n, 0, n - 1, n / 2), tectonics::grid::east), east_face, east_column, east_row);
    int crossed = 0;
    for (CellId c = 0; c < count; ++c) {
        ASSERT_TRUE(std::isfinite(result[c].x) && std::isfinite(result[c].y) && std::isfinite(result[c].z));
        if (owners[c] != 1) continue;
        const mesh::Vec3 change = result[c] - rigid[c];
        EXPECT_NEAR(mesh::dot(result[c], mesh.centre(c)), 0, 1e-6 * (1 + mesh::length(result[c])));
        if (mesh::length(change) == 0) continue;
        int face = 0, column = 0, row = 0;
        mesh::cubeindices(n, c, face, column, row);
        if (face != east_face || std::abs(row - n / 2) > n / 4) continue;
        // Beyond the edge, the push points along the ocean plate's motion (into the continent).
        EXPECT_GT(mesh::dot(change, alongright(c)), .9 * mesh::length(change));
        ++crossed;
    }
    EXPECT_GT(crossed, 0);
}
