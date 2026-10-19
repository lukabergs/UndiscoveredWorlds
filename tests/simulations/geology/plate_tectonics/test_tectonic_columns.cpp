// The isostatic columns on the cube (W2-E, E3): the column laws, the column model's processes
// and its transport, the rock uplift rate and the strain, on the equatorial band of the cube
// (cube_fixture.hpp) or the whole sphere; the lithosphere runs on small cubes.
#include "gtest/gtest.h"
#include "cube_fixture.hpp"
#include "simulations/geology/plate_tectonics/tectonic_columns.hpp"
#include "simulations/geology/plate_tectonics/tectonic_initialization.hpp"
#include "simulations/geology/plate_tectonics/tectonic_margins.hpp"
#include "core/mesh/cube_geometry.hpp"
#include "core/mesh/mesh_geometry.hpp"
#include "simulations/geology/plate_tectonics/tectonic_thermal.hpp"
#include "simulations/geology/plate_tectonics/lithosphere.hpp"
#include "validation/tectonics/tectonic_snapshots.hpp"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <functional>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

using namespace tectonics::columns;
using cubefixture::Band;
using cubefixture::CellId;
using cubefixture::Vec3;
using tectonics::grid::TectonicGrid;
namespace st = tectonics::sphere_transport;

namespace {
constexpr double pi = 3.14159265358979323846;
constexpr double earth_km = cubefixture::earth_km;

// Two plates split by the meridian plane at longitude 45 / -135 degrees: plate 0 west of
// 45 degrees (band columns 0 .. 2N - 1), plate 1 east of it.
std::vector<uint32_t> halves(const TectonicGrid& grid) {
    const Vec3 normal{std::cos(-pi / 4), std::sin(-pi / 4), 0};
    std::vector<uint32_t> owners(grid.cellcount());
    for (CellId c = 0; c < grid.cellcount(); ++c) owners[c] = mesh::dot(grid.position(c), normal) > 0 ? 0 : 1;
    return owners;
}
// The two halves turning about the z axis by +-cells band columns per update: they converge at
// 45 degrees and part at -135 degrees.
st::Motion converging(const Band& band, double cells) {
    return cubefixture::plates(halves(*band.grid), {band.eastward(cells), band.eastward(-cells)});
}

double area_weighted(const Model& model, double Cell::* member) {
    double sum = 0;
    for (CellId c = 0; c < model.cells().size(); ++c) sum += model.cells()[c].*member * model.metric().area_m2(c);
    return sum;
}

// The budget's continental area terms: stretched + arc + continentalized - collided - oceanized, m².
double continental_area_terms(const Budget& b) {
    return b.stretched_continental_area_m2 + b.arc_continental_area_m2 + b.continentalized_oceanic_area_m2
        - b.collided_continental_area_m2 - b.oceanized_continental_area_m2;
}

std::vector<Displacement> zero(const Model& model) { return std::vector<Displacement>(model.cells().size()); }
} // namespace

TEST(TectonicColumns, AiryRootsAndWaterLoadingHaveAnalyticalElevations) {
    for (double z : {-5000., -500., 0., 1500., 8000., 12000.}) {
        const auto c = seed(z, true, 100);
        EXPECT_NEAR(elevation_m(c), z, 1e-10);
        EXPECT_GT(c.continent, z);
    }
    const auto ridge = seed(-100, false, 0);
    EXPECT_DOUBLE_EQ(ridge.ocean, 7000);
    EXPECT_DOUBLE_EQ(elevation_m(ridge), -2500);
    EXPECT_NEAR(elevation_m(seed(-100, false, 80)), -2500 - tectonics::thermal::subsidence_m(80), 1e-10);
}

// The lithospheric root (continental_root_kg_m3): seed and elevation stay inverses,
// the slope per km of crust is (500 - root) / 3300 above sea and / 2270 below, and a
// root of 0 reproduces the crust-only expressions bit for bit.
TEST(TectonicColumns, LithosphericRootLowersTheSlopeAndZeroKeepsTheCrustOnlyColumn) {
    constexpr double root = 246;
    // -3500 m, not -5000: at this root continental crust thins to zero at -3.9 km.
    for (double z : {-3500., -500., 0., 1500., 8000., 12000.}) {
        const auto c = seed(z, true, 100, root);
        EXPECT_NEAR(elevation_m(c, 1, root), z, 1e-10);
    }
    Cell thick, thin;
    thick.continental_area = thin.continental_area = 1;
    thick.continent = 36000; thin.continent = 34000;
    EXPECT_NEAR(elevation_m(thick, 1, root), 77.0, .05);
    EXPECT_NEAR(elevation_m(thin, 1, root), -111.9, .05);
    EXPECT_NEAR(elevation_m(thick, 1, root), 254000. / 3300, 1e-9);
    EXPECT_NEAR(elevation_m(thin, 1, root), -254000. / 2270, 1e-9);
    // The crust-only expressions as they stood before the root was added.
    constexpr double ocean_reference = (mantle_density - oceanic_density) * ridge_crust_m
        + (mantle_density - water_density) * ridge_depth_m;
    constexpr double continental_mantle_load =
        (mantle_density - continental_density) * reference_continent_m - ocean_reference;
    const auto legacy_elevation = [&](const Cell& c, double thermal_strength) {
        const double ocean_age = c.ocean > 0 ? std::max(0.0,
            (c.oceanic_age_volume-c.thermal_rejuvenation_volume)/c.ocean) : 0;
        const double load = ocean_reference + continental_mantle_load * c.continental_area
            + (mantle_density - water_density) * thermal_strength
                * tectonics::thermal::subsidence_m(ocean_age) * c.oceanic_area;
        const double buoyancy = (mantle_density - continental_density) * c.continent
            + (mantle_density - oceanic_density) * c.ocean
            + (mantle_density - sediment_density) * c.sediment - load;
        return buoyancy / (buoyancy >= 0 ? mantle_density : mantle_density - water_density) - c.flexure_m;
    };
    Cell mixed;
    mixed.continent = 0.3 * 41234.5; mixed.ocean = 0.7 * 6800; mixed.sediment = 812.25;
    mixed.continental_area = .3; mixed.oceanic_area = .7;
    mixed.continental_age_volume = mixed.continent * 310; mixed.oceanic_age_volume = mixed.ocean * 47.5;
    mixed.thermal_rejuvenation_volume = mixed.ocean * 3; mixed.flexure_m = 123.4;
    Cell ghost; // Continental traces without continental area.
    ghost.continent = 1e-3; ghost.ocean = 7100; ghost.sediment = 20;
    ghost.oceanic_area = 1; ghost.oceanic_age_volume = ghost.ocean * 80;
    for (const auto& c : {mixed, ghost}) for (double thermal : {0., .5, 1.}) {
        EXPECT_EQ(elevation_m(c, thermal, 0), legacy_elevation(c, thermal));
        EXPECT_EQ(elevation_m(c, thermal, 0), elevation_m(c, thermal));
    }
    for (double z : {-5000., -500., 0., 1500., 12000.})
        EXPECT_EQ(seed(z, true, 100, 0).continent,
                  reference_continent_m + z * (z >= 0 ? mantle_density : mantle_density - water_density)
                      / (mantle_density - continental_density));
    EXPECT_THROW(seed(0, true, 100, 500), std::runtime_error);
    EXPECT_THROW(seed(0, true, 100, -1), std::runtime_error);
    const auto grid = cubefixture::cube(2);
    Model model;
    EXPECT_THROW(model.initialize(grid, std::vector<Cell>(24, seed(0, true, 100)), 1, 0, {}, false, 500),
                 std::runtime_error);
    EXPECT_THROW(model.initialize(grid, std::vector<Cell>(23, seed(0, true, 100))), std::runtime_error);
    EXPECT_THROW(model.initialize(nullptr, {}), std::runtime_error);
}

// The root's delamination cap (continental_root_max_thickness_m): above the cap the
// root stops growing and the crust-only slope 500 / 3300 returns; seed inverts it
// piecewise; cap 0 is the uncapped expression bit for bit. With a delamination width
// the root's load falls to zero over the width instead; width 0 is the frozen cap.
TEST(TectonicColumns, RootCapRestoresTheCrustOnlySlopeAboveIt) {
    constexpr double root = 246, cap = 50000;
    const double cap_z = (cap - reference_continent_m) * (500 - root) / mantle_density;
    EXPECT_NEAR(cap_z, 1155, .5);
    for (double z : {-3500., 0., 500., 1155., cap_z, 2000., 5000., 8000.}) {
        const auto c = seed(z, true, 100, root, cap);
        EXPECT_NEAR(elevation_m(c, 1, root, cap), z, 1e-9) << z;
    }
    EXPECT_NEAR(seed(cap_z, true, 100, root, cap).continent, cap, 1e-9);
    EXPECT_NEAR(seed(8000, true, 100, root, cap).continent, cap + (8000 - cap_z) * 3300 / 500, 1e-9);
    const auto at = [&](double continent) {
        Cell c; c.continental_area = 1; c.continent = continent;
        return elevation_m(c, 1, root, cap);
    };
    EXPECT_NEAR((at(41000) - at(40000)) / 1000, .077, .0005);
    EXPECT_NEAR((at(61000) - at(60000)) / 1000, .1515, .0005);
    EXPECT_NEAR((at(41000) - at(40000)) / 1000, 254. / 3300, 1e-12);
    EXPECT_NEAR((at(61000) - at(60000)) / 1000, 500. / 3300, 1e-12);
    // Mixed cells: the cap is area-weighted like the reference.
    Cell half; half.continental_area = .5; half.oceanic_area = .5;
    half.continent = 30000; half.ocean = 3500; half.oceanic_age_volume = half.ocean * 20;
    // The root acts on `rooted` metres of crust (the whole continent when uncapped).
    const auto expression = [](const Cell& c, double thermal, double r, double rooted) {
        constexpr double ocean_reference = (mantle_density - oceanic_density) * ridge_crust_m
            + (mantle_density - water_density) * ridge_depth_m;
        constexpr double continental_mantle_load =
            (mantle_density - continental_density) * reference_continent_m - ocean_reference;
        const double ocean_age = c.ocean > 0 ? std::max(0.0,
            (c.oceanic_age_volume-c.thermal_rejuvenation_volume)/c.ocean) : 0;
        const double load = ocean_reference + continental_mantle_load * c.continental_area
            + (mantle_density - water_density) * thermal
                * tectonics::thermal::subsidence_m(ocean_age) * c.oceanic_area;
        double buoyancy = (mantle_density - continental_density) * c.continent
            + (mantle_density - oceanic_density) * c.ocean
            + (mantle_density - sediment_density) * c.sediment - load;
        buoyancy -= r * (rooted - reference_continent_m * c.continental_area);
        return buoyancy / (buoyancy >= 0 ? mantle_density : mantle_density - water_density) - c.flexure_m;
    };
    const auto uncapped = [&](const Cell& c, double thermal, double r) { return expression(c, thermal, r, c.continent); };
    EXPECT_EQ(elevation_m(half, 1, root, cap), expression(half, 1, root, cap * half.continental_area));
    EXPECT_GT(elevation_m(half, 1, root, cap), elevation_m(half, 1, root));
    // Cap 0: the uncapped expression, at root 246 and at root 0.
    Cell mixed;
    mixed.continent = 0.3 * 61234.5; mixed.ocean = 0.7 * 6800; mixed.sediment = 812.25;
    mixed.continental_area = .3; mixed.oceanic_area = .7;
    mixed.continental_age_volume = mixed.continent * 310; mixed.oceanic_age_volume = mixed.ocean * 47.5;
    mixed.thermal_rejuvenation_volume = mixed.ocean * 3; mixed.flexure_m = 123.4;
    Cell thick; thick.continental_area = 1; thick.continent = 70000;
    for (const auto& c : {mixed, thick, half}) for (double r : {0., root}) for (double thermal : {0., 1.}) {
        EXPECT_EQ(elevation_m(c, thermal, r, 0), uncapped(c, thermal, r));
        EXPECT_EQ(elevation_m(c, thermal, r, 0), elevation_m(c, thermal, r));
    }
    for (double r : {0., root}) for (double z : {-3500., 0., 1500., 8000.})
        EXPECT_EQ(seed(z, true, 100, r, 0).continent,
                  reference_continent_m + z * (z >= 0 ? mantle_density : mantle_density - water_density)
                      / (mantle_density - continental_density - r));
    // Validation: 0 or above 35 km.
    EXPECT_NO_THROW(seed(0, true, 100, root, 0));
    EXPECT_THROW(seed(0, true, 100, root, 30000), std::runtime_error);
    EXPECT_THROW(seed(0, true, 100, root, reference_continent_m), std::runtime_error);
    EXPECT_THROW(seed(0, true, 100, root, std::numeric_limits<double>::quiet_NaN()), std::runtime_error);
    const auto grid = cubefixture::cube(2);
    Model model;
    const std::vector<Cell> flat(grid->cellcount(), seed(0, true, 100));
    EXPECT_NO_THROW(model.initialize(grid, flat, 1, 0, {}, false, root, 0));
    EXPECT_NO_THROW(model.initialize(grid, flat, 1, 0, {}, false, root, cap));
    EXPECT_THROW(model.initialize(grid, flat, 1, 0, {}, false, root, 30000), std::runtime_error);

    // Delamination ramp (continental_root_delamination_width_m): the root's load is full
    // up to the cap, falls linearly to zero at cap + width, and the crust-only Airy height
    // holds above; continuous, strictly increasing, inverted by seed().
    constexpr double width = 10000;
    const auto ramped = [&](double continent, double w, double r) {
        Cell c; c.continental_area = 1; c.continent = continent;
        return elevation_m(c, 1, r, cap, w);
    };
    const auto crust_only = [&](double continent) {
        Cell c; c.continental_area = 1; c.continent = continent;
        return elevation_m(c, 1, 0);
    };
    // Continuity at both ends: the frozen cap's height at the cap, no root at cap + width.
    EXPECT_NEAR(ramped(cap, width, root), at(cap), 1e-9);
    EXPECT_NEAR(ramped(cap, width, root), cap_z, 1e-9);
    EXPECT_NEAR(ramped(cap + width, width, root), crust_only(cap + width), 1e-9);
    EXPECT_EQ(ramped(cap + width + 1, width, root), crust_only(cap + width + 1));
    EXPECT_NEAR(ramped(60000, width, root), 25000. * 500 / 3300, 1e-9); // 3788 m
    EXPECT_NEAR(ramped(70000, width, root), 35000. * 500 / 3300, 1e-9); // 5303 m
    EXPECT_NEAR(ramped(cap + width / 2, width, root), (500. * 20000 - root * .5 * 20000) / 3300, 1e-9);
    const double mid_z = ramped(cap + width / 2, width, root), top_z = ramped(cap + width, width, root);
    for (double z : {0., 500., cap_z, mid_z, top_z, 3000., 5000., 8000.}) {
        const auto c = seed(z, true, 100, root, cap, width);
        EXPECT_NEAR(elevation_m(c, 1, root, cap, width), z, 1e-9) << z;
    }
    EXPECT_NEAR(seed(cap_z, true, 100, root, cap, width).continent, cap, 1e-6);
    EXPECT_NEAR(seed(mid_z, true, 100, root, cap, width).continent, cap + width / 2, 1e-6);
    EXPECT_NEAR(seed(top_z, true, 100, root, cap, width).continent, cap + width, 1e-6);
    EXPECT_NEAR(seed(-3500, true, 100, root, cap, width).continent, seed(-3500, true, 100, root, cap).continent, 1e-9);
    // Strictly increasing through the ramp for every width and root: d(buoyancy)/dH =
    // 500 - root (1 - (2H - cap - 35 km) / width) > 500 - root > 0 on (cap, cap + width).
    for (double w : {1., 1000., width, 1e6}) for (double r : {root, 499.}) {
        const double step = (w + 2000) / 2000;
        double previous = ramped(cap - 1000, w, r);
        for (double h = cap - 1000 + step; h <= cap + w + 1000; h += step) {
            const double z = ramped(h, w, r);
            ASSERT_GT(z, previous) << "width " << w << " root " << r << " h " << h;
            previous = z;
        }
        for (double z : {0., 2000., 5000., 8000.})
            EXPECT_NEAR(elevation_m(seed(z, true, 100, r, cap, w), 1, r, cap, w), z, 1e-8) << w << ' ' << r << ' ' << z;
    }
    // Mixed cells by area: the per-area thickness picks the factor.
    Cell ramp_half = half; ramp_half.continent = .5 * (cap + width / 2); // factor 1/2
    EXPECT_NEAR(elevation_m(ramp_half, 1, root, cap, width), expression(ramp_half, 1, root * .5, ramp_half.continent), 1e-9);
    EXPECT_NEAR(elevation_m(half, 1, root, cap, width), uncapped(half, 1, 0), 1e-9); // 60 km per area: no root
    EXPECT_NEAR(elevation_m(mixed, 1, root, cap, width), uncapped(mixed, 1, 0), 1e-9); // 61.2 km per area
    Cell ghost; ghost.continent = 1e-3; ghost.ocean = 7100; ghost.oceanic_area = 1; ghost.oceanic_age_volume = ghost.ocean * 80;
    for (double thermal : {0., 1.})
        EXPECT_EQ(elevation_m(ghost, thermal, root, cap, width), elevation_m(ghost, thermal, root, cap));
    // Width 0: the frozen cap's expression bit for bit (mixed, 70 km and half cells).
    for (const auto& c : {mixed, thick, half}) for (double thermal : {0., 1.}) {
        EXPECT_EQ(elevation_m(c, thermal, root, cap, 0),
                  expression(c, thermal, root, std::min(c.continent, cap * c.continental_area)));
        EXPECT_EQ(elevation_m(c, thermal, root, cap, 0), elevation_m(c, thermal, root, cap));
    }
    for (double z : {-3500., 0., 1500., 8000.})
        EXPECT_EQ(seed(z, true, 100, root, cap, 0).continent, seed(z, true, 100, root, cap).continent);
    // Validation: 0, or finite and positive with a cap.
    EXPECT_THROW(seed(0, true, 100, root, cap, -1), std::runtime_error);
    EXPECT_THROW(seed(0, true, 100, root, 0, width), std::runtime_error);
    EXPECT_THROW(seed(0, true, 100, root, cap, std::numeric_limits<double>::quiet_NaN()), std::runtime_error);
    EXPECT_THROW(seed(0, true, 100, root, cap, std::numeric_limits<double>::infinity()), std::runtime_error);
    EXPECT_NO_THROW(model.initialize(grid, flat, 1, 0, {}, false, root, cap, width));
    EXPECT_THROW(model.initialize(grid, flat, 1, 0, {}, false, root, 0, width), std::runtime_error);
}

TEST(TectonicColumns, StationaryCoolingIsIndependentOfTimestepAndDoesNotRemoveCrust) {
    const auto grid = cubefixture::cube(4);
    Model a, b;
    const std::vector<Cell> initial(grid->cellcount(), seed(-5000, false, 20));
    a.initialize(grid, initial); b.initialize(grid, initial);
    const auto still = cubefixture::still(*grid);
    for (int i = 0; i < 20; ++i) a.advance(still, 1, 0, 1);
    for (int i = 0; i < 80; ++i) b.advance(still, .25, 0, 1);
    for (CellId c = 0; c < grid->cellcount(); ++c) {
        EXPECT_NEAR(a.cells()[c].ocean, 7000, 1e-9);
        EXPECT_NEAR(age_myr(a.cells()[c]), 40, 1e-9);
        EXPECT_NEAR(elevation_m(a.cells()[c]), elevation_m(b.cells()[c]), 1e-9);
    }
    EXPECT_NEAR(a.budget().relative_mass_residual, 0, 1e-14);
}

// Two continental halves of a 100 km band (N 8) close at 45 degrees and part at -135: the
// continental area changes by the budget's area terms only, the middle collides and the far
// side opens (stretched without margins), mass is conserved and the collision builds
// mountains (a root, no 4 km ceiling).
TEST(TectonicColumns, ConvergenceThickensContinentsWithoutInventingOrDeletingMass) {
    const Band band(8, cubefixture::radiusfor_km(8, 100000));
    Model world;
    world.initialize(band.grid, std::vector<Cell>(band.cellcount(), seed(0, true, 100)));
    const auto motion = converging(band, .05);
    for (int i = 0; i < 80; ++i) {
        const double area = area_weighted(world, &Cell::continental_area), terms = continental_area_terms(world.budget());
        world.advance(motion, .25, 0, 1);
        EXPECT_NEAR(area_weighted(world, &Cell::continental_area) - area, continental_area_terms(world.budget()) - terms,
                    1e-9 * band.spacing_m() * band.spacing_m());
    }
    EXPECT_GT(world.budget().collided_continental_area_m2, 0);
    EXPECT_GT(world.budget().stretched_continental_area_m2, 0);
    EXPECT_DOUBLE_EQ(world.budget().arc_continental_area_m2, 0);
    EXPECT_DOUBLE_EQ(world.budget().oceanized_continental_area_m2, 0);
    EXPECT_DOUBLE_EQ(world.budget().continentalized_oceanic_area_m2, 0);
    double peak = 0;
    for (const auto& c : world.cells()) peak = std::max(peak, elevation_m(c));
    RecordProperty("collision_peak_m", std::to_string(peak));
    EXPECT_GT(peak, 5000);
    EXPECT_LT(std::abs(world.budget().relative_mass_residual), 1e-12);
    EXPECT_DOUBLE_EQ(world.budget().subducted_ocean_kg, 0);
    EXPECT_DOUBLE_EQ(world.budget().generated_ocean_kg, 0);
}

TEST(TectonicColumns, OceanOpeningAndSubductionAreExplicitMassSourcesAndSinks) {
    const Band band(8, cubefixture::radiusfor_km(8, 100000));
    Model world;
    world.initialize(band.grid, std::vector<Cell>(band.cellcount(), seed(-4000, false, 80)));
    const auto motion = converging(band, .05);
    for (int i = 0; i < 20; ++i) world.advance(motion, 1, 0, 1);
    EXPECT_GT(world.budget().generated_ocean_kg, 0);
    EXPECT_GT(world.budget().subducted_ocean_kg, 0);
    EXPECT_LT(std::abs(world.budget().relative_mass_residual), 1e-12);
    // Crust of the ridge's thickness everywhere: the per-plate transport keeps a value per
    // covered area between its donors'.
    for (const auto& c : world.cells()) EXPECT_NEAR(c.ocean, 7000, 1e-6);
}

TEST(TectonicColumns, ErodedRockReappearsAsSedimentWithEqualMass) {
    const Band band(4, cubefixture::radiusfor_km(4, 100000));
    const auto cells = cubefixture::layout(band, seed(0, true, 100), [](int x, int) {
        return x % 16 == 8 ? seed(6000, true, 100) : seed(0, true, 100); });
    Model world;
    world.initialize(band.grid, cells);
    const auto still = cubefixture::still(*band.grid);
    for (int i = 0; i < 20; ++i) world.advance(still, 1, 1, 1);
    EXPECT_GT(world.budget().eroded_kg, 0);
    EXPECT_NEAR(world.budget().deposited_kg / world.budget().eroded_kg, 1, 1e-12);
    EXPECT_LT(std::abs(world.budget().relative_mass_residual), 1e-12);
    EXPECT_LT(elevation_m(world.cells()[band.cell(8, 1)]), 6000);
    EXPECT_GT(world.cells()[band.cell(7, 1)].sediment, 0);
}

// A rigid rotation carries an anomaly of older crust across a face edge, the date line and the
// north pole: the continental volume and the age moment are conserved (the age moment grows by
// the update's time), the anomaly's excess moves with the rotation, and a displacement beyond
// 0.95 of a cell's smallest edge is refused.
TEST(TectonicColumns, FractionalTransportCrossesFaceEdgesTheDateLineAndThePoleAndPreservesVolumeAndAge) {
    Model uninitialized;
    EXPECT_THROW(uninitialized.advance(st::Motion{}, 1, 0, 1), std::runtime_error);
    constexpr int n = 8;
    const auto grid = cubefixture::cube(n);
    const auto& mesh = *grid->mesh();
    struct Case { const char* name; Vec3 start; Vec3 axis; };
    const double column = pi / (2.0 * n);
    const Case cases[] = {
        {"face_edge", mesh::cubepoint(0, pi / 4 - column / 2, 0), {0, 0, 1}},
        {"date_line", mesh::unitvector(0, pi - column / 2), {0, 0, 1}},
        {"pole", mesh::unitvector(pi / 2 - column, -pi / 2), {1, 0, 0}},
    };
    for (const auto& test : cases) {
        auto cells = std::vector<Cell>(grid->cellcount(), seed(0, true, 100));
        const CellId anomaly = mesh.locate(test.start);
        cells[anomaly] = seed(0, true, 250);  // Uniform elevation isolates advection from root flow.
        Model world;
        world.initialize(grid, cells);
        const double volume = area_weighted(world, &Cell::continent);
        const double moment = area_weighted(world, &Cell::continental_age_volume);
        const Vec3 rotation = .25 * column * test.axis;
        world.advance(cubefixture::rotating(*grid, rotation), 1, 0, 0);
        EXPECT_NEAR(area_weighted(world, &Cell::continent) / volume, 1, 1e-12) << test.name;
        EXPECT_NEAR(area_weighted(world, &Cell::continental_age_volume) / (moment + volume), 1, 1e-12) << test.name;
        EXPECT_LT(std::abs(world.budget().relative_mass_residual), 1e-12) << test.name;
        // The excess age (over the background's 101 Myr) has moved by the rotation.
        Vec3 centre;
        double excess = 0;
        for (CellId c = 0; c < grid->cellcount(); ++c) {
            const double e = (age_myr(world.cells()[c]) - 101) * world.cells()[c].continent * grid->area_m2(c);
            centre = centre + e * grid->position(c);
            excess += e;
        }
        EXPECT_NEAR(excess / (150 * cells[anomaly].continent * grid->area_m2(anomaly)), 1, 1e-6) << test.name;
        const double angle = mesh::length(rotation);
        const Vec3 k = mesh::normalise(rotation), p = grid->position(anomaly);
        const Vec3 moved = std::cos(angle) * p + std::sin(angle) * mesh::cross(k, p) + (mesh::dot(k, p) * (1 - std::cos(angle))) * k;
        const double error = mesh::greatcircle_rad(mesh::normalise(centre), moved) / column;
        RecordProperty(std::string(test.name) + "_centre_error_columns", std::to_string(error));
        EXPECT_LT(error, .1) << test.name;
        // The next cell along the motion (a column on) got about a quarter of the anomaly.
        const CellId receiver = mesh.locate(mesh::normalise(p + (4 * angle) * mesh::cross(k, p)));
        EXPECT_NE(receiver, anomaly) << test.name;
        EXPECT_GT(age_myr(world.cells()[receiver]), 101 + .15 * 150) << test.name;
    }
    Model world;
    world.initialize(grid, std::vector<Cell>(grid->cellcount(), seed(0, true, 100)));
    EXPECT_THROW(world.advance(cubefixture::rotating(*grid, {0, 0, 1.1 * column}), 1, 0, 1), std::runtime_error);
}

TEST(TectonicColumns, RootsSpreadWithoutErosionAndPreserveMass) {
    const Band band(4, cubefixture::radiusfor_km(4, 100000));
    const auto cells = cubefixture::layout(band, seed(0, true, 100), [](int x, int) {
        return x % 16 == 8 ? seed(12000, true, 100) : seed(0, true, 100); });
    Model world;
    world.initialize(band.grid, cells);
    const auto still = cubefixture::still(*band.grid);
    for (int i = 0; i < 20; ++i) world.advance(still, 1, 0, 1);
    EXPECT_LT(elevation_m(world.cells()[band.cell(8, 1)]), 6000);
    EXPECT_GT(elevation_m(world.cells()[band.cell(7, 1)]), 0);
    EXPECT_GT(world.budget().redistributed_continent_kg, 0);
    EXPECT_DOUBLE_EQ(world.budget().eroded_kg, 0);
    EXPECT_LT(std::abs(world.budget().relative_mass_residual), 1e-12);
}

TEST(TectonicColumns, GentleInteriorReliefSurvivesRootFlowOver200Myr) {
    const Band band(8, cubefixture::radiusfor_km(8, 156250));
    const auto relief = [](int x) { return 300 * std::sin(2 * pi * x / 8); };
    const auto initial = cubefixture::layout(band, seed(700, true, 100), [&](int x, int) {
        return seed(700 + relief(x), true, 100); });
    Model world; world.initialize(band.grid, initial);
    const auto still = cubefixture::still(*band.grid);
    for (int i = 0; i < 200; ++i) world.advance(still, 1, 0, 1);
    double amplitude = 0, norm = 0;
    for (int y = 0; y < band.rows(); ++y) for (int x = 0; x < band.columns(); ++x) {
        amplitude += (elevation_m(world.cells()[band.cell(x, y)]) - 700) * relief(x);
        norm += relief(x) * relief(x);
    }
    const double retained = amplitude / norm;
    RecordProperty("root_flow_relief_retained", retained);
    EXPECT_GT(retained, .95);
    EXPECT_LT(std::abs(world.budget().relative_mass_residual), 1e-12);
}

TEST(TectonicColumns, VanishingContinentalTracesHaveZeroRootConductance) {
    auto tiny = seed(0, true, 100);
    tiny.continent = 1e-80;
    tiny.continental_age_volume = tiny.continent * 100;
    const auto grid = cubefixture::cube(2, cubefixture::radiusfor_km(2, 156250));
    Model world; world.initialize(grid, std::vector<Cell>(grid->cellcount(), tiny));
    for (int i = 0; i < 5; ++i) world.advance(cubefixture::still(*grid), 1, 0, 0);
    for (const auto& c : world.cells()) {
        EXPECT_DOUBLE_EQ(c.continent, tiny.continent);
        EXPECT_NEAR(age_myr(c), 105, 1e-12);
    }
    EXPECT_DOUBLE_EQ(world.budget().redistributed_continent_kg, 0);
}

TEST(TectonicColumns, RootFlowAndErosionConvergeWhenTimestepIsReduced) {
    const Band band(4, cubefixture::radiusfor_km(4, 100000));
    const auto cells = cubefixture::layout(band, seed(0, true, 100), [](int x, int) {
        return x % 16 == 8 ? seed(6000, true, 100) : seed(0, true, 100); });
    const auto still = cubefixture::still(*band.grid);
    const auto run = [&](double dt) {
        Model world; world.initialize(band.grid, cells);
        for (int i = 0; i < static_cast<int>(10 / dt); ++i) world.advance(still, dt, 1, 1);
        EXPECT_LT(std::abs(world.budget().relative_mass_residual), 1e-12);
        return world.cells();
    };
    // Stay below the adaptive root-flow stability step; halving a large outer
    // step need not halve its internally selected effective timestep.
    const auto coarse = run(.125), medium = run(.0625), fine = run(.015625);
    double coarse_error = 0, medium_error = 0;
    for (size_t i = 0; i < cells.size(); ++i) {
        coarse_error += std::abs(elevation_m(coarse[i]) - elevation_m(fine[i]));
        medium_error += std::abs(elevation_m(medium[i]) - elevation_m(fine[i]));
    }
    EXPECT_LT(medium_error, coarse_error * .8);
}

TEST(TectonicColumns, InitializationKeepsLowContinentalLandAboveWater) {
    tectonics::scenario::Scenario s;
    s.seed = 1; s.sphere_face_cells = 8; s.plate_count = 4;
    s.evolution_duration_myr = 1;
    lithosphere world(s);
    const auto datum = world.getSeaLevelMeters();
    std::vector<uint16_t> lowland(world.getGrid().cellcount(), datum + 200);
    world.importMetricHeightMap(lowland.data(), datum);
    for (const auto& c : world.getColumnModel().cells()) {
        EXPECT_DOUBLE_EQ(c.continental_area, 1);
        EXPECT_NEAR(elevation_m(c), 200, .1);
    }
    // With an Earth-like lithospheric root the same lowland stands on a 37.6 km column.
    s.continental_root_kg_m3 = 246;
    lithosphere rooted(s);
    rooted.importMetricHeightMap(lowland.data(), datum);
    for (const auto& c : rooted.getColumnModel().cells()) {
        EXPECT_DOUBLE_EQ(c.continental_area, 1);
        EXPECT_NEAR(elevation_m(c, 1, 246), 200, .1);
        EXPECT_NEAR(c.continent, 35000 + 200 * 3300. / 254, 2);
    }
}

// The continental fraction seeds shelves: the continental area (true areas) reaches the
// requested share, the land share stays.
TEST(TectonicColumns, ContinentalFractionSeedsShelvesToTheRequestedArea) {
    tectonics::scenario::Scenario s;
    s.seed = 1; s.sphere_face_cells = 32; s.plate_count = 8; s.sea_level = .65f;
    s.evolution_duration_myr = 1; s.boundary_conditioned_terrain = true;
    const auto share = [](lithosphere& world, bool columns) {
        double continental = 0, total = 0;
        const auto& cells = world.getColumnModel().cells();
        for (CellId c = 0; c < cells.size(); ++c) {
            const double area = world.getGrid().area_m2(c);
            continental += (columns ? cells[c].continental_area : (cells[c].continental_area >= .5 ? 1.0 : 0.0)) * area;
            total += area;
        }
        return continental / total;
    };
    const auto land = [](lithosphere& world) {
        double above = 0, total = 0;
        const auto& cells = world.getColumnModel().cells();
        for (CellId c = 0; c < cells.size(); ++c) {
            above += (elevation_m(cells[c]) > 0 ? 1.0 : 0.0) * world.getGrid().area_m2(c);
            total += world.getGrid().area_m2(c);
        }
        return above / total;
    };
    lithosphere legacy(s);
    s.continental_fraction = .40;
    lithosphere shelved(s);
    EXPECT_GT(share(legacy, false), .5);
    EXPECT_NEAR(share(shelved, false), .40, .01);
    // The land share is unchanged: the shelf lies below sea level.
    EXPECT_NEAR(land(shelved), land(legacy), .01);
    // With the baseline's 2000 m margin taper the slope continues to -3500 m, so the
    // continental area summed over the columns, not only the cells over one half, is f.
    s.initial_margin_width_m = 2000;
    lithosphere tapered(s);
    EXPECT_NEAR(share(tapered, true), .40, .01);
    s.initial_margin_width_m = 0;
    s.continental_fraction = .35; // 1 - sea_level: no shelf, a 2.5 km cliff at the coast
    EXPECT_THROW(lithosphere{s}, std::runtime_error);
    s.continental_fraction = .30;
    EXPECT_THROW(lithosphere{s}, std::runtime_error);
}

TEST(TectonicColumns, Integrated200MyrReplayClosesMassBudgetAndPreservesContinents) {
    tectonics::scenario::Scenario s;
    s.seed = 1; s.sphere_face_cells = 16; s.plate_count = 8;
    s.evolution_duration_myr = 200;
    lithosphere a(s), b(s);
    double initial_continent = 0;
    double initial_max_age = 0;
    for (const auto& c : a.getColumnModel().cells()) {
        initial_continent += c.continent;
        initial_max_age = std::max(initial_max_age, age_myr(c));
    }
    for (int i = 0; i < 200; ++i) {
        a.update(); b.update();
        ASSERT_LT(std::abs(a.getColumnModel().budget().relative_mass_residual), 1e-10) << i;
        for (const auto& c : a.getColumnModel().cells())
            ASSERT_LE(age_myr(c), initial_max_age + i + 1 + 1e-6) << "update " << i;
    }
    EXPECT_TRUE(a.isFinished());
    EXPECT_DOUBLE_EQ(a.getTimeMyr(), 200);
    EXPECT_GT(initial_continent, 0);
    for (size_t i = 0; i < a.getColumnModel().cells().size(); ++i) {
        EXPECT_DOUBLE_EQ(a.getColumnModel().cells()[i].continent, b.getColumnModel().cells()[i].continent);
        EXPECT_FLOAT_EQ(a.getTopography()[i], b.getTopography()[i]);
    }
}

TEST(TectonicColumns, ResolvedRiftsCreateIndependentPlatesAndReplayWithoutLosingMaterial) {
    // Stronger forcing is a controlled failure fixture, not the benchmark recipe: two plates and
    // a 40 cm/yr mantle flow (the flow varies over about a radius: a plate needs strong flow
    // to fail within 100 Myr).
    tectonics::scenario::Scenario s;
    s.seed=1; s.sphere_face_cells=16; s.plate_count=2;  // large plates: the tension grows with their area
    s.evolution_duration_myr=100; s.delta_time_myr=.5;
    s.mantle_flow_speed_cm_yr=40;
    lithosphere a(s), b(s);
    for (int step=0;step<200;++step) {
        a.update(); b.update();
        ASSERT_LT(std::abs(a.getColumnModel().budget().relative_mass_residual),1e-10);
        ASSERT_EQ(a.getTopologyDiagnostics().stable_ids,b.getTopologyDiagnostics().stable_ids);
    }
    const auto& events=a.getTopologyDiagnostics().events;
    RecordProperty("max_tension_n_m",std::to_string(a.getTopologyDiagnostics().max_tension_n_m));
    RecordProperty("max_rift_damage",std::to_string(a.getTopologyDiagnostics().max_rift_damage));
    ASSERT_GE(events.size(),1U);
    EXPECT_EQ(events[0].kind,tectonics::topology::EventKind::Split);
    EXPECT_GT(a.getTopologyDiagnostics().stable_ids.size(),s.plate_count);
    EXPECT_EQ(events.size(),b.getTopologyDiagnostics().events.size());
    for (size_t i=0;i<events.size();++i) {
        EXPECT_DOUBLE_EQ(events[i].time_myr,b.getTopologyDiagnostics().events[i].time_myr);
        EXPECT_EQ(events[i].result,b.getTopologyDiagnostics().events[i].result);
    }
    for (uint32_t i=0;i<a.getGrid().cellcount();++i) {
        EXPECT_FLOAT_EQ(a.getTopography()[i],b.getTopography()[i]);
        EXPECT_EQ(a.getPlatesMap()[i],b.getPlatesMap()[i]);
        EXPECT_DOUBLE_EQ(a.getColumnModel().cells()[i].continental_age_volume,
                         b.getColumnModel().cells()[i].continental_age_volume);
    }
    RecordProperty("split_events",static_cast<int>(events.size()));
}

TEST(TectonicColumns, ExcessMountainForceMatchesAiryColumnUnitsAndVanishesOnUnthickenedCrust) {
    EXPECT_DOUBLE_EQ(mountain_force_n_m(seed(0,true,100)),0);
    EXPECT_DOUBLE_EQ(mountain_force_n_m(seed(-4000,false,100)),0);
    const auto mountain=seed(5000,true,100);
    const double expected=.5*9.80665*2800*(1-2800./3300)*(68000.*68000-35000.*35000);
    EXPECT_NEAR(mountain_force_n_m(mountain),expected,expected*1e-14);
    EXPECT_GT(mountain_force_n_m(seed(10000,true,100)),mountain_force_n_m(mountain));
    // A lithospheric root scales the force with the elevation per metre of crust.
    EXPECT_NEAR(mountain_force_n_m(mountain,false,246),expected*.508,expected*1e-14);
    auto covered=mountain; covered.sediment=1500;
    EXPECT_NEAR(mountain_force_n_m(covered,true,246),.508*mountain_force_n_m(covered,true),
                mountain_force_n_m(covered,true)*1e-14);
    EXPECT_EQ(mountain_force_n_m(covered,true,0),mountain_force_n_m(covered,true));
}

// Every field of two snapshots is equal.
static bool same_fields(const tectonics::contract::Snapshot& a, const tectonics::contract::Snapshot& b) {
    return a.continental_crust_m == b.continental_crust_m && a.oceanic_crust_m == b.oceanic_crust_m &&
        a.sediment_m == b.sediment_m && a.rift_damage == b.rift_damage && a.plume_influence == b.plume_influence &&
        a.heightmap == b.heightmap && a.plate_id == b.plate_id && a.crust_age_steps == b.crust_age_steps &&
        a.crust_age_myr == b.crust_age_myr && a.crust_thickness == b.crust_thickness && a.crust_class == b.crust_class &&
        a.uplift_tendency == b.uplift_tendency && a.subsidence_tendency == b.subsidence_tendency &&
        a.accumulated_strain == b.accumulated_strain && a.boundary_type == b.boundary_type &&
        a.boundary_distance == b.boundary_distance && a.boundary_segment_id == b.boundary_segment_id &&
        a.nearest_boundary_id == b.nearest_boundary_id && a.deforming_region_id == b.deforming_region_id &&
        a.deforming_region_type == b.deforming_region_type && a.deformation_rate == b.deformation_rate &&
        a.deformation_velocity_x == b.deformation_velocity_x && a.deformation_velocity_y == b.deformation_velocity_y &&
        a.convergence_score == b.convergence_score && a.divergence_score == b.divergence_score &&
        a.shear_score == b.shear_score && a.geologic_regime == b.geologic_regime;
}

TEST(TectonicColumns, PhysicalLayersSurviveBundleReadAndAffectComparison) {
    tectonics::scenario::Scenario s;
    s.seed = 1; s.sphere_face_cells = 8; s.plate_count = 4;
    s.evolution_duration_myr = 1;
    lithosphere world(s); world.update();
    auto snapshot = tectonics::contract::capture_snapshot(world);
    EXPECT_EQ(snapshot.width, 8u);
    EXPECT_EQ(snapshot.height, 48u);
    snapshot.column_budget.noise_redistributed_kg=123456;
    snapshot.column_budget.redistributed_sediment_kg=654321;
    snapshot.topology_diagnostics.events = {{tectonics::topology::EventKind::Split,1,0,0,4,128},
        {tectonics::topology::EventKind::Fusion,2,1,2,5,200}};
    const auto folder = std::filesystem::temp_directory_path() / "uw-isostatic-column-test";
    tectonics::tooling::write_snapshot_bundle(folder, snapshot, {});
    auto loaded = tectonics::tooling::load_snapshot_bundle(folder).snapshot;
    EXPECT_EQ(loaded.continental_crust_m, snapshot.continental_crust_m);
    EXPECT_EQ(loaded.oceanic_crust_m, snapshot.oceanic_crust_m);
    EXPECT_EQ(loaded.sediment_m, snapshot.sediment_m);
    EXPECT_DOUBLE_EQ(loaded.column_budget.mass_kg, snapshot.column_budget.mass_kg);
    EXPECT_DOUBLE_EQ(loaded.column_budget.noise_redistributed_kg,123456);
    EXPECT_DOUBLE_EQ(loaded.column_budget.redistributed_sediment_kg,654321);
    EXPECT_TRUE(loaded.topology_diagnostics.enabled);
    EXPECT_EQ(loaded.topology_diagnostics.reassigned_cells, snapshot.topology_diagnostics.reassigned_cells);
    EXPECT_EQ(loaded.topology_diagnostics.max_detached_components, snapshot.topology_diagnostics.max_detached_components);
    EXPECT_EQ(loaded.topology_diagnostics.stable_ids, snapshot.topology_diagnostics.stable_ids);
    ASSERT_EQ(loaded.topology_diagnostics.events.size(),2U);
    EXPECT_EQ(loaded.topology_diagnostics.events[0].kind,tectonics::topology::EventKind::Split);
    EXPECT_EQ(loaded.topology_diagnostics.events[1].kind,tectonics::topology::EventKind::Fusion);
    EXPECT_EQ(loaded.topology_diagnostics.events[1].result,5U);
    EXPECT_DOUBLE_EQ(loaded.topology_diagnostics.events[1].time_myr,2);
    EXPECT_TRUE(same_fields(snapshot, loaded));
    loaded.continental_crust_m[0] += 10;
    EXPECT_FALSE(same_fields(snapshot, loaded));
}

TEST(TectonicColumns, SedimentLoadContributesToLayeredGpeAndConservativeRootFlow) {
    for (double z:{-1000.,0.,1000.,5000.}) {
        const auto c=seed(z,true,40);
        EXPECT_NEAR(mountain_force_n_m(c,true),mountain_force_n_m(c),1);
    }
    auto column=seed(0,true,40);column.sediment=27500;
    EXPECT_NEAR(elevation_m(column,0),5000,1e-9);
    EXPECT_DOUBLE_EQ(mountain_force_n_m(column),0);
    EXPECT_NEAR(mountain_force_n_m(column,true),.5*9.80665*1.15875e12,1);
    const Band band(8, cubefixture::radiusfor_km(8, 50000));
    const auto cells = cubefixture::layout(band, seed(0,true,40), [&](int x, int) {
        auto c = seed(0,true,40);
        c.sediment = 50000*std::pow(.5+.5*std::cos(2*pi*x/band.columns()),4);
        return c; });
    Model loaded,old;loaded.initialize(band.grid,cells);old.initialize(band.grid,cells);
    const auto stopped = cubefixture::still(*band.grid);
    for (int i=0;i<10;++i) {
        loaded.advance(stopped,1,0,0,0,nullptr,{}, {},true);
        old.advance(stopped,1,0,0,0);
    }
    double initial_s=0,final_s=0,initial_c=0,final_c=0,initial_age=0,final_age=0,peak=-1e9,old_peak=-1e9;
    for (CellId i=0;i<cells.size();++i) {
        const auto& c=loaded.cells()[i];
        const double area=band.grid->area_m2(i);
        initial_s+=cells[i].sediment*area;final_s+=c.sediment*area;
        initial_c+=cells[i].continent*area;final_c+=c.continent*area;
        initial_age+=cells[i].continental_age_volume*area;final_age+=c.continental_age_volume*area;
        EXPECT_GE(c.sediment,0);EXPECT_GE(c.continent,0);
        peak=std::max(peak,elevation_m(c,0));old_peak=std::max(old_peak,elevation_m(old.cells()[i],0));
    }
    EXPECT_LT(peak,old_peak);
    EXPECT_NEAR(initial_s,final_s,initial_s*1e-12);
    EXPECT_NEAR(initial_c,final_c,initial_c*1e-12);
    EXPECT_NEAR(initial_age+10*initial_c,final_age,final_age*1e-12);
    EXPECT_GT(loaded.budget().redistributed_sediment_kg,0);
    EXPECT_DOUBLE_EQ(loaded.budget().eroded_kg,0);
    EXPECT_LT(std::abs(loaded.budget().relative_mass_residual),1e-12);
}

TEST(TectonicColumns, SedimentLoadedFlowConvergesWithSmallerTimesteps) {
    const Band band(4, cubefixture::radiusfor_km(4, 100000));
    const auto cells = cubefixture::layout(band, seed(0,true,40), [](int x, int) {
        auto c = seed(0,true,40); c.sediment = x%16==8?25000:1000; return c; });
    const auto stopped = cubefixture::still(*band.grid);
    const auto run=[&](double dt) {
        Model model;model.initialize(band.grid,cells);
        for (int i=0;i<static_cast<int>(10/dt);++i)
            model.advance(stopped,dt,0,0,0,nullptr,{}, {},true);
        EXPECT_LT(std::abs(model.budget().relative_mass_residual),1e-12);
        return model.cells();
    };
    const auto coarse=run(.125),medium=run(.0625),fine=run(.015625);
    double coarse_error=0,medium_error=0;
    for (size_t i=0;i<cells.size();++i) {
        coarse_error+=std::abs(elevation_m(coarse[i],0)-elevation_m(fine[i],0));
        medium_error+=std::abs(elevation_m(medium[i],0)-elevation_m(fine[i],0));
    }
    EXPECT_LT(medium_error,.8*coarse_error);
}

// The rock uplift rate (erosion plan D1) counts the tectonic sub-steps only: the
// diffusion and the rivers change the columns, never the rate; a plume's basalt on a
// continental column raises it by 400/3300 of its thickness (Airy, land), so 330 m
// over 1 Myr is a step rate of 40 m/Myr, kept as an exponential mean with tau = 5 Myr.
TEST(TectonicColumns, RockUpliftRateCountsTectonicChangeOnly) {
    tectonics::margins::Parameters still; // No root flow: the roots would carry the change sideways.
    still.root_mobility = 0;
    {
        const Band band(4, cubefixture::radiusfor_km(4, 100000));
        const auto cells = cubefixture::layout(band, seed(500, true, 100), [](int x, int) {
            return x % 16 == 8 ? seed(3000, true, 100) : seed(500, true, 100); });
        auto rivers = still;
        rivers.enabled = true; // routeSediment with the default river erosion; no forcing, no thermal change
        Model world;
        world.initialize(band.grid, cells);
        for (int i = 0; i < 3; ++i) world.advance(cubefixture::still(*band.grid), 1, 1, 1, 1, nullptr, rivers);
        EXPECT_GT(world.budget().eroded_kg, 0);
        EXPECT_GT(world.budget().river_eroded_kg, 0);
        EXPECT_LT(elevation_m(world.cells()[band.cell(8, 1)]), 3000);
        ASSERT_EQ(world.upliftRate().size(), cells.size());
        for (double rate : world.upliftRate()) EXPECT_NEAR(rate, 0, 1e-11);
    }
    const auto grid = cubefixture::cube(4, cubefixture::radiusfor_km(4, 100000));
    const CellId lifted = 9;
    tectonics::plumes::Source plume;
    plume.basalt_m.assign(grid->cellcount(), 0);
    plume.influence.assign(grid->cellcount(), 0);
    plume.basalt_m[lifted] = 330;
    Model world;
    world.initialize(grid, std::vector<Cell>(grid->cellcount(), seed(500, true, 100)));
    const auto at_rest = cubefixture::still(*grid);
    world.advance(at_rest, 1, 0, 1, 1, nullptr, still, {}, false, &plume, 5);
    const double weight = -std::expm1(-1.0 / 5); // 1 - e^(-dt/tau)
    EXPECT_NEAR(world.upliftRate()[lifted], 40 * weight, 1e-11);
    EXPECT_NEAR(world.upliftRate()[lifted], 7.2507698768807245, 1e-11);
    for (CellId i = 0; i < grid->cellcount(); ++i) if (i != lifted) EXPECT_NEAR(world.upliftRate()[i], 0, 1e-11);
    // A step without tectonic change decays the mean by e^(-dt/tau).
    world.advance(at_rest, 1, 0, 1, 1, nullptr, still, {}, false, nullptr, 5);
    EXPECT_NEAR(world.upliftRate()[lifted], 40 * weight * std::exp(-1.0 / 5), 1e-11);
    EXPECT_THROW(world.advance(at_rest, 1, 0, 1, 1, nullptr, still, {}, false, nullptr, 0),
        std::runtime_error);
}

namespace {
// The largest |rate| of a model, and the largest Eulerian rate (fixed-cell elevation change
// over dt) of the same step, over the band's cells.
struct RateMeasure { double measured = 0, eulerian = 0; };
RateMeasure step_rates(Model& world, const Band& band, const st::Motion& motion, const tectonics::margins::Parameters& p,
                       const tectonics::plumes::Source* plume = nullptr) {
    std::vector<double> before(world.cells().size());
    for (size_t i = 0; i < before.size(); ++i) before[i] = elevation_m(world.cells()[i]);
    world.advance(motion, 1, 0, 1, 1, nullptr, p, {}, false, plume, 5);
    const double weight = -std::expm1(-1.0 / 5);
    RateMeasure m;
    for (int y = 0; y < band.rows(); ++y) for (int x = 0; x < band.columns(); ++x) {
        const CellId c = band.cell(x, y);
        m.measured = std::max(m.measured, std::abs(world.upliftRate()[c]));
        m.eulerian = std::max(m.eulerian, std::abs(elevation_m(world.cells()[c]) - before[c]) * weight);
    }
    return m;
}
}

// The measurement is in the plate's frame (erosion plan H2): plate motion over a slope is not
// uplift. Land columns 500 + 100 x m (x the band's column, a 1500 m jump at the band's start)
// move east by 0.25 column with nothing else acting. The elevation before the step is carried
// with the rock by the same transport as the columns; the per-plate normalisation scales each
// field by its own ripple factor (|scale - 1| about 1e-5), so the rate is zero to that factor's
// share of the carried values instead of to rounding. The fixed-cell (Eulerian) rate is
// -v . grad z, 25 m/Myr (r = 4.53).
TEST(TectonicColumns, RockUpliftRateIsZeroForPlateMotionOverASlope) {
    tectonics::margins::Parameters still;
    still.root_mobility = 0;
    const Band band(4, cubefixture::radiusfor_km(4, 100000));
    const auto cells = cubefixture::layout(band, seed(500, true, 100), [](int x, int) {
        return seed(500 + 100. * x, true, 100); });
    Model world;
    world.initialize(band.grid, cells);
    const auto m = step_rates(world, band, cubefixture::rotating(*band.grid, band.eastward(.25)), still);
    RecordProperty("largest_rate_m_myr", std::to_string(m.measured));
    RecordProperty("largest_eulerian_rate_m_myr", std::to_string(m.eulerian));
    EXPECT_GT(m.eulerian, 4);
    EXPECT_LT(m.measured, .05 * m.eulerian);
}

// A tectonic change on a moving plate is measured as the change alone. On the same ramp moving
// east by 0.25 column, a plume adds 330 m of basalt under column 9 of every row: 40 m on a land
// column, so the step rate there is 40 m/Myr and r = 40 (1 - e^(-1/5)) = 7.25, not the Eulerian
// 40 - 25 = 15 m/Myr; every other cell stays near 0.
TEST(TectonicColumns, RockUpliftRateOnAMovingPlateIsTheTectonicChange) {
    tectonics::margins::Parameters still;
    still.root_mobility = 0;
    const Band band(4, cubefixture::radiusfor_km(4, 100000));
    const auto cells = cubefixture::layout(band, seed(500, true, 100), [](int x, int) {
        return seed(500 + 100. * x, true, 100); });
    tectonics::plumes::Source plume;
    plume.basalt_m.assign(cells.size(), 0);
    plume.influence.assign(cells.size(), 0);
    for (int row = 0; row < band.rows(); ++row) plume.basalt_m[band.cell(9, row)] = 330;
    Model world;
    world.initialize(band.grid, cells);
    world.advance(cubefixture::rotating(*band.grid, band.eastward(.25)), 1, 0, 1, 1, nullptr, still, {}, false, &plume, 5);
    const double expected = 40 * -std::expm1(-1.0 / 5);
    double largest_error = 0;
    for (int y = 0; y < band.rows(); ++y) for (int x = 0; x < band.columns(); ++x) {
        const double rate = world.upliftRate()[band.cell(x, y)];
        largest_error = std::max(largest_error, std::abs(rate - (x == 9 ? expected : 0)));
    }
    RecordProperty("largest_rate_error_m_myr", std::to_string(largest_error));
    EXPECT_LT(largest_error, .05 * expected);
}

// Rigid rotation of a curved surface (no other process): land columns
// z = 1000 + 500 sin(2 pi x / 16) + 300 cos(2 pi y / 4) on the band rotate about a tilted
// axis, at subdivision 1 and at the baseline's 4. The carried elevation moves with the engine's
// own transport, so the mean is far below the fixed-cell rate.
TEST(TectonicColumns, RockUpliftRateIsZeroForPlateMotionOverACurvedSurface) {
    tectonics::margins::Parameters still;
    still.root_mobility = 0;
    const Band band(4, cubefixture::radiusfor_km(4, 100000));
    const auto cells = cubefixture::layout(band, seed(1000, true, 100), [](int x, int y) {
        return seed(1000 + 500 * std::sin(2 * pi * x / 16) + 300 * std::cos(2 * pi * y / 4), true, 100); });
    const double column = pi / 8;
    for (uint32_t subdivision : {1u, 4u}) {
        Model world;
        world.initialize(band.grid, cells, subdivision);
        const auto m = step_rates(world, band, cubefixture::rotating(*band.grid, {.06 * column, -.08 * column, .16 * column}), still);
        RecordProperty("largest_rate_s" + std::to_string(subdivision), std::to_string(m.measured));
        RecordProperty("largest_eulerian_s" + std::to_string(subdivision), std::to_string(m.eulerian));
        EXPECT_GT(m.eulerian, 10) << subdivision;
        EXPECT_LT(m.measured, .05 * m.eulerian) << subdivision;
    }
}

// The running mean follows the rock across a face edge and over the pole. Uniform land 500 m
// rotates by a quarter column per 1 Myr step; a plume adds 330 m of basalt (40 m) under one cell
// in the first step only: r1 = 40 w there (w = 1 - e^(-1/5)). The second step has no tectonic
// change in the rock's frame, so its mean is e^(-1/5) r1 spread where the rock went: the
// integral of the rate keeps that share of r1 times the cell's area, its centre moves with the
// rotation, and cells away from the path keep about 0.
TEST(TectonicColumns, RockUpliftRateMeanFollowsTheRockAcrossAFaceEdgeAndThePole) {
    tectonics::margins::Parameters still;
    still.root_mobility = 0;
    constexpr int n = 8;
    const double column = pi / (2.0 * n);
    const auto grid = cubefixture::cube(n, cubefixture::radiusfor_km(n, 100000));
    const auto& mesh = *grid->mesh();
    struct Case { const char* name; Vec3 start; Vec3 axis; };
    const Case cases[] = {
        {"face_edge", mesh::cubepoint(0, pi / 4 - column / 2, 0), {0, 0, 1}},
        {"pole", mesh::unitvector(pi / 2 - column / 2, -pi / 2), {1, 0, 0}},
    };
    for (uint32_t subdivision : {1u, 4u}) for (const auto& test : cases) {
        const CellId source = mesh.locate(test.start);
        tectonics::plumes::Source plume;
        plume.basalt_m.assign(grid->cellcount(), 0);
        plume.influence.assign(grid->cellcount(), 0);
        plume.basalt_m[source] = 330;
        Model world;
        world.initialize(grid, std::vector<Cell>(grid->cellcount(), seed(500, true, 100)), subdivision);
        const Vec3 rotation = .25 * column * test.axis;
        const auto motion = cubefixture::rotating(*grid, rotation);
        world.advance(motion, 1, 0, 1, 1, nullptr, still, {}, false, &plume, 5);
        const double r1 = 40 * -std::expm1(-1.0 / 5);
        EXPECT_NEAR(world.upliftRate()[source], r1, .05 * r1) << test.name << " s" << subdivision;
        world.advance(motion, 1, 0, 1, 1, nullptr, still, {}, false, nullptr, 5);
        const double carried = r1 * std::exp(-1.0 / 5);
        Vec3 centre;
        double integral = 0, distant = 0;
        const Vec3 p = grid->position(source);
        for (CellId c = 0; c < grid->cellcount(); ++c) {
            const double rate = world.upliftRate()[c];
            integral += rate * grid->area_m2(c);
            centre = centre + (rate * grid->area_m2(c)) * grid->position(c);
            if (mesh::greatcircle_rad(p, grid->position(c)) > 3 * column) distant = std::max(distant, std::abs(rate));
        }
        const std::string label = std::string(test.name) + "_s" + std::to_string(subdivision);
        EXPECT_NEAR(integral / (carried * grid->area_m2(source)), 1, .02) << label;
        const double angle = mesh::length(rotation);
        const Vec3 k = test.axis;
        const Vec3 moved = std::cos(angle) * p + std::sin(angle) * mesh::cross(k, p) + (mesh::dot(k, p) * (1 - std::cos(angle))) * k;
        // The rate of the first step sat in the source cell (measured where the rock is at the
        // step's end): the second step's rotation moves it from the source cell's centre.
        const double error = mesh::greatcircle_rad(mesh::normalise(centre), moved) / column;
        RecordProperty(label + "_centre_error_columns", std::to_string(error));
        EXPECT_LT(error, .1) << label;
        // Away from the path only the transport's per-plate total-restoring scale (|scale - 1| of
        // order 1e-7 on the shifted elevation) reads as a rate.
        RecordProperty(label + "_distant_over_carried", std::to_string(distant / carried));
        EXPECT_LT(distant, 1e-4 * carried) << label;
    }
}

// Diverging plates at subdivision 4: land z = 1000 + 500 sin(2 pi x / 16) (range R = 1000 m)
// on the band, the halves part at -135 degrees (and close at 45) by 0.2 column per 1 Myr step;
// the opening stretches the boundary columns (real subsidence). With a window of 1e-3 Myr the
// mean is the step rate, so the carried source elevation is z_after - rate dt. Over six steps
// it stays within [low - 2R, high + 2R] of the elevations before each step; off the
// boundaries the rate stays small (the transport of one varying column field).
TEST(TectonicColumns, RockUpliftRateAtDivergentBoundariesIsBounded) {
    tectonics::margins::Parameters still;
    still.root_mobility = 0;
    const Band band(4, cubefixture::radiusfor_km(4, 100000));
    const auto cells = cubefixture::layout(band, seed(1000, true, 100), [](int x, int) {
        return seed(1000 + 500 * std::sin(2 * pi * x / 16), true, 100); });
    // Plate 0 west of 45 degrees moves west, plate 1 east: they part at 45 and close at -135.
    const auto motion = cubefixture::plates(halves(*band.grid), {band.eastward(-.2), band.eastward(.2)});
    Model world;
    world.initialize(band.grid, cells, 4);
    double interior = 0, lowest = 1e300, highest = -1e300, boundary_rate = 0;
    for (int step = 0; step < 6; ++step) {
        std::vector<double> before(cells.size());
        for (size_t i = 0; i < cells.size(); ++i) before[i] = elevation_m(world.cells()[i]);
        const auto [low, high] = std::minmax_element(before.begin(), before.end());
        const double range = *high - *low;
        world.advance(motion, 1, 0, 1, 1, nullptr, still, {}, false, nullptr, 1e-3);
        for (int y = 0; y < band.rows(); ++y) for (int x = 0; x < band.columns(); ++x) {
            const CellId i = band.cell(x, y);
            const double rate = world.upliftRate()[i];
            const double source = elevation_m(world.cells()[i]) - rate;
            EXPECT_GE(source, *low - 2 * range) << "step " << step << " column " << x;
            EXPECT_LE(source, *high + 2 * range) << "step " << step << " column " << x;
            lowest = std::min(lowest, source - *low);
            highest = std::max(highest, source - *high);
            // The boundaries lie between columns 7 and 8 (45 degrees) and 15 and 0 (-135).
            if (x == 0 || x == 7 || x == 8 || x == 15) boundary_rate = std::max(boundary_rate, std::abs(rate));
            else if ((x >= 2 && x <= 5) || (x >= 10 && x <= 13)) interior = std::max(interior, std::abs(rate));
        }
    }
    RecordProperty("interior_rate_m_myr", std::to_string(interior));
    RecordProperty("boundary_rate_m_myr", std::to_string(boundary_rate));
    EXPECT_LT(interior, .05 * boundary_rate);
    EXPECT_GT(boundary_rate, 0);
    EXPECT_GT(lowest, -2000); EXPECT_LT(highest, 1000);
}

// Thermal subsidence of an ageing ocean floor is tectonic: a 20 Myr old ocean column
// (thermal strength 1, nothing else acting) ages to 21 Myr in one step and sinks by
// subsidence_m(21) - subsidence_m(20) = 39.50061759809455 m (Parsons & Sclater law
// of tectonic_thermal.hpp), so r = -39.50061759809455 (1 - e^(-1/5)) = -7.160247204961215.
TEST(TectonicColumns, RockUpliftRateCountsThermalSubsidence) {
    const auto grid = cubefixture::cube(4, cubefixture::radiusfor_km(4, 100000));
    Model world;
    world.initialize(grid, std::vector<Cell>(grid->cellcount(), seed(-5000, false, 20)));
    world.advance(cubefixture::still(*grid), 1, 0, 1, 1, nullptr, {}, {}, false, nullptr, 5);
    const double sinking = tectonics::thermal::subsidence_m(21) - tectonics::thermal::subsidence_m(20);
    EXPECT_NEAR(sinking, 39.50061759809455, 1e-9);
    for (double rate : world.upliftRate()) {
        EXPECT_NEAR(rate, -sinking * -std::expm1(-1.0 / 5), 1e-9);
        EXPECT_NEAR(rate, -7.160247204961215, 1e-9);
    }
}

// Cumulative shortening (multiscale plan S5.1): the export of a strain tensor. The
// shortening is minus the least principal component times the cell width; the
// direction is that component's axis in [0, pi) from x toward y. Simple shear xy = g/2
// has components +-g/2 and shortens along (1, -1) (3 pi / 4) for g > 0, along (1, 1)
// for g < 0; extension alone shortens nothing; equal components have direction 0.
TEST(TectonicColumns, ShorteningAndDirectionOfExplicitStrains) {
    EXPECT_NEAR(shortening_m({-.1, 0, 0}, 1000), 100, 1e-12);
    EXPECT_EQ(shortening_direction({-.1, 0, 0}), 0);
    EXPECT_NEAR(shortening_m({0, -.1, 0}, 1000), 100, 1e-12);
    EXPECT_NEAR(shortening_direction({0, -.1, 0}), pi / 2, 1e-15);
    EXPECT_NEAR(shortening_m({0, 0, .05}, 1000), 50, 1e-12);
    EXPECT_NEAR(shortening_direction({0, 0, .05}), 3 * pi / 4, 1e-15);
    EXPECT_NEAR(shortening_direction({0, 0, -.05}), pi / 4, 1e-15);
    EXPECT_EQ(shortening_m({.1, .2, 0}, 1000), 0);
    EXPECT_NEAR(shortening_direction({.1, .2, 0}), 0, 1e-15); // the least extension's axis
    EXPECT_NEAR(shortening_m({-.1, -.1, 0}, 1000), 100, 1e-12);
    EXPECT_EQ(shortening_direction({-.1, -.1, 0}), 0);
    EXPECT_EQ(shortening_m({}, 1000), 0);
    EXPECT_EQ(shortening_direction({}), 0);
}

namespace {
// Land columns at rest except for the motion: no erosion, no roots, no forcing; a band of N 8
// with 100 km columns at the equator.
struct Strained {
    Band band{8, cubefixture::radiusfor_km(8, 100000)};
    Model world;
    explicit Strained(uint32_t subdivision, double elevation = 500, bool continental = true) {
        world.initialize(band.grid, std::vector<Cell>(band.cellcount(), seed(elevation, continental, 100)), subdivision);
    }
    void advance(const st::Motion& motion, const tectonics::margins::Forcing* forcing = nullptr) {
        tectonics::margins::Parameters still;
        still.root_mobility = 0;
        world.advance(motion, 1, 0, 1, 1, forcing, still, {}, false, nullptr, 5);
    }
    // Offsets east (u, column widths at the equator, times the cosine of the cell's latitude:
    // a uniform u is a rotation about the axis, no strain) and north (v, column widths) at the
    // band's cells, zero elsewhere.
    st::Motion field(const std::function<double(int, int)>& u, const std::function<double(int, int)>& v,
                     const std::vector<uint32_t>& owners = {}) const {
        std::vector<Vec3> offset(band.cellcount());
        for (int y = 0; y < band.rows(); ++y) for (int x = 0; x < band.columns(); ++x) {
            const CellId c = band.cell(x, y);
            const Vec3 r = band.grid->position(c);
            const double parallel = std::sqrt(r.x * r.x + r.y * r.y);
            offset[c] = (u(x, y) * band.spacing_m() * parallel) * cubefixture::east(*band.grid, c)
                + (v(x, y) * band.spacing_m()) * cubefixture::north(*band.grid, c);
        }
        return cubefixture::offsets(*band.grid, offset, owners);
    }
    double shortening(int x, int y) const { return shortening_m(world.strain()[band.cell(x, y)], band.spacing_m()); }
    const Strain& strain(int x, int y) const { return world.strain()[band.cell(x, y)]; }
};
// Periodic triangle profile of period 2 n: v at 0, -v at n, slope -v / (n / 2) between.
double triangle(int x, int n, double v) {
    const int p = ((x % (2 * n)) + 2 * n) % (2 * n);
    return p <= n ? v - 2 * v * p / n : -v + 2 * v * (p - n) / n;
}
// The two rows that straddle the equator.
constexpr int equator[2] = {3, 4};
}

// A uniform convergence along the band: u_x a periodic triangle of period 16 columns (0.2 at
// x = 4, -0.2 at x = 12, slope -0.05 column per column between), so each update shortens the
// rock of x = 5..11 by about 0.05 along x near the equator: 5 km per 100 km column, direction
// about 0 (the cells' axes are their east neighbours'). Rock arriving from the kinks mixes into
// x = 5 and 11, so x = 6..10 after one update, x = 7..9 after two (10 km). Subdivisions 1 and 4.
TEST(TectonicColumns, ShorteningOfAUniformConvergence) {
    for (uint32_t subdivision : {1u, 4u}) {
        Strained s(subdivision);
        const auto u = s.field([](int x, int) { return triangle(x - 4, 8, .2); }, [](int, int) { return 0.0; });
        s.advance(u);
        for (const int y : equator) for (int x = 6; x <= 10; ++x) {
            EXPECT_NEAR(s.shortening(x, y), 5000, 300) << "subdivision " << subdivision << " x " << x;
            EXPECT_NEAR(std::sin(shortening_direction(s.strain(x, y))), 0, .05);  // an axis: 0 or pi
        }
        for (const int y : equator) for (int x : {14, 15, 0, 1, 2})
            EXPECT_LT(s.shortening(x, y), 50) << "x " << x;
        s.advance(u);
        for (const int y : equator) for (int x = 7; x <= 9; ++x)
            EXPECT_NEAR(s.shortening(x, y), 10000, 600) << "subdivision " << subdivision << " x " << x;
    }
}

// The tensor moves with the rock. u_x = 0.1 for x <= 4, 0 at x = 5, -0.1 for 6 <= x <= 10
// (rest 0) strains column 5; a second update rotates everything by a quarter column east (no
// strain, no spin): the total over the equatorial rows is kept and its centre moves east by
// the rotation, the direction staying about 0.
TEST(TectonicColumns, ShorteningIsCarriedWithTheRock) {
    Strained s(1);
    s.advance(s.field([](int x, int) { return x <= 4 ? .1 : x == 5 ? 0.0 : x <= 10 ? -.1 : 0.0; },
                      [](int, int) { return 0.0; }));
    const auto rows = [&] {
        double total = 0, moment = 0;
        for (const int y : equator) for (int x = 0; x < 16; ++x) {
            const double amount = s.shortening(x, y);
            total += amount; moment += amount * x;
        }
        return std::pair{total, moment / total};
    };
    const auto [first, before] = rows();
    EXPECT_GT(s.shortening(5, 3), 8000);
    s.advance(cubefixture::rotating(*s.band.grid, s.band.eastward(.25)));
    const auto [total, centroid] = rows();
    EXPECT_NEAR(total / first, 1, .02);
    EXPECT_NEAR(centroid, before + .25, .03);
    for (int x = 0; x < 16; ++x) for (const int y : equator)
        if (s.shortening(x, y) > 1000) EXPECT_NEAR(std::sin(shortening_direction(s.strain(x, y))), 0, .1);
}

// The tensor turns with the plate. After one update of the uniform convergence along x, a
// rigid rotation by w = 0.05 rad about the outward axis at (longitude 0, latitude 0) turned
// from east toward south (about -x) adds no strain and turns the axis by w there: the cells
// next to longitude 0 on the equatorial rows keep their shortening at direction about w.
TEST(TectonicColumns, ShorteningTurnsWithAPlateRotation) {
    const double w = .05;
    for (uint32_t subdivision : {1u, 4u}) {
        Strained s(subdivision);
        // Longitude 0 lies between band columns 11 and 12 (N 8: face 0 is columns 8..15).
        s.advance(s.field([](int x, int) { return triangle(x - 8, 8, .2); }, [](int, int) { return 0.0; }));
        std::vector<double> before, direction;
        for (const int y : equator) for (int x : {11, 12}) {
            before.push_back(s.shortening(x, y));
            direction.push_back(shortening_direction(s.strain(x, y)));
        }
        s.advance(cubefixture::rotating(*s.band.grid, {-w, 0, 0}));
        size_t k = 0;
        for (const int y : equator) for (int x : {11, 12}) {
            EXPECT_NEAR(s.shortening(x, y) / before[k], 1, .1) << "subdivision " << subdivision << " x " << x;
            // The axis turned by w (directions are axes, modulo pi).
            const double turned = std::remainder(shortening_direction(s.strain(x, y)) - direction[k], pi);
            EXPECT_NEAR(turned, w, .02) << "subdivision " << subdivision << " x " << x;
            ++k;
        }
    }
}

// Plate edges: a jump of velocity to another plate is not strain. Plate 0 (x < 8) moves east by
// 0.2 column and slows to 0.1 at x = 6 and 0 at x = 7; plate 1 (x >= 8) moves east by 0.1. With
// the plates, x = 6 and 7 take the in-plate difference 0.1 (10 km), plate 1 and the rest of
// plate 0 nothing, though the boundary opens between x = 7 and 8. As one plate, x = 7's
// gradient is the central difference 0: only x = 6 is strained (and part of its rock reaches
// x = 7). Rigid plates converging (rotations) shorten nothing.
TEST(TectonicColumns, ShorteningIgnoresJumpsBetweenPlates) {
    Strained owned(1), single(1), converging(1);
    std::vector<uint32_t> owners(owned.band.cellcount(), 1);
    for (int y = 0; y < owned.band.rows(); ++y) for (int x = 0; x < 8; ++x) owners[owned.band.cell(x, y)] = 0;
    const auto u = [](int x, int) { return x < 6 ? .2 : x == 6 ? .1 : x == 7 ? 0.0 : .1; };
    const auto none = [](int, int) { return 0.0; };
    // The plates' common parts as their rotations (0.2 and 0.1 column), the rest as offsets.
    auto owned_motion = owned.field([&](int x, int y) { return u(x, y) - (x < 8 ? .2 : .1); }, none, owners);
    owned_motion.rotation = {owned.band.eastward(.2), owned.band.eastward(.1)};
    owned.advance(owned_motion);
    single.advance(single.field(u, none));
    converging.advance(cubefixture::plates(owners, {converging.band.eastward(.2), converging.band.eastward(-.2)}));
    for (const int y : equator) {
        EXPECT_NEAR(owned.shortening(7, y), 10000, 800) << "row " << y;
        EXPECT_NEAR(owned.shortening(6, y), 10000 * 9. / 11, 1500) << "row " << y;
        for (int x : {0, 1, 2, 3, 9, 10, 11, 12, 13, 14}) EXPECT_LT(owned.shortening(x, y), 100) << "x " << x;
        EXPECT_LT(single.shortening(7, y), 2500) << "row " << y;
        EXPECT_GT(single.shortening(6, y), 5000) << "row " << y;
        for (int x = 0; x < 32; ++x) EXPECT_LT(converging.shortening(x, y), 1) << "x " << x;
    }
}

// New crust has no strain. Ocean floor is first shortened along the rows (the triangle profile
// of period 8 rows, one plate) across every column; then the two halves part for 50 updates at
// subdivision 4 (0.2 column per update, a ridge at 45 degrees): the crust the ridge creates
// carries none and the rock that stayed is diluted by it, so the ridge columns keep almost none
// of the 5 km.
TEST(TectonicColumns, ShorteningOfNewCrustAtADivergentBoundaryIsZero) {
    Strained s(4, -4000, false);
    const auto owners = halves(*s.band.grid);
    s.advance(s.field([](int, int) { return 0.0; }, [](int, int y) { return -triangle(y - 2, 4, .1); }));
    double squeezed = 0;
    for (int x = 0; x < s.band.columns(); ++x) squeezed = std::max(squeezed, s.shortening(x, 4));
    ASSERT_GT(squeezed, 2000);
    const auto apart = cubefixture::plates(owners, {s.band.eastward(-.2), s.band.eastward(.2)});
    tectonics::margins::Forcing plates;
    plates.thermal_age_myr.assign(owners.size(), -1);
    plates.trench_m.assign(owners.size(), 0);
    plates.arc_receiver.assign(owners.size(), std::numeric_limits<uint32_t>::max());
    plates.owners = owners;
    for (int step = 0; step < 50; ++step) s.advance(apart, &plates);
    // The ridge at 45 degrees lies between band columns 15 and 16 (N 8).
    double ridge = 0;
    for (int y = 0; y < s.band.rows(); ++y) for (int x : {15, 16}) ridge = std::max(ridge, s.shortening(x, y));
    RecordProperty("largest_ridge_shortening_m", std::to_string(ridge));
    EXPECT_LT(ridge, .02 * squeezed);
}

// Over several updates of a varied displacement (two plates, subdivision 4) the shortening
// stays finite and non-negative and the direction in [0, pi).
TEST(TectonicColumns, ShorteningIsNonNegativeAndItsDirectionInRange) {
    Strained s(4);
    std::vector<uint32_t> owners(s.band.cellcount(), 0);
    for (int y = 0; y < s.band.rows(); ++y) for (int x = 0; x < s.band.columns(); ++x)
        owners[s.band.cell(x, y)] = (x + y) % 32 < 16 ? 0 : 1;
    const auto u = s.field([](int x, int y) { return .3 * std::sin(2 * pi * x / 32 + y); },
                           [](int x, int y) { return .2 * std::cos(2 * pi * y / 8 - x / 3.0); }, owners);
    double largest = 0;
    for (int step = 0; step < 6; ++step) {
        s.advance(u);
        for (const auto& strain : s.world.strain()) {
            const double amount = shortening_m(strain, 100000), direction = shortening_direction(strain);
            ASSERT_TRUE(std::isfinite(amount));
            EXPECT_GE(amount, 0);
            EXPECT_GE(direction, 0);
            EXPECT_LT(direction, pi);
            largest = std::max(largest, amount);
        }
    }
    EXPECT_GT(largest, 0);
}

// The column processes on the cube (W2-E, E2B): finite volumes over the cell table.
namespace {
double total_continent_m3(const Metric& metric, const std::vector<Cell>& cells) {
    double total = 0;
    for (CellId c = 0; c < cells.size(); ++c) total += cells[c].continent * metric.area_m2(c);
    return total;
}

// One explicit diffusion step of a Gaussian bump (sigma 600 km, 1000 m on a 500 m continent)
// centred at a cube point C: the spread of its mass, the change of the integral of
// f = |P - C|² R² (the chord squared) times the mass per area, over the continuum's
// rate kappa rho_s ∫ (z - 500) ∇²f dA with ∇²f = 4 cos(angle to C). 1 for a consistent
// discretization (the gross face transfer, by contrast, scales with each cell's aspect).
double bumpspread(const std::shared_ptr<const TectonicGrid>& grid, const Metric& metric, const mesh::Vec3& centre) {
    const auto& mesh = *grid->mesh();
    const double radius_m = mesh.radius_m();
    std::vector<Cell> cells(grid->cellcount());
    std::vector<double> f(cells.size()), elevation(cells.size());
    double predicted = 0, before = 0;
    for (CellId c = 0; c < cells.size(); ++c) {
        const double angle = mesh::greatcircle_rad(centre, mesh.centre(c));
        const double distance_km = angle * radius_m / 1000;
        cells[c] = seed(500 + 1000 * std::exp(-.5 * distance_km * distance_km / (600.0 * 600.0)), true, 100);
        elevation[c] = elevation_m(cells[c]);
        f[c] = 2 * radius_m * radius_m * (1 - std::cos(angle));
        predicted += transport_diffusivity * sediment_density * 4 * std::cos(angle) * (elevation[c] - 500) * metric.area_m2(c);
        before += f[c] * mass_per_area(cells[c]) * metric.area_m2(c);
    }
    Budget budget;
    diffuse_surface(metric, cells, elevation, surface_diffusion_number(metric, 1, 1), 1, budget);
    EXPECT_NEAR(budget.deposited_kg / budget.eroded_kg, 1, 1e-12);
    double after = 0;
    for (CellId c = 0; c < cells.size(); ++c) after += f[c] * mass_per_area(cells[c]) * metric.area_m2(c);
    return (after - before) / predicted;
}
} // namespace

TEST(TectonicColumnsCube, DiffusionOfABumpIsEqualAtAFaceCentreEdgeAndCorner) {
    constexpr double quarter = 0.78539816339744830962;
    for (const int n : {32, 64}) {
        const auto grid = TectonicGrid::cube(n, earth_km);
        const Metric metric(*grid);
        // Each centre a cell vertex (N even): four cells at the face centre and the edge, three at the corner.
        const double centre = bumpspread(grid, metric, mesh::cubepoint(0, 0, 0));
        const double edge = bumpspread(grid, metric, mesh::cubepoint(0, quarter, 0));
        const double corner = bumpspread(grid, metric, mesh::cubepoint(0, quarter, quarter));
        const std::string suffix = "_n" + std::to_string(n);
        RecordProperty("diffusion_face_centre" + suffix, std::to_string(centre));
        RecordProperty("diffusion_face_edge" + suffix, std::to_string(edge));
        RecordProperty("diffusion_cube_corner" + suffix, std::to_string(corner));
        RecordProperty("reference_area_over_mean" + suffix, std::to_string(metric.reference_area_m2()
            / (grid->totalarea_m2() / static_cast<double>(grid->cellcount()))));
        RecordProperty("correction_factor" + suffix, std::to_string(metric.correctionfactor()));
        // With the non-orthogonal correction (without it the two-point weights spread 7 % fast
        // over a face edge and 11 % slow over a corner).
        EXPECT_NEAR(centre, 1, .02) << n;
        EXPECT_NEAR(edge / centre, 1, .02) << n;
        EXPECT_NEAR(corner / centre, 1, .02) << n;
    }
}

// The explicit limits of root flow, collapse and diffusion hold with the corrected differences:
// a checkerboard of thickness or elevation over a cube corner (the most skewed cells), stepped
// at the limits the processes choose, decays in amplitude and keeps every column positive.
TEST(TectonicColumnsCube, CheckerboardAtTheExplicitLimitsDecaysOverACubeCorner) {
    constexpr int n = 16;
    const auto grid = TectonicGrid::cube(n, cubefixture::radiusfor_km(n, 100000));
    const Metric metric(*grid);
    const auto& mesh = *grid->mesh();
    const Vec3 corner = mesh::normalise({1, 1, 1});
    EXPECT_GT(metric.correctionfactor(), 1);
    const auto checker = [&](double low, double high) {
        std::vector<Cell> cells(grid->cellcount());
        for (CellId c = 0; c < cells.size(); ++c) {
            int face = 0, i = 0, j = 0;
            mesh::cubeindices(n, c, face, i, j);
            const bool close = mesh::greatcircle_rad(corner, mesh.centre(c)) * grid->radius_m() < 6e5;
            cells[c] = seed(close && (i + j) % 2 ? high : low, true, 100);
        }
        return cells;
    };
    const auto amplitude = [&](const std::vector<Cell>& cells, double Cell::* member) {
        double lo = 1e300, hi = -1e300;
        for (CellId c = 0; c < cells.size(); ++c) {
            if (mesh::greatcircle_rad(corner, mesh.centre(c)) * grid->radius_m() >= 6e5) continue;
            lo = std::min(lo, cells[c].*member); hi = std::max(hi, cells[c].*member);
        }
        return hi - lo;
    };
    // Diffusion at its limit (0.2 / correction factor).
    {
        auto cells = checker(0, 3000);
        std::vector<double> elevation(cells.size());
        for (size_t c = 0; c < cells.size(); ++c) elevation[c] = elevation_m(cells[c]);
        const double start = *std::max_element(elevation.begin(), elevation.end()) - *std::min_element(elevation.begin(), elevation.end());
        Budget budget;
        for (int step = 0; step < 20; ++step) {
            for (size_t c = 0; c < cells.size(); ++c) elevation[c] = elevation_m(cells[c]);
            diffuse_surface(metric, cells, elevation, .2 / metric.correctionfactor(), 1, budget);
            for (const auto& c : cells) { ASSERT_GE(c.continent, 0); ASSERT_GE(c.sediment, 0); }
        }
        for (size_t c = 0; c < cells.size(); ++c) elevation[c] = elevation_m(cells[c]);
        const double end = *std::max_element(elevation.begin(), elevation.end()) - *std::min_element(elevation.begin(), elevation.end());
        EXPECT_LT(end, start);
    }
    // Root flow and collapse choose their own sub-steps.
    {
        auto cells = checker(0, 8000);
        const double start = amplitude(cells, &Cell::continent), volume = total_continent_m3(metric, cells);
        Budget budget;
        for (int step = 0; step < 5; ++step) spread_roots(metric, cells, 1, Isostasy{}, false, budget);
        EXPECT_LT(amplitude(cells, &Cell::continent), start);
        EXPECT_NEAR(total_continent_m3(metric, cells) / volume, 1, 1e-12);
        for (const auto& c : cells) EXPECT_GT(c.continent, 0);
    }
    {
        auto cells = checker(0, 9000);  // 35 km against about 95 km of crust
        tectonics::margins::Parameters p;
        p.enabled = true;
        p.continental_max_thickness_m = 40000;
        p.orogenic_collapse_km2_myr = 30000;
        const double start = amplitude(cells, &Cell::continent), volume = total_continent_m3(metric, cells);
        Budget budget;
        for (int step = 0; step < 5; ++step) collapse_orogens(metric, cells, 1, p, budget);
        EXPECT_LT(amplitude(cells, &Cell::continent), start);
        EXPECT_NEAR(total_continent_m3(metric, cells) / volume, 1, 1e-12);
        for (const auto& c : cells) EXPECT_GT(c.continent, 0);
    }
}

TEST(TectonicColumnsCube, CollapseAndRootFlowConserveVolumeAcrossFaceEdgesAndCorners) {
    const auto grid = TectonicGrid::cube(16, earth_km);
    const Metric metric(*grid);
    const auto& mesh = *grid->mesh();
    constexpr double quarter = 0.78539816339744830962;
    // An orogen of 80 km crust over a cube corner and along a face edge, on a continent at sea level.
    const mesh::Vec3 corner = mesh::cubepoint(0, quarter, quarter), edge = mesh::cubepoint(0, quarter, 0);
    std::vector<Cell> cells(grid->cellcount());
    for (CellId c = 0; c < cells.size(); ++c) {
        const double closest = std::min(mesh::greatcircle_rad(corner, mesh.centre(c)), mesh::greatcircle_rad(edge, mesh.centre(c)));
        cells[c] = seed(0, true, 100);
        if (closest * earth_km < 700) {
            cells[c].continent = 80000;
            cells[c].continental_age_volume = 80000 * 300;
        }
    }
    tectonics::margins::Parameters p;
    p.enabled = true;
    p.continental_max_thickness_m = 60000;
    p.orogenic_collapse_km2_myr = 30000;
    const double volume = total_continent_m3(metric, cells);
    double age = 0;
    for (CellId c = 0; c < cells.size(); ++c) age += cells[c].continental_age_volume * metric.area_m2(c);
    Budget budget;
    for (int step = 0; step < 10; ++step) collapse_orogens(metric, cells, 1, p, budget);
    EXPECT_GT(budget.collapsed_continent_kg, 0);
    EXPECT_NEAR(total_continent_m3(metric, cells) / volume, 1, 1e-12);
    double age_after = 0, thickest = 0;
    for (CellId c = 0; c < cells.size(); ++c) {
        age_after += cells[c].continental_age_volume * metric.area_m2(c);
        thickest = std::max(thickest, cells[c].continent);
    }
    EXPECT_NEAR(age_after / age, 1, 1e-12);
    EXPECT_LT(thickest, 80000);
    const Isostasy isostasy;
    for (int step = 0; step < 5; ++step) spread_roots(metric, cells, 1, isostasy, false, budget);
    EXPECT_GT(budget.redistributed_continent_kg, 0);
    EXPECT_NEAR(total_continent_m3(metric, cells) / volume, 1, 1e-12);
}

TEST(TectonicColumnsCube, PlumeBasaltConservesItsVolumeOnTheSphere) {
    const auto grid = TectonicGrid::cube(32, earth_km);
    tectonics::plumes::Parameters p;
    p.count = 6;
    p.head_volume_km3 = 2e6;
    const auto source = tectonics::plumes::source(*grid, 1, 0, 2, p);
    const auto plumes = tectonics::plumes::placesphere(earth_km, 1, p);
    ASSERT_EQ(plumes.size(), 6u);
    for (size_t a = 0; a < plumes.size(); ++a) for (size_t b = 0; b < a; ++b)
        EXPECT_GE(mesh::greatcircle_rad(plumes[a].position, plumes[b].position) * earth_km, p.min_spacing_km);
    double expected = 0;
    for (const auto& plume : plumes) {
        const auto smooth = [](double t) { t = std::clamp(t, 0.0, 1.0); return t * t * (3 - 2 * t); };
        expected += p.tail_flux_km3_yr * 1e15 * std::max(0.0, 2 - plume.birth_myr)
            + p.head_volume_km3 * 1e9 * (smooth((2 - plume.birth_myr) / p.head_duration_myr) - smooth(-plume.birth_myr / p.head_duration_myr));
    }
    double volume = 0;
    for (CellId c = 0; c < grid->cellcount(); ++c) volume += source.basalt_m[c] * grid->area_m2(c);
    EXPECT_NEAR(volume / expected, 1, 1e-12);
    for (const auto& plume : plumes) {
        EXPECT_GT(source.influence[grid->mesh()->locate(plume.position)], .8);
    }
}

TEST(TectonicColumnsCube, InitialTerrainRunsOnTheCube) {
    const auto grid = TectonicGrid::cube(16, earth_km);
    const Metric metric(*grid);
    const auto& mesh = *grid->mesh();
    const size_t count = grid->cellcount();
    std::vector<double> elevations(count), ages(count, 50);
    std::vector<uint32_t> owners(count);
    std::vector<mesh::Vec3> velocity(count);
    const mesh::Vec3 pole{0, 0, 1};
    for (CellId c = 0; c < count; ++c) {
        const mesh::Vec3 r = mesh.centre(c);
        owners[c] = r.x > 0 ? 1 : 0;
        elevations[c] = r.y > .3 ? 400 : -4000;
        // Each plate turns about the pole, the two apart (spreading along x = 0).
        velocity[c] = (owners[c] ? 20000.0 : -20000.0) * mesh::cross(pole, r);
    }
    const auto cells = tectonics::initialization::terrain(metric, 1, elevations, ages, owners, velocity, 1, 0, 1);
    EXPECT_THROW(tectonics::initialization::terrain(metric, 1, elevations, ages, owners,
        std::vector<Displacement>(count - 1), 1, 0, 1), std::runtime_error);
    ASSERT_EQ(cells.size(), count);
    for (const auto& c : cells) EXPECT_TRUE(std::isfinite(elevation_m(c)));
}

// Explicit rivers and waves on the cube (routeSediment): a continental cap over a cube corner
// drains into the surrounding ocean; eroded and deposited mass agree in true areas, sediment
// reaches the ocean beside the coast and the mass budget closes. Implicit rivers (the landscape
// kernels on the cube's metric) do the same.
TEST(TectonicColumnsCube, RiversCloseTheirBudgetOverACubeCorner) {
    constexpr int n = 16;
    const auto grid = TectonicGrid::cube(n, cubefixture::radiusfor_km(n, 150000));
    const auto& mesh = *grid->mesh();
    const Vec3 corner = mesh::normalise({1, 1, 1});
    std::vector<Cell> cells(grid->cellcount());
    std::vector<uint8_t> coast(cells.size(), 0);
    for (CellId c = 0; c < cells.size(); ++c) {
        const double distance = mesh::greatcircle_rad(corner, mesh.centre(c)) * grid->radius_m() / 150000;
        cells[c] = distance < 6 ? seed(2000 - 300 * distance, true, 40) : seed(-4000, false, 40);
        coast[c] = distance >= 6 && distance < 7.5;
    }
    for (const bool implicit : {false, true}) {
        tectonics::margins::Parameters p;
        p.enabled = true;
        p.river_implicit = implicit;
        p.wave_erosion_m_myr = implicit ? 0 : 20;
        Model model;
        model.initialize(grid, cells);
        model.advance(cubefixture::still(*grid), 1, 0, 1, 0, nullptr, p);
        const auto& b = model.budget();
        EXPECT_GT(b.river_eroded_kg, 0) << implicit;
        EXPECT_NEAR(b.eroded_kg, b.deposited_kg, b.eroded_kg * 1e-12) << implicit;
        EXPECT_LT(std::abs(b.relative_mass_residual), 1e-12) << implicit;
        double offshore = 0;
        for (CellId c = 0; c < cells.size(); ++c) if (coast[c]) offshore += model.cells()[c].sediment * grid->area_m2(c);
        EXPECT_GT(offshore, 0) << implicit;
        RecordProperty(implicit ? "implicit_river_kg" : "explicit_river_kg", std::to_string(b.river_eroded_kg));
    }
}

// Arc melt is emplaced as a Gaussian over the receiver plate's cells within 300 km (great
// circle) of its receiver: an ocean plate closing on a continental one (margins on, the forcing
// of margins::forcing) adds arc basalt only on the receivers' plate near a receiver, and the
// added volume is the budget's generated arc mass.
TEST(TectonicColumnsCube, ArcMeltSpreadsWithinThreeHundredKilometresOfItsReceiver) {
    const Band band(16, cubefixture::radiusfor_km(16, 100000));
    const auto owners = halves(*band.grid);
    std::vector<Cell> cells(band.cellcount());
    for (CellId c = 0; c < cells.size(); ++c) cells[c] = owners[c] == 0 ? seed(-4000, false, 80) : seed(500, true, 100);
    tectonics::margins::Parameters p;
    p.enabled = true;
    p.river_erosion = 0;
    const auto motion = cubefixture::plates(owners, {band.eastward(.3), {}});
    std::vector<uint32_t> stable(owners.begin(), owners.end());
    const auto rigid = Model::displacement(*band.grid, motion);
    Model world, control;
    world.initialize(band.grid, cells);
    control.initialize(band.grid, cells);
    tectonics::margins::Polarity polarity;
    for (int step = 0; step < 6; ++step) {
        const auto forcing = tectonics::margins::forcing(world.metric(), 1, stable, world.cells(), rigid, p, &polarity);
        auto silent = p;
        silent.arc_melt_ratio = 0;
        control = world;
        world.advance(motion, 1, 0, 1, 0, &forcing, p);
        control.advance(motion, 1, 0, 1, 0, &forcing, silent);
        if (world.budget().generated_arc_kg == 0) continue;
        double added = 0;
        const auto& mesh = *band.grid->mesh();
        std::vector<CellId> receivers;
        for (CellId c = 0; c < cells.size(); ++c)
            if (forcing.arc_receiver[c] != std::numeric_limits<uint32_t>::max()) receivers.push_back(forcing.arc_receiver[c]);
        std::sort(receivers.begin(), receivers.end());
        receivers.erase(std::unique(receivers.begin(), receivers.end()), receivers.end());
        for (CellId c = 0; c < cells.size(); ++c) {
            const double extra = world.cells()[c].ocean - control.cells()[c].ocean;
            if (std::abs(extra) < 1e-9) continue;
            added += extra * band.grid->area_m2(c) * oceanic_density;
            double nearest = 1e300;
            for (const CellId r : receivers) nearest = std::min(nearest, mesh::greatcircle_m(mesh.centre(r), mesh.centre(c), mesh.radius_m()));
            EXPECT_LE(nearest, 300000 + 1) << c;
            EXPECT_EQ(owners[c], 1u) << c;
        }
        const double generated = world.budget().generated_arc_kg - control.budget().generated_arc_kg;
        EXPECT_NEAR(added / generated, 1, 1e-9);
        return;
    }
    FAIL() << "no arc melt in six updates";
}
