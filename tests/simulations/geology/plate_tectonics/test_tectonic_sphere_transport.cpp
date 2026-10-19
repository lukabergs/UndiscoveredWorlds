// The rock transport on the cube (W2-E, step E2D): the subcell lattice, conservation and
// positivity under plate rotations, crossings of a face edge, a cube corner and a pole, the
// coverage of a rigid plate, and (disabled, run once) its cost.
#include "gtest/gtest.h"

#include "core/mesh/cell_grid.hpp"
#include "core/mesh/cube_geometry.hpp"
#include "simulations/geology/plate_tectonics/tectonic_sphere_transport.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace st = tectonics::sphere_transport;
using mesh::Vec3;

namespace {
constexpr double pi = 3.14159265358979323846;
constexpr double radius_m = 6371e3;

// splitmix64, seed 1: portable pseudo-random numbers for fixtures.
struct Random {
    std::uint64_t state = 1;
    std::uint64_t next() {
        std::uint64_t z = (state += 0x9e3779b97f4a7c15ULL);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        return z ^ (z >> 31);
    }
    double uniform() { return static_cast<double>(next() >> 11) * (1.0 / 9007199254740992.0); }
    Vec3 direction() {
        const double z = 2 * uniform() - 1, lon = 2 * pi * uniform(), r = std::sqrt(1 - z * z);
        return {r * std::cos(lon), r * std::sin(lon), z};
    }
};

double total(const st::Lattice& l, const std::vector<double>& v) {
    double sum = 0;
    for (st::CellId c = 0; c < l.subcellcount(); ++c) sum += v[c] * l.area_m2(c);
    return sum;
}

Vec3 masscentre(const st::Lattice& l, const std::vector<double>& v) {
    Vec3 sum;
    for (st::CellId c = 0; c < l.subcellcount(); ++c) sum = sum + (v[c] * l.area_m2(c)) * l.centre(c);
    return mesh::normalise(sum);
}

Vec3 rotate(const Vec3& w, const Vec3& v) {
    const double angle = mesh::length(w);
    const Vec3 k = (1 / angle) * w;
    return std::cos(angle) * v + std::sin(angle) * mesh::cross(k, v) + (mesh::dot(k, v) * (1 - std::cos(angle))) * k;
}

st::Motion rigid(const st::Lattice& l, const Vec3& rotation) {
    st::Motion m;
    m.rotation = {rotation};
    m.plate.assign(l.cellcount(), 0);
    return m;
}

// A smooth cap of mass around p (a Gaussian of sigma subcell spacings, zero below 1e-8).
std::vector<double> cap(const st::Lattice& l, const Vec3& p, double sigma_cells) {
    const double sigma = sigma_cells * mesh::cubestep(l.subfacecells());
    std::vector<double> v(l.subcellcount(), 0);
    for (st::CellId c = 0; c < l.subcellcount(); ++c) {
        const double angle = std::acos(std::clamp(mesh::dot(l.centre(c), p), -1.0, 1.0));
        const double value = 1000 * std::exp(-0.5 * angle * angle / (sigma * sigma));
        v[c] = value > 1e-8 ? value : 0;
    }
    return v;
}

// Moves a cap from start to end along the great circle in a number of updates; checks mass,
// positivity, the mass centre against the rotated start and the share on the end's face.
void crossing(const Vec3& start, const Vec3& end, const char* name) {
    const auto lattice = st::Lattice::make(8, 2, radius_m);
    const st::Lattice& l = *lattice;
    const double angle = std::acos(mesh::dot(start, end));
    const int updates = 12;
    const Vec3 rotation = (angle / updates) * mesh::normalise(mesh::cross(start, end));
    const st::Stencil stencil(lattice, rigid(l, rotation));
    auto field = cap(l, start, 1.5);
    const double mass = total(l, field);
    for (int t = 0; t < updates; ++t) field = stencil.remap(field);
    EXPECT_GE(*std::min_element(field.begin(), field.end()), 0) << name;
    EXPECT_NEAR(total(l, field), mass, mass * 1e-12) << name;
    const Vec3 expected = rotate(updates * rotation, start);
    const double error_cells = std::acos(std::clamp(mesh::dot(masscentre(l, field), expected), -1.0, 1.0))
        / mesh::cubestep(l.subfacecells());
    const int endface = mesh::cubeposition(end).face;
    double onface = 0;
    for (st::CellId c = 0; c < l.subcellcount(); ++c)
        if (static_cast<int>(c / (l.subfacecells() * l.subfacecells())) == endface) onface += field[c] * l.area_m2(c);
    ::testing::Test::RecordProperty(std::string(name) + "_centre_error_cells", std::to_string(error_cells));
    std::printf("%s: substeps %u, mass centre %.4f subcells from the rotated start, %.4f of the mass on face %d\n",
                name, stencil.substeps(), error_cells, onface / mass, endface);
    EXPECT_LT(error_cells, 0.5) << name;
    EXPECT_GT(onface / mass, 0.6) << name;
}
} // namespace

TEST(TectonicSphereTransport, LatticeNestsInItsParentsWithExactAreas) {
    const auto lattice = st::Lattice::make(8, 2, radius_m);
    const st::Lattice& l = *lattice;
    ASSERT_EQ(l.subcellcount(), 6U * 16 * 16);
    const auto parents = mesh::makecubegrid(8, radius_m);
    std::vector<double> children(l.cellcount(), 0);
    double sum = 0;
    for (st::CellId c = 0; c < l.subcellcount(); ++c) {
        int f = 0, i = 0, j = 0, pf = 0, pi_ = 0, pj = 0;
        mesh::cubeindices(16, c, f, i, j);
        mesh::cubeindices(8, l.parent(c), pf, pi_, pj);
        ASSERT_EQ(pf, f);
        ASSERT_EQ(pi_, i / 2);
        ASSERT_EQ(pj, j / 2);
        children[l.parent(c)] += l.area_m2(c);
        sum += l.area_m2(c);
        for (const auto p : l.profile(c)) ASSERT_LT(p, l.subcellcount());
    }
    EXPECT_NEAR(sum, 4 * pi * radius_m * radius_m, 4 * pi * radius_m * radius_m * 1e-12);
    for (st::CellId p = 0; p < l.cellcount(); ++p)
        EXPECT_NEAR(children[p], parents->area_m2(p), parents->area_m2(p) * 1e-12) << p;
    int gaps = 0;
    for (st::CellId c = 0; c < l.subcellcount(); ++c)
        for (int k = 1; k <= 8; ++k) gaps += l.neighbour(c, k) == st::nocell;
    EXPECT_EQ(gaps, 24);  // one diagonal per corner subcell
}

TEST(TectonicSphereTransport, RotationAboutAnArbitraryAxisConservesMassAndMoment) {
    const auto lattice = st::Lattice::make(8, 2, radius_m);
    const st::Lattice& l = *lattice;
    // About two subcells per update at the rotation's equator: several substeps.
    const Vec3 rotation = (2.2 * mesh::cubestep(16)) * mesh::normalise({0.3, -0.5, 0.8});
    const st::Stencil stencil(lattice, rigid(l, rotation));
    EXPECT_GE(stencil.substeps(), 3U);
    Random random;
    std::vector<double> rock(l.subcellcount()), moment(l.subcellcount());
    for (std::size_t c = 0; c < rock.size(); ++c) {
        rock[c] = c % 7 == 0 ? 0 : 100 + 27000 * random.uniform();
        moment[c] = rock[c] * (random.uniform() < 0.5 ? 100 : 300);
    }
    const double mass = total(l, rock), age = total(l, moment);
    for (int t = 0; t < 40; ++t) {
        rock = stencil.remap(rock, {&moment});
        for (std::size_t c = 0; c < rock.size(); ++c) {
            ASSERT_GE(rock[c], 0) << t << ", subcell " << c;
            if (rock[c] > 0) {
                ASSERT_GE(moment[c] / rock[c], 100 - 1e-9) << t << ", subcell " << c;
                ASSERT_LE(moment[c] / rock[c], 300 + 1e-9) << t << ", subcell " << c;
            }
        }
    }
    EXPECT_NEAR(total(l, rock), mass, mass * 1e-12);
    EXPECT_NEAR(total(l, moment), age, age * 1e-12);
}

// The cube's replacement of TectonicTransport.DiagonalVariableMotionPreservesPositiveMassAcrossBothSeams:
// four plates alternating by parent cell, tangent offsets, zeros in the field.
TEST(TectonicSphereTransport, VariableMotionPreservesPositiveMassAcrossFaceEdges) {
    const auto lattice = st::Lattice::make(8, 2, radius_m);
    const st::Lattice& l = *lattice;
    Random random;
    st::Motion motion;
    const double step = mesh::cubestep(16);
    for (int p = 0; p < 4; ++p) motion.rotation.push_back(((0.4 + 0.5 * p) * step) * random.direction());
    motion.plate.resize(l.cellcount());
    motion.offset_m.resize(l.cellcount());
    for (st::CellId c = 0; c < l.cellcount(); ++c) {
        motion.plate[c] = c % 4;
        motion.offset_m[c] = (0.9 * step * radius_m * random.uniform()) * random.direction();
    }
    motion.plate[5] = st::noplate;
    const st::Stencil stencil(lattice, motion);
    std::vector<double> field(l.subcellcount());
    for (std::size_t c = 0; c < field.size(); ++c) field[c] = c % 5 == 0 ? 0 : static_cast<double>((c * 17) % 31);
    const double mass = total(l, field);
    for (int t = 0; t < 100; ++t) field = stencil.remap(field);
    EXPECT_GE(*std::min_element(field.begin(), field.end()), 0);
    EXPECT_NEAR(total(l, field), mass, mass * 1e-12);
    std::vector<double> negative(l.subcellcount(), 1);
    negative[3] = -1;
    EXPECT_THROW(stencil.remap(negative), std::invalid_argument);
    EXPECT_THROW(stencil.remap({}), std::invalid_argument);
    motion.plate[0] = 4;
    EXPECT_THROW(st::Stencil(lattice, motion), std::invalid_argument);
}

TEST(TectonicSphereTransport, StationaryCellsKeepTheirValuesExactly) {
    const auto lattice = st::Lattice::make(8, 2, radius_m);
    const st::Lattice& l = *lattice;
    st::Motion motion = rigid(l, {0, 0, 0});
    const st::Stencil stencil(lattice, motion);
    EXPECT_EQ(stencil.substeps(), 1U);
    EXPECT_EQ(stencil.maximumdisplacement(), 0);
    std::vector<double> field(l.subcellcount());
    for (std::size_t c = 0; c < field.size(); ++c) field[c] = static_cast<double>((c * 13) % 29);
    const auto moved = stencil.remap(field);
    for (std::size_t c = 0; c < field.size(); ++c) EXPECT_NEAR(moved[c], field[c], 1e-12 * (1 + field[c])) << c;
}

TEST(TectonicSphereTransport, CapCrossesAFaceEdge) {
    crossing(mesh::unitvector(5 * pi / 180, 33 * pi / 180), mesh::unitvector(-5 * pi / 180, 57 * pi / 180), "edge");
}

TEST(TectonicSphereTransport, CapCrossesACubeCorner) {
    const Vec3 corner = mesh::normalise({1, 1, 1});
    const Vec3 e = mesh::normalise({1.4, -0.6, -0.8});
    crossing(mesh::normalise(corner - 0.3 * e), mesh::normalise(corner + 0.3 * e), "corner");
}

TEST(TectonicSphereTransport, CapCrossesAPole) {
    crossing(mesh::unitvector(72 * pi / 180, 10 * pi / 180), mesh::unitvector(72 * pi / 180, -170 * pi / 180), "pole");
}

// The maxima of |value - expected| by subcell class: face interiors, face-edge subcells (corners
// excluded) and the 24 corner subcells.
namespace {
struct Classes {
    double interior = 0, edge = 0, corner = 0;
    double all() const { return std::max({interior, edge, corner}); }
};
Classes classify(const st::Lattice& l, const std::vector<double>& values, double expected) {
    Classes result;
    const int m = l.subfacecells();
    int corners = 0;
    for (st::CellId c = 0; c < l.subcellcount(); ++c) {
        int f = 0, i = 0, j = 0;
        mesh::cubeindices(m, c, f, i, j);
        const bool bi = i == 0 || i == m - 1, bj = j == 0 || j == m - 1;
        const double error = std::abs(values[c] - expected);
        if (bi && bj) { result.corner = std::max(result.corner, error); ++corners; }
        else if (bi || bj) result.edge = std::max(result.edge, error);
        else result.interior = std::max(result.interior, error);
    }
    EXPECT_EQ(corners, 24);
    return result;
}
} // namespace

// One plate over the whole sphere, rotating about two axes for 100 updates: the transported
// uniform field stays uniform and the plate's coverage 1 (interiors, face edges, corners), and
// the plate's amounts of a uniform and a varying field are conserved. The raw kernel's coverage
// after one update is recorded for comparison.
TEST(TectonicSphereTransport, RigidPlateKeepsAUniformFieldAndItsCoverage) {
    for (const auto [n, s] : {std::pair{8, 2}, std::pair{32, 4}}) {
        const auto lattice = st::Lattice::make(n, s, radius_m);
        const st::Lattice& l = *lattice;
        for (const Vec3 axis : {mesh::normalise({0.3, -0.5, 0.8}), mesh::normalise({-0.7, 0.2, 0.1})}) {
            const Vec3 rotation = (0.6 * mesh::cubestep(l.subfacecells())) * axis;
            const st::Stencil stencil(lattice, rigid(l, rotation));
            const Classes raw = classify(l, stencil.remap(std::vector<double>(l.subcellcount(), 1)), 1);
            std::vector<double> uniform(l.subcellcount(), 2500), varying(l.subcellcount());
            Random random;
            for (auto& v : varying) v = random.uniform() < 0.1 ? 0 : 100 + 27000 * random.uniform();
            const double uniformmass = total(l, uniform), varyingmass = total(l, varying);
            for (int t = 1; t <= 100; ++t) {
                const auto layers = stencil.transport({&uniform, &varying});
                ASSERT_EQ(layers.size(), 1U);
                const st::PlateLayer& layer = layers[0];
                ASSERT_EQ(layer.cells.size(), l.subcellcount());
                std::vector<double> coverage(l.subcellcount(), 0);
                double um = 0, vm = 0;
                for (std::size_t k = 0; k < layer.cells.size(); ++k) {
                    const st::CellId c = layer.cells[k];
                    coverage[c] = layer.coverage[k];
                    uniform[c] = layer.values[0][k];
                    varying[c] = layer.values[1][k];
                    um += uniform[c] * coverage[c] * l.area_m2(c);
                    vm += varying[c] * coverage[c] * l.area_m2(c);
                }
                ASSERT_GE(*std::min_element(varying.begin(), varying.end()), 0);
                if (t != 1 && t != 20 && t != 100) continue;
                const Classes value = classify(l, uniform, 2500), cover = classify(l, coverage, 1);
                std::printf("N %d s %d axis (%.2f, %.2f, %.2f) after %3d: max |value / 2500 - 1| interior %.1e, edge "
                            "%.1e, corner %.1e; max |coverage - 1| interior %.1e, edge %.1e, corner %.1e; mass error "
                            "uniform %.1e, varying %.1e (raw kernel, one update: interior %.1e, edge %.1e, corner %.1e)\n",
                            n, s, axis.x, axis.y, axis.z, t, value.interior / 2500, value.edge / 2500,
                            value.corner / 2500, cover.interior, cover.edge, cover.corner,
                            std::abs(um / uniformmass - 1), std::abs(vm / varyingmass - 1), raw.interior, raw.edge,
                            raw.corner);
                const std::string key = "N" + std::to_string(n) + "_s" + std::to_string(s) + "_t" + std::to_string(t);
                RecordProperty(key + "_value_error", std::to_string(value.all() / 2500));
                RecordProperty(key + "_coverage_error", std::to_string(cover.all()));
                EXPECT_LT(value.all() / 2500, 1e-9) << t;
                EXPECT_LT(cover.all(), 1e-9) << t;
                EXPECT_NEAR(um, uniformmass, uniformmass * 1e-9) << t;
                EXPECT_NEAR(vm, varyingmass, varyingmass * 1e-9) << t;
            }
        }
    }
}

// Two plates split by the plane y = 0 turning in opposite senses about Z: they converge along
// longitude 0 and diverge along 180. The summed coverage exceeds 1 only along the converging
// contact, falls short only along the diverging one and is 1 within 1e-9 everywhere else.
TEST(TectonicSphereTransport, ConvergingPlatesOverlapOnlyAlongTheirContact) {
    for (const auto [n, s] : {std::pair{8, 2}, std::pair{32, 4}}) {
        const auto lattice = st::Lattice::make(n, s, radius_m);
        const st::Lattice& l = *lattice;
        const auto parents = mesh::makecubegrid(n, radius_m);
        const double step = mesh::cubestep(l.subfacecells());
        st::Motion motion;
        motion.rotation = {Vec3{0, 0, 0.6 * step}, Vec3{0, 0, -0.6 * step}};
        motion.plate.resize(l.cellcount());
        for (st::CellId c = 0; c < l.cellcount(); ++c) motion.plate[c] = parents->centre(c).y < 0 ? 0U : 1U;
        const st::Stencil stencil(lattice, motion);
        std::vector<double> rock(l.subcellcount(), 30000);
        const auto layers = stencil.transport({&rock});
        ASSERT_EQ(layers.size(), 2U);
        std::vector<double> sum(l.subcellcount(), 0);
        for (const auto& layer : layers)
            for (std::size_t k = 0; k < layer.cells.size(); ++k) {
                sum[layer.cells[k]] += layer.coverage[k];
                EXPECT_NEAR(layer.values[0][k], 30000, 30000 * 1e-9);
            }
        // The angle to the half meridian at longitude 0 (x >= 0) or 180 (x <= 0).
        const auto distance = [](const Vec3& p, bool east) {
            return (east ? p.x >= 0 : p.x <= 0) ? std::asin(std::min(1.0, std::abs(p.y))) : std::acos(std::min(1.0, std::abs(p.z)));
        };
        const double band = 6 * step;
        int overlap = 0, gap = 0, stray = 0;
        double elsewhere = 0;
        for (st::CellId c = 0; c < l.subcellcount(); ++c) {
            const Vec3& p = l.centre(c);
            const bool converging = distance(p, true) <= band, diverging = distance(p, false) <= band;
            if (sum[c] > 1 + 1e-6) { ++overlap; stray += !converging; }
            else if (sum[c] < 1 - 1e-6) { ++gap; stray += !diverging; }
            if (!converging && !diverging) elsewhere = std::max(elsewhere, std::abs(sum[c] - 1));
        }
        std::printf("N %d s %d: %d overlapping and %d gap subcells, %d outside their contact band; max |coverage - 1| "
                    "elsewhere %.1e\n", n, s, overlap, gap, stray, elsewhere);
        EXPECT_GT(overlap, 0);
        EXPECT_GT(gap, 0);
        EXPECT_EQ(stray, 0);
        EXPECT_LT(elsewhere, 1e-9);
        // Plate 0's front lies inside plate 1 by the rotation's angle: its rock beyond y = 0 is
        // the lune of the converging contact, 2 theta R^2 (rock moved by the motion of the plate
        // it enters would stop at the contact and leave none).
        double front = 0;
        for (const auto& layer : layers)
            if (layer.plate == 0)
                for (std::size_t k = 0; k < layer.cells.size(); ++k)
                    if (l.centre(layer.cells[k]).y > 0) front += layer.coverage[k] * l.area_m2(layer.cells[k]);
        const double lune = 2 * 0.6 * step * radius_m * radius_m;
        RecordProperty("front_over_lune_n" + std::to_string(n), std::to_string(front / lune));
        EXPECT_NEAR(front / lune, 1, 0.02);
    }
}

namespace {
// A motion scaled so that its largest index displacement is target lattice cells per update.
st::Motion scaled(const std::shared_ptr<const st::Lattice>& lattice, st::Motion motion, double target) {
    const double largest = st::Stencil(lattice, motion).maximumdisplacement();
    const double factor = target / largest;
    for (auto& w : motion.rotation) w = factor * w;
    for (auto& t : motion.offset_m) t = factor * t;
    return motion;
}

// A layer's coverage or field on the whole lattice (zero where the plate is absent).
std::vector<double> dense(const st::Lattice& l, const st::PlateLayer& layer, int field) {
    std::vector<double> v(l.subcellcount(), 0);
    for (std::size_t k = 0; k < layer.cells.size(); ++k)
        v[layer.cells[k]] = field < 0 ? layer.coverage[k] : layer.values[static_cast<std::size_t>(field)][k];
    return v;
}

// The minimum and maximum of v over each subcell's rings-step 8-neighbourhood, restricted to
// the subcells where use is set (+-infinity where none is).
void neighbourhood(const st::Lattice& l, const std::vector<double>& v, const std::vector<std::uint8_t>& use, int rings,
                   std::vector<double>& lo, std::vector<double>& hi) {
    const double inf = std::numeric_limits<double>::infinity();
    lo.assign(l.subcellcount(), inf);
    hi.assign(l.subcellcount(), -inf);
    for (st::CellId c = 0; c < l.subcellcount(); ++c)
        if (use[c]) lo[c] = hi[c] = v[c];
    for (int r = 0; r < rings; ++r) {
        auto nlo = lo, nhi = hi;
        for (st::CellId c = 0; c < l.subcellcount(); ++c)
            for (int k = 1; k <= 8; ++k) {
                const st::CellId m = l.neighbour(c, k);
                if (m == st::nocell) continue;
                nlo[c] = std::min(nlo[c], lo[m]);
                nhi[c] = std::max(nhi[c], hi[m]);
            }
        lo.swap(nlo);
        hi.swap(nhi);
    }
}

// A smooth field on the sphere.
double smooth(const Vec3& p) {
    const Vec3 a = mesh::normalise({0.2, 0.9, -0.4}), b = mesh::normalise({-0.6, 0.1, 0.8});
    return 1000 + 600 * mesh::dot(p, a) * mesh::dot(p, a) + 300 * mesh::dot(p, b) + 150 * std::sin(3 * p.z);
}

bool cornercell(int m, int i, int j) { return (i == 0 || i == m - 1) && (j == 0 || j == m - 1); }
bool edgecell(int m, int i, int j) { return !cornercell(m, i, j) && (i == 0 || i == m - 1 || j == 0 || j == m - 1); }
} // namespace

// Two plates split by y = 0 turning oppositely about Z at about 1.6 and 3.2 subcells per update
// (2 and 4 substeps), carrying a varying field: overlap only along the converging contact, gaps
// only along the diverging one, and every value within its plate's donors (the values of the
// plate's cells within substeps + 2 rings) before the plate's total-restoring scale (reported):
// exactly where the plate's coverage is fractional, within a stated margin inside the plate (the
// field's own profile keeps smooth extrema).
TEST(TectonicSphereTransport, ConvergingPlatesKeepVaryingValuesWithinTheirDonors) {
    for (const auto [n, s] : {std::pair{8, 2}, std::pair{32, 4}})
        for (const double target : {1.6, 3.2}) {
            const auto lattice = st::Lattice::make(n, s, radius_m);
            const st::Lattice& l = *lattice;
            const auto parents = mesh::makecubegrid(n, radius_m);
            st::Motion motion;
            motion.rotation = {Vec3{0, 0, 1e-3}, Vec3{0, 0, -1e-3}};
            motion.plate.resize(l.cellcount());
            for (st::CellId c = 0; c < l.cellcount(); ++c) motion.plate[c] = parents->centre(c).y < 0 ? 0U : 1U;
            motion = scaled(lattice, motion, target);
            const st::Stencil stencil(lattice, motion);
            EXPECT_EQ(stencil.substeps(), target < 2 ? 2U : 4U);
            std::vector<double> field(l.subcellcount());
            for (st::CellId c = 0; c < l.subcellcount(); ++c) field[c] = smooth(l.centre(c));
            const double mass = total(l, field);
            const auto layers = stencil.transport({&field});
            ASSERT_EQ(layers.size(), 2U);
            std::vector<double> sum(l.subcellcount(), 0);
            double amount = 0, edgeexcess = 0, innerexcess = 0, scaling = 0;
            for (std::uint32_t p = 0; p < 2; ++p) {
                const auto& layer = layers[p];
                scaling = std::max(scaling, std::abs(layer.scale[0] - 1));
                std::vector<std::uint8_t> own(l.subcellcount(), 0);
                for (st::CellId c = 0; c < l.subcellcount(); ++c) own[c] = motion.plate[l.parent(c)] == p;
                std::vector<double> lo, hi;
                neighbourhood(l, field, own, static_cast<int>(stencil.substeps()) + 2, lo, hi);
                for (std::size_t k = 0; k < layer.cells.size(); ++k) {
                    const st::CellId c = layer.cells[k];
                    const double v = layer.values[0][k], cov = layer.coverage[k];
                    sum[c] += cov;
                    amount += v * cov * l.area_m2(c);
                    // Before the plate's total-restoring scale.
                    const double unscaled = v / layer.scale[0];
                    const double excess = std::max(lo[c] - unscaled, unscaled - hi[c]) / hi[c];
                    if (cov < 1 - 1e-9) edgeexcess = std::max(edgeexcess, excess);
                    else innerexcess = std::max(innerexcess, excess);
                }
            }
            const double step = mesh::cubestep(l.subfacecells());
            const auto distance = [](const Vec3& p, bool east) {
                return (east ? p.x >= 0 : p.x <= 0) ? std::asin(std::min(1.0, std::abs(p.y)))
                                                     : std::acos(std::min(1.0, std::abs(p.z)));
            };
            const double band = (target + 4) * step;
            int overlap = 0, gap = 0, stray = 0;
            double elsewhere = 0;
            for (st::CellId c = 0; c < l.subcellcount(); ++c) {
                const Vec3& p = l.centre(c);
                const bool converging = distance(p, true) <= band, diverging = distance(p, false) <= band;
                if (sum[c] > 1 + 1e-6) { ++overlap; stray += !converging; }
                else if (sum[c] < 1 - 1e-6) { ++gap; stray += !diverging; }
                if (!converging && !diverging) elsewhere = std::max(elsewhere, std::abs(sum[c] - 1));
            }
            std::printf("converging N %d s %d, %.1f subcells (%u substeps): %d overlapping, %d gap subcells, %d stray; "
                        "|coverage - 1| elsewhere %.1e; values beyond their donors before the scale (relative): "
                        "fractional coverage %.1e, inside %.1e; |scale - 1| %.1e; amount error %.1e\n", n, s, target,
                        stencil.substeps(), overlap, gap, stray, elsewhere, edgeexcess, innerexcess, scaling,
                        std::abs(amount / mass - 1));
            EXPECT_GT(overlap, 0);
            EXPECT_GT(gap, 0);
            EXPECT_EQ(stray, 0);
            EXPECT_LT(elsewhere, 1e-9);
            EXPECT_LE(edgeexcess, 1e-12);
            // Inside, each field's own profile keeps smooth extrema: about twice the largest
            // measured excess (N 8: 6.6e-4, N 32: 1.3e-5; 2026-10-19).
            EXPECT_LE(innerexcess, n == 8 ? 1.5e-3 : 3e-5);
            EXPECT_NEAR(amount, mass, mass * 1e-9);
        }
}

// One rigid plate over the sphere at about 1.6 and 3.2 subcells per update: the transported
// smooth field against its exact pull-back (the field at the back-rotated centre), by class.
TEST(TectonicSphereTransport, RigidPlateFollowsTheExactPullBack) {
    for (const auto [n, s] : {std::pair{8, 2}, std::pair{32, 4}})
        for (const double target : {1.6, 3.2}) {
            const auto lattice = st::Lattice::make(n, s, radius_m);
            const st::Lattice& l = *lattice;
            const st::Motion motion = scaled(lattice, rigid(l, 1e-3 * mesh::normalise({0.3, -0.5, 0.8})), target);
            const st::Stencil stencil(lattice, motion);
            std::vector<double> field(l.subcellcount());
            for (st::CellId c = 0; c < l.subcellcount(); ++c) field[c] = smooth(l.centre(c));
            const auto layers = stencil.transport({&field});
            const auto moved = dense(l, layers[0], 0);
            const int m = l.subfacecells();
            double interior = 0, edge = 0, corner = 0;
            for (st::CellId c = 0; c < l.subcellcount(); ++c) {
                const double exact = smooth(rotate(-1.0 * motion.rotation[0], l.centre(c)));
                const double error = std::abs(moved[c] - exact) / exact;
                int f = 0, i = 0, j = 0;
                mesh::cubeindices(m, c, f, i, j);
                if (cornercell(m, i, j)) corner = std::max(corner, error);
                else if (edgecell(m, i, j)) edge = std::max(edge, error);
                else interior = std::max(interior, error);
            }
            std::printf("pull-back N %d s %d, %.1f subcells (%u substeps): max relative error interior %.2e, edge "
                        "%.2e, corner %.2e\n", n, s, target, stencil.substeps(), interior, edge, corner);
            // Bounds: about twice the largest measured value of 2026-10-19 (E2D rework) over both
            // speeds (N 8: 1.9e-2, 2.6e-2, 2.0e-2; N 32: 2.8e-3, 2.7e-3, 3.2e-3).
            const bool coarse = n == 8;
            EXPECT_LT(interior, coarse ? 4e-2 : 6e-3);
            EXPECT_LT(edge, coarse ? 6e-2 : 6e-3);
            EXPECT_LT(corner, coarse ? 5e-2 : 7e-3);
        }
}

// One plate with tangent offsets converging on the equator within 20 degrees of it, at about 1.6
// and 3.2 subcells per update: coverage above 1 only in the shortening band, below 1 only in it
// or at its edges (the outermost moving rows leave the still ones; moving along meridians toward
// the equator also stretches the parallels, about 0.5 % per update at 12 degrees, seen along the
// parent columns' edges), 1 outside it, and the amounts conserved.
TEST(TectonicSphereTransport, ConvergingOffsetsShortenOnlyTheirBand) {
    for (const auto [n, s] : {std::pair{8, 2}, std::pair{32, 4}})
        for (const double target : {1.6, 3.2}) {
            const auto lattice = st::Lattice::make(n, s, radius_m);
            const st::Lattice& l = *lattice;
            const auto parents = mesh::makecubegrid(n, radius_m);
            const double edge = 20 * pi / 180;
            st::Motion motion = rigid(l, {0, 0, 0});
            motion.offset_m.resize(l.cellcount());
            for (st::CellId c = 0; c < l.cellcount(); ++c) {
                const double lat = parents->latitude_rad(c);
                motion.offset_m[c] = std::abs(lat) < edge ? (-1000 * std::sin(lat) / std::sin(edge)) * parents->north(c) : Vec3{};
            }
            motion = scaled(lattice, motion, target);
            const st::Stencil stencil(lattice, motion);
            EXPECT_EQ(stencil.substeps(), target < 2 ? 2U : 4U);
            std::vector<double> field(l.subcellcount());
            for (st::CellId c = 0; c < l.subcellcount(); ++c) field[c] = smooth(l.centre(c));
            const double mass = total(l, field);
            const auto layers = stencil.transport({&field});
            ASSERT_EQ(layers.size(), 1U);
            const auto& layer = layers[0];
            // The band's edge follows the parent cells: half a parent beyond 20 degrees at most.
            const double margin = (target + s + 4) * mesh::cubestep(l.subfacecells());
            int shortening = 0, extension = 0, stray = 0;
            double amount = 0, elsewhere = 0;
            for (std::size_t k = 0; k < layer.cells.size(); ++k) {
                const st::CellId c = layer.cells[k];
                const double lat = std::abs(mesh::latitudeof(l.centre(c))), cov = layer.coverage[k];
                amount += layer.values[0][k] * cov * l.area_m2(c);
                if (cov > 1 + 1e-6) { ++shortening; stray += lat > edge + margin; }
                else if (cov < 1 - 1e-6) { ++extension; stray += lat > edge + margin; }
                if (lat > edge + margin) elsewhere = std::max(elsewhere, std::abs(cov - 1));
            }
            EXPECT_EQ(layer.cells.size(), l.subcellcount());
            std::printf("offsets N %d s %d, %.1f subcells (%u substeps): %d shortening, %d extension subcells, %d "
                        "stray; |coverage - 1| outside the band %.1e; amount error %.1e\n", n, s, target,
                        stencil.substeps(), shortening, extension, stray, elsewhere, std::abs(amount / mass - 1));
            EXPECT_GT(shortening, 0);
            EXPECT_EQ(stray, 0);
            EXPECT_LT(elsewhere, 1e-9);
            EXPECT_NEAR(amount, mass, mass * 1e-9);
        }
}

// A rotating plate shortened by converging offsets for 100 updates, carrying rock and an age
// moment (rock x age, ages 100 and 300): the moment total's drift per update and the ages'
// range. The next update starts from value x coverage on every subcell.
TEST(TectonicSphereTransport, AgeMomentDriftOver100Updates) {
    const auto lattice = st::Lattice::make(8, 2, radius_m);
    const st::Lattice& l = *lattice;
    const auto parents = mesh::makecubegrid(8, radius_m);
    const double edge = 20 * pi / 180;
    st::Motion motion = rigid(l, 1e-3 * mesh::normalise({0.3, -0.5, 0.8}));
    motion.offset_m.resize(l.cellcount());
    for (st::CellId c = 0; c < l.cellcount(); ++c) {
        const double lat = parents->latitude_rad(c);
        motion.offset_m[c] = std::abs(lat) < edge ? (-2000 * std::sin(lat) / std::sin(edge)) * parents->north(c) : Vec3{};
    }
    motion = scaled(lattice, motion, 1.6);
    const st::Stencil stencil(lattice, motion);
    Random random;
    std::vector<double> rock(l.subcellcount()), moment(l.subcellcount());
    for (std::size_t c = 0; c < rock.size(); ++c) {
        rock[c] = 100 + 27000 * random.uniform();
        moment[c] = rock[c] * (random.uniform() < 0.5 ? 100 : 300);
    }
    double rockdrift = 0, momentdrift = 0, lowest = 1e30, highest = 0;
    for (int t = 0; t < 100; ++t) {
        const double rockmass = total(l, rock), momentmass = total(l, moment);
        const auto layers = stencil.transport({&rock, &moment}, {-1, 0});
        ASSERT_EQ(layers.size(), 1U);
        const auto& layer = layers[0];
        std::fill(rock.begin(), rock.end(), 0.0);
        std::fill(moment.begin(), moment.end(), 0.0);
        for (std::size_t k = 0; k < layer.cells.size(); ++k) {
            const st::CellId c = layer.cells[k];
            rock[c] = layer.values[0][k] * layer.coverage[k];
            moment[c] = layer.values[1][k] * layer.coverage[k];
            if (rock[c] > 0) {
                lowest = std::min(lowest, moment[c] / rock[c]);
                highest = std::max(highest, moment[c] / rock[c]);
            }
        }
        rockdrift = std::max(rockdrift, std::abs(total(l, rock) / rockmass - 1));
        momentdrift = std::max(momentdrift, std::abs(total(l, moment) / momentmass - 1));
    }
    std::printf("100 updates: largest drift per update of the rock total %.1e, of the age moment total %.1e; ages %.12g "
                "to %.12g\n", rockdrift, momentdrift, lowest, highest);
    RecordProperty("moment_drift_per_update", std::to_string(momentdrift));
    EXPECT_LT(rockdrift, 1e-9);
    EXPECT_LT(momentdrift, 1e-9);
    // Each moment restores its own total, so an age moves by the ratio of its moment's scale to
    // its rock's (measured 2.6e-4 of the range over 100 updates at N 8): bounds 1e-3.
    EXPECT_GE(lowest, 100 * (1 - 1e-3));
    EXPECT_LE(highest, 300 * (1 + 1e-3));
}

// Run once: tectonics_tests --gtest_also_run_disabled_tests --gtest_filter=*Timing*
TEST(TectonicSphereTransport, DISABLED_Timing) {
    using clock = std::chrono::steady_clock;
    const auto seconds = [](clock::duration d) { return std::chrono::duration<double>(d).count(); };
    const auto lattice = st::Lattice::make(64, 4, radius_m);
    const st::Lattice& l = *lattice;
    Random random;
    // Twelve plates by the nearest of twelve random poles, about 1.2 subcells per update.
    std::vector<Vec3> seeds(12);
    st::Motion motion;
    for (auto& p : seeds) {
        p = random.direction();
        motion.rotation.push_back((1.2 * mesh::cubestep(256)) * random.direction());
    }
    const auto parents = mesh::makecubegrid(64, radius_m);
    motion.plate.resize(l.cellcount());
    for (st::CellId c = 0; c < l.cellcount(); ++c) {
        std::uint32_t best = 0;
        for (std::uint32_t p = 1; p < seeds.size(); ++p)
            if (mesh::dot(parents->centre(c), seeds[p]) > mesh::dot(parents->centre(c), seeds[best])) best = p;
        motion.plate[c] = best;
    }
    std::vector<double> field(l.subcellcount());
    for (auto& v : field) v = 1000 + 30000 * random.uniform();
    std::vector<double> second(field), third(field), fourth(field);
    double build = 1e30, sphere = 1e30, one = 1e30, four = 1e30;
    std::uint32_t substeps = 0;
    std::size_t stencilbytes = 0;
    for (int repeat = 0; repeat < 3; ++repeat) {
        auto t0 = clock::now();
        const st::Stencil stencil(lattice, motion);
        auto t1 = clock::now();
        const auto moved = stencil.remap(field);
        auto t2 = clock::now();
        const auto layers = stencil.transport({&field});
        auto t3 = clock::now();
        const auto layers4 = stencil.transport({&field, &second, &third, &fourth});
        auto t4 = clock::now();
        substeps = stencil.substeps();
        stencilbytes = stencil.bytes();
        build = std::min(build, seconds(t1 - t0));
        sphere = std::min(sphere, seconds(t2 - t1) / substeps);
        one = std::min(one, seconds(t3 - t2) / substeps);
        four = std::min(four, seconds(t4 - t3) / substeps / 4);
        ASSERT_EQ(moved.size(), field.size());
        ASSERT_EQ(layers.size(), 12U);
        ASSERT_EQ(layers4.size(), 12U);
    }
    const double nsphere = 1e9 * sphere / l.subcellcount(), nbuild = 1e9 * build / l.subcellcount();
    const double none = 1e9 * one / l.subcellcount(), nfour = 1e9 * four / l.subcellcount();
    std::printf("cube N 64 s 4 (%zu subcells, 12 plates, %u substeps): stencil with the rigid references %.1f ns per "
                "subcell; raw remap %.1f ns per subcell and substep; per-plate transport, one field (and the coverage) "
                "%.1f, four fields %.1f per field, ns per subcell and substep\n", l.subcellcount(), substeps, nbuild, nsphere,
                none, nfour);
    std::printf("allocated: lattice %.1f MB, stencil %.1f MB (transport's working arrays: 2 x (fields + 1) x %.1f MB "
                "and a %.1f MB mark, freed on return)\n", l.bytes() / 1e6, stencilbytes / 1e6,
                l.subcellcount() * sizeof(double) / 1e6, l.subcellcount() / 1e6);
    RecordProperty("lattice_bytes", std::to_string(l.bytes()));
    RecordProperty("stencil_bytes", std::to_string(stencilbytes));
    RecordProperty("sphere_stencil_ns", std::to_string(nbuild));
    RecordProperty("sphere_remap_ns", std::to_string(nsphere));
    RecordProperty("sphere_transport_one_field_ns", std::to_string(none));
    RecordProperty("sphere_transport_per_field_of_four_ns", std::to_string(nfour));
}
