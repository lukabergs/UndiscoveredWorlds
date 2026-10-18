// The rock transport on the cube (W2-E, step E2D): the subcell lattice, conservation and
// positivity under plate rotations, crossings of a face edge, a cube corner and a pole, the
// coverage of a rigid plate, and (disabled, run once) the cost against the torus kernel.
#include "gtest/gtest.h"

#include "core/mesh/cell_grid.hpp"
#include "core/mesh/cube_geometry.hpp"
#include "simulations/geology/plate_tectonics/tectonic_sphere_transport.hpp"
#include "simulations/geology/plate_tectonics/tectonic_transport.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
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

// A uniform field on one rigid plate: max |coverage - 1| in face interiors, on face-edge
// subcells (corners excluded) and on the 24 corner subcells, after one update and after 20.
TEST(TectonicSphereTransport, RigidPlateCoverage) {
    for (const auto [n, s] : {std::pair{8, 2}, std::pair{32, 4}}) {
        const auto lattice = st::Lattice::make(n, s, radius_m);
        const st::Lattice& l = *lattice;
        const int m = l.subfacecells();
        const Vec3 rotation = (0.6 * mesh::cubestep(m)) * mesh::normalise({0.3, -0.5, 0.8});
        const st::Stencil stencil(lattice, rigid(l, rotation));
        std::vector<double> field(l.subcellcount(), 1);
        // A subcell's upstream path stays farther than 4 subcells from every cube corner.
        std::vector<std::uint8_t> clear(l.subcellcount(), 1);
        const auto nearcorner = [&](const Vec3& p) {
            const double nearest = (std::abs(p.x) + std::abs(p.y) + std::abs(p.z)) / std::sqrt(3.0);
            return std::acos(std::min(1.0, nearest)) <= 4 * mesh::cubestep(m);
        };
        for (st::CellId c = 0; c < l.subcellcount(); ++c) clear[c] = !nearcorner(l.centre(c));
        for (int t = 1; t <= 100; ++t) {
            field = stencil.remap(field);
            for (st::CellId c = 0; c < l.subcellcount(); ++c)
                if (clear[c] && nearcorner(rotate(-static_cast<double>(t) * rotation, l.centre(c)))) clear[c] = 0;
            if (t != 1 && t != 20 && t != 100) continue;
            double interior = 0, edge = 0, corner = 0, away = 0;
            int corners = 0;
            for (st::CellId c = 0; c < l.subcellcount(); ++c) {
                int f = 0, i = 0, j = 0;
                mesh::cubeindices(m, c, f, i, j);
                const bool bi = i == 0 || i == m - 1, bj = j == 0 || j == m - 1;
                const double error = std::abs(field[c] - 1);
                if (bi && bj) { corner = std::max(corner, error); ++corners; }
                else if (bi || bj) edge = std::max(edge, error);
                else interior = std::max(interior, error);
                if (clear[c]) away = std::max(away, error);
            }
            ASSERT_EQ(corners, 24);
            const std::string key = "N" + std::to_string(n) + "_s" + std::to_string(s) + "_t" + std::to_string(t);
            RecordProperty(key + "_interior", std::to_string(interior));
            RecordProperty(key + "_edge", std::to_string(edge));
            RecordProperty(key + "_corner", std::to_string(corner));
            RecordProperty(key + "_away_from_corners", std::to_string(away));
            std::printf("coverage N %d s %d after %d update(s): max |coverage - 1| interior %.3e, edge %.3e, "
                        "corner %.3e, away from the corners' wakes %.3e\n", n, s, t, interior, edge, corner, away);
            // Regression bounds over the measured values (E2D, 2026-10-19), not quality claims.
            // The corner subcells keep about a third: three face frames meet at 120 degrees there
            // and the splat assumes 90; the defect is carried downstream as a wake (the later
            // interior and edge maxima). Elsewhere the index displacement of a rigid rotation
            // varies in space and the bilinear splat ripples by about its gradient per update
            // (first order: smaller at finer s N), and the ripples grow over updates ("away").
            if (t == 1) {
                EXPECT_LT(interior, n == 8 ? 0.025 : 0.006);
                EXPECT_LT(edge, n == 8 ? 0.05 : 0.008);
            }
            EXPECT_LT(corner, 0.5);
            EXPECT_LT(away, 0.1);
        }
    }
}

// Run once: tectonics_tests --gtest_also_run_disabled_tests --gtest_filter=*Timing*
TEST(TectonicSphereTransport, DISABLED_TimingAgainstTheTorusKernel) {
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
    double build = 1e30, sphere = 1e30;
    std::uint32_t substeps = 0;
    for (int repeat = 0; repeat < 3; ++repeat) {
        auto t0 = clock::now();
        const st::Stencil stencil(lattice, motion);
        auto t1 = clock::now();
        const auto moved = stencil.remap(field);
        auto t2 = clock::now();
        substeps = stencil.substeps();
        build = std::min(build, seconds(t1 - t0));
        sphere = std::min(sphere, seconds(t2 - t1) / substeps);
        ASSERT_EQ(moved.size(), field.size());
    }
    // The torus kernel at 256 x 128 with subdivision 4 (the baseline's): 1024 x 512 subcells.
    const std::uint32_t w = 1024, h = 512;
    std::vector<double> torus(static_cast<std::size_t>(w) * h);
    std::vector<tectonics::transport::Displacement> v(torus.size());
    for (std::size_t i = 0; i < torus.size(); ++i) {
        torus[i] = 1000 + 30000 * random.uniform();
        v[i] = {0.9 * (2 * random.uniform() - 1), 0.9 * (2 * random.uniform() - 1)};
    }
    double kernel = 1e30;
    for (int repeat = 0; repeat < 3; ++repeat) {
        auto t0 = clock::now();
        const auto moved = tectonics::transport::remap(w, h, torus, v);
        kernel = std::min(kernel, seconds(clock::now() - t0));
        ASSERT_EQ(moved.size(), torus.size());
    }
    const double nsphere = 1e9 * sphere / l.subcellcount(), nbuild = 1e9 * build / l.subcellcount();
    const double ntorus = 1e9 * kernel / torus.size();
    std::printf("cube N 64 s 4 (%zu subcells, %u substeps): stencil %.1f ns per subcell, remap %.1f ns per subcell and "
                "substep; torus 1024 x 512: remap %.1f ns per subcell and substep\n",
                l.subcellcount(), substeps, nbuild, nsphere, ntorus);
    RecordProperty("sphere_stencil_ns", std::to_string(nbuild));
    RecordProperty("sphere_remap_ns", std::to_string(nsphere));
    RecordProperty("torus_remap_ns", std::to_string(ntorus));
}
