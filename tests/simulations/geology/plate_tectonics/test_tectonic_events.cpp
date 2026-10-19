// Plate topology events on the cube (W2-E, E3): rifts that split a plate under sustained
// extension, fusion of locked continental plates, and their independence of where on the cube
// the plate lies (a face edge and corner, the pole).
#include "gtest/gtest.h"
#include "cube_fixture.hpp"
#include "simulations/geology/plate_tectonics/tectonic_topology.hpp"
#include "simulations/geology/plate_tectonics/tectonic_columns.hpp"
#include <algorithm>
#include <cmath>
#include <deque>
#include <set>
#include <string>

namespace {
using namespace tectonics::topology;
using cubefixture::CellId;
using cubefixture::Vec3;
using tectonics::grid::TectonicGrid;
constexpr double pi = cubefixture::pi;

// exp([w]x) v.
Vec3 rotate(const Vec3& w, const Vec3& v) {
    const double angle = mesh::length(w);
    if (angle == 0) return v;
    const Vec3 k = (1 / angle) * w;
    return std::cos(angle) * v + std::sin(angle) * mesh::cross(k, v) + (mesh::dot(k, v) * (1 - std::cos(angle))) * k;
}

// An island plate (0) of 80 x 60 degrees in longitude and latitude around the point that the
// rotation `placement` carries (longitude 0, latitude 0) to, inside plate 1, on an N 16 Earth;
// the mantle beneath it spreads east-west, 2e-10 m/s per 250 km away from its centre.
struct Island {
    std::shared_ptr<const TectonicGrid> grid = cubefixture::cube(16);
    std::vector<uint32_t> owners;
    std::vector<Vec3> mantle;
    explicit Island(const Vec3& placement = {}, double sign = 1, int kind = 0) {
        owners.assign(grid->cellcount(), 1);
        mantle.assign(grid->cellcount(), Vec3{});
        for (CellId c = 0; c < grid->cellcount(); ++c) {
            const Vec3 local = rotate(-1.0 * placement, grid->position(c));
            const double lon = mesh::longitudeof(local), lat = mesh::latitudeof(local);
            if (std::abs(lon) >= 40 * pi / 180 || std::abs(lat) >= 30 * pi / 180) continue;
            owners[c] = 0;
            const double delta = sign * lon * grid->radius_m() / 250000 * 2e-10;
            const Vec3 east = mesh::eastof(lat, lon), north = mesh::northof(lat, lon);
            const Vec3 v = kind == 0 ? delta * east : kind == 1 ? delta * north : 1e-9 * (east + north);
            mantle[c] = rotate(placement, v);
        }
    }
    uint32_t cells() const { return static_cast<uint32_t>(std::count(owners.begin(), owners.end(), 0u)); }
};

void check_connected(const TectonicGrid& grid, const std::vector<uint32_t>& owners) {
    std::vector<bool> visited(owners.size()); std::set<uint32_t> seen;
    for (CellId i = 0; i < owners.size(); ++i) if (!visited[i]) {
        EXPECT_TRUE(seen.insert(owners[i]).second) << "detached owner " << owners[i];
        std::deque<CellId> queue{i}; visited[i] = true;
        while (!queue.empty()) {
            const auto a = queue.front(); queue.pop_front();
            for (auto b : grid.neighbours(a))
                if (!visited[b] && owners[b] == owners[a]) { visited[b] = true; queue.push_back(b); }
        }
    }
}

const std::vector<Vec3> at_rest(size_t plates) { return std::vector<Vec3>(plates); }
}

TEST(TectonicEvents, SustainedExtensionSplitsConnectedDomainsAndRecordsAncestry) {
    Island island;
    auto owners = island.owners;
    const auto initial = owners;
    Model model; model.initialize(island.grid, 2, owners);
    std::vector<Change> changes;
    for (int t = 1; t <= 100 && changes.empty(); ++t) {
        owners = model.advance(at_rest(2));
        changes = model.evolve(owners, island.mantle, at_rest(2), std::vector<double>(owners.size(), 1), 1, t);
        if (t < 20) EXPECT_TRUE(changes.empty());
    }
    ASSERT_EQ(changes.size(), 1U); EXPECT_EQ(changes[0].kind, EventKind::Split);
    ASSERT_EQ(model.diagnostics().events.size(), 1U);
    const auto& event = model.diagnostics().events[0];
    EXPECT_EQ(event.parent_a, 0U); EXPECT_EQ(event.result, 2U); EXPECT_GE(event.time_myr, 20);
    EXPECT_EQ(model.diagnostics().stable_ids, (std::vector<uint32_t>{0, 1, 2}));
    for (size_t i = 0; i < owners.size(); ++i) if (initial[i] == 1) EXPECT_EQ(owners[i], 1U);
    check_connected(*island.grid, owners);
    // Once born, daughters are independent bodies with continuous coverage.
    const double column = pi / 32;
    owners = model.advance({{0, 0, -.1 * column}, {}, {0, 0, .1 * column}});
    check_connected(*island.grid, owners);
    EXPECT_EQ(model.diagnostics().events.size(), 1U);
}

TEST(TectonicEvents, ProgressiveDamagePrecedesConnectedSplitAndCreatesExtension) {
    Island island;
    auto owners = island.owners;
    Model model; model.initialize(island.grid, 2, owners);
    Parameters p; p.progressive_rifting = true; p.yield_n_m = 5e11; p.damage_time_myr = 30;
    std::vector<double> strength(owners.size(), 1);
    for (CellId c = 0; c < owners.size(); ++c)
        if (owners[c] == 0 && std::abs(cubefixture::longitude_deg(*island.grid, c)) < 6) strength[c] = .5;
    bool damaged = false, extending = false; std::vector<Change> changes;
    for (int t = 1; t <= 200 && changes.empty(); ++t) {
        owners = model.advance(at_rest(2));
        changes = model.evolve(owners, island.mantle, at_rest(2), std::vector<double>(owners.size(), 1), 1, t, p, strength);
        const auto damage = model.damageField(owners);
        damaged |= *std::max_element(damage.begin(), damage.end()) > 0;
        for (auto v : model.riftDisplacement(owners, 1, p)) extending |= mesh::length(v) > 0;
    }
    EXPECT_TRUE(damaged); EXPECT_TRUE(extending);
    ASSERT_EQ(changes.size(), 1U); EXPECT_EQ(changes[0].kind, EventKind::Split);
    check_connected(*island.grid, owners);
}

TEST(TectonicEvents, CompressionTransformShearAndUniformFlowDoNotSplit) {
    for (int kind = 0; kind < 3; ++kind) {
        Island island({}, kind == 0 ? -1 : 1, kind);
        auto owners = island.owners;
        Model model; model.initialize(island.grid, 2, owners);
        for (int t = 1; t <= 100; ++t)
            EXPECT_TRUE(model.evolve(owners, island.mantle, at_rest(2), std::vector<double>(owners.size(), 1), 1, t).empty());
        EXPECT_TRUE(model.diagnostics().events.empty()) << kind;
    }
}

namespace {
// Two plates split by y = 0 turning about z: +-speed at the equator (m/s); they close at
// longitude 0 and part at 180 degrees.
std::vector<Vec3> about_z(const TectonicGrid& grid, double speed) {
    return {{0, 0, speed / grid.radius_m()}, {0, 0, -speed / grid.radius_m()}};
}
}

TEST(TectonicEvents, FusionRequiresPersistentContinentalLockingAndPreservesIdentityAfterCompaction) {
    const auto grid = cubefixture::cube(16);
    auto owners = cubefixture::hemispheres(*grid);
    Model model; model.initialize(grid, 2, owners);
    const auto rotation = about_z(*grid, 1e-11);
    const std::vector<Vec3> mantle(owners.size());
    const std::vector<double> continent(owners.size(), 1);
    for (int t = 1; t < 30; ++t) EXPECT_TRUE(model.evolve(owners, mantle, rotation, continent, 1, t).empty());
    const auto change = model.evolve(owners, mantle, rotation, continent, 1, 30);
    ASSERT_EQ(change.size(), 1U); EXPECT_EQ(change[0].kind, EventKind::Fusion);
    check_connected(*grid, owners);
    EXPECT_EQ(*std::max_element(owners.begin(), owners.end()), 0U);
    ASSERT_EQ(model.diagnostics().events.size(), 1U);
    EXPECT_EQ(model.diagnostics().events[0].parent_a, 0U);
    EXPECT_EQ(model.diagnostics().events[0].parent_b, 1U);
    EXPECT_EQ(model.diagnostics().events[0].result, 2U);
    model.removePlate(1);
    EXPECT_EQ(model.diagnostics().stable_ids, (std::vector<uint32_t>{2}));
    EXPECT_EQ(model.advance({{.01, -.01, 0}}), owners);
}

TEST(TectonicEvents, OceanContactsFastSlipAndStationaryNeighborsDoNotFuse) {
    for (int kind = 0; kind < 4; ++kind) {
        const auto grid = cubefixture::cube(16);
        auto owners = cubefixture::hemispheres(*grid);
        Model model; model.initialize(grid, 2, owners);
        const double r = grid->radius_m();
        // Ocean closing slowly; continents closing fast; continents sliding along the boundary
        // (turning about y, its plane's normal); continents at rest.
        const std::vector<Vec3> rotation = kind == 0 ? about_z(*grid, 1e-11)
            : kind == 1 ? about_z(*grid, 1e-9)
            : kind == 2 ? std::vector<Vec3>{{0, 1e-11 / r, 0}, {0, -1e-11 / r, 0}} : at_rest(2);
        for (int t = 1; t <= 100; ++t)
            EXPECT_TRUE(model.evolve(owners, std::vector<Vec3>(owners.size()), rotation,
                std::vector<double>(owners.size(), kind == 0 ? 0 : 1), 1, t).empty()) << kind;
    }
}

TEST(TectonicEvents, ShortLockingEpisodeDoesNotCarryOverAcrossRenewedSlip) {
    const auto grid = cubefixture::cube(16);
    auto owners = cubefixture::hemispheres(*grid);
    Model model; model.initialize(grid, 2, owners);
    for (int t = 1; t <= 50; ++t) {
        const double speed = t == 25 ? 1e-9 : 1e-11;
        EXPECT_TRUE(model.evolve(owners, std::vector<Vec3>(owners.size()), about_z(*grid, speed / 2),
            std::vector<double>(owners.size(), 1), 1, t).empty());
    }
}

TEST(TectonicEvents, TimestepRefinementPreservesSustainedFailureTime) {
    double times[2]{};
    for (int k = 0; k < 2; ++k) {
        Island island;
        auto owners = island.owners;
        Model model; model.initialize(island.grid, 2, owners);
        const double dt = k == 0 ? 1 : .25;
        for (double t = dt; t <= 100; t += dt)
            if (!model.evolve(owners, island.mantle, at_rest(2), std::vector<double>(owners.size(), 1), dt, t).empty()) {
                times[k] = t; break;
            }
        EXPECT_GT(times[k], 0);
    }
    EXPECT_NEAR(times[0], times[1], 1);
}

// The same island across a cube face edge and corner (centred on the corner (1, 1, 1)) and over
// the north pole fails at about the same time and splits off about the same area as at
// (longitude 0, latitude 0), a face's centre.
TEST(TectonicEvents, RiftFailureIsAboutUnchangedWhenThePlateCrossesAFaceCornerOrThePole) {
    const Vec3 x{1, 0, 0}, corner = mesh::normalise({1, 1, 1}), pole{0, 0, 1};
    const auto turning = [&](const Vec3& to) {
        const Vec3 axis = mesh::cross(x, to);
        const double sine = mesh::length(axis);
        return sine == 0 ? Vec3{} : (std::atan2(sine, mesh::dot(x, to)) / sine) * axis;
    };
    const Vec3 placements[] = {{}, turning(corner), turning(pole)};
    const char* names[] = {"face_centre", "corner", "pole"};
    double times[3]{}, shares[3]{};
    for (int k = 0; k < 3; ++k) {
        Island island(placements[k]);
        auto owners = island.owners;
        const double total = island.cells();
        Model model; model.initialize(island.grid, 2, owners);
        for (int t = 1; t <= 100; ++t)
            if (!model.evolve(owners, island.mantle, at_rest(2), std::vector<double>(owners.size(), 1), 1, t).empty()) {
                times[k] = t; shares[k] = model.diagnostics().events[0].cells / total; break;
            }
        RecordProperty(std::string(names[k]) + "_max_tension_n_m", std::to_string(model.diagnostics().max_tension_n_m));
        EXPECT_GT(times[k], 0) << names[k]; check_connected(*island.grid, owners);
        RecordProperty(std::string(names[k]) + "_failure_myr", std::to_string(times[k]));
        RecordProperty(std::string(names[k]) + "_daughter_share", std::to_string(shares[k]));
    }
    for (int k = 1; k < 3; ++k) {
        // The failure time goes as 1 / (tension / yield - 1): measured 40 Myr at the face
        // centre and the pole, 46 at the corner (tension 3.31e12 against 3.52e12 N/m).
        EXPECT_NEAR(times[k], times[0], .2 * times[0]) << names[k];
        EXPECT_NEAR(shares[k], shares[0], .15) << names[k];
    }
}

// Detached fragments of a plate are given to an adjacent plate, a plate that wraps a whole
// circle of the sphere (the band's rows 3 to 7 across all four side faces, with the south
// polar face) stays connected.
TEST(TectonicTopology, DetachedFragmentIsRepairedButAPlateRoundTheSphereIsConnected) {
    const auto grid = cubefixture::cube(8);
    const cubefixture::Band band(8);
    std::vector<uint32_t> owners(grid->cellcount(), 1);
    for (int y = 3; y < band.rows(); ++y) for (int x = 0; x < band.columns(); ++x) owners[band.cell(x, y)] = 0;
    for (CellId c = 5 * 64; c < 6 * 64; ++c) owners[c] = 0;
    // A one-cell fragment of plate 0 on the north polar face.
    const CellId fragment = mesh::cubeid(8, 4, 4, 4);
    owners[fragment] = 0;
    Model model; model.initialize(grid, 2, owners);
    const auto moved = model.advance(at_rest(2));
    EXPECT_EQ(moved[fragment], 1u);
    for (int y = 3; y <= 4; ++y) for (int x = 0; x < band.columns(); ++x) EXPECT_EQ(moved[band.cell(x, y)], 0u);
    EXPECT_EQ(model.diagnostics().reassigned_cells, 1u);
    check_connected(*grid, moved);
}
