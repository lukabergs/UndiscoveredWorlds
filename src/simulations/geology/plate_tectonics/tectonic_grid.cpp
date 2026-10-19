#include "simulations/geology/plate_tectonics/tectonic_grid.hpp"

#include "core/mesh/cube_geometry.hpp"
#include "core/mesh/mesh_geometry.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <tuple>

namespace tectonics::grid {
namespace {

constexpr double kPi = 3.14159265358979323846;
// Mesh slots of west, east, north, south (D8 codes 7 W, 3 E, 1 N, 5 S).
constexpr int kSlot[directioncount] = {7, 3, 1, 5};

std::uint8_t direction_of_slot(int slot) {
    switch (slot) {
    case 7: return west;
    case 3: return east;
    case 1: return north;
    case 5: return south;
    default: throw std::logic_error("TectonicGrid: a cardinal step reached a diagonal slot");
    }
}

// The unit tangent of v at the unit vector r.
Vec3 tangent(const Vec3& v, const Vec3& r) { return mesh::normalise(v - mesh::dot(v, r) * r); }

} // namespace

std::shared_ptr<const TectonicGrid> TectonicGrid::cube(std::shared_ptr<const mesh::CellGrid> mesh) {
    if (!mesh || mesh->kind() != mesh::GridKind::cube || mesh->issubgrid())
        throw std::invalid_argument("TectonicGrid::cube needs a full cube mesh");
    std::shared_ptr<TectonicGrid> grid(new TectonicGrid());
    grid->facecells_ = mesh->facecells();
    grid->cellcount_ = mesh->cellcount();
    grid->totalarea_m2_ = mesh->totalarea_m2();
    grid->cellspacing_m_ = std::sqrt(grid->totalarea_m2_ / static_cast<double>(grid->cellcount_));
    grid->neighbours_.resize(4 * grid->cellcount_);
    grid->continuing_.resize(4 * grid->cellcount_);
    grid->area_.resize(grid->cellcount_);
    grid->smallestedge_.resize(grid->cellcount_);
    for (CellId c = 0; c < grid->cellcount_; ++c) {
        grid->area_[c] = mesh->area_m2(c);
        double smallest = std::numeric_limits<double>::infinity();
        for (std::uint8_t k = 0; k < directioncount; ++k) {
            const std::int32_t n = mesh->neighbour(c, kSlot[k]);
            if (n < 0) throw std::logic_error("TectonicGrid::cube: a cell without a cardinal neighbour");
            const int back = mesh->reverse(c, kSlot[k]);
            // Straight on at n is the opposite of the slot by which n reaches c.
            grid->neighbours_[4 * static_cast<std::size_t>(c) + k] = static_cast<CellId>(n);
            grid->continuing_[4 * static_cast<std::size_t>(c) + k] = direction_of_slot((back + 3) % 8 + 1);
            smallest = std::min(smallest, mesh->facelength_m(c, kSlot[k]));
        }
        grid->smallestedge_[c] = smallest;
    }
    grid->mesh_ = std::move(mesh);
    grid->buildedges();
    return grid;
}

std::shared_ptr<const TectonicGrid> TectonicGrid::cube(int facecells, double radius_km) {
    if (facecells < 1 || !std::isfinite(radius_km) || radius_km <= 0)
        throw std::invalid_argument("TectonicGrid::cube needs N >= 1 and a positive radius");
    return cube(mesh::makecubegrid(facecells, radius_km * 1000.0));
}

void TectonicGrid::buildedges() {
    edges_.clear();
    edges_.reserve(2 * cellcount_);
    edgeof_.assign(4 * cellcount_, UINT32_MAX);
    const int n = facecells_;
    const std::uint8_t order[directioncount] = {east, south, west, north};
    for (CellId c = 0; c < cellcount_; ++c) {
        const Vec3 pa = mesh_->centre(c);
        int face = 0, column = 0, row = 0;
        mesh::cubeindices(n, c, face, column, row);
        for (const std::uint8_t k : order) {
            const CellId other = neighbour(c, k);
            if (other <= c) continue;
            const Vec3 pb = mesh_->centre(other);
            // The shared edge's vertices in c's face coordinates.
            Vec3 p0, p1;
            if (k == east || k == west) {
                const double a = mesh::cubeedgea(n, k == east ? column + 1 : column);
                p0 = mesh::cubepoint(face, a, mesh::cubeedgeb(n, row));
                p1 = mesh::cubepoint(face, a, mesh::cubeedgeb(n, row + 1));
            } else {
                const double b = mesh::cubeedgeb(n, k == south ? row + 1 : row);
                p0 = mesh::cubepoint(face, mesh::cubeedgea(n, column), b);
                p1 = mesh::cubepoint(face, mesh::cubeedgea(n, column + 1), b);
            }
            Edge edge;
            edge.a = c;
            edge.b = other;
            edge.direction = k;
            edge.axis = (k == east || k == west) ? 0 : 1;
            edge.normal = mesh::normalise(mesh::cross(p0, p1));
            if (mesh::dot(edge.normal, pb - pa) < 0) edge.normal = -1.0 * edge.normal;
            edge.centreline = tangent(pb - pa, mesh::normalise(pa + pb));
            edge.midpoint = mesh::normalise(p0 + p1);
            edge.length_m = mesh_->facelength_m(c, kSlot[k]);
            edge.distance_m = mesh_->distance_m(c, kSlot[k]);
            const auto id = static_cast<std::uint32_t>(edges_.size());
            edgeof_[4 * static_cast<std::size_t>(c) + k] = id;
            for (std::uint8_t back = 0; back < directioncount; ++back)
                if (neighbour(other, back) == c) edgeof_[4 * static_cast<std::size_t>(other) + back] = id;
            edges_.push_back(edge);
        }
    }
    // Each adjacent pair shares one edge: 2 per cell on a closed quadrilateral mesh.
    if (edges_.size() != 2 * cellcount_)
        throw std::logic_error("TectonicGrid::cube: the edge table is not 2 edges per cell");
    for (const std::uint32_t e : edgeof_)
        if (e == UINT32_MAX) throw std::logic_error("TectonicGrid::cube: a cell face without an edge");
}

std::vector<std::vector<CellId>> TectonicGrid::rings(CellId c, int radius) const {
    std::vector<std::vector<CellId>> result;
    if (radius < 0 || c >= cellcount_) return result;
    std::vector<std::uint8_t> seen(cellcount_, 0);
    result.push_back({c});
    seen[c] = 1;
    for (int d = 1; d <= radius; ++d) {
        std::vector<CellId> next;
        for (const CellId cell : result.back())
            for (const CellId n : neighbours(cell))
                if (!seen[n]) { seen[n] = 1; next.push_back(n); }
        if (next.empty()) break;
        std::sort(next.begin(), next.end());
        result.push_back(std::move(next));
    }
    return result;
}

std::vector<CellId> TectonicGrid::withinkm(CellId c, double radius_km) const {
    if (c >= cellcount_ || !(radius_km >= 0)) return {};
    return mesh::cellswithin(*mesh_, mesh_->centre(c), radius_km * 1000.0);
}

Distances TectonicGrid::dijkstra(const std::vector<CellId>& sources, const std::vector<std::uint8_t>& passable,
                                 double limit_m) const {
    Distances result;
    result.distance_m.assign(cellcount_, std::numeric_limits<double>::infinity());
    result.source.assign(cellcount_, UINT32_MAX);
    if (!passable.empty() && passable.size() != cellcount_)
        throw std::invalid_argument("TectonicGrid::dijkstra: the passable mask has the wrong size");
    using Entry = std::tuple<double, std::uint32_t, CellId>;
    std::priority_queue<Entry, std::vector<Entry>, std::greater<Entry>> queue;
    for (std::uint32_t s = 0; s < sources.size(); ++s) {
        const CellId c = sources[s];
        if (c >= cellcount_ || (!passable.empty() && !passable[c])) continue;
        if (result.distance_m[c] == 0 && result.source[c] <= s) continue;
        result.distance_m[c] = 0;
        result.source[c] = s;
        queue.emplace(0.0, s, c);
    }
    std::vector<std::uint8_t> settled(cellcount_, 0);
    while (!queue.empty()) {
        const auto [distance, source, c] = queue.top();
        queue.pop();
        if (settled[c] || distance != result.distance_m[c] || source != result.source[c]) continue;
        settled[c] = 1;
        for (std::uint8_t k = 0; k < directioncount; ++k) {
            const CellId n = neighbour(c, k);
            if (settled[n] || (!passable.empty() && !passable[n])) continue;
            const double next = distance + edges_[edge(c, k)].distance_m;
            if (next > limit_m) continue;
            if (next < result.distance_m[n] || (next == result.distance_m[n] && source < result.source[n])) {
                result.distance_m[n] = next;
                result.source[n] = source;
                queue.emplace(next, source, n);
            }
        }
    }
    return result;
}

float TectonicGrid::quantile(const std::vector<float>& values, float q) const {
    if (values.empty()) return 0.0f;
    if (values.size() != cellcount_)
        throw std::invalid_argument("TectonicGrid::quantile needs one value per cell");
    const float clamped = std::max(0.0f, std::min(1.0f, q));
    std::vector<CellId> order(values.size());
    std::iota(order.begin(), order.end(), 0U);
    std::stable_sort(order.begin(), order.end(), [&](CellId l, CellId r) { return values[l] < values[r]; });
    double total = 0;
    for (const CellId c : order) total += area_[c];
    const double target = static_cast<double>(clamped) * (total - area_[order.back()]);
    double below = 0;
    CellId chosen = order.front();
    for (const CellId c : order) {
        if (below > target) break;
        chosen = c;
        below += area_[c];
    }
    return values[chosen];
}

std::vector<Vec3> TectonicGrid::centroids(const std::vector<std::uint32_t>& owners, std::uint32_t count) const {
    if (owners.size() != cellcount_) throw std::invalid_argument("TectonicGrid::centroids needs one owner per cell");
    std::vector<Vec3> sum(count);
    std::vector<double> weight(count, 0.0);
    for (CellId c = 0; c < cellcount_; ++c) {
        const std::uint32_t owner = owners[c];
        if (owner >= count) continue;
        sum[owner] = sum[owner] + area_[c] * mesh_->centre(c);
        weight[owner] += area_[c];
    }
    std::vector<Vec3> result(count);
    std::vector<std::uint32_t> vanishing;
    for (std::uint32_t p = 0; p < count; ++p) {
        if (!(weight[p] > 0)) continue;
        if (mesh::length(sum[p]) >= 1e-6 * weight[p]) result[p] = mesh::normalise(sum[p]);
        else vanishing.push_back(p);
    }
    if (!vanishing.empty()) {
        // Rare: collect the members of the plates whose sums vanish.
        std::vector<std::vector<Vec3>> members(count);
        for (CellId c = 0; c < cellcount_; ++c)
            if (owners[c] < count && std::binary_search(vanishing.begin(), vanishing.end(), owners[c]))
                members[owners[c]].push_back(mesh_->centre(c));
        for (const std::uint32_t p : vanishing)
            result[p] = boundary_cube::meandirection(sum[p], weight[p], members[p]);
    }
    return result;
}

std::array<double, 2> TectonicGrid::tangentoffset_m(const Vec3& origin, CellId c) const {
    const Vec3 o = mesh::normalise(origin);
    const Vec3 p = mesh_->centre(c);
    const double along = mesh::dot(p, o);
    if (!(along > 0)) return {std::nan(""), std::nan("")};
    const mesh::TangentFrame frame = mesh::tangentframe(o, mesh_->radius_m());
    const Vec3 q = (1.0 / along) * p - o;
    return {mesh::dot(q, frame.east) * mesh_->radius_m(), mesh::dot(q, frame.north) * mesh_->radius_m()};
}

} // namespace tectonics::grid

namespace tectonics::boundary_cube {

namespace {
constexpr double kDegrees = 57.295779513082320876798154814105;
} // namespace

double cellspacing_m(const grid::TectonicGrid& grid) { return grid.cellspacing_m(); }

std::array<float, 2> encodecentroid(const mesh::Vec3& unit) {
    return {static_cast<float>(mesh::longitudeof(unit) * kDegrees),
            static_cast<float>(mesh::latitudeof(unit) * kDegrees)};
}

mesh::Vec3 decodecentroid(float x, float y) {
    return mesh::unitvector(static_cast<double>(y) / kDegrees, static_cast<double>(x) / kDegrees);
}

std::array<double, 3> velocity3(const grid::TectonicGrid& grid, const mesh::Vec3* rotation,
                                uint32_t plate, uint32_t plate_count, const mesh::Vec3& p) {
    if (rotation == nullptr || plate >= plate_count) return {0.0, 0.0, 0.0};
    const double scale = grid.radius_m() / grid.cellspacing_m();
    const mesh::Vec3 v = scale * mesh::cross(rotation[plate], p);
    return {v.x, v.y, v.z};
}

std::array<double, 2> eastnorth(const mesh::Vec3& p, const mesh::Vec3& vector) {
    const double latitude = mesh::latitudeof(p), longitude = mesh::longitudeof(p);
    return {mesh::dot(vector, mesh::eastof(latitude, longitude)),
            mesh::dot(vector, mesh::northof(latitude, longitude))};
}

std::array<float, 2> velocity(const grid::TectonicGrid& grid, const mesh::Vec3* rotation,
                              uint32_t plate, uint32_t plate_count, const mesh::Vec3& p) {
    const auto v = velocity3(grid, rotation, plate, plate_count, p);
    const auto en = eastnorth(p, {v[0], v[1], v[2]});
    return {static_cast<float>(en[0]), static_cast<float>(en[1])};
}

mesh::Vec3 meandirection(const mesh::Vec3& sum, double weight, const std::vector<mesh::Vec3>& members) {
    const double length = mesh::length(sum);
    if (length > 0 && length >= 1e-6 * weight) return (1.0 / length) * sum;
    if (members.empty()) return {};
    std::size_t best = 0;
    double nearest = -std::numeric_limits<double>::infinity();
    for (std::size_t k = 0; k < members.size(); ++k) {
        const double along = mesh::dot(members[k], sum);
        if (along > nearest) { nearest = along; best = k; }
    }
    return members[best];
}

} // namespace tectonics::boundary_cube
