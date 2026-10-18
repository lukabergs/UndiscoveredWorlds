#include "simulations/geology/plate_tectonics/tectonic_grid.hpp"

#include "core/mesh/mesh_geometry.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <tuple>

namespace tectonics::grid {
namespace {

constexpr double kTwoPi = 6.283185307179586476925286766559;
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

} // namespace

std::shared_ptr<const TectonicGrid> TectonicGrid::planar(std::uint32_t width, std::uint32_t height,
                                                         double width_km) {
    if (width == 0 || height == 0 || !std::isfinite(width_km) || width_km <= 0)
        throw std::invalid_argument("TectonicGrid::planar needs a positive size and width");
    std::shared_ptr<TectonicGrid> grid(new TectonicGrid());
    grid->kind_ = Kind::planar;
    grid->width_ = width;
    grid->height_ = height;
    grid->cellcount_ = static_cast<std::size_t>(width) * height;
    // The engine's cell edge: the same expression as the torus (world_width_km * 1000.0 / width).
    grid->cell_m_ = width_km * 1000.0 / width;
    grid->cell_area_m2_ = grid->cell_m_ * grid->cell_m_;
    grid->totalarea_m2_ = grid->cell_area_m2_ * static_cast<double>(grid->cellcount_);
    grid->neighbours_.resize(4 * grid->cellcount_);
    grid->continuing_.resize(4 * grid->cellcount_);
    for (std::uint32_t y = 0; y < height; ++y)
        for (std::uint32_t x = 0; x < width; ++x) {
            const std::size_t base = 4 * (static_cast<std::size_t>(y) * width + x);
            grid->neighbours_[base + west] = y * width + (x + width - 1) % width;
            grid->neighbours_[base + east] = y * width + (x + 1) % width;
            grid->neighbours_[base + north] = ((y + height - 1) % height) * width + x;
            grid->neighbours_[base + south] = ((y + 1) % height) * width + x;
            for (std::uint8_t k = 0; k < directioncount; ++k) grid->continuing_[base + k] = k;
        }
    grid->buildedges();
    return grid;
}

std::shared_ptr<const TectonicGrid> TectonicGrid::cube(std::shared_ptr<const mesh::CellGrid> mesh) {
    if (!mesh || mesh->kind() != mesh::GridKind::cube || mesh->issubgrid())
        throw std::invalid_argument("TectonicGrid::cube needs a full cube mesh");
    std::shared_ptr<TectonicGrid> grid(new TectonicGrid());
    grid->kind_ = Kind::cube;
    grid->facecells_ = mesh->facecells();
    grid->cellcount_ = mesh->cellcount();
    grid->totalarea_m2_ = mesh->totalarea_m2();
    grid->neighbours_.resize(4 * grid->cellcount_);
    grid->continuing_.resize(4 * grid->cellcount_);
    grid->area_.resize(grid->cellcount_);
    for (CellId c = 0; c < grid->cellcount_; ++c) {
        grid->area_[c] = mesh->area_m2(c);
        for (std::uint8_t k = 0; k < directioncount; ++k) {
            const std::int32_t n = mesh->neighbour(c, kSlot[k]);
            if (n < 0) throw std::logic_error("TectonicGrid::cube: a cell without a cardinal neighbour");
            const int back = mesh->reverse(c, kSlot[k]);
            // Straight on at n is the opposite of the slot by which n reaches c.
            grid->neighbours_[4 * static_cast<std::size_t>(c) + k] = static_cast<CellId>(n);
            grid->continuing_[4 * static_cast<std::size_t>(c) + k] = direction_of_slot((back + 3) % 8 + 1);
        }
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
    if (kind_ == Kind::planar) {
        edges_.resize(2 * cellcount_);
        for (CellId c = 0; c < cellcount_; ++c) {
            Edge& across_x = edges_[2 * static_cast<std::size_t>(c)];
            across_x = {c, neighbour(c, east), east, 0, {1.0, 0.0, 0.0}, cell_m_, cell_m_};
            Edge& across_y = edges_[2 * static_cast<std::size_t>(c) + 1];
            across_y = {c, neighbour(c, south), south, 1, {0.0, 1.0, 0.0}, cell_m_, cell_m_};
        }
        return;
    }
    edges_.reserve(2 * cellcount_);
    const std::uint8_t order[directioncount] = {east, south, west, north};
    for (CellId c = 0; c < cellcount_; ++c) {
        const Vec3 pa = mesh_->centre(c);
        for (const std::uint8_t k : order) {
            const CellId n = neighbour(c, k);
            if (n <= c) continue;
            const Vec3 pb = mesh_->centre(n);
            const Vec3 middle = mesh::normalise(pa + pb);
            Vec3 tangent = pb - pa;
            tangent = tangent - mesh::dot(tangent, middle) * middle;
            Edge edge;
            edge.a = c;
            edge.b = n;
            edge.direction = k;
            edge.axis = (k == east || k == west) ? 0 : 1;
            edge.normal = mesh::normalise(tangent);
            edge.length_m = mesh_->facelength_m(c, kSlot[k]);
            edge.distance_m = mesh_->distance_m(c, kSlot[k]);
            edges_.push_back(edge);
        }
    }
    // Each adjacent pair shares one edge: 2 per cell on a closed quadrilateral mesh.
    if (edges_.size() != 2 * cellcount_)
        throw std::logic_error("TectonicGrid::cube: the edge table is not 2 edges per cell");
}

Vec3 TectonicGrid::position(CellId c) const {
    if (kind_ == Kind::cube) return mesh_->centre(c);
    return {static_cast<double>(c % width_), static_cast<double>(c / width_), 0.0};
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
    std::vector<CellId> result;
    if (c >= cellcount_ || !(radius_km >= 0)) return result;
    if (kind_ == Kind::cube) return mesh::cellswithin(*mesh_, mesh_->centre(c), radius_km * 1000.0);
    const double radius_m = radius_km * 1000.0;
    const auto reach = static_cast<std::int64_t>(std::ceil(radius_m / cell_m_));
    const std::int64_t reach_x = std::min<std::int64_t>(reach, width_), reach_y = std::min<std::int64_t>(reach, height_);
    const std::int64_t x0 = c % width_, y0 = c / width_;
    const Vec3 origin = position(c);
    for (std::int64_t dy = -reach_y; dy <= reach_y; ++dy)
        for (std::int64_t dx = -reach_x; dx <= reach_x; ++dx) {
            const auto x = static_cast<std::uint32_t>(((x0 + dx) % width_ + width_) % width_);
            const auto y = static_cast<std::uint32_t>(((y0 + dy) % height_ + height_) % height_);
            const CellId n = y * width_ + x;
            const auto offset = tangentoffset_m(origin, n);
            if (std::hypot(offset[0], offset[1]) <= radius_m) result.push_back(n);
        }
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
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
            const double step = kind_ == Kind::cube ? mesh_->distance_m(c, kSlot[k]) : cell_m_;
            const double next = distance + step;
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
    const float clamped = std::max(0.0f, std::min(1.0f, q));
    if (kind_ == Kind::planar) {
        // Equal areas: today's quantile_threshold, float arithmetic included.
        std::vector<float> sorted = values;
        std::sort(sorted.begin(), sorted.end());
        const size_t index = static_cast<size_t>(std::floor(clamped * static_cast<float>(sorted.size() - 1)));
        return sorted[index];
    }
    if (values.size() != cellcount_)
        throw std::invalid_argument("TectonicGrid::quantile needs one value per cell");
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
    if (kind_ == Kind::cube) {
        for (CellId c = 0; c < cellcount_; ++c) {
            const std::uint32_t owner = owners[c];
            if (owner >= count) continue;
            sum[owner] = sum[owner] + area_[c] * mesh_->centre(c);
            weight[owner] += area_[c];
        }
        for (std::uint32_t p = 0; p < count; ++p) {
            const double length = mesh::length(sum[p]);
            sum[p] = weight[p] > 0 && length > 0 ? (1.0 / length) * sum[p] : Vec3{};
        }
        return sum;
    }
    // Circular means per axis: sum stores (cos x, sin x, cos y) and the sin y sum separately.
    std::vector<double> sin_y(count, 0.0);
    for (CellId c = 0; c < cellcount_; ++c) {
        const std::uint32_t owner = owners[c];
        if (owner >= count) continue;
        const double ax = kTwoPi * static_cast<double>(c % width_) / width_;
        const double ay = kTwoPi * static_cast<double>(c / width_) / height_;
        sum[owner].x += std::cos(ax);
        sum[owner].y += std::sin(ax);
        sum[owner].z += std::cos(ay);
        sin_y[owner] += std::sin(ay);
        weight[owner] += 1;
    }
    std::vector<Vec3> result(count);
    for (std::uint32_t p = 0; p < count; ++p) {
        if (weight[p] <= 0) continue;
        double x = std::atan2(sum[p].y, sum[p].x) / kTwoPi * width_;
        double y = std::atan2(sin_y[p], sum[p].z) / kTwoPi * height_;
        if (x < 0) x += width_;
        if (y < 0) y += height_;
        result[p] = {x >= width_ ? 0.0 : x, y >= height_ ? 0.0 : y, 0.0};
    }
    return result;
}

std::array<double, 2> TectonicGrid::tangentoffset_m(const Vec3& origin, CellId c) const {
    if (kind_ == Kind::planar) {
        const Vec3 p = position(c);
        return {std::remainder(p.x - origin.x, static_cast<double>(width_)) * cell_m_,
                std::remainder(p.y - origin.y, static_cast<double>(height_)) * cell_m_};
    }
    const Vec3 o = mesh::normalise(origin);
    const Vec3 p = mesh_->centre(c);
    const double along = mesh::dot(p, o);
    if (!(along > 0)) return {std::nan(""), std::nan("")};
    const mesh::TangentFrame frame = mesh::tangentframe(o, mesh_->radius_m());
    const Vec3 q = (1.0 / along) * p - o;
    return {mesh::dot(q, frame.east) * mesh_->radius_m(), mesh::dot(q, frame.north) * mesh_->radius_m()};
}

} // namespace tectonics::grid
