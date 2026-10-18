#include "simulations/geology/plate_tectonics/tectonic_sphere_transport.hpp"

#include "core/mesh/cube_geometry.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace tectonics::sphere_transport {
namespace {

bool finite(const Vec3& v) { return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z); }

// The slot of a column and row offset in -1..1 (mesh D8 codes; 0 the cell itself).
constexpr int slots[3][3] = {{8, 1, 2}, {7, 0, 3}, {6, 5, 4}};
int slotof(int columnoffset, int rowoffset) { return slots[rowoffset + 1][columnoffset + 1]; }
constexpr int oppositeslot(int slot) { return (slot + 3) % 8 + 1; }

// exp([w]x) v (Rodrigues).
Vec3 rotate(const Vec3& w, const Vec3& v) {
    const double angle = mesh::length(w);
    if (angle == 0) return v;
    const Vec3 k = (1 / angle) * w;
    const double c = std::cos(angle), s = std::sin(angle);
    return c * v + s * mesh::cross(k, v) + (mesh::dot(k, v) * (1 - c)) * k;
}

// The torus kernel's piecewise-parabolic profile (tectonic_transport.cpp).
struct Profile { double linear, quadratic; };
Profile profile(double ll, double l, double q, double r, double rr) {
    double lo = std::clamp((7 * (l + q) - ll - r) / 12, std::min(l, q), std::max(l, q));
    double hi = std::clamp((7 * (q + r) - l - rr) / 12, std::min(q, r), std::max(q, r));
    const double curvature = l - 2 * q + r;
    const bool smooth_extremum = curvature * (ll - 2 * l + q) > 0 && curvature * (q - 2 * r + rr) > 0;
    if (!smooth_extremum) {
        const double delta = hi - lo, six = 6 * q - 3 * (lo + hi);
        if ((hi - q) * (q - lo) <= 0) lo = hi = q;
        else if (delta * six > delta * delta) lo = 3 * q - 2 * hi;
        else if (delta * six < -delta * delta) hi = 3 * q - 2 * lo;
    }
    return {hi - lo, 3 * (lo + hi - 2 * q)};
}
double minimum(const Profile& p) {
    double result = -.5 * std::abs(p.linear) + p.quadratic / 6;
    if (p.quadratic > 0) {
        const double x = -.5 * p.linear / p.quadratic;
        if (std::abs(x) < .5) result = std::min(result, p.linear * x + p.quadratic * (x * x - 1. / 12));
    }
    return result;
}

} // namespace

std::shared_ptr<const Lattice> Lattice::make(int facecells, int subdivision, double radius_m) {
    if (facecells < 1 || subdivision < 1 || facecells * subdivision < 2 || facecells * subdivision > 4096
        || !std::isfinite(radius_m) || radius_m <= 0)
        throw std::invalid_argument("Invalid sphere transport lattice");
    std::shared_ptr<Lattice> lattice(new Lattice());
    Lattice& l = *lattice;
    const int n = facecells * subdivision;
    l.facecells_ = facecells;
    l.subdivision_ = subdivision;
    l.subfacecells_ = n;
    l.cellcount_ = static_cast<std::size_t>(mesh::cubefacecount) * facecells * facecells;
    l.subcellcount_ = static_cast<std::size_t>(mesh::cubefacecount) * n * n;
    l.radius_m_ = radius_m;
    l.parent_.resize(l.subcellcount_);
    l.centre_.resize(l.subcellcount_);
    l.area_.resize(l.subcellcount_);
    l.neighbour_.resize(l.subcellcount_);
    l.profile_.resize(l.subcellcount_);
    const double r2 = radius_m * radius_m;
    for (int f = 0; f < mesh::cubefacecount; ++f)
        for (int j = 0; j < n; ++j)
            for (int i = 0; i < n; ++i) {
                const CellId c = mesh::cubeid(n, f, i, j);
                l.parent_[c] = mesh::cubeid(facecells, f, i / subdivision, j / subdivision);
                l.centre_[c] = mesh::cubepoint(f, mesh::cubecentrea(n, i), mesh::cubecentreb(n, j));
                l.area_[c] = r2 * mesh::cubeunitarea(mesh::cubeedgea(n, i), mesh::cubeedgea(n, i + 1),
                                                     mesh::cubeedgeb(n, j + 1), mesh::cubeedgeb(n, j));
                for (int k = 1; k <= mesh::slotcount; ++k) {
                    const std::int32_t m = mesh::cubeneighbour(n, c, k);
                    l.neighbour_[c][k - 1] = m < 0 ? nocell : static_cast<CellId>(m);
                }
            }
    // Two steps along each face axis; the second continues straight in the first neighbour's
    // frame (opposite the slot by which it reaches the cell).
    const auto ahead = [&](CellId c, int slot, CellId& first) {
        first = l.neighbour_[c][slot - 1];
        for (int back = 1; back <= mesh::slotcount; back += 2)
            if (l.neighbour_[first][back - 1] == c) return l.neighbour_[first][oppositeslot(back) - 1];
        throw std::logic_error("Sphere transport lattice: no reverse slot");
    };
    for (CellId c = 0; c < l.subcellcount_; ++c) {
        auto& p = l.profile_[c];
        p[0] = ahead(c, 7, p[1]);
        p[3] = ahead(c, 3, p[2]);
        p[4] = ahead(c, 1, p[5]);
        p[7] = ahead(c, 5, p[6]);
    }
    return lattice;
}

CellId Lattice::neighbour(CellId subcell, int slot) const {
    if (slot == 0) return subcell;
    if (slot < 1 || slot > mesh::slotcount) return nocell;
    return neighbour_[subcell][slot - 1];
}

Stencil::Stencil(std::shared_ptr<const Lattice> lattice, const Motion& motion) : lattice_(std::move(lattice)) {
    if (!lattice_) throw std::invalid_argument("Sphere transport: no lattice");
    const Lattice& l = *lattice_;
    if (motion.plate.size() != l.cellcount() || (!motion.offset_m.empty() && motion.offset_m.size() != l.cellcount()))
        throw std::invalid_argument("Sphere transport: motion size");
    for (const auto& w : motion.rotation)
        if (!finite(w)) throw std::invalid_argument("Sphere transport: rotation not finite");
    for (const auto& t : motion.offset_m)
        if (!finite(t)) throw std::invalid_argument("Sphere transport: offset not finite");
    for (const auto p : motion.plate)
        if (p != noplate && p >= motion.rotation.size()) throw std::invalid_argument("Sphere transport: plate index");

    const int n = l.subfacecells();
    const double step = mesh::cubestep(n), radius = l.radius_m();
    displacement_.assign(l.subcellcount(), {0, 0});
    maximum_ = 0;
    for (CellId c = 0; c < l.subcellcount(); ++c) {
        const CellId p = l.parent(c);
        Vec3 r = l.centre(c);
        bool moved = false;
        if (!motion.offset_m.empty()) {
            const Vec3& t = motion.offset_m[p];
            const Vec3 tangent = t - mesh::dot(t, r) * r;
            const double length = mesh::length(tangent);
            if (length > 0) {
                const double angle = length / radius;
                r = mesh::normalise(std::cos(angle) * r + (std::sin(angle) / length) * tangent);
                moved = true;
            }
        }
        const std::uint32_t plate = motion.plate[p];
        if (plate != noplate && mesh::length(motion.rotation[plate]) > 0) {
            r = rotate(motion.rotation[plate], r);
            moved = true;
        }
        if (!moved) continue;
        int f = 0, i = 0, j = 0;
        mesh::cubeindices(n, c, f, i, j);
        const mesh::CubeFace& face = mesh::cubeface(f);
        const double d = mesh::dot(r, face.centre);
        if (!(d > 0)) throw std::invalid_argument("Sphere transport: destination a quarter turn from the face");
        const double a = std::atan(mesh::dot(r, face.right) / d), b = std::atan(mesh::dot(r, face.up) / d);
        const std::array<double, 2> delta{(a - mesh::cubecentrea(n, i)) / step, (mesh::cubecentreb(n, j) - b) / step};
        displacement_[c] = delta;
        maximum_ = std::max({maximum_, std::abs(delta[0]), std::abs(delta[1])});
    }
    if (!std::isfinite(maximum_)) throw std::invalid_argument("Sphere transport: displacement not finite");
    substeps_ = static_cast<std::uint32_t>(std::max(1.0, std::ceil(maximum_ / maximumsubstep)));

    entries_.resize(l.subcellcount());
    for (CellId c = 0; c < l.subcellcount(); ++c) {
        Entry& e = entries_[c];
        const double tx = displacement_[c][0] / substeps_, ty = displacement_[c][1] / substeps_;
        const int ix = static_cast<int>(std::floor(tx)), iy = static_cast<int>(std::floor(ty));
        const double fx = tx - ix, fy = ty - iy;
        e.wx = {1 - fx, fx};
        e.cx = {-.5 * fx, .5 * (1 - fx)};
        e.qx = {fx * (fx - .5) / 3, (1 - fx) * (.5 - fx) / 3};
        e.wy = {1 - fy, fy};
        e.cy = {-.5 * fy, .5 * (1 - fy)};
        e.qy = {fy * (fy - .5) / 3, (1 - fy) * (.5 - fy) / 3};
        for (int k = 0; k < 4; ++k) e.target[k] = l.neighbour(c, slotof(ix + (k & 1), iy + (k >> 1)));
    }
}

std::array<double, 2> Stencil::displacement(CellId subcell) const { return displacement_[subcell]; }

void Stencil::substep(const std::vector<double>& values, const std::vector<std::vector<double>>& moments,
                      std::vector<double>& next, std::vector<std::vector<double>>& nextmoments) const {
    const Lattice& l = *lattice_;
    const std::size_t count = l.subcellcount();
    next.assign(count, 0);
    nextmoments.resize(moments.size());
    for (auto& m : nextmoments) m.assign(count, 0);
    for (CellId c = 0; c < count; ++c) {
        const double q = values[c];
        if (q == 0) continue;  // the limiter flattens both profiles: nothing moves
        const auto& s = l.profile(c);
        auto px = profile(values[s[0]], values[s[1]], q, values[s[2]], values[s[3]]);
        auto py = profile(values[s[4]], values[s[5]], q, values[s[6]], values[s[7]]);
        const double extent = -minimum(px) - minimum(py);
        if (extent > q) {
            const double factor = q / extent;
            px.linear *= factor; px.quadratic *= factor;
            py.linear *= factor; py.quadratic *= factor;
        }
        const Entry& e = entries_[c];
        std::array<double, 4> amount{};
        std::array<double, 4> weight{};
        int missing = -1;
        for (int k = 0; k < 4; ++k) {
            const int dx = k & 1, dy = k >> 1;
            const double average = q + px.linear * e.cx[dx] + px.quadratic * e.qx[dx]
                + py.linear * e.cy[dy] + py.quadratic * e.qy[dy];
            weight[k] = e.wx[dx] * e.wy[dy];
            amount[k] = weight[k] * std::max(0.0, average) * l.area_m2(c);  // roundoff at zero only
            if (e.target[k] == nocell) missing = k;
        }
        if (missing >= 0) {
            // A cube corner: the diagonal's share goes to the other three by their weights.
            // The others weigh at least 1 - 0.95^2: a substep moves at most 0.95 cells per axis.
            double others = 0;
            for (int k = 0; k < 4; ++k)
                if (k != missing) others += weight[k];
            const double gap = amount[missing];
            amount[missing] = 0;
            for (int k = 0; k < 4; ++k)
                if (k != missing) amount[k] += gap * (weight[k] / others);
        }
        for (int k = 0; k < 4; ++k) {
            if (amount[k] == 0) continue;
            const CellId t = e.target[k];
            next[t] += amount[k];
            for (std::size_t m = 0; m < moments.size(); ++m) nextmoments[m][t] += amount[k] * (moments[m][c] / q);
        }
    }
    for (CellId c = 0; c < count; ++c) {
        const double area = l.area_m2(c);
        next[c] /= area;
        for (auto& m : nextmoments) m[c] /= area;
    }
}

std::vector<double> Stencil::remap(const std::vector<double>& values) const { return remap(values, {}); }

std::vector<double> Stencil::remap(const std::vector<double>& values,
                                   const std::vector<std::vector<double>*>& moments) const {
    const std::size_t count = lattice_->subcellcount();
    if (values.size() != count) throw std::invalid_argument("Sphere transport: field size");
    for (const double v : values)
        if (!std::isfinite(v) || v < 0) throw std::invalid_argument("Sphere transport: value not finite or negative");
    std::vector<std::vector<double>> current(moments.size()), moved;
    for (std::size_t m = 0; m < moments.size(); ++m) {
        if (!moments[m] || moments[m]->size() != count) throw std::invalid_argument("Sphere transport: moment size");
        for (std::size_t c = 0; c < count; ++c) {
            const double v = (*moments[m])[c];
            if (!std::isfinite(v) || v < 0 || (values[c] == 0 && v != 0))
                throw std::invalid_argument("Sphere transport: invalid material moment");
        }
        current[m] = *moments[m];
    }
    std::vector<double> field = values, next;
    for (std::uint32_t step = 0; step < substeps_; ++step) {
        substep(field, current, next, moved);
        field.swap(next);
        current.swap(moved);
    }
    for (std::size_t m = 0; m < moments.size(); ++m) *moments[m] = std::move(current[m]);
    return field;
}

} // namespace tectonics::sphere_transport
