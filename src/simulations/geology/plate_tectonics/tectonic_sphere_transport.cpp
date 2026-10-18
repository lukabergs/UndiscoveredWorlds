#include "simulations/geology/plate_tectonics/tectonic_sphere_transport.hpp"

#include "core/mesh/cube_geometry.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace tectonics::sphere_transport {
namespace {

using detail::Entry;

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

// The index displacement from a subcell's centre to r, in extended coordinates of its face.
std::array<double, 2> indexdisplacement(const Lattice& l, CellId c, const Vec3& r) {
    const int n = l.subfacecells();
    int f = 0, i = 0, j = 0;
    mesh::cubeindices(n, c, f, i, j);
    const mesh::CubeFace& face = mesh::cubeface(f);
    const double d = mesh::dot(r, face.centre);
    if (!(d > 0)) throw std::invalid_argument("Sphere transport: destination a quarter turn from the face");
    const double a = std::atan(mesh::dot(r, face.right) / d), b = std::atan(mesh::dot(r, face.up) / d);
    const double step = mesh::cubestep(n);
    return {(a - mesh::cubecentrea(n, i)) / step, (mesh::cubecentreb(n, j) - b) / step};
}

// A rotation's index displacement at a subcell ({0, 0} for no rotation).
std::array<double, 2> rotationdisplacement(const Lattice& l, CellId c, const Vec3& w) {
    if (!(mesh::length(w) > 0)) return {0, 0};
    return indexdisplacement(l, c, rotate(w, l.centre(c)));
}

Entry makeentry(const Lattice& l, CellId c, double tx, double ty) {
    Entry e;
    const int ix = static_cast<int>(std::floor(tx)), iy = static_cast<int>(std::floor(ty));
    const double fx = tx - ix, fy = ty - iy;
    e.wx = {1 - fx, fx};
    e.cx = {-.5 * fx, .5 * (1 - fx)};
    e.qx = {fx * (fx - .5) / 3, (1 - fx) * (.5 - fx) / 3};
    e.wy = {1 - fy, fy};
    e.cy = {-.5 * fy, .5 * (1 - fy)};
    e.qy = {fy * (fy - .5) / 3, (1 - fy) * (.5 - fy) / 3};
    for (int k = 0; k < 4; ++k) e.target[k] = l.neighbour(c, slotof(ix + (k & 1), iy + (k >> 1)));
    return e;
}

// The seeds (ascending) and every subcell within rings 8-neighbour steps, ascending. mark is
// all zero on entry and on return.
std::vector<CellId> dilate(const Lattice& l, const std::vector<CellId>& seeds, std::uint32_t rings,
                           std::vector<std::uint8_t>& mark) {
    std::vector<CellId> result = seeds, front = seeds, next;
    for (const CellId c : seeds) mark[c] = 1;
    for (std::uint32_t r = 0; r < rings && !front.empty(); ++r) {
        next.clear();
        for (const CellId c : front)
            for (int k = 1; k <= mesh::slotcount; ++k) {
                const CellId m = l.neighbour(c, k);
                if (m != nocell && !mark[m]) {
                    mark[m] = 1;
                    next.push_back(m);
                }
            }
        result.insert(result.end(), next.begin(), next.end());
        front.swap(next);
    }
    for (const CellId c : result) mark[c] = 0;
    std::sort(result.begin(), result.end());
    return result;
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

// Which fields move on their own and which follow a rock field's transfer fractions.
struct Layout {
    std::vector<int> independent;
    std::vector<std::vector<int>> moments;  // per field
};

// One substep over the sources of region (ascending ids; entryof(k, region[k]) their entries).
// current: dense fields, zero outside region; next: dense, zero on entry and on return. A
// destination outside inside (when given) is dropped, or an error when strict. The result,
// divided by the destinations' areas, replaces current on region.
template <class EntryOf>
void substep(const Lattice& l, const std::vector<CellId>& region, const EntryOf& entryof, const Layout& layout,
             const std::vector<std::uint8_t>* inside, bool strict, std::vector<std::vector<double>>& current,
             std::vector<std::vector<double>>& next) {
    for (std::size_t k = 0; k < region.size(); ++k) {
        const CellId c = region[k];
        const Entry& e = entryof(k, c);
        const auto& s = l.profile(c);
        for (const int f : layout.independent) {
            const std::vector<double>& values = current[f];
            const double q = values[c];
            if (q == 0) continue;  // the limiter flattens both profiles: nothing moves
            auto px = profile(values[s[0]], values[s[1]], q, values[s[2]], values[s[3]]);
            auto py = profile(values[s[4]], values[s[5]], q, values[s[6]], values[s[7]]);
            const double extent = -minimum(px) - minimum(py);
            if (extent > q) {
                const double factor = q / extent;
                px.linear *= factor; px.quadratic *= factor;
                py.linear *= factor; py.quadratic *= factor;
            }
            std::array<double, 4> amount{};
            std::array<double, 4> weight{};
            int missing = -1;
            for (int t = 0; t < 4; ++t) {
                const int dx = t & 1, dy = t >> 1;
                const double average = q + px.linear * e.cx[dx] + px.quadratic * e.qx[dx]
                    + py.linear * e.cy[dy] + py.quadratic * e.qy[dy];
                weight[t] = e.wx[dx] * e.wy[dy];
                amount[t] = weight[t] * std::max(0.0, average) * l.area_m2(c);  // roundoff at zero only
                if (e.target[t] == nocell) missing = t;
            }
            if (missing >= 0) {
                // A cube corner: the diagonal's share goes to the other three by their weights,
                // which weigh at least 1 - 0.95^2 (a substep moves at most 0.95 cells per axis).
                double others = 0;
                for (int t = 0; t < 4; ++t)
                    if (t != missing) others += weight[t];
                const double gap = amount[missing];
                amount[missing] = 0;
                for (int t = 0; t < 4; ++t)
                    if (t != missing) amount[t] += gap * (weight[t] / others);
            }
            for (int t = 0; t < 4; ++t) {
                if (amount[t] == 0) continue;
                const CellId target = e.target[t];
                if (inside && !(*inside)[target]) {
                    if (strict) throw std::logic_error("Sphere transport: rock beyond its plate's reach");
                    continue;
                }
                next[f][target] += amount[t];
                for (const int m : layout.moments[f]) next[m][target] += amount[t] * (current[m][c] / q);
            }
        }
    }
    for (const CellId c : region) {
        const double area = l.area_m2(c);
        for (std::size_t f = 0; f < current.size(); ++f) {
            current[f][c] = next[f][c] / area;
            next[f][c] = 0;
        }
    }
}

// Sizes, finite non-negative values, moments of independent fields and zero where their rock is.
Layout validate(std::size_t count, const std::vector<const std::vector<double>*>& fields, const std::vector<int>& momentof) {
    if (!momentof.empty() && momentof.size() != fields.size()) throw std::invalid_argument("Sphere transport: moment layout");
    Layout layout;
    layout.moments.resize(fields.size());
    for (std::size_t f = 0; f < fields.size(); ++f) {
        if (!fields[f] || fields[f]->size() != count) throw std::invalid_argument("Sphere transport: field size");
        for (const double v : *fields[f])
            if (!std::isfinite(v) || v < 0) throw std::invalid_argument("Sphere transport: value not finite or negative");
        const int rock = momentof.empty() ? -1 : momentof[f];
        if (rock < 0) { layout.independent.push_back(static_cast<int>(f)); continue; }
        if (static_cast<std::size_t>(rock) >= fields.size() || rock == static_cast<int>(f) || momentof[rock] >= 0)
            throw std::invalid_argument("Sphere transport: a moment must follow a field of its own");
        layout.moments[rock].push_back(static_cast<int>(f));
        for (std::size_t c = 0; c < count; ++c)
            if ((*fields[rock])[c] == 0 && (*fields[f])[c] != 0)
                throw std::invalid_argument("Sphere transport: invalid material moment");
    }
    return layout;
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

    const double radius = l.radius_m();
    const std::size_t count = l.subcellcount();
    const auto rotationof = [&](std::uint32_t plate) { return plate == noplate ? Vec3{} : motion.rotation[plate]; };
    displacement_.assign(count, {0, 0});
    maximum_ = 0;
    for (CellId c = 0; c < count; ++c) {
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
        const Vec3 w = rotationof(motion.plate[p]);
        if (mesh::length(w) > 0) {
            r = rotate(w, r);
            moved = true;
        }
        if (!moved) continue;
        displacement_[c] = indexdisplacement(l, c, r);
        maximum_ = std::max({maximum_, std::abs(displacement_[c][0]), std::abs(displacement_[c][1])});
    }
    if (!std::isfinite(maximum_)) throw std::invalid_argument("Sphere transport: displacement not finite");

    // The plates (ascending, noplate last) and their cells.
    std::vector<std::uint32_t> present;
    for (const auto p : motion.plate) present.push_back(p);
    std::sort(present.begin(), present.end());
    present.erase(std::unique(present.begin(), present.end()), present.end());
    groups_.resize(present.size());
    for (std::size_t g = 0; g < present.size(); ++g) groups_[g].plate = present[g];
    for (CellId c = 0; c < count; ++c) {
        const auto g = std::lower_bound(present.begin(), present.end(), motion.plate[l.parent(c)]) - present.begin();
        groups_[static_cast<std::size_t>(g)].cells.push_back(c);
    }

    // Substeps over the subcells and every plate's rotation on its reference domain (the reach
    // dilated by three rings per substep: a substep's profile reads two rings and splats one,
    // so the reference is exact on the reach). A larger count widens the domains: repeat.
    std::vector<std::uint8_t> mark(count, 0);
    std::vector<std::vector<CellId>> domains(groups_.size());
    std::vector<std::vector<std::array<double, 2>>> rigid(groups_.size());
    substeps_ = static_cast<std::uint32_t>(std::max(1.0, std::ceil(maximum_ / maximumsubstep)));
    for (int attempt = 0;; ++attempt) {
        if (attempt == 8) throw std::invalid_argument("Sphere transport: the substeps do not settle");
        double largest = maximum_;
        for (std::size_t g = 0; g < groups_.size(); ++g) {
            Group& group = groups_[g];
            group.reach = dilate(l, group.cells, substeps_, mark);
            domains[g] = dilate(l, group.reach, 3 * substeps_, mark);
            rigid[g].resize(domains[g].size());
            const Vec3 w = rotationof(group.plate);
            for (std::size_t k = 0; k < domains[g].size(); ++k) {
                rigid[g][k] = rotationdisplacement(l, domains[g][k], w);
                largest = std::max({largest, std::abs(rigid[g][k][0]), std::abs(rigid[g][k][1])});
            }
        }
        if (!std::isfinite(largest)) throw std::invalid_argument("Sphere transport: displacement not finite");
        const auto needed = static_cast<std::uint32_t>(std::max(1.0, std::ceil(largest / maximumsubstep)));
        if (needed <= substeps_) { maximum_ = largest; break; }
        substeps_ = needed;
    }

    entries_.resize(count);
    for (CellId c = 0; c < count; ++c)
        entries_[c] = makeentry(l, c, displacement_[c][0] / substeps_, displacement_[c][1] / substeps_);

    // The rigid references: 1 on the domain moved by the plate's rotation alone, read on the reach.
    Layout one;
    one.independent = {0};
    one.moments.resize(1);
    std::vector<std::vector<double>> current(1, std::vector<double>(count, 0)), next(1, std::vector<double>(count, 0));
    std::vector<Entry> entries;
    for (std::size_t g = 0; g < groups_.size(); ++g) {
        Group& group = groups_[g];
        const std::vector<CellId>& domain = domains[g];
        entries.resize(domain.size());
        for (std::size_t k = 0; k < domain.size(); ++k)
            entries[k] = makeentry(l, domain[k], rigid[g][k][0] / substeps_, rigid[g][k][1] / substeps_);
        for (const CellId c : domain) { current[0][c] = 1; mark[c] = 1; }
        const auto entryof = [&](std::size_t k, CellId) -> const Entry& { return entries[k]; };
        for (std::uint32_t step = 0; step < substeps_; ++step)
            substep(l, domain, entryof, one, &mark, false, current, next);
        group.reference.resize(group.reach.size());
        for (std::size_t k = 0; k < group.reach.size(); ++k) group.reference[k] = current[0][group.reach[k]];
        for (const CellId c : domain) { current[0][c] = 0; mark[c] = 0; }
    }
}

std::array<double, 2> Stencil::displacement(CellId subcell) const { return displacement_[subcell]; }

std::vector<PlateLayer> Stencil::transport(const std::vector<const std::vector<double>*>& fields,
                                           const std::vector<int>& momentof) const {
    const Lattice& l = *lattice_;
    const std::size_t count = l.subcellcount();
    const Layout given = validate(count, fields, momentof);
    // Field 0 is the plate's coverage; the given fields follow.
    Layout layout;
    layout.independent.push_back(0);
    for (const int f : given.independent) layout.independent.push_back(f + 1);
    layout.moments.resize(fields.size() + 1);
    for (std::size_t f = 0; f < fields.size(); ++f)
        for (const int m : given.moments[f]) layout.moments[f + 1].push_back(m + 1);

    std::vector<std::vector<double>> current(fields.size() + 1, std::vector<double>(count, 0));
    std::vector<std::vector<double>> next(fields.size() + 1, std::vector<double>(count, 0));
    std::vector<std::uint8_t> mark(count, 0);
    const auto entryof = [&](std::size_t, CellId c) -> const Entry& { return entries_[c]; };
    std::vector<PlateLayer> layers;
    layers.reserve(groups_.size());
    for (const Group& group : groups_) {
        for (const CellId c : group.cells) {
            current[0][c] = 1;
            for (std::size_t f = 0; f < fields.size(); ++f) current[f + 1][c] = (*fields[f])[c];
        }
        for (const CellId c : group.reach) mark[c] = 1;
        for (std::uint32_t step = 0; step < substeps_; ++step)
            substep(l, group.reach, entryof, layout, &mark, true, current, next);

        PlateLayer layer;
        layer.plate = group.plate;
        layer.values.resize(fields.size());
        std::vector<double> transported(fields.size(), 0), kept(fields.size(), 0);
        for (std::size_t k = 0; k < group.reach.size(); ++k) {
            const CellId c = group.reach[k];
            const double area = l.area_m2(c);
            for (std::size_t f = 0; f < fields.size(); ++f) transported[f] += current[f + 1][c] * area;
            const double covered = current[0][c];
            if (!(covered > 0)) continue;
            if (!(group.reference[k] > 0)) throw std::logic_error("Sphere transport: coverage beyond the rigid reference");
            const double coverage = covered / group.reference[k];
            layer.cells.push_back(c);
            layer.coverage.push_back(coverage);
            for (std::size_t f = 0; f < fields.size(); ++f) {
                const double value = current[f + 1][c] / covered;
                layer.values[f].push_back(value);
                kept[f] += value * coverage * area;
            }
        }
        // Per field of its own, the scale that restores the plate's transported amount; its
        // moments take the same scale, so their specific values stay those transported.
        for (const int f : given.independent) {
            if (!(kept[f] > 0)) continue;
            const double scale = transported[f] / kept[f];
            for (auto& v : layer.values[f]) v *= scale;
            for (const int m : given.moments[f])
                for (auto& v : layer.values[m]) v *= scale;
        }
        for (const CellId c : group.reach) {
            mark[c] = 0;
            for (auto& field : current) field[c] = 0;
        }
        layers.push_back(std::move(layer));
    }
    return layers;
}

std::vector<double> Stencil::remap(const std::vector<double>& values) const { return remap(values, {}); }

std::vector<double> Stencil::remap(const std::vector<double>& values,
                                   const std::vector<std::vector<double>*>& moments) const {
    const Lattice& l = *lattice_;
    const std::size_t count = l.subcellcount();
    std::vector<const std::vector<double>*> fields{&values};
    std::vector<int> momentof{-1};
    for (const auto* m : moments) {
        fields.push_back(m);
        momentof.push_back(0);
    }
    const Layout layout = validate(count, fields, momentof);
    std::vector<std::vector<double>> current, next(fields.size(), std::vector<double>(count, 0));
    for (const auto* f : fields) current.push_back(*f);
    std::vector<CellId> all(count);
    std::iota(all.begin(), all.end(), CellId{0});
    const auto entryof = [&](std::size_t, CellId c) -> const Entry& { return entries_[c]; };
    for (std::uint32_t step = 0; step < substeps_; ++step) substep(l, all, entryof, layout, nullptr, false, current, next);
    for (std::size_t m = 0; m < moments.size(); ++m) *moments[m] = std::move(current[m + 1]);
    return std::move(current[0]);
}

} // namespace tectonics::sphere_transport
