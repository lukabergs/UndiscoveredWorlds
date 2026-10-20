#include "simulations/geology/plate_tectonics/tectonic_sphere_transport.hpp"

#include "core/mesh/cube_geometry.hpp"
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable : 4996)  // getenv for UW_THREADS; the engine builds with /WX
#endif
#include "core/parallel_rows.hpp"
#if defined(_MSC_VER)
#pragma warning(pop)
#endif

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


// A source's amounts of one field into its four targets by the field's own profile (the corner's
// gap shared); false when the source holds none of the field.
bool ownamounts(const Lattice& l, CellId c, const Entry& e, const std::vector<double>& values,
                std::array<double, 4>& amount) {
    const double q = values[c];
    if (q == 0) return false;  // the limiter flattens both profiles: nothing moves
    const auto& s = l.profile(c);
    auto px = profile(values[s[0]], values[s[1]], q, values[s[2]], values[s[3]]);
    auto py = profile(values[s[4]], values[s[5]], q, values[s[6]], values[s[7]]);
    const double extent = -minimum(px) - minimum(py);
    if (extent > q) {
        const double factor = q / extent;
        px.linear *= factor; px.quadratic *= factor;
        py.linear *= factor; py.quadratic *= factor;
    }
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
        // A cube corner: the diagonal's share goes to the other three by their weights, which
        // weigh at least 1 - 0.95^2 (a substep moves at most 0.95 cells per axis).
        double others = 0;
        for (int t = 0; t < 4; ++t)
            if (t != missing) others += weight[t];
        const double gap = amount[missing];
        amount[missing] = 0;
        for (int t = 0; t < 4; ++t)
            if (t != missing) amount[t] += gap * (weight[t] / others);
    }
    return true;
}

// How a source moves its fields in one substep of a plate, and at a margin the reference's
// amounts into its four targets (per source of the region; reused between substeps).
enum class Move : std::uint8_t { none, margin, own };
struct Sources {
    std::vector<Move> move;
    std::vector<std::array<double, 4>> amount;
};

// One substep over the sources of region (ascending ids; entryof(k, region[k]) their entries).
// current: dense fields, zero outside region; next: dense, zero on entry and on return. A
// destination outside inside (when given) is dropped, or an error when strict. With a plate's
// reference at the substep's start (dense) and full (1 where the plate's coverage equals it),
// field 0 is the coverage, and a source whose own and profile and neighbour subcells are not all
// full moves every field by the reference's transfer fractions (each field's share of the
// reference); the others move each field of its own by its profile. The result, divided by
// the destinations' areas, replaces current on region.
// Each independent field and its moments run as one task (parallelforrows): a field's amounts
// are added in source order whatever the thread count, so the result is bit-identical to a
// serial scatter.
template <class EntryOf>
void substep(const Lattice& l, const std::vector<CellId>& region, const EntryOf& entryof, const Layout& layout,
             const std::vector<std::uint8_t>* inside, bool strict, const std::vector<double>* reference,
             const std::vector<std::uint8_t>* full, std::vector<std::vector<double>>& current,
             std::vector<std::vector<double>>& next, Sources& sources) {
    const bool plates = reference != nullptr;
    const auto allowed = [&](CellId target) {
        if (!inside || (*inside)[target]) return true;
        if (strict) throw std::logic_error("Sphere transport: rock beyond its plate's reach");
        return false;
    };
    const int count = static_cast<int>(region.size());
    if (plates) {
        sources.move.resize(region.size());
        sources.amount.resize(region.size());
        parallelforrows(0, count - 1, [&](int first, int last) {
            for (int k = first; k <= last; ++k) {
                const CellId c = region[static_cast<std::size_t>(k)];
                Move move = Move::own;
                if (current[0][c] == 0) move = Move::none;  // the plate holds nothing here
                else {
                    const auto& s = l.profile(c);
                    bool intact = (*full)[c] && (*full)[s[0]] && (*full)[s[3]] && (*full)[s[4]] && (*full)[s[7]];
                    for (int slot = 1; slot <= mesh::slotcount && intact; ++slot) {
                        const CellId n = l.neighbour(c, slot);
                        intact = n == nocell || (*full)[n];
                    }
                    const auto source = static_cast<std::size_t>(k);
                    if (!intact && (*reference)[c] > 0
                        && ownamounts(l, c, entryof(source, c), *reference, sources.amount[source]))
                        move = Move::margin;
                }
                sources.move[static_cast<std::size_t>(k)] = move;
            }
        }, 4096);
    }
    const int tasks = static_cast<int>(layout.independent.size());
    parallelforrows(0, tasks - 1, [&](int firsttask, int lasttask) {
        std::array<double, 4> amount{};
        for (int task = firsttask; task <= lasttask; ++task) {
            const int f = layout.independent[static_cast<std::size_t>(task)];
            const std::vector<int>& moments = layout.moments[static_cast<std::size_t>(f)];
            for (std::size_t k = 0; k < region.size(); ++k) {
                const CellId c = region[k];
                const Entry& e = entryof(k, c);
                if (plates) {
                    if (sources.move[k] == Move::none) continue;
                    if (sources.move[k] == Move::margin) {
                        // At a margin every field moves by the reference's transfer fractions:
                        // coverage never exceeds the reference under rigid motion and a value is
                        // a convex combination of its donors'.
                        const std::array<double, 4>& a = sources.amount[k];
                        const double r = (*reference)[c];
                        const auto spread = [&](int g) {
                            const double share = current[g][c] / r;
                            if (share == 0) return;
                            for (int t = 0; t < 4; ++t)
                                if (a[t] != 0 && allowed(e.target[t])) next[g][e.target[t]] += a[t] * share;
                        };
                        spread(f);
                        for (const int m : moments) spread(m);
                        continue;
                    }
                }
                if (!ownamounts(l, c, e, current[f], amount)) continue;
                const double q = current[f][c];
                for (int t = 0; t < 4; ++t) {
                    if (amount[t] == 0 || !allowed(e.target[t])) continue;
                    const CellId target = e.target[t];
                    next[f][target] += amount[t];
                    for (const int m : moments) next[m][target] += amount[t] * (current[m][c] / q);
                }
            }
        }
    }, 1);
    parallelforrows(0, static_cast<int>(current.size()) - 1, [&](int firstfield, int lastfield) {
        for (int f = firstfield; f <= lastfield; ++f) {
            std::vector<double>& now = current[static_cast<std::size_t>(f)];
            std::vector<double>& after = next[static_cast<std::size_t>(f)];
            for (const CellId c : region) {
                now[c] = after[c] / l.area_m2(c);
                after[c] = 0;
            }
        }
    }, 1);
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

    // Each plate's entries on its reach: on its own cells the subcell's motion (the rotation and
    // the offset), elsewhere the plate's rotation alone (the domain's displacement).
    for (std::size_t g = 0; g < groups_.size(); ++g) {
        Group& group = groups_[g];
        const std::vector<CellId>& domain = domains[g];
        group.entries.resize(group.reach.size());
        std::size_t j = 0;
        for (std::size_t k = 0; k < group.reach.size(); ++k) {
            const CellId c = group.reach[k];
            while (domain[j] != c) ++j;
            const auto& d = motion.plate[l.parent(c)] == group.plate ? displacement_[c] : rigid[g][j];
            group.entries[k] = makeentry(l, c, d[0] / substeps_, d[1] / substeps_);
        }
    }

    // The rigid references: 1 on the domain moved by the plate's rotation alone, read on the reach.
    Layout one;
    one.independent = {0};
    one.moments.resize(1);
    std::vector<std::vector<double>> current(1, std::vector<double>(count, 0)), next(1, std::vector<double>(count, 0));
    std::vector<Entry> entries;
    Sources sources;
    for (std::size_t g = 0; g < groups_.size(); ++g) {
        Group& group = groups_[g];
        const std::vector<CellId>& domain = domains[g];
        entries.resize(domain.size());
        for (std::size_t k = 0; k < domain.size(); ++k)
            entries[k] = makeentry(l, domain[k], rigid[g][k][0] / substeps_, rigid[g][k][1] / substeps_);
        for (const CellId c : domain) { current[0][c] = 1; mark[c] = 1; }
        const auto entryof = [&](std::size_t k, CellId) -> const Entry& { return entries[k]; };
        group.domain = domain;
        group.history.clear();
        for (std::uint32_t step = 0; step < substeps_; ++step) {
            substep(l, domain, entryof, one, &mark, false, nullptr, nullptr, current, next, sources);
            if (step + 1 < substeps_)
                for (const CellId c : domain) group.history.push_back(current[0][c]);
        }
        group.reference.resize(group.reach.size());
        for (std::size_t k = 0; k < group.reach.size(); ++k) group.reference[k] = current[0][group.reach[k]];
        for (const CellId c : domain) { current[0][c] = 0; mark[c] = 0; }
    }
}

std::array<double, 2> Stencil::displacement(CellId subcell) const { return displacement_[subcell]; }

namespace {
template <class T> std::size_t allocated(const std::vector<T>& v) { return v.capacity() * sizeof(T); }
} // namespace

std::size_t Lattice::bytes() const {
    return sizeof(*this) + allocated(parent_) + allocated(centre_) + allocated(area_) + allocated(neighbour_)
        + allocated(profile_);
}

std::size_t Stencil::bytes() const {
    std::size_t total = sizeof(*this) + allocated(displacement_) + allocated(groups_);
    for (const Group& g : groups_)
        total += allocated(g.cells) + allocated(g.reach) + allocated(g.domain) + allocated(g.entries)
            + allocated(g.reference) + allocated(g.history);
    return total;
}

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
    std::vector<std::uint8_t> mark(count, 0), full(count, 0);
    std::vector<double> reference(count, 0);
    Sources sources;
    std::vector<PlateLayer> layers;
    layers.reserve(groups_.size());
    for (const Group& group : groups_) {
        for (const CellId c : group.cells) {
            current[0][c] = 1;
            for (std::size_t f = 0; f < fields.size(); ++f) current[f + 1][c] = (*fields[f])[c];
        }
        for (const CellId c : group.reach) mark[c] = 1;
        const auto entryof = [&](std::size_t k, CellId) -> const Entry& { return group.entries[k]; };
        const std::size_t span = group.domain.size();
        for (std::uint32_t step = 0; step < substeps_; ++step) {
            for (std::size_t j = 0; j < span; ++j)
                reference[group.domain[j]] = step == 0 ? 1 : group.history[(step - 1) * span + j];
            for (const CellId c : group.reach) {
                const double cov = current[0][c], r = reference[c];
                full[c] = cov > 0 && cov == r;  // bit-equal: an intact interior
            }
            substep(l, group.reach, entryof, layout, &mark, true, &reference, &full, current, next, sources);
        }
        for (const CellId c : group.domain) reference[c] = 0;
        for (const CellId c : group.reach) full[c] = 0;

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
        // Per field, the scale that restores the plate's transported amount (the coverage
        // divided by the reference changes it by the reference's ripple).
        // Moments take their own scale: with their rock's, an age moment's total drifted by up
        // to 2.7e-4 per update (E2D rework, 2026-10-19, AgeMomentDriftOver100Updates).
        layer.scale.assign(fields.size(), 1.0);
        for (std::size_t f = 0; f < fields.size(); ++f) {
            if (!(kept[f] > 0)) continue;
            const double scale = transported[f] / kept[f];
            layer.scale[f] = scale;
            for (auto& v : layer.values[f]) v *= scale;
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
    std::vector<Entry> entries(count);
    for (CellId c = 0; c < count; ++c)
        entries[c] = makeentry(l, c, displacement_[c][0] / substeps_, displacement_[c][1] / substeps_);
    const auto entryof = [&](std::size_t, CellId c) -> const Entry& { return entries[c]; };
    Sources sources;
    for (std::uint32_t step = 0; step < substeps_; ++step)
        substep(l, all, entryof, layout, nullptr, false, nullptr, nullptr, current, next, sources);
    for (std::size_t m = 0; m < moments.size(); ++m) *moments[m] = std::move(current[m + 1]);
    return std::move(current[0]);
}

} // namespace tectonics::sphere_transport
