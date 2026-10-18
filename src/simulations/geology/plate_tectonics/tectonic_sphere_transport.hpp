#ifndef TECTONIC_SPHERE_TRANSPORT_HPP
#define TECTONIC_SPHERE_TRANSPORT_HPP

// The rock transport on the cube (W2-E, step E2D; docs/architecture/terrain-roadmap.md,
// "W2-E"). Not wired into the lithosphere yet (E3).
//
// The subcell lattice is the cube of s N cells per face edge (uw-cube-1 ids at s N,
// docs/architecture/sphere-mesh.md section 2): subcell (f, I, J) has the parent cell
// (f, I / s, J / s) on the same face and its exact area. Once per update a stencil is built:
// each subcell centre r is moved by the tangent offset of its parent cell (the exponential
// map on the sphere) and then by the rotation exp([w]x) of its parent's plate, the
// destination is located in extended face coordinates of the source face, and the index
// displacement D gives ceil(max |D| / 0.95) substeps (max over subcells and both axes). Per
// substep the kernel is the torus kernel's (tectonic_transport.cpp): the piecewise-parabolic
// profiles along the face axes (neighbours ahead, across face edges in the neighbour's frame),
// the limiter, a 2 x 2 splat in the source face's index space with the targets beyond the face
// from the neighbour table; at a cube corner the missing diagonal's amount goes to the other
// three targets in proportion to their weights. Amounts are value x source subcell area and the
// result is the amount divided by the destination's area, so the integral of value x area is
// conserved. The stencil serves every field and substep; the scatter is serial in id order.

#include "core/mesh/cell_grid.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

namespace tectonics::sphere_transport {

using CellId = mesh::CellId;
using mesh::Vec3;
constexpr CellId nocell = mesh::nocell;
constexpr std::uint32_t noplate = 0xffffffffu;  // a cell without a plate: no rotation
constexpr double maximumsubstep = 0.95;         // lattice cells per substep and axis

class Lattice {
public:
    // N >= 1 parent cells per face edge, s >= 1 subcells per parent edge, s N >= 2.
    static std::shared_ptr<const Lattice> make(int facecells, int subdivision, double radius_m);

    int facecells() const { return facecells_; }        // N
    int subdivision() const { return subdivision_; }    // s
    int subfacecells() const { return subfacecells_; }  // s N
    std::size_t cellcount() const { return cellcount_; }        // 6 N^2 parents
    std::size_t subcellcount() const { return subcellcount_; }  // 6 (s N)^2
    double radius_m() const { return radius_m_; }

    CellId parent(CellId subcell) const { return parent_[subcell]; }
    const Vec3& centre(CellId subcell) const { return centre_[subcell]; }  // unit vector
    double area_m2(CellId subcell) const { return area_[subcell]; }
    // Slot 0 (the subcell itself) to 8 (mesh D8 codes in the face frame); nocell for the
    // diagonal of a corner subcell that would cross two face edges.
    CellId neighbour(CellId subcell, int slot) const;
    // The profile cells along the face axes: west 2, west 1, east 1, east 2, north 2, north 1,
    // south 1, south 2 (the second step continues straight in the first neighbour's frame).
    const std::array<CellId, 8>& profile(CellId subcell) const { return profile_[subcell]; }

private:
    Lattice() = default;
    int facecells_ = 0, subdivision_ = 0, subfacecells_ = 0;
    std::size_t cellcount_ = 0, subcellcount_ = 0;
    double radius_m_ = 0;
    std::vector<CellId> parent_;
    std::vector<Vec3> centre_;
    std::vector<double> area_;
    std::vector<std::array<CellId, 8>> neighbour_;  // slots 1..8
    std::vector<std::array<CellId, 8>> profile_;
};

// One update's motion. Rotation vectors are the update's rotations (axis times angle in
// radians, the plate's angular velocity times the step). Offsets are tangent vectors in metres
// at the parent cell's centre (the deformation and rift offsets), projected onto each
// subcell's tangent plane; empty means none.
struct Motion {
    std::vector<Vec3> rotation;             // per plate
    std::vector<std::uint32_t> plate;       // per parent cell: an index into rotation, or noplate
    std::vector<Vec3> offset_m;             // per parent cell, or empty
};

class Stencil {
public:
    // Throws std::invalid_argument on a size mismatch, a non-finite input, a plate index out
    // of range, or a destination a quarter turn or more from the source face's centre.
    Stencil(std::shared_ptr<const Lattice> lattice, const Motion& motion);

    const Lattice& lattice() const { return *lattice_; }
    std::uint32_t substeps() const { return substeps_; }
    double maximumdisplacement() const { return maximum_; }  // lattice cells per update
    // The index displacement of a subcell over the whole update (columns east, rows south in
    // the source face's frame).
    std::array<double, 2> displacement(CellId subcell) const;

    // A density (amount per area) on the subcells through every substep. Values must be
    // finite and non-negative (std::invalid_argument otherwise).
    std::vector<double> remap(const std::vector<double>& values) const;
    // The same, with material moments (e.g. thickness x age) moved by the rock's transfer
    // fractions and replaced in place: both integrals are conserved and each moment's specific
    // value stays within its donors' range. A moment must be finite, non-negative and zero
    // where the rock is zero.
    std::vector<double> remap(const std::vector<double>& values,
                              const std::vector<std::vector<double>*>& moments) const;

private:
    struct Entry {
        std::array<CellId, 4> target;  // (0, 0), (1, 0), (0, 1), (1, 1); nocell: a corner's gap
        std::array<double, 2> wx, cx, qx, wy, cy, qy;  // weight, linear and quadratic terms per axis
    };
    void substep(const std::vector<double>& values, const std::vector<std::vector<double>>& moments,
                 std::vector<double>& next, std::vector<std::vector<double>>& nextmoments) const;

    std::shared_ptr<const Lattice> lattice_;
    std::uint32_t substeps_ = 1;
    double maximum_ = 0;
    std::vector<std::array<double, 2>> displacement_;
    std::vector<Entry> entries_;
};

} // namespace tectonics::sphere_transport
#endif
