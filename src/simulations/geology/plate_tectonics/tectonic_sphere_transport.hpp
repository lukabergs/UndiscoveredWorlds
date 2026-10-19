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
//
// The splat alone does not keep a rigid plate uniform: the index displacement of a rotation
// varies over a face (a ripple of about its gradient per update) and three face frames meet at
// 120 degrees at a cube corner (about a third there). So the engine moves each plate
// separately (transport, the E2D decision of 2026-10-19). A plate's rock moves with the plate's
// own motion over its whole reach and every substep: its rotation, plus the tangent offset of
// the cells the plate owns (offsets are the plate's deformation; elsewhere the rotation alone).
// A cell's amount of a field is value x coverage x area. The plate's rigid reference is the
// transport of 1 everywhere under the plate's rotation alone (kept after every substep). A
// source subcell whose coverage equals the reference on itself, its eight neighbours and its
// profile cells (an intact interior) moves the coverage and each field by its own profile, as
// the torus does (resolved relief is kept; a smooth extremum may exceed its donors slightly);
// any other source (a margin, a deforming or overlapping region) moves every field by the
// reference's transfer fractions, each in proportion to its share of the reference there, so
// a value is a convex combination of donor values and a rigid plate's coverage never exceeds
// its reference. The transported value is the transported amount over the transported
// coverage, and the coverage is divided by the reference: a rigid interior keeps coverage 1
// and a uniform field uniform exactly. Summed over plates, coverage above 1 is the overlap of
// converging plates and below 1 a gap; within one plate, coverage above 1 is shortening by
// its offsets (and below 1 extension), not overlap.

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
    std::size_t bytes() const;  // the tables' allocated sizes

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

// One plate's rock after an update, on the subcells its coverage reaches. A cell's amount of
// a field is value x coverage x area.
struct PlateLayer {
    std::uint32_t plate = 0;          // an index into Motion::rotation, or noplate
    std::vector<CellId> cells;        // ascending: transported coverage > 0
    // Transported coverage over the plate's rigid reference: 1 in a rigid interior, the covered
    // share at a margin, above 1 where the plate's offsets shorten it. Summed over plates, above
    // 1 is overlap (convergence), below 1 a gap (divergence).
    std::vector<double> coverage;
    // Per field, aligned with cells: the transported amount over the transported coverage (a
    // value per covered area), scaled per plate and field so that sum value x coverage x area
    // equals the plate's transported amount of the field.
    std::vector<std::vector<double>> values;
    // Per field, that scale (1 where the reference's ripple cancels; a value's excess beyond its
    // donors at a margin is at most |scale - 1|).
    std::vector<double> scale;
};

namespace detail {
struct Entry {
    std::array<CellId, 4> target;  // (0, 0), (1, 0), (0, 1), (1, 1); nocell: a corner's gap
    std::array<double, 2> wx, cx, qx, wy, cy, qy;  // weight, linear and quadratic terms per axis
};
} // namespace detail

class Stencil {
public:
    // Throws std::invalid_argument on a size mismatch, a non-finite input, a plate index out
    // of range, or a destination a quarter turn or more from the source face's centre.
    Stencil(std::shared_ptr<const Lattice> lattice, const Motion& motion);

    const Lattice& lattice() const { return *lattice_; }
    std::uint32_t substeps() const { return substeps_; }
    // Lattice cells per update, over the subcells, the plates' reaches and their references.
    double maximumdisplacement() const { return maximum_; }
    // The index displacement of a subcell over the whole update under its own plate's motion
    // (columns east, rows south in the source face's frame).
    std::array<double, 2> displacement(CellId subcell) const;
    std::size_t bytes() const;  // the members' allocated sizes

    // The engine's transport: every plate's rock moved separately (the header comment). fields:
    // densities on the subcells (amount per area), each subcell's rock belonging to its
    // parent's plate. momentof: empty, or per field the index of the rock field whose transfer
    // fractions it uses inside a plate (a material moment such as thickness x age; its specific
    // value stays within its donors' range) or -1 for a field of its own. Layers in ascending
    // plate order, noplate last; a plate without cells has none. Values finite and non-negative,
    // a moment zero where its rock is (std::invalid_argument otherwise).
    std::vector<PlateLayer> transport(const std::vector<const std::vector<double>*>& fields,
                                      const std::vector<int>& momentof = {}) const;

    // The raw kernel over the whole lattice, each subcell moving with its own plate's motion and
    // no per-plate normalisation (the crossing and conservation checks): a density through every
    // substep; the integral of value x area is conserved.
    std::vector<double> remap(const std::vector<double>& values) const;
    // The same, with material moments moved by the rock's transfer fractions and replaced in
    // place: both integrals are conserved and each moment's specific value stays within its
    // donors' range.
    std::vector<double> remap(const std::vector<double>& values,
                              const std::vector<std::vector<double>*>& moments) const;

private:
    // A plate's cells, its reach (the subcells its rock can arrive at: the cells dilated by one
    // ring per substep), the plate's entries on the reach, its rigid reference there, and the
    // reference's domain (the reach dilated by three rings per substep) with the reference after
    // each substep but the last on it.
    struct Group {
        std::uint32_t plate = 0;
        std::vector<CellId> cells, reach, domain;
        std::vector<detail::Entry> entries;
        std::vector<double> reference, history;
    };

    std::shared_ptr<const Lattice> lattice_;
    std::uint32_t substeps_ = 1;
    double maximum_ = 0;
    std::vector<std::array<double, 2>> displacement_;
    std::vector<Group> groups_;
};

} // namespace tectonics::sphere_transport
#endif
