#include "simulations/geology/plate_tectonics/tectonic_plumes.hpp"
#include "simulations/geology/plate_tectonics/tectonic_grid.hpp"
#include "core/mesh/mesh_geometry.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace tectonics::plumes {
namespace {
// Stateless seed-salted hash, as the ocean relief: no engine random state is drawn.
double hash(uint32_t index, uint32_t channel, uint32_t seed) {
    uint32_t value = seed ^ 0x5bd1e995U ^ (index * 0x9e3779b9U) ^ (channel * 0x85ebca6bU);
    value ^= value >> 16; value *= 0x7feb352dU;
    value ^= value >> 15; value *= 0x846ca68bU; value ^= value >> 16;
    return value / 4294967295.0;
}
double smooth(double t) { t = std::clamp(t, 0.0, 1.0); return t*t*(3-2*t); }
constexpr double km3_per_myr_to_m3 = 1e15; // 1 km3/yr = 1e9 m3/yr = 1e15 m3/Myr

// Normalized Gaussian weights of the cells around one point on the torus, so a
// volume spread with them is conserved exactly at any resolution.
void spread(std::vector<double>& field, uint32_t w, uint32_t h, double cell_m, const Plume& plume,
            double sigma_km, double volume_m3) {
    if (volume_m3 <= 0) return;
    const double cell_km = cell_m / 1000, width_km = w*cell_km, height_km = h*cell_km;
    const int reach = static_cast<int>(std::ceil(3*sigma_km/cell_km)) + 1;
    const int cx = static_cast<int>(std::floor(plume.x_km/cell_km)), cy = static_cast<int>(std::floor(plume.y_km/cell_km));
    std::vector<std::pair<size_t, double>> weights;
    double total = 0;
    for (int dy = -reach; dy <= reach; ++dy) for (int dx = -reach; dx <= reach; ++dx) {
        const int x = ((cx+dx) % static_cast<int>(w) + static_cast<int>(w)) % static_cast<int>(w);
        const int y = ((cy+dy) % static_cast<int>(h) + static_cast<int>(h)) % static_cast<int>(h);
        const double ex = std::remainder((x+.5)*cell_km - plume.x_km, width_km);
        const double ey = std::remainder((y+.5)*cell_km - plume.y_km, height_km);
        const double weight = std::exp(-.5*(ex*ex+ey*ey)/(sigma_km*sigma_km));
        if (weight < 1e-9) continue;
        weights.push_back({static_cast<size_t>(y)*w + x, weight});
        total += weight;
    }
    if (total <= 0) { // narrower than a cell: all of it into the plume's cell
        const auto x = static_cast<size_t>(std::clamp(cx, 0, static_cast<int>(w) - 1));
        const auto y = static_cast<size_t>(std::clamp(cy, 0, static_cast<int>(h) - 1));
        field[y*w + x] += volume_m3/(cell_m*cell_m);
        return;
    }
    for (const auto& [cell, weight] : weights) field[cell] += volume_m3*weight/total/(cell_m*cell_m);
}
// The sphere: Gaussian weights by great-circle distance over the cells within reach of the
// plume, the volume shared by weight x true area, so the basalt thickness is the Gaussian
// and the volume is conserved exactly.
void spreadsphere(std::vector<double>& field, const grid::TectonicGrid& grid, const Plume& plume,
                  double sigma_km, double volume_m3) {
    if (volume_m3 <= 0) return;
    const mesh::CellGrid& mesh = *grid.mesh();
    const double sigma_m = sigma_km*1000;
    const double spacing_m = std::sqrt(grid.totalarea_m2()/static_cast<double>(grid.cellcount()));
    std::vector<std::pair<grid::CellId, double>> weights;
    double total = 0;
    for (const auto c : mesh::cellswithin(mesh, plume.position, 3*sigma_m + 2*spacing_m)) {
        const double distance = mesh::greatcircle_m(plume.position, mesh.centre(c), mesh.radius_m());
        const double weight = std::exp(-.5*distance*distance/(sigma_m*sigma_m));
        if (weight < 1e-9) continue;
        weights.push_back({c, weight});
        total += weight*grid.area_m2(c);
    }
    if (total <= 0) { // narrower than a cell: all of it into the plume's cell
        const grid::CellId c = mesh.locate(plume.position);
        field[c] += volume_m3/grid.area_m2(c);
        return;
    }
    for (const auto& [cell, weight] : weights) field[cell] += volume_m3*weight/total;
}
}

void validate(const Parameters& p) {
    const auto nonnegative = [](double v) { return std::isfinite(v) && v >= 0; };
    const auto positive = [](double v) { return std::isfinite(v) && v > 0; };
    if (p.count > 64 || !nonnegative(p.tail_flux_km3_yr) || !positive(p.tail_radius_km) ||
        !nonnegative(p.head_volume_km3) || !positive(p.head_radius_km) || !positive(p.head_duration_myr) ||
        !nonnegative(p.birth_window_myr) || !nonnegative(p.swell_age_myr) || !positive(p.swell_radius_km) ||
        !std::isfinite(p.weakening) || p.weakening < 0 || p.weakening > 1 || !positive(p.min_spacing_km))
        throw std::runtime_error("Invalid mantle plume parameters: count <= 64, positive radii, spacing and head "
                                 "duration, nonnegative fluxes, volumes and ages, weakening in 0..1");
}

std::vector<Plume> place(double width_km, double height_km, long seed, const Parameters& p) {
    std::vector<Plume> plumes;
    const auto key = static_cast<uint32_t>(seed);
    // Ordered rejection: candidate k is accepted when it keeps the spacing to every
    // accepted plume. The sequence depends on the seed only.
    for (uint32_t k = 0; plumes.size() < p.count && k < 256*std::max(1U, p.count); ++k) {
        const Plume candidate{hash(k, 0, key)*width_km, hash(k, 1, key)*height_km,
                              hash(k, 2, key)*p.birth_window_myr};
        bool spaced = true;
        for (const auto& other : plumes) {
            const double dx = std::remainder(candidate.x_km - other.x_km, width_km);
            const double dy = std::remainder(candidate.y_km - other.y_km, height_km);
            if (std::hypot(dx, dy) < p.min_spacing_km) { spaced = false; break; }
        }
        if (spaced) plumes.push_back(candidate);
    }
    return plumes;
}

std::vector<Plume> placesphere(double radius_km, long seed, const Parameters& p) {
    std::vector<Plume> plumes;
    const auto key = static_cast<uint32_t>(seed);
    for (uint32_t k = 0; plumes.size() < p.count && k < 256*std::max(1U, p.count); ++k) {
        // hash() reaches 1 inclusive; uniformpoint takes [0, 1) (1 is the same point as 0 in longitude).
        const double u = std::min(hash(k, 0, key), std::nextafter(1.0, 0.0));
        Plume candidate{0, 0, hash(k, 2, key)*p.birth_window_myr, mesh::uniformpoint(u, hash(k, 1, key))};
        bool spaced = true;
        for (const auto& other : plumes)
            if (mesh::greatcircle_rad(candidate.position, other.position)*radius_km < p.min_spacing_km) { spaced = false; break; }
        if (spaced) plumes.push_back(candidate);
    }
    return plumes;
}

Source source(uint32_t w, uint32_t h, double cell_m, long seed, double t0, double t1, const Parameters& p,
              const grid::TectonicGrid* grid) {
    validate(p);
    const bool planar = !grid || grid->kind() == grid::TectonicGrid::Kind::planar;
    if ((planar && (w == 0 || h == 0 || !(cell_m > 0))) || !std::isfinite(t0) || !(t1 > t0))
        throw std::runtime_error("Invalid mantle plume source grid or interval");
    Source result;
    const size_t count = planar ? static_cast<size_t>(w)*h : grid->cellcount();
    result.basalt_m.assign(count, 0.0);
    result.influence.assign(count, 0.0);
    result.swell_age_myr = p.swell_age_myr;
    const auto head_volume = [&](const Plume& plume) {
        // Smoothstep head: the cumulative difference integrates to the head volume.
        return p.head_volume_km3*1e9*(smooth((t1 - plume.birth_myr)/p.head_duration_myr)
                                    - smooth((t0 - plume.birth_myr)/p.head_duration_myr));
    };
    if (!planar) {
        const mesh::CellGrid& mesh = *grid->mesh();
        const double radius_km = mesh.radius_m()/1000;
        for (const auto& plume : placesphere(radius_km, seed, p)) {
            const double active = std::max(0.0, t1 - std::max(t0, plume.birth_myr));
            spreadsphere(result.basalt_m, *grid, plume, p.tail_radius_km, p.tail_flux_km3_yr*km3_per_myr_to_m3*active);
            spreadsphere(result.basalt_m, *grid, plume, p.head_radius_km, head_volume(plume));
            if (t1 <= plume.birth_myr) continue;
            for (grid::CellId c = 0; c < count; ++c) {
                const double distance_km = mesh::greatcircle_rad(plume.position, mesh.centre(c))*radius_km;
                auto& influence = result.influence[c];
                influence = std::max(influence, std::exp(-.5*distance_km*distance_km/(p.swell_radius_km*p.swell_radius_km)));
            }
        }
        return result;
    }
    const double cell_km = cell_m/1000, width_km = w*cell_km, height_km = h*cell_km;
    for (const auto& plume : place(width_km, height_km, seed, p)) {
        const double active = std::max(0.0, t1 - std::max(t0, plume.birth_myr));
        spread(result.basalt_m, w, h, cell_m, plume, p.tail_radius_km, p.tail_flux_km3_yr*km3_per_myr_to_m3*active);
        spread(result.basalt_m, w, h, cell_m, plume, p.head_radius_km, head_volume(plume));
        if (t1 <= plume.birth_myr) continue;
        for (uint32_t y = 0; y < h; ++y) for (uint32_t x = 0; x < w; ++x) {
            const double ex = std::remainder((x+.5)*cell_km - plume.x_km, width_km);
            const double ey = std::remainder((y+.5)*cell_km - plume.y_km, height_km);
            auto& influence = result.influence[static_cast<size_t>(y)*w + x];
            influence = std::max(influence, std::exp(-.5*(ex*ex+ey*ey)/(p.swell_radius_km*p.swell_radius_km)));
        }
    }
    return result;
}
}
