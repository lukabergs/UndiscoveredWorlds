#include "simulations/geology/plate_tectonics/tectonic_columns.hpp"
#include "core/mesh/mesh_geometry.hpp"
#include "simulations/geology/plate_tectonics/tectonic_grid.hpp"
#include "simulations/geology/plate_tectonics/tectonic_thermal.hpp"
#include "simulations/landscape/depressions.hpp"
#include "simulations/landscape/flow_graph.hpp"
#include "simulations/landscape/sphere_metric.hpp"
#include "simulations/landscape/stream_power.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <iomanip>
#include <numeric>
#include <limits>

namespace tectonics::columns {
namespace {
constexpr double ocean_reference = (mantle_density - oceanic_density) * ridge_crust_m
    + (mantle_density - water_density) * ridge_depth_m;
// Continental mantle has a different reference buoyancy than a young ridge.
constexpr double continental_mantle_load =
    (mantle_density - continental_density) * reference_continent_m - ocean_reference;
void validate(const Cell& c) {
    int field = 0;
    for (double v : {c.continent, c.ocean, c.sediment, c.continental_area,
                     c.oceanic_area, c.continental_age_volume, c.oceanic_age_volume,
                     c.thermal_rejuvenation_volume, c.flexure_m}) {
        if (!std::isfinite(v) || v < 0) {
            std::ostringstream message;
            message << "Invalid crust column field " << field << ": " << std::setprecision(17) << v;
            throw std::runtime_error(message.str());
        }
        ++field;
    }
    if (std::abs(c.continental_area + c.oceanic_area - 1) > 1e-8)
        throw std::runtime_error("Crust area fractions must sum to one");
}
// Removes up to `mass` (kg/m²) of rock from a column: sediment first, then
// continental, then oceanic crust. Age volumes and the ocean's thermal
// rejuvenation scale with the volume left. Returns the mass removed; the part
// taken from oceanic crust is added to *ocean_mass.
double removeRock(Cell& c, double mass, double* ocean_mass = nullptr) {
    double remaining=mass,removed=0;
    const auto erode=[&](double& volume,double density,double* age) {
        const double taken=std::min(volume,std::max(0.0,remaining)/density);
        if (age) *age*=volume>0?(volume-taken)/volume:0;
        volume-=taken;remaining-=taken*density;removed+=taken*density;
        return taken*density;
    };
    erode(c.sediment,sediment_density,nullptr);
    erode(c.continent,continental_density,&c.continental_age_volume);
    const double deficit=c.ocean>0?c.thermal_rejuvenation_volume/c.ocean:0;
    const double ocean=erode(c.ocean,oceanic_density,&c.oceanic_age_volume);
    if (ocean_mass) *ocean_mass+=ocean;
    c.thermal_rejuvenation_volume=c.ocean*deficit;
    return removed;
}
// Relaxation time of the continental area cap (continental_min_thickness_m), Myr.
constexpr double area_cap_response_myr=5;
// Relaxation time of the upper area cap (continental_max_thickness_m), Myr: gravitational
// spreading of overthickened crust, slower than the lower cap's cleanup.
constexpr double area_upper_cap_response_myr=10;
// Oceanic area below which the upper cap converts the rest of a cell at once (share of the cell).
constexpr double upper_cap_residual_area=1e-6;
// Largest orogenic collapse coefficient per face and sub-step (D dt / cell²): at most
// 0.4 of a cell's excess leaves through its four faces, so no source drops below the collapse
// threshold (the maximum thickness, or the onset when set; with the creep below the onset or
// the margin cap the lowest threshold is the maximum, and both regimes share the limit).
constexpr double collapse_step_limit=.1;
// Response time of the continental growth controller to the reservoir's error, Myr.
constexpr double growth_response_myr=50;
// The unit tangent of v at the unit vector r.
mesh::Vec3 tangent(const mesh::Vec3& v, const mesh::Vec3& r) { return mesh::normalise(v - mesh::dot(v, r) * r); }
}

Metric::Metric(const grid::TectonicGrid* grid, uint32_t width, uint32_t height, double cell_m) {
    if (!grid) {
        owned_ = grid::TectonicGrid::planar(width, height, width * cell_m / 1000.0);
        grid = owned_.get();
    }
    grid_ = grid;
    planar_ = grid->kind() == grid::TectonicGrid::Kind::planar;
    if (planar_) {
        if (grid->width() != width || grid->height() != height || !std::isfinite(cell_m) || cell_m <= 0)
            throw std::invalid_argument("Metric: the planar scaffold does not match the grid");
        width_cells_ = width; height_cells_ = height;
        cell_m_ = cell_m;
        cell_area_m2_ = cell_m * cell_m;
        reference_area_m2_ = cell_area_m2_;
        return;
    }
    const size_t count = grid->cellcount();
    const mesh::CellGrid& mesh = *grid->mesh();
    area_.resize(count); width_.resize(count); scale_.resize(count);
    edges_.assign(4 * count, UINT32_MAX);
    const auto& table = grid->edges();
    weight_.resize(table.size());
    for (size_t e = 0; e < table.size(); ++e) {
        const auto& edge = table[e];
        weight_[e] = edge.length_m / edge.distance_m;
        edges_[4 * static_cast<size_t>(edge.a) + edge.direction] = static_cast<uint32_t>(e);
        for (uint8_t k = 0; k < grid::directioncount; ++k)
            if (grid->neighbour(edge.b, k) == edge.a) edges_[4 * static_cast<size_t>(edge.b) + k] = static_cast<uint32_t>(e);
    }
    reference_area_m2_ = std::numeric_limits<double>::infinity();
    east_.resize(count); south_.resize(count);
    for (grid::CellId c = 0; c < count; ++c) {
        area_[c] = grid->area_m2(c);
        width_[c] = std::sqrt(area_[c]);
        for (uint8_t k = 0; k < grid::directioncount; ++k) {
            const uint32_t e = edges_[4 * static_cast<size_t>(c) + k];
            if (e == UINT32_MAX) throw std::logic_error("Metric: a cell face without an edge");
            reference_area_m2_ = std::min(reference_area_m2_, area_[c] / weight_[e]);
        }
        // The frame: x toward the east neighbour (central difference), y orthogonal toward the south one.
        const mesh::Vec3 r = mesh.centre(c);
        east_[c] = tangent(mesh.centre(grid->neighbour(c, grid::east)) - mesh.centre(grid->neighbour(c, grid::west)), r);
        const mesh::Vec3 down = mesh.centre(grid->neighbour(c, grid::south)) - mesh.centre(grid->neighbour(c, grid::north));
        south_[c] = mesh::normalise(down - mesh::dot(down, east_[c]) * east_[c] - mesh::dot(down, r) * r);
    }
    for (grid::CellId c = 0; c < count; ++c) scale_[c] = reference_area_m2_ / area_[c];
}

size_t Metric::edge(grid::CellId c, uint8_t k) const {
    if (!planar_) return edges_[4 * static_cast<size_t>(c) + k];
    // Planar: edge 2 i is (i, east of i), 2 i + 1 is (i, south of i).
    switch (k) {
    case grid::east: return 2 * static_cast<size_t>(c);
    case grid::south: return 2 * static_cast<size_t>(c) + 1;
    case grid::west: return 2 * static_cast<size_t>(grid_->neighbour(c, grid::west));
    default: return 2 * static_cast<size_t>(grid_->neighbour(c, grid::north)) + 1;
    }
}

std::array<Metric::Neighbour, 8> Metric::eight(grid::CellId c) const {
    std::array<Neighbour, 8> result;
    if (planar_) {
        const grid::CellId n = grid_->neighbour(c, grid::north), s = grid_->neighbour(c, grid::south);
        const grid::CellId cells[8] = {grid_->neighbour(n, grid::west), n, grid_->neighbour(n, grid::east),
            grid_->neighbour(c, grid::west), grid_->neighbour(c, grid::east),
            grid_->neighbour(s, grid::west), s, grid_->neighbour(s, grid::east)};
        for (int k = 0; k < 8; ++k)
            result[k] = {cells[k], cell_m_ * (k == 0 || k == 2 || k == 5 || k == 7 ? std::sqrt(2.) : 1.)};
        return result;
    }
    constexpr int slots[8] = {8, 1, 2, 7, 3, 6, 5, 4};
    const mesh::CellGrid& mesh = *grid_->mesh();
    for (int k = 0; k < 8; ++k) {
        const std::int32_t n = mesh.neighbour(c, slots[k]);
        if (n >= 0) result[k] = {static_cast<grid::CellId>(n), mesh.distance_m(c, slots[k])};
    }
    return result;
}

mesh::Vec3 Metric::vector(grid::CellId c, const Displacement& d) const {
    if (planar_) return {d.x, d.y, 0.0};
    return d.x * east_[c] + d.y * south_[c];
}

Displacement Metric::local(grid::CellId c, const mesh::Vec3& v) const {
    if (planar_) return {v.x, v.y};
    return {mesh::dot(v, east_[c]), mesh::dot(v, south_[c])};
}

Displacement Metric::inframe(grid::CellId to, grid::CellId from, const Displacement& d) const {
    if (planar_ || to == from) return d;
    return local(to, vector(from, d));
}

double Metric::closing(size_t e, const Displacement& a, const Displacement& b) const {
    const auto& edge = grid_->edges()[e];
    if (planar_) return edge.axis == 0 ? a.x - b.x : a.y - b.y;
    return mesh::dot(vector(edge.a, a) - vector(edge.b, b), edge.normal);
}

std::array<double, 2> Metric::unitposition(grid::CellId c) const {
    if (planar_) return {(c % width_cells_ + .5) / width_cells_, (c / width_cells_ + .5) / height_cells_};
    constexpr double pi = 3.14159265358979323846;
    const mesh::CellGrid& mesh = *grid_->mesh();
    const double u = (mesh.longitude_rad(c) + pi) / (2 * pi), v = (pi / 2 - mesh.latitude_rad(c)) / pi;
    return {std::min(u, std::nextafter(1.0, 0.0)), std::clamp(v, 0.0, std::nextafter(1.0, 0.0))};
}

std::array<double, 2> Metric::extent_km() const {
    if (planar_) return {width_cells_ * cell_m_ / 1000, height_cells_ * cell_m_ / 1000};
    constexpr double pi = 3.14159265358979323846;
    const double radius_km = grid_->mesh()->radius_m() / 1000;
    return {2 * pi * radius_km, pi * radius_km};
}
double shortening_m(const Strain& s, double cell_m) {
    const double least = (s.xx + s.yy) / 2 - std::hypot((s.xx - s.yy) / 2, s.xy);
    return std::max(0.0, -least) * cell_m;
}

double shortening_direction(const Strain& s) {
    constexpr double pi = 3.14159265358979323846;
    if (std::hypot((s.xx - s.yy) / 2, s.xy) == 0) return 0;
    // The largest component's axis is atan2(2 xy, xx - yy) / 2; the least is normal to it.
    double angle = std::atan2(2 * s.xy, s.xx - s.yy) / 2 + pi / 2;
    if (angle >= pi) angle -= pi;
    return angle < 0 ? angle + pi : angle;
}

double mass_per_area(const Cell& c) {
    return continental_density * c.continent + oceanic_density * c.ocean
        + sediment_density * c.sediment;
}
double felsic_mass_per_area(const Cell& c) {
    return continental_density * c.continent + sediment_density * c.sediment;
}
double age_myr(const Cell& c) {
    return (c.continent + c.ocean) > 0
        ? (c.continental_age_volume + c.oceanic_age_volume) / (c.continent + c.ocean) : 0;
}
double mountain_force_n_m(const Cell& c, bool sediment_feedback, double root_kg_m3) {
    // Airy-compensated continental column, relative to unthickened crust:
    // delta GPE = g*rho_c*(1-rho_c/rho_m)*(H²-H0²)/2, in N/m.
    // Effective mixed-cell column: no division by vanishing area fractions.
    // This reduced crust-only estimate omits thermal mantle GPE and flexure.
    // A lithospheric root scales the elevation per metre of crust by
    // (rho_m - rho_c - root) / (rho_m - rho_c); the GPE g rho_c z (H + H0) / 2 scales
    // with the elevation, so both branches take that factor (the same crust-only
    // estimate; the root's own GPE is not resolved).
    const double root_scale = (mantle_density - continental_density - root_kg_m3)
        / (mantle_density - continental_density);
    const double reference = reference_continent_m * c.continental_area;
    if (sediment_feedback) {
        // Two-layer Airy GPE for continental basement plus its sediment cover,
        // relative to the same unthickened reference column, N/m; exact for the
        // crust-only column (root 0). With a root, the (500 - root) / 500 factor is a
        // model choice: the crust-only estimate scaled with the elevation; a real dense
        // root lowers the GPE further. Mixed-cell sediment participation follows the
        // continental area.
        const double s=c.sediment*c.continental_area;
        const double z=((mantle_density-continental_density)*(c.continent-reference)
            +(mantle_density-sediment_density)*s)/mantle_density;
        const double base=z-c.continent-s,interface_z=z-s;
        const double gpe=sediment_density*z*z+(continental_density-sediment_density)*interface_z*interface_z
            +(mantle_density-continental_density)*(base*base-reference*reference);
        return .5*9.80665*std::max(0.0,gpe)*root_scale;
    }
    return .5 * 9.80665 * continental_density * (1 - continental_density / mantle_density)
        * std::max(0.0, c.continent * c.continent - reference * reference) * root_scale;
}
double elevation_m(const Cell& c, double thermal_strength, double root_kg_m3, double root_max_thickness_m,
        double root_delamination_width_m) {
    const double ocean_age = c.ocean > 0 ? std::max(0.0,
        (c.oceanic_age_volume-c.thermal_rejuvenation_volume)/c.ocean) : 0;
    const double load = ocean_reference + continental_mantle_load * c.continental_area
        + (mantle_density - water_density) * thermal_strength
            * thermal::subsidence_m(ocean_age) * c.oceanic_area;
    double buoyancy = (mantle_density - continental_density) * c.continent
        + (mantle_density - oceanic_density) * c.ocean
        + (mantle_density - sediment_density) * c.sediment - load;
    // Lithospheric root: crust thicker than the (area-weighted) reference rides on a
    // denser mantle root, thinner crust on a lighter one. Exact zero when root is 0.
    if (root_max_thickness_m > 0 && root_delamination_width_m > 0) {
        // Delamination ramp: the whole root's load scales by 1 up to the cap, falls
        // linearly to 0 at cap + width (per-area crust thickness h; mixed cells by area
        // like the reference), and is gone above, where the crust-only Airy height
        // returns. Continuous with the frozen cap at h = cap. No continental area: no
        // root load (the frozen cap's min(continent, 0) as well).
        double factor = 0;
        if (c.continental_area > 0) {
            const double h = c.continent / c.continental_area;
            factor = h <= root_max_thickness_m ? 1
                : h >= root_max_thickness_m + root_delamination_width_m ? 0
                : 1 - (h - root_max_thickness_m) / root_delamination_width_m;
        }
        buoyancy -= root_kg_m3 * factor * (c.continent - reference_continent_m * c.continental_area);
    } else {
        // Frozen cap (root_max_thickness_m > 0, width 0): the root stops growing above that
        // crust thickness (area-weighted like the reference), where the crust-only slope
        // returns.
        const double rooted_crust = root_max_thickness_m > 0
            ? std::min(c.continent, root_max_thickness_m * c.continental_area) : c.continent;
        buoyancy -= root_kg_m3 * (rooted_crust - reference_continent_m * c.continental_area);
    }
    // Hydrostatic compensation includes water loading below sea level.
    return buoyancy / (buoyancy >= 0 ? mantle_density : mantle_density - water_density) - c.flexure_m;
}
Cell seed(double z, bool continental, double age, double root_kg_m3, double root_max_thickness_m,
        double root_delamination_width_m) {
    if (!std::isfinite(z) || !std::isfinite(age) || age < 0 || !std::isfinite(root_kg_m3) || root_kg_m3 < 0
        || root_kg_m3 >= mantle_density - continental_density || !valid_root_max_thickness(root_max_thickness_m)
        || !valid_root_delamination_width(root_max_thickness_m, root_delamination_width_m))
        throw std::runtime_error("Invalid initial column elevation, age or root load");
    Cell c;
    if (continental) {
        c.continental_area = 1;
        // Above the height of a column of the cap's thickness the root no longer grows:
        // the crust-only slope from there (the inverse of elevation_m).
        const double cap_z = root_max_thickness_m > 0 ? (root_max_thickness_m - reference_continent_m)
            * (mantle_density - continental_density - root_kg_m3) / mantle_density : 0;
        // Delamination ramp (the inverse of elevation_m's ramp branch; root 0 has none).
        const bool ramp = root_max_thickness_m > 0 && root_delamination_width_m > 0 && root_kg_m3 > 0;
        // Height of a column at cap + width: the crust-only Airy height, root gone.
        const double top_z = ramp ? (root_max_thickness_m + root_delamination_width_m - reference_continent_m)
            * (mantle_density - continental_density) / mantle_density : 0;
        if (ramp && z >= top_z)
            c.continent = reference_continent_m + z * mantle_density / (mantle_density - continental_density);
        else if (ramp && z > cap_z) {
            // x = H - 35 km, k = cap - 35 km, w = width: rho_m z = (rho_m - rho_c) x
            // - root (k + w - x) x / w, a quadratic a x^2 + b x - rho_m z = 0 with a > 0
            // and z > 0: one positive root (the column's thickness), in the stable form.
            const double k = root_max_thickness_m - reference_continent_m, w = root_delamination_width_m;
            const double a = root_kg_m3 / w;
            const double b = mantle_density - continental_density - root_kg_m3 * (k + w) / w;
            const double q = mantle_density * z;
            const double root_d = std::sqrt(b * b + 4 * a * q);
            c.continent = reference_continent_m + (b >= 0 ? 2 * q / (b + root_d) : (root_d - b) / (2 * a));
        } else if (root_max_thickness_m > 0 && z > cap_z)
            c.continent = root_max_thickness_m + (z - cap_z) * mantle_density / (mantle_density - continental_density);
        else
            c.continent = reference_continent_m
                + z * (z >= 0 ? mantle_density : mantle_density - water_density)
                    / (mantle_density - continental_density - root_kg_m3);
        if (!(c.continent > 0))
            throw std::runtime_error("Continental column has no crust: initial elevation below the root's crust-zero depth");
        c.continental_age_volume = c.continent * age;
    } else {
        // Ocean depth follows age/cooling; do not infer implausibly thick basalt
        // to reproduce the old independently generated height noise.
        c.oceanic_area = 1;
        c.ocean = ridge_crust_m;
        c.oceanic_age_volume = c.ocean * age;
    }
    validate(c);
    return c;
}
void Model::initialize(uint32_t width, uint32_t height, double cell_m, std::vector<Cell> cells,
        uint32_t subdivision, long seed, const ocean_relief::Parameters& relief, bool preserve_initial,
        double root_kg_m3, double root_max_thickness_m, double root_delamination_width_m,
        const grid::TectonicGrid* grid) {
    // Transport, reconciliation, implicit rivers and the newborn crust's noise are planar until E3.
    grid::requireplanar(grid, "columns");
    if (width < 2 || height < 2 || cells.size() != static_cast<size_t>(width) * height
        || !std::isfinite(cell_m) || cell_m <= 0 || (subdivision != 1 && subdivision != 2 && subdivision != 4))
        throw std::runtime_error("Invalid column grid");
    if (!std::isfinite(root_kg_m3) || root_kg_m3 < 0 || root_kg_m3 >= mantle_density - continental_density)
        throw std::runtime_error("Continental root load must be in 0..500 kg/m^3 (500 excluded)");
    if (!valid_root_max_thickness(root_max_thickness_m))
        throw std::runtime_error("Continental root maximum thickness must be 0 (no cap) or above 35000 m");
    if (!valid_root_delamination_width(root_max_thickness_m, root_delamination_width_m))
        throw std::runtime_error("Continental root delamination width must be 0, or finite and positive "
                                 "with a root maximum thickness");
    metric_ = Metric(grid, width, height, cell_m);
    width_ = width; height_ = height; cell_m_ = cell_m; cells_ = std::move(cells);
    seed_=seed;root_kg_m3_=root_kg_m3;root_max_thickness_m_=root_max_thickness_m;time_myr_=0;
    root_delamination_width_m_=root_delamination_width_m;
    newborn_ocean_crust_ = ocean_relief::crust(width,height,cell_m,seed,relief);
    for (size_t i=0; i<cells_.size(); ++i) {
        auto& c=cells_[i];
        const double ratio=preserve_initial?1:newborn_ocean_crust_[i]/ridge_crust_m;
        c.ocean*=ratio; c.oceanic_age_volume*=ratio; c.thermal_rejuvenation_volume*=ratio;
    }
    budget_ = {}; budget_.enabled = true;
    uplift_rate_m_myr_.assign(cells_.size(), 0.0);
    strain_.assign(cells_.size(), Strain{});
    for (grid::CellId i = 0; i < cells_.size(); ++i) {
        const auto& c = cells_[i];
        validate(c);
        budget_.initial_mass_kg += metric_.mass(mass_per_area(c), i);
        budget_.initial_felsic_kg += metric_.mass(felsic_mass_per_area(c), i);
        budget_.initial_continent_kg += metric_.mass(c.continent * continental_density, i);
        budget_.initial_continental_area_m2 += metric_.mass(c.continental_area, i);
    }
    budget_.mass_kg = budget_.initial_mass_kg;
    budget_.felsic_kg = budget_.initial_felsic_kg;
    budget_.continent_kg = budget_.initial_continent_kg;
    subdivision_ = subdivision;
    transport_cells_.resize(static_cast<size_t>(width)*height*subdivision*subdivision);
    for (uint32_t y=0; y<height*subdivision; ++y) for (uint32_t x=0; x<width*subdivision; ++x)
        transport_cells_[static_cast<size_t>(y)*width*subdivision+x] = cells_[(y/subdivision)*width+x/subdivision];
}

std::vector<Cell> Model::transportColumns(const std::vector<Displacement>& velocity, double dt,
                                          std::vector<std::vector<double>>& carried) {
    const uint32_t w=width_*subdivision_, h=height_*subdivision_;
    double max_motion=0;
    for (auto v : velocity) max_motion=std::max({max_motion,std::abs(v.x),std::abs(v.y)});
    const auto steps=static_cast<uint32_t>(std::max(1.0,std::ceil(max_motion*subdivision_/.95)));
    std::vector<Displacement> fine_velocity(transport_cells_.size());
    for (uint32_t y=0; y<h; ++y) for (uint32_t x=0; x<w; ++x) {
        const auto v=velocity[(y/subdivision_)*width_+x/subdivision_];
        fine_velocity[static_cast<size_t>(y)*w+x]={v.x*subdivision_/steps,v.y*subdivision_/steps};
    }
    for (auto& c : transport_cells_) {
        c.continental_age_volume += c.continent*dt;
        c.oceanic_age_volume += c.ocean*dt;
    }
    for (uint32_t step=0; step<steps; ++step) {
        std::vector<Cell> moved_cells(transport_cells_.size());
        for (auto member : {&Cell::continent,&Cell::ocean,&Cell::sediment,&Cell::continental_area,&Cell::oceanic_area,&Cell::flexure_m}) {
            std::vector<double> values(transport_cells_.size()), moved;
            for (size_t i=0; i<values.size(); ++i) values[i]=transport_cells_[i].*member;
            if (member==&Cell::continent || member==&Cell::ocean) {
                const auto moment=member==&Cell::continent?&Cell::continental_age_volume:&Cell::oceanic_age_volume;
                std::vector<double> age(values.size()), moved_age;
                for (size_t i=0; i<values.size(); ++i) age[i]=transport_cells_[i].*moment;
                moved=transport::remap_with_moment(w,h,values,age,fine_velocity,moved_age);
                for (size_t i=0; i<values.size(); ++i) moved_cells[i].*moment=moved_age[i];
                if (member==&Cell::ocean) {
                    for (size_t i=0;i<values.size();++i) age[i]=transport_cells_[i].thermal_rejuvenation_volume;
                    transport::remap_with_moment(w,h,values,age,fine_velocity,moved_age);
                    for (size_t i=0;i<values.size();++i) moved_cells[i].thermal_rejuvenation_volume=moved_age[i];
                }
            } else moved=transport::remap(w,h,values,fine_velocity);
            for (size_t i=0; i<values.size(); ++i) moved_cells[i].*member=moved[i];
        }
        transport_cells_=std::move(moved_cells);
        // Measurement fields on the subcells move by the same remap and sub-steps.
        for (auto& values : carried) values=transport::remap(w,h,values,fine_velocity);
    }
    std::vector<Cell> coarse(cells_.size());
    const double fraction=1.0/(subdivision_*subdivision_);
    for (uint32_t y=0; y<h; ++y) for (uint32_t x=0; x<w; ++x) {
        auto& c=coarse[(y/subdivision_)*width_+x/subdivision_];
        const auto& f=transport_cells_[static_cast<size_t>(y)*w+x];
        for (auto member : {&Cell::continent,&Cell::ocean,&Cell::sediment,&Cell::continental_area,
                &Cell::oceanic_area,&Cell::continental_age_volume,&Cell::oceanic_age_volume,
                &Cell::thermal_rejuvenation_volume,&Cell::flexure_m})
            c.*member += f.*member*fraction;
    }
    return coarse;
}

void Model::reconcileTransport(const std::vector<Cell>& evolved) {
    if (subdivision_==1) { transport_cells_=evolved; return; }
    const uint32_t w=width_*subdivision_;
    const double fraction=1.0/(subdivision_*subdivision_);
    // Coarse physical sources/sinks update persistent subcells conservatively:
    // remove proportionally, add a uniform layer. Existing rock relief is never
    // reconstructed from the averaged map. Fractions and specific ages are
    // mixed at the physics-grid scale; we make no finer-age-resolution claim.
    for (uint32_t y=0; y<height_; ++y) for (uint32_t x=0; x<width_; ++x) {
        const auto& target=evolved[y*width_+x];
        Cell mean;
        for (uint32_t sy=0; sy<subdivision_; ++sy) for (uint32_t sx=0; sx<subdivision_; ++sx) {
            const auto& f=transport_cells_[(y*subdivision_+sy)*w+x*subdivision_+sx];
            for (auto member : {&Cell::continent,&Cell::ocean,&Cell::sediment}) mean.*member += f.*member*fraction;
        }
        for (uint32_t sy=0; sy<subdivision_; ++sy) for (uint32_t sx=0; sx<subdivision_; ++sx) {
            auto& f=transport_cells_[(y*subdivision_+sy)*w+x*subdivision_+sx];
            for (auto member : {&Cell::continent,&Cell::ocean,&Cell::sediment}) {
                const double desired=target.*member, old=mean.*member;
                f.*member = desired < old ? f.*member*(desired/old) : f.*member+(desired-old);
            }
            f.continental_area=target.continental_area; f.oceanic_area=target.oceanic_area;
            f.continental_age_volume=f.continent*(target.continent>0?target.continental_age_volume/target.continent:0);
            f.oceanic_age_volume=f.ocean*(target.ocean>0?target.oceanic_age_volume/target.ocean:0);
            f.thermal_rejuvenation_volume=f.ocean*(target.ocean>0?target.thermal_rejuvenation_volume/target.ocean:0);
            f.flexure_m=target.flexure_m;
        }
    }
}
void spread_roots(const Metric& metric, std::vector<Cell>& cells, double dt, const Isostasy& isostasy,
                  bool sediment_feedback, Budget& budget) {
    // Effective n=3 power-law creep in a lower-crust channel. Unlike globally
    // Newtonian weak crust, small pressure gradients do not rapidly erase quiet
    // terrain. eta=1e20 Pa s at 5 MPa is an explicit rheological assumption;
    // temperature/composition are not resolved. Integrating the cubic shear
    // profile gives 3/5 times the Newtonian flux times (wall stress / 5 MPa)^2.
    // Finite volumes on the cell table: a face's flux per length over its centre
    // distance, times its length / distance over the reference area (mobility per
    // reference area), each side's thickness by its reference / own area.
    constexpr double viscosity = 1e20, reference_stress_pa = 5e6, seconds_per_myr = 3.15576e13;
    constexpr double mobility_factor = continental_density * 9.80665 * seconds_per_myr / (12 * viscosity);
    if (cells.size() != metric.cellcount()) throw std::runtime_error("Root flow cell count mismatch");
    const auto& table = metric.grid().edges();
    const double reference_area = metric.reference_area_m2();
    struct Face { size_t a, b; double mobility; };
    std::vector<Face> faces;
    faces.reserve(table.size());
    std::vector<double> rates(cells.size()), z(cells.size()), channel(cells.size());
    std::vector<double> outgoing(cells.size()), incoming(cells.size()), incoming_age(cells.size());
    std::vector<double> mobile(cells.size()),out_sediment(cells.size()),in_sediment(cells.size());
    double remaining = dt;
    for (int step = 0; remaining > dt * 1e-12; ++step) {
        if (step > 10000) throw std::runtime_error("Root flow requires a smaller outer timestep");
        faces.clear();
        for (auto* buffer : {&rates, &outgoing, &incoming, &incoming_age}) std::fill(buffer->begin(), buffer->end(), 0);
        std::fill(out_sediment.begin(),out_sediment.end(),0);std::fill(in_sediment.begin(),in_sediment.end(),0);
        for (size_t i = 0; i < cells.size(); ++i) {
            z[i] = isostasy.elevation(cells[i]);
            mobile[i]=cells[i].continent+(sediment_feedback?cells[i].sediment*cells[i].continental_area:0);
            channel[i] = std::pow(.5 * mobile[i], 5);
        }
        for (size_t e = 0; e < table.size(); ++e) {
            const size_t a = table[e].a, b = table[e].b;
            const double ka = channel[a], kb = channel[b];
            // Vanishing continental traces can underflow h^5 to zero.
            // Zero conductance means zero flow, never a 0/0 mobility.
            if (ka == 0 || kb == 0) continue;
            const double small = std::min(ka, kb), large = std::max(ka, kb);
            const double harmonic = 2 * small / (1 + small / large);
            const double pressure_gradient = continental_density * 9.80665
                * std::abs(z[a] - z[b]) / metric.distance_m(e);
            const double mobility = mobility_factor * .6 * harmonic
                * std::pow(pressure_gradient / (2 * reference_stress_pa), 2) * metric.weight(e) / reference_area;
            faces.push_back({a, b, mobility});
            // Cubic flux has three times its secant slope. Thickness also
            // changes mobility; a factor two safety margin covers that here.
            rates[a] += 6 * mobility * metric.scale(static_cast<grid::CellId>(a));
            rates[b] += 6 * mobility * metric.scale(static_cast<grid::CellId>(b));
        }
        const double sensitivity = (mantle_density - (sediment_feedback?sediment_density:continental_density)) / (mantle_density - water_density);
        const double max_rate = *std::max_element(rates.begin(), rates.end()) * sensitivity;
        const double sub_dt = max_rate > 0 ? std::min(remaining, .4 / max_rate) : remaining;
        for (const auto& face : faces) {
            const size_t from = z[face.a] > z[face.b] ? face.a : face.b;
            const size_t to = from == face.a ? face.b : face.a;
            const double from_scale = metric.scale(static_cast<grid::CellId>(from));
            const double to_scale = metric.scale(static_cast<grid::CellId>(to));
            // The quarter-volume bound keeps all four outgoing faces positive;
            // it limits flux, never elevation, and preserves the transferred mass.
            const double volume = std::min(.25 * mobile[from] / from_scale,
                face.mobility * std::abs(z[from] - z[to]) * sub_dt);
            const double sediment=sediment_feedback?volume*(cells[from].sediment*cells[from].continental_area/mobile[from]):0;
            const double continent=volume-sediment;
            const double age = cells[from].continent > 0
                ? cells[from].continental_age_volume / cells[from].continent : 0;
            outgoing[from] += continent * from_scale; incoming[to] += continent * to_scale;
            incoming_age[to] += continent * to_scale * age;
            out_sediment[from]+=sediment*from_scale;in_sediment[to]+=sediment*to_scale;
            budget.redistributed_continent_kg += metric.referencemass(continent * continental_density);
            budget.redistributed_sediment_kg += metric.referencemass(sediment * sediment_density);
        }
        for (size_t i = 0; i < cells.size(); ++i) {
            const double retained = cells[i].continent > 0
                ? 1 - std::min(1.0, outgoing[i] / cells[i].continent) : 0;
            cells[i].continent = cells[i].continent * retained + incoming[i];
            cells[i].continental_age_volume = cells[i].continental_age_volume * retained + incoming_age[i];
            if (sediment_feedback) cells[i].sediment=cells[i].sediment-out_sediment[i]+in_sediment[i];
        }
        remaining -= sub_dt;
    }
}
void collapse_orogens(const Metric& metric, std::vector<Cell>& cells, double dt, const margins::Parameters& p,
                      Budget& budget) {
    // Linear diffusion of the cell-mean continent above the maximum thickness, toward the
    // lower side of each face (the cell table's edges): v = kappa max(0, c_a - max(c_b,
    // h_max)). Unlike root flow it reaches interior cells and pure ocean. Continental age
    // moves with the rock; sediment stays. Simultaneous update per sub-step. A positive
    // orogenic_collapse_onset_thickness_m replaces h_max here (not in the area cap): crust
    // is held by its strength below about 60-70 km and spreads only above it.
    // With the onset, two options (both off: the single threshold above):
    // orogenic_collapse_margin_cap: a face whose receiver has oceanic area (more than
    // upper_cap_residual_area) keeps the legacy threshold, the maximum thickness: thick crust
    // spreads into its weak margins and forelands. The receiver alone decides: a face into a
    // cell without oceanic area waits for the onset, whatever the source (a margin feeding
    // an interior cell too). orogenic_collapse_below_onset_km2_myr: weak creep D_low (min(c_a, onset) -
    // max(c_b, cap))+ plus the full D (c_a - max(c_b, onset))+. With both, a margin receiver
    // takes the full collapse above the cap (the creep is subsumed) and interior faces the
    // two regimes. The sub-step count follows the larger coefficient and the reference area
    // (finite volumes: kappa per reference area times the face's length / distance, each
    // side's thickness by its reference / own area).
    if (cells.size() != metric.cellcount()) throw std::runtime_error("Orogenic collapse cell count mismatch");
    const double cap=p.continental_max_thickness_m, onset=p.orogenic_collapse_onset_thickness_m;
    const double h_max=onset>0 ? onset : cap;
    const double diffusivity=p.orogenic_collapse_km2_myr*1e6;
    const double low_diffusivity=onset>0 ? p.orogenic_collapse_below_onset_km2_myr*1e6 : 0;
    const bool margin_cap=onset>0 && p.orogenic_collapse_margin_cap;
    const bool two_regimes=margin_cap || low_diffusivity>0;
    if (diffusivity==0 && low_diffusivity==0) return;
    const double reference_area=metric.reference_area_m2();
    const double total=std::max(diffusivity,low_diffusivity)*dt/reference_area;
    const double steps=std::max(1.0,std::ceil(total/collapse_step_limit));
    if (steps>10000) throw std::runtime_error("Orogenic collapse requires a smaller outer timestep");
    const double kappa=diffusivity*dt/reference_area/steps;
    const double kappa_low=low_diffusivity*dt/reference_area/steps;
    const auto& table=metric.grid().edges();
    std::vector<double> outgoing(cells.size()),incoming(cells.size()),incoming_age(cells.size());
    for (int step=0;step<static_cast<int>(steps);++step) {
        for (auto* buffer : {&outgoing,&incoming,&incoming_age}) std::fill(buffer->begin(),buffer->end(),0);
        for (size_t e=0;e<table.size();++e) {
            const size_t a=table[e].a,b=table[e].b;
            const double ca=cells[a].continent,cb=cells[b].continent;
            size_t from=cells.size();
            double volume=0;
            if (!two_regimes) {
                from=ca>std::max(cb,h_max)?a:cb>std::max(ca,h_max)?b:cells.size();
                if (from==cells.size()) continue;
                volume=kappa*metric.weight(e)*(cells[from].continent-std::max(cells[(from==a?b:a)].continent,h_max));
            } else {
                if (ca==cb) continue;
                from=ca>cb?a:b;
                const double c_from=cells[from].continent,c_to=cells[from==a?b:a].continent;
                if (margin_cap && cells[from==a?b:a].oceanic_area>upper_cap_residual_area)
                    volume=kappa*std::max(0.0,c_from-std::max(c_to,cap));
                else
                    volume=kappa_low*std::max(0.0,std::min(c_from,onset)-std::max(c_to,cap))
                        +kappa*std::max(0.0,c_from-std::max(c_to,onset));
                volume*=metric.weight(e);
                if (!(volume>0)) continue;
            }
            const size_t to=from==a?b:a;
            const double moved_from=volume*metric.scale(static_cast<grid::CellId>(from));
            const double moved_to=volume*metric.scale(static_cast<grid::CellId>(to));
            outgoing[from]+=moved_from;incoming[to]+=moved_to;
            incoming_age[to]+=moved_to*cells[from].continental_age_volume/cells[from].continent;
            budget.collapsed_continent_kg+=metric.referencemass(volume*continental_density);
        }
        for (size_t i=0;i<cells.size();++i) {
            if (outgoing[i]==0 && incoming[i]==0) continue;
            const double retained=cells[i].continent>0?1-outgoing[i]/cells[i].continent:0;
            cells[i].continent=cells[i].continent*retained+incoming[i];
            cells[i].continental_age_volume=cells[i].continental_age_volume*retained+incoming_age[i];
        }
    }
}
double surface_diffusion_number(const Metric& metric, double erosion_strength, double dt) {
    return transport_diffusivity * erosion_strength * dt / metric.reference_area_m2();
}
void diffuse_surface(const Metric& metric, std::vector<Cell>& cells, const std::vector<double>& elevation,
                     double diffusion, double submarine_erosion_ratio, Budget& budget) {
    // Conservative downhill transport. Erode the top sediment first, then
    // continental/oceanic rock; deposit the same mass at the lower neighbour.
    // Density conversion includes sediment expansion; no sediment is deleted.
    // Finite volumes: each face's drop weighted by its length / distance, the source's
    // removal by its reference / own area, each deposit by the source / receiver areas.
    if (cells.size() != metric.cellcount() || elevation.size() != cells.size())
        throw std::runtime_error("Surface diffusion cell count mismatch");
    const auto& grid = metric.grid();
    std::vector<double> deposits(cells.size(), 0);
    for (grid::CellId i = 0; i < cells.size(); ++i) {
        const auto neighbors = grid.neighbours(i);
        double weighted[grid::directioncount];
        double drop_sum = 0;
        for (uint8_t k = 0; k < grid::directioncount; ++k) {
            weighted[k] = metric.weight(metric.edge(i, k)) * std::max(0.0, elevation[i] - elevation[neighbors[k]]);
            drop_sum += weighted[k];
        }
        if (drop_sum == 0 || diffusion == 0) continue;
        auto& c = cells[i];
        // Reduced submarine redistribution: land weathering is not applied at
        // full strength in deep water. Smoothly transition over the upper 1 km
        // to avoid a mobility jump when a cell crosses sea level. This depth and
        // the ratio are coarse modelling choices, not an ocean-current model.
        const double depth_fraction=std::clamp(-elevation[i]/1000,0.0,1.0);
        const double submerged=depth_fraction*depth_fraction*(3-2*depth_fraction);
        const double mobility=1-(1-submarine_erosion_ratio)*submerged;
        const double requested = diffusion * mobility * (drop_sum * metric.scale(i)) * sediment_density;
        double remaining = std::min(requested, mass_per_area(c)), removed = 0;
        const auto erode = [&](double& volume, double density, double* age_volume) {
            const double taken = std::min(volume, std::max(0.0, remaining) / density);
            const double age = age_volume && volume > 0 ? *age_volume / volume : 0;
            volume -= taken;
            if (age_volume) *age_volume = volume * age;
            remaining -= taken * density; removed += taken * density;
            return taken * density;
        };
        erode(c.sediment, sediment_density, nullptr);
        erode(c.continent, continental_density, &c.continental_age_volume);
        const double thermal_age_reduction=c.ocean>0?c.thermal_rejuvenation_volume/c.ocean:0;
        budget.eroded_ocean_kg += erode(c.ocean, oceanic_density, &c.oceanic_age_volume) * metric.area_m2(i);
        c.thermal_rejuvenation_volume=c.ocean*thermal_age_reduction;
        for (uint8_t k = 0; k < grid::directioncount; ++k) {
            const auto j = neighbors[k];
            deposits[j] += removed * weighted[k] / drop_sum / sediment_density * (metric.area_m2(i) / metric.area_m2(j));
        }
        budget.eroded_kg += removed * metric.area_m2(i);
    }
    for (grid::CellId i = 0; i < cells.size(); ++i) {
        cells[i].sediment += deposits[i];
        budget.deposited_kg += deposits[i] * sediment_density * metric.area_m2(i);
    }
}


// Relative tolerance of the implicit kernel's eroded against deposited thickness:
// about 10 times the worst float rounding measured (2.4e-8 in the margins tests,
// 4.9e-9 over the 50 updates of the minimal 128 x 64 tectonics run).
constexpr double implicit_budget_tolerance=2.5e-7;

void Model::routeRiversImplicit(std::vector<Cell>& cells, double dt, double thermal_strength,
                                const margins::Parameters& p, const margins::Forcing* forcing) {
    const size_t count=cells.size();
    const auto metric=landscape::torusmetric(static_cast<int>(width_),static_cast<int>(height_),cell_m_);
    // With a uniform runoff r = 1 m/yr on land the kernel's discharge is
    // Q = A_m2 * 1 m/yr / 1e9 = A_km2 * 1e-3 km^3/yr, so K = river_erosion / sqrt(1e-3)
    // gives K sqrt(Q) = river_erosion sqrt(A_km2): the explicit law (m/Myr per km of
    // sqrt(catchment km^2)) times the slope to the receiver (sea level at an outlet).
    const double erodibility_scale=p.river_erosion/std::sqrt(1e-3);
    std::vector<float> z(count),runoff(count),erodibility(count),uplift(count,0.0f),erosion(count),deposition(count);
    std::vector<uint8_t> sea(count);
    for (size_t i=0;i<count;++i) {
        z[i]=static_cast<float>(elevation_m(cells[i],thermal_strength,root_kg_m3_, root_max_thickness_m_, root_delamination_width_m_));
        sea[i]=z[i]<=0?1:0;
        runoff[i]=sea[i]?0.0f:1.0f;
        const double factor=forcing && !forcing->erosion_multiplier.empty()?forcing->erosion_multiplier[i]:1;
        erodibility[i]=static_cast<float>(erodibility_scale*factor);
    }
    // minimumdepth_m also chooses receivers (a cell raised by less than it takes the
    // steepest descent on the filled surface); 0 routes every flooded cell along the
    // flood order through its spill.
    const auto depressions=landscape::filldepressions(metric,z,sea,0.0f);
    auto graph=landscape::buildflowgraph(metric,z,sea,0.0f,depressions);
    landscape::accumulateflow(metric,runoff,graph);
    landscape::StreamPowerParameters law;
    law.step_myr=dt;
    law.baselevel_m=0;
    law.areaexponent=.5;
    law.slopeexponent=1;
    law.depositioncoefficient=0;
    law.hillslopediffusivity_km2_myr=0; // The engine diffuses separately.
    // A scalar travel length: the kernel has no per-cell settling, so the
    // forcing's sediment_multiplier does not act on implicit rivers.
    law.sedimenttravel_km=p.sediment_travel_km;
    law.deltatop_m=std::numeric_limits<double>::infinity();
    auto surface=z;
    landscape::stepstreampower(metric,graph,erodibility,uplift,law,surface,erosion,deposition);
    // Erosion thickness becomes rock mass at the sediment density, as in the
    // explicit path. A column that cannot supply it gives what it has, and the
    // deposits are scaled so that deposited mass equals removed mass exactly.
    const double area=cell_m_*cell_m_;
    // The kernel's own budget (equal cell areas on the torus): eroded and
    // deposited thickness agree up to the rounding of each cell's value to float.
    double eroded_m=0,deposited_m=0;
    for (size_t i=0;i<count;++i) { eroded_m+=erosion[i];deposited_m+=deposition[i]; }
    if (std::abs(eroded_m-deposited_m)>implicit_budget_tolerance*std::max(eroded_m,deposited_m))
        throw std::runtime_error("Implicit river kernel budget does not close");
    double removed_total=0,requested_total=0;
    for (size_t i=0;i<count;++i) {
        requested_total+=static_cast<double>(deposition[i])*sediment_density;
        if (erosion[i]<=0) continue;
        const double wanted=static_cast<double>(erosion[i])*sediment_density;
        auto& c=cells[i];
        double ocean=0;
        const double removed=removeRock(c,wanted,&ocean);
        budget_.eroded_ocean_kg+=ocean*area;
        if (removed<wanted && c.sediment==0 && c.continent==0 && c.ocean==0) ++budget_.river_unsupplied_cells;
        removed_total+=removed;
        budget_.river_eroded_kg+=removed*area;
        budget_.eroded_kg+=removed*area;
    }
    if (removed_total==0) return;
    if (!(requested_total>0)) throw std::runtime_error("Implicit river erosion without deposition");
    const double scale=removed_total/requested_total;
    for (size_t i=0;i<count;++i) {
        if (deposition[i]<=0) continue;
        const double deposit=static_cast<double>(deposition[i])*sediment_density*scale;
        cells[i].sediment+=deposit/sediment_density;
        budget_.deposited_kg+=deposit*area;
    }
}

void Model::routeSediment(std::vector<Cell>& cells, double dt, double thermal_strength,
                         const margins::Parameters& parameters, const margins::Forcing* forcing) {
    if (!parameters.enabled || (parameters.river_erosion==0 && parameters.wave_erosion_m_myr==0
        && parameters.wave_retreat_km_myr==0)) return;
    // Implicit rivers first; the explicit loop below then carries only wave
    // attack and its load (river term zero), and is skipped without waves.
    margins::Parameters waves_only;
    if (parameters.river_implicit) {
        if (parameters.river_erosion>0) routeRiversImplicit(cells,dt,thermal_strength,parameters,forcing);
        if (parameters.wave_erosion_m_myr==0 && parameters.wave_retreat_km_myr==0) return;
        waves_only=parameters;
        waves_only.river_erosion=0;
    }
    const margins::Parameters& p=parameters.river_implicit?waves_only:parameters;
    const size_t count=cells.size();
    std::vector<double> z(count),catchment(count),load(count),slope(count),distance(count),exposure(count);
    std::vector<double> upstream_drop(count,std::numeric_limits<double>::infinity());
    std::vector<size_t> order(count),downstream(count);
    std::iota(order.begin(),order.end(),0);
    for (size_t i=0;i<count;++i) { z[i]=elevation_m(cells[i],thermal_strength,root_kg_m3_, root_max_thickness_m_, root_delamination_width_m_); downstream[i]=i; }
    for (size_t i=0;i<count;++i) {
        const auto cell=static_cast<grid::CellId>(i);
        catchment[i]=z[i]>0?metric_.area_m2(cell)/1e6:0; // Uniform runoff proxy, km².
        const auto neighbours=metric_.eight(cell);
        int present=0;
        for (const auto& n : neighbours) present+=n.cell!=grid::nocell;
        for (const auto& n : neighbours) {
            if (n.cell==grid::nocell) continue;
            const size_t j=n.cell;
            if (z[i]>0 && z[j]<=0) exposure[i]+=1./present;
            const double length=n.distance_m;
            const double gradient=(z[i]-z[j])/length;
            if (gradient>slope[i]) { slope[i]=gradient;downstream[i]=j;distance[i]=length; }
        }
    }
    std::sort(order.begin(),order.end(),[&](size_t a,size_t b) { return z[a]!=z[b]?z[a]>z[b]:a<b; });
    for (auto i : order) {
        const auto j=downstream[i];
        auto& c=cells[i];
        double removed=0;
        // A river discharges at sea level, not at the adjacent abyssal floor.
        // Ocean routing still uses seabed slopes for sediment transport.
        const double river_slope=p.river_outlet_sea_level && z[i]>0 && j!=i
            ?std::max(0.0,z[i]-std::max(0.0,z[j]))/distance[i]:slope[i];
        const double capacity=p.alluvial_transport*p.river_erosion*std::sqrt(catchment[i])*river_slope*dt*sediment_density;
        if (z[i]>0 && j!=i) {
            // Stream power exports material from catchments without applying a
            // smoothing filter to the ocean floor. Limit incision to a fraction
            // of the downstream drop for explicit-step stability, not a height cap.
            const double factor=forcing && !forcing->erosion_multiplier.empty()?forcing->erosion_multiplier[i]:1;
            // At this resolution beaches/cliffs are subgrid. Exposure is the
            // fraction of neighbouring ocean; exponential attenuation protects
            // high inland terrain. The same routing conserves eroded rock mass.
            // Wave attack at the foot of a high coast can retreat a cliff even
            // when its cell-mean elevation is far above wave reach. Retreat
            // length * exposed relief / cell width gives removed rock thickness.
            const double waves=exposure[i]*(p.wave_erosion_m_myr*std::exp(-z[i]/p.wave_height_scale_m)
                +p.wave_retreat_km_myr*1000*z[i]/metric_.width_m(static_cast<grid::CellId>(i)));
            const double incision=std::min(.25*(z[i]-z[j]),
                factor*(p.river_erosion*std::sqrt(catchment[i])*river_slope+waves)*dt);
            double ocean=0;
            removed=removeRock(c,incision*sediment_density,&ocean);
            budget_.eroded_ocean_kg+=metric_.mass(ocean,static_cast<grid::CellId>(i));
            budget_.river_eroded_kg+=metric_.mass(removed,static_cast<grid::CellId>(i));
            budget_.eroded_kg+=metric_.mass(removed,static_cast<grid::CellId>(i));
        }
        const double sediment=load[i]+removed;
        // Land transfers river load downstream; lakes retain it. Submarine
        // settling deposits a wedge progressively, without eroding quiet ridges.
        const double settling=forcing && !forcing->sediment_multiplier.empty()?forcing->sediment_multiplier[i]:1;
        const double retained=j==i?1:z[i]>0?0:-std::expm1(-distance[i]/(p.sediment_travel_km*1000*settling));
        double deposit=sediment*retained;
        if (p.alluvial_transport>0 && z[i]>0 && j!=i) {
            const double settling_fraction=-std::expm1(-distance[i]/(p.sediment_travel_km*1000*settling));
            deposit=std::max(0.0,sediment-capacity)*settling_fraction;
            // Do not reverse the drainage slope within one explicit update.
            // Excess load continues downstream; no mass or relief is clipped.
            const double room=std::isfinite(upstream_drop[i])?.25*upstream_drop[i]:0;
            deposit=std::min(deposit,room*mantle_density/(mantle_density-sediment_density)*sediment_density);
        }
        c.sediment+=deposit/sediment_density;
        budget_.deposited_kg+=metric_.mass(deposit,static_cast<grid::CellId>(i));
        if (j!=i) {
            // The load per area of the source, on the receiver's area.
            load[j]+=(sediment-deposit)*(metric_.area_m2(static_cast<grid::CellId>(i))/metric_.area_m2(static_cast<grid::CellId>(j)));
            catchment[j]+=catchment[i];
            upstream_drop[j]=std::min(upstream_drop[j],z[i]-z[j]);
        }
    }
}

void Model::advance(const std::vector<Displacement>& velocity, double dt,
                    double erosion_strength, double thermal_strength, double submarine_erosion_ratio,
                    const margins::Forcing* forcing, const margins::Parameters& margin_parameters,
                    const noise_pass::Parameters& noise, bool sediment_feedback,
                    const plumes::Source* plume, double uplift_rate_window_myr) {
    margins::validate(margin_parameters);
    if (plume && (plume->basalt_m.size()!=cells_.size() || plume->influence.size()!=cells_.size()))
        throw std::runtime_error("Mantle plume source shape mismatch");
    noise_pass::validate(noise);
    if (forcing && (forcing->thermal_age_myr.size()!=cells_.size() || forcing->trench_m.size()!=cells_.size()
        || forcing->arc_receiver.size()!=cells_.size() || forcing->owners.size()!=cells_.size()))
        throw std::runtime_error("Margin forcing shape mismatch");
    if (forcing) for (const auto* values:{&forcing->arc_multiplier,&forcing->erosion_multiplier,&forcing->sediment_multiplier}) {
        if (!values->empty() && values->size()!=cells_.size()) throw std::runtime_error("Process variation shape mismatch");
        for (double value:*values) if (!std::isfinite(value) || value<=0)
            throw std::runtime_error("Invalid process variation multiplier");
    }
    if (!budget_.enabled || velocity.size() != cells_.size() || !std::isfinite(dt) || dt <= 0
        || !std::isfinite(erosion_strength) || erosion_strength < 0
        || !std::isfinite(thermal_strength) || thermal_strength < 0
        || !std::isfinite(submarine_erosion_ratio) || submarine_erosion_ratio < 0 || submarine_erosion_ratio > 1
        || !std::isfinite(uplift_rate_window_myr) || uplift_rate_window_myr <= 0)
        throw std::runtime_error("Invalid column update parameters");
    const double diffusion = surface_diffusion_number(metric_, erosion_strength, dt);
    if (diffusion > .2) throw std::runtime_error("Sediment transport requires a smaller timestep");
    for (const auto& v : velocity)
        if (!std::isfinite(v.x) || !std::isfinite(v.y) || std::max(std::abs(v.x), std::abs(v.y)) > .95)
            throw std::runtime_error("Column transport requires displacement <= 0.95 cell");
    const auto index = [&](int x, int y) {
        return static_cast<size_t>((y + static_cast<int>(height_)) % static_cast<int>(height_)) * width_
            + (x + static_cast<int>(width_)) % static_cast<int>(width_);
    };
    // Elevation before the tectonic sub-steps and the running mean of the uplift
    // rate, carried with the rock for the measurement below: on the subcells (the
    // subcell's own elevation, the coarse mean), shifted to positive values because
    // the remap rejects negative ones, with a field of ones whose remap is the
    // arriving rock's coverage. In the first sub-step the positivity limiter cannot
    // act (the shift exceeds 4 times the range, the reconstruction's largest
    // excursion over both axes is 3 times); from the second sub-step on the fields
    // are amounts (value x coverage, coverage 0 to about 2 at plate boundaries), the
    // positivity and shape limiters can act, and moved / coverage is a bounded
    // approximation there.
    std::vector<double> elevation_before(cells_.size());
    for (size_t i = 0; i < cells_.size(); ++i) elevation_before[i] = elevation_m(cells_[i], thermal_strength, root_kg_m3_, root_max_thickness_m_, root_delamination_width_m_);
    const uint32_t fine_width = width_ * subdivision_;
    const auto parent = [&](size_t f) {
        return static_cast<size_t>(f / fine_width / subdivision_) * width_ + (f % fine_width) / subdivision_;
    };
    // Strain of this update on the rock before it moves (strain()): rotate the tensor by
    // the local spin, add sym(grad u); derivatives by one-sided differences within the plate,
    // the neighbour's displacement on the cell's frame over the edge's centre distance
    // (planar: cells per cell, the identity).
    const auto& grid = metric_.grid();
    const auto same_plate = [&](size_t a, size_t b) { return !forcing || forcing->owners[a] == forcing->owners[b]; };
    const auto derivative = [&](grid::CellId i, uint8_t back_k, uint8_t ahead_k, double Displacement::* component) {
        const grid::CellId behind = grid.neighbour(i, back_k), ahead = grid.neighbour(i, ahead_k);
        const bool back = same_plate(i, behind), front = same_plate(i, ahead);
        const double db = (velocity[i].*component - metric_.inframe(i, behind, velocity[behind]).*component)
            / metric_.displacementdistance(metric_.edge(i, back_k));
        const double df = (metric_.inframe(i, ahead, velocity[ahead]).*component - velocity[i].*component)
            / metric_.displacementdistance(metric_.edge(i, ahead_k));
        if (back && front) return std::abs(db) < std::abs(df) ? db : std::abs(df) < std::abs(db) ? df : (db + df) / 2;
        return back ? db : front ? df : 0.0;
    };
    for (grid::CellId i = 0; i < cells_.size(); ++i) {
        const double dux_dx = derivative(i, grid::west, grid::east, &Displacement::x);
        const double duy_dx = derivative(i, grid::west, grid::east, &Displacement::y);
        const double dux_dy = derivative(i, grid::north, grid::south, &Displacement::x);
        const double duy_dy = derivative(i, grid::north, grid::south, &Displacement::y);
        const double spin = (duy_dx - dux_dy) / 2, c = std::cos(spin), s = std::sin(spin);
        const Strain old = strain_[i];
        strain_[i].xx = c * c * old.xx - 2 * c * s * old.xy + s * s * old.yy + dux_dx;
        strain_[i].yy = s * s * old.xx + 2 * c * s * old.xy + c * c * old.yy + duy_dy;
        strain_[i].xy = c * s * (old.xx - old.yy) + (c * c - s * s) * old.xy + (dux_dy + duy_dx) / 2;
    }
    // Carried with the rock: 0 elevation, 1 uplift rate, 2 coverage, 3 to 5 strain xx, yy, xy.
    std::vector<std::vector<double>> carried(6, std::vector<double>(transport_cells_.size(), 1.0));
    for (size_t f = 0; f < transport_cells_.size(); ++f) {
        carried[0][f] = elevation_m(transport_cells_[f], thermal_strength, root_kg_m3_, root_max_thickness_m_, root_delamination_width_m_);
        carried[1][f] = uplift_rate_m_myr_[parent(f)];
        carried[3][f] = strain_[parent(f)].xx;
        carried[4][f] = strain_[parent(f)].yy;
        carried[5][f] = strain_[parent(f)].xy;
    }
    for (int k : {0, 1, 3, 4, 5}) {
        const auto [low, high] = std::minmax_element(carried[k].begin(), carried[k].end());
        const double shift = 4 * (*high - *low) - *low + 1;
        for (double& value : carried[k]) value += shift;
    }
    const auto coarse_mean = [&](const std::vector<double>& fine) {
        std::vector<double> mean(cells_.size(), 0.0);
        const double fraction = 1.0 / (subdivision_ * subdivision_);
        for (size_t f = 0; f < fine.size(); ++f) mean[parent(f)] += fine[f] * fraction;
        return mean;
    };
    const auto elevation_stay = coarse_mean(carried[0]), rate_stay = coarse_mean(carried[1]);
    const auto xx_stay = coarse_mean(carried[3]), yy_stay = coarse_mean(carried[4]), xy_stay = coarse_mean(carried[5]);
    // True cell areas from the cell table (planar: cell_m², the torus' area).
    const auto area_of = [&](size_t i) { return metric_.area_m2(static_cast<grid::CellId>(i)); };
    // The reservoir before this update, for the growth controller (serial sum): felsic,
    // or the continent alone.
    const bool continent_only = margin_parameters.continental_growth_continent_only;
    double reservoir_start = 0;
    if (margin_parameters.continental_growth_control)
        for (size_t i = 0; i < cells_.size(); ++i) {
            const auto& c = cells_[i];
            reservoir_start += (continent_only ? c.continent * continental_density : felsic_mass_per_area(c)) * area_of(i);
        }
    auto next = transportColumns(velocity, dt, carried);
    std::vector<double> arc_source(next.size());
    std::vector<double> created_area(next.size(), 0.0); // New ocean floor of this update, share of the cell.
    const double min_thickness = margin_parameters.continental_min_thickness_m;
    double reservoir_sinks = 0; // This update's subducted sediment and continent, kg.
    double continent_sinks = 0; // This update's subducted continent, kg.
    for (size_t i=0; i<next.size(); ++i) {
        auto& c=next[i];
        const double area = area_of(i);
        double coverage = c.continental_area + c.oceanic_area;
        if (coverage > 1 && c.oceanic_area > 0) {
            // Overlapping oceanic area sinks; continental volume is retained
            // as a thickened root. The mixed ocean slab has a mean age here.
            const double lost_area = std::min(c.oceanic_area, coverage - 1);
            const double fraction = lost_area / c.oceanic_area;
            const double sediment_loss=c.sediment*(lost_area/coverage)*margin_parameters.sediment_subduction_fraction;
            c.sediment-=sediment_loss;
            budget_.subducted_sediment_kg+=sediment_loss*sediment_density*area;
            reservoir_sinks+=sediment_loss*sediment_density*area;
            if (margin_parameters.subduction_erosion_ratio>0 && c.continent>0) {
                // Subduction erosion: the slab drags continental crust down with it.
                const double eroded=std::min(c.continent,margin_parameters.subduction_erosion_ratio*c.ocean*fraction);
                c.continental_age_volume*=(c.continent-eroded)/c.continent;
                c.continent-=eroded;
                budget_.subducted_continent_kg+=eroded*continental_density*area;
                reservoir_sinks+=eroded*continental_density*area;
                continent_sinks+=eroded*continental_density*area;
            }
            budget_.subducted_ocean_kg += c.ocean * fraction * oceanic_density * area;
            const double melt=forcing && !forcing->arc_multiplier.empty()?forcing->arc_multiplier[i]:1;
            arc_source[i]=c.ocean*fraction*margin_parameters.arc_melt_ratio*melt;
            c.ocean *= 1 - fraction;
            c.oceanic_age_volume *= 1 - fraction;
            c.thermal_rejuvenation_volume *= 1-fraction;
            c.oceanic_area -= lost_area;
            coverage -= lost_area;
        }
        if (coverage < 1 && coverage > 0) {
            // Stretch surviving continental columns; only the oceanic share of
            // opening creates new basalt. This avoids cutting instant deep gaps
            // through thick continents. New basalt has age zero.
            double ocean_share=c.oceanic_area/coverage;
            if (margin_parameters.enabled && c.continental_area>0) {
                const double thickness=c.continent/c.continental_area;
                // From the onset down to the continental minimum thickness (0 by default).
                const double onset=margin_parameters.rift_breakup_thickness_m;
                const double thin=min_thickness>0?std::clamp((onset-thickness)/(onset-min_thickness),0.0,1.0)
                    :std::clamp(1-thickness/onset,0.0,1.0);
                ocean_share+=(c.continental_area/coverage)*thin*thin*(3-2*thin);
            }
            const double opening=1-coverage;
            const double generated_area = opening*ocean_share;
            created_area[i] = generated_area;
            // Without margins the division by coverage below grows the continental area by the same share.
            budget_.stretched_continental_area_m2 += (opening-generated_area)*area;
            const double generated = newborn_ocean_crust_[i] * generated_area;
            c.ocean += generated;
            budget_.generated_ocean_kg += generated * oceanic_density * area;
            if (margin_parameters.enabled) {
                // An exhausted continental column cannot remain a 100%
                // continental mantle reference forever. Fill only the newly
                // opened area with ocean floor as stretched crust breaks up.
                // Remaining continental rock is retained; no height is clamped.
                c.oceanic_area+=generated_area;
                c.continental_area+=opening-generated_area;
                coverage=1;
            }
        }
        if (!(coverage > 0)) throw std::runtime_error("Unresolved empty column after remap");
        if (coverage > 1) {
            // Collision: the continental area the overlap squeezes out of the column.
            const double before = c.continental_area;
            budget_.collided_continental_area_m2 += (before - before / coverage) * area;
        }
        c.continental_area /= coverage;
        c.oceanic_area /= coverage;
        c.flexure_m /= coverage;
        if (margin_parameters.enabled && min_thickness>0 && c.continental_area*min_thickness>c.continent) {
            // Continental area beyond continent / minimum thickness (stretched or
            // eroded crust, area kept apart from its volume) relaxes into oceanic
            // area; no rock moves.
            const double freed=(c.continental_area-c.continent/min_thickness)*(-std::expm1(-dt/area_cap_response_myr));
            c.continental_area-=freed;
            c.oceanic_area+=freed;
            budget_.oceanized_continental_area_m2+=freed*area;
        }
        // Advected thermal memory advances at a configurable rate while crust
        // chronology still advances by dt. Newborn basalt starts at zero age.
        // Reducing this rate delays cooling; it does not change initial relief
        // or the final thermal contraction budget, nor depend on trench proximity.
        c.thermal_rejuvenation_volume=std::clamp(c.thermal_rejuvenation_volume
            +(1-margin_parameters.thermal_age_rate)*dt*c.ocean,0.0,c.oceanic_age_volume);
        if (forcing && margin_parameters.enabled) {
            const double response=-std::expm1(-dt/5); // 5 Myr transient; no instantaneous terrain stamp.
            const double target=forcing->thermal_age_myr[i]>=0
                ? forcing->thermal_age_myr[i]*margin_parameters.thermal_age_rate : -1;
            if (target>=0) {
                const double deficit=std::max(0.0,c.oceanic_age_volume-target*c.ocean);
                c.thermal_rejuvenation_volume+=response*std::max(0.0,deficit-c.thermal_rejuvenation_volume);
            }
            const double bending=-std::expm1(-dt/margin_parameters.trench_response_myr);
            c.flexure_m+=bending*(forcing->trench_m[i]-c.flexure_m);
        }
    }
    // Share of arc melt emplaced as continent: fixed, or set by the growth controller
    // to hold the reservoir (felsic, or the continent alone) at M0 (1 + g t): this update's sinks, the target's
    // growth and the reservoir's error over growth_response_myr, over this update's melt.
    double arc_share=margin_parameters.arc_continent_fraction;
    if (margin_parameters.continental_growth_control) {
        double melt=0;
        if (forcing) for (size_t i=0;i<next.size();++i)
            if (arc_source[i]>0 && forcing->arc_receiver[i]!=std::numeric_limits<uint32_t>::max())
                melt+=arc_source[i]*oceanic_density*area_of(i);
        const double initial=continent_only?budget_.initial_continent_kg:budget_.initial_felsic_kg;
        const double growth=margin_parameters.continental_growth_per_gyr;
        const double sinks=continent_only?continent_sinks:reservoir_sinks;
        const double target=initial*(1+growth*time_myr_/1000);
        const double need=sinks+initial*growth*dt/1000+(target-reservoir_start)*dt/growth_response_myr;
        arc_share=melt>0?std::clamp(need/melt,0.0,1.0):0;
        // Book only this update's unmet physical demand; the error term would recount a lasting deficit.
        budget_.continental_growth_shortfall_kg+=std::max(0.0,sinks+initial*growth*dt/1000-melt);
    }
    budget_.arc_continent_fraction=arc_share;
    if (forcing && margin_parameters.enabled) for (size_t i=0;i<next.size();++i) {
        const auto receiver=forcing->arc_receiver[i];
        if (arc_source[i]<=0 || receiver==std::numeric_limits<uint32_t>::max()) continue;
        if (receiver>=next.size()) throw std::runtime_error("Invalid volcanic arc receiver");
        // Resolve a 100 km magma emplacement width over native cell centres.
        // The source is actual removed ocean slab; new magma is separately budgeted.
        std::vector<std::pair<size_t,double>> outlets;
        double total=0;
        if (metric_.planar()) {
            // The torus: a box of ceil(300 km / cell) cells, equal areas.
            const double area=area_of(receiver);
            const int radius=static_cast<int>(std::ceil(300000/cell_m_));
            for (int dy=-radius;dy<=radius;++dy) for (int dx=-radius;dx<=radius;++dx) {
                const auto j=index(static_cast<int>(receiver%width_)+dx,static_cast<int>(receiver/width_)+dy);
                if (forcing->owners[j]!=forcing->owners[receiver]) continue;
                const double weight=std::exp(-.5*(dx*dx+dy*dy)*area/1e10);
                outlets.push_back({j,weight});total+=weight;
            }
        } else {
            // The cells within 300 km (great circle), each weight times its area: the share
            // of the melt's volume, so the emplaced thickness is the Gaussian.
            const auto& mesh=*grid.mesh();
            const mesh::Vec3 centre=mesh.centre(receiver);
            for (const auto j : grid.withinkm(receiver,300)) {
                if (forcing->owners[j]!=forcing->owners[receiver]) continue;
                const double distance=mesh::greatcircle_m(centre,mesh.centre(j),mesh.radius_m());
                const double weight=std::exp(-.5*distance*distance/1e10)*area_of(j);
                outlets.push_back({j,weight});total+=weight;
            }
        }
        for (const auto [j,weight] : outlets) {
            // The source's thickness on the outlet's area (planar: 1).
            const double melt=arc_source[i]*(area_of(i)/area_of(j))*weight/total;
            if (arc_share>0) {
                // The continental share keeps its mass (volume x 2900/2800), age 0. An
                // island arc (mostly oceanic outlet) also gains continental area.
                const double continent=melt*arc_share*oceanic_density/continental_density;
                next[j].ocean+=melt-melt*arc_share;
                next[j].continent+=continent;
                if (next[j].oceanic_area>.5) {
                    const double converted=std::min(next[j].oceanic_area,continent/min_thickness);
                    next[j].oceanic_area-=converted;
                    next[j].continental_area+=converted;
                    budget_.arc_continental_area_m2+=converted*area_of(j);
                }
            } else next[j].ocean+=melt;
        }
        const double mass=arc_source[i]*oceanic_density*area_of(i);
        if (arc_share>0) {
            const double continental_mass=mass*arc_share;
            budget_.generated_ocean_kg+=mass-continental_mass;
            budget_.generated_continent_kg+=continental_mass;
        } else budget_.generated_ocean_kg+=mass;
        budget_.generated_arc_kg+=mass;
    }
    if (plume) for (size_t i=0;i<next.size();++i) {
        // New plume basalt has age zero, like arc magma, and is budgeted as generated ocean crust.
        const double added=plume->basalt_m[i];
        if (!std::isfinite(added) || added<0) throw std::runtime_error("Invalid mantle plume basalt");
        next[i].ocean+=added;
        const double mass=added*oceanic_density*area_of(i);
        budget_.generated_ocean_kg+=mass;budget_.generated_plume_kg+=mass;
        // Reheat ocean floor above the plume toward the swell age, with the ridge
        // heating response; the plate carries the warm, shallow rock away.
        if (plume->swell_age_myr>0 && plume->influence[i]>0 && next[i].ocean>0) {
            const double response=-std::expm1(-dt/5)*plume->influence[i];
            const double deficit=std::max(0.0,next[i].oceanic_age_volume-plume->swell_age_myr*next[i].ocean);
            next[i].thermal_rejuvenation_volume=std::clamp(next[i].thermal_rejuvenation_volume
                +response*std::max(0.0,deficit-next[i].thermal_rejuvenation_volume),0.0,next[i].oceanic_age_volume);
        }
    }
    const Isostasy isostasy{thermal_strength, root_kg_m3_, root_max_thickness_m_, root_delamination_width_m_};
    spread_roots(metric_, next, dt*margin_parameters.root_mobility, isostasy, sediment_feedback, budget_);
    collapse_orogens(metric_, next, dt, margin_parameters, budget_);
    // Upper area cap: oceanic area under continent beyond continental area x the maximum
    // thickness (collisions, arc continent, root inflow, collapse) relaxes into continental
    // area; no rock moves. Before the uplift rate is read, so the rate records it as
    // tectonic: collapse and spreading lower an overthickened margin.
    const double max_thickness=margin_parameters.continental_max_thickness_m;
    if (margin_parameters.enabled && max_thickness>0) for (size_t i=0;i<next.size();++i) {
        auto& c=next[i];
        if (!(c.oceanic_area>0 && c.continent>c.continental_area*max_thickness)) continue;
        double gained=std::min(c.oceanic_area,c.continent/max_thickness-c.continental_area)
            *(-std::expm1(-dt/area_upper_cap_response_myr));
        // A vanishing oceanic area converts fully: its basalt stays in the column as crust the
        // continent overrode, never left under a near-zero area for a rounding overlap to subduct.
        if (c.oceanic_area-gained<upper_cap_residual_area) gained=c.oceanic_area;
        c.continental_area+=gained;
        c.oceanic_area-=gained;
        budget_.continentalized_oceanic_area_m2+=gained*area_of(i);
    }
    // Conservative downhill transport (diffuse_surface) after the uplift rate below.
    std::vector<double> elevation(next.size());
    for (size_t i = 0; i < next.size(); ++i) elevation[i] = elevation_m(next[i], thermal_strength, root_kg_m3_, root_max_thickness_m_, root_delamination_width_m_);
    // Rock uplift rate: the tectonic elevation change of this step (before the
    // surface processes below), kept as an exponential mean. Measurement only.
    // In the plate's frame (semi-Lagrangian: the elevation before the step and the
    // previous mean are those of the rock now in the cell, carried by the engine's
    // own remap on the subcells with the material's sub-steps and displacements), so
    // plate motion over relief is not uplift and the mean follows the rock. The
    // carried value is the arriving rock's mean (divided by its coverage), compared
    // with the value that stayed, so a cell without motion is exact. Exact for a
    // translation when a single column field varies and the elevation is linear in
    // it (land): each column field gets its own limited reconstruction, so with
    // several varying fields the carried elevation differs from the columns' by the
    // limiters' differences. Limits: the remap smears the carried mean by at most
    // a(1 - a) subcells^2 of variance per axis and sub-step (a the fractional
    // subcell displacement; a flattened peak), a(1 - a) n / s^2 cells^2 per step
    // with n sub-steps at subdivision s, and (1 - w) / w = 9.5 times that in steady
    // state at a 10 Myr window and 1 Myr steps (5 cm/yr at 128 columns: a = 0.64,
    // n = 1, s = 4, sigma about 0.37 cells); where plates converge the rock of both
    // sides mixes; at a divergent boundary the cell keeps the rock that stayed (its
    // own side); new crust has no history of its own.
    const double uplift_weight = -std::expm1(-dt / uplift_rate_window_myr);
    const auto coverage = coarse_mean(carried[2]);
    const auto elevation_moved = coarse_mean(carried[0]), rate_moved = coarse_mean(carried[1]);
    for (size_t i = 0; i < next.size(); ++i) {
        const double source_elevation = elevation_before[i] + (elevation_moved[i] / coverage[i] - elevation_stay[i]);
        const double source_rate = uplift_rate_m_myr_[i] + (rate_moved[i] / coverage[i] - rate_stay[i]);
        const double rate = (elevation[i] - source_elevation) / dt;
        uplift_rate_m_myr_[i] = source_rate + (rate - source_rate) * uplift_weight;
    }
    // The strain of the rock now in the cell, carried as the mean above; the ocean floor
    // created in the opening has none, so the cell keeps the arriving rock's tensor on
    // the rest of its area only (rigid plates parting leave the ridge and its floor at 0).
    const auto xx_moved = coarse_mean(carried[3]), yy_moved = coarse_mean(carried[4]), xy_moved = coarse_mean(carried[5]);
    for (size_t i = 0; i < next.size(); ++i) {
        const double kept = 1 - created_area[i];
        strain_[i].xx = (strain_[i].xx + xx_moved[i] / coverage[i] - xx_stay[i]) * kept;
        strain_[i].yy = (strain_[i].yy + yy_moved[i] / coverage[i] - yy_stay[i]) * kept;
        strain_[i].xy = (strain_[i].xy + xy_moved[i] / coverage[i] - xy_stay[i]) * kept;
    }
    diffuse_surface(metric_, next, elevation, diffusion, submarine_erosion_ratio, budget_);
    routeSediment(next,dt,thermal_strength,margin_parameters,forcing);
    if (noise.amplitude_m>0)
        budget_.noise_redistributed_kg+=noise_pass::apply(width_,height_,cell_m_,next,seed_,time_myr_,time_myr_+dt,noise,thermal_strength,root_kg_m3_, root_max_thickness_m_, root_delamination_width_m_,&metric_);
    budget_.mass_kg = 0;
    budget_.felsic_kg = 0;
    budget_.continent_kg = 0;
    for (size_t i=0;i<next.size();++i) {
        validate(next[i]);
        const double area = area_of(i);
        budget_.mass_kg += mass_per_area(next[i]) * area;
        budget_.felsic_kg += felsic_mass_per_area(next[i]) * area;
        budget_.continent_kg += next[i].continent * continental_density * area;
    }
    budget_.relative_mass_residual = (budget_.mass_kg - budget_.initial_mass_kg
        - budget_.generated_ocean_kg - budget_.generated_continent_kg + budget_.subducted_ocean_kg
        + budget_.subducted_sediment_kg + budget_.subducted_continent_kg) / std::max(1.0, budget_.initial_mass_kg);
    // The reservoir's own budget: arc continent and oceanic crust eroded into sediment
    // add to it, subducted sediment and subduction erosion take from it.
    budget_.felsic_relative_residual = (budget_.felsic_kg - budget_.initial_felsic_kg
        - budget_.generated_continent_kg - budget_.eroded_ocean_kg + budget_.subducted_sediment_kg
        + budget_.subducted_continent_kg) / std::max(1.0, budget_.initial_felsic_kg);
    reconcileTransport(next);
    cells_ = std::move(next);
    time_myr_+=dt;
}
} // namespace tectonics::columns
