#include "simulations/geology/plate_tectonics/tectonic_columns.hpp"
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
// rejuvenation scale with the volume left. Returns the mass removed.
double removeRock(Cell& c, double mass) {
    double remaining=mass,removed=0;
    const auto erode=[&](double& volume,double density,double* age) {
        const double taken=std::min(volume,std::max(0.0,remaining)/density);
        if (age) *age*=volume>0?(volume-taken)/volume:0;
        volume-=taken;remaining-=taken*density;removed+=taken*density;
    };
    erode(c.sediment,sediment_density,nullptr);
    erode(c.continent,continental_density,&c.continental_age_volume);
    const double deficit=c.ocean>0?c.thermal_rejuvenation_volume/c.ocean:0;
    erode(c.ocean,oceanic_density,&c.oceanic_age_volume);
    c.thermal_rejuvenation_volume=c.ocean*deficit;
    return removed;
}
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
double age_myr(const Cell& c) {
    return (c.continent + c.ocean) > 0
        ? (c.continental_age_volume + c.oceanic_age_volume) / (c.continent + c.ocean) : 0;
}
double mountain_force_n_m(const Cell& c, bool sediment_feedback) {
    // Airy-compensated continental column, relative to unthickened crust:
    // delta GPE = g*rho_c*(1-rho_c/rho_m)*(H²-H0²)/2, in N/m.
    // Effective mixed-cell column: no division by vanishing area fractions.
    // This reduced crust-only estimate omits thermal mantle GPE and flexure.
    const double reference = reference_continent_m * c.continental_area;
    if (sediment_feedback) {
        // Exact two-layer Airy GPE for continental basement plus its sediment
        // cover, relative to the same unthickened reference column. Units N/m.
        // Mixed-cell sediment participation follows the continental area.
        const double s=c.sediment*c.continental_area;
        const double z=((mantle_density-continental_density)*(c.continent-reference)
            +(mantle_density-sediment_density)*s)/mantle_density;
        const double base=z-c.continent-s,interface_z=z-s;
        const double gpe=sediment_density*z*z+(continental_density-sediment_density)*interface_z*interface_z
            +(mantle_density-continental_density)*(base*base-reference*reference);
        return .5*9.80665*std::max(0.0,gpe);
    }
    return .5 * 9.80665 * continental_density * (1 - continental_density / mantle_density)
        * std::max(0.0, c.continent * c.continent - reference * reference);
}
double elevation_m(const Cell& c, double thermal_strength) {
    const double ocean_age = c.ocean > 0 ? std::max(0.0,
        (c.oceanic_age_volume-c.thermal_rejuvenation_volume)/c.ocean) : 0;
    const double load = ocean_reference + continental_mantle_load * c.continental_area
        + (mantle_density - water_density) * thermal_strength
            * thermal::subsidence_m(ocean_age) * c.oceanic_area;
    const double buoyancy = (mantle_density - continental_density) * c.continent
        + (mantle_density - oceanic_density) * c.ocean
        + (mantle_density - sediment_density) * c.sediment - load;
    // Hydrostatic compensation includes water loading below sea level.
    return buoyancy / (buoyancy >= 0 ? mantle_density : mantle_density - water_density) - c.flexure_m;
}
Cell seed(double z, bool continental, double age) {
    if (!std::isfinite(z) || !std::isfinite(age) || age < 0)
        throw std::runtime_error("Invalid initial column elevation or age");
    Cell c;
    if (continental) {
        c.continental_area = 1;
        c.continent = reference_continent_m
            + z * (z >= 0 ? mantle_density : mantle_density - water_density)
                / (mantle_density - continental_density);
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
        uint32_t subdivision, long seed, const ocean_relief::Parameters& relief, bool preserve_initial) {
    if (width < 2 || height < 2 || cells.size() != static_cast<size_t>(width) * height
        || !std::isfinite(cell_m) || cell_m <= 0 || (subdivision != 1 && subdivision != 2 && subdivision != 4))
        throw std::runtime_error("Invalid column grid");
    width_ = width; height_ = height; cell_m_ = cell_m; cells_ = std::move(cells);
    seed_=seed;time_myr_=0;
    newborn_ocean_crust_ = ocean_relief::crust(width,height,cell_m,seed,relief);
    for (size_t i=0; i<cells_.size(); ++i) {
        auto& c=cells_[i];
        const double ratio=preserve_initial?1:newborn_ocean_crust_[i]/ridge_crust_m;
        c.ocean*=ratio; c.oceanic_age_volume*=ratio; c.thermal_rejuvenation_volume*=ratio;
    }
    budget_ = {}; budget_.enabled = true;
    uplift_rate_m_myr_.assign(cells_.size(), 0.0);
    strain_.assign(cells_.size(), Strain{});
    for (const auto& c : cells_) { validate(c); budget_.initial_mass_kg += mass_per_area(c) * cell_m * cell_m; }
    budget_.mass_kg = budget_.initial_mass_kg;
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
void Model::spreadRoots(std::vector<Cell>& cells, double dt, double thermal_strength, bool sediment_feedback) {
    // Effective n=3 power-law creep in a lower-crust channel. Unlike globally
    // Newtonian weak crust, small pressure gradients do not rapidly erase quiet
    // terrain. eta=1e20 Pa s at 5 MPa is an explicit rheological assumption;
    // temperature/composition are not resolved. Integrating the cubic shear
    // profile gives 3/5 times the Newtonian flux times (wall stress / 5 MPa)^2.
    constexpr double viscosity = 1e20, reference_stress_pa = 5e6, seconds_per_myr = 3.15576e13;
    constexpr double mobility_factor = continental_density * 9.80665 * seconds_per_myr / (12 * viscosity);
    const auto index = [&](uint32_t x, uint32_t y) { return static_cast<size_t>(y) * width_ + x; };
    struct Edge { size_t a, b; double mobility; };
    std::vector<Edge> edges;
    edges.reserve(2 * cells.size());
    std::vector<double> rates(cells.size()), z(cells.size()), channel(cells.size());
    std::vector<double> outgoing(cells.size()), incoming(cells.size()), incoming_age(cells.size());
    std::vector<double> mobile(cells.size()),out_sediment(cells.size()),in_sediment(cells.size());
    double remaining = dt;
    for (int step = 0; remaining > dt * 1e-12; ++step) {
        if (step > 10000) throw std::runtime_error("Root flow requires a smaller outer timestep");
        edges.clear();
        for (auto* buffer : {&rates, &outgoing, &incoming, &incoming_age}) std::fill(buffer->begin(), buffer->end(), 0);
        std::fill(out_sediment.begin(),out_sediment.end(),0);std::fill(in_sediment.begin(),in_sediment.end(),0);
        for (size_t i = 0; i < cells.size(); ++i) {
            z[i] = elevation_m(cells[i], thermal_strength);
            mobile[i]=cells[i].continent+(sediment_feedback?cells[i].sediment*cells[i].continental_area:0);
            channel[i] = std::pow(.5 * mobile[i], 5);
        }
        for (uint32_t y = 0; y < height_; ++y) for (uint32_t x = 0; x < width_; ++x) {
            const auto a = index(x, y);
            for (auto b : {index((x + 1) % width_, y), index(x, (y + 1) % height_)}) {
                const double ka = channel[a], kb = channel[b];
                // Vanishing continental traces can underflow h^5 to zero.
                // Zero conductance means zero flow, never a 0/0 mobility.
                if (ka == 0 || kb == 0) continue;
                const double small = std::min(ka, kb), large = std::max(ka, kb);
                const double harmonic = 2 * small / (1 + small / large);
                const double pressure_gradient = continental_density * 9.80665
                    * std::abs(z[a] - z[b]) / cell_m_;
                const double mobility = mobility_factor * .6 * harmonic
                    * std::pow(pressure_gradient / (2 * reference_stress_pa), 2) / (cell_m_ * cell_m_);
                edges.push_back({a, b, mobility});
                // Cubic flux has three times its secant slope. Thickness also
                // changes mobility; a factor two safety margin covers that here.
                rates[a] += 6 * mobility; rates[b] += 6 * mobility;
            }
        }
        const double sensitivity = (mantle_density - (sediment_feedback?sediment_density:continental_density)) / (mantle_density - water_density);
        const double max_rate = *std::max_element(rates.begin(), rates.end()) * sensitivity;
        const double sub_dt = max_rate > 0 ? std::min(remaining, .4 / max_rate) : remaining;
        for (const auto& edge : edges) {
            const size_t from = z[edge.a] > z[edge.b] ? edge.a : edge.b;
            const size_t to = from == edge.a ? edge.b : edge.a;
            // The quarter-volume bound keeps all four outgoing faces positive;
            // it limits flux, never elevation, and preserves the transferred mass.
            const double volume = std::min(.25 * mobile[from],
                edge.mobility * std::abs(z[from] - z[to]) * sub_dt);
            const double sediment=sediment_feedback?volume*(cells[from].sediment*cells[from].continental_area/mobile[from]):0;
            const double continent=volume-sediment;
            const double age = cells[from].continent > 0
                ? cells[from].continental_age_volume / cells[from].continent : 0;
            outgoing[from] += continent; incoming[to] += continent;
            incoming_age[to] += continent * age;
            out_sediment[from]+=sediment;in_sediment[to]+=sediment;
            budget_.redistributed_continent_kg += continent * continental_density * cell_m_ * cell_m_;
            budget_.redistributed_sediment_kg += sediment * sediment_density * cell_m_ * cell_m_;
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
        z[i]=static_cast<float>(elevation_m(cells[i],thermal_strength));
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
        const double removed=removeRock(c,wanted);
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
    for (size_t i=0;i<count;++i) { z[i]=elevation_m(cells[i],thermal_strength); downstream[i]=i; }
    const auto index=[&](int x,int y) {
        return static_cast<size_t>((y+static_cast<int>(height_))%static_cast<int>(height_))*width_
            +(x+static_cast<int>(width_))%static_cast<int>(width_);
    };
    for (size_t i=0;i<count;++i) {
        const int x=static_cast<int>(i%width_),y=static_cast<int>(i/width_);
        catchment[i]=z[i]>0?cell_m_*cell_m_/1e6:0; // Uniform runoff proxy, km².
        for (int dy=-1;dy<=1;++dy) for (int dx=-1;dx<=1;++dx) {
            if (!dx && !dy) continue;
            const auto j=index(x+dx,y+dy);
            if (z[i]>0 && z[j]<=0) exposure[i]+=1./8;
            const double length=cell_m_*(dx && dy?std::sqrt(2.):1.);
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
                +p.wave_retreat_km_myr*1000*z[i]/cell_m_);
            const double incision=std::min(.25*(z[i]-z[j]),
                factor*(p.river_erosion*std::sqrt(catchment[i])*river_slope+waves)*dt);
            removed=removeRock(c,incision*sediment_density);
            budget_.river_eroded_kg+=removed*cell_m_*cell_m_;
            budget_.eroded_kg+=removed*cell_m_*cell_m_;
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
        budget_.deposited_kg+=deposit*cell_m_*cell_m_;
        if (j!=i) {
            load[j]+=sediment-deposit;catchment[j]+=catchment[i];
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
    const double diffusion = transport_diffusivity * erosion_strength * dt / (cell_m_ * cell_m_);
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
    for (size_t i = 0; i < cells_.size(); ++i) elevation_before[i] = elevation_m(cells_[i], thermal_strength);
    const uint32_t fine_width = width_ * subdivision_;
    const auto parent = [&](size_t f) {
        return static_cast<size_t>(f / fine_width / subdivision_) * width_ + (f % fine_width) / subdivision_;
    };
    // Strain of this update on the rock before it moves (strain()): rotate the tensor by
    // the local spin, add sym(grad u); derivatives by one-sided differences within the plate.
    const auto same_plate = [&](size_t a, size_t b) { return !forcing || forcing->owners[a] == forcing->owners[b]; };
    const auto derivative = [&](size_t i, size_t behind, size_t ahead, double Displacement::* component) {
        const bool back = same_plate(i, behind), front = same_plate(i, ahead);
        const double db = velocity[i].*component - velocity[behind].*component;
        const double df = velocity[ahead].*component - velocity[i].*component;
        if (back && front) return std::abs(db) < std::abs(df) ? db : std::abs(df) < std::abs(db) ? df : (db + df) / 2;
        return back ? db : front ? df : 0.0;
    };
    for (uint32_t y = 0; y < height_; ++y) for (uint32_t x = 0; x < width_; ++x) {
        const int ix = static_cast<int>(x), iy = static_cast<int>(y);
        const auto i = index(ix, iy), west = index(ix - 1, iy), east = index(ix + 1, iy);
        const auto north = index(ix, iy - 1), south = index(ix, iy + 1);
        const double dux_dx = derivative(i, west, east, &Displacement::x);
        const double duy_dx = derivative(i, west, east, &Displacement::y);
        const double dux_dy = derivative(i, north, south, &Displacement::x);
        const double duy_dy = derivative(i, north, south, &Displacement::y);
        const double spin = (duy_dx - dux_dy) / 2, c = std::cos(spin), s = std::sin(spin);
        const Strain old = strain_[i];
        strain_[i].xx = c * c * old.xx - 2 * c * s * old.xy + s * s * old.yy + dux_dx;
        strain_[i].yy = s * s * old.xx + 2 * c * s * old.xy + c * c * old.yy + duy_dy;
        strain_[i].xy = c * s * (old.xx - old.yy) + (c * c - s * s) * old.xy + (dux_dy + duy_dx) / 2;
    }
    // Carried with the rock: 0 elevation, 1 uplift rate, 2 coverage, 3 to 5 strain xx, yy, xy.
    std::vector<std::vector<double>> carried(6, std::vector<double>(transport_cells_.size(), 1.0));
    for (size_t f = 0; f < transport_cells_.size(); ++f) {
        carried[0][f] = elevation_m(transport_cells_[f], thermal_strength);
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
    auto next = transportColumns(velocity, dt, carried);
    const double area = cell_m_ * cell_m_;
    std::vector<double> arc_source(next.size());
    std::vector<double> created_area(next.size(), 0.0); // New ocean floor of this update, share of the cell.
    for (size_t i=0; i<next.size(); ++i) {
        auto& c=next[i];
        double coverage = c.continental_area + c.oceanic_area;
        if (coverage > 1 && c.oceanic_area > 0) {
            // Overlapping oceanic area sinks; continental volume is retained
            // as a thickened root. The mixed ocean slab has a mean age here.
            const double lost_area = std::min(c.oceanic_area, coverage - 1);
            const double fraction = lost_area / c.oceanic_area;
            const double sediment_loss=c.sediment*(lost_area/coverage)*margin_parameters.sediment_subduction_fraction;
            c.sediment-=sediment_loss;
            budget_.subducted_sediment_kg+=sediment_loss*sediment_density*area;
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
                const double thin=std::clamp(1-thickness/margin_parameters.rift_breakup_thickness_m,0.0,1.0);
                ocean_share+=(c.continental_area/coverage)*thin*thin*(3-2*thin);
            }
            const double opening=1-coverage;
            const double generated_area = opening*ocean_share;
            created_area[i] = generated_area;
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
        c.continental_area /= coverage;
        c.oceanic_area /= coverage;
        c.flexure_m /= coverage;
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
    if (forcing && margin_parameters.enabled) for (size_t i=0;i<next.size();++i) {
        const auto receiver=forcing->arc_receiver[i];
        if (arc_source[i]<=0 || receiver==std::numeric_limits<uint32_t>::max()) continue;
        if (receiver>=next.size()) throw std::runtime_error("Invalid volcanic arc receiver");
        // Resolve a 100 km magma emplacement width over native cell centres.
        // The source is actual removed ocean slab; new magma is separately budgeted.
        const int radius=static_cast<int>(std::ceil(300000/cell_m_));
        std::vector<std::pair<size_t,double>> outlets;
        double total=0;
        for (int dy=-radius;dy<=radius;++dy) for (int dx=-radius;dx<=radius;++dx) {
            const auto j=index(static_cast<int>(receiver%width_)+dx,static_cast<int>(receiver/width_)+dy);
            if (forcing->owners[j]!=forcing->owners[receiver]) continue;
            const double weight=std::exp(-.5*(dx*dx+dy*dy)*area/1e10);
            outlets.push_back({j,weight});total+=weight;
        }
        for (const auto [j,weight] : outlets) next[j].ocean+=arc_source[i]*weight/total;
        const double mass=arc_source[i]*oceanic_density*area;
        budget_.generated_ocean_kg+=mass;budget_.generated_arc_kg+=mass;
    }
    if (plume) for (size_t i=0;i<next.size();++i) {
        // New plume basalt has age zero, like arc magma, and is budgeted as generated ocean crust.
        const double added=plume->basalt_m[i];
        if (!std::isfinite(added) || added<0) throw std::runtime_error("Invalid mantle plume basalt");
        next[i].ocean+=added;
        const double mass=added*oceanic_density*area;
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
    spreadRoots(next, dt*margin_parameters.root_mobility, thermal_strength,sediment_feedback);
    // Conservative downhill transport. Erode the top sediment first, then
    // continental/oceanic rock; deposit the same mass at the lower neighbour.
    // Density conversion includes sediment expansion; no sediment is deleted.
    std::vector<double> elevation(next.size()), deposits(next.size(), 0);
    for (size_t i = 0; i < next.size(); ++i) elevation[i] = elevation_m(next[i], thermal_strength);
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
    for (uint32_t y = 0; y < height_; ++y) for (uint32_t x = 0; x < width_; ++x) {
        const auto i = index(x, y);
        const size_t neighbors[] = {index(static_cast<int>(x) - 1, y), index(x + 1, y),
            index(x, static_cast<int>(y) - 1), index(x, y + 1)};
        double drop_sum = 0;
        for (auto j : neighbors) drop_sum += std::max(0.0, elevation[i] - elevation[j]);
        if (drop_sum == 0 || diffusion == 0) continue;
        auto& c = next[i];
        // Reduced submarine redistribution: land weathering is not applied at
        // full strength in deep water. Smoothly transition over the upper 1 km
        // to avoid a mobility jump when a cell crosses sea level. This depth and
        // the ratio are coarse modelling choices, not an ocean-current model.
        const double depth_fraction=std::clamp(-elevation[i]/1000,0.0,1.0);
        const double submerged=depth_fraction*depth_fraction*(3-2*depth_fraction);
        const double mobility=1-(1-submarine_erosion_ratio)*submerged;
        const double requested = diffusion * mobility * drop_sum * sediment_density;
        double remaining = std::min(requested, mass_per_area(c)), removed = 0;
        const auto erode = [&](double& volume, double density, double* age_volume) {
            const double taken = std::min(volume, std::max(0.0, remaining) / density);
            const double age = age_volume && volume > 0 ? *age_volume / volume : 0;
            volume -= taken;
            if (age_volume) *age_volume = volume * age;
            remaining -= taken * density; removed += taken * density;
        };
        erode(c.sediment, sediment_density, nullptr);
        erode(c.continent, continental_density, &c.continental_age_volume);
        const double thermal_age_reduction=c.ocean>0?c.thermal_rejuvenation_volume/c.ocean:0;
        erode(c.ocean, oceanic_density, &c.oceanic_age_volume);
        c.thermal_rejuvenation_volume=c.ocean*thermal_age_reduction;
        for (auto j : neighbors)
            deposits[j] += removed * std::max(0.0, elevation[i] - elevation[j]) / drop_sum / sediment_density;
        budget_.eroded_kg += removed * area;
    }
    for (size_t i = 0; i < next.size(); ++i) {
        next[i].sediment += deposits[i];
        budget_.deposited_kg += deposits[i] * sediment_density * area;
    }
    routeSediment(next,dt,thermal_strength,margin_parameters,forcing);
    if (noise.amplitude_m>0)
        budget_.noise_redistributed_kg+=noise_pass::apply(width_,height_,cell_m_,next,seed_,time_myr_,time_myr_+dt,noise,thermal_strength);
    budget_.mass_kg = 0;
    for (size_t i=0;i<next.size();++i) {
        validate(next[i]);
        budget_.mass_kg += mass_per_area(next[i]) * area;
    }
    budget_.relative_mass_residual = (budget_.mass_kg - budget_.initial_mass_kg
        - budget_.generated_ocean_kg + budget_.subducted_ocean_kg + budget_.subducted_sediment_kg) / std::max(1.0, budget_.initial_mass_kg);
    reconcileTransport(next);
    cells_ = std::move(next);
    time_myr_+=dt;
}
} // namespace tectonics::columns
