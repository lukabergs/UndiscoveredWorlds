#include "simulations/geology/plate_tectonics/tectonic_noise_pass.hpp"
#include "simulations/geology/plate_tectonics/tectonic_columns.hpp"
#include "simulations/geology/plate_tectonics/tectonic_ocean_relief.hpp"
#include <algorithm>
#include <cmath>
#include <optional>
#include <stdexcept>

namespace tectonics::noise_pass {
void validate(const Parameters& p) {
    if (!std::isfinite(p.amplitude_m) || p.amplitude_m<0 ||
        !std::isfinite(p.wavelength_km) || p.wavelength_km<=0 ||
        !std::isfinite(p.period_myr) || p.period_myr<=0 || !std::isfinite(p.coast_taper_m) || p.coast_taper_m<0)
        throw std::runtime_error("Invalid experimental noise pass");
}
double apply(uint32_t w,uint32_t h,double cell_m,std::vector<columns::Cell>& cells,
             long seed,double start,double end,const Parameters& p,double thermal_strength,double root_kg_m3,
             double root_max_thickness_m, double root_delamination_width_m, const columns::Metric* metric) {
    using namespace columns;
    validate(p);
    const bool planar=!metric || metric->planar();
    if ((planar && (w<2 || h<2 || cells.size()!=static_cast<size_t>(w)*h ||
        !std::isfinite(cell_m) || cell_m<=0)) || (!planar && cells.size()!=metric->cellcount()) ||
        !std::isfinite(start) || !std::isfinite(end) || start<0 || end<start ||
        !std::isfinite(thermal_strength) || thermal_strength<0 ||
        !std::isfinite(root_kg_m3) || root_kg_m3<0 || root_kg_m3>=mantle_density-continental_density ||
        !valid_root_max_thickness(root_max_thickness_m) ||
        !valid_root_delamination_width(root_max_thickness_m,root_delamination_width_m))
        throw std::runtime_error("Invalid noise pass grid or time");
    if (p.amplitude_m==0 || start==end) return 0;
    if (!planar) throw std::logic_error("noise pass: the noise is not on the cube yet (E4)");
    std::optional<Metric> own;
    if (!metric) metric=&own.emplace(nullptr,w,h,cell_m);
    constexpr double tau=6.28318530717958647692;
    const auto kx=static_cast<uint32_t>(std::clamp(std::round(w*cell_m/(p.wavelength_km*1000)),1.,static_cast<double>(w)));
    const auto ky=static_cast<uint32_t>(std::clamp(std::round(h*cell_m/(p.wavelength_km*1000)),1.,static_cast<double>(h)));
    const auto ramp=[](double t) { return -std::expm1(-t/30.); };
    const double steady=.6*(ramp(end)-ramp(start));
    const double varying=.4*(ramp(end)*std::sin(tau*end/p.period_myr)
        -ramp(start)*std::sin(tau*start/p.period_myr));
    std::vector<double> potential(cells.size());
    std::vector<double> coast_weight(cells.size(),1);
    std::vector<bool> land(cells.size());
    for (uint32_t y=0;y<h;++y) for (uint32_t x=0;x<w;++x) {
        const double u=(x+.5)/w,v=(y+.5)/h;
        potential[y*w+x]=p.amplitude_m*(steady*ocean_relief::noise(u,v,kx,ky,seed)
            +varying*ocean_relief::noise(u,v,kx,ky,static_cast<uint32_t>(seed)^0x9e3779b9U));
        if (p.coast_taper_m>0) {
            const double z=elevation_m(cells[y*w+x],thermal_strength,root_kg_m3,root_max_thickness_m,root_delamination_width_m);
            const double t=std::min(1.,std::abs(z)/p.coast_taper_m);
            coast_weight[y*w+x]=t*t*(3-2*t);land[y*w+x]=z>=0;
        }
    }
    // Face transfers over the cell table: each face's potential difference times its
    // length / distance, each side's thickness by its reference / own area (planar: 1).
    const auto& table=metric->grid().edges();
    double transferred=0;
    for (const auto member:{&Cell::continent,&Cell::ocean}) {
        const bool continental=member==&Cell::continent;
        const double density=continental?continental_density:oceanic_density;
        const auto moment=continental?&Cell::continental_age_volume:&Cell::oceanic_age_volume;
        std::vector<double> volume(cells.size()),age(cells.size()),thermal(cells.size());
        for (size_t e=0;e<table.size();++e) {
            const size_t a=table[e].a,b=table[e].b;
            // Do not move a rock component across its support boundary.
            if (cells[a].*member<=0 || cells[b].*member<=0) continue;
            if (p.coast_taper_m>0 && land[a]!=land[b]) continue;
            // The potential's difference, with the cube's non-orthogonal correction.
            const double rise=metric->planar()?potential[b]-potential[a]:-metric->difference(e,potential);
            const double q=.125*rise*mantle_density/(mantle_density-density)
                *std::min(coast_weight[a],coast_weight[b])*metric->weight(e);
            const auto from=q>0?a:b,to=q>0?b:a;
            const double from_scale=metric->scale(static_cast<grid::CellId>(from));
            const double to_scale=metric->scale(static_cast<grid::CellId>(to));
            const double original=cells[from].*member;
            // Four incident faces: positivity bound, never a height clamp.
            const double moved=std::min(.24*original/from_scale,std::abs(q));
            volume[from]-=moved*from_scale;volume[to]+=moved*to_scale;
            const double moved_age=moved*(cells[from].*moment)/original;
            age[from]-=moved_age*from_scale;age[to]+=moved_age*to_scale;
            if (!continental) {
                const double deficit=moved*cells[from].thermal_rejuvenation_volume/original;
                thermal[from]-=deficit*from_scale;thermal[to]+=deficit*to_scale;
            }
            transferred+=metric->referencemass(moved*density);
        }
        for (size_t i=0;i<cells.size();++i) {
            cells[i].*member+=volume[i];cells[i].*moment+=age[i];
            if (!continental) cells[i].thermal_rejuvenation_volume+=thermal[i];
        }
    }
    return transferred;
}
}
