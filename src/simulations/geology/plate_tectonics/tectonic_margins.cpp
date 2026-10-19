#include "simulations/geology/plate_tectonics/tectonic_margins.hpp"
#include "simulations/geology/plate_tectonics/tectonic_grid.hpp"
#include "simulations/geology/plate_tectonics/tectonic_columns.hpp"
#include "simulations/geology/plate_tectonics/tectonic_variation.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <stdexcept>
#include <tuple>

namespace tectonics::margins {
void validate(const Parameters& p) {
    variation::validate(p.process_variation,p.process_wavelength_km,p.process_period_myr);
    if (!std::isfinite(p.sediment_subduction_fraction) || p.sediment_subduction_fraction<0 || p.sediment_subduction_fraction>1)
        throw std::runtime_error("Invalid sediment subduction fraction");
    for (double value : {p.trench_load_n_m,p.arc_distance_km,p.arc_melt_ratio,p.river_erosion,
            p.thermal_age_rate,p.wave_erosion_m_myr,p.root_mobility,p.alluvial_transport,p.wave_retreat_km_myr})
        if (!std::isfinite(value) || value<0) throw std::runtime_error("Invalid margin process parameter");
    for (double value : {p.flexural_length_km,p.trench_response_myr,p.sediment_travel_km,p.rift_breakup_thickness_m,p.wave_height_scale_m})
        if (!std::isfinite(value) || value<=0) throw std::runtime_error("Invalid margin process length");
    if (p.river_implicit && p.alluvial_transport>0)
        throw std::runtime_error("Implicit rivers have no alluvial transport capacity: set alluvial_transport to 0");
    if (!std::isfinite(p.continental_min_thickness_m) || p.continental_min_thickness_m<0
        || p.continental_min_thickness_m>=p.rift_breakup_thickness_m)
        throw std::runtime_error("Continental minimum thickness must be 0 or below the rift breakup onset");
    if (!std::isfinite(p.continental_max_thickness_m) || p.continental_max_thickness_m<0
        || (p.continental_max_thickness_m>0 && (p.continental_max_thickness_m<=p.rift_breakup_thickness_m
            || p.continental_max_thickness_m<=p.continental_min_thickness_m)))
        throw std::runtime_error("Continental maximum thickness must be 0 or above the rift breakup onset and the minimum thickness");
    if (!std::isfinite(p.orogenic_collapse_km2_myr) || p.orogenic_collapse_km2_myr<0)
        throw std::runtime_error("Invalid orogenic collapse diffusivity");
    if (p.orogenic_collapse_km2_myr>0 && p.continental_max_thickness_m<=0)
        throw std::runtime_error("Orogenic collapse requires a continental maximum thickness");
    if (!std::isfinite(p.orogenic_collapse_onset_thickness_m) || p.orogenic_collapse_onset_thickness_m<0
        || (p.orogenic_collapse_onset_thickness_m>0 && (p.continental_max_thickness_m<=0
            || p.orogenic_collapse_onset_thickness_m<p.continental_max_thickness_m)))
        throw std::runtime_error("Orogenic collapse onset must be 0 or at least the continental maximum thickness");
    if (!std::isfinite(p.orogenic_collapse_below_onset_km2_myr) || p.orogenic_collapse_below_onset_km2_myr<0)
        throw std::runtime_error("Invalid orogenic collapse diffusivity below the onset");
    if ((p.orogenic_collapse_margin_cap || p.orogenic_collapse_below_onset_km2_myr>0) && p.orogenic_collapse_onset_thickness_m<=0)
        throw std::runtime_error("Orogenic collapse margin cap and creep below the onset require an onset thickness");
    if ((p.orogenic_collapse_onset_thickness_m>0 || p.orogenic_collapse_margin_cap || p.orogenic_collapse_below_onset_km2_myr>0)
        && p.orogenic_collapse_km2_myr<=0)
        throw std::runtime_error("Orogenic collapse onset, margin cap and creep below the onset require orogenic collapse");
    if (p.orogenic_collapse_below_onset_km2_myr>p.orogenic_collapse_km2_myr)
        throw std::runtime_error("Orogenic collapse creep below the onset must not exceed the collapse diffusivity");
    if (!std::isfinite(p.arc_continent_fraction) || p.arc_continent_fraction<0 || p.arc_continent_fraction>1)
        throw std::runtime_error("Invalid arc continent fraction");
    if (!std::isfinite(p.subduction_erosion_ratio) || p.subduction_erosion_ratio<0)
        throw std::runtime_error("Invalid subduction erosion ratio");
    if (!std::isfinite(p.continental_growth_per_gyr) || p.continental_growth_per_gyr<-1)
        throw std::runtime_error("Continental growth per Gyr must be finite and at least -1");
    if ((p.arc_continent_fraction>0 || p.continental_growth_control) && p.continental_min_thickness_m<=0)
        throw std::runtime_error("Arc continent and continental growth require a continental minimum thickness");
    if (!p.enabled && (p.continental_min_thickness_m>0 || p.continental_max_thickness_m>0 || p.arc_continent_fraction>0
        || p.subduction_erosion_ratio>0 || p.continental_growth_control))
        throw std::runtime_error("Continental growth processes require active margins");
}
namespace {
bool finite(const mesh::Vec3& v) { return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z); }
} // namespace

Forcing forcing(const columns::Metric& metric,double dt,const std::vector<uint32_t>& owners,
        const std::vector<columns::Cell>& cells,const std::vector<mesh::Vec3>& rigid,const Parameters& p,
        Polarity* polarity,long seed,double time_myr) {
    validate(p);
    const size_t count=metric.cellcount();
    if (count==0 || owners.size()!=count || cells.size()!=count ||
        rigid.size()!=count || !std::isfinite(dt) || dt<=0 || !std::isfinite(time_myr))
        throw std::runtime_error("Invalid margin forcing grid");
    constexpr auto none=std::numeric_limits<uint32_t>::max();
    Forcing result{std::vector<double>(count,-1),std::vector<double>(count),std::vector<uint32_t>(count,none),owners,{},{},{}};
    if (!p.enabled) return result;
    // The process variation at each cell's centre on the sphere.
    const double radius_km=metric.grid().radius_m()/1000;
    const auto factor=[&](size_t i,uint32_t channel) {
        return variation::multiplier(metric.grid().position(static_cast<grid::CellId>(i)),radius_km,time_myr,
            seed,channel,p.process_variation,p.process_wavelength_km,p.process_period_myr);
    };
    if (p.process_variation>0) {
        result.arc_multiplier.resize(count);result.erosion_multiplier.resize(count);result.sediment_multiplier.resize(count);
        for (size_t i=0;i<count;++i) {
            result.arc_multiplier[i]=factor(i,1);result.erosion_multiplier[i]=factor(i,2);
            result.sediment_multiplier[i]=factor(i,3);
        }
    }
    const auto& table=metric.grid().edges();
    struct Source { double rate=0, load=0; uint32_t receiver=std::numeric_limits<uint32_t>::max(); };
    std::vector<Source> ridge(count),trench(count);
    std::vector<Source> overlap(count);
    const auto ocean_age=[&](size_t i) { return cells[i].ocean>0?cells[i].oceanic_age_volume/cells[i].ocean:0; };
    for (size_t i=0;i<count;++i)
        if (!finite(rigid[i]))
            throw std::runtime_error("Invalid margin velocity");
    // The arc walk's reach: the steps continue while the distance walked is below the arc distance.
    const double arc_m=p.arc_distance_km*1000;
    for (size_t e=0;e<table.size();++e) {
        const auto a=table[e].a,b=table[e].b;
        if (owners[a]==owners[b]) continue;
        // Normal relative motion in m/Myr. Sliding alone creates no arc/trench.
        const double closing=metric.closing(e,rigid[a],rigid[b])/dt;
        if (closing<0) {
            for (auto i : {a,b}) if (cells[i].oceanic_area>.5)
                ridge[i].rate=std::max(ridge[i].rate,-closing*.5);
            continue;
        }
        if (closing==0 || (cells[a].continental_area>=.5 && cells[b].continental_area>=.5)) continue;
        // Continental crust overrides oceanic crust; for ocean/ocean contact,
        // older (colder) lithosphere descends. Ties use stable cell ordering.
        auto lower=cells[a].continental_area>=.5?b:cells[b].continental_area>=.5?a:
            (ocean_age(a)>=ocean_age(b)?a:b);
        if (polarity && cells[a].continental_area<.5 && cells[b].continental_area<.5) {
            const std::pair<uint32_t,uint32_t> pair=std::minmax(owners[a],owners[b]);
            const auto [entry,inserted]=polarity->emplace(pair,owners[lower]);
            (void)inserted;
            lower=owners[a]==entry->second?a:b;
        }
        const auto upper=lower==a?b:a;
        uint32_t receiver=upper;
        // Walk inward on the overriding plate, straight on along the face normal (across
        // cube face edges by the cell table's ahead), never across another plate.
        grid::Step at=lower==a?metric.grid().ahead(a,table[e].direction)
                              :grid::Step{a,grid::opposite(table[e].direction)};
        double walked=0;
        while (walked<arc_m) {
            const grid::Step next=metric.grid().ahead(at.cell,at.direction);
            if (owners[next.cell]!=owners[upper]) break;
            walked+=metric.distance_m(metric.edge(at.cell,at.direction));
            receiver=next.cell;
            at=next;
        }
        const double maturity=1-std::exp(-ocean_age(lower)/20);
        const double load=p.trench_load_n_m*(closing/(closing+10000))*maturity*factor(lower,0);
        if (load>trench[lower].load) trench[lower]={closing,load,receiver};
        // Eulerian overlap is often received by the upper boundary cell as
        // the descending slab moves into it. Route that removed basalt to
        // the same arc without incorrectly bending the overriding plate.
        if (load>overlap[upper].load) overlap[upper]={closing,load,receiver};
    }
    // Nearest-source propagation over the cell table's centre distances is confined to its
    // plate, crosses face edges and has deterministic ties. It does not
    // add dozens of raster-face forces at a jagged boundary onto the same cell.
    const auto spread=[&](const std::vector<Source>& sources,bool heating) {
        using Entry=std::tuple<double,uint32_t,uint32_t>;
        std::priority_queue<Entry,std::vector<Entry>,std::greater<Entry>> queue;
        std::vector<double> distance(count,std::numeric_limits<double>::infinity());
        std::vector<uint32_t> origin(count,none);
        for (uint32_t i=0;i<count;++i) if (sources[i].rate>0) {
            distance[i]=metric.width_m(i)*.5; origin[i]=i; queue.emplace(distance[i],i,i);
        }
        const double extent_m=heating?600000:4*p.flexural_length_km*1000;
        while (!queue.empty()) {
            const auto [d,s,i]=queue.top(); queue.pop();
            if (d!=distance[i] || s!=origin[i]) continue;
            for (uint8_t k=0;k<grid::directioncount;++k) {
                const auto j=metric.grid().neighbour(i,k);
                if (owners[j]!=owners[i] || cells[j].oceanic_area<=.5) continue;
                const double next=d+metric.distance_m(metric.edge(i,k));
                if (next>extent_m) continue;
                if (next<distance[j] || (next==distance[j] && s<origin[j])) {
                    distance[j]=next;origin[j]=s;queue.emplace(next,s,j);
                }
            }
        }
        for (size_t i=0;i<count;++i) if (origin[i]!=none) {
            const auto& source=sources[origin[i]];
            if (heating) result.thermal_age_myr[i]=distance[i]/source.rate;
            else {
                const double length=p.flexural_length_km*1000,u=distance[i]/length;
                // Downward lobe of an elastic line-load response. The unresolved
                // outer bulge is omitted. Water-loaded restoring force is rho*g.
                result.trench_m[i]=source.load/(2*(columns::mantle_density-columns::water_density)*9.80665*length)
                    *std::max(0.0,std::exp(-u)*(std::cos(u)+std::sin(u)));
                result.arc_receiver[i]=source.receiver;
            }
        }
    };
    spread(ridge,true); spread(trench,false);
    for (size_t i=0;i<count;++i) if (overlap[i].receiver!=none)
        result.arc_receiver[i]=overlap[i].receiver;
    return result;
}
}
