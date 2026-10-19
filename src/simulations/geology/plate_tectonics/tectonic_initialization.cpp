#include "simulations/geology/plate_tectonics/tectonic_initialization.hpp"
#include "simulations/geology/plate_tectonics/tectonic_grid.hpp"
#include "simulations/geology/plate_tectonics/tectonic_thermal.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <stdexcept>
#include <tuple>
#include <type_traits>

namespace tectonics::initialization {
namespace {
bool finite(const transport::Displacement& d) { return std::isfinite(d.x) && std::isfinite(d.y); }
bool finite(const mesh::Vec3& v) { return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z); }

template <class V>
std::vector<columns::Cell> terrainon(uint32_t w,uint32_t h,double cell_m,double dt,
        const std::vector<double>& elevations,const std::vector<double>& ages,
        const std::vector<uint32_t>& owners,const std::vector<V>& velocity,
        double thermal_strength, double margin_width_m, double boundary_influence, double root_kg_m3,
        double root_max_thickness_m, double root_delamination_width_m, const grid::TectonicGrid* grid) {
    using namespace columns;
    const bool planar=!grid || grid->kind()==grid::TectonicGrid::Kind::planar;
    const size_t count=planar?static_cast<size_t>(w)*h:grid->cellcount();
    if ((planar && (w<2 || h<2 || !std::isfinite(cell_m) || cell_m<=0)) || elevations.size()!=count ||
        ages.size()!=count || owners.size()!=count || velocity.size()!=count || !std::isfinite(dt) || dt<=0 ||
        !std::isfinite(thermal_strength) || thermal_strength<0 ||
        !std::isfinite(margin_width_m) || margin_width_m<0 || margin_width_m>4000 ||
        !std::isfinite(boundary_influence) || boundary_influence<0 || boundary_influence>1 ||
        !std::isfinite(root_kg_m3) || root_kg_m3<0 || root_kg_m3>=mantle_density-continental_density ||
        !valid_root_max_thickness(root_max_thickness_m) ||
        !valid_root_delamination_width(root_max_thickness_m,root_delamination_width_m))
        throw std::runtime_error("Invalid initial terrain grid");
    for (size_t i=0;i<count;++i) if (!std::isfinite(elevations[i]) || !std::isfinite(ages[i]) || ages[i]<0 ||
        !finite(velocity[i]))
        throw std::runtime_error("Invalid initial terrain cell");
    // The cell table (W2-E): openings on its edges through the displacement adapter, the
    // spreading distance by Dijkstra in metres.
    const Metric metric(grid,w,h,cell_m);
    metric.requiredisplacement(std::is_same_v<V,mesh::Vec3>,"initialization");
    const auto& table=metric.grid().edges();
    constexpr uint32_t none=std::numeric_limits<uint32_t>::max();
    std::vector<double> distance(count,std::numeric_limits<double>::infinity()),rate(count);
    std::vector<uint32_t> source(count,none);
    using Entry=std::tuple<double,uint32_t,uint32_t>;
    std::priority_queue<Entry,std::vector<Entry>,std::greater<Entry>> queue;
    for (size_t e=0;e<table.size();++e) {
        const auto a=table[e].a,b=table[e].b;
        if (owners[a]==owners[b]) continue;
        const double opening=-metric.closing(e,velocity[a],velocity[b])*metric.displacementunit_m()/dt;
        if (opening<=0) continue;
        for (auto i:{a,b}) rate[i]=std::max(rate[i],opening*.5);
    }
    for (uint32_t i=0;i<count;++i) if (rate[i]>0) {
        distance[i]=metric.width_m(i)*.5;source[i]=i;queue.emplace(distance[i],i,i);
    }
    constexpr double extent=800000;
    while (!queue.empty()) {
        const auto [d,s,i]=queue.top();queue.pop();
        if (d!=distance[i] || s!=source[i]) continue;
        for (uint8_t k=0;k<grid::directioncount;++k) {
            const auto j=metric.grid().neighbour(i,k);
            const double next=d+metric.distance_m(metric.edge(i,k));
            if (owners[j]!=owners[i] || next>extent) continue;
            if (next<distance[j] || (next==distance[j] && s<source[j])) {
                distance[j]=next;source[j]=s;queue.emplace(next,s,j);
            }
        }
    }
    // The long-wave background: planar the 5 x 5 box mean; on the cube the area-weighted
    // mean of the two-step ring over the eight-neighbour table (the box's equivalent).
    const auto background=[&](uint32_t i) {
        double mean=0;
        if (metric.planar()) {
            const auto index=[&](int x,int y) { return static_cast<uint32_t>((y%static_cast<int>(h)+h)%h)*w
                +static_cast<uint32_t>((x%static_cast<int>(w)+w)%w); };
            const int x=i%w,y=i/w;
            for (int dy=-2;dy<=2;++dy) for (int dx=-2;dx<=2;++dx) mean+=elevations[index(x+dx,y+dy)]/25;
            return mean;
        }
        std::vector<grid::CellId> ring{i};
        for (int step=0;step<2;++step) {
            const auto reached=ring;
            for (const auto c : reached) for (const auto& n : metric.eight(c))
                if (n.cell!=grid::nocell) ring.push_back(n.cell);
            std::sort(ring.begin(),ring.end());
            ring.erase(std::unique(ring.begin(),ring.end()),ring.end());
        }
        double area=0;
        for (const auto c : ring) { mean+=elevations[c]*metric.area_m2(c); area+=metric.area_m2(c); }
        return mean/area;
    };
    std::vector<Cell> result(count);
    for (uint32_t i=0;i<count;++i) {
        const double mean=background(i);
        // Keep the seeded fine relief; modify only the long-wave background.
        const double detail=elevations[i]-mean;
        double z=elevations[i];
        if (source[i]!=none) {
            const double t=std::min(1.0,distance[i]/extent),weight=(1-t)*(1-t)*(1+2*t);
            const double age=distance[i]/rate[source[i]];
            const double target=elevations[i]>-2000?-200:
                -ridge_depth_m-thermal_strength*thermal::subsidence_m(age);
            // An initial spreading heat source can lift existing ocean floor,
            // never deepen a local depression while trying to retain its detail.
            z+=boundary_influence*weight*(elevations[i]>-2000?target-mean:std::max(0.0,target-mean));
        }
        const double transition=margin_width_m>0?std::clamp(.5+(-z-2500)/margin_width_m,0.0,1.0):0;
        const double ocean_fraction=margin_width_m>0?transition*transition*(3-2*transition):(z>-2500?0:1);
        if (ocean_fraction==0) {
            // Submerged continental basement represents shelves and young rifts.
            result[i]=seed(z,true,ages[i],root_kg_m3,root_max_thickness_m,root_delamination_width_m);
        } else {
            auto c=seed(z,false,0,root_kg_m3,root_max_thickness_m,root_delamination_width_m);
            // Carry high-frequency relief in physical basalt thickness, so it
            // does not disappear as thermal ages approach equilibrium.
            const double anomaly=3500*std::tanh(detail*(mantle_density-water_density)
                /((mantle_density-oceanic_density)*3500));
            c.ocean+=anomaly;
            const double target=std::max(0.0,(elevation_m(c,0,root_kg_m3,root_max_thickness_m,root_delamination_width_m)-z)/std::max(thermal_strength,1e-30));
            double lo=0,hi=1000;
            for (int iteration=0;iteration<60;++iteration) {
                const double mid=(lo+hi)*.5;
                if (thermal::subsidence_m(mid)<target) lo=mid;else hi=mid;
            }
            const double age=thermal_strength>0?(lo+hi)*.5:0;
            c.oceanic_age_volume=c.ocean*age;
            // Residual beyond the invertible cooling range belongs to crust
            // thickness. This does not clip or remap the requested terrain.
            c.ocean+=(z-elevation_m(c,thermal_strength,root_kg_m3,root_max_thickness_m,root_delamination_width_m))*(mantle_density-water_density)
                /(mantle_density-oceanic_density);
            if (c.ocean<=0) throw std::runtime_error("Initial bathymetry requires nonpositive ocean crust");
            c.oceanic_age_volume=c.ocean*age;
            if (ocean_fraction<1) {
                // Both end-member columns reproduce the same submerged height.
                // Their mixture therefore preserves it exactly, while thermal
                // subsidence tapers with oceanic area instead of jumping at -2500 m.
                const auto continent=seed(z,true,ages[i],root_kg_m3,root_max_thickness_m,root_delamination_width_m);
                c.ocean*=ocean_fraction;c.oceanic_age_volume*=ocean_fraction;
                c.oceanic_area=ocean_fraction;c.continental_area=1-ocean_fraction;
                c.continent=continent.continent*c.continental_area;
                c.continental_age_volume=continent.continental_age_volume*c.continental_area;
            }
            result[i]=c;
        }
    }
    return result;
}
} // namespace

std::vector<columns::Cell> terrain(uint32_t w,uint32_t h,double cell_m,double dt,
        const std::vector<double>& elevations,const std::vector<double>& ages,
        const std::vector<uint32_t>& owners,const std::vector<transport::Displacement>& velocity,
        double thermal_strength, double margin_width_m, double boundary_influence, double root_kg_m3,
        double root_max_thickness_m, double root_delamination_width_m, const grid::TectonicGrid* grid) {
    return terrainon(w,h,cell_m,dt,elevations,ages,owners,velocity,thermal_strength,margin_width_m,
        boundary_influence,root_kg_m3,root_max_thickness_m,root_delamination_width_m,grid);
}

std::vector<columns::Cell> terrain(const grid::TectonicGrid& grid,double dt,
        const std::vector<double>& elevations,const std::vector<double>& ages,
        const std::vector<uint32_t>& owners,const std::vector<mesh::Vec3>& velocity_m,
        double thermal_strength, double margin_width_m, double boundary_influence, double root_kg_m3,
        double root_max_thickness_m, double root_delamination_width_m) {
    return terrainon(0,0,0,dt,elevations,ages,owners,velocity_m,thermal_strength,margin_width_m,
        boundary_influence,root_kg_m3,root_max_thickness_m,root_delamination_width_m,&grid);
}
}
