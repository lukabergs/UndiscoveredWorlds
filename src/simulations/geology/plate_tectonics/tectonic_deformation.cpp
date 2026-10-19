#include "simulations/geology/plate_tectonics/tectonic_deformation.hpp"
#include "simulations/geology/plate_tectonics/tectonic_columns.hpp"
#include "simulations/geology/plate_tectonics/tectonic_grid.hpp"
#include "simulations/geology/plate_tectonics/tectonic_variation.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <stdexcept>
#include <tuple>

namespace tectonics::deformation {
double taper(double distance, double peak) {
    if (!std::isfinite(distance) || distance < 0 || !std::isfinite(peak) || peak <= 0)
        throw std::runtime_error("Invalid deformation profile distance");
    const double x=std::min(1.0,distance/(2*peak));
    return (1-x)*(1-x)*(1+2*x);
}

namespace {
bool finite(const mesh::Vec3& v) { return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z); }
void addscaled(mesh::Vec3& to, double k, const mesh::Vec3& v) { to.x += k*v.x; to.y += k*v.y; to.z += k*v.z; }
void adddifference(mesh::Vec3& to, double k, const mesh::Vec3& a, const mesh::Vec3& b) {
    to.x += k*(a.x-b.x); to.y += k*(a.y-b.y); to.z += k*(a.z-b.z);
}
void divide(mesh::Vec3& v, uint32_t n) { v.x/=n; v.y/=n; v.z/=n; }
} // namespace

std::vector<mesh::Vec3> distribute(const columns::Metric& metric, const std::vector<uint32_t>& owners,
        const std::vector<double>& continent, const std::vector<mesh::Vec3>& rigid, long seed,
        const Parameters& p) {
    using V = mesh::Vec3;
    variation::validate(p.process_variation,p.process_wavelength_km,p.process_period_myr);
    const size_t size=metric.cellcount();
    if (size==0 || owners.size()!=size ||
        continent.size()!=size || rigid.size()!=size || !std::isfinite(p.peak_distance_km) || p.peak_distance_km<0 ||
        !std::isfinite(p.distance_variation) || p.distance_variation<0 || p.distance_variation>.5 ||
        !std::isfinite(p.subduction_coupling) || p.subduction_coupling<0 || p.subduction_coupling>1 || !std::isfinite(p.time_myr))
        throw std::runtime_error("Invalid distributed deformation inputs");
    auto result=rigid;
    for (size_t i=0;i<size;++i) if (!finite(rigid[i]) ||
        !std::isfinite(continent[i]) || continent[i]<0 || continent[i]>1+1e-8)
        throw std::runtime_error("Invalid deformation cell");
    if (p.peak_distance_km==0) return result;
    // The cell table (W2-E): contacts on its edges, displacements through the adapter
    // (a neighbour's on the cell's own frame), distances in metres.
    const auto& table=metric.grid().edges();
    std::vector<V> correction(size);
    std::vector<uint32_t> counts(size,0), source(size,std::numeric_limits<uint32_t>::max());
    std::vector<double> distance(size,std::numeric_limits<double>::infinity()), peak(size,0);
    // Only currently convergent contacts activate a deforming zone. Continental
    // collision shares motion equally; ocean/continent contact transmits only
    // the configured fraction to the upper plate, leaving room for subduction.
    for (size_t e=0;e<table.size();++e) {
        const auto a=table[e].a,b=table[e].b;
        if (owners[a]==owners[b]) continue;
        const double closing=metric.closing(e,rigid[a],rigid[b]);
        if (closing<=0) continue;
        for (const auto pair : {std::pair<uint32_t,uint32_t>{a,b},{b,a}}) {
            const auto i=pair.first,j=pair.second;
            if (continent[i]<.5 && continent[j]<.5) continue;
            const double coupling=continent[i]>=.5
                ? (continent[j]>=.5?.5:p.subduction_coupling) : 1-p.subduction_coupling;
            if (coupling==0) continue;
            adddifference(correction[i],coupling,metric.inframe(i,rigid[j]),rigid[i]);
            ++counts[i];
        }
    }
    using Entry=std::tuple<double,uint32_t,uint32_t>; // distance, source, cell: deterministic ties.
    std::priority_queue<Entry,std::vector<Entry>,std::greater<Entry>> queue;
    constexpr double pi=3.14159265358979323846;
    const auto bits=static_cast<uint32_t>(seed)*1664525U+1013904223U;
    const double phase=2*pi*bits/4294967296.0;
    // E4: replace. The peak's variation on the 4N x 2N lat-lon proxy.
    const auto extent_km=metric.extent_km();
    for (uint32_t i=0;i<size;++i) if (counts[i]) {
        divide(correction[i],counts[i]);
        const auto position=metric.unitposition(i);
        const double x=position[0],y=position[1];
        const double variation=.5*(std::sin(2*pi*(4*x+2*y)+phase)+std::cos(2*pi*(2*x-3*y)-phase));
        peak[i]=p.peak_distance_km*1000*(1+p.distance_variation*variation);
        peak[i]*=tectonics::variation::multiplier(x,y,extent_km[0],extent_km[1],p.time_myr,seed,4,
            p.process_variation,p.process_wavelength_km,p.process_period_myr);
        distance[i]=.5*metric.width_m(i); source[i]=i;
        queue.emplace(distance[i],i,i);
    }
    const double extent=2*p.peak_distance_km*1000*(1+p.distance_variation)*(1+p.process_variation);
    // The eight neighbours (NW, N, NE, W, E, SW, S, SE); a diagonal's two cardinals.
    constexpr int sides[8][2]={{3,1},{-1,-1},{4,1},{-1,-1},{-1,-1},{3,6},{-1,-1},{4,6}};
    const auto same=[&](uint32_t i,uint32_t j) {
        return owners[j]==owners[i] && (continent[j]>=.5)==(continent[i]>=.5);
    };
    while (!queue.empty()) {
        const auto [d,s,i]=queue.top(); queue.pop();
        if (d!=distance[i] || s!=source[i]) continue;
        const auto neighbours=metric.eight(i);
        for (int k=0;k<8;++k) {
            const auto j=neighbours[k].cell;
            if (j==grid::nocell || !same(i,j)) continue;
            // No diagonal shortcut through water or another plate.
            if (sides[k][0]>=0 && (!same(i,neighbours[sides[k][0]].cell) || !same(i,neighbours[sides[k][1]].cell))) continue;
            const double next=d+neighbours[k].distance_m;
            if (next>extent) continue;
            if (next<distance[j] || (next==distance[j] && s<source[j])) {
                distance[j]=next; source[j]=s; queue.emplace(next,s,j);
            }
        }
    }
    for (size_t i=0;i<size;++i) if (source[i]!=std::numeric_limits<uint32_t>::max()) {
        const auto s=source[i]; const double f=taper(distance[i],peak[s]);
        addscaled(result[i],f,metric.inframe(static_cast<grid::CellId>(i),correction[s]));
    }
    return result;
}
}
