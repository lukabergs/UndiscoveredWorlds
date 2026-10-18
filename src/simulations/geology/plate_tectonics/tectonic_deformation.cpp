#include "simulations/geology/plate_tectonics/tectonic_deformation.hpp"
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

std::vector<transport::Displacement> distribute(uint32_t w, uint32_t h, double cell_m,
        const std::vector<uint32_t>& owners, const std::vector<double>& continent,
        const std::vector<transport::Displacement>& rigid, long seed, const Parameters& p,
        const grid::TectonicGrid* /*grid*/) {
    variation::validate(p.process_variation,p.process_wavelength_km,p.process_period_myr);
    const size_t size=static_cast<size_t>(w)*h;
    if (w<2 || h<2 || owners.size()!=size || continent.size()!=size || rigid.size()!=size ||
        !std::isfinite(cell_m) || cell_m<=0 || !std::isfinite(p.peak_distance_km) || p.peak_distance_km<0 ||
        !std::isfinite(p.distance_variation) || p.distance_variation<0 || p.distance_variation>.5 ||
        !std::isfinite(p.subduction_coupling) || p.subduction_coupling<0 || p.subduction_coupling>1 || !std::isfinite(p.time_myr))
        throw std::runtime_error("Invalid distributed deformation inputs");
    auto result=rigid;
    for (size_t i=0;i<size;++i) if (!std::isfinite(rigid[i].x) || !std::isfinite(rigid[i].y) ||
        !std::isfinite(continent[i]) || continent[i]<0 || continent[i]>1+1e-8)
        throw std::runtime_error("Invalid deformation cell");
    if (p.peak_distance_km==0) return result;
    const auto index=[&](int x,int y) {
        return static_cast<uint32_t>((y+static_cast<int>(h))%static_cast<int>(h))*w
            + static_cast<uint32_t>((x+static_cast<int>(w))%static_cast<int>(w));
    };
    std::vector<transport::Displacement> correction(size);
    std::vector<uint32_t> counts(size,0), source(size,std::numeric_limits<uint32_t>::max());
    std::vector<double> distance(size,std::numeric_limits<double>::infinity()), peak(size,0);
    // Only currently convergent contacts activate a deforming zone. Continental
    // collision shares motion equally; ocean/continent contact transmits only
    // the configured fraction to the upper plate, leaving room for subduction.
    for (uint32_t y=0;y<h;++y) for (uint32_t x=0;x<w;++x) {
        const auto a=index(x,y);
        const uint32_t ns[]={index(x+1,y),index(x,y+1)};
        for (size_t axis=0;axis<2;++axis) {
            const auto b=ns[axis];
            if (owners[a]==owners[b]) continue;
            const double closing=axis==0?rigid[a].x-rigid[b].x:rigid[a].y-rigid[b].y;
            if (closing<=0) continue;
            for (const auto pair : {std::pair<uint32_t,uint32_t>{a,b},{b,a}}) {
                const auto i=pair.first,j=pair.second;
                if (continent[i]<.5 && continent[j]<.5) continue;
                const double coupling=continent[i]>=.5
                    ? (continent[j]>=.5?.5:p.subduction_coupling) : 1-p.subduction_coupling;
                if (coupling==0) continue;
                correction[i].x += coupling*(rigid[j].x-rigid[i].x);
                correction[i].y += coupling*(rigid[j].y-rigid[i].y);
                ++counts[i];
            }
        }
    }
    using Entry=std::tuple<double,uint32_t,uint32_t>; // distance, source, cell: deterministic ties.
    std::priority_queue<Entry,std::vector<Entry>,std::greater<Entry>> queue;
    constexpr double pi=3.14159265358979323846;
    const auto bits=static_cast<uint32_t>(seed)*1664525U+1013904223U;
    const double phase=2*pi*bits/4294967296.0;
    for (uint32_t i=0;i<size;++i) if (counts[i]) {
        correction[i].x/=counts[i]; correction[i].y/=counts[i];
        const double x=(i%w+.5)/w,y=(i/w+.5)/h;
        const double variation=.5*(std::sin(2*pi*(4*x+2*y)+phase)+std::cos(2*pi*(2*x-3*y)-phase));
        peak[i]=p.peak_distance_km*1000*(1+p.distance_variation*variation);
        peak[i]*=tectonics::variation::multiplier(x,y,w*cell_m/1000,h*cell_m/1000,p.time_myr,seed,4,
            p.process_variation,p.process_wavelength_km,p.process_period_myr);
        distance[i]=.5*cell_m; source[i]=i;
        queue.emplace(distance[i],i,i);
    }
    const double extent=2*p.peak_distance_km*1000*(1+p.distance_variation)*(1+p.process_variation);
    while (!queue.empty()) {
        const auto [d,s,i]=queue.top(); queue.pop();
        if (d!=distance[i] || s!=source[i]) continue;
        const auto x=i%w,y=i/w;
        for (int dy=-1;dy<=1;++dy) for (int dx=-1;dx<=1;++dx) {
            if (dx==0 && dy==0) continue;
            const auto j=index(static_cast<int>(x)+dx,static_cast<int>(y)+dy);
            if (owners[j]!=owners[i] || (continent[j]>=.5)!=(continent[i]>=.5)) continue;
            // No diagonal shortcut through water or another plate.
            if (dx && dy && (owners[index(static_cast<int>(x)+dx,y)]!=owners[i] ||
                owners[index(x,static_cast<int>(y)+dy)]!=owners[i] ||
                (continent[index(static_cast<int>(x)+dx,y)]>=.5)!=(continent[i]>=.5) ||
                (continent[index(x,static_cast<int>(y)+dy)]>=.5)!=(continent[i]>=.5))) continue;
            const double next=d+cell_m*(dx && dy?std::sqrt(2.0):1.0);
            if (next>extent) continue;
            if (next<distance[j] || (next==distance[j] && s<source[j])) {
                distance[j]=next; source[j]=s; queue.emplace(next,s,j);
            }
        }
    }
    for (size_t i=0;i<size;++i) if (source[i]!=std::numeric_limits<uint32_t>::max()) {
        const auto s=source[i]; const double f=taper(distance[i],peak[s]);
        result[i].x+=f*correction[s].x; result[i].y+=f*correction[s].y;
    }
    return result;
}
}
