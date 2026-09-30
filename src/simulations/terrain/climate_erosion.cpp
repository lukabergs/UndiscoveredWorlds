#include "climate_erosion.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace climateerosion {
Result evolve(uint32_t w,uint32_t h,const std::vector<double>& initial,const Forcing& f,const Parameters& p) {
    const size_t size=static_cast<size_t>(w)*h;
    if (w!=2*h || h<2 || initial.size()!=size || p.duration_myr<=0 || p.step_myr<=0 || p.radius_m<=0 || p.sediment_travel_km<=0)
        throw std::runtime_error("Invalid climate erosion geometry or time");
    for (double v:{p.duration_myr,p.step_myr,p.radius_m,p.sediment_travel_km,p.river_coefficient,
            p.wave_m_myr,p.ice_m_myr,p.wind_m_myr,p.hillslope_m2_myr})
        if (!std::isfinite(v) || v<0) throw std::runtime_error("Invalid climate erosion parameter");
    for (const auto* field:{&initial,&f.rain_mm_year,&f.temperature_c,&f.wind_m_s,&f.snow_fraction,&f.ice_fraction}) {
        if (field->size()!=size) throw std::runtime_error("Climate erosion field shape mismatch");
        for (double v:*field) if (!std::isfinite(v)) throw std::runtime_error("Nonfinite climate erosion input");
    }
    constexpr double pi=3.14159265358979323846;
    const double dy=pi*p.radius_m/h;
    std::vector<double> area(size),dx(size),runoff(size),load(size),slope(size),length(size),delta(size);
    std::vector<size_t> order(size),down(size);
    for (uint32_t y=0;y<h;++y) for (uint32_t x=0;x<w;++x) {
        const auto i=static_cast<size_t>(y)*w+x;
        const double north=pi/2-pi*y/h,south=pi/2-pi*(y+1)/h;
        area[i]=p.radius_m*p.radius_m*(2*pi/w)*(std::sin(north)-std::sin(south));
        dx[i]=2*pi*p.radius_m/w*std::cos((north+south)/2);
    }
    const auto neighbor=[&](size_t i,int ox,int oy) {
        const int y=static_cast<int>(i/w)+oy;
        if (y<0 || y>=static_cast<int>(h)) return i; // Closed polar faces, periodic longitude.
        return static_cast<size_t>(y)*w+(static_cast<int>(i%w)+ox+static_cast<int>(w))%w;
    };
    Result result;result.elevation_m=initial;result.erosion_m.assign(size,0);result.deposition_m.assign(size,0);
    for (double time=0;time<p.duration_myr;) {
        const double dt=std::min(p.step_myr,p.duration_myr-time);time+=dt;
        auto& z=result.elevation_m;
        std::fill(load.begin(),load.end(),0);std::fill(slope.begin(),slope.end(),0);std::fill(delta.begin(),delta.end(),0);
        std::iota(order.begin(),order.end(),0);
        for (size_t i=0;i<size;++i) {
            down[i]=i;
            runoff[i]=z[i]>0?area[i]/1e6*std::max(0.0,f.rain_mm_year[i])/1000*
                (1-.7*std::clamp(f.ice_fraction[i],0.0,1.0))*(1-.4*std::clamp(f.snow_fraction[i],0.0,1.0)):0;
            for (int oy=-1;oy<=1;++oy) for (int ox=-1;ox<=1;++ox) {
                if (!ox && !oy) continue;
                const auto j=neighbor(i,ox,oy);if (j==i) continue;
                const double distance=std::hypot(ox*.5*(dx[i]+dx[j]),oy*dy);
                const double gradient=(z[i]-z[j])/distance;
                if (gradient>slope[i]) {slope[i]=gradient;length[i]=distance;down[i]=j;}
            }
        }
        std::sort(order.begin(),order.end(),[&](size_t a,size_t b){return z[a]!=z[b]?z[a]>z[b]:a<b;});
        for (auto i:order) {
            const auto j=down[i];double removed=0;
            if (z[i]>0 && j!=i) {
                int ocean_faces=0;
                for (const auto offset:{std::pair<int,int>{-1,0},{1,0},{0,-1},{0,1}})
                    ocean_faces+=z[neighbor(i,offset.first,offset.second)]<=0;
                const double rain=std::max(0.0,f.rain_mm_year[i]);
                const double wind=std::max(0.0,f.wind_m_s[i]);
                const double frost=std::exp(-std::pow(f.temperature_c[i]/8,2));
                const double river=p.river_coefficient*std::sqrt(runoff[i])*slope[i]*(1+.5*frost);
                // Warm-based ice potential only: a coarse slope-driven abrasion
                // closure; permanent cold polar ice is protected from full erosion.
                const double warm_bed=1/(1+std::exp(-(f.temperature_c[i]+12)/3));
                const double glacier=p.ice_m_myr*std::clamp(f.ice_fraction[i],0.0,1.0)*warm_bed*
                    std::sqrt(slope[i]/.01);
                const double aeolian=p.wind_m_myr*std::pow(wind/8,3)*std::exp(-rain/400)*
                    (1-std::clamp(f.snow_fraction[i],0.0,1.0))*std::min(1.0,slope[i]/.005);
                const double waves=p.wave_m_myr*.25*ocean_faces*std::exp(-z[i]/200)*(.5+.5*std::min(3.0,wind/8));
                const double incision=std::min(.2*(z[i]-z[j]),dt*(river+glacier+aeolian+waves));
                removed=incision*area[i];delta[i]-=incision;result.erosion_m[i]+=incision;result.eroded_m3+=removed;
            }
            const double sediment=load[i]+removed;
            const double retained=j==i?1:z[i]>0?0:-std::expm1(-length[i]/(p.sediment_travel_km*1000));
            const double deposited=sediment*retained;
            delta[i]+=deposited/area[i];result.deposition_m[i]+=deposited/area[i];result.deposited_m3+=deposited;
            if (j!=i) {load[j]+=sediment-deposited;runoff[j]+=runoff[i];}
        }
        // Conservative hillslope transport only on land, with a per-face flux
        // limiter for steep/polar cells; no global smoothing of ocean detail.
        for (size_t i=0;i<size;++i) for (const auto offset:{std::pair<int,int>{1,0},{0,1}}) {
            const auto j=neighbor(i,offset.first,offset.second);
            if (j==i || z[i]<=0 || z[j]<=0) continue;
            const size_t from=z[i]>z[j]?i:j,to=from==i?j:i;
            const double distance=offset.first?.5*(dx[i]+dx[j]):dy;
            const double face=offset.first?dy:.5*(dx[i]+dx[j]);
            const double volume=std::min(.05*(z[from]-z[to])*std::min(area[i],area[j]),
                p.hillslope_m2_myr*dt*(z[from]-z[to])/distance*face);
            delta[from]-=volume/area[from];delta[to]+=volume/area[to];
            result.erosion_m[from]+=volume/area[from];result.deposition_m[to]+=volume/area[to];
            result.eroded_m3+=volume;result.deposited_m3+=volume;
        }
        for (size_t i=0;i<size;++i) {
            z[i]+=delta[i];
            if (!std::isfinite(z[i])) throw std::runtime_error("Nonfinite climate erosion result");
        }
    }
    double residual=0;
    for (size_t i=0;i<size;++i) residual+=(result.elevation_m[i]-initial[i])*area[i];
    result.relative_volume_residual=residual/std::max(1.0,result.eroded_m3);
    return result;
}
}
