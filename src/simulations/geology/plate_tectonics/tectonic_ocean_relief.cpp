#include "simulations/geology/plate_tectonics/tectonic_ocean_relief.hpp"
#include "simulations/geology/plate_tectonics/tectonic_columns.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace tectonics::ocean_relief {
namespace {
double hash(uint32_t x, uint32_t y, uint32_t seed) {
    uint32_t value = seed ^ (x * 0x9e3779b9U) ^ (y * 0x85ebca6bU);
    value ^= value >> 16; value *= 0x7feb352dU;
    value ^= value >> 15; value *= 0x846ca68bU; value ^= value >> 16;
    return value / 4294967295.0;
}
double smooth(double t) { return t*t*t*(t*(t*6-15)+10); }
}
void validate(const Parameters& p) {
    if (!std::isfinite(p.crust_variation_m) || p.crust_variation_m < 0 ||
        p.crust_variation_m >= columns::ridge_crust_m ||
        !std::isfinite(p.wavelength_km) || p.wavelength_km <= 0 ||
        !std::isfinite(p.seamount_crust_m) || p.seamount_crust_m < 0 ||
        !std::isfinite(p.seamount_radius_km) || p.seamount_radius_km <= 0 ||
        !std::isfinite(p.seamount_spacing_km) || p.seamount_spacing_km < 6*p.seamount_radius_km ||
        !std::isfinite(p.seamount_roughness) || p.seamount_roughness<0 || p.seamount_roughness>1)
        throw std::runtime_error("Invalid ocean relief parameters: require positive scales, nonnegative amplitudes and spacing >= 6 radii");
}
double noise(double u, double v, uint32_t nx, uint32_t ny, uint32_t seed) {
    if (!std::isfinite(u) || !std::isfinite(v) || nx == 0 || ny == 0)
        throw std::runtime_error("Invalid ocean noise coordinates");
    const double x=(u-std::floor(u))*nx, y=(v-std::floor(v))*ny;
    const auto ix=static_cast<uint32_t>(std::floor(x)), iy=static_cast<uint32_t>(std::floor(y));
    const double tx=smooth(x-ix), ty=smooth(y-iy);
    const double a=hash(ix%nx,iy%ny,seed)*(1-tx)+hash((ix+1)%nx,iy%ny,seed)*tx;
    const double b=hash(ix%nx,(iy+1)%ny,seed)*(1-tx)+hash((ix+1)%nx,(iy+1)%ny,seed)*tx;
    return 2*(a*(1-ty)+b*ty)-1;
}
std::vector<double> crust(uint32_t w, uint32_t h, double cell_m, long seed, const Parameters& p) {
    validate(p);
    if (w<2 || h<2 || !std::isfinite(cell_m) || cell_m<=0)
        throw std::runtime_error("Invalid ocean relief grid");
    // Omit unresolved octaves; never alias fine noise into isolated pixels.
    const auto nx=static_cast<uint32_t>(std::clamp(std::round(w*cell_m/(p.wavelength_km*1000)),1.0,w/2.0));
    const auto ny=static_cast<uint32_t>(std::clamp(std::round(h*cell_m/(p.wavelength_km*1000)),1.0,h/2.0));
    const double width_km=w*cell_m/1000, height_km=h*cell_m/1000;
    const auto sx=static_cast<int>(std::max(1.0,std::round(width_km/p.seamount_spacing_km)));
    const auto sy=static_cast<int>(std::max(1.0,std::round(height_km/p.seamount_spacing_km)));
    std::vector<double> result(static_cast<size_t>(w)*h);
    for (uint32_t y=0;y<h;++y) for (uint32_t x=0;x<w;++x) {
        double value=0,weight=1,total=0;
        for (uint32_t octave=0;octave<3;++octave) {
            if ((nx<<octave)>w/2 || (ny<<octave)>h/2) break;
            value+=weight*noise((x+.5)/w,(y+.5)/h,nx<<octave,ny<<octave,
                               static_cast<uint32_t>(seed)^0xc2b2ae35U^(octave*0x9e3779b9U));
            total+=weight; weight*=.5;
        }
        double seamount=0;
        if (p.seamount_crust_m>0) {
            const double px=(x+.5)*cell_m/1000, py=(y+.5)*cell_m/1000;
            const int tx=static_cast<int>(px/width_km*sx),ty=static_cast<int>(py/height_km*sy);
            for (int dy=-1;dy<=1;++dy) for (int dx=-1;dx<=1;++dx) {
                const auto ix=static_cast<uint32_t>((tx+dx+sx)%sx),iy=static_cast<uint32_t>((ty+dy+sy)%sy);
                const auto key=static_cast<uint32_t>(seed)^0x27d4eb2fU;
                const double cx=(tx+dx+.2+.6*hash(ix,iy,key))*width_km/sx;
                const double cy=(ty+dy+.2+.6*hash(ix,iy,key^0x165667b1U))*height_km/sy;
                const double angle=6.283185307179586*hash(ix,iy,key^0xd3a2646cU);
                const double aspect=.75+.5*hash(ix,iy,key^0xfd7046c5U);
                const double rx=((px-cx)*std::cos(angle)+(py-cy)*std::sin(angle))/(p.seamount_radius_km*aspect);
                const double ry=(-(px-cx)*std::sin(angle)+(py-cy)*std::cos(angle))/(p.seamount_radius_km/aspect);
                const double radius=std::hypot(rx,ry);
                if (radius<3) seamount+=std::exp(-.5*radius*radius)*(1-smooth(std::clamp(radius-2,0.0,1.0)))
                    *(.8+.4*hash(ix,iy,key^0xb55a4f09U));
            }
        }
        // Distort the initial/birth profile, not evolved terrain. The same
        // positive rock layer subsequently advects; no per-frame noise state.
        seamount*=1+p.seamount_roughness*noise((x+.5)/w,(y+.5)/h,nx,ny,
                                             static_cast<uint32_t>(seed)^0x94d049bbU);
        result[static_cast<size_t>(y)*w+x]=columns::ridge_crust_m+p.crust_variation_m*value/total
            +p.seamount_crust_m*seamount;
    }
    return result;
}
}
