#include "simulations/geology/plate_tectonics/tectonic_ocean_relief.hpp"
#include "simulations/geology/plate_tectonics/tectonic_columns.hpp"
#include "simulations/geology/plate_tectonics/tectonic_sphere_noise.hpp"
#include "core/mesh/mesh_geometry.hpp"
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
std::vector<double> crust(const mesh::CellGrid& grid, long seed, const Parameters& p) {
    validate(p);
    const size_t count=grid.cellcount();
    if (count==0 || !std::isfinite(grid.radius_m()) || grid.radius_m()<=0)
        throw std::runtime_error("Invalid ocean relief grid");
    const double radius_km=grid.radius_m()/1000;
    const double spacing_km=std::sqrt(grid.totalarea_m2()/static_cast<double>(count))/1000;
    // Omit unresolved octaves (shorter than two cell spacings); never alias fine noise into
    // isolated cells.
    const auto key=static_cast<uint32_t>(seed);
    const sphere_noise::Octaves relief(key^0xc2b2ae35U,p.wavelength_km,radius_km,3,.5,
                                       sphere_noise::value_deviation,2*spacing_km);
    const sphere_noise::Octaves roughness(key^0x94d049bbU,p.wavelength_km,radius_km,1,.5,
                                          sphere_noise::value_deviation,2*spacing_km);
    std::vector<double> seamount(count,0.0);
    if (p.seamount_crust_m>0) {
        // One seamount per cell of a coarse cube of cells about the spacing wide.
        const double face_km=std::sqrt(grid.totalarea_m2()/6)/1000;
        const int n=std::max(1,static_cast<int>(std::lround(face_km/p.seamount_spacing_km)));
        const auto coarse=mesh::makecubegrid(n,grid.radius_m());
        const auto salt=key^0x27d4eb2fU;
        const double reach_m=3*p.seamount_radius_km*1.25*1000;
        for (mesh::CellId k=0;k<coarse->cellcount();++k) {
            const double width_km=std::sqrt(coarse->area_m2(k))/1000;
            const auto cell=mesh::tangentframe(coarse->centre(k),grid.radius_m());
            const mesh::Vec3 centre=mesh::fromtangent(cell,(.6*hash(k,0,salt)-.3)*width_km,
                                                      (.6*hash(k,1,salt)-.3)*width_km);
            const auto frame=mesh::tangentframe(centre,grid.radius_m());
            const double angle=6.283185307179586*hash(k,2,salt);
            const double aspect=.75+.5*hash(k,3,salt);
            const double height=.8+.4*hash(k,4,salt);
            for (const auto c : mesh::cellswithin(grid,centre,reach_m)) {
                // The cell centre in the seamount's tangent frame by great-circle distance.
                const mesh::Vec3 q=grid.centre(c);
                const double distance_km=mesh::greatcircle_rad(centre,q)*radius_km;
                const mesh::Vec3 t=q-mesh::dot(q,centre)*centre;
                const double length=mesh::length(t);
                const double px=length>0?distance_km*mesh::dot(t,frame.east)/length:0;
                const double py=length>0?distance_km*mesh::dot(t,frame.north)/length:0;
                const double rx=(px*std::cos(angle)+py*std::sin(angle))/(p.seamount_radius_km*aspect);
                const double ry=(-px*std::sin(angle)+py*std::cos(angle))/(p.seamount_radius_km/aspect);
                const double radius=std::hypot(rx,ry);
                if (radius<3) seamount[c]+=std::exp(-.5*radius*radius)*(1-smooth(std::clamp(radius-2,0.0,1.0)))*height;
            }
        }
    }
    std::vector<double> result(count);
    for (mesh::CellId c=0;c<count;++c) {
        const mesh::Vec3 q=grid.centre(c);
        const double variation=p.crust_variation_m>0?p.crust_variation_m*relief(q):0;
        // Distort the initial/birth profile, not evolved terrain. The same
        // positive rock layer subsequently advects; no per-frame noise state.
        const double raised=seamount[c]>0?seamount[c]*(1+p.seamount_roughness*roughness(q)):0;
        result[c]=columns::ridge_crust_m+variation+p.seamount_crust_m*raised;
    }
    return result;
}
}
