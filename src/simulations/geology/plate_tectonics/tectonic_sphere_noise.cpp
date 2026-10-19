#include "simulations/geology/plate_tectonics/tectonic_sphere_noise.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace tectonics::sphere_noise {
namespace {
landscape::FractalNoise settings(std::uint64_t seed, double wavelength_km, double radius_km, int octaves,
                                 double persistence, double shortest_km) {
    if (!std::isfinite(wavelength_km) || wavelength_km <= 0 || !std::isfinite(radius_km) || radius_km <= 0 ||
        octaves < 1 || !std::isfinite(persistence) || persistence <= 0 || !std::isfinite(shortest_km) || shortest_km < 0)
        throw std::invalid_argument("sphere_noise::Octaves: positive wavelength, radius, octaves and persistence");
    landscape::FractalNoise noise;
    noise.seed = seed;
    noise.radius_km = radius_km;
    noise.wavelength_km = std::max(wavelength_km, shortest_km);
    noise.octaves = octaves;
    noise.persistence = persistence;
    noise.shortest_km = shortest_km;
    // The project's isotropic choice (multi-scale plan S3.16), fixed here rather than read from
    // landscape.noise: unrotated Perlin octaves share one lattice with an axis on the pole axis,
    // whose orientation shows at the seeding's wavelengths (up to several times the sphere).
    noise.kernel = landscape::NoiseKernel::opensimplex2s;
    noise.rotateoctaves = 1;
    return noise;
}
}

Octaves::Octaves(std::uint64_t seed, double wavelength_km, double radius_km, int octaves, double persistence,
                 double deviation, double shortest_km)
    : field_(settings(seed, wavelength_km, radius_km, octaves, persistence, shortest_km)) {
    if (!std::isfinite(deviation) || deviation <= 0)
        throw std::invalid_argument("sphere_noise::Octaves: positive deviation");
    // The kept octaves' weights, as FractalNoiseField keeps them (wavelength at least shortest_km).
    const double longest = std::max(wavelength_km, shortest_km);
    double weights = 0;
    for (int o = 0; o < octaves; ++o)
        if (longest / std::ldexp(1.0, o) >= shortest_km) weights += std::pow(persistence, o);
    factor_ = deviation / (field_.deviation() * weights);
}

double Octaves::operator()(const mesh::Vec3& p) const {
    return std::clamp(factor_ * field_.sum(p.x, p.y, p.z), -1.0, 1.0);
}

mesh::Vec3 warp(const mesh::Vec3& p, const Octaves& x, const Octaves& y, const Octaves& z, double amplitude_km,
                double radius_km) {
    const mesh::Vec3 d = amplitude_km * mesh::Vec3{x(p), y(p), z(p)};
    const mesh::Vec3 t = d - mesh::dot(d, p) * p;
    const double length = mesh::length(t);
    if (length == 0) return p;
    const double angle = length / radius_km;
    return mesh::normalise(std::cos(angle) * p + (std::sin(angle) / length) * t);
}

std::array<mesh::Vec3, 3> frame(std::uint64_t seed) {
    const auto m = landscape::octaverotation(seed);
    return {mesh::Vec3{m[0], m[1], m[2]}, mesh::Vec3{m[3], m[4], m[5]}, mesh::Vec3{m[6], m[7], m[8]}};
}
}
