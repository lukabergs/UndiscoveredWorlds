#include "simulations/landscape/relief_roughness.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

#include "simulations/landscape/spherical_noise.hpp"

namespace landscape
{
namespace
{
// The octave sum is divided by its own standard deviation (the kernel's one-octave
// deviation, spherical_noise.hpp) and by `peakdeviations`, so that amplitude_m is the value
// reached at that many standard deviations: with every octave resolved, the
// largest magnitude on a grid of 64 to 1024 columns was 0.9 to 1.2 x amplitude_m
// (C3 measurement, before the band limit below).
constexpr double peakdeviations = 3.5;
}

void addreliefroughness(const SphereMetric& metric, const ReliefRoughness& roughness, const std::vector<uint8_t>& sea,
                        float sealevel, std::vector<float>& elevation)
{
    if (roughness.amplitude_m == 0.0)
        return;
    // The noise lives on the sphere; a torus metric has no position on it (erosion plan D2a).
    if (metric.wrapsrows)
        throw std::invalid_argument("addreliefroughness: a torus metric is not supported");
    const std::size_t count = metric.size();
    if (sea.size() != count || elevation.size() != count)
        throw std::invalid_argument("addreliefroughness: input sizes differ from the metric");
    if (!(roughness.wavelength_km > 0.0) || roughness.octaves < 1)
        throw std::invalid_argument("addreliefroughness: wavelength_km must be positive and octaves at least 1");

    // Octave k samples the noise at the cell centre's position on the sphere in
    // units of its wavelength (wavelength_km / 2^k is the lattice spacing), with
    // its own seed and lattice offset; its weight is persistence^k. The weighted
    // sum is divided by the standard deviation of the full octave set,
    // octavedeviation x sqrt(sum of the squared weights). An octave whose
    // wavelength is below twice the north-south cell size is not resolved by the
    // grid (it would act as per-cell noise whose slope grows with the width): its
    // weight is 0, the normalisation unchanged.
    FractalNoise noise;
    noise.seed = roughness.seed;
    noise.radius_km = metric.radius_m / 1000.0;
    noise.wavelength_km = roughness.wavelength_km;
    noise.octaves = roughness.octaves;
    noise.persistence = roughness.persistence;
    noise.shortest_km = 2.0 * metric.dy_m / 1000.0;
    const FractalNoiseField field(noise);
    if (field.resolved() == 0)
        return;
    // The operations and their order of step C3 (products are not associative), with
    // the kernel's own one-octave deviation (octavedeviation for Perlin, so
    // bit-identical there; multi-scale plan S3.16).
    const double factor =
        roughness.amplitude_m / (peakdeviations * field.deviation() * std::sqrt(field.squaredweights()));

    const double pi = 3.14159265358979323846;
    const double floor = static_cast<double>(sealevel) + 1.0;
    // The noise position of a cell: its centre on a cube or a subgrid; on a full
    // latitude-longitude grid (or a gridless metric of rows) the row-centre latitude
    // and the longitude 2 pi (x + 0.5) / W, which runs from 0 at the west edge of
    // column 0 (the grid's own longitudes run from -pi: the roughness of a lat-lon
    // world is kept as it was).
    const bool rows = !metric.grid || (metric.grid->kind() == mesh::GridKind::latlon && !metric.grid->issubgrid());
    for (std::size_t cell = 0; cell < count; cell++)
    {
        if (sea[cell] != 0)
            continue;
        double px = 0.0, py = 0.0, pz = 0.0;
        if (rows)
        {
            const double latitude = metric.rowlatitude_rad[static_cast<std::size_t>(metric.row(cell))];
            const double longitude = 2.0 * pi * (static_cast<double>(metric.column(cell)) + 0.5) / metric.width;
            unitvector(latitude, longitude, px, py, pz);
        }
        else
        {
            const mesh::Vec3 centre = metric.grid->centre(static_cast<mesh::CellId>(cell));
            px = centre.x;
            py = centre.y;
            pz = centre.z;
        }
        const double sum = field.sum(px, py, pz);
        // A lowering tapers towards the coast: it is scaled by
        // min(1, (old - floor) / amplitude), so low land is not flattened onto
        // floor = sealevel + 1 m. The floor stays as a backstop for the rare
        // perturbation beyond the amplitude: a land cell is not lowered to it or
        // below (one already there is not lowered at all), so it stays land.
        const double old = static_cast<double>(elevation[cell]);
        double delta = factor * sum;
        if (delta < 0.0)
            delta *= std::clamp((old - floor) / std::fabs(roughness.amplitude_m), 0.0, 1.0);
        const double lowest = std::min(old, floor);
        elevation[cell] = static_cast<float>(std::max(old + delta, lowest));
    }
}
}
