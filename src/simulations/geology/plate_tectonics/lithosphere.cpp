/******************************************************************************
 *  plate-tectonics, a plate tectonics simulation library
 *  Copyright (C) 2012-2013 Lauri Viitanen
 *  Copyright (C) 2014-2015 Federico Tomassetti, Bret Curtis
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, see http://www.gnu.org/licenses/
 *****************************************************************************/

#include "simulations/geology/plate_tectonics/lithosphere.hpp"
#include "simulations/geology/plate_tectonics/tectonic_sphere_noise.hpp"
#include "simulations/geology/plate_tectonics/tectonic_boundary_graph.hpp"
#include "simulations/geology/plate_tectonics/tectonic_deforming_regions.hpp"
#include "simulations/geology/plate_tectonics/tectonic_junction_coupling.hpp"
#include "simulations/geology/plate_tectonics/tectonic_junctions.hpp"
#include "simulations/geology/plate_tectonics/tectonic_provenance.hpp"
#include "simulations/geology/plate_tectonics/tectonic_thermal.hpp"
#include "simulations/geology/plate_tectonics/tectonic_initialization.hpp"
#include "simulations/geology/plate_tectonics/tectonic_deformation.hpp"
#include "core/mesh/mesh_geometry.hpp"

#include <algorithm>
#include <array>
#include <cfloat>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <iostream>
#include <limits>
#include <map>
#include <tuple>
#include <utility>
#include <vector>


using namespace std;

static constexpr uint32_t MAX_BUOYANCY_AGE = 20;





namespace {

// The torus engine's seeding noise on the sphere (E4): a field keeps the torus' octaves and
// persistence; its frequency f (the torus' scale x noise scale x the coordinate factor it was
// read at) sets the wavelength unit_km / f (scenario.initial_noise_wavelength_km).
tectonics::sphere_noise::Octaves seeding_noise(uint32_t seed, double unit_km, double radius_km, int octaves,
                                               double persistence, double frequency) {
    return tectonics::sphere_noise::Octaves(seed, unit_km / frequency, radius_km, octaves, persistence);
}

float clamp_unit(float value) {
    return std::max(0.0f, std::min(1.0f, value));
}

float interpolate(float start, float end, float alpha) {
    return start + (end - start) * clamp_unit(alpha);
}

float smoothstep(float edge0, float edge1, float value) {
    if (std::fabs(edge1 - edge0) <= FLT_EPSILON) {
        return value >= edge1 ? 1.0f : 0.0f;
    }
    const float t = clamp_unit((value - edge0) / (edge1 - edge0));
    return t * t * (3.0f - 2.0f * t);
}

float logistic(float value, float sharpness) {
    return 1.0f / (1.0f + std::exp(-value * sharpness));
}

float remap_to_range(float value, float min_value, float max_value, float target_min,
                     float target_max) {
    if (max_value - min_value <= FLT_EPSILON) {
        return 0.5f * (target_min + target_max);
    }
    return target_min +
           ((value - min_value) / (max_value - min_value)) * (target_max - target_min);
}

float internal_to_physical_meters(float internal, uint16_t sea_level_m) {
    if (internal <= TopographyCodec::kOceanicBase) {
        return 0.0f;
    }

    if (TopographyCodec::is_oceanic_internal(internal)) {
        if (sea_level_m <= 1U) {
            return 0.0f;
        }

        const float oceanic_ceiling =
            std::nextafter(TopographyCodec::kContinentalBase, TopographyCodec::kOceanicBase);
        const float ratio =
            (std::min(oceanic_ceiling, internal) - TopographyCodec::kOceanicBase) /
            (oceanic_ceiling - TopographyCodec::kOceanicBase);
        return ratio * static_cast<float>(sea_level_m - 1U);
    }

    return static_cast<float>(sea_level_m) +
           (internal - TopographyCodec::kContinentalBase) *
               static_cast<float>(TopographyCodec::kMaxHeightMeters - sea_level_m);
}

// Inverse of the native physical conversion, without uint16 rounding or an
// upper image-format clamp. Relief processes must not use encoded columns as
// distances: their zero is the codec floor, not sea level.
float elevation_to_internal(float elevation_m, uint16_t sea_level_m) {
    if (elevation_m >= 0.0f) {
        return TopographyCodec::kContinentalBase + elevation_m /
            std::max(1.0f, static_cast<float>(TopographyCodec::kMaxHeightMeters - sea_level_m));
    }
    const float ceiling = std::nextafter(TopographyCodec::kContinentalBase,
                                        TopographyCodec::kOceanicBase);
    const float absolute_m = std::max(0.0f, elevation_m + sea_level_m);
    return TopographyCodec::kOceanicBase + std::min(1.0f,
        absolute_m / std::max(1.0f, static_cast<float>(sea_level_m) - 1.0f)) *
        (ceiling - TopographyCodec::kOceanicBase);
}

tectonics::contract::CrustClass classify_crust_class(float total_thickness,
                                                  float continental_thickness) {
    using tectonics::contract::CrustClass;

    if (total_thickness <= FLT_EPSILON) {
        return CrustClass::None;
    }

    const float continental_ratio = continental_thickness / total_thickness;
    if (continental_ratio >= 0.70f) {
        return CrustClass::Continental;
    }
    if (continental_ratio <= 0.30f) {
        return CrustClass::Oceanic;
    }
    return CrustClass::Transitional;
}

tectonics::contract::CrustClass crust_class_from_material(uint8_t material_index) {
    using tectonics::contract::CrustClass;
    return tectonics::material::from_index(material_index) == tectonics::material::Type::Basalt
        ? CrustClass::Oceanic
        : CrustClass::Continental;
}

bool is_continental_material(uint8_t material_index) {
    return crust_class_from_material(material_index) == tectonics::contract::CrustClass::Continental;
}

uint16_t resolve_crust_type_boundary_m(const tectonics::scenario::Scenario& scenario,
                                       uint16_t sea_level_m) {
    if (scenario.crust_type_boundary_m != tectonics::scenario::kAutoCrustTypeBoundaryMeters) {
        return static_cast<uint16_t>(scenario.crust_type_boundary_m);
    }

    const uint32_t threshold =
        static_cast<uint32_t>(sea_level_m) +
        static_cast<uint32_t>(tectonics::scenario::kDefaultCrustTypeBoundaryOffsetMeters);
    return static_cast<uint16_t>(std::min<uint32_t>(TopographyCodec::kMaxHeightMeters, threshold));
}

namespace Phase {
enum : int {
    motion, margins, deformation, rift_displacement, plumes, columns, topology_advance, topology_evolve,
    ownership_publish, provenance, boundary_graph, junctions, deforming_regions, junction_coupling,
    field_stack, creation
};
}
static_assert(Phase::creation + 1 == tectonics::PhaseTimings::count, "one name per phase");

// Adds the wall-clock time of its scope to one phase.
class PhaseTimer {
  public:
    PhaseTimer(tectonics::PhaseTimings& timings, int phase)
        : timings_(timings), phase_(phase), start_(std::chrono::steady_clock::now()) {}
    PhaseTimer(const PhaseTimer&) = delete;
    PhaseTimer& operator=(const PhaseTimer&) = delete;
    ~PhaseTimer() {
        timings_.seconds[static_cast<size_t>(phase_)] +=
            std::chrono::duration<double>(std::chrono::steady_clock::now() - start_).count();
    }

  private:
    tectonics::PhaseTimings& timings_;
    int phase_;
    std::chrono::steady_clock::time_point start_;
};

} // namespace

const char* tectonics::PhaseTimings::name(int phase) {
    static const char* const names[count] = {
        "motion", "margins", "deformation", "rift_displacement", "plumes", "columns",
        "topology_advance", "topology_evolve", "ownership_publish", "provenance", "boundary_graph",
        "junctions", "deforming_regions", "junction_coupling", "field_stack", "creation"};
    return phase >= 0 && phase < count ? names[phase] : "";
}

tectonics::PlateState tectonics::PlateState::fromSeed(uint32_t seed, double spacing_m, double radius_m, double dt_s) {
    // A seed plate's motion: speed 1 cell spacing per update, a rotation about a pole drawn
    // from the plate's own sequence.
    SimpleRandom random(seed);
    const double u = std::min(static_cast<double>(random.next_double()), std::nextafter(1.0, 0.0));
    const double v = std::min(static_cast<double>(random.next_double()), std::nextafter(1.0, 0.0));
    PlateState state;
    state.rotation_rad_s = (spacing_m / (radius_m * dt_s)) * mesh::uniformpoint(u, v);
    return state;
}

namespace {
// The cube's cells (6 N^2); 0 for an invalid N (makegrid throws first).
size_t cubecells(const tectonics::scenario::Scenario& params) {
    return params.sphere_face_cells >= 2 && params.sphere_face_cells <= 1024
        ? 6 * static_cast<size_t>(params.sphere_face_cells) * params.sphere_face_cells : 0;
}
std::shared_ptr<const tectonics::grid::TectonicGrid> makegrid(const tectonics::scenario::Scenario& params) {
    if (params.sphere_face_cells < 2 || params.sphere_face_cells > 1024)
        throw runtime_error("sphere_face_cells must be 2 to 1024");
    if (!std::isfinite(params.radius_km) || !(params.radius_km > 0))
        throw runtime_error("radius_km must be finite and positive");
    return tectonics::grid::TectonicGrid::cube(static_cast<int>(params.sphere_face_cells), params.radius_km);
}
tectonics::dynamics::MantleFlowShape mantleshape(const tectonics::scenario::Scenario& params) {
    tectonics::dynamics::MantleFlowShape shape;
    shape.wavelength_km = params.mantle_flow_wavelength_km;
    shape.octaves = static_cast<int>(params.mantle_flow_octaves);
    shape.wobble_rad = params.mantle_flow_wobble;
    shape.drift_cycles = params.mantle_flow_drift;
    shape.shear_ratio = params.mantle_flow_shear_ratio;
    return shape;
}
} // namespace

lithosphere::lithosphere(const tectonics::scenario::Scenario& params) noexcept(false)
    : scenario(params),
      grid_(makegrid(params)),
      hmap(cubecells(params), 0.0f),
      prev_hmap(cubecells(params), 0.0f),
      display_hmap(cubecells(params), 0.0f),
      prev_display_hmap(cubecells(params), 0.0f),
      initial_hmap(cubecells(params), 0.0f),
      imap(cubecells(params), 0U),
      prev_imap(cubecells(params), 0U),
      amap(cubecells(params), 0U),
      owner_crust_map(cubecells(params), 0.0f),
      owner_material_map(cubecells(params),
                         tectonics::material::to_index(tectonics::material::kDefaultType)),
      iter_count(0), time_origin_step(0), max_plates(params.plate_count), num_plates(0),
      erosion_strength(std::max(0.0f, params.erosion_strength)),
      sea_level_m(TopographyCodec::default_import_sea_level_m()),
      initial_ocean_depth_m(params.initial_ocean_depth_m),
      initial_land_relief_m(params.initial_land_relief_m),
      convergence_map(cubecells(params), 0U),
      divergence_map(cubecells(params), 0U),
      shear_map(cubecells(params), 0U),
      geologic_regime_map(cubecells(params),
                          static_cast<uint8_t>(tectonics::contract::GeologicRegime::Stable)),
      crust_age_myr_map(cubecells(params), 0.0f),
      crust_thickness_map(cubecells(params),
                          0.0f),
      crust_class_map(cubecells(params),
                      static_cast<uint8_t>(tectonics::contract::CrustClass::None)),
      uplift_tendency_map(cubecells(params),
                          0.0f),
      rock_uplift_rate_map(cubecells(params),
                           0.0f),
      shortening_map(cubecells(params), 0.0f),
      shortening_direction_map(cubecells(params), 0.0f),
      subsidence_tendency_map(
          cubecells(params), 0.0f),
      accumulated_strain_map(
          cubecells(params), 0.0f),
      boundary_type_map(cubecells(params),
                        static_cast<uint8_t>(tectonics::contract::BoundaryType::None)),
      boundary_distance_map(
          cubecells(params), 0U),
      boundary_segment_id_map(
          cubecells(params),
          tectonics::contract::kNoBoundaryId),
      nearest_boundary_id_map(
          cubecells(params),
          tectonics::contract::kNoBoundaryId),
      deforming_region_id_map(
          cubecells(params),
          tectonics::contract::kNoDeformingRegionId),
      deforming_region_type_map(
          cubecells(params),
          static_cast<uint8_t>(tectonics::contract::DeformingRegionType::None)),
      deformation_rate_map(
          cubecells(params), 0.0f),
      deformation_velocity_x_map(
          cubecells(params), 0.0f),
      deformation_velocity_y_map(
          cubecells(params), 0.0f),
      junction_rift_influence_map(
          cubecells(params), 0.0f),
      junction_transform_influence_map(
          cubecells(params), 0.0f),
      junction_subduction_influence_map(
          cubecells(params), 0.0f),
      surface_crust_type_map(
          cubecells(params),
          static_cast<uint8_t>(tectonics::contract::CrustClass::Oceanic)),
      material_map(cubecells(params),
                   tectonics::material::to_index(tectonics::material::kDefaultType)),
      continental_lock_map(
          cubecells(params), 0U),
      rift_tension_x_map(cubecells(params), 0.0f),
      rift_tension_y_map(cubecells(params), 0.0f),
      rift_tension_map(cubecells(params), 0.0f),
      _randsource(params.seed), _steps(0) {
    if (!std::isfinite(params.delta_time_myr) || params.delta_time_myr <= 0.0) {
        throw runtime_error("delta_time_myr must be finite and > 0");
    }
    if (!std::isfinite(params.uplift_rate_window_myr) || params.uplift_rate_window_myr <= 0.0) {
        throw runtime_error("uplift_rate_window_myr must be finite and > 0");
    }
    if (!std::isfinite(params.evolution_duration_myr) || params.evolution_duration_myr <= 0.0 ||
        params.evolution_duration_myr / params.delta_time_myr >
            static_cast<double>(std::numeric_limits<int>::max() - 10000)) {
        throw runtime_error("evolution_duration_myr must be finite, positive and fit the update counter");
    }
    if (!std::isfinite(params.mantle_flow_speed_cm_yr) || params.mantle_flow_speed_cm_yr < 0.0) {
        throw runtime_error("mantle_flow_speed_cm_yr must be nonnegative and finite");
    }
    mantle_flow = std::make_shared<const tectonics::dynamics::SphereMantleFlow>(params.seed,
        params.mantle_flow_speed_cm_yr, mantleshape(params), params.radius_km);
    tectonics::ocean_relief::validate({params.ocean_crust_variation_m, params.ocean_noise_wavelength_km,
        params.seamount_crust_m,params.seamount_radius_km,params.seamount_spacing_km,params.seamount_roughness});
    if (!std::isfinite(params.submarine_erosion_ratio) || params.submarine_erosion_ratio<0 || params.submarine_erosion_ratio>1)
        throw runtime_error("Submarine erosion ratio must be in 0..1");
    tectonics::margins::validate(params.margin_processes);
    tectonics::noise_pass::validate(params.experimental_noise);
    tectonics::plumes::validate(params.plume_processes);
    if (params.plume_processes.weakening>0 && !params.topology_processes.progressive_rifting)
        throw runtime_error("Plume weakening acts through progressive rifting; enable it or set weakening to 0");
    if (!std::isfinite(params.mountain_force_range_km) || params.mountain_force_range_km<0)
        throw runtime_error("Inland mountain force requires isostatic columns and a nonnegative range");
    if (!std::isfinite(params.initial_margin_width_m) || params.initial_margin_width_m<0 || params.initial_margin_width_m>4000 ||
        (params.initial_margin_width_m>0 && !params.boundary_conditioned_terrain))
        throw runtime_error("Initial margin width requires boundary-conditioned terrain and 0..4000 m");
    if (!std::isfinite(params.initial_noise_wavelength_km) || params.initial_noise_wavelength_km <= 0 ||
        !std::isfinite(params.initial_noise_warp_km) || params.initial_noise_warp_km < 0)
        throw runtime_error("Initial noise wavelength must be positive and its warp nonnegative (km)");
    if (params.margin_processes.process_variation>0 && !params.margin_processes.enabled)
        throw runtime_error("Process variation requires isostatic columns and active margins");
    if (!std::isfinite(params.continental_fraction) || params.continental_fraction<0 || params.continental_fraction>1 ||
        (params.continental_fraction>0 && (params.continental_fraction<1-static_cast<double>(params.sea_level)+.01-1e-6 ||
            !params.boundary_conditioned_terrain || params.initial_ocean_depth_m<=2500+params.initial_margin_width_m/2)))
        throw runtime_error("Continental fraction must be 0, or in 1 - sea_level + 0.01..1 (a shelf of at least 1 % "
                            "of the surface) with boundary-conditioned terrain and an initial ocean depth over "
                            "2500 m + initial_margin_width_m / 2 (the slope foot)");
    if (!std::isfinite(params.continental_root_kg_m3) || params.continental_root_kg_m3<0 ||
        params.continental_root_kg_m3>=tectonics::columns::mantle_density-tectonics::columns::continental_density)
        throw runtime_error("Continental root load must be in 0..500 kg/m^3 (500 excluded)");
    if (!tectonics::columns::valid_root_max_thickness(params.continental_root_max_thickness_m))
        throw runtime_error("Continental root maximum thickness must be 0 (no cap) or above 35000 m");
    if (!tectonics::columns::valid_root_delamination_width(params.continental_root_max_thickness_m,
            params.continental_root_delamination_width_m))
        throw runtime_error("Continental root delamination width must be 0, or finite and positive with a "
                            "continental root maximum thickness");
    {
        // The slope foot (2500 m + width / 2 deep) must lie above the depth where a seeded
        // continental column has no crust left: 35 km x (500 - root) / (rho_m - rho_w).
        const double crust_zero_depth=tectonics::columns::reference_continent_m
            *(tectonics::columns::mantle_density-tectonics::columns::continental_density-params.continental_root_kg_m3)
            /(tectonics::columns::mantle_density-tectonics::columns::water_density);
        if (2500+params.initial_margin_width_m/2>=crust_zero_depth)
            throw runtime_error("Initial margin width must be below "+std::to_string(std::max(0.0,2*(crust_zero_depth-2500)))
                +" m at this continental root load (slope foot 2500 m + width / 2 above the crust-zero depth "
                +std::to_string(crust_zero_depth)+" m)");
    }
    if ((params.transport_subdivision != 1 && params.transport_subdivision != 2 && params.transport_subdivision != 4) ||
        !std::isfinite(params.deformation_peak_km) || params.deformation_peak_km < 0 ||
        !std::isfinite(params.deformation_variation) || params.deformation_variation < 0 || params.deformation_variation > .5 ||
        !std::isfinite(params.deformation_variation_wavelength_km) || params.deformation_variation_wavelength_km <= 0 ||
        !std::isfinite(params.subduction_coupling) || params.subduction_coupling < 0 || params.subduction_coupling > 1 ||
        !std::isfinite(params.coastal_gpe_strength) || params.coastal_gpe_strength < 0 || params.coastal_gpe_strength > 1)
        throw runtime_error("Invalid transport subdivision or distributed deformation parameters");
    if (!std::isfinite(params.thermal_subsidence_strength) || params.thermal_subsidence_strength < 0.0f) {
        throw runtime_error("thermal_subsidence_strength must be finite and >= 0");
    }
    if (!std::isfinite(params.boundary_force_ramp_cm_yr) || params.boundary_force_ramp_cm_yr < 0)
        throw runtime_error("Boundary force activation speed must be finite and nonnegative");
    if (params.consistent_motion && params.boundary_force_ramp_cm_yr<=0)
        throw runtime_error("Consistent motion requires force balance, isostatic columns and a positive force ramp");
    if (!std::isfinite(params.boundary_memory_myr) || params.boundary_memory_myr<=0)
        throw runtime_error("Boundary memory must be finite and positive");
    if (params.crust_type_boundary_m < tectonics::scenario::kAutoCrustTypeBoundaryMeters ||
        params.crust_type_boundary_m > static_cast<int32_t>(TopographyCodec::kMaxHeightMeters)) {
        throw runtime_error("crust_type_boundary_m is out of range");
    }



    PhaseTimer timer(phase_timings, Phase::creation);
    seedInitialTopography(params.sea_level, params.sea_level_m);
    seedInitialCrustTypes();
    seedMaterialMap();
    createPlates();
}

lithosphere::~lithosphere() noexcept = default;

void lithosphere::growPlates(std::vector<std::vector<uint32_t>>& borders) {
    // "Grow" plates from their origins until surface is fully populated:
    // a random border cell of each plate claims its free cardinal neighbours.
    using tectonics::grid::Direction;
    constexpr std::uint8_t claim_order[] = {Direction::north, Direction::south, Direction::west,
                                            Direction::east};
    uint32_t max_border = 1;
    while (max_border) {
        max_border = 0;
        for (uint32_t i = 0; i < num_plates; ++i) {
            std::vector<uint32_t>& border = borders[i];
            const uint32_t N = (uint32_t)border.size();
            max_border = max_border > N ? max_border : N;

            if (N == 0) {
                continue;
            }
            const uint32_t j = _randsource.next() % N;
            const uint32_t p = border[j];
            for (const std::uint8_t k : claim_order) {
                const uint32_t n = grid_->neighbour(p, k);
                if (imap[n] >= num_plates) {
                    imap[n] = i;
                    border.push_back(n);
                }
            }

            // Overwrite processed point with unprocessed one.
            border[j] = border.back();
            border.pop_back();
        }
    }
}

void lithosphere::createPlates() {
    try {
        const uint32_t map_area = cellCount();
        num_plates = max_plates;
        _steps = 0;
        if (num_plates > map_area) {
            throw runtime_error("plate_count cannot exceed map area");
        }
        std::vector<std::vector<uint32_t>> borders(num_plates);

        // Initialize "Free plate center position" lookup table.
        // This way two plate centers will never be identical.
        std::vector<uint32_t> free_centers(map_area);
        for (uint32_t i = 0; i < map_area; ++i)
            free_centers[i] = i;

        // Select N plate centers from the global map.
        for (uint32_t i = 0; i < num_plates; ++i) {
            // Randomly select an unused plate origin.
            const uint32_t remaining = map_area - i;
            const uint32_t random_index = (uint32_t)_randsource.next() % remaining;
            const uint32_t p = free_centers[random_index];
            borders[i].reserve(8);
            borders[i].push_back(p); // ...and mark it as border.

            // Overwrite used entry with last unused entry in array.
            free_centers[random_index] = free_centers[remaining - 1U];
        }

        std::fill(imap.begin(), imap.end(), 0xFFFFFFFFU);
        growPlates(borders);

        // check all the points of the map are owned
        for (uint32_t i = 0; i < map_area; i++) {
            ASSERT(imap[i] < num_plates, "A point was not assigned to any plate");
        }

        // Each plate draws its motion from a seed of the world's sequence; it records the
        // crust and material of its cells.
        plates.assign(num_plates, {});
        const double dt_s = scenario.delta_time_myr * tectonics::dynamics::kSecondsPerMyr;
        for (uint32_t i = 0; i < num_plates; ++i)
            plates[i] = tectonics::PlateState::fromSeed(_randsource.next(), grid_->cellspacing_m(), grid_->radius_m(), dt_s);
        std::vector<uint32_t> occupied(num_plates, 0U), continental(num_plates, 0U);
        for (uint32_t k = 0; k < map_area; ++k) {
            owner_crust_map[k] = hmap[k];
            owner_material_map[k] = material_map[k];
            if (hmap[k] > 0.0f) {
                ++occupied[imap[k]];
                continental[imap[k]] += is_continental_material(material_map[k]) ? 1U : 0U;
            }
        }
        for (uint32_t i = 0; i < num_plates; ++i) {
            plates[i].creation_continental_share = occupied[i] > 0U
                ? static_cast<float>(continental[i]) / static_cast<float>(occupied[i]) : 0.0f;
            plates[i].continental_share = plates[i].creation_continental_share;
        }
        plate_indices_found.assign(num_plates, 0U);

        if (iter_count == 0U) {
            const double age_horizon_steps =
                std::ceil(260.0 / std::max(0.001, scenario.delta_time_myr));
            iter_count = std::max<uint32_t>(
                num_plates + MAX_BUOYANCY_AGE,
                static_cast<uint32_t>(std::min(
                    age_horizon_steps + static_cast<double>(MAX_BUOYANCY_AGE),
                    static_cast<double>(std::numeric_limits<uint32_t>::max() - 1U))));
        }
        const bool initial_creation = time_origin_step == 0U;
        if (initial_creation) {
            time_origin_step = iter_count;
        }


        updateDerivedMaps();
        // Repartitioning plates is not the creation of new crust. Restarts
        // restore these birth timestamps into the new owners below.
        if (initial_creation) {
            seedInitialCrustAges();
        }
        if (initial_creation) initializeColumns();
        updateDerivedMaps();

    } catch (const exception& e) {
        string msg = "Problem during createPlates: ";
        msg = msg + e.what();
        throw runtime_error(msg.c_str());
    }
}

uint32_t lithosphere::getPlateCount() const throw() {
    return num_plates;
}

uint32_t lithosphere::getProvenancePlateCount() const noexcept {
    return static_cast<uint32_t>(plate_kinematics.size());
}

const uint32_t* lithosphere::getAgeMap() const throw() {
    return amap.data();
}

const float* lithosphere::getCrustAgeMyrMap() const noexcept {
    return crust_age_myr_map.data();
}

const float* lithosphere::getCrustThicknessMap() const noexcept {
    return crust_thickness_map.data();
}

const uint8_t* lithosphere::getCrustClassMap() const noexcept {
    return crust_class_map.data();
}

const uint8_t* lithosphere::getSurfaceCrustTypeMap() const noexcept {
    return surface_crust_type_map.data();
}

const float* lithosphere::getUpliftTendencyMap() const noexcept {
    return uplift_tendency_map.data();
}

const float* lithosphere::getRockUpliftRateMap() const noexcept {
    return rock_uplift_rate_map.data();
}

const float* lithosphere::getShorteningMap() const noexcept {
    return shortening_map.data();
}

const float* lithosphere::getShorteningDirectionMap() const noexcept {
    return shortening_direction_map.data();
}

const float* lithosphere::getSubsidenceTendencyMap() const noexcept {
    return subsidence_tendency_map.data();
}

const float* lithosphere::getAccumulatedStrainMap() const noexcept {
    return accumulated_strain_map.data();
}

const float* lithosphere::getTensionMap() const noexcept {
    return rift_tension_map.data();
}

const uint8_t* lithosphere::getBoundaryTypeMap() const noexcept {
    return boundary_type_map.data();
}

const uint16_t* lithosphere::getBoundaryDistanceMap() const noexcept {
    return boundary_distance_map.data();
}

const uint32_t* lithosphere::getBoundarySegmentIdMap() const noexcept {
    return boundary_segment_id_map.data();
}

const uint32_t* lithosphere::getNearestBoundaryIdMap() const noexcept {
    return nearest_boundary_id_map.data();
}

const tectonics::contract::BoundarySegment* lithosphere::getBoundarySegments() const noexcept {
    return boundary_segments.data();
}

uint32_t lithosphere::getBoundarySegmentCount() const noexcept {
    return static_cast<uint32_t>(boundary_segments.size());
}

const uint32_t* lithosphere::getDeformingRegionIdMap() const noexcept {
    return deforming_region_id_map.data();
}

const uint8_t* lithosphere::getDeformingRegionTypeMap() const noexcept {
    return deforming_region_type_map.data();
}

const float* lithosphere::getDeformationRateMap() const noexcept {
    return deformation_rate_map.data();
}

const float* lithosphere::getDeformationVelocityXMap() const noexcept {
    return deformation_velocity_x_map.data();
}

const float* lithosphere::getDeformationVelocityYMap() const noexcept {
    return deformation_velocity_y_map.data();
}

const tectonics::contract::DeformingRegion* lithosphere::getDeformingRegions() const noexcept {
    return deforming_regions.data();
}

uint32_t lithosphere::getDeformingRegionCount() const noexcept {
    return static_cast<uint32_t>(deforming_regions.size());
}

const tectonics::contract::TectonicJunction* lithosphere::getTectonicJunctions() const noexcept {
    return tectonic_junctions.data();
}

uint32_t lithosphere::getTectonicJunctionCount() const noexcept {
    return static_cast<uint32_t>(tectonic_junctions.size());
}

const uint8_t* lithosphere::getMaterialMap() const noexcept {
    return material_map.data();
}

float* lithosphere::getTopography() const throw() {
    // The C interface hands out a mutable pointer; the engine never reads it back.
    return const_cast<float*>(display_hmap.data());
}

const uint8_t* lithosphere::getConvergenceMap() const noexcept {
    return convergence_map.data();
}

const uint8_t* lithosphere::getDivergenceMap() const noexcept {
    return divergence_map.data();
}

const uint8_t* lithosphere::getShearMap() const noexcept {
    return shear_map.data();
}

const uint8_t* lithosphere::getGeologicRegimeMap() const noexcept {
    return geologic_regime_map.data();
}

const tectonics::contract::PlateKinematics* lithosphere::getPlateKinematics() const noexcept {
    return plate_kinematics.data();
}

double lithosphere::getTimeMyr() const noexcept {
    const uint32_t elapsed_steps =
        iter_count >= time_origin_step ? iter_count - time_origin_step : 0U;
    return tectonics::scenario::steps_to_myr(elapsed_steps, scenario.delta_time_myr);
}

bool lithosphere::isOceanic(float value) const noexcept {
    return TopographyCodec::is_oceanic_internal(value);
}

void lithosphere::resetSimulationState() {
    std::fill(amap.begin(), amap.end(), 0U);
    plume_influence_map.clear();
    std::fill(imap.begin(), imap.end(), 0xFFFFFFFFU);
    std::fill(prev_imap.begin(), prev_imap.end(), 0xFFFFFFFFU);
    prev_hmap = hmap;
    prev_display_hmap = display_hmap;
    std::fill(convergence_map.begin(), convergence_map.end(), static_cast<uint8_t>(0));
    std::fill(divergence_map.begin(), divergence_map.end(), static_cast<uint8_t>(0));
    std::fill(shear_map.begin(), shear_map.end(), static_cast<uint8_t>(0));
    std::fill(geologic_regime_map.begin(), geologic_regime_map.end(),
              static_cast<uint8_t>(tectonics::contract::GeologicRegime::Stable));
    std::fill(crust_age_myr_map.begin(), crust_age_myr_map.end(), 0.0f);
    std::fill(crust_thickness_map.begin(), crust_thickness_map.end(), 0.0f);
    std::fill(crust_class_map.begin(), crust_class_map.end(),
              static_cast<uint8_t>(tectonics::contract::CrustClass::None));
    std::fill(uplift_tendency_map.begin(), uplift_tendency_map.end(), 0.0f);
    std::fill(rock_uplift_rate_map.begin(), rock_uplift_rate_map.end(), 0.0f);
    std::fill(shortening_map.begin(), shortening_map.end(), 0.0f);
    std::fill(shortening_direction_map.begin(), shortening_direction_map.end(), 0.0f);
    std::fill(subsidence_tendency_map.begin(), subsidence_tendency_map.end(), 0.0f);
    std::fill(accumulated_strain_map.begin(), accumulated_strain_map.end(), 0.0f);
    std::fill(boundary_type_map.begin(), boundary_type_map.end(),
              static_cast<uint8_t>(tectonics::contract::BoundaryType::None));
    std::fill(boundary_distance_map.begin(), boundary_distance_map.end(), static_cast<uint16_t>(0));
    std::fill(boundary_segment_id_map.begin(), boundary_segment_id_map.end(),
              tectonics::contract::kNoBoundaryId);
    std::fill(nearest_boundary_id_map.begin(), nearest_boundary_id_map.end(),
              tectonics::contract::kNoBoundaryId);
    std::fill(deforming_region_id_map.begin(), deforming_region_id_map.end(),
              tectonics::contract::kNoDeformingRegionId);
    std::fill(deforming_region_type_map.begin(), deforming_region_type_map.end(),
              static_cast<uint8_t>(tectonics::contract::DeformingRegionType::None));
    std::fill(deformation_rate_map.begin(), deformation_rate_map.end(), 0.0f);
    std::fill(deformation_velocity_x_map.begin(), deformation_velocity_x_map.end(), 0.0f);
    std::fill(deformation_velocity_y_map.begin(), deformation_velocity_y_map.end(), 0.0f);
    std::fill(junction_rift_influence_map.begin(), junction_rift_influence_map.end(), 0.0f);
    std::fill(junction_transform_influence_map.begin(), junction_transform_influence_map.end(), 0.0f);
    std::fill(junction_subduction_influence_map.begin(), junction_subduction_influence_map.end(), 0.0f);
    std::fill(continental_lock_map.begin(), continental_lock_map.end(), static_cast<uint8_t>(0));
    std::fill(rift_tension_x_map.begin(), rift_tension_x_map.end(), 0.0f);
    std::fill(rift_tension_y_map.begin(), rift_tension_y_map.end(), 0.0f);
    std::fill(rift_tension_map.begin(), rift_tension_map.end(), 0.0f);
    plate_kinematics.clear();
    boundary_segments.clear();
    deforming_regions.clear();
    tectonic_junctions.clear();
    next_boundary_segment_id = tectonics::contract::kNoBoundaryId + 1U;

    plates.clear();
    num_plates = 0;

    iter_count = 0;
    time_origin_step = 0;
    _steps = 0;
    finished = false;
    motion_diagnostics = {};
    column_model = {};
}

void lithosphere::updateCrustAgeMyrMap() {
    const uint32_t map_area = cellCount();
    for (uint32_t i = 0; i < map_area; ++i) {
        const uint32_t age_steps = iter_count >= amap[i] ? iter_count - amap[i] : 0U;
        crust_age_myr_map[i] = static_cast<float>(
            tectonics::scenario::steps_to_myr(age_steps, scenario.delta_time_myr));
    }
}

void lithosphere::updateTectonicProvenance() {
    updateCrustAgeMyrMap();
    plate_kinematics.resize(num_plates);
    // Identity, area and centroid of every plate from the ownership map.
    const std::vector<tectonics::grid::Vec3> centroids = grid_->centroids(imap, num_plates);
    for (tectonics::PlateState& state : plates) state.area_m2 = 0;
    for (uint32_t index = 0; index < cellCount(); ++index)
        if (imap[index] < num_plates) plates[imap[index]].area_m2 += grid_->area_m2(index);
    const auto& stable_ids = topology_model.diagnostics().stable_ids;
    plate_rotation_update = plateRotations();
    for (uint32_t plate_index = 0; plate_index < num_plates; ++plate_index) {
        tectonics::PlateState& state = plates[plate_index];
        state.centroid = centroids[plate_index];
        state.stable_id = plate_index < stable_ids.size() ? stable_ids[plate_index] : plate_index;
        // Contract schema 8: omega in rad/s, the centroid as longitude and latitude in degrees,
        // the velocity omega x R c at the centroid (east, north) in cm/yr.
        tectonics::contract::PlateKinematics& kinematics = plate_kinematics[plate_index];
        kinematics.rotation_rad_s[0] = state.rotation_rad_s.x;
        kinematics.rotation_rad_s[1] = state.rotation_rad_s.y;
        kinematics.rotation_rad_s[2] = state.rotation_rad_s.z;
        const auto centre = tectonics::boundary_cube::encodecentroid(state.centroid);
        kinematics.centroid_longitude_deg = centre[0];
        kinematics.centroid_latitude_deg = centre[1];
        const mesh::Vec3 velocity_m_s = mesh::cross(state.rotation_rad_s, grid_->radius_m() * state.centroid);
        const auto velocity = tectonics::boundary_cube::eastnorth(state.centroid, velocity_m_s);
        const double cm_yr = 100.0 * tectonics::dynamics::kSecondsPerMyr / 1.0e6;
        kinematics.velocity_east_cm_yr = static_cast<float>(velocity[0] * cm_yr);
        kinematics.velocity_north_cm_yr = static_cast<float>(velocity[1] * cm_yr);
    }

    // Velocities in cell spacings per update to cm/yr: cellspacing_m / dt.
    const tectonics::provenance::Inputs inputs{
        display_hmap.data(),
        imap.data(),
        num_plates,
        crust_class_map.data(),
        motion_diagnostics.enabled
            ? grid_->cellspacing_m() * 100.0 / (scenario.delta_time_myr * 1e6)
            : 0.0,
        grid_.get(),
        plate_rotation_update.data(),
    };
    tectonics::provenance::compute_maps(inputs, convergence_map.data(), divergence_map.data(),
                                     shear_map.data(), geologic_regime_map.data());
}

void lithosphere::updateBoundaryDistanceMap() {
    const uint32_t map_area = cellCount();
    // Farther than any cell can be from a boundary (half a great circle is 2 N steps).
    const uint16_t fallback_distance = static_cast<uint16_t>(
        std::min<uint32_t>(std::numeric_limits<uint16_t>::max(), 6U * static_cast<uint32_t>(grid_->facecells())));

    std::fill(boundary_distance_map.begin(), boundary_distance_map.end(), fallback_distance);
    std::fill(nearest_boundary_id_map.begin(), nearest_boundary_id_map.end(),
              tectonics::contract::kNoBoundaryId);

    std::deque<uint32_t> frontier;
    frontier.clear();
    for (uint32_t index = 0; index < map_area; ++index) {
        if (boundary_segment_id_map[index] != tectonics::contract::kNoBoundaryId) {
            boundary_distance_map[index] = 0U;
            nearest_boundary_id_map[index] = boundary_segment_id_map[index];
            frontier.push_back(index);
        }
    }

    if (frontier.empty()) {
        return;
    }

    while (!frontier.empty()) {
        const uint32_t index = frontier.front();
        frontier.pop_front();

        const auto neighbors = grid_->neighbours(index);
        const uint16_t next_distance = static_cast<uint16_t>(boundary_distance_map[index] + 1U);
        const uint32_t boundary_id = nearest_boundary_id_map[index];
        for (uint32_t neighbor : neighbors) {
            if (next_distance < boundary_distance_map[neighbor] ||
                (next_distance == boundary_distance_map[neighbor] &&
                 boundary_id != tectonics::contract::kNoBoundaryId &&
                 (nearest_boundary_id_map[neighbor] == tectonics::contract::kNoBoundaryId ||
                  boundary_id < nearest_boundary_id_map[neighbor]))) {
                boundary_distance_map[neighbor] = next_distance;
                nearest_boundary_id_map[neighbor] = boundary_id;
                frontier.push_back(neighbor);
            }
        }
    }
}

void lithosphere::updateBoundaryGraphMaps() {
    const uint32_t map_area = cellCount();
    const std::vector<uint8_t> raw_geologic_regime = geologic_regime_map;

    const tectonics::boundary_graph::Inputs inputs{
        imap.data(),
        convergence_map.data(),
        divergence_map.data(),
        shear_map.data(),
        raw_geologic_regime.data(),
        crust_class_map.data(),
        boundary_segment_id_map.data(),
        nearest_boundary_id_map.data(),
        boundary_segments.data(),
        static_cast<uint32_t>(boundary_segments.size()),
        next_boundary_segment_id,
        scenario.delta_time_myr,
        grid_.get(),
    };
    tectonics::boundary_graph::Outputs outputs = tectonics::boundary_graph::build(inputs);
    boundary_segments = std::move(outputs.segments);
    boundary_segment_id_map = std::move(outputs.boundary_segment_id);
    boundary_type_map = std::move(outputs.boundary_type);
    geologic_regime_map = std::move(outputs.boundary_regime);
    next_boundary_segment_id = outputs.next_boundary_id;

    updateBoundaryDistanceMap();

    std::vector<const tectonics::contract::BoundarySegment*> segments_by_id(
        next_boundary_segment_id, nullptr);
    for (const tectonics::contract::BoundarySegment& segment : boundary_segments) {
        if (segment.id < segments_by_id.size()) {
            segments_by_id[segment.id] = &segment;
        }
    }

    for (uint32_t index = 0; index < map_area; ++index) {
        if (boundary_segment_id_map[index] != tectonics::contract::kNoBoundaryId) {
            continue;
        }

        const uint32_t nearest_boundary_id = nearest_boundary_id_map[index];
        const auto raw_regime =
            static_cast<tectonics::contract::GeologicRegime>(raw_geologic_regime[index]);
        if (nearest_boundary_id == tectonics::contract::kNoBoundaryId ||
            nearest_boundary_id >= segments_by_id.size() ||
            segments_by_id[nearest_boundary_id] == nullptr) {
            geologic_regime_map[index] =
                static_cast<uint8_t>(raw_regime == tectonics::contract::GeologicRegime::PassiveMargin
                                         ? raw_regime
                                         : tectonics::contract::GeologicRegime::Stable);
            continue;
        }

        const tectonics::contract::BoundarySegment& segment =
            *segments_by_id[nearest_boundary_id];
        const uint16_t distance = boundary_distance_map[index];
        auto resolved_regime = tectonics::contract::GeologicRegime::Stable;
        if (distance <= 1U) {
            switch (segment.boundary_type) {
            case tectonics::contract::BoundaryType::Convergent:
                if (raw_regime == tectonics::contract::GeologicRegime::TrenchAdjacent ||
                    segment.geologic_regime == tectonics::contract::GeologicRegime::TrenchAdjacent) {
                    resolved_regime = tectonics::contract::GeologicRegime::TrenchAdjacent;
                } else {
                    resolved_regime = segment.geologic_regime;
                }
                break;
            case tectonics::contract::BoundaryType::Divergent:
                resolved_regime = segment.geologic_regime;
                break;
            case tectonics::contract::BoundaryType::Transform:
                resolved_regime = tectonics::contract::GeologicRegime::Transform;
                break;
            case tectonics::contract::BoundaryType::PassiveMargin:
                resolved_regime = tectonics::contract::GeologicRegime::PassiveMargin;
                break;
            case tectonics::contract::BoundaryType::None:
            default:
                resolved_regime = tectonics::contract::GeologicRegime::Stable;
                break;
            }
        } else if (distance <= 3U &&
                   segment.boundary_type == tectonics::contract::BoundaryType::PassiveMargin) {
            resolved_regime = tectonics::contract::GeologicRegime::PassiveMargin;
        } else if (distance <= 2U && raw_regime != tectonics::contract::GeologicRegime::Stable) {
            resolved_regime = raw_regime;
        }

        // Nearby collision belts cannot change an oceanic cell into continental crust.
        if (resolved_regime == tectonics::contract::GeologicRegime::ContinentCollision &&
            crust_class_map[index] == static_cast<uint8_t>(tectonics::contract::CrustClass::Oceanic)) {
            resolved_regime = tectonics::contract::GeologicRegime::ConvergentArc;
        }
        geologic_regime_map[index] = static_cast<uint8_t>(resolved_regime);
    }
}

void lithosphere::updateTectonicJunctions() {
    const tectonics::tectonic_junctions::Inputs inputs{
        imap.data(),
        boundary_segment_id_map.data(),
        boundary_segments.data(),
        static_cast<uint32_t>(boundary_segments.size()),
        static_cast<uint32_t>(plate_rotation_update.size()),
        grid_.get(),
        plate_rotation_update.data(),
    };
    tectonics::tectonic_junctions::Outputs outputs = tectonics::tectonic_junctions::build(inputs);
    tectonic_junctions = std::move(outputs.junctions);
}

void lithosphere::updateTectonicJunctionCouplingMaps() {
    const uint32_t map_area = cellCount();
    const tectonics::tectonic_junction_coupling::Inputs inputs{
        tectonic_junctions.data(),
        static_cast<uint32_t>(tectonic_junctions.size()),
        grid_.get(),
    };
    tectonics::tectonic_junction_coupling::Outputs outputs =
        tectonics::tectonic_junction_coupling::build(inputs);
    junction_rift_influence_map = std::move(outputs.rift_influence);
    junction_transform_influence_map = std::move(outputs.transform_influence);
    junction_subduction_influence_map = std::move(outputs.subduction_influence);

    if (outputs.rift_tension_x.size() != map_area || outputs.rift_tension_y.size() != map_area ||
        rift_tension_x_map.size() != map_area || rift_tension_y_map.size() != map_area) {
        return;
    }

    constexpr float kJunctionRiftTensionGain = 18.0f * 0.35f;
    constexpr float kMaxTension = 100.0f * 2.5f;
    for (uint32_t index = 0; index < map_area; ++index) {
        rift_tension_x_map[index] += outputs.rift_tension_x[index] * kJunctionRiftTensionGain;
        rift_tension_y_map[index] += outputs.rift_tension_y[index] * kJunctionRiftTensionGain;
        const float magnitude =
            std::sqrt(rift_tension_x_map[index] * rift_tension_x_map[index] +
                      rift_tension_y_map[index] * rift_tension_y_map[index]);
        if (magnitude > kMaxTension) {
            const float scale = kMaxTension / magnitude;
            rift_tension_x_map[index] *= scale;
            rift_tension_y_map[index] *= scale;
        }
    }
    updateTensionMagnitudeMap();
}

void lithosphere::updateDeformingRegionMaps() {
    const tectonics::deforming_regions::Inputs inputs{
        imap.data(),
        crust_class_map.data(),
        convergence_map.data(),
        divergence_map.data(),
        shear_map.data(),
        boundary_type_map.data(),
        geologic_regime_map.data(),
        boundary_distance_map.data(),
        nearest_boundary_id_map.data(),
        boundary_segments.data(),
        static_cast<uint32_t>(boundary_segments.size()),
        static_cast<uint32_t>(plate_rotation_update.size()),
        grid_.get(),
        plate_rotation_update.data(),
    };
    tectonics::deforming_regions::Outputs outputs =
        tectonics::deforming_regions::build(inputs);
    deforming_regions = std::move(outputs.regions);
    deforming_region_id_map = std::move(outputs.deforming_region_id);
    deforming_region_type_map = std::move(outputs.deforming_region_type);
    deformation_rate_map = std::move(outputs.deformation_rate);
    deformation_velocity_x_map = std::move(outputs.deformation_velocity_x);
    deformation_velocity_y_map = std::move(outputs.deformation_velocity_y);
}

void lithosphere::updateFieldStackMaps() {
    if (num_plates == 0) {
        return;
    }

    const uint32_t map_area = cellCount();
    std::vector<float> continental_thickness(map_area, 0.0f);
    std::vector<float> oceanic_thickness(map_area, 0.0f);

    std::fill(crust_thickness_map.begin(), crust_thickness_map.end(), 0.0f);

    // The owner's recorded crust, split by its recorded material; each plate's continental
    // share counts its cells with crust (above FLT_EPSILON) of continental material.
    std::vector<uint32_t> occupied(num_plates, 0U), continental_cells(num_plates, 0U);
    for (uint32_t index = 0; index < map_area; ++index) {
        const float crust = owner_crust_map[index];
        const uint32_t owner = imap[index];
        const bool continental = is_continental_material(owner_material_map[index]);
        if (owner < num_plates && crust > FLT_EPSILON) {
            ++occupied[owner];
            continental_cells[owner] += continental ? 1U : 0U;
        }
        if (crust <= 2.0f * FLT_EPSILON) {
            continue;
        }
        crust_thickness_map[index] += crust;
        if (continental) {
            continental_thickness[index] += crust;
        } else {
            oceanic_thickness[index] += crust;
        }
    }
    for (uint32_t plate_index = 0; plate_index < num_plates; ++plate_index)
        plates[plate_index].continental_share = occupied[plate_index] > 0U
            ? static_cast<float>(continental_cells[plate_index]) / static_cast<float>(occupied[plate_index])
            : plates[plate_index].creation_continental_share;

    for (uint32_t index = 0; index < map_area; ++index) {
        crust_thickness_map[index] =
            std::max(crust_thickness_map[index], std::max(0.0f, display_hmap[index]));
        crust_class_map[index] = static_cast<uint8_t>(classify_crust_class(
            crust_thickness_map[index], continental_thickness[index]));
        if (imap[index] >= num_plates || hmap[index] <= 2.0f * FLT_EPSILON) {
            surface_crust_type_map[index] =
                static_cast<uint8_t>(tectonics::contract::CrustClass::None);
        } else {
            surface_crust_type_map[index] =
                static_cast<uint8_t>(crust_class_from_material(material_map[index]));
        }
    }

    {
        PhaseTimer timer(phase_timings, Phase::provenance);
        updateTectonicProvenance();
    }
    {
        PhaseTimer timer(phase_timings, Phase::boundary_graph);
        updateBoundaryGraphMaps();
    }
    {
        PhaseTimer timer(phase_timings, Phase::junctions);
        updateTectonicJunctions();
    }
    {
        PhaseTimer timer(phase_timings, Phase::deforming_regions);
        updateDeformingRegionMaps();
    }
    {
        PhaseTimer timer(phase_timings, Phase::junction_coupling);
        updateTectonicJunctionCouplingMaps();
    }
    PhaseTimer field_timer(phase_timings, Phase::field_stack);

    constexpr float kSurfaceDeltaScale = 0.08f;
    constexpr float kYoungCrustMyr = 35.0f;
    constexpr float kOldCrustMyr = 80.0f;

    for (uint32_t index = 0; index < map_area; ++index) {
        const float total_thickness = crust_thickness_map[index];
        const float continental_ratio =
            total_thickness > FLT_EPSILON ? continental_thickness[index] / total_thickness : 0.0f;
        const float oceanic_ratio =
            total_thickness > FLT_EPSILON ? oceanic_thickness[index] / total_thickness : 0.0f;

        const float convergence = static_cast<float>(convergence_map[index]) / 100.0f;
        const float divergence = static_cast<float>(divergence_map[index]) / 100.0f;
        const float shear = static_cast<float>(shear_map[index]) / 100.0f;
        const auto regime =
            static_cast<tectonics::contract::GeologicRegime>(geologic_regime_map[index]);
        const auto boundary_type =
            static_cast<tectonics::contract::BoundaryType>(boundary_type_map[index]);
        const auto deforming_region_type = static_cast<tectonics::contract::DeformingRegionType>(
            deforming_region_type_map[index]);
        const float deformation_rate = deformation_rate_map[index];
        const float deformation_velocity = std::sqrt(
            deformation_velocity_x_map[index] * deformation_velocity_x_map[index] +
            deformation_velocity_y_map[index] * deformation_velocity_y_map[index]);
        const float junction_rift = junction_rift_influence_map[index];
        const float junction_transform = junction_transform_influence_map[index];
        const float junction_subduction = junction_subduction_influence_map[index];

        const float delta_surface = display_hmap[index] - prev_display_hmap[index];
        const float positive_surface = clamp_unit(delta_surface / kSurfaceDeltaScale);
        const float negative_surface = clamp_unit(-delta_surface / kSurfaceDeltaScale);
        const float age_old = clamp_unit(crust_age_myr_map[index] / kOldCrustMyr);
        const float age_young = 1.0f - clamp_unit(crust_age_myr_map[index] / kYoungCrustMyr);
        if (deforming_region_type ==
            tectonics::contract::DeformingRegionType::ContinentalRift) {
            crust_thickness_map[index] =
                std::max(display_hmap[index],
                         crust_thickness_map[index] * (1.0f - 0.12f * deformation_rate));
        } else if (deforming_region_type ==
                   tectonics::contract::DeformingRegionType::DiffuseCollision) {
            crust_thickness_map[index] +=
                deformation_rate * (0.12f + 0.18f * continental_ratio);
        }
        if (junction_rift > 0.0f &&
            (static_cast<tectonics::contract::CrustClass>(crust_class_map[index]) ==
                 tectonics::contract::CrustClass::Continental ||
             static_cast<tectonics::contract::CrustClass>(crust_class_map[index]) ==
                 tectonics::contract::CrustClass::Transitional)) {
            crust_thickness_map[index] =
                std::max(display_hmap[index], crust_thickness_map[index] * (1.0f - 0.045f * junction_rift));
        }

        float uplift =
            0.50f * convergence + 0.10f * shear + 0.15f * positive_surface +
            0.10f * continental_ratio + 0.08f * age_young;
        float subsidence =
            0.24f * divergence + 0.18f * negative_surface +
            0.18f * (oceanic_ratio * age_old) + 0.10f * (1.0f - continental_ratio) +
            0.06f * (boundary_type == tectonics::contract::BoundaryType::PassiveMargin ? 1.0f : 0.0f);

        switch (regime) {
        case tectonics::contract::GeologicRegime::ConvergentArc:
            uplift += 0.24f;
            break;
        case tectonics::contract::GeologicRegime::ContinentCollision:
            uplift += 0.30f;
            break;
        case tectonics::contract::GeologicRegime::TrenchAdjacent:
            uplift += 0.10f;
            subsidence += 0.16f;
            break;
        case tectonics::contract::GeologicRegime::DivergentRift:
            subsidence += 0.18f;
            break;
        case tectonics::contract::GeologicRegime::PassiveMargin:
            subsidence += 0.20f;
            break;
        case tectonics::contract::GeologicRegime::MidOceanRidge:
            uplift += 0.12f;
            subsidence -= 0.08f;
            break;
        case tectonics::contract::GeologicRegime::Transform:
            uplift += 0.05f;
            break;
        case tectonics::contract::GeologicRegime::Stable:
        default:
            break;
        }

        switch (deforming_region_type) {
        case tectonics::contract::DeformingRegionType::ContinentalRift:
            subsidence += 0.28f * deformation_rate;
            subsidence += 0.05f * deformation_velocity;
            break;
        case tectonics::contract::DeformingRegionType::DiffuseCollision:
            uplift += 0.30f * deformation_rate;
            uplift += 0.04f * deformation_velocity;
            subsidence -= 0.05f * deformation_rate;
            break;
        case tectonics::contract::DeformingRegionType::None:
        default:
            break;
        }
        subsidence += 0.20f * junction_rift;
        uplift += 0.05f * junction_rift;
        uplift += 0.04f * junction_transform;
        subsidence += 0.02f * junction_transform;
        uplift += 0.09f * junction_subduction;
        subsidence += 0.12f * junction_subduction;

        uplift_tendency_map[index] = clamp_unit(uplift);
        subsidence_tendency_map[index] = clamp_unit(subsidence);

        float strain_input =
            0.26f * convergence + 0.22f * divergence + 0.34f * shear +
            0.10f * std::max(positive_surface, negative_surface);
        switch (boundary_type) {
        case tectonics::contract::BoundaryType::Convergent:
            strain_input += 0.10f;
            break;
        case tectonics::contract::BoundaryType::Divergent:
            strain_input += 0.08f;
            break;
        case tectonics::contract::BoundaryType::Transform:
            strain_input += 0.12f;
            break;
        case tectonics::contract::BoundaryType::PassiveMargin:
            strain_input += 0.03f;
            break;
        case tectonics::contract::BoundaryType::None:
        default:
            break;
        }
        switch (deforming_region_type) {
        case tectonics::contract::DeformingRegionType::ContinentalRift:
            strain_input += 0.14f * deformation_rate;
            break;
        case tectonics::contract::DeformingRegionType::DiffuseCollision:
            strain_input += 0.18f * deformation_rate;
            break;
        case tectonics::contract::DeformingRegionType::None:
        default:
            break;
        }
        accumulated_strain_map[index] =
            clamp_unit(accumulated_strain_map[index] * 0.965f + clamp_unit(strain_input) * 0.22f);
    }
}

void lithosphere::updateDerivedMaps() {
    if (num_plates == 0) updateTectonicProvenance();
    updateFieldStackMaps();
}

void lithosphere::initializeHeightMapFromFloat(const float* heightmap) {
    if (heightmap == nullptr) {
        throw invalid_argument("Height map cannot be null");
    }

    const uint32_t map_area = cellCount();
    for (uint32_t i = 0; i < map_area; ++i) {
        hmap[i] = heightmap[i];
    }
    prev_hmap = hmap;
    display_hmap = hmap;
    prev_display_hmap = display_hmap;
    initial_hmap = hmap;
}

void lithosphere::initializeHeightMapFromMetric(const uint16_t* heightmap_m, uint16_t new_sea_level_m) {
    if (heightmap_m == nullptr) {
        throw invalid_argument("Metric height map cannot be null");
    }

    sea_level_m = new_sea_level_m;
    const uint32_t map_area = cellCount();
    for (uint32_t i = 0; i < map_area; ++i) {
        hmap[i] = TopographyCodec::meters_to_internal(heightmap_m[i], sea_level_m);
    }
    prev_hmap = hmap;
    display_hmap = hmap;
    prev_display_hmap = display_hmap;
    initial_hmap = hmap;
}

void lithosphere::seedInitialTopography(float ocean_coverage, int32_t sea_level_m_override) {
    const uint32_t map_area = cellCount();
    const bool has_override = sea_level_m_override != TopographyCodec::kNoSeaLevelOverride;
    if (sea_level_m_override < TopographyCodec::kNoSeaLevelOverride ||
        sea_level_m_override > TopographyCodec::kMaxHeightMeters) {
        throw invalid_argument("sea_level_m is out of range");
    }
    sea_level_m = has_override ? static_cast<uint16_t>(sea_level_m_override)
                               : TopographyCodec::default_import_sea_level_m();
    if (initial_ocean_depth_m == 0U || initial_ocean_depth_m > sea_level_m ||
        initial_land_relief_m > TopographyCodec::kMaxHeightMeters - sea_level_m) {
        throw invalid_argument("Initial depth/relief must fit around the sea-level datum");
    }
    SimpleRandom topo_rng(_randsource);

    // The torus engine's fields on the sphere (seeding_noise): octaves, persistence and
    // frequency, the last with the coordinate factor each field was read at.
    const double radius_km = grid_->radius_m() / 1000;
    const double unit_km = scenario.initial_noise_wavelength_km;
    const auto field = [&](int octaves, double persistence, double frequency) {
        return seeding_noise(topo_rng.next(), unit_km, radius_km, octaves, persistence, frequency);
    };
    const auto warp_x_noise = field(4, 0.58, 0.45 * 1.6);
    const auto warp_y_noise = field(4, 0.58, 0.45 * 1.6);
    const auto continent_macro_noise = field(5, 0.57, 0.16 * 0.95);
    const auto continent_detail_noise = field(6, 0.59, 0.42 * 1.85 * 1.15);
    const auto coast_breakup_noise = field(6, 0.57, 0.62 * 2.2 * 1.55);
    const auto shelf_noise = field(6, 0.61, 0.85 * 2.9 * 1.4);
    const auto sea_base_noise = field(7, 0.63, 1.15 * 3.9 * 1.15);
    const auto sea_detail_noise = field(7, 0.61, 1.95 * 5.6 * 2.3);
    const auto sea_micro_noise = field(5, 0.48, 3.95 * 10.2 * 4.8);
    const auto land_noise = field(6, 0.56, 0.95 * 3.6 * 1.55);
    const auto foothill_noise = field(6, 0.53, 1.55 * 4.7 * 2.15);
    const auto mountain_noise = field(5, 0.47, 2.45 * 6.3 * 3.05);
    const auto ridge_noise = field(5, 0.50, 2.6 * 7.4 * 2.1);
    const auto trench_noise = field(5, 0.48, 3.8 * 9.8 * 2.35);
    // The warp's third component (the torus had two axes).
    const auto warp_z_noise = field(4, 0.58, 0.45 * 1.6);
    const auto at = [](const tectonics::sphere_noise::Octaves& noise, const mesh::Vec3& q) {
        return static_cast<float>(noise(q));
    };

    std::vector<mesh::Vec3> warped(map_area);
    std::vector<float> continentality(map_area);
    std::vector<float> raw_height(map_area);
    std::vector<uint8_t> is_land(map_area, 0);

    // Each cell's centre moved by the warp, then the fields read there.
    for (uint32_t index = 0; index < map_area; ++index) {
        const mesh::Vec3 q = tectonics::sphere_noise::warp(grid_->position(index), warp_x_noise, warp_y_noise,
                                                           warp_z_noise, scenario.initial_noise_warp_km, radius_km);
        warped[index] = q;
        const float macro = at(continent_macro_noise, q);
        const float detail = at(continent_detail_noise, q);
        const float coast_breakup = at(coast_breakup_noise, q);
        continentality[index] = macro * 0.70f + detail * 0.22f + coast_breakup * 0.08f;
    }

    const float continent_threshold = grid_->quantile(continentality, ocean_coverage);
    // A continental fraction above the land share (scenario.continental_fraction):
    // the ocean side between the land and the continental edge is shelf and slope,
    // and the deep floor's profile starts at the continental edge.
    const bool shelf_mode = scenario.continental_fraction > 0;
    const float shelf_edge = shelf_mode
        ? grid_->quantile(continentality, static_cast<float>(1 - scenario.continental_fraction))
        : continent_threshold;
    const auto is_shelf = [&](uint32_t index) {
        return shelf_mode && continentality[index] < continent_threshold && continentality[index] >= shelf_edge;
    };
    // The deep floor's profile is measured from the continental edge.
    const auto ocean_edge = [&](uint32_t index) {
        return shelf_mode && continentality[index] < shelf_edge ? shelf_edge : continent_threshold;
    };
    float raw_land_min = FLT_MAX;
    float raw_land_max = -FLT_MAX;
    float raw_ocean_min = FLT_MAX;
    float raw_ocean_max = -FLT_MAX;

    for (uint32_t index = 0; index < map_area; ++index) {
        {
            const mesh::Vec3& q = warped[index];
            const float detail = at(continent_detail_noise, q);
            const float shelf = at(shelf_noise, q);
            const float sea_base = at(sea_base_noise, q);
            const float sea_detail = at(sea_detail_noise, q);
            const float sea_micro = at(sea_micro_noise, q);
            const float land_detail = at(land_noise, q);
            const float foothills = at(foothill_noise, q);
            const float mountains = 1.0f - std::abs(at(mountain_noise, q));
            const float ridge = 1.0f - std::abs(at(ridge_noise, q));
            const float trench = 1.0f - std::abs(at(trench_noise, q));
            const float distance_to_coast = continentality[index] - ocean_edge(index);
            const float coastal_transition =
                smoothstep(0.010f, 0.110f, std::abs(distance_to_coast));
            const float shoreline_profile =
                0.09f * std::tanh(distance_to_coast * 18.0f) + 0.07f * detail + 0.05f * shelf +
                0.04f * land_detail + 0.06f * ridge - 0.06f * trench + 0.03f * foothills;

            if (distance_to_coast >= 0.0f) {
                is_land[index] = 1U;
                const float inland_bias = logistic(distance_to_coast + 0.05f * land_detail, 8.0f);
                const float coastal_bias = 1.0f - inland_bias;
                const float inland_profile =
                    0.34f * inland_bias - 0.08f * coastal_bias + 0.18f * detail +
                    0.19f * land_detail + 0.15f * foothills + 0.17f * mountains +
                    0.10f * ridge - 0.05f * trench + 0.04f * shelf;
                raw_height[index] =
                    interpolate(shoreline_profile, inland_profile, coastal_transition);
                raw_land_min = std::min(raw_land_min, raw_height[index]);
                raw_land_max = std::max(raw_land_max, raw_height[index]);
            } else {
                const float shelf_weight = logistic(distance_to_coast + 0.10f * shelf + 0.05f, 11.0f);
                const float abyss_weight = 1.0f - shelf_weight;
                const float ocean_profile =
                    -0.62f * abyss_weight - 0.16f * shelf_weight + 0.17f * sea_base +
                    0.15f * sea_detail + 0.09f * sea_micro + 0.08f * shelf + 0.17f * ridge -
                    0.27f * trench + 0.04f * detail;
                raw_height[index] =
                    interpolate(shoreline_profile, ocean_profile, coastal_transition);
                if (!is_shelf(index)) {
                    raw_ocean_min = std::min(raw_ocean_min, raw_height[index]);
                    raw_ocean_max = std::max(raw_ocean_max, raw_height[index]);
                }
            }
        }
    }

    // Scale physical offsets, then encode once. Both profiles meet the sea datum
    // at the continentality threshold; independent rescaling must not create a
    // multi-kilometre jump between the nearest ocean and land samples.
    std::vector<float> profile(map_area);
    float land_peak = 0.0f;
    float ocean_peak = 0.0f;
    for (uint32_t i = 0; i < map_area; ++i) {
        if (is_shelf(i)) continue;
        const float coast_weight = smoothstep(
            0.0f, 0.110f, std::abs(continentality[i] - ocean_edge(i)));
        const float amplitude = is_land[i]
            ? remap_to_range(raw_height[i], raw_land_min, raw_land_max, 0.0f, 1.0f)
            : remap_to_range(raw_height[i], raw_ocean_min, raw_ocean_max, 1.0f, 0.0f);
        profile[i] = coast_weight * amplitude;
        float& peak = is_land[i] ? land_peak : ocean_peak;
        peak = std::max(peak, profile[i]);
    }
    // Shelf from the coast to -200 m, then a slope crossing -2500 m at the
    // continental edge over the last continentality of about 0.01 before it.
    // Columns are continental above -2500 m with a taper of the initial margin
    // width centred there (tectonic_initialization.cpp), so the slope continues
    // at the same gradient to -(2500 + margin / 2) m and the deep floor starts
    // there: the taper's area above and below the edge balance.
    constexpr double shelf_depth_m = 200, slope_edge_m = 2500, slope_width = 0.010;
    const double slope_foot_m = slope_edge_m + scenario.initial_margin_width_m / 2;
    const double shelf_width = static_cast<double>(continent_threshold) - shelf_edge;
    const double slope = std::min(slope_width, shelf_width / 2);
    const auto slope_depth = [&](double d) {
        return shelf_depth_m + (slope_edge_m - shelf_depth_m) * (d - (shelf_width - slope)) / slope;
    };
    std::vector<uint16_t> seeded_heightmap(map_area);
    for (uint32_t i = 0; i < map_area; ++i) {
        const float peak = is_land[i] ? land_peak : ocean_peak;
        const float fraction = peak > FLT_EPSILON ? profile[i] / peak : 0.0f;
        int offset_m = is_land[i]
            ? static_cast<int>(std::lround(initial_land_relief_m * fraction))
            : -std::max(1, static_cast<int>(std::lround(initial_ocean_depth_m * fraction)));
        if (shelf_mode && !is_land[i]) {
            if (is_shelf(i)) {
                const double d = static_cast<double>(continent_threshold) - continentality[i];
                const double depth = d <= shelf_width - slope
                    ? shelf_depth_m * d / (shelf_width - slope)
                    : slope_depth(d);
                offset_m = -std::max(1, static_cast<int>(std::lround(depth)));
            } else {
                // The slope past the continental edge down to its foot, then the deep
                // floor below the foot, never at it (margin 0: the foot is the edge).
                const double depth = slope_depth(static_cast<double>(continent_threshold) - continentality[i]);
                offset_m = slope_foot_m > slope_edge_m && depth < slope_foot_m
                    ? -static_cast<int>(std::lround(depth))
                    : -static_cast<int>(std::lround(slope_foot_m)) - std::max(1, static_cast<int>(std::lround(
                        (initial_ocean_depth_m - slope_foot_m) * fraction)));
            }
        }
        seeded_heightmap[i] = static_cast<uint16_t>(sea_level_m + offset_m);
    }
    initializeHeightMapFromMetric(seeded_heightmap.data(), sea_level_m);
}

void lithosphere::fillEnclosedInitialOceans() {
    using tectonics::contract::CrustClass;

    const uint32_t map_area = cellCount();
    if (map_area == 0U) {
        return;
    }

    std::vector<uint8_t> visited(map_area, 0U);
    std::vector<uint32_t> largest_component;
    largest_component.reserve(map_area / 2U);

    for (uint32_t index = 0; index < map_area; ++index) {
        if (visited[index] != 0U ||
            surface_crust_type_map[index] != static_cast<uint8_t>(CrustClass::Oceanic)) {
            continue;
        }

        std::vector<uint32_t> component;
        component.reserve(256U);
        std::deque<uint32_t> frontier;
        frontier.push_back(index);
        visited[index] = 1U;

        while (!frontier.empty()) {
            const uint32_t current = frontier.front();
            frontier.pop_front();
            component.push_back(current);

            for (const uint32_t neighbor : grid_->neighbours(current)) {
                if (visited[neighbor] != 0U ||
                    surface_crust_type_map[neighbor] !=
                        static_cast<uint8_t>(CrustClass::Oceanic)) {
                    continue;
                }
                visited[neighbor] = 1U;
                frontier.push_back(neighbor);
            }
        }

        if (component.size() > largest_component.size()) {
            largest_component.swap(component);
        }
    }

    std::vector<uint8_t> keep_ocean(map_area, 0U);
    for (uint32_t index : largest_component) {
        keep_ocean[index] = 1U;
    }

    for (uint32_t index = 0; index < map_area; ++index) {
        if (surface_crust_type_map[index] == static_cast<uint8_t>(CrustClass::Oceanic) &&
            keep_ocean[index] == 0U) {
            surface_crust_type_map[index] = static_cast<uint8_t>(CrustClass::Continental);
        }
    }
}

void lithosphere::seedInitialCrustTypes() {
    using tectonics::contract::CrustClass;

    const uint32_t map_area = cellCount();
    const uint16_t threshold_m = resolve_crust_type_boundary_m(scenario, sea_level_m);
    for (uint32_t index = 0; index < map_area; ++index) {
        const uint16_t height_m = TopographyCodec::internal_to_meters(hmap[index], sea_level_m);
        surface_crust_type_map[index] =
            height_m >= threshold_m ? static_cast<uint8_t>(CrustClass::Continental)
                                    : static_cast<uint8_t>(CrustClass::Oceanic);
    }

    fillEnclosedInitialOceans();
}

void lithosphere::seedMaterialMap() {
    const uint32_t map_area = cellCount();
    if (map_area == 0U) {
        material_map.clear();
        return;
    }

    SimpleRandom material_rng(static_cast<uint32_t>(scenario.seed) ^ 0x6A09E667U);
    // The torus engine's fields on the sphere (seeding_noise); the warp 2.4 / 18 of the relief's.
    const double radius_km = grid_->radius_m() / 1000;
    const double unit_km = scenario.initial_noise_wavelength_km;
    const auto field = [&](int octaves, double persistence, double frequency) {
        return seeding_noise(material_rng.next(), unit_km, radius_km, octaves, persistence, frequency);
    };
    const auto warp_x_noise = field(4, 0.58, 0.28 * 1.4);
    const auto warp_y_noise = field(4, 0.58, 0.28 * 1.4);
    const auto material_macro_noise = field(5, 0.57, 0.22 * 1.2);
    const auto material_detail_noise = field(6, 0.61, 0.75 * 2.4 * 1.8);
    const auto warp_z_noise = field(4, 0.58, 0.28 * 1.4);
    const double warp_km = scenario.initial_noise_warp_km * 2.4 / 18;

    for (uint32_t index = 0; index < map_area; ++index) {
        {
            const mesh::Vec3 q = tectonics::sphere_noise::warp(grid_->position(index), warp_x_noise, warp_y_noise,
                                                               warp_z_noise, warp_km, radius_km);
            const float macro = static_cast<float>(material_macro_noise(q));
            const float detail = static_cast<float>(material_detail_noise(q));
            const float composite = 0.68f * macro + 0.32f * detail;

            uint8_t material_index =
                tectonics::material::to_index(tectonics::material::Type::Granite);
            if (surface_crust_type_map[index] ==
                static_cast<uint8_t>(tectonics::contract::CrustClass::Oceanic)) {
                material_index = tectonics::material::to_index(tectonics::material::Type::Basalt);
            } else if (composite < -0.08f) {
                material_index =
                    tectonics::material::to_index(tectonics::material::Type::Sedimentary);
            } else if (composite > 0.48f) {
                material_index =
                    tectonics::material::to_index(tectonics::material::Type::Metamorphic);
            }

            material_map[index] = material_index;
        }
    }
}

void lithosphere::seedInitialCrustAges() {
    using tectonics::contract::BoundaryType;
    using tectonics::contract::CrustClass;
    using tectonics::contract::GeologicRegime;

    const uint32_t map_area = cellCount();
    if (map_area == 0U || iter_count == 0U) {
        return;
    }

    const double step_myr = std::max(0.001, scenario.delta_time_myr);
    std::vector<int32_t> ocean_distance(map_area, -1);
    std::deque<uint32_t> frontier;

    const auto is_oceanic_cell = [&](uint32_t index) {
        return index < map_area && imap[index] < num_plates && hmap[index] > FLT_EPSILON &&
               crust_class_map[index] == static_cast<uint8_t>(CrustClass::Oceanic);
    };
    const auto push_seed = [&](uint32_t index) {
        if (!is_oceanic_cell(index) || ocean_distance[index] >= 0) {
            return;
        }
        ocean_distance[index] = 0;
        frontier.push_back(index);
    };

    for (uint32_t index = 0; index < map_area; ++index) {
        if (!is_oceanic_cell(index)) {
            continue;
        }

        const auto boundary_type = static_cast<BoundaryType>(boundary_type_map[index]);
        const auto regime = static_cast<GeologicRegime>(geologic_regime_map[index]);
        if (boundary_type == BoundaryType::Divergent || regime == GeologicRegime::DivergentRift ||
            regime == GeologicRegime::MidOceanRidge) {
            push_seed(index);
        }
    }

    if (frontier.empty()) {
        for (uint32_t index = 0; index < map_area; ++index) {
            if (!is_oceanic_cell(index)) {
                continue;
            }
            if (static_cast<BoundaryType>(boundary_type_map[index]) != BoundaryType::None) {
                push_seed(index);
            }
        }
    }

    if (frontier.empty()) {
        for (uint32_t index = 0; index < map_area; ++index) {
            push_seed(index);
        }
    }

    while (!frontier.empty()) {
        const uint32_t current = frontier.front();
        frontier.pop_front();
        const int32_t next_distance = ocean_distance[current] + 1;
        for (const uint32_t neighbor : grid_->neighbours(current)) {
            if (!is_oceanic_cell(neighbor) || ocean_distance[neighbor] >= 0) {
                continue;
            }
            ocean_distance[neighbor] = next_distance;
            frontier.push_back(neighbor);
        }
    }

    // The torus engine's fields on the sphere (seeding_noise), each seeded by its own salt.
    const double radius_km = grid_->radius_m() / 1000;
    const double unit_km = scenario.initial_noise_wavelength_km;
    const auto field = [&](uint32_t salt, int octaves, double persistence, double frequency) {
        SimpleRandom source(static_cast<uint32_t>(scenario.seed) ^ salt);
        return seeding_noise(source.next(), unit_km, radius_km, octaves, persistence, frequency);
    };
    const auto continental_macro_noise = field(0x3C6EF372U, 5, 0.57, 0.26 * 1.3 * 0.75);
    const auto continental_detail_noise = field(0xBB67AE85U, 6, 0.60, 0.72 * 2.5 * 1.7);
    const auto ocean_age_noise = field(0xA54FF53AU, 4, 0.58, 0.34 * 1.6);
    const auto ocean_detail_noise = field(0x510E527FU, 5, 0.61, 0.92 * 3.1 * 1.8);

    uint32_t max_age_steps = 0U;
    std::vector<uint32_t> timestamps(map_area, 0U);

    for (uint32_t index = 0; index < map_area; ++index) {
        {
            const mesh::Vec3 q = grid_->position(index);
            if (imap[index] >= num_plates || hmap[index] <= FLT_EPSILON) {
                continue;
            }

            const auto crust_class = static_cast<CrustClass>(crust_class_map[index]);
            float age_myr = 0.0f;
            if (crust_class == CrustClass::Oceanic) {
                const int32_t distance = std::max(0, ocean_distance[index]);
                const float macro = static_cast<float>(ocean_age_noise(q));
                const float detail = static_cast<float>(ocean_detail_noise(q));
                age_myr = static_cast<float>(distance) * 0.55f;
                if (distance > 0) {
                    age_myr += 2.0f * macro + 1.5f * detail;
                }
                age_myr = std::max(0.0f, std::min(140.0f, age_myr));
            } else {
                const float macro = 0.5f + 0.5f * static_cast<float>(continental_macro_noise(q));
                const float detail = 0.5f + 0.5f * static_cast<float>(continental_detail_noise(q));
                const float normalized = clamp_unit(0.68f * macro + 0.32f * detail);
                const float min_age = crust_class == CrustClass::Transitional ? 55.0f : 110.0f;
                const float span = crust_class == CrustClass::Transitional ? 55.0f : 120.0f;
                age_myr = min_age + span * normalized;
            }

            const uint32_t age_steps =
                static_cast<uint32_t>(std::llround(age_myr / step_myr));
            max_age_steps = std::max(max_age_steps, age_steps);
            timestamps[index] = age_steps;
        }
    }

    if (max_age_steps + MAX_BUOYANCY_AGE > iter_count) {
        iter_count = max_age_steps + MAX_BUOYANCY_AGE;
        time_origin_step = iter_count;
    }

    for (uint32_t index = 0; index < map_area; ++index) {
        {
            const uint32_t owner = imap[index];
            if (owner >= num_plates || hmap[index] <= FLT_EPSILON) {
                amap[index] = 0U;
                continue;
            }

            const uint32_t timestamp =
                timestamps[index] <= iter_count ? iter_count - timestamps[index] : 0U;
            amap[index] = timestamp;
        }
    }
}

float lithosphere::internalToPhysicalMeters(float internal_height) const noexcept {
    return internal_to_physical_meters(internal_height, sea_level_m);
}

float lithosphere::metersToInternalDelta(float delta_m) const noexcept {
    if (delta_m <= 0.0f) {
        return 0.0f;
    }

    const float land_range =
        static_cast<float>(std::max<uint16_t>(1U, TopographyCodec::kMaxHeightMeters - sea_level_m));
    return delta_m / land_range;
}

void lithosphere::importNormalizedHeightMap(const float* normalized_map, float sea_level) {
    if (normalized_map == nullptr) {
        throw invalid_argument("Normalized height map cannot be null");
    }

    const uint32_t map_area = cellCount();
    const float sea_threshold =
        TopographyCodec::infer_normalized_sea_threshold(normalized_map, map_area, sea_level);
    sea_level_m = static_cast<uint16_t>(std::lround(
        sea_threshold * static_cast<float>(TopographyCodec::kMaxHeightMeters)));

    for (uint32_t i = 0; i < map_area; ++i) {
        hmap[i] = TopographyCodec::normalized_to_internal(normalized_map[i], sea_threshold);
    }
    prev_hmap = hmap;
    display_hmap = hmap;
    prev_display_hmap = display_hmap;
    initial_hmap = hmap;
    seedInitialCrustTypes();
    seedMaterialMap();

    resetSimulationState();
    createPlates();
}

void lithosphere::importRawHeightMap(const float* normalized_map) {
    if (normalized_map == nullptr) {
        throw invalid_argument("Normalized height map cannot be null");
    }

    const uint32_t map_area = cellCount();
    std::vector<uint16_t> metric_heightmap(map_area);
    for (uint32_t i = 0; i < map_area; ++i) {
        metric_heightmap[i] = static_cast<uint16_t>(std::lround(
            TopographyCodec::clamp_normalized(normalized_map[i]) *
            static_cast<float>(TopographyCodec::kMaxHeightMeters)));
    }

    initializeHeightMapFromMetric(metric_heightmap.data(), TopographyCodec::default_import_sea_level_m());
    seedInitialCrustTypes();
    seedMaterialMap();
    resetSimulationState();
    createPlates();
}

void lithosphere::importMetricHeightMap(const uint16_t* heightmap_m, uint16_t new_sea_level_m) {
    initializeHeightMapFromMetric(heightmap_m, new_sea_level_m);
    seedInitialCrustTypes();
    seedMaterialMap();
    resetSimulationState();
    createPlates();
}

bool lithosphere::isFinished() const {
    return finished || getPlateCount() == 0;
}

// At least two plates are at same location.
// Move some crust from the SMALLER plate onto LARGER one.
// Update height and plate index maps.
// Doing it plate by plate is much faster than doing it index wise:
// Each plate's map's memory area is accessed sequentially and only
// once as opposed to calculating "num_plates" indices within plate
// maps in order to find out which plate(s) own current location.
void lithosphere::updateTensionMagnitudeMap() noexcept {
    const uint32_t map_area = cellCount();
    if (rift_tension_map.size() != map_area || rift_tension_x_map.size() != map_area ||
        rift_tension_y_map.size() != map_area) {
        return;
    }

    for (uint32_t index = 0; index < map_area; ++index) {
        const float tx = rift_tension_x_map[index];
        const float ty = rift_tension_y_map[index];
        rift_tension_map[index] = std::sqrt(tx * tx + ty * ty);
    }
}

// Remove empty plates from the system.
void lithosphere::removeEmptyPlates() {
    for (uint32_t i = 0; i < num_plates; ++i) {
        if (num_plates == 1)
            puts("ONLY ONE PLATE LEFT!");
        else if (plate_indices_found[i] == 0) {
            topology_model.removePlate(i);
            plates[i] = plates[num_plates - 1];
            plates.pop_back();
            plate_indices_found[i] = plate_indices_found[num_plates - 1];
            plate_indices_found.pop_back();

            // Life is seldom as simple as seems at first.
            // Replace the moved plate's index in the index map
            // to match its current position in the array!
            for (uint32_t j = 0; j < cellCount(); ++j)
                if (imap[j] == num_plates - 1)
                    imap[j] = i;

            --num_plates;
            --i;
        }
    }
}

std::vector<tectonics::dynamics::MantleSample> lithosphere::mantleSamples(double time_myr) const {
    if (mantle_cache_time == time_myr && mantle_cache.size() == cellCount()) return mantle_cache;
    const double radius = grid_->radius_m();
    std::vector<tectonics::dynamics::MantleSample> samples(cellCount());
    for (uint32_t c = 0; c < cellCount(); ++c) samples[c] = mantle_flow->at(grid_->position(c), time_myr, radius);
    mantle_cache = samples;
    mantle_cache_time = time_myr;
    return samples;
}

std::vector<mesh::Vec3> lithosphere::plateRotations() const {
    const double dt_s = scenario.delta_time_myr * tectonics::dynamics::kSecondsPerMyr;
    std::vector<mesh::Vec3> rotation(num_plates);
    for (uint32_t p = 0; p < num_plates; ++p) rotation[p] = dt_s * plates[p].rotation_rad_s;
    return rotation;
}

std::vector<mesh::Vec3> lithosphere::rigidDisplacement() const {
    tectonics::sphere_transport::Motion motion;
    motion.rotation = plateRotations();
    motion.plate.assign(imap.begin(), imap.end());
    return tectonics::columns::Model::displacement(*grid_, motion);
}

void lithosphere::solvePlateMotion() {
    using namespace tectonics::dynamics;
    PhaseTimer timer(phase_timings, Phase::motion);
    const double radius = grid_->radius_m();
    const double dt_s = scenario.delta_time_myr * kSecondsPerMyr;
    const std::vector<tectonics::grid::Edge>& edges = grid_->edges();
    std::vector<SphereContact> contacts;
    std::vector<size_t> contact_faces; // edge ids
    std::map<std::tuple<uint32_t,uint32_t,uint32_t>,std::array<double,3>> activity;
    for (size_t face=0;face<motion_history.size() && face<edges.size();++face) {
        const auto& sample=motion_history[face];
        if (sample.a==sample.b || sample.a==UINT32_MAX) continue;
        auto& sum=activity[{sample.a,sample.b,static_cast<uint32_t>(edges[face].axis)}];
        sum[0]+=sample.ridge;sum[1]+=sample.slab;sum[2]+=1;
    }
    std::vector<uint32_t> pressure_owners;
    std::vector<double> pressure;
    if (scenario.mountain_force_range_km>0) {
        pressure_owners = imap;
        pressure.reserve(pressure_owners.size());
        for (const auto& c:column_model.cells()) pressure.push_back(tectonics::columns::mountain_force_n_m(c,scenario.sediment_load_feedback,
            scenario.continental_root_kg_m3));
    }
    // Plate drag moments and their mantle forcing.
    const auto mantle = mantleSamples(getTimeMyr());
    const std::vector<SpherePlate> bodies = sphere_plates(*grid_, imap, num_plates, mantle);
    for (const auto& body : bodies)
        if (body.area_m2 <= 0.0) throw runtime_error("Force balance found an empty plate owner");
    // One contact per edge between two plates, its normal from the owner of edge.a.
    for (size_t face = 0; face < edges.size(); ++face) {
        const tectonics::grid::Edge& edge = edges[face];
        const uint32_t index = edge.a, neighbor = edge.b, owner = imap[index], other = imap[neighbor];
        if (owner >= num_plates || other == owner || other >= num_plates) continue;
        contacts.push_back(sphere_contact(*grid_, edge, imap));
        SphereContact& contact = contacts.back();
        contact.continental_a = is_continental_material(material_map[index]);
        contact.continental_b = is_continental_material(material_map[neighbor]);
        contact.age_a_myr = (iter_count - std::min(iter_count, amap[index])) * scenario.delta_time_myr;
        contact.age_b_myr = (iter_count - std::min(iter_count, amap[neighbor])) * scenario.delta_time_myr;
        contact.closing_m_s = sphere_closing_m_s(contact, plates[owner].rotation_rad_s, plates[other].rotation_rad_s, radius);
        contact.activation_speed_m_s = scenario.boundary_force_ramp_cm_yr / 3.15576e9;
        contact_faces.push_back(face);
        if (scenario.consistent_motion) {
            const auto& stable=topology_model.diagnostics().stable_ids;
            const auto found=subduction_polarity.find(std::minmax(stable[owner],stable[other]));
            if (found!=subduction_polarity.end())
                contact.descending_side=found->second==stable[owner]?0:1;
            else {
                const auto& ca=column_model.cells()[index];
                const auto& cb=column_model.cells()[neighbor];
                const double aa=ca.ocean>0?ca.oceanic_age_volume/ca.ocean:0;
                const double ab=cb.ocean>0?cb.oceanic_age_volume/cb.ocean:0;
                contact.descending_side=aa>=ab?0:1;
                if (!contact.continental_a && !contact.continental_b)
                    subduction_polarity.emplace(std::minmax(stable[owner],stable[other]),
                        contact.descending_side==0?stable[owner]:stable[other]);
            }
        }
        {
            const auto& cells = column_model.cells();
            if (scenario.thermal_force_feedback) {
                const auto thermal_age=[](const tectonics::columns::Cell& c) {
                    return c.ocean>0?std::max(0.0,
                        (c.oceanic_age_volume-c.thermal_rejuvenation_volume)/c.ocean):0.0;
                };
                contact.age_a_myr=thermal_age(cells[index]);
                contact.age_b_myr=thermal_age(cells[neighbor]);
                contact.slab_pull_scale=scenario.thermal_subsidence_strength;
            }
            contact.mountain_force_n_m = .5 * (tectonics::columns::mountain_force_n_m(cells[index],scenario.sediment_load_feedback,
                    scenario.continental_root_kg_m3)
                + tectonics::columns::mountain_force_n_m(cells[neighbor],scenario.sediment_load_feedback,
                    scenario.continental_root_kg_m3));
            // The inland strips walk away from the contact on each side: from edge.a
            // opposite the edge's direction, from edge.b straight on in its own frame.
            if (!pressure.empty())
                contact.mountain_force_n_m=.5*(inland_pressure(*grid_,pressure_owners,pressure,
                    index,tectonics::grid::opposite(edge.direction),scenario.mountain_force_range_km*1000)
                    +inland_pressure(*grid_,pressure_owners,pressure,
                    neighbor,grid_->ahead(index,edge.direction).direction,scenario.mountain_force_range_km*1000));
            // The same excess crustal GPE can resist compression at an
            // active continental margin. Previously the force solver
            // discarded this pressure whenever the other side was oceanic.
            if (!contact.continental_a || !contact.continental_b)
                contact.mountain_force_n_m *= (contact.continental_a || contact.continental_b)
                    ? scenario.coastal_gpe_strength : 0;
        }
    }
    if (scenario.consistent_motion) {
        // Initial activity follows the prescribed mantle forcing, not the seeded rotations:
        // the closing of each plate's motion under its basal drag alone.
        std::vector<Vec3> free(num_plates);
        for (uint32_t p = 0; p < num_plates; ++p) free[p] = drag_rotation(bodies[p], radius);
        for (size_t k=0;k<contacts.size();++k) {
            auto& edge=contacts[k];
            const auto& stable=topology_model.diagnostics().stable_ids;
            const uint32_t a=stable[edge.a],b=stable[edge.b],axis=edges[contact_faces[k]].axis;
            const auto found=activity.find({a,b,axis});
            if (contact_faces[k]<motion_history.size() && motion_history[contact_faces[k]].a==a && motion_history[contact_faces[k]].b==b) {
                edge.ridge_activity=motion_history[contact_faces[k]].ridge;
                edge.slab_activity=motion_history[contact_faces[k]].slab;
            } else if (found!=activity.end()) {
                edge.ridge_activity=found->second[0]/found->second[2];
                edge.slab_activity=found->second[1]/found->second[2];
            } else {
                // Newly exposed faces inherit their pair.
                const double closing=sphere_closing_m_s(edge,free[edge.a],free[edge.b],radius);
                const double strength=std::abs(closing)/std::hypot(closing,edge.activation_speed_m_s);
                edge.ridge_activity=closing<0?strength:0;edge.slab_activity=closing>0?strength:0;
            }
            edge.activity_response=scenario.delta_time_myr/(scenario.boundary_memory_myr+scenario.delta_time_myr);
        }
    }
    std::vector<Vec3> old_rotation(num_plates);
    for (uint32_t p=0;p<num_plates;++p) old_rotation[p]=plates[p].rotation_rad_s;
    const auto solution = scenario.consistent_motion
        ? solve_sphere_consistent(bodies,contacts,old_rotation,radius) : solve_sphere(bodies, contacts, radius);
    const auto previous_diagnostics=motion_diagnostics;
    motion_diagnostics = solution.diagnostics;
    motion_diagnostics.max_velocity_consistency_cm_yr=std::max(previous_diagnostics.max_velocity_consistency_cm_yr,
        motion_diagnostics.velocity_consistency_cm_yr);
    motion_diagnostics.persistent_face_steps=previous_diagnostics.persistent_face_steps;
    motion_diagnostics.reversing_face_steps=previous_diagnostics.reversing_face_steps;
    motion_diagnostics.strong_face_steps=previous_diagnostics.strong_face_steps;
    motion_diagnostics.strong_reversing_face_steps=previous_diagnostics.strong_reversing_face_steps;
    sphere_measures(bodies, contacts, solution.rotation_rad_s, radius, motion_diagnostics);
    if (getTimeMyr()>motion_history_time) {
        // The motion history of every edge: its plate pair and closing speed (cm/yr).
        motion_history.resize(edges.size());
        const auto& stable=topology_model.diagnostics().stable_ids;
        for (size_t face=0;face<edges.size();++face) {
            const auto& edge=edges[face];
            const auto a=imap[edge.a],b=imap[edge.b];
            auto& sample=motion_history[face];
            const double closing=mesh::dot(mesh::cross(solution.rotation_rad_s[a]-solution.rotation_rad_s[b],
                radius*edge.midpoint),edge.normal)*3.15576e9;
            if (a!=b && sample.a==stable[a] && sample.b==stable[b]) {
                ++motion_diagnostics.persistent_face_steps;
                const bool reversal=sample.closing*closing<0;
                motion_diagnostics.reversing_face_steps+=reversal;
                if (std::abs(sample.closing)>.5 && std::abs(closing)>.5) {
                    ++motion_diagnostics.strong_face_steps;
                    motion_diagnostics.strong_reversing_face_steps+=reversal;
                }
            }
            sample={stable[a],stable[b],closing,0,0};
        }
        if (scenario.consistent_motion) for (size_t k=0;k<contacts.size();++k) {
            const auto& edge=contacts[k];auto& sample=motion_history[contact_faces[k]];
            const double strength=std::abs(sample.closing)/std::hypot(sample.closing,scenario.boundary_force_ramp_cm_yr);
            sample.ridge=edge.ridge_activity+edge.activity_response*((sample.closing<0?strength:0)-edge.ridge_activity);
            sample.slab=edge.slab_activity+edge.activity_response*((sample.closing>0?strength:0)-edge.slab_activity);
        }
        motion_history_time=getTimeMyr();
    }
    for (uint32_t p = 0; p < num_plates; ++p) {
        const Vec3 w = solution.rotation_rad_s[p];
        if (!std::isfinite(w.x) || !std::isfinite(w.y) || !std::isfinite(w.z))
            throw runtime_error("Force balance produced a non-finite plate rotation");
        plates[p].rotation_rad_s = w;
    }
    // The transport cannot resolve a cell skipping its neighbour. Fail explicitly rather than
    // hiding excessive velocities behind a speed clamp: each cell moves at most 0.95 of its
    // smallest edge per update.
    double largest = 0;
    for (uint32_t c = 0; c < cellCount(); ++c) {
        const double moved = mesh::length(mesh::cross(plates[imap[c]].rotation_rad_s, radius * grid_->position(c))) * dt_s;
        largest = std::max(largest, moved / grid_->smallestedge_m(c));
    }
    if (largest > tectonics::columns::displacement_limit)
        throw runtime_error("Force-balance displacement exceeds 0.95 of a cell's smallest edge; reduce delta_time_myr");
}

void lithosphere::initializeColumns() {
    subduction_polarity.clear();
    motion_history.clear();motion_history_time=-1;
    updateCrustAgeMyrMap();
    std::vector<tectonics::columns::Cell> cells(cellCount());
    for (uint32_t i = 0; i < cellCount(); ++i)
        cells[i] = tectonics::columns::seed(internalToPhysicalMeters(hmap[i]) - sea_level_m,
            hmap[i] >= CONTINENTAL_BASE || is_continental_material(material_map[i]), crust_age_myr_map[i],
            scenario.continental_root_kg_m3, scenario.continental_root_max_thickness_m, scenario.continental_root_delamination_width_m);
    const tectonics::ocean_relief::Parameters relief{scenario.ocean_crust_variation_m, scenario.ocean_noise_wavelength_km,
        scenario.seamount_crust_m,scenario.seamount_radius_km,scenario.seamount_spacing_km,scenario.seamount_roughness};
    column_model.initialize(grid_, std::move(cells), scenario.transport_subdivision, scenario.seed, relief,
        false, scenario.continental_root_kg_m3, scenario.continental_root_max_thickness_m, scenario.continental_root_delamination_width_m);
    topology_model.initialize(grid_, num_plates, imap);
    if (scenario.boundary_conditioned_terrain) {
        // Use the initial force-balanced motion when available, so terrain is
        // conditioned on opening contacts that will actually drive the first step.
        solvePlateMotion();
        std::vector<double> elevations(cellCount()),ages(cellCount());
        for (uint32_t i=0;i<cellCount();++i) {
            elevations[i]=internalToPhysicalMeters(hmap[i])-sea_level_m;
            ages[i]=crust_age_myr_map[i];
        }
        auto conditioned=tectonics::initialization::terrain(column_model.metric(),scenario.delta_time_myr,
            elevations,ages,imap,rigidDisplacement(),scenario.thermal_subsidence_strength,scenario.initial_margin_width_m,
            scenario.initial_boundary_influence,scenario.continental_root_kg_m3,
            scenario.continental_root_max_thickness_m, scenario.continental_root_delamination_width_m);
        column_model.initialize(grid_,std::move(conditioned),scenario.transport_subdivision,
            scenario.seed,relief,true,
            scenario.continental_root_kg_m3,scenario.continental_root_max_thickness_m, scenario.continental_root_delamination_width_m);
    }
    publishColumns();
    // Conditioning changes buoyancy/GPE. Publish motion for that actual initial
    // surface, instead of retaining the solve from the preconditioned surface.
    if (scenario.consistent_motion) solvePlateMotion();
    prev_hmap = hmap; prev_display_hmap = display_hmap;
}

void lithosphere::publishColumns() {
    const auto& cells = column_model.cells();
    for (uint32_t i = 0; i < cellCount(); ++i) {
        const auto& c = cells[i];
        const double z = tectonics::columns::elevation_m(c, scenario.thermal_subsidence_strength,
            scenario.continental_root_kg_m3, scenario.continental_root_max_thickness_m, scenario.continental_root_delamination_width_m);
        if (!std::isfinite(z) || z <= -sea_level_m)
            throw runtime_error("Isostatic elevation cannot be represented by the native height datum");
        hmap[i] = elevation_to_internal(static_cast<float>(z), sea_level_m);
        const auto material = c.continental_area >= .5 ? tectonics::material::Type::Granite
                                                     : tectonics::material::Type::Basalt;
        material_map[i] = tectonics::material::to_index(material);
        const double age_steps = tectonics::columns::age_myr(c) / scenario.delta_time_myr;
        amap[i] = iter_count - std::min(iter_count, static_cast<uint32_t>(std::round(age_steps)));
        // The owner's record: the crust, and the material unless recorded continental crust
        // would turn oceanic (it stays continental until the cell changes owner).
        float crust = hmap[i];
        if (!std::isfinite(crust) || crust < 0) crust = 0;
        if (crust > 0.0f && !(owner_crust_map[i] > FLT_EPSILON && is_continental_material(owner_material_map[i]) &&
                              !is_continental_material(material_map[i])))
            owner_material_map[i] = material_map[i];
        owner_crust_map[i] = crust;
    }
    display_hmap = hmap;
    const auto& uplift_rate = column_model.upliftRate();
    for (size_t i = 0; i < rock_uplift_rate_map.size(); ++i)
        rock_uplift_rate_map[i] = i < uplift_rate.size() ? static_cast<float>(uplift_rate[i]) : 0.0f;
    const auto& strain = column_model.strain();
    const auto& metric = column_model.metric();
    const double cell_m = grid_->cellspacing_m();
    constexpr double pi = 3.14159265358979323846;
    for (size_t i = 0; i < shortening_map.size(); ++i) {
        const auto s = i < strain.size() ? strain[i] : tectonics::columns::Strain{};
        shortening_map[i] = static_cast<float>(tectonics::columns::shortening_m(s, cell_m));
        // The axis from the cell's stencil axes to geographic east toward south.
        const auto c = static_cast<tectonics::grid::CellId>(i);
        const double local = tectonics::columns::shortening_direction(s);
        const mesh::Vec3 axis = std::cos(local) * metric.axis(c, 0) + std::sin(local) * metric.axis(c, 1);
        const mesh::Vec3 p = grid_->position(c);
        const double latitude = mesh::latitudeof(p), longitude = mesh::longitudeof(p);
        double direction = std::atan2(-mesh::dot(axis, mesh::northof(latitude, longitude)),
                                      mesh::dot(axis, mesh::eastof(latitude, longitude)));
        if (direction < 0) direction += pi;
        if (direction >= pi) direction -= pi;
        // An axis just below pi rounds to the float pi, the same axis as 0; equal principal
        // components have no axis (0).
        const auto rounded = static_cast<float>(direction);
        shortening_direction_map[i] = rounded >= static_cast<float>(pi) || std::hypot((s.xx - s.yy) / 2, s.xy) == 0
            ? 0.0f : rounded;
    }
}

void lithosphere::updateColumns() {
    prev_hmap = hmap; prev_display_hmap = display_hmap; prev_imap = imap;
    solvePlateMotion();
    const std::vector<mesh::Vec3> rigid = rigidDisplacement();
    std::vector<uint32_t> material_owners(rigid.size());
    std::vector<double> continental_fraction(rigid.size());
    for (uint32_t i=0; i<rigid.size(); ++i) {
        material_owners[i]=imap[i]; continental_fraction[i]=column_model.cells()[i].continental_area;
    }
    auto stable_owners=material_owners;
    for (auto& owner:stable_owners) owner=topology_model.diagnostics().stable_ids[owner];
    tectonics::margins::Forcing margin_forcing;
    {
        PhaseTimer timer(phase_timings, Phase::margins);
        margin_forcing=tectonics::margins::forcing(column_model.metric(),scenario.delta_time_myr,
            stable_owners,column_model.cells(),rigid,scenario.margin_processes,&subduction_polarity,
            scenario.seed,getTimeMyr());
    }
    std::vector<mesh::Vec3> displacement;
    {
        PhaseTimer timer(phase_timings, Phase::deformation);
        displacement=tectonics::deformation::distribute(column_model.metric(),material_owners,continental_fraction,
            rigid,scenario.seed,{scenario.deformation_peak_km,scenario.deformation_variation,scenario.subduction_coupling,
            scenario.margin_processes.process_variation,getTimeMyr(),scenario.margin_processes.process_wavelength_km,
            scenario.margin_processes.process_period_myr,scenario.deformation_variation_wavelength_km});
    }
    if (scenario.topology_processes.progressive_rifting) {
        PhaseTimer timer(phase_timings, Phase::rift_displacement);
        const auto rift=topology_model.riftDisplacement(material_owners,scenario.delta_time_myr,scenario.topology_processes);
        for (size_t i=0;i<displacement.size();++i) displacement[i]=displacement[i]+rift[i];
    }
    // The update's motion: the plates' rotations, each cell's plate and its offsets beyond the
    // rigid motion (the distributed deformation and the rift opening).
    tectonics::sphere_transport::Motion motion;
    motion.rotation = plateRotations();
    motion.plate.assign(imap.begin(), imap.end());
    motion.offset_m.resize(displacement.size());
    for (size_t i=0;i<displacement.size();++i) motion.offset_m[i]=displacement[i]-rigid[i];
    {
        // Published: (east, north) in cell spacings per update.
        const auto moved = tectonics::columns::Model::displacement(*grid_, motion);
        material_displacement.resize(moved.size());
        const double spacing = grid_->cellspacing_m();
        for (size_t i=0;i<moved.size();++i) {
            const auto en = tectonics::boundary_cube::eastnorth(grid_->position(static_cast<tectonics::grid::CellId>(i)), moved[i]);
            material_displacement[i] = {en[0] / spacing, en[1] / spacing};
        }
    }
    tectonics::plumes::Source plume;
    if (scenario.plume_processes.count>0) {
        PhaseTimer timer(phase_timings, Phase::plumes);
        plume=tectonics::plumes::source(*grid_,scenario.seed,getTimeMyr(),
            getTimeMyr()+scenario.delta_time_myr,scenario.plume_processes);
    }
    {
        PhaseTimer timer(phase_timings, Phase::columns);
        column_model.advance(motion, scenario.delta_time_myr,
            scenario.erosion_strength, scenario.thermal_subsidence_strength, scenario.submarine_erosion_ratio,
            &margin_forcing,scenario.margin_processes,scenario.experimental_noise,scenario.sediment_load_feedback,
            plume.basalt_m.empty()?nullptr:&plume,scenario.uplift_rate_window_myr);
    }
    plume_influence_map.resize(plume.influence.size());
    std::transform(plume.influence.begin(),plume.influence.end(),plume_influence_map.begin(),
        [](double value){return static_cast<float>(value);});
    std::vector<uint32_t> owners;
    {
        PhaseTimer timer(phase_timings, Phase::topology_advance);
        owners = topology_model.advance(motion.rotation);
    }
    std::vector<tectonics::topology::Change> changes;
    {
        PhaseTimer timer(phase_timings, Phase::topology_evolve);
        const auto samples = mantleSamples(getTimeMyr() + scenario.delta_time_myr);
        std::vector<mesh::Vec3> mantle(owners.size()), rotation(num_plates);
        std::vector<double> continent(owners.size());
        std::vector<double> strength(owners.size());
        for (uint32_t i = 0; i < owners.size(); ++i) {
            mantle[i] = samples[i].velocity_m_s;
            continent[i] = column_model.cells()[i].continental_area;
            const auto& c=column_model.cells()[i];
            const double thickness=c.continental_area>1e-6?c.continent/c.continental_area:35000;
            const double thermal_age=c.ocean>0?std::max(0.0,(c.oceanic_age_volume-c.thermal_rejuvenation_volume)/c.ocean):60;
            const double weakening=plume_influence_map.empty()?0:
                scenario.plume_processes.weakening*plume_influence_map[i];
            strength[i]=std::clamp((continent[i]*std::pow(thickness/35000,2)+
                (1-continent[i])*std::sqrt((thermal_age+5)/65))*(1-weakening),.15,4.0);
        }
        for (uint32_t p = 0; p < num_plates; ++p) rotation[p] = plates[p].rotation_rad_s;
        changes = topology_model.evolve(owners, mantle, rotation, continent,
            scenario.delta_time_myr, getTimeMyr() + scenario.delta_time_myr,scenario.topology_processes,strength);
    }
    {
    PhaseTimer ownership_timer(phase_timings, Phase::ownership_publish);
    std::vector<uint8_t> newborn(num_plates, 0U);
    for (const auto& change : changes) if (change.kind == tectonics::topology::EventKind::Split) {
        if (change.b != plates.size()) throw runtime_error("A split plate's index is not the next plate");
        // The newborn inherits motion until the next force solve, with no kick.
        tectonics::PlateState daughter = tectonics::PlateState::fromSeed(static_cast<uint32_t>(scenario.seed),
            grid_->cellspacing_m(), grid_->radius_m(), scenario.delta_time_myr * tectonics::dynamics::kSecondsPerMyr);
        daughter.rotation_rad_s = plates[change.a].rotation_rad_s;
        plates.push_back(daughter);
        newborn.push_back(1U);
        ++num_plates;
    }
    // A newborn's creation share counts its cells with crust (above 0), as at creation.
    std::vector<uint32_t> born_occupied(num_plates, 0U), born_continental(num_plates, 0U);
    plate_indices_found.assign(num_plates, 0U);
    for (uint32_t i = 0; i < cellCount(); ++i) {
        // A new owner records no crust at the cell, except a newborn plate, which starts with
        // the crust and material of its cells. Material columns are untouched.
        if (imap[i] != owners[i]) {
            if (newborn[owners[i]]) {
                owner_crust_map[i] = hmap[i];
                owner_material_map[i] = material_map[i];
                if (hmap[i] > 0.0f) {
                    ++born_occupied[owners[i]];
                    born_continental[owners[i]] += is_continental_material(material_map[i]) ? 1U : 0U;
                }
            } else {
                owner_crust_map[i] = 0.0f;
            }
        }
        imap[i] = owners[i];
        ++plate_indices_found[imap[i]];
    }
    for (uint32_t p = 0; p < num_plates; ++p)
        if (newborn[p])
            plates[p].creation_continental_share = born_occupied[p] > 0U
                ? static_cast<float>(born_continental[p]) / static_cast<float>(born_occupied[p]) : 0.0f;
    ++iter_count;
    publishColumns();
    removeEmptyPlates();
    }
    if (!changes.empty()) solvePlateMotion();
    updateDerivedMaps();
    // Preserve fractional material ages rather than rounding them to updates.
    for (size_t i = 0; i < cellCount(); ++i)
        crust_age_myr_map[i] = static_cast<float>(tectonics::columns::age_myr(column_model.cells()[i]));
    finished = iter_count - time_origin_step >= static_cast<uint32_t>(
        std::ceil(scenario.evolution_duration_myr / scenario.delta_time_myr));
}

uint32_t lithosphere::getWidth() const {
    return static_cast<uint32_t>(grid_->facecells());
}

uint32_t lithosphere::getHeight() const {
    return 6U * static_cast<uint32_t>(grid_->facecells());
}

uint32_t* lithosphere::getPlatesMap() const throw() {
    // The C interface hands out a mutable pointer; the engine never reads it back.
    return const_cast<uint32_t*>(imap.data());
}

const tectonics::PlateState& lithosphere::getPlate(uint32_t index) const {
    ASSERT(index < num_plates, "invalid plate index");
    return plates[index];
}

void lithosphere::update() {
    if (isFinished()) return;
    ++_steps;
    updateColumns();
    ++phase_timings.updates;
}
