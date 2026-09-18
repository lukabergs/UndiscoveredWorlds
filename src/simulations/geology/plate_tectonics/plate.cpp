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

#define _USE_MATH_DEFINES // Winblow$.
#include <cfloat>         // FT_EPSILON
#ifdef __MINGW32__ // this is to avoid a problem with the hypot function which is messed up by
                   // Python...
    #undef __STRICT_ANSI__
#endif
#include <assert.h>
#include <cmath>     // sin, cos
#include <cstdlib>   // rand
#include <stdexcept> // std::invalid_argument
#include <vector>

#include "heightmap.hpp"
#include "material.hpp"
#include "plate.hpp"
#include "plate_functions.hpp"
#include "rectangle.hpp"
#include "simplexnoise.hpp"
#include "utils.hpp"

using namespace std;

namespace {

static const float kRiverErosionStrength = 0.055f;
static const float kRainfallRedistributionStrength = 0.14f;
static const float kWindRedistributionStrength = 0.08f;
static const float kSlopeNormalization = 0.35f;
static const float kWindRandomVariation = 0.45f;
static const float kUpliftCompensationStrength = 1.35f;
static const float kUpliftDecaySteps = 18.0f;
static const float kMinRasterRotationAngle = 0.015f;
static const float kContinentalBase = 1.0f;
static const float kOceanicPlateMaxContinentalAreaRatio = 0.10f;

bool material_index_is_continental(uint8_t material_index) {
    return platec::material::from_index(material_index) != platec::material::Type::Basalt;
}

float plate_buoyancy_from_metrics(PlateType type, float continental_area_ratio,
                                  float continental_mass_ratio, float mean_crust) {
    return type == PlateType::Continental
        ? 1.10f + 0.85f * continental_area_ratio + 0.35f * continental_mass_ratio +
              0.10f * mean_crust
        : 0.25f + 0.20f * continental_area_ratio + 0.55f * mean_crust;
}

float clamp_unit(float value) {
    return std::max(0.0f, std::min(1.0f, value));
}

float wrapped_delta(float point, float center, uint32_t world_size) {
    float delta = point - center;
    const float period = static_cast<float>(world_size);
    const float half_period = period * 0.5f;
    if (delta > half_period) {
        delta -= period;
    } else if (delta < -half_period) {
        delta += period;
    }
    return delta;
}

float wrap_coordinate(float value, uint32_t world_size) {
    const float period = static_cast<float>(world_size);
    while (value < 0.0f) {
        value += period;
    }
    while (value >= period) {
        value -= period;
    }
    return value;
}

uint32_t wrap_index(int value, uint32_t size) {
    ASSERT(size > 0U, "Wrap size must be positive");
    int wrapped = value % static_cast<int>(size);
    if (wrapped < 0) {
        wrapped += static_cast<int>(size);
    }
    return static_cast<uint32_t>(wrapped);
}

float rainfall_from_latitude(uint32_t world_y, uint32_t world_height) {
    if (world_height == 0) {
        return 0.5f;
    }
    const float latitude =
        (static_cast<float>(world_y) + 0.5f) / static_cast<float>(world_height);
    return std::sin(latitude * PI);
}

float normalized_slope(float center, float west, float east, float north, float south) {
    const float west_height = west > 0.0f ? west : center;
    const float east_height = east > 0.0f ? east : center;
    const float north_height = north > 0.0f ? north : center;
    const float south_height = south > 0.0f ? south : center;
    const float dx = 0.5f * (east_height - west_height);
    const float dy = 0.5f * (south_height - north_height);
    const float slope = std::sqrt(dx * dx + dy * dy);
    return slope / (slope + kSlopeNormalization);
}

float uplift_compensation(float crust, float lower_bound, uint32_t current_iteration,
                          uint32_t timestamp) {
    if (crust <= lower_bound) {
        return 0.0f;
    }
    const float relief = crust - lower_bound;
    const float relief_factor = relief / (relief + kSlopeNormalization);
    const uint32_t crust_age =
        current_iteration > timestamp ? current_iteration - timestamp : 0U;
    const float youth = std::exp(-static_cast<float>(crust_age) / kUpliftDecaySteps);
    return 1.0f + kUpliftCompensationStrength * relief_factor * youth;
}

float fracture_edge_exposure(const HeightMap& map, uint32_t width, uint32_t height, uint32_t x,
                             uint32_t y, bool wrap_x, bool wrap_y) {
    const uint32_t index = y * width + x;
    if (map[index] <= 0.0f) {
        return 0.0f;
    }

    uint32_t empty_neighbors = 0U;
    uint32_t total_neighbors = 0U;
    for (int dy = -1; dy <= 1; ++dy) {
        for (int dx = -1; dx <= 1; ++dx) {
            if (dx == 0 && dy == 0) {
                continue;
            }

            ++total_neighbors;
            int nx = static_cast<int>(x) + dx;
            int ny = static_cast<int>(y) + dy;
            const bool outside_x = nx < 0 || nx >= static_cast<int>(width);
            const bool outside_y = ny < 0 || ny >= static_cast<int>(height);

            if ((outside_x && !wrap_x) || (outside_y && !wrap_y)) {
                ++empty_neighbors;
                continue;
            }

            if (wrap_x) {
                nx = static_cast<int>(wrap_index(nx, width));
            }
            if (wrap_y) {
                ny = static_cast<int>(wrap_index(ny, height));
            }

            if (map[static_cast<uint32_t>(ny) * width + static_cast<uint32_t>(nx)] <= 0.0f) {
                ++empty_neighbors;
            }
        }
    }

    return total_neighbors > 0U
               ? static_cast<float>(empty_neighbors) / static_cast<float>(total_neighbors)
               : 0.0f;
}

float fracture_warp_noise(float x, float y, float z0, float z1) {
    const float macro =
        scaled_octave_noise_3d(4.0f, 0.56f, 0.17f, -1.0f, 1.0f, x, y, z0);
    const float detail =
        scaled_octave_noise_3d(2.0f, 0.50f, 0.43f, -1.0f, 1.0f, x + 19.0f, y + 31.0f, z1);
    return 0.68f * macro + 0.32f * detail;
}

float fracture_warp_noise_wrapped(float x, float y, float linear_x, float linear_y,
                                  uint32_t width, uint32_t height, bool wrap_x, bool wrap_y,
                                  float z0, float z1) {
    if (!wrap_x && !wrap_y) {
        return fracture_warp_noise(linear_x, linear_y, z0, z1);
    }

    if (wrap_x && wrap_y) {
        const float theta_x =
            ((x + z0 * 0.031f) / static_cast<float>(width)) * 2.0f * static_cast<float>(PI);
        const float theta_y =
            ((y + z1 * 0.037f) / static_cast<float>(height)) * 2.0f * static_cast<float>(PI);
        const float macro =
            scaled_octave_noise_4d(4.0f, 0.56f, 0.17f, -1.0f, 1.0f, cosf(theta_x), sinf(theta_x),
                                   cosf(theta_y), sinf(theta_y));
        const float detail = scaled_octave_noise_4d(
            2.0f, 0.50f, 0.43f, -1.0f, 1.0f, cosf(theta_x + 0.61f), sinf(theta_x + 0.61f),
            cosf(theta_y + 1.17f), sinf(theta_y + 1.17f));
        return 0.68f * macro + 0.32f * detail;
    }

    if (wrap_x) {
        const float theta_x =
            ((x + z0 * 0.031f) / static_cast<float>(width)) * 2.0f * static_cast<float>(PI);
        const float macro =
            scaled_octave_noise_4d(4.0f, 0.56f, 0.17f, -1.0f, 1.0f, cosf(theta_x), sinf(theta_x),
                                   linear_y, z0);
        const float detail =
            scaled_octave_noise_4d(2.0f, 0.50f, 0.43f, -1.0f, 1.0f, cosf(theta_x + 0.61f),
                                   sinf(theta_x + 0.61f), linear_y + 31.0f, z1);
        return 0.68f * macro + 0.32f * detail;
    }

    const float theta_y =
        ((y + z1 * 0.037f) / static_cast<float>(height)) * 2.0f * static_cast<float>(PI);
    const float macro =
        scaled_octave_noise_4d(4.0f, 0.56f, 0.17f, -1.0f, 1.0f, linear_x, cosf(theta_y),
                               sinf(theta_y), z0);
    const float detail =
        scaled_octave_noise_4d(2.0f, 0.50f, 0.43f, -1.0f, 1.0f, linear_x + 19.0f,
                               cosf(theta_y + 1.17f), sinf(theta_y + 1.17f), z1);
    return 0.68f * macro + 0.32f * detail;
}

float subduction_wavelet(float x, float trench_scale, float arc_scale) {
    const float scale = x < 0.0f ? trench_scale : arc_scale;
    if (scale <= FLT_EPSILON) {
        return 0.0f;
    }
    return (x / scale) * std::exp(-(x * x) / (2.0f * scale * scale));
}

} // namespace

plate::plate(long seed, float* m, uint32_t w, uint32_t h, uint32_t _x, uint32_t _y,
             uint32_t plate_age, WorldDimension worldDimension, const uint8_t* local_material_map,
             float erosion_strength, float crust_rotation_strength, float rotation_strength,
             float movement_energy)
    : _worldDimension(worldDimension), _randsource(seed), map(m, w, h), age_map(w, h),
      material_map(w, h), submerged_map(w, h), submerged_age_map(w, h),
      submerged_since_map(w, h), submerged_material_map(w, h),
      _bounds(nullptr), _mass(MassBuilder(m, Dimension(w, h)).build()),
      _movement(_randsource, worldDimension, rotation_strength), _plate_type(PlateType::Oceanic),
      _buoyancy(0.0f), _submerged_mass(0.0f), _continental_area_ratio(0.0f),
      _continental_mass_ratio(0.0f), _mean_crust(0.0f),
      _erosion_strength(erosion_strength < 0.0f ? 0.0f : erosion_strength),
      _crust_rotation_strength(crust_rotation_strength < 0.0f ? 0.0f : crust_rotation_strength),
      _movement_energy(movement_energy < 0.0f ? 0.0f : movement_energy),
      _pending_crust_rotation(0.0f),
      _segments(nullptr), _mySegmentCreator(nullptr) {
    const uint32_t plate_area = w * h;
    uint32_t occupied_area = 0;
    uint32_t continental_area = 0;
    float occupied_mass = 0.0f;
    float continental_mass = 0.0f;

    _bounds = new Bounds(worldDimension, FloatPoint(static_cast<float>(_x), static_cast<float>(_y)),
                         Dimension(w, h));
    submerged_map.set_all(0.0f);
    submerged_age_map.set_all(0U);
    submerged_since_map.set_all(0U);
    submerged_material_map.set_all(platec::material::to_index(platec::material::kDefaultType));

    uint32_t k;
    for (uint32_t y = k = 0; y < _bounds->height(); ++y) {
        for (uint32_t x = 0; x < _bounds->width(); ++x, ++k) {
            if (m[k] > 0.0f) {
                ++occupied_area;
                occupied_mass += m[k];
                const uint8_t local_material =
                    local_material_map != nullptr
                        ? local_material_map[k]
                        : platec::material::to_index(platec::material::kDefaultType);
                if (material_index_is_continental(local_material)) {
                    ++continental_area;
                    continental_mass += m[k];
                }
            }

            // Set the age of ALL points in this plate to same
            // value. The right thing to do would be to simulate
            // the generation of new oceanic crust as if the plate
            // had been moving to its current direction until all
            // plate's (oceanic) crust receive an age.
            age_map.set(x, y, plate_age & -(m[k] > 0));
            material_map.set(
                x, y,
                local_material_map != nullptr
                    ? local_material_map[k]
                    : platec::material::to_index(platec::material::kDefaultType));
        }
    }

    _continental_area_ratio =
        occupied_area > 0 ? static_cast<float>(continental_area) / static_cast<float>(occupied_area)
                          : 0.0f;
    _continental_mass_ratio =
        occupied_mass > FLT_EPSILON ? continental_mass / occupied_mass : 0.0f;
    _plate_type = _continental_area_ratio <= kOceanicPlateMaxContinentalAreaRatio
        ? PlateType::Oceanic
        : PlateType::Continental;

    _mean_crust =
        occupied_area > 0 ? occupied_mass / static_cast<float>(occupied_area) : 0.0f;
    _buoyancy = plate_buoyancy_from_metrics(_plate_type, _continental_area_ratio,
                                            _continental_mass_ratio, _mean_crust);

    Segments* segments = new Segments(plate_area);
    _segments = segments;
    _mySegmentCreator = new MySegmentCreator(*_bounds, _segments, map, _worldDimension);
    segments->setSegmentCreator(_mySegmentCreator);
    segments->setBounds(_bounds);
}

plate::~plate() {
    delete _mySegmentCreator;
    delete _segments;
    delete _bounds;
}

uint32_t plate::addCollision(uint32_t wx, uint32_t wy) {
    ISegmentData& seg = getContinentAt(wx, wy);
    seg.incCollCount();
    return seg.area();
}

void plate::addCrustByCollision(uint32_t x, uint32_t y, float z, uint32_t time,
                                ContinentId activeContinent, int material_index) {
    // Add crust. Extend plate if necessary.
    setCrust(x, y, getCrust(x, y) + z, time, material_index);

    uint32_t index = _bounds->getMapIndex(&x, &y);
    if (index == BAD_INDEX) {
        ASSERT(false, "BAD map index found");
        return;
    }
    _segments->setId(index, activeContinent);

    ISegmentData& data = (*_segments)[activeContinent];
    data.incArea();
    data.enlarge_to_contain(x, y);
}

void plate::addCrustBySubduction(uint32_t x, uint32_t y, float z, uint32_t t, float dx, float dy,
                                 int material_index) {
    const uint32_t world_x = x;
    const uint32_t world_y = y;
    const uint32_t impact_index = _bounds->getValidMapIndex(&x, &y);

    auto deposit_at = [&](uint32_t index, float amount) {
        const float old_height = map[index];
        if (old_height <= 0.0f || amount <= 0.0f) {
            return false;
        }

        const float combined = old_height + amount;
        age_map[index] = combined > FLT_EPSILON
            ? static_cast<uint32_t>(
                  (old_height * static_cast<float>(age_map[index]) +
                   amount * static_cast<float>(t)) /
                  combined)
            : t;
        map[index] = combined;
        if (material_index >= 0) {
            const uint8_t next_material = static_cast<uint8_t>(material_index);
            if (!(material_index_is_continental(material_map[index]) &&
                  !material_index_is_continental(next_material))) {
                material_map[index] = next_material;
            }
        }
        _mass.incMass(amount);
        return true;
    };

    Platec::FloatVector local_velocity = surfaceVelocityAt(world_x, world_y);
    const bool has_guided_direction =
        std::isfinite(dx) && std::isfinite(dy) &&
        (std::fabs(dx) > FLT_EPSILON || std::fabs(dy) > FLT_EPSILON);
    if (!has_guided_direction) {
        dx = local_velocity.x();
        dy = local_velocity.y();
    } else {
        const float dot = local_velocity.dotProduct(Platec::FloatVector(dx, dy));
        if (dot > 0.0f && local_velocity.normalize() > 0.0f) {
            dx -= local_velocity.x();
            dy -= local_velocity.y();
        }
    }

    float direction_length = std::sqrt(dx * dx + dy * dy);
    if (!std::isfinite(direction_length) || direction_length <= FLT_EPSILON) {
        deposit_at(impact_index, z);
        return;
    }

    dx /= direction_length;
    dy /= direction_length;
    const float lateral_x = -dy;
    const float lateral_y = dx;

    struct DepositTarget {
        uint32_t index;
        float weight;
    };
    std::vector<DepositTarget> targets;
    targets.reserve(32);
    float weight_sum = 0.0f;

    for (int step = 1; step <= 20; ++step) {
        const float along = static_cast<float>(step);
        const float arc_wave = std::max(0.0f, subduction_wavelet(along - 3.0f, 2.0f, 5.0f));
        if (arc_wave <= FLT_EPSILON) {
            continue;
        }

        const int lateral_limit = std::max(1, step / 5);
        for (int lateral_step = -lateral_limit; lateral_step <= lateral_limit; ++lateral_step) {
            const float lateral = static_cast<float>(lateral_step);
            const float lateral_weight =
                std::exp(-(lateral * lateral) / (2.0f * 1.7f * 1.7f));
            const float fx = static_cast<float>(x) + dx * along + lateral_x * lateral;
            const float fy = static_cast<float>(y) + dy * along + lateral_y * lateral;
            if (!_bounds->isInLimits(fx, fy)) {
                continue;
            }

            const uint32_t candidate_index =
                _bounds->index(static_cast<uint32_t>(fx), static_cast<uint32_t>(fy));
            if (map[candidate_index] <= 0.0f) {
                continue;
            }

            const float weight = arc_wave * lateral_weight;
            targets.push_back({candidate_index, weight});
            weight_sum += weight;
        }
    }

    if (weight_sum <= FLT_EPSILON) {
        deposit_at(impact_index, z);
        return;
    }

    for (const DepositTarget& target : targets) {
        deposit_at(target.index, z * (target.weight / weight_sum));
    }
}

float plate::aggregateCrust(plate* p, uint32_t wx, uint32_t wy) {
    uint32_t lx = wx, ly = wy;
    const uint32_t index = _bounds->getMapIndex(&lx, &ly);
    if (index == BAD_INDEX) {
        ASSERT(false, "BAD map index found");
        return 0;
    }

    const ContinentId seg_id = _segments->id(index);
    if (seg_id >= _segments->size()) {
        return 0;
    }

    // This check forces the caller to do things in proper order!
    //
    // Usually continents collide at several locations simultaneously.
    // Thus if this segment that is being merged now is removed from
    // segmentation bookkeeping, then the next point of collision that is
    // processed during the same iteration step would cause the test
    // below to be true and system would experience a premature abort.
    //
    // Therefore, segmentation bookkeeping is left intact. It doesn't
    // cause significant problems because all crust is cleared and empty
    // points are not processed at all. (Test on (seg_id >= seg_data.size()) removed)

    // One continent may have many points of collision. If one of them
    // causes continent to aggregate then all successive collisions and
    // attempts of aggregation would necessarily change nothing at all,
    // because the continent was removed from this plate earlier!
    if ((*_segments)[seg_id].isEmpty()) {
        return 0; // Do not process empty continents.
    }

    ContinentId activeContinent = p->selectCollisionSegment(wx, wy);
    if (activeContinent >= p->_segments->size()) {
        return 0;
    }

    // Wrap coordinates around world edges to safeguard subtractions.
    wx += _worldDimension.getWidth();
    wy += _worldDimension.getHeight();

    // Aggregating segment [%u, %u]x[%u, %u] vs. [%u, %u]@[%u, %u]\n",
    //      seg_data[seg_id].x0, seg_data[seg_id].y0,
    //      seg_data[seg_id].x1, seg_data[seg_id].y1,
    //      _dimension.getWidth(), _dimension.getHeight(), lx, ly);

    float old_mass = _mass.getMass();

    // Add all of the collided continent's crust to destination plate.
    for (uint32_t y = (*_segments)[seg_id].getTop(); y <= (*_segments)[seg_id].getBottom(); ++y) {
        for (uint32_t x = (*_segments)[seg_id].getLeft(); x <= (*_segments)[seg_id].getRight();
             ++x) {
            const uint32_t i = y * _bounds->width() + x;
            if ((_segments->id(i) == seg_id) && (map[i] > 0)) {
                p->addCrustByCollision(wx + x - lx, wy + y - ly, map[i], age_map[i],
                                       activeContinent, material_map[i]);

                _mass.incMass(-1.0f * map[i]);
                map[i] = 0.0f;
            }
        }
    }

    (*_segments)[seg_id].markNonExistent(); // Mark segment as non-existent
    return old_mass - _mass.getMass();
}

void plate::applyFriction(float deformed_mass) {
    // Remove the energy that deformation consumed from plate's kinetic
    // energy: F - dF = ma - dF => a = dF/m.
    if (!_mass.null()) {
        const float friction_scale = 1.0f / std::max(0.05f, _movement_energy);
        _movement.applyFriction(deformed_mass, _mass.getMass(), friction_scale);
    }
}

void plate::collide(plate& p, float coll_mass) {
    if (!_mass.null() && coll_mass > 0) {
        _movement.collide(*this, p, coll_mass);
    }
}

void plate::calculateCrust(uint32_t x, uint32_t y, uint32_t index, float& w_crust, float& e_crust,
                           float& n_crust, float& s_crust, uint32_t& w, uint32_t& e, uint32_t& n,
                           uint32_t& s) {
    ::calculateCrust(x, y, index, w_crust, e_crust, n_crust, s_crust, w, e, n, s, _worldDimension,
                     map, _bounds->width(), _bounds->height());
}

void plate::findRiverSources(float lower_bound, vector<uint32_t>* sources) {
    const uint32_t bounds_height = _bounds->height();
    const uint32_t bounds_width = _bounds->width();

    // Find all tops.
    for (uint32_t y = 0; y < bounds_height; ++y) {
        const uint32_t y_width = y * bounds_width;
        for (uint32_t x = 0; x < bounds_width; ++x) {
            const uint32_t index = y_width + x;

            if (map[index] < lower_bound) {
                continue;
            }

            float w_crust, e_crust, n_crust, s_crust;
            uint32_t w, e, n, s;
            calculateCrust(x, y, index, w_crust, e_crust, n_crust, s_crust, w, e, n, s);

            // This location is either at the edge of the plate or it is not the
            // tallest of its neightbours. Don't start a river from here.
            if (w_crust * e_crust * n_crust * s_crust == 0) {
                continue;
            }

            sources->push_back(index);
        }
    }
}

void plate::flowRivers(float lower_bound, const std::vector<float>& river_strength,
                       vector<uint32_t>* sources, HeightMap& tmp) {
    const uint32_t bounds_area = _bounds->area();
    vector<uint32_t> sinks_data;
    vector<uint32_t>* sinks = &sinks_data;

    static vector<bool> s_flowDone;
    if (s_flowDone.size() < bounds_area) {
        s_flowDone.resize(bounds_area);
    }
    fill(s_flowDone.begin(), s_flowDone.begin() + bounds_area, false);

    // From each top, start flowing water along the steepest slope.
    while (!sources->empty()) {
        while (!sources->empty()) {
            const uint32_t index = sources->back();
            const uint32_t y = index / _bounds->width();
            const uint32_t x = index - y * _bounds->width();

            sources->pop_back();

            if (map[index] < lower_bound) {
                continue;
            }

            float w_crust, e_crust, n_crust, s_crust;
            uint32_t w, e, n, s;
            calculateCrust(x, y, index, w_crust, e_crust, n_crust, s_crust, w, e, n, s);

            // If this is the lowest part of its neighbourhood, stop.
            if (w_crust + e_crust + n_crust + s_crust == 0) {
                continue;
            }

            w_crust += (w_crust == 0) * map[index];
            e_crust += (e_crust == 0) * map[index];
            n_crust += (n_crust == 0) * map[index];
            s_crust += (s_crust == 0) * map[index];

            // Find lowest neighbour.
            float lowest_crust = w_crust;
            uint32_t dest = index - 1;

            if (e_crust < lowest_crust) {
                lowest_crust = e_crust;
                dest = index + 1;
            }

            if (n_crust < lowest_crust) {
                lowest_crust = n_crust;
                dest = index - _bounds->width();
            }

            if (s_crust < lowest_crust) {
                lowest_crust = s_crust;
                dest = index + _bounds->width();
            }

            // if it's not handled yet, add it as new sink.
            if (dest < _bounds->area() && !s_flowDone[dest]) {
                sinks->push_back(dest);
                s_flowDone[dest] = true;
            }

            // Erode this location with the water flow.
            const float strength =
                index < river_strength.size() ? river_strength[index] : _erosion_strength;
            tmp[index] -= (tmp[index] - lower_bound) * kRiverErosionStrength * strength;
        }

        vector<uint32_t>* v_tmp = sources;
        sources = sinks;
        sinks = v_tmp;
        sinks->clear();
    }
}

void plate::erode(float lower_bound, uint32_t current_iteration) {
    if (_erosion_strength <= 0.0f) {
        return;
    }

    const uint32_t bounds_width = _bounds->width();
    const uint32_t bounds_height = _bounds->height();
    std::vector<float> river_strength(_bounds->area(), 0.0f);
    std::vector<float> redistribution_strength(_bounds->area(), 0.0f);

    for (uint32_t y = 0; y < bounds_height; ++y) {
        const uint32_t world_y = _worldDimension.yMod(_bounds->topAsUint() + y);
        const float rainfall = rainfall_from_latitude(world_y, _worldDimension.getHeight());
        for (uint32_t x = 0; x < bounds_width; ++x) {
            const uint32_t index = y * bounds_width + x;
            if (map[index] < lower_bound) {
                continue;
            }

            float w_crust, e_crust, n_crust, s_crust;
            uint32_t w, e, n, s;
            calculateCrust(x, y, index, w_crust, e_crust, n_crust, s_crust, w, e, n, s);

            const float slope = normalized_slope(map[index], w_crust, e_crust, n_crust, s_crust);
            const float wind_variation =
                1.0f + (_randsource.next_float() * 2.0f - 1.0f) * kWindRandomVariation;
            const float wind = clamp_unit(slope * wind_variation);
            const float uplift =
                uplift_compensation(map[index], lower_bound, current_iteration, age_map[index]);

            river_strength[index] = std::min(
                2.5f, _erosion_strength * (0.30f + 0.70f * rainfall) * uplift);
            redistribution_strength[index] = std::min(
                0.45f,
                _erosion_strength *
                    (kRainfallRedistributionStrength * rainfall +
                     kWindRedistributionStrength * wind) *
                    uplift);
        }
    }

    vector<uint32_t> sources_data;
    vector<uint32_t>* sources = &sources_data;
    const float erosion_floor = lower_bound;

    HeightMap tmpHm(map);
    findRiverSources(erosion_floor, sources);
    flowRivers(erosion_floor, river_strength, sources, tmpHm);
    for (uint32_t i = 0; i < _bounds->area(); ++i) {
        if (tmpHm[i] < 0.0f) {
            tmpHm[i] = 0.0f;
        }
    }

    map = tmpHm;
    tmpHm.set_all(0.0f);
    for (uint32_t y = 0; y < _bounds->height(); ++y) {
        for (uint32_t x = 0; x < _bounds->width(); ++x) {
            const uint32_t index = y * _bounds->width() + x;
            tmpHm[index] += map[index]; // Careful not to overwrite earlier amounts.

            if (map[index] < erosion_floor)
                continue;

            const float local_strength = redistribution_strength[index];
            if (local_strength <= FLT_EPSILON) {
                continue;
            }

            float w_crust, e_crust, n_crust, s_crust;
            uint32_t w, e, n, s;
            calculateCrust(x, y, index, w_crust, e_crust, n_crust, s_crust, w, e, n, s);

            // This location has no neighbours (ARTIFACT!) or it is the lowest
            // part of its area. In either case the work here is done.
            if (w_crust + e_crust + n_crust + s_crust == 0)
                continue;

            // The steeper the slope, the more water flows along it.
            // The more downhill (sources), the more water flows to here.
            // 1+1+10 = 12, avg = 4, stdev = sqrt((3*3+3*3+6*6)/3) = 4.2, var = 18,
            //  1*1+1*1+10*10 = 102, 102/4.2=24
            // 1+4+7 = 12, avg = 4, stdev = sqrt((3*3+0*0+3*3)/3) = 2.4, var = 6,
            //  1*1+4*4+7*7 = 66, 66/2.4 = 27
            // 4+4+4 = 12, avg = 4, stdev = sqrt((0*0+0*0+0*0)/3) = 0, var = 0,
            //  4*4+4*4+4*4 = 48, 48/0 = inf -> 48
            // If there's a source slope of height X then it will always cause
            // water erosion of amount Y. Then again from one spot only so much
            // water can flow.
            // Thus, the calculated non-linear flow value for this location is
            // multiplied by the "water erosion" constant.
            // The result is max(result, 1.0). New height of this location could
            // be e.g. h_lowest + (1 - 1 / result) * (h_0 - h_lowest).

            // Calculate the difference in height between this point and its
            // nbours that are lower than this point.
            float w_diff = map[index] - w_crust;
            float e_diff = map[index] - e_crust;
            float n_diff = map[index] - n_crust;
            float s_diff = map[index] - s_crust;

            float min_diff = w_diff;
            min_diff -= (min_diff - e_diff) * (e_diff < min_diff);
            min_diff -= (min_diff - n_diff) * (n_diff < min_diff);
            min_diff -= (min_diff - s_diff) * (s_diff < min_diff);

            // Calculate the sum of difference between lower neighbours and
            // the TALLEST lower neighbour.
            float diff_sum =
                (w_diff - min_diff) * (w_crust > 0) + (e_diff - min_diff) * (e_crust > 0) +
                (n_diff - min_diff) * (n_crust > 0) + (s_diff - min_diff) * (s_crust > 0);

            // Erosion difference sum is negative!
            ASSERT(diff_sum >= 0, "Difference sum must be positive");

            if (diff_sum < min_diff) {
                // There's NOT enough room in neighbours to contain all the
                // crust from this peak so that it would be as tall as its
                // tallest lower neighbour. Thus first step is make ALL
                // lower neighbours and this point equally tall.
                tmpHm[w] += (w_diff - min_diff) * (w_crust > 0) * local_strength;
                tmpHm[e] += (e_diff - min_diff) * (e_crust > 0) * local_strength;
                tmpHm[n] += (n_diff - min_diff) * (n_crust > 0) * local_strength;
                tmpHm[s] += (s_diff - min_diff) * (s_crust > 0) * local_strength;
                tmpHm[index] -= min_diff * local_strength;

                min_diff = (min_diff - diff_sum) * local_strength;

                // Spread the remaining crust equally among all lower nbours.
                min_diff /= 1 + (w_crust > 0) + (e_crust > 0) + (n_crust > 0) + (s_crust > 0);

                tmpHm[w] += min_diff * (w_crust > 0);
                tmpHm[e] += min_diff * (e_crust > 0);
                tmpHm[n] += min_diff * (n_crust > 0);
                tmpHm[s] += min_diff * (s_crust > 0);
                tmpHm[index] += min_diff;
            } else {
                float unit = (min_diff * local_strength) / diff_sum;

                // Remove all crust from this location making it as tall as
                // its tallest lower neighbour.
                tmpHm[index] -= min_diff * local_strength;

                // Spread all removed crust among all other lower neighbours.
                tmpHm[w] += unit * (w_diff - min_diff) * (w_crust > 0);
                tmpHm[e] += unit * (e_diff - min_diff) * (e_crust > 0);
                tmpHm[n] += unit * (n_diff - min_diff) * (n_crust > 0);
                tmpHm[s] += unit * (s_diff - min_diff) * (s_crust > 0);
            }
        }
    }

    // Clamp all heightmap values to prevent negative mass from floating point errors
    // This is a safety measure for Issue #30
    for (uint32_t i = 0; i < _bounds->area(); ++i) {
        if (tmpHm[i] < 0.0f) {
            tmpHm[i] = 0.0f;
        }
    }

    map = tmpHm;
    _mass = MassBuilder(map.raw_data(), Dimension(_bounds->width(), _bounds->height())).build();
}

void plate::getCollisionInfo(uint32_t wx, uint32_t wy, uint32_t* count, float* ratio) const {
    uint32_t lx = wx;
    uint32_t ly = wy;
    const uint32_t index = _bounds->getMapIndex(&lx, &ly);
    if (index == BAD_INDEX || map[index] <= 0.0f) {
        *count = 0;
        *ratio = 0.0f;
        return;
    }

    const ISegmentData& seg = getContinentAt(wx, wy);

    *count = seg.collCount();
    *ratio = (float)seg.collCount() / (float)(1 + seg.area()); // +1 avoids DIV with zero.
}

uint32_t plate::getContinentArea(uint32_t wx, uint32_t wy) const {
    const uint32_t index = _bounds->getMapIndex(&wx, &wy);
    if (index == BAD_INDEX) {
        ASSERT(false, "BAD map index found");
        return 0;
    }
    ASSERT(_segments->id(index) < _segments->size(), "Segment index invalid");
    if (_segments->id(index) >= _segments->size()) {
        return 0;
    }
    return (*_segments)[_segments->id(index)].area();
}

float plate::getCrust(uint32_t x, uint32_t y) const {
    const uint32_t index = _bounds->getMapIndex(&x, &y);
    return index != BAD_INDEX ? map[index] : 0;
}

uint32_t plate::getCrustTimestamp(uint32_t x, uint32_t y) const {
    const uint32_t index = _bounds->getMapIndex(&x, &y);
    return index != BAD_INDEX ? age_map[index] : 0;
}

void plate::setCrustTimestamp(uint32_t x, uint32_t y, uint32_t t) {
    const uint32_t index = _bounds->getMapIndex(&x, &y);
    if (index == BAD_INDEX || map[index] <= FLT_EPSILON) {
        return;
    }
    age_map[index] = t;
}

uint8_t plate::getMaterial(uint32_t x, uint32_t y) const {
    const uint32_t index = _bounds->getMapIndex(&x, &y);
    return index != BAD_INDEX ? material_map[index]
                              : platec::material::to_index(platec::material::kDefaultType);
}

float plate::getSubmergedCrust(uint32_t x, uint32_t y) const {
    const uint32_t index = _bounds->getMapIndex(&x, &y);
    return index != BAD_INDEX ? submerged_map[index] : 0.0f;
}

uint32_t plate::getSubmergedSince(uint32_t x, uint32_t y) const {
    const uint32_t index = _bounds->getMapIndex(&x, &y);
    return index != BAD_INDEX ? submerged_since_map[index] : 0U;
}

uint8_t plate::getSubmergedMaterial(uint32_t x, uint32_t y) const {
    const uint32_t index = _bounds->getMapIndex(&x, &y);
    return index != BAD_INDEX ? submerged_material_map[index]
                              : platec::material::to_index(platec::material::kDefaultType);
}

void plate::submergeCrust(uint32_t x, uint32_t y, float z, uint32_t crust_age,
                          uint32_t hidden_since, int material_index) {
    if (z <= 0.0f) {
        return;
    }

    const uint32_t index = _bounds->getValidMapIndex(&x, &y);
    const float previous_hidden = submerged_map[index];
    if (z > previous_hidden + FLT_EPSILON) {
        submerged_map[index] = z;
        submerged_age_map[index] = crust_age;
        submerged_since_map[index] = hidden_since;
        if (material_index >= 0) {
            submerged_material_map[index] = static_cast<uint8_t>(material_index);
        }
        _submerged_mass += z - previous_hidden;
    }
}

float plate::resurfaceSubmergedCrust(uint32_t x, uint32_t y, float fraction, float min_amount,
                                     uint32_t current_iteration) {
    (void)current_iteration;
    const uint32_t index = _bounds->getMapIndex(&x, &y);
    if (index == BAD_INDEX || submerged_map[index] <= FLT_EPSILON) {
        return 0.0f;
    }

    const float resurfaced = std::min(
        submerged_map[index], std::max(min_amount, submerged_map[index] * std::max(0.0f, fraction)));
    if (resurfaced <= FLT_EPSILON) {
        return 0.0f;
    }

    if (map[index] > FLT_EPSILON) {
        const float combined = map[index] + resurfaced;
        age_map[index] = static_cast<uint32_t>(
            (map[index] * static_cast<float>(age_map[index]) +
             resurfaced * static_cast<float>(submerged_age_map[index])) /
            combined);
        map[index] = combined;
    } else {
        map[index] = resurfaced;
        age_map[index] = submerged_age_map[index];
        material_map[index] = submerged_material_map[index];
    }

    _mass.incMass(resurfaced);
    submerged_map[index] -= resurfaced;
    _submerged_mass -= resurfaced;
    if (submerged_map[index] <= FLT_EPSILON) {
        submerged_map[index] = 0.0f;
        submerged_age_map[index] = 0U;
        submerged_since_map[index] = 0U;
        submerged_material_map[index] =
            platec::material::to_index(platec::material::kDefaultType);
    }

    return resurfaced;
}

FloatPoint plate::worldMassCenter() const {
    if (_mass.null()) {
        return FloatPoint(static_cast<float>(_bounds->leftAsUint()),
                          static_cast<float>(_bounds->topAsUint()));
    }

    const float world_x =
        wrap_coordinate(static_cast<float>(_bounds->leftAsUint()) + _mass.getCx(),
                        _worldDimension.getWidth());
    const float world_y =
        wrap_coordinate(static_cast<float>(_bounds->topAsUint()) + _mass.getCy(),
                        _worldDimension.getHeight());
    return FloatPoint(world_x, world_y);
}

Platec::FloatVector plate::surfaceVelocityAt(uint32_t x, uint32_t y) const {
    const Platec::FloatVector linear_velocity = _movement.velocityVector();
    if (_mass.null()) {
        return linear_velocity;
    }

    const FloatPoint center = worldMassCenter();
    const float dx =
        wrapped_delta(static_cast<float>(x), center.getX(), _worldDimension.getWidth());
    const float dy =
        wrapped_delta(static_cast<float>(y), center.getY(), _worldDimension.getHeight());
    const float angular_velocity = _movement.rotationAngle();

    return Platec::FloatVector(linear_velocity.x() - dy * angular_velocity,
                               linear_velocity.y() + dx * angular_velocity);
}

float plate::continentalAreaRatio() const noexcept {
    uint32_t occupied_area = 0U;
    uint32_t continental_area = 0U;
    const uint32_t area = map.area();
    for (uint32_t index = 0; index < area; ++index) {
        if (map[index] <= FLT_EPSILON) {
            continue;
        }
        ++occupied_area;
        if (material_index_is_continental(material_map[index])) {
            ++continental_area;
        }
    }

    return occupied_area > 0U
        ? static_cast<float>(continental_area) / static_cast<float>(occupied_area)
        : _continental_area_ratio;
}

PlateType plate::plateType() const noexcept {
    return isOceanicPlate() ? PlateType::Oceanic : PlateType::Continental;
}

bool plate::isOceanicPlate() const noexcept {
    return continentalAreaRatio() <= kOceanicPlateMaxContinentalAreaRatio;
}

bool plate::isContinentalPlate() const noexcept {
    return !isOceanicPlate();
}

float plate::continentalityScore() const noexcept {
    return 0.65f * _continental_mass_ratio + 0.35f * continentalAreaRatio() +
           0.05f * _mean_crust;
}

void plate::forcePlateType(PlateType type) noexcept {
    _plate_type = type;
    _buoyancy =
        plate_buoyancy_from_metrics(_plate_type, _continental_area_ratio, _continental_mass_ratio,
                                    _mean_crust);
}

void plate::getMap(const float** c, const uint32_t** t, const uint8_t** m) const {
    if (c) {
        *c = map.raw_data();
    }
    if (t) {
        *t = age_map.raw_data();
    }
    if (m) {
        *m = material_map.raw_data();
    }
}

void plate::ensureRotationPadding(uint32_t margin) {
    if (margin == 0) {
        return;
    }

    const uint32_t width = _bounds->width();
    const uint32_t height = _bounds->height();
    if (width >= _worldDimension.getWidth() && height >= _worldDimension.getHeight()) {
        return;
    }

    bool need_left = false;
    bool need_right = false;
    bool need_top = false;
    bool need_bottom = false;

    for (uint32_t y = 0; y < height && !(need_left && need_right); ++y) {
        for (uint32_t x = 0; x < margin && x < width; ++x) {
            if (map[y * width + x] > 0.0f) {
                need_left = true;
                break;
            }
        }
        for (uint32_t x = 0; x < margin && x < width; ++x) {
            if (map[y * width + (width - 1 - x)] > 0.0f) {
                need_right = true;
                break;
            }
        }
    }

    for (uint32_t x = 0; x < width && !(need_top && need_bottom); ++x) {
        for (uint32_t y = 0; y < margin && y < height; ++y) {
            if (map[y * width + x] > 0.0f) {
                need_top = true;
                break;
            }
        }
        for (uint32_t y = 0; y < margin && y < height; ++y) {
            if (map[(height - 1 - y) * width + x] > 0.0f) {
                need_bottom = true;
                break;
            }
        }
    }

    uint32_t available_x = _worldDimension.getWidth() - width;
    uint32_t available_y = _worldDimension.getHeight() - height;
    uint32_t pad_left = need_left && available_x > 0 ? 1U : 0U;
    uint32_t pad_right = need_right && available_x > pad_left ? 1U : 0U;
    uint32_t pad_top = need_top && available_y > 0 ? 1U : 0U;
    uint32_t pad_bottom = need_bottom && available_y > pad_top ? 1U : 0U;

    if (pad_left == 0 && pad_right == 0 && pad_top == 0 && pad_bottom == 0) {
        return;
    }

    _bounds->shift(-1.0f * static_cast<float>(pad_left), -1.0f * static_cast<float>(pad_top));
    _bounds->grow(static_cast<int>(pad_left + pad_right), static_cast<int>(pad_top + pad_bottom));

    HeightMap padded_map(_bounds->width(), _bounds->height());
    AgeMap padded_age(_bounds->width(), _bounds->height());
    MaterialMap padded_material(_bounds->width(), _bounds->height());
    HeightMap padded_submerged(_bounds->width(), _bounds->height());
    AgeMap padded_submerged_age(_bounds->width(), _bounds->height());
    AgeMap padded_submerged_since(_bounds->width(), _bounds->height());
    MaterialMap padded_submerged_material(_bounds->width(), _bounds->height());
    uint32_t* padded_segments = new uint32_t[_bounds->area()];
    padded_map.set_all(0.0f);
    padded_age.set_all(0);
    padded_material.set_all(platec::material::to_index(platec::material::kDefaultType));
    padded_submerged.set_all(0.0f);
    padded_submerged_age.set_all(0U);
    padded_submerged_since.set_all(0U);
    padded_submerged_material.set_all(platec::material::to_index(platec::material::kDefaultType));
    memset(padded_segments, 255, _bounds->area() * sizeof(uint32_t));

    for (uint32_t y = 0; y < height; ++y) {
        const uint32_t dest_index = (pad_top + y) * _bounds->width() + pad_left;
        const uint32_t src_index = y * width;
        memcpy(&padded_map[dest_index], &map[src_index], width * sizeof(float));
        memcpy(&padded_age[dest_index], &age_map[src_index], width * sizeof(uint32_t));
        memcpy(&padded_material[dest_index], &material_map[src_index], width * sizeof(uint8_t));
        memcpy(&padded_submerged[dest_index], &submerged_map[src_index], width * sizeof(float));
        memcpy(&padded_submerged_age[dest_index], &submerged_age_map[src_index],
               width * sizeof(uint32_t));
        memcpy(&padded_submerged_since[dest_index], &submerged_since_map[src_index],
               width * sizeof(uint32_t));
        memcpy(&padded_submerged_material[dest_index], &submerged_material_map[src_index],
               width * sizeof(uint8_t));
        memcpy(&padded_segments[dest_index], &_segments->id(src_index), width * sizeof(uint32_t));
    }

    map = padded_map;
    age_map = padded_age;
    material_map = padded_material;
    submerged_map = padded_submerged;
    submerged_age_map = padded_submerged_age;
    submerged_since_map = padded_submerged_since;
    submerged_material_map = padded_submerged_material;
    _segments->reassign(_bounds->area(), padded_segments);
    _segments->shift(pad_left, pad_top);
}

void plate::rotateCrust(float angle) {
    if (fabsf(angle) < 0.0001f || _mass.null()) {
        return;
    }

    ensureRotationPadding(2);

    const uint32_t width = _bounds->width();
    const uint32_t height = _bounds->height();
    const bool wrap_x = width == _worldDimension.getWidth();
    const bool wrap_y = height == _worldDimension.getHeight();
    HeightMap rotated(width, height);
    AgeMap rotated_age(width, height);
    MaterialMap rotated_material(width, height);
    HeightMap rotated_submerged(width, height);
    AgeMap rotated_submerged_age(width, height);
    AgeMap rotated_submerged_since(width, height);
    MaterialMap rotated_submerged_material(width, height);
    rotated.set_all(0.0f);
    rotated_age.set_all(0);
    rotated_material.set_all(platec::material::to_index(platec::material::kDefaultType));
    rotated_submerged.set_all(0.0f);
    rotated_submerged_age.set_all(0U);
    rotated_submerged_since.set_all(0U);
    rotated_submerged_material.set_all(platec::material::to_index(platec::material::kDefaultType));

    float center_x = _mass.getCx();
    float center_y = _mass.getCy();
    if (!_bounds->isInLimits(center_x, center_y)) {
        center_x = (static_cast<float>(width) - 1.0f) * 0.5f;
        center_y = (static_cast<float>(height) - 1.0f) * 0.5f;
    }

    const float c = cosf(angle);
    const float s = sinf(angle);
    MassBuilder mass_builder;
    for (uint32_t y = 0; y < height; ++y) {
        for (uint32_t x = 0; x < width; ++x) {
            const float dx = static_cast<float>(x) - center_x;
            const float dy = static_cast<float>(y) - center_y;
            const float source_x = dx * c + dy * s + center_x;
            const float source_y = -dx * s + dy * c + center_y;

            if ((!wrap_x && (source_x < 0.0f || source_x > static_cast<float>(width - 1U))) ||
                (!wrap_y && (source_y < 0.0f || source_y > static_cast<float>(height - 1U)))) {
                continue;
            }

            const int rounded_x = static_cast<int>(std::lround(
                wrap_x ? wrap_coordinate(source_x, width) : source_x));
            const int rounded_y = static_cast<int>(std::lround(
                wrap_y ? wrap_coordinate(source_y, height) : source_y));
            const uint32_t sx = wrap_x
                ? wrap_index(rounded_x, width)
                : static_cast<uint32_t>(
                      std::max(0, std::min(static_cast<int>(width - 1U), rounded_x)));
            const uint32_t sy = wrap_y
                ? wrap_index(rounded_y, height)
                : static_cast<uint32_t>(
                      std::max(0, std::min(static_cast<int>(height - 1U), rounded_y)));
            const uint32_t source_index = sy * width + sx;
            const uint32_t dest_index = y * width + x;

            const float crust = map[source_index];
            if (crust > 0.00001f) {
                rotated[dest_index] = crust;
                rotated_age[dest_index] = age_map[source_index];
                rotated_material[dest_index] = material_map[source_index];
                mass_builder.addPoint(x, y, crust);
            }

            const float submerged_crust = submerged_map[source_index];
            if (submerged_crust > 0.00001f) {
                rotated_submerged[dest_index] = submerged_crust;
                rotated_submerged_age[dest_index] = submerged_age_map[source_index];
                rotated_submerged_since[dest_index] = submerged_since_map[source_index];
                rotated_submerged_material[dest_index] = submerged_material_map[source_index];
            }
        }
    }

    map = rotated;
    age_map = rotated_age;
    material_map = rotated_material;
    submerged_map = rotated_submerged;
    submerged_age_map = rotated_submerged_age;
    submerged_since_map = rotated_submerged_since;
    submerged_material_map = rotated_submerged_material;
    _mass = mass_builder.build();
}

void plate::move() {
    _movement.move();
    _pending_crust_rotation += _movement.rotationAngle() * _crust_rotation_strength;
    if (fabsf(_pending_crust_rotation) >= kMinRasterRotationAngle) {
        rotateCrust(_pending_crust_rotation);
        _pending_crust_rotation = 0.0f;
    }

    // Location modulations into range [0..world width/height[ are a have to!
    // If left undone SOMETHING WILL BREAK DOWN SOMEWHERE in the code!

    _bounds->shift(_movement.velocityOnX(), _movement.velocityOnY());
}

void plate::resetSegments() {
    ASSERT(_bounds->area() == _segments->area(), "Segments doesn't have the expected area");
    _segments->reset();
}

void plate::setCrust(uint32_t x, uint32_t y, float z, uint32_t t, int material_index) {
    if (!std::isfinite(z) || z < 0) { // Do not accept invalid or negative values.
        z = 0;
    }

    uint32_t _x = x;
    uint32_t _y = y;
    uint32_t index = _bounds->getMapIndex(&_x, &_y);

    if (index == BAD_INDEX) {
        // Extending plate for nothing!
        if (z <= 0.0f) {
            return;
        }

        const uint32_t ilft = _bounds->leftAsUint();
        const uint32_t itop = _bounds->topAsUint();
        const uint32_t irgt = _bounds->rightAsUintNonInclusive();
        const uint32_t ibtm = _bounds->bottomAsUintNonInclusive();

        _worldDimension.normalize(x, y);

        // Calculate distance of new point from plate edges.
        const uint32_t _lft = ilft - x;
        const uint32_t _rgt = (_worldDimension.getWidth() & -(x < ilft)) + x - irgt;
        const uint32_t _top = itop - y;
        const uint32_t _btm = (_worldDimension.getHeight() & -(y < itop)) + y - ibtm;

        // Set larger of horizontal/vertical distance to zero.
        // A valid distance is NEVER larger than world's side's length!
        uint32_t d_lft = _lft & -(_lft < _rgt) & -(_lft < _worldDimension.getWidth());
        uint32_t d_rgt = _rgt & -(_rgt <= _lft) & -(_rgt < _worldDimension.getWidth());
        uint32_t d_top = _top & -(_top < _btm) & -(_top < _worldDimension.getHeight());
        uint32_t d_btm = _btm & -(_btm <= _top) & -(_btm < _worldDimension.getHeight());

        // Scale all changes to multiple of 8.
        d_lft = ((d_lft > 0) + (d_lft >> 3)) << 3;
        d_rgt = ((d_rgt > 0) + (d_rgt >> 3)) << 3;
        d_top = ((d_top > 0) + (d_top >> 3)) << 3;
        d_btm = ((d_btm > 0) + (d_btm >> 3)) << 3;

        // Make sure plate doesn't grow bigger than the system it's in!
        if (_bounds->width() + d_lft + d_rgt > _worldDimension.getWidth()) {
            d_lft = 0;
            d_rgt = _worldDimension.getWidth() - _bounds->width();
        }

        if (_bounds->height() + d_top + d_btm > _worldDimension.getHeight()) {
            d_top = 0;
            d_btm = _worldDimension.getHeight() - _bounds->height();
        }

        // Index out of bounds, but nowhere to grow!
        if (d_lft + d_rgt + d_top + d_btm == 0) {
            ASSERT(false, "Invalid plate growth deltas");
            return;
        }

        const uint32_t old_width = _bounds->width();
        const uint32_t old_height = _bounds->height();

        _bounds->shift(-1.0f * d_lft, -1.0f * d_top);
        _bounds->grow(d_lft + d_rgt, d_top + d_btm);

        HeightMap tmph = HeightMap(_bounds->width(), _bounds->height());
        AgeMap tmpa = AgeMap(_bounds->width(), _bounds->height());
        MaterialMap tmpm = MaterialMap(_bounds->width(), _bounds->height());
        HeightMap tmpsub = HeightMap(_bounds->width(), _bounds->height());
        AgeMap tmpsub_age = AgeMap(_bounds->width(), _bounds->height());
        AgeMap tmpsub_since = AgeMap(_bounds->width(), _bounds->height());
        MaterialMap tmpsub_mat = MaterialMap(_bounds->width(), _bounds->height());
        uint32_t* tmps = new uint32_t[_bounds->area()];
        tmph.set_all(0);
        tmpa.set_all(0);
        tmpm.set_all(platec::material::to_index(platec::material::kDefaultType));
        tmpsub.set_all(0.0f);
        tmpsub_age.set_all(0U);
        tmpsub_since.set_all(0U);
        tmpsub_mat.set_all(platec::material::to_index(platec::material::kDefaultType));
        memset(tmps, 255, _bounds->area() * sizeof(uint32_t));

        // copy old plate into new.
        for (uint32_t j = 0; j < old_height; ++j) {
            const uint32_t dest_i = (d_top + j) * _bounds->width() + d_lft;
            const uint32_t src_i = j * old_width;
            memcpy(&tmph[dest_i], &map[src_i], old_width * sizeof(float));
            memcpy(&tmpa[dest_i], &age_map[src_i], old_width * sizeof(uint32_t));
            memcpy(&tmpm[dest_i], &material_map[src_i], old_width * sizeof(uint8_t));
            memcpy(&tmpsub[dest_i], &submerged_map[src_i], old_width * sizeof(float));
            memcpy(&tmpsub_age[dest_i], &submerged_age_map[src_i], old_width * sizeof(uint32_t));
            memcpy(&tmpsub_since[dest_i], &submerged_since_map[src_i], old_width * sizeof(uint32_t));
            memcpy(&tmpsub_mat[dest_i], &submerged_material_map[src_i], old_width * sizeof(uint8_t));
            memcpy(&tmps[dest_i], &_segments->id(src_i), old_width * sizeof(uint32_t));
        }

        map = tmph;
        age_map = tmpa;
        material_map = tmpm;
        submerged_map = tmpsub;
        submerged_age_map = tmpsub_age;
        submerged_since_map = tmpsub_since;
        submerged_material_map = tmpsub_mat;
        _segments->reassign(_bounds->area(), tmps);

        // Shift all segment data to match new coordinates.
        _segments->shift(d_lft, d_top);

        _x = x, _y = y;
        index = _bounds->getMapIndex(&_x, &_y);
        if (index == BAD_INDEX) {
            ASSERT(false, "BAD map index found");
            return;
        }
    }

    // Clamp to prevent floating point precision errors (Issue #30)
    if (z < 0.0f) {
        z = 0.0f;
    }

    // Update crust age without evaluating the empty->empty 0/0 case.
    const float old_height = map[index];
    uint32_t next_age = 0U;
    if (z > 0.0f) {
        if (old_height > 0.0f) {
            const float total_height = old_height + z;
            next_age = total_height > FLT_EPSILON
                ? static_cast<uint32_t>(
                      (old_height * static_cast<float>(age_map[index]) + z * static_cast<float>(t)) /
                      total_height)
                : t;
        } else {
            next_age = t;
        }
    }

    _mass.incMass(-1.0f * old_height);
    _mass.incMass(z); // Update mass counter.
    map[index] = z;   // Set new crust height to desired location.
    age_map[index] = next_age;
    if (material_index >= 0 && z > 0.0f) {
        const uint8_t next_material = static_cast<uint8_t>(material_index);
        if (!(old_height > FLT_EPSILON && material_index_is_continental(material_map[index]) &&
              !material_index_is_continental(next_material))) {
            material_map[index] = next_material;
        }
    }
}

ContinentId plate::selectCollisionSegment(uint32_t coll_x, uint32_t coll_y) {
    uint32_t index = _bounds->getMapIndex(&coll_x, &coll_y);
    if (index == BAD_INDEX) {
        ASSERT(false, "BAD map index found");
        return _segments->size();
    }
    ContinentId activeContinent = _segments->id(index);
    return activeContinent;
}

///////////////////////////////////////////////////////////////////////////////
/// Private methods ///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

uint32_t plate::createSegment(uint32_t x, uint32_t y) throw() {
    return _mySegmentCreator->createSegment(x, y);
}

ISegmentData& plate::getContinentAt(int x, int y) {
    return (*_segments)[_segments->getContinentAt(x, y)];
}

const ISegmentData& plate::getContinentAt(int x, int y) const {
    return (*_segments)[_segments->getContinentAt(x, y)];
}
