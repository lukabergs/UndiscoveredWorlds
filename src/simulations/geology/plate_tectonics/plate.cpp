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

#include "simulations/geology/plate_tectonics/heightmap.hpp"
#include "simulations/geology/plate_tectonics/material.hpp"
#include "simulations/geology/plate_tectonics/plate.hpp"
#include "simulations/geology/plate_tectonics/rectangle.hpp"
#include "simulations/geology/plate_tectonics/utils.hpp"

using namespace std;

namespace {

static const float kOceanicPlateMaxContinentalAreaRatio = 0.10f;

bool material_index_is_continental(uint8_t material_index) {
    return tectonics::material::from_index(material_index) != tectonics::material::Type::Basalt;
}

// Identity only; never used for randomness or ordering of simulation work.
float plate_buoyancy_from_metrics(PlateType type, float continental_area_ratio,
                                  float continental_mass_ratio, float mean_crust) {
    return type == PlateType::Continental
        ? 1.10f + 0.85f * continental_area_ratio + 0.35f * continental_mass_ratio +
              0.10f * mean_crust
        : 0.25f + 0.20f * continental_area_ratio + 0.55f * mean_crust;
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

} // namespace

plate::plate(long seed, float* m, uint32_t w, uint32_t h, uint32_t _x, uint32_t _y,
             uint32_t plate_age, WorldDimension worldDimension, const uint8_t* local_material_map,
             float rotation_strength)
    : _worldDimension(worldDimension), _randsource(seed), map(m, w, h), age_map(w, h),
      material_map(w, h),
      _bounds(nullptr), _mass(MassBuilder(m, Dimension(w, h)).build()),
      _movement(_randsource, worldDimension, rotation_strength), _plate_type(PlateType::Oceanic),
      _buoyancy(0.0f), _continental_area_ratio(0.0f),
      _continental_mass_ratio(0.0f), _mean_crust(0.0f) {
    uint32_t occupied_area = 0;
    uint32_t continental_area = 0;
    float occupied_mass = 0.0f;
    float continental_mass = 0.0f;

    _bounds = new Bounds(worldDimension, FloatPoint(static_cast<float>(_x), static_cast<float>(_y)),
                         Dimension(w, h));

    uint32_t k;
    for (uint32_t y = k = 0; y < _bounds->height(); ++y) {
        for (uint32_t x = 0; x < _bounds->width(); ++x, ++k) {
            if (m[k] > 0.0f) {
                ++occupied_area;
                occupied_mass += m[k];
                const uint8_t local_material =
                    local_material_map != nullptr
                        ? local_material_map[k]
                        : tectonics::material::to_index(tectonics::material::kDefaultType);
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
                    : tectonics::material::to_index(tectonics::material::kDefaultType));
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

}

plate::~plate() {
    delete _bounds;
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
                              : tectonics::material::to_index(tectonics::material::kDefaultType);
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

Tectonics::FloatVector plate::surfaceVelocityAt(uint32_t x, uint32_t y) const {
    const Tectonics::FloatVector linear_velocity = _movement.velocityVector();
    if (_mass.null()) {
        return linear_velocity;
    }

    const FloatPoint center = worldMassCenter();
    const float dx =
        wrapped_delta(static_cast<float>(x), center.getX(), _worldDimension.getWidth());
    const float dy =
        wrapped_delta(static_cast<float>(y), center.getY(), _worldDimension.getHeight());
    const float angular_velocity = _movement.rotationAngle();

    return Tectonics::FloatVector(linear_velocity.x() - dy * angular_velocity,
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
        tmph.set_all(0);
        tmpa.set_all(0);
        tmpm.set_all(tectonics::material::to_index(tectonics::material::kDefaultType));

        // copy old plate into new.
        for (uint32_t j = 0; j < old_height; ++j) {
            const uint32_t dest_i = (d_top + j) * _bounds->width() + d_lft;
            const uint32_t src_i = j * old_width;
            memcpy(&tmph[dest_i], &map[src_i], old_width * sizeof(float));
            memcpy(&tmpa[dest_i], &age_map[src_i], old_width * sizeof(uint32_t));
            memcpy(&tmpm[dest_i], &material_map[src_i], old_width * sizeof(uint8_t));
        }

        map = tmph;
        age_map = tmpa;
        material_map = tmpm;


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

///////////////////////////////////////////////////////////////////////////////
/// Private methods ///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

