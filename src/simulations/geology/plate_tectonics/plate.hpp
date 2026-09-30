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

#ifndef PLATE_HPP
#define PLATE_HPP

#include "bounds.hpp"
#include "heightmap.hpp"
#include "mass.hpp"
#include "movement.hpp"
#include "rectangle.hpp"
#include "simplerandom.hpp"
#include "utils.hpp"
#include <cmath> // sin, cos
#include <cstdint>
#include <vector>

enum class PlateType {
    Oceanic,
    Continental,
};

class IPlate : public IMass, public IMovement {
  public:
    ~IPlate() override = default;
    virtual FloatPoint worldMassCenter() const = 0;
};

// Per-plate geometry and kinematic cache. Material evolution belongs to columns::Model.
class plate : public IPlate {
  public:
    /// Initializes plate with the supplied height map.
    ///
    /// @param  m              Pointer to array to height map of terrain.
    /// @param  w              Width of height map in pixels.
    /// @param  h              Height of height map in pixels.
    /// @param  _x             X of height map's left-top corner on world map.
    /// @param  _y             Y of height map's left-top corner on world map.
    /// @param  worldDimension Dimension of world map's either side in pixels.
    plate(long seed, float* m, uint32_t w, uint32_t h, uint32_t _x, uint32_t _y, uint32_t plate_age,
          WorldDimension worldDimension, const uint8_t* local_material_map = nullptr,
          float rotation_strength = 1.0f);

    ~plate() override;

    /// Get the amount of plate's crustal material at some location.
    ///
    /// @param  x   Offset on the global world map along X axis.
    /// @param  y   Offset on the global world map along Y axis.
    /// @return     Amount of crust at requested location.
    float getCrust(uint32_t x, uint32_t y) const;

    /// Get the timestamp of plate's crustal material at some location.
    ///
    /// @param  x   Offset on the global world map along X axis.
    /// @param  y   Offset on the global world map along Y axis.
    /// @return     Timestamp of creation of crust at the location.
    ///                     Zero is returned if location contains no crust.
    uint32_t getCrustTimestamp(uint32_t x, uint32_t y) const;
    void setCrustTimestamp(uint32_t x, uint32_t y, uint32_t t);
    uint8_t getMaterial(uint32_t x, uint32_t y) const;
    /// Get pointers to plate's data.
    ///
    /// @param  c   Adress of crust height map is stored here.
    /// @param  t   Adress of crust timestamp map is stored here.
    void getMap(const float** c, const uint32_t** t, const uint8_t** m = nullptr) const;

    /// Set the amount of plate's crustal material at some location.
    ///
    /// If amount of crust to be set is negative, it'll be set to zero.
    ///
    /// @param  x   Offset on the global world map along X axis.
    /// @param  y   Offset on the global world map along Y axis.
    /// @param  z   Amount of crust at given location.
    /// @param  t   Time of creation of new crust.
    void setCrust(uint32_t x, uint32_t y, float z, uint32_t t, int material_index = -1);

    float getMass() const throw() override {
        return _mass.getMass();
    }
    float getMomentum() const throw() {
        return _movement.momentum(_mass);
    }
    uint32_t getHeight() const throw() {
        return _bounds->height();
    }
    uint32_t getLeftAsUint() const throw() {
        return _bounds->leftAsUint();
    }
    uint32_t getTopAsUint() const throw() {
        return _bounds->topAsUint();
    }
    float getVelocity() const throw() {
        return _movement.getVelocity();
    }
    float getAngularVelocity() const throw() {
        return _movement.rotationAngle();
    }
    Tectonics::FloatVector linearVelocityVector() const {
        return _movement.velocityVector();
    }
    FloatPoint worldMassCenter() const override;
    Tectonics::FloatVector surfaceVelocityAt(uint32_t x, uint32_t y) const;
    void setKinematics(const Tectonics::FloatVector& linear_velocity, float angular_velocity) {
        _movement.setKinematics(linear_velocity, angular_velocity);
    }
    PlateType plateType() const noexcept;
    bool isOceanicPlate() const noexcept;
    bool isContinentalPlate() const noexcept;
    float continentalAreaRatio() const noexcept;
    float buoyancy() const noexcept {
        return _buoyancy;
    }
    float continentalityScore() const noexcept;
    void forcePlateType(PlateType type) noexcept;

    Tectonics::FloatVector velocityUnitVector() const override {
        return _movement.velocityUnitVector();
    }

    /// @Deprecated, use velocityUnitVector instead
    float getVelX() const throw() {
        return _movement.velX();
    }
    /// @Deprecated, use velocityUnitVector instead
    float getVelY() const throw() {
        return _movement.velY();
    }

    uint32_t getWidth() const throw() {
        return _bounds->width();
    }
    bool isEmpty() const throw() {
        return _mass.null();
    }
    float getCx() const {
        return _mass.getCx();
    }
    float getCy() const {
        return _mass.getCy();
    }
    FloatPoint massCenter() const override {
        return _mass.massCenter();
    }

    void decImpulse(const Tectonics::FloatVector& delta) override {
        _movement.decDx(delta.x());
        _movement.decDy(delta.y());
    }

    // @Deprecated, use decImpulse instead
    void decDx(float delta) {
        _movement.decDx(delta);
    }
    void decDy(float delta) {
        _movement.decDy(delta);
    }

    // Visible for testing
    void injectBounds(IBounds* bounds) {
        delete _bounds;
        _bounds = bounds;
    }

  private:
    const WorldDimension _worldDimension;
    SimpleRandom _randsource;
    HeightMap map;  ///< Bitmap of plate's structure/height.
    AgeMap age_map; ///< Bitmap of plate's soil's age: timestamp of creation.
    MaterialMap material_map; ///< Bitmap of plate material indices.
    IBounds* _bounds;
    Mass _mass;
    Movement _movement;
    PlateType _plate_type;
    float _buoyancy;
    float _continental_area_ratio;
    float _continental_mass_ratio;
    float _mean_crust;
};

#endif
