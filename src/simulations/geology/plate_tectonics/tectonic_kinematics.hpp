#ifndef TECTONIC_KINEMATICS_HPP
#define TECTONIC_KINEMATICS_HPP

#include "tectonic_contract.hpp"
#include <cmath>

namespace platec::kinematics {

struct Velocity { float x; float y; };

inline float wrapped_offset(float point, float center, uint32_t period) {
    float delta = point - center;
    const float size = static_cast<float>(period);
    if (delta > size * 0.5f) delta -= size;
    else if (delta < -size * 0.5f) delta += size;
    return delta;
}

// Native planar velocities are cells/update; angular velocity is radians/update.
// Evaluate both plates at the same boundary point before taking their difference.
inline Velocity local_velocity(const contract::PlateKinematics& plate, float x, float y,
                               uint32_t width, uint32_t height) {
    const float dx = wrapped_offset(x, plate.mass_center_x, width);
    const float dy = wrapped_offset(y, plate.mass_center_y, height);
    return {plate.linear_velocity_x - dy * plate.angular_velocity,
            plate.linear_velocity_y + dx * plate.angular_velocity};
}

} // namespace platec::kinematics
#endif
