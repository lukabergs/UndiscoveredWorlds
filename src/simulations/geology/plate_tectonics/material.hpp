#ifndef MATERIAL_HPP
#define MATERIAL_HPP

#include <array>
#include <cstdint>

namespace platec::material {

enum class Type : uint8_t {
    Granite = 0,
    Basalt = 1,
    Sedimentary = 2,
    Metamorphic = 3,
};

struct Properties {
    const char* name;
    uint8_t r;
    uint8_t g;
    uint8_t b;
    float density_kg_m3;
    float tensile_strength_pa;
};

inline constexpr std::array<Properties, 4> kProperties = {{
    {"granite", 178, 176, 182, 2700.0f, 2.0e8f},
    {"basalt", 72, 80, 96, 3000.0f, 1.7e8f},
    {"sedimentary", 206, 176, 122, 2350.0f, 1.0e8f},
    {"metamorphic", 112, 124, 138, 2850.0f, 2.2e8f},
}};

inline constexpr Type kDefaultType = Type::Granite;

inline constexpr uint8_t to_index(Type type) {
    return static_cast<uint8_t>(type);
}

inline constexpr Type from_index(uint8_t index) {
    return index < kProperties.size() ? static_cast<Type>(index) : kDefaultType;
}

inline constexpr const Properties& properties(Type type) {
    return kProperties[to_index(type)];
}

inline constexpr const Properties& properties_from_index(uint8_t index) {
    return properties(from_index(index));
}

} // namespace platec::material

#endif
