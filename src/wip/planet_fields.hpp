#pragma once

// Read access to the world catalogue (world_fields.hpp) over planet's current
// storage. Used by run bundles and by the scheduler's access audit; it never
// changes the world.

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "world/field_encoding.hpp"

class planet;

using PlanetFieldData = worldmodel::FieldData;

// Copies one field in row-major order with its catalogue element type.
PlanetFieldData readplanetfield(const planet& world, worldmodel::FieldId id);

// Order-sensitive 64-bit hash of a field's shape and contents.
std::uint64_t hashplanetfield(const planet& world, worldmodel::FieldId id);

// Entity tables as text: CSV, or JSON for the mineral budget, or hexadecimal
// for the serialized weather state. The extension says which.
std::string planettabletext(const planet& world, worldmodel::TableId id, std::string* extension = nullptr);
std::uint64_t hashplanettable(const planet& world, worldmodel::TableId id);

std::uint64_t fnv1a64(const void* data, std::size_t size, std::uint64_t hash = 1469598103934665603ULL);
